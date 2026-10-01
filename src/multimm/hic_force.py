"""
hic_force.py  —  Hi-C contact-probability force for MultiMM / OpenMM
=====================================================================

Organised in three clearly-separated layers:

  Layer 1 — Pre-processing   (pure numpy, no OpenMM)
      symmetrize_and_clean()
      diagonal_normalize()
      resize_matrix()

  Layer 2 — OpenMM force builder
      build_crossentropy_force()  → sparse CustomBondForce

  Layer 3 — Public entry point
      build_hic_force()  — cleans → resizes → normalises → builds force

All layers emit structured log messages through the module logger
(see logger.py for setup).

Method
------
Rather than converting contacts into fixed target distances, the force
models the contact *probability* directly as a monotone-decreasing
function of 3-D distance, P_ij(r), and minimises the binary
cross-entropy (negative log-likelihood) between the model probability
and the observed, row-normalised Hi-C contact frequency c_ij ∈ [0, 1]:

    L = -Σ_{i<j} [ c_ij · log P_ij(r_ij) + (1 - c_ij) · log(1 - P_ij(r_ij)) ]

The force is F_i = -∇_i L.  Because the energy is built from P_ij
directly, the residual (c_ij - P_ij) is self-regulating: it vanishes
once the simulated structure reproduces the data, in either direction
(push apart if too close, pull together if too far).

P_ij(r) is pluggable (`kernel=` argument).  Five kernels are provided;
`gaussian` is the default:

  * "gaussian"    P = exp(-r² / (2σ²))                         (default)
  * "power_law"   P = 1 / (1 + (r / r_c)^α)                    (sigmoid; also "sigmoid")
  * "exponential" P = exp(-r / r_c)
  * "erfc"        P = ½ erfc((r - r_c) / (√2 σ_s))             (soft step)
  * "rouse"       P = erfc(r / √(2 s b²))                      (s = |i-j| in beads,
                                                                  b = Kuhn length;
                                                                  separation-aware)

Soft, proportional weighting
-----------------------------
A hard distance cutoff is not the only thing kept soft here: every
pair's contribution to the loss (and hence its force) is additionally
scaled by a weight w_ij = c_ij^β (β = `weight_power`, default 1 ⇒ force
strictly proportional to the observed contact strength).  This means
pairs with c_ij close to the inclusion `threshold` (i.e. almost no
evidence of contact) exert a correspondingly tiny force — both the
attractive and the repulsive branch of the loss are damped — instead of
acting as a hard constraint the moment they clear the threshold.  Only
pairs with strong contact evidence behave close to a hard restraint.
`threshold` therefore mostly controls sparsity (how many bonds are
built, for performance), not force "hardness" — that is governed
continuously by c_ij itself.

Usage in model.py
------------------
    from hic_force import build_hic_force

    force = build_hic_force(args.hic_matrix, N_beads=self.N,
                             r_comp=self.r_comp, kernel="gaussian")
    self.system.addForce(force)
"""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np
from scipy.ndimage import zoom
from scipy.special import erfc as _erfc

try:
    import openmm as mm
except ImportError:                     # legacy simtk namespace
    from simtk import openmm as mm

# ── module-level logger ───────────────────────────────────────────────────────
log = logging.getLogger(__name__)       # hic_force

from .logger import log_table  # noqa: E402 — after log is set up


# ═════════════════════════════════════════════════════════════════════════════
# Layer 1 — Pre-processing
# ═════════════════════════════════════════════════════════════════════════════

def symmetrize_and_clean(H: np.ndarray) -> np.ndarray:
    """
    Symmetrise *H* and remove pathological values.

    Steps
    -----
    1. Cast to float64.
    2. Average H and H^T  →  enforces exact symmetry.
    3. Replace NaN / ±Inf with 0.
    4. Clip negatives to 0 (counts cannot be negative).
    5. Zero the main diagonal (self-contacts are uninformative).

    Parameters
    ----------
    H : (N, N) array_like
        Raw Hi-C counts or pre-normalised matrix.

    Returns
    -------
    H_clean : (N, N) ndarray, float64, symmetric, non-negative
    """
    H = np.array(H, dtype=np.float64)

    if H.ndim != 2 or H.shape[0] != H.shape[1]:
        raise ValueError(
            f"Hi-C matrix must be square 2-D ndarray, got shape {H.shape}"
        )

    N = H.shape[0]
    log.info("symmetrize_and_clean: input shape (%d, %d)", N, N)

    H = 0.5 * (H + H.T)
    bad = ~np.isfinite(H)
    if bad.any():
        log.warning("  %d non-finite entries replaced with 0", bad.sum())
    H = np.where(bad, 0.0, H)
    H = np.maximum(H, 0.0)
    np.fill_diagonal(H, 0.0)

    log.debug("  non-zero entries: %d / %d", (H > 0).sum(), N * N)
    return H


def diagonal_normalize(H: np.ndarray, n_iter: int = 50) -> np.ndarray:
    """
    Knight-Ruiz iterative row/column balancing.

    Iterates   H ← D⁻¹ H D⁻¹   until marginals are approximately 1.
    Rows/columns with zero sum are excluded from scaling (masked out).

    Parameters
    ----------
    H      : (N, N) ndarray — cleaned, symmetric, non-negative
    n_iter : number of balancing iterations (default 50)

    Returns
    -------
    H_bal : (N, N) ndarray — doubly stochastic up to a global constant
    """
    log.info("diagonal_normalize: %d iterations", n_iter)
    H = H.copy()

    for it in range(n_iter):
        row_sums = H.sum(axis=1)
        scale    = np.where(row_sums > 0, 1.0 / np.sqrt(row_sums), 0.0)
        H        = scale[:, None] * H * scale[None, :]
        H        = np.where(np.isfinite(H), H, 0.0)

    H = np.maximum(H, 0.0)
    np.fill_diagonal(H, 0.0)

    residual = np.abs(H.sum(axis=1) - 1.0)
    log.debug("  balancing residual  max=%.3e  mean=%.3e",
              residual.max(), residual.mean())
    return H


def resize_matrix(H: np.ndarray, N_target: int) -> np.ndarray:
    """
    Bilinear resampling of a square Hi-C matrix to a new resolution.

    Parameters
    ----------
    H        : (N_src, N_src) ndarray
    N_target : desired side length

    Returns
    -------
    H_resized : (N_target, N_target) ndarray, symmetric
    """
    N_src = H.shape[0]
    if N_src == N_target:
        return H.copy()

    log.info("resize_matrix: %d → %d beads (factor %.3f)",
             N_src, N_target, N_target / N_src)

    H_res = zoom(H, N_target / N_src, order=1)

    # zoom can produce N ± 1 due to floating-point rounding
    if H_res.shape[0] != N_target:
        tmp   = np.zeros((N_target, N_target), dtype=H_res.dtype)
        s     = min(H_res.shape[0], N_target)
        tmp[:s, :s] = H_res[:s, :s]
        H_res = tmp

    H_res = 0.5 * (H_res + H_res.T)
    np.fill_diagonal(H_res, 0.0)
    return H_res


# ═════════════════════════════════════════════════════════════════════════════
# Layer 2 — Contact-probability kernels
# ═════════════════════════════════════════════════════════════════════════════

#: Names accepted by `kernel=`.  "sigmoid" is kept as an alias of "power_law"
#: for backward compatibility with earlier MultiMM configs.
VALID_KERNELS = ("gaussian", "power_law", "sigmoid", "exponential", "erfc", "rouse")


def _kernel_spec(
    kernel      : str,
    r_comp      : float,
    alpha       : float,
    sigma       : float,
    sigma_s     : float,
    kuhn_length : float,
):
    """
    Return the OpenMM expression fragment for P(r), the extra global
    parameters it needs, whether it needs a per-bond genomic-separation
    parameter, and a plain-numpy callable P(r, sep) used only for
    calibration logging (never executed inside OpenMM).

    Returns
    -------
    p_expr      : str   — OpenMM/Lepton fragment "P = ...;"
    globals_    : dict  — {name: value} additional global parameters
    needs_sep   : bool  — True only for "rouse" (needs per-bond msd_ij)
    p_func      : callable(r, sep=None) -> P, for diagnostics only
    ref_r       : float — characteristic length scale, for diagnostics only
    """
    kernel = kernel.lower()

    if kernel in ("power_law", "sigmoid"):
        p_expr = "P = 1 / (1 + (r / hic_rc)^hic_alpha);"
        globals_ = {"hic_rc": r_comp, "hic_alpha": alpha}
        needs_sep = False

        def p_func(r, sep=None):
            return 1.0 / (1.0 + (r / r_comp) ** alpha)

        ref_r = r_comp

    elif kernel == "gaussian":
        p_expr = "P = exp(-(r^2) / (2 * hic_sigma^2));"
        globals_ = {"hic_sigma": sigma}
        needs_sep = False

        def p_func(r, sep=None):
            return np.exp(-(r ** 2) / (2.0 * sigma ** 2))

        ref_r = sigma

    elif kernel == "exponential":
        p_expr = "P = exp(-r / hic_rc);"
        globals_ = {"hic_rc": r_comp}
        needs_sep = False

        def p_func(r, sep=None):
            return np.exp(-r / r_comp)

        ref_r = r_comp

    elif kernel == "erfc":
        p_expr = "P = 0.5 * erfc((r - hic_rc) / (sqrt(2) * hic_sigmas));"
        globals_ = {"hic_rc": r_comp, "hic_sigmas": sigma_s}
        needs_sep = False

        def p_func(r, sep=None):
            return 0.5 * _erfc((r - r_comp) / (np.sqrt(2.0) * sigma_s))

        ref_r = r_comp

    elif kernel == "rouse":
        # msd_ij = <r^2(s)> = s * b^2 is precomputed per-bond (s differs per
        # pair), so the OpenMM expression only needs the per-bond parameter.
        p_expr = "P = erfc(r / sqrt(2 * msd_ij));"
        globals_ = {}
        needs_sep = True

        def p_func(r, sep=1):
            msd = np.maximum(sep, 1) * kuhn_length ** 2
            return _erfc(r / np.sqrt(2.0 * msd))

        ref_r = kuhn_length

    else:
        raise ValueError(
            f"Unknown Hi-C kernel '{kernel}'. Valid options: {VALID_KERNELS}"
        )

    return p_expr, globals_, needs_sep, p_func, ref_r


def build_crossentropy_force(
    C            : np.ndarray,
    r_comp       : float,
    kernel       : str             = "gaussian",
    threshold    : float           = 0.01,
    alpha        : float           = 3.0,
    sigma        : Optional[float] = None,
    sigma_s      : Optional[float] = None,
    kuhn_length  : Optional[float] = None,
    k_scale      : float           = 1.0,
    weight_power : float           = 1.0,
    force_group  : int             = 1,
) -> mm.CustomBondForce:
    """
    Construct a sparse CustomBondForce from the binary cross-entropy loss
    between observed contact frequencies and a model contact probability
    P(r), with a pluggable distance→probability kernel.

    Physics
    -------
    Per-pair cross-entropy energy (weighted, see below):

        U_ij(r) = -k_scale · w_ij · [c_ij·log(P+ε) + (1-c_ij)·log(1-P+ε)]

    where P = P(r) comes from the chosen `kernel` (see module docstring
    for the five options) and the force is F_i = -∇_i U summed over all
    bonds, which OpenMM derives automatically from the symbolic
    expression.

    Soft / proportional weighting
    ------------------------------
    w_ij = c_ij ** weight_power.  With the default weight_power=1, the
    force on a pair is directly proportional to how strong its observed
    contact is: pairs just above `threshold` (almost no contact
    evidence) contribute a correspondingly tiny force, while only
    well-supported contacts (c_ij close to 1) behave like a firm
    restraint. This keeps the force "relaxed" everywhere the data don't
    support a contact, rather than imposing a hard constraint on every
    pair that merely clears the sparsity threshold.

    Only pairs with  C[i, j] ≥ threshold  are built into bonds at all
    (for performance — O(M) bonds, M ≪ N²); within that set, the actual
    force strength still scales continuously with c_ij via w_ij.

    Parameters
    ----------
    C            : (N_beads, N_beads) ndarray — normalised contact matrix
                   with values in [0, 1].  Typically output of
                   diagonal_normalize() divided by its maximum.
    r_comp       : characteristic contact distance [nm], used as r_c for
                   "power_law"/"sigmoid"/"exponential"/"erfc".  Ignored
                   by "gaussian" (use `sigma`) and "rouse" (use
                   `kuhn_length`), except as their fallback default.
    kernel       : one of VALID_KERNELS. Default "gaussian".
    threshold    : minimum C[i,j] to add a bond (sparsity control only).
    alpha        : "power_law" sigmoid steepness (typical: 2–4).
    sigma        : "gaussian" width [nm]. Defaults to r_comp.
    sigma_s      : "erfc" softening width [nm]. Defaults to 0.3 * r_comp.
    kuhn_length  : "rouse" Kuhn (statistical segment) length [nm].
                   Defaults to r_comp.
    k_scale      : global energy scale [kJ/mol].
    weight_power : exponent β in w_ij = c_ij^β (default 1.0 ⇒ force
                   strictly proportional to observed contact strength).
    force_group  : OpenMM force-group index.

    Returns
    -------
    force : mm.CustomBondForce
    """
    kernel = kernel.lower()
    if kernel not in VALID_KERNELS:
        raise ValueError(f"Unknown Hi-C kernel '{kernel}'. Valid options: {VALID_KERNELS}")

    sigma       = r_comp if sigma is None else sigma
    sigma_s     = 0.3 * r_comp if sigma_s is None else sigma_s
    kuhn_length = r_comp if kuhn_length is None else kuhn_length

    N_beads = C.shape[0]
    log.info(
        "build_crossentropy_force: N=%d, kernel=%s, r_comp=%.3f nm, "
        "threshold=%.4f, weight_power=%.2f",
        N_beads, kernel, r_comp, threshold, weight_power,
    )

    # normalise to [0, 1]
    C_max = C.max()
    if C_max > 0:
        C = C / C_max
    else:
        log.warning("  contact matrix is all-zero — no bonds will be added")

    # upper-triangle pairs above threshold (exclude diagonal by construction)
    rows, cols = np.triu_indices(N_beads, k=1)
    mask       = C[rows, cols] >= threshold
    rows, cols = rows[mask], cols[mask]
    c_vals     = C[rows, cols]
    sep_vals   = (cols - rows).astype(np.float64)   # genomic separation, beads

    n_bonds = len(rows)
    log.info("  pairs above threshold: %d  (%.2f%% of upper triangle)",
             n_bonds, 100.0 * n_bonds / (N_beads * (N_beads - 1) / 2))

    if n_bonds == 0:
        log.warning("  no bonds added — consider lowering threshold or "
                    "checking matrix normalisation")

    p_expr, kernel_globals, needs_sep, p_func, ref_r = _kernel_spec(
        kernel, r_comp, alpha, sigma, sigma_s, kuhn_length,
    )

    # ── calibration diagnostics (approximate — numeric, kernel-agnostic) ─────
    w_mean = float(np.mean(c_vals ** weight_power)) if n_bonds else 0.0
    dr     = max(ref_r * 1e-3, 1e-6)
    sep_ref = float(np.median(sep_vals)) if n_bonds else 1.0
    P_lo   = p_func(ref_r - dr, sep_ref)
    P_hi   = p_func(ref_r + dr, sep_ref)
    P_ref  = p_func(ref_r, sep_ref)
    dP_dr  = (P_hi - P_lo) / (2.0 * dr)
    denom  = max(P_ref * (1.0 - P_ref), 1e-6)
    f_typical = k_scale * w_mean * abs(dP_dr) / denom      # kJ/mol/nm (approx)
    bonds_per_bead = 2.0 * n_bonds / N_beads
    log.info(
        "  calibration (kernel=%s): k_scale=%.2f kJ/mol, <c_ij>=%.3f, "
        "<w_ij>=%.3f → F/bond@r≈ref≈%.1f kJ/mol/nm, "
        "%.1f bonds/bead → total~%.0f kJ/mol/nm per bead",
        kernel, k_scale, float(c_vals.mean()) if n_bonds else 0.0, w_mean,
        f_typical, bonds_per_bead, f_typical * bonds_per_bead,
    )
    if f_typical * bonds_per_bead < 25.0:
        log.warning(
            "  Hi-C cross-entropy force may be too weak to drive folding "
            "(total force/bead < 10 kT/nm near the kernel's characteristic "
            "scale).  Consider increasing HIC_K_SCALE (current %.2f kJ/mol); "
            "5–20 kJ/mol is typical.",
            k_scale,
        )
    elif k_scale > 30.0:
        log.warning(
            "  HIC_K_SCALE=%.2f kJ/mol is very high.  High values make the "
            "well around strongly-supported contacts very stiff, shrinking "
            "thermal fluctuations and potentially freezing MD sampling.  "
            "5–20 kJ/mol is typical; weight_power (currently %.2f) already "
            "softens weak contacts independently of k_scale.",
            k_scale, weight_power,
        )

    # OpenMM expression  — note: OpenMM's parser uses `x^y`, not `pow(x, y)`
    # Order matches existing MultiMM style: main energy term first, then
    # the auxiliary variables it depends on (w_ij, then P).
    expression = (
        "-hic_k * w_ij * ( c_ij * log(P + hic_eps) + (1 - c_ij) * log(1 - P + hic_eps) );"
        "w_ij = c_ij^hic_wpow;"
        f"{p_expr}"
    )

    force = mm.CustomBondForce(expression)
    force.addGlobalParameter("hic_k",    k_scale)
    force.addGlobalParameter("hic_eps",  1e-8)
    force.addGlobalParameter("hic_wpow", weight_power)
    for name, val in kernel_globals.items():
        force.addGlobalParameter(name, val)

    force.addPerBondParameter("c_ij")
    if needs_sep:
        force.addPerBondParameter("msd_ij")   # <r^2(s)> = s * b^2, precomputed

    if needs_sep:
        msd_vals = sep_vals * (kuhn_length ** 2)
        for idx in range(n_bonds):
            force.addBond(int(rows[idx]), int(cols[idx]),
                          [float(c_vals[idx]), float(msd_vals[idx])])
    else:
        for idx in range(n_bonds):
            force.addBond(int(rows[idx]), int(cols[idx]), [float(c_vals[idx])])

    force.setForceGroup(force_group)
    log.info("  Hi-C contact-probability force ready  (%d bonds, kernel=%s)",
             n_bonds, kernel)
    return force


# ═════════════════════════════════════════════════════════════════════════════
# Layer 3 — Public entry point
# ═════════════════════════════════════════════════════════════════════════════

def build_hic_force(
    H_raw            : np.ndarray,
    N_beads          : int,
    r_comp           : float,
    kernel           : str             = "gaussian",
    threshold        : float           = 0.01,
    alpha            : float           = 3.0,
    sigma            : Optional[float] = None,
    sigma_s          : Optional[float] = None,
    kuhn_length      : Optional[float] = None,
    k_scale          : float           = 1.0,
    weight_power     : float           = 1.0,
    force_group      : int             = 1,
    already_balanced : bool            = False,
) -> mm.CustomBondForce:
    """
    Full pipeline: raw Hi-C array → sparse contact-probability force.

    Pipeline
    --------
    H_raw
      │  symmetrize_and_clean()
      │  resize_matrix()              ← to N_beads × N_beads
      │  diagonal_normalize()         ← skip if already_balanced=True
      └─ build_crossentropy_force()   → sparse CustomBondForce

    Parameters
    ----------
    H_raw            : (M, M) ndarray — raw Hi-C counts (any resolution).
                       Resampled to N_beads × N_beads if M ≠ N_beads.
    N_beads          : number of simulation beads.
    r_comp           : characteristic contact distance [nm]. Used
                       directly by "power_law"/"sigmoid"/"exponential"/
                       "erfc" as r_c, and as the fallback default for
                       `sigma` / `kuhn_length` when those are not given.
    kernel           : distance→probability model. One of VALID_KERNELS.
                       Default "gaussian".
    threshold        : minimum normalised c_ij to add a bond (sparsity
                       only — see build_crossentropy_force doc for how
                       force magnitude still scales with c_ij).
    alpha            : "power_law" sigmoid steepness (2–4 typical).
    sigma            : "gaussian" width [nm] (defaults to r_comp).
    sigma_s          : "erfc" softening width [nm] (defaults to 0.3*r_comp).
    kuhn_length      : "rouse" Kuhn length [nm] (defaults to r_comp).
    k_scale          : global energy scale [kJ/mol].
                       Recommended range: 5–20 kJ/mol.
                       Values > 30 kJ/mol risk freezing MD sampling.
    weight_power     : β in w_ij = c_ij^β (default 1.0 — force
                       proportional to observed contact strength).
    force_group      : OpenMM force-group index.
    already_balanced : True → skip diagonal_normalize().

    Returns
    -------
    force : mm.CustomBondForce

    Examples
    --------
    >>> force = build_hic_force(hic_array, N_beads=500, r_comp=0.15,
    ...                         kernel="gaussian", k_scale=10.0)
    >>> system.addForce(force)
    """
    log_table(
        [
            ("N beads",       N_beads),
            ("Kernel",        kernel),
            ("r_comp",        f"{r_comp:.3f} nm"),
            ("threshold",     threshold),
            ("alpha",         alpha if kernel in ("power_law", "sigmoid") else "n/a"),
            ("k_scale",       f"{k_scale:.3f} kJ/mol"),
            ("weight_power",  weight_power),
            ("Force group",   force_group),
            ("Input shape",   str(H_raw.shape)),
        ],
        title="Hi-C Force — build (contact-probability cross-entropy)",
        log_fn=log.info,
    )

    # ── Layer 1: pre-processing ──────────────────────────────────────────────
    H = symmetrize_and_clean(H_raw)

    if H.shape[0] != N_beads:
        H = resize_matrix(H, N_beads)

    if not already_balanced:
        H = diagonal_normalize(H)
    else:
        log.info("diagonal_normalize: skipped (already_balanced=True)")

    # ── Layer 2: build force ─────────────────────────────────────────────────
    force = build_crossentropy_force(
        H, r_comp,
        kernel=kernel,
        threshold=threshold,
        alpha=alpha,
        sigma=sigma,
        sigma_s=sigma_s,
        kuhn_length=kuhn_length,
        k_scale=k_scale,
        weight_power=weight_power,
        force_group=force_group,
    )

    log.info("build_hic_force: done (kernel=%s, group %d)", kernel, force_group)
    log.info("═" * 60)
    return force
