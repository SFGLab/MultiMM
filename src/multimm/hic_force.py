"""
hic_force.py — Hi-C contact-probability force for MultiMM / OpenMM

Three layers: pre-processing (symmetrize_and_clean, diagonal_normalize,
resize_matrix; pure numpy), the OpenMM force builder
(build_crossentropy_force -> sparse CustomBondForce), and the public
entry point (build_hic_force: cleans -> resizes -> normalises -> builds).

Method: models contact *probability* as a monotone-decreasing function of
3-D distance, P_ij(r), and minimises binary cross-entropy between P_ij and
the observed contact frequency c_ij:

    L = -Σ_{i<j} [ c_ij · log P_ij(r_ij) + (1 - c_ij) · log(1 - P_ij(r_ij)) ]

F_i = -∇_i L, so the residual (c_ij - P_ij) is self-regulating and vanishes
once the structure matches the data. P_ij(r) is pluggable (`kernel=`):

  * "gaussian"    P = exp(-r² / (2σ²))                         (default)
  * "power_law"   P = 1 / (1 + (r / r_c)^α)                    (sigmoid; also "sigmoid")
  * "exponential" P = exp(-r / r_c)
  * "erfc"        P = ½ erfc((r - r_c) / (√2 σ_s))             (soft step)
  * "rouse"       P = erfc(r / √(2 s b²))                      (s = separation in beads,
                                                                  b = Kuhn length)

Each pair's loss contribution is also scaled by w_ij = c_ij^β
(`weight_power`, default 1), so weak-evidence pairs exert proportionally
small force instead of a hard constraint once past `threshold` (which
mostly controls sparsity/performance, not force hardness).

Usage (model.py):
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


def get_kernel_p_func(
    kernel      : str,
    r_comp      : float,
    alpha       : float           = 3.0,
    sigma       : Optional[float] = None,
    sigma_s     : Optional[float] = None,
    kuhn_length : Optional[float] = None,
):
    """
    Public helper so other modules (notably `validation.py`) can compute
    the exact same distance→contact-probability kernel P(r) used to build
    the Hi-C force, instead of a generic/unrelated proxy such as 1/(d+ε).

    Resolves the same None-fallback defaults as `build_crossentropy_force`
    (sigma → r_comp, sigma_s → 0.3*r_comp, kuhn_length → r_comp), so a
    caller only needs the same parameters that were used to build the
    force itself.

    Returns
    -------
    p_func    : callable(r, sep=None) -> P, vectorised over numpy arrays
                (r and, for "rouse", sep may be full (N, N) matrices).
    needs_sep : bool — True only for "rouse" (needs a genomic-separation
                matrix sep[i, j] = |i - j| passed as the second argument).
    """
    kernel = kernel.lower()
    if kernel not in VALID_KERNELS:
        raise ValueError(f"Unknown Hi-C kernel '{kernel}'. Valid options: {VALID_KERNELS}")

    sigma       = r_comp if sigma is None else sigma
    sigma_s     = 0.3 * r_comp if sigma_s is None else sigma_s
    kuhn_length = r_comp if kuhn_length is None else kuhn_length

    _, _, needs_sep, p_func, _ = _kernel_spec(kernel, r_comp, alpha, sigma, sigma_s, kuhn_length)
    return p_func, needs_sep


def auto_contact_scale(
    coords: np.ndarray,
    percentile: float = 10.0,
    n_samples: int = 20000,
    seed: int = 0,
) -> float:
    """Estimate a contact length scale from a structure's own pairwise-distance
    distribution, instead of reusing the Hi-C force's microscopic
    ``r_comp``/``sigma`` (tuned for bead-overlap, so it underflows to ~0 for
    any pair more than a few bead-spacings apart and gives an uninformative,
    near-diagonal-only contact map).

    Samples random bead pairs from the 3-D structure and returns the given
    percentile (default 10th, i.e. a "frequently-contacting" scale) of their
    pairwise distances, keeping the kernel informative across the whole
    distance range rather than saturating to zero almost everywhere.

    coords: (N, 3) structure coordinates (nm). percentile: lower biases
    tighter/diagonal-focused, higher spreads contacts further out.
    n_samples/seed: sampling size and RNG seed. Returns the calibrated
    length scale in the same units as coords.
    """
    coords = np.asarray(coords)
    N = coords.shape[0]
    if N < 2:
        return 1.0
    rng = np.random.default_rng(seed)
    # Oversample, then drop any accidental i == j pairs and trim to n_samples.
    draw = max(n_samples * 2, 64)
    i_idx = rng.integers(0, N, size=draw)
    j_idx = rng.integers(0, N, size=draw)
    keep = i_idx != j_idx
    i_idx, j_idx = i_idx[keep], j_idx[keep]
    n_use = min(n_samples, i_idx.size)
    if n_use == 0:
        return 1.0
    i_idx, j_idx = i_idx[:n_use], j_idx[:n_use]
    d = np.linalg.norm(coords[i_idx] - coords[j_idx], axis=1)
    d = d[np.isfinite(d) & (d > 0)]
    if d.size == 0:
        return 1.0
    return float(np.percentile(d, percentile))


def _oe_normalize_matrix(H: np.ndarray) -> np.ndarray:
    """Observed/Expected enrichment, floored at background: divide each
    diagonal by its mean, then subtract 1 and clip at 0.

    Used only when `HIC_FORCE_OE=True`, as the cross-entropy loss's c_ij
    target. OE=1 means "exactly the expected distance-decay background" —
    neither enriched nor depleted — so pairs at or below that baseline
    (OE <= 1, matching "towards -1" on the displayed log2(O/E) scale) get
    c_ij = 0 exactly, which makes their force weight w_ij = c_ij^weight_power
    exactly 0 too: the cross-entropy term vanishes regardless of current
    distance, i.e. genuinely no attraction *or* repulsion ("loose"), rather
    than a weak two-sided pull toward some nonzero background intensity.
    Only the EXCESS enrichment above background (OE > 1) becomes a positive
    c_ij, so only genuinely enriched pairs are pulled together — up to
    c_ij = 1 (the most-enriched pair after re-normalising by max in
    `build_crossentropy_force`), whose target P -> 1 drives it toward the
    smallest distance the kernel/excluded-volume force allow.
    """
    N = H.shape[0]
    oe = np.zeros_like(H, dtype=np.float64)
    for k in range(N):
        diag = np.diag(H, k)
        mean_k = diag.mean()
        oe_diag = diag / mean_k if mean_k > 0 else diag
        idx = np.arange(N - k)
        oe[idx, idx + k] = oe_diag
        oe[idx + k, idx] = oe_diag
    return np.maximum(oe - 1.0, 0.0)


def preprocess_hic_matrix(
    H_raw: np.ndarray,
    N_beads: int,
    already_balanced: bool = False,
    oe_normalize: bool = False,
) -> np.ndarray:
    """Run the same clean -> resize -> balance -> (optional) OE-floor
    pipeline `build_hic_force` uses, without building the OpenMM force —
    so diagnostics can see exactly the c_ij values the force actually
    targets (see `_oe_normalize_matrix` for what `oe_normalize` does).
    """
    H = symmetrize_and_clean(H_raw)
    if H.shape[0] != N_beads:
        H = resize_matrix(H, N_beads)
    if not already_balanced:
        H = diagonal_normalize(H)
    else:
        log.info("diagonal_normalize: skipped (already_balanced=True)")
    if oe_normalize:
        H = _oe_normalize_matrix(H)
        log.info("  OE normalisation applied — c_ij now targets enrichment "
                 "above background, floored at 0")
    return H


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
    return_controller: bool        = False,
    noise_seed   : int             = 0,
):
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

    Soft/proportional weighting: w_ij = c_ij ** weight_power (default 1), so
    force scales continuously with observed contact strength rather than
    acting as a hard constraint. Only pairs with C[i,j] >= threshold are
    built into bonds at all (sparsity/performance), but within that set the
    force strength still scales with c_ij via w_ij.

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
    return_controller : if True, also return a `HiCNoiseController` that
                   can periodically redraw each bond's c_ij with fresh
                   noise around its original value (see that class).
    noise_seed   : RNG seed for the returned controller.

    Returns
    -------
    force : mm.CustomBondForce, or (force, HiCNoiseController) if
            return_controller=True.
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

    if not return_controller:
        return force

    controller = HiCNoiseController(
        force, rows, cols, c_vals,
        needs_sep=needs_sep, msd_vals=msd_vals if needs_sep else None,
        seed=noise_seed,
    )
    return force, controller


class HiCNoiseController:
    """Periodically redraws each Hi-C bond's c_ij around its ORIGINAL,
    data-derived value with fresh Gaussian noise, instead of letting it
    drift — every perturbation is anchored back to the real contact
    strength, not the previous (already-noisy) one, so the noise never
    runs away; it only lets different contacts take turns pulling
    strongest, nudging the structure to explore nearby configurations
    instead of settling into exactly one attractor.

    Call :meth:`resample` periodically (e.g. once per saved MD frame) with
    the running `Context` and the desired noise intensity (std-dev of the
    noise added to c_ij, same [0, 1] scale; 0 disables it).
    """

    def __init__(self, force, rows, cols, base_c_vals, needs_sep=False,
                 msd_vals=None, seed=0):
        self.force = force
        self.rows = np.asarray(rows)
        self.cols = np.asarray(cols)
        self.base_c_vals = np.asarray(base_c_vals, dtype=np.float64)
        self.needs_sep = needs_sep
        self.msd_vals = msd_vals
        self.rng = np.random.default_rng(seed)

    def resample(self, context, intensity: float) -> None:
        if intensity <= 0 or self.base_c_vals.size == 0:
            return
        noisy = np.clip(
            self.base_c_vals + self.rng.normal(0.0, intensity, size=self.base_c_vals.shape),
            0.0, 1.0,
        )
        if self.needs_sep:
            for idx in range(noisy.shape[0]):
                self.force.setBondParameters(
                    idx, int(self.rows[idx]), int(self.cols[idx]),
                    [float(noisy[idx]), float(self.msd_vals[idx])],
                )
        else:
            for idx in range(noisy.shape[0]):
                self.force.setBondParameters(
                    idx, int(self.rows[idx]), int(self.cols[idx]), [float(noisy[idx])],
                )
        self.force.updateParametersInContext(context)


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
    oe_normalize     : bool            = False,
    return_controller: bool            = False,
    noise_seed       : int             = 0,
):
    """
    Full pipeline: raw Hi-C array → sparse contact-probability force.

    Pipeline
    --------
    H_raw
      │  symmetrize_and_clean()
      │  resize_matrix()              ← to N_beads × N_beads
      │  diagonal_normalize()         ← skip if already_balanced=True
      │  OE normalisation             ← only if oe_normalize=True
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
    oe_normalize     : True → target enrichment-above-background (OE - 1,
                       floored at 0 — see `_oe_normalize_matrix`) instead of
                       raw contact frequency, so pairs at/below the expected
                       distance-decay baseline exert exactly zero force
                       ("loose"), and only genuinely enriched pairs attract,
                       scaling up to the smallest allowed distance for the
                       most-enriched pair. Default False.
    return_controller : if True, also return a `HiCNoiseController` for
                       periodically re-noising c_ij around its original
                       value (see that class's docstring).
    noise_seed       : RNG seed for the returned controller.

    Returns
    -------
    force : mm.CustomBondForce, or (force, HiCNoiseController) if
            return_controller=True.

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
            ("OE normalize",  oe_normalize),
            ("Input shape",   str(H_raw.shape)),
        ],
        title="Hi-C Force — build (contact-probability cross-entropy)",
        log_fn=log.info,
    )

    # ── Layer 1: pre-processing ──────────────────────────────────────────────
    H = preprocess_hic_matrix(
        H_raw, N_beads, already_balanced=already_balanced, oe_normalize=oe_normalize,
    )

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
        return_controller=return_controller,
        noise_seed=noise_seed,
    )

    log.info("build_hic_force: done (kernel=%s, group %d)", kernel, force_group)
    log.info("═" * 60)
    return force
