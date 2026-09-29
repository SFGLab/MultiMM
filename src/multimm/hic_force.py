"""
hic_force.py  —  Hi-C cross-entropy force for MultiMM / OpenMM
===============================================================

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

Usage in model.py
-----------------
    from hic_force import build_hic_force

    force = build_hic_force(args.hic_matrix, N_beads=self.N,
                            r_comp=self.r_comp)
    self.system.addForce(force)
"""

from __future__ import annotations

import logging

import numpy as np
from scipy.ndimage import zoom

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
# Layer 2 — Cross-entropy force builder
# ═════════════════════════════════════════════════════════════════════════════

def build_crossentropy_force(
    C           : np.ndarray,
    r_comp      : float,
    threshold   : float = 0.01,
    alpha       : float = 3.0,
    k_scale     : float = 1.0,
    force_group : int   = 1,
) -> mm.CustomBondForce:
    """
    Construct a sparse CustomBondForce from the binary cross-entropy loss.

    Physics
    -------
    Contact probability model:

        P(r) = 1 / (1 + (r / r_c)^α)

    Per-pair cross-entropy energy:

        U_{ij}(r) = −k_scale · [c_ij · log(P + ε) + (1 − c_ij) · log(1 − P + ε)]

    Resulting force on bead i (self-regulating residual):

        F_ij = (α / r) · (P − c_ij)  in direction r̂_{ij}

    When P > c_ij the beads are too close → push apart;
    when P < c_ij they are too far        → pull together.

    Only pairs with  C[i, j] ≥ threshold  are included, keeping the
    bond list sparse (O(M) bonds, M ≪ N²).

    Parameters
    ----------
    C           : (N_beads, N_beads) ndarray — normalised contact matrix
                  with values in [0, 1].  Typically output of
                  diagonal_normalize() divided by its maximum.
    r_comp      : contact radius r_c [nm]
    threshold   : minimum C[i,j] to add a bond
    alpha       : sigmoid steepness (typical: 2–4)
    k_scale     : global energy scale [kJ/mol]
    force_group : OpenMM force-group index

    Returns
    -------
    force : mm.CustomBondForce
    """
    N_beads = C.shape[0]
    log.info("build_crossentropy_force: N=%d, r_c=%.3f nm, α=%.1f, "
             "threshold=%.4f", N_beads, r_comp, alpha, threshold)

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

    n_bonds = len(rows)
    log.info("  pairs above threshold: %d  (%.2f%% of upper triangle)",
             n_bonds, 100.0 * n_bonds / (N_beads * (N_beads - 1) / 2))

    if n_bonds == 0:
        log.warning("  no bonds added — consider lowering threshold or "
                    "checking matrix normalisation")

    # ── calibration diagnostics ──────────────────────────────────────────────
    f_typical = k_scale * alpha * float(c_vals.mean()) / 1.0   # kJ/mol/nm
    bonds_per_bead = 2.0 * n_bonds / N_beads
    log.info(
        "  calibration: k_scale=%.2f kJ/mol, α=%.1f, "
        "<c_ij>=%.3f → F/bond@1nm≈%.1f kJ/mol/nm, "
        "%.1f bonds/bead → total~%.0f kJ/mol/nm per bead",
        k_scale, alpha, float(c_vals.mean()),
        f_typical, bonds_per_bead, f_typical * bonds_per_bead,
    )
    if f_typical * bonds_per_bead < 25.0:
        log.warning(
            "  Hi-C cross-entropy force may be too weak to drive folding "
            "(total force/bead < 10 kT/nm at 1 nm).  Consider increasing "
            "HIC_K_SCALE (current %.2f kJ/mol); 5–20 kJ/mol is typical.",
            k_scale,
        )
    elif k_scale > 30.0:
        log.warning(
            "  HIC_K_SCALE=%.2f kJ/mol is very high for crossentropy mode.  "
            "The effective spring constant at equilibrium is ~%.0f kJ/mol/nm², "
            "shrinking thermal fluctuations to <0.1 Å and freezing MD sampling.  "
            "For crossentropy, use HIC_K_SCALE in the 5–20 kJ/mol range.",
            k_scale,
            k_scale * alpha ** 2 * 0.25 / (r_comp ** 2),
        )

    # OpenMM expression  — note: OpenMM's parser uses `x^y`, not `pow(x, y)`
    expression = (
        "-hic_k * ( c_ij * log(P + hic_eps) + (1 - c_ij) * log(1 - P + hic_eps) );"
        "P = 1 / (1 + (r / hic_rc)^hic_alpha)"
    )

    force = mm.CustomBondForce(expression)
    force.addGlobalParameter("hic_k",     k_scale)
    force.addGlobalParameter("hic_rc",    r_comp)
    force.addGlobalParameter("hic_alpha", alpha)
    force.addGlobalParameter("hic_eps",   1e-8)
    force.addPerBondParameter("c_ij")

    for idx in range(n_bonds):
        force.addBond(int(rows[idx]), int(cols[idx]), [float(c_vals[idx])])

    force.setForceGroup(force_group)
    log.info("  CrossEntropy force ready  (%d bonds)", n_bonds)
    return force


# ═════════════════════════════════════════════════════════════════════════════
# Layer 3 — Public entry point
# ═════════════════════════════════════════════════════════════════════════════

def build_hic_force(
    H_raw            : np.ndarray,
    N_beads          : int,
    r_comp           : float,
    threshold        : float = 0.01,
    alpha            : float = 3.0,
    k_scale          : float = 1.0,
    force_group      : int   = 1,
    already_balanced : bool  = False,
) -> mm.CustomBondForce:
    """
    Full pipeline: raw Hi-C array → sparse CrossEntropy force.

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
    r_comp           : sigmoid contact radius r_c [nm].
    threshold        : minimum normalised c_ij to add a bond.
    alpha            : sigmoid steepness (2–4 typical).
    k_scale          : global energy scale [kJ/mol].
                       Recommended range: 5–20 kJ/mol.
                       Values > 30 kJ/mol freeze MD sampling (runtime warning).
    force_group      : OpenMM force-group index.
    already_balanced : True → skip diagonal_normalize().

    Returns
    -------
    force : mm.CustomBondForce

    Examples
    --------
    >>> force = build_hic_force(hic_array, N_beads=500, r_comp=0.15,
    ...                         k_scale=10.0)
    >>> system.addForce(force)
    """
    log_table(
        [
            ("N beads",       N_beads),
            ("r_comp",        f"{r_comp:.3f} nm"),
            ("threshold",     threshold),
            ("alpha",         alpha),
            ("k_scale",       f"{k_scale:.3f} kJ/mol"),
            ("Force group",   force_group),
            ("Input shape",   str(H_raw.shape)),
        ],
        title="Hi-C Force — build (crossentropy)",
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
        threshold=threshold,
        alpha=alpha,
        k_scale=k_scale,
        force_group=force_group,
    )

    log.info("build_hic_force: done (crossentropy, group %d)", force_group)
    log.info("═" * 60)
    return force
