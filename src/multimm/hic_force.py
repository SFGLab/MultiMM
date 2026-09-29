"""
hic_force.py  —  Hi-C-guided force fields for MultiMM / OpenMM
===============================================================

Organised in four clearly-separated layers:

  Layer 1 — Pre-processing   (pure numpy, no OpenMM)
      symmetrize_and_clean()
      diagonal_normalize()
      compute_oe_matrix()
      resize_matrix()

  Layer 2 — Decomposition    (pure numpy, no OpenMM)
      svd_decompose()

  Layer 3 — OpenMM force builders   (return mm.Force objects)
      build_svd_force()
      build_crossentropy_force()

  Layer 4 — Public entry point
      build_hic_force()  — cleans → resizes → normalises → decomposes
                           → dispatches to the right builder

All layers emit structured log messages through the module logger
(see logger.py for setup).

Usage in model.py
-----------------
    from hic_force import build_hic_force

    force = build_hic_force(args.hic_matrix, N_beads=self.N,
                            r_comp=self.r_beads * 3, mode='svd')
    self.system.addForce(force)
"""

from __future__ import annotations

import logging
from typing import Literal, Tuple

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


def compute_oe_matrix(H: np.ndarray, pseudocount: float = 1e-6) -> np.ndarray:
    """
    Compute the log Observed / Expected (O/E) matrix.

        OE[i, j] = log( H[i, j] / E[i, j] )

    where E[i, j] is the genome-wide mean contact frequency at genomic
    separation |i − j|, capturing the expected polymer distance-decay.

    Values in ℝ:  positive → attraction,  negative → repulsion.
    The main diagonal is set to 0.

    Parameters
    ----------
    H           : (N, N) ndarray — balanced contact matrix
    pseudocount : small constant added before the log to avoid log(0)

    Returns
    -------
    OE : (N, N) ndarray, float64
    """
    N = H.shape[0]
    log.info("compute_oe_matrix: N=%d, pseudocount=%.1e", N, pseudocount)

    H = H.copy() + pseudocount

    # vectorised per-diagonal means (expected contact at each separation)
    expected = np.ones(N, dtype=np.float64)
    for k in range(1, N):
        diag_vals  = np.diagonal(H, offset=k)
        expected[k] = diag_vals.mean() if diag_vals.size else 1.0

    idx = np.arange(N)
    sep = np.abs(idx[:, None] - idx[None, :])      # (N, N) separation matrix
    E   = expected[sep]                             # (N, N) expected matrix

    OE  = np.log(H / E)
    np.fill_diagonal(OE, 0.0)

    log.debug("  OE  min=%.3f  max=%.3f  mean=%.3f",
              OE.min(), OE.max(), OE.mean())
    return OE


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
# Layer 2 — SVD decomposition  (pure numpy)
# ═════════════════════════════════════════════════════════════════════════════

def svd_decompose(
    H_oe : np.ndarray,
    K    : int = 10,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Truncated eigendecomposition of the symmetric O/E matrix.

    Because H_oe is symmetric, SVD reduces to eigendecomposition:

        H_oe  ≈  Σ_{k=1}^{K}  λ_k · v_k · v_k^T

    The K components with the **largest |λ_k|** are retained.

    Per-particle scalars that encode both sign and magnitude:

        a_k(i) = sign(λ_k) · √|λ_k| · v_k(i)

    so that  a_k(i) · a_k(j) = sign(λ_k) · |λ_k| · v_k(i) · v_k(j),
    and the sum Σ_k a_k(i)·a_k(j) reconstructs H_oe[i,j] exactly for K=N.

    Parameters
    ----------
    H_oe : (N, N) ndarray — O/E matrix (symmetric)
    K    : number of eigenvectors to retain

    Returns
    -------
    lam  : (K,) eigenvalues, sorted by |λ| descending
    vecs : (N, K) eigenvectors (columns)
    A    : (N, K) per-particle parameter matrix  a_k(i)
    """
    N = H_oe.shape[0]
    K = min(K, N - 1)
    log.info("svd_decompose: N=%d, K=%d", N, K)

    eigvals, eigvecs = np.linalg.eigh(H_oe)       # ascending order
    order            = np.argsort(np.abs(eigvals))[::-1][:K]
    lam              = eigvals[order]              # (K,)
    vecs             = eigvecs[:, order]           # (N, K)

    signs  = np.sign(lam)
    scales = np.sqrt(np.abs(lam))
    A      = vecs * (signs * scales)[None, :]     # (N, K)

    # Normalise so that the maximum self-dot-product (= diagonal of rank-K approx)
    # equals 1.  This makes k_scale directly control the energy ceiling in kJ/mol.
    diag_vals = np.sum(A ** 2, axis=1)           # (N,)  A[i]·A[i] per bead
    max_diag  = diag_vals.max()
    explained = np.abs(lam).sum() / np.abs(eigvals).sum() * 100
    if max_diag > 0:
        A /= np.sqrt(max_diag)
        log_table(
            [
                ("N",                  N),
                ("K retained",         K),
                ("Variance explained", f"{explained:.1f}%"),
                ("A normalised by √",  f"{max_diag:.4g}"),
                ("Max self-dot-prod",  "1.0"),
                ("λ[:5]",             str(np.round(lam[:5], 4))),
            ],
            title="SVD decomposition",
            log_fn=log.info,
        )
    else:
        log.warning("  A matrix is all-zero after decomposition — Hi-C force will have no effect")
    return lam, vecs, A


# ─────────────────────────────────────────────────────────────────────────────
# Per-component σ helpers  (used by multi-scale SVD)
# ─────────────────────────────────────────────────────────────────────────────

def compute_sigma_eigenvalue(
    eigvals  : np.ndarray,
    sigma_max: float,
    beta     : float = 0.5,
) -> np.ndarray:
    """
    Assign per-component length scales from eigenvalue magnitudes.

        σ_k = σ_max · (|λ_k| / |λ_1|)^β

    Larger eigenvalues → larger spatial scale.  β = 0.5 is a sensible
    default; reduce it to compress the range, increase it to spread it.

    Parameters
    ----------
    eigvals   : (K,) eigenvalues sorted by |λ| descending (from svd_decompose)
    sigma_max : length scale for the dominant component [nm]
    beta      : power-law exponent (default 0.5)

    Returns
    -------
    sigmas : (K,) per-component length scales [nm], descending
    """
    ratios = np.abs(eigvals) / (np.abs(eigvals[0]) + 1e-12)
    sigmas = sigma_max * np.power(ratios, beta)
    log.debug("compute_sigma_eigenvalue: σ range [%.3f, %.3f] nm",
              sigmas[-1], sigmas[0])
    return sigmas


def compute_sigma_autocorr(
    vecs  : np.ndarray,
    r_bead: float,
    nu    : float = 1.0 / 3.0,
) -> np.ndarray:
    """
    Assign per-component length scales from eigenvector autocorrelation lengths.

    The genomic autocorrelation length of v_k is the weighted mean separation:

        ξ_k = Σ_s  s · |ρ_k(s)|  /  Σ_s |ρ_k(s)|

    where ρ_k(s) = mean( v_k(i) · v_k(i+s) ) over all i.

    Converted to 3D via polymer scaling:

        σ_k = r_bead · ξ_k^ν

    with ν = 1/3 for a fractal/compact globule (default) or 1/2 for an
    ideal chain.

    Parameters
    ----------
    vecs   : (N, K) eigenvectors, columns from svd_decompose()
    r_bead : bead radius / unit length [nm]
    nu     : polymer scaling exponent

    Returns
    -------
    sigmas : (K,) per-component length scales [nm]
    """
    N, K   = vecs.shape
    sigmas = np.zeros(K)
    seps   = np.arange(1, N, dtype=np.float64)   # separations 1 … N-1

    for k in range(K):
        v       = vecs[:, k]
        # vectorised: autocorrelation at each separation s in one pass
        rho     = np.array([np.mean(v[:N - s] * v[s:]) for s in range(1, N)])
        abs_rho = np.abs(rho)
        total   = abs_rho.sum()
        xi_k    = (seps * abs_rho).sum() / (total + 1e-12)
        sigmas[k] = r_bead * (xi_k ** nu)

    log.debug("compute_sigma_autocorr: σ range [%.3f, %.3f] nm",
              sigmas.min(), sigmas.max())
    return sigmas


# ═════════════════════════════════════════════════════════════════════════════
# Layer 3 — OpenMM force builders
# ═════════════════════════════════════════════════════════════════════════════

def build_svd_force(
    A           : np.ndarray,
    r_comp      : float,
    k_scale     : float = 1.0,
    force_group : int   = 1,
) -> mm.CustomNonbondedForce:
    """
    Construct a CustomNonbondedForce from the SVD parameter matrix *A*.

    Physics
    -------
    Given the per-particle matrix A  (shape N × K, from svd_decompose),
    the pairwise energy is:

        U(i, j; r) = −k_scale · [Σ_{k=1}^{K} a_k(i) · a_k(j)] · g(r)

    where g(r) = exp(−r² / (2σ²)),  σ = r_comp / 3,  cutoff at 3σ.

    Same-type beads (A or B compartment) accumulate a positive dot product
    → attractive;  opposite types → repulsive.  This is the rank-K
    generalisation of the existing compartment-block force in MultiMM.

    Parameters
    ----------
    A           : (N_beads, K) ndarray — output of svd_decompose()
    r_comp      : compartment length scale [nm]; sets σ = r_comp / 3
    k_scale     : global energy scale [kJ/mol]
    force_group : OpenMM force-group index

    Returns
    -------
    force : mm.CustomNonbondedForce
    """
    N_beads, K = A.shape
    sigma  = r_comp / 3.0
    cutoff = 3.0 * sigma

    log.info("build_svd_force: N=%d, K=%d, σ=%.3f nm, cutoff=%.3f nm",
             N_beads, K, sigma, cutoff)

    # energy expression: inner product of K-dim vectors, Gaussian envelope
    dot_terms  = " + ".join(f"h{k+1}1*h{k+1}2" for k in range(K))
    expression = f"-hic_k * ({dot_terms}) * exp(-r^2 / (2*hic_sigma^2))"

    force = mm.CustomNonbondedForce(expression)
    force.addGlobalParameter("hic_k",     k_scale)
    force.addGlobalParameter("hic_sigma", sigma)

    for k in range(K):
        force.addPerParticleParameter(f"h{k+1}")

    for i in range(N_beads):
        force.addParticle(list(A[i]))             # K per-particle values

    # NoCutoff: the Gaussian exp(-r²/2σ²) already decays to ~0 beyond the cutoff
    # distance, so a hard cutoff is unnecessary and would conflict with other
    # CustomNonbondedForce instances that use NoCutoff (required by OpenCL/CUDA).
    force.setNonbondedMethod(mm.CustomNonbondedForce.NoCutoff)
    force.setForceGroup(force_group)

    log.info("  SVD force ready  (%d particles, %d parameters/particle)",
             N_beads, K)
    return force


def build_svd_force_multiscale(
    A           : np.ndarray,
    sigmas      : np.ndarray,
    k_scale     : float = 1.0,
    force_group : int   = 1,
) -> mm.CustomNonbondedForce:
    """
    Multi-scale SVD force: each eigenvector gets its own Gaussian σ_k.

    Physics
    -------
        U(i, j; r) = −k_scale · Σ_k  a_k(i)·a_k(j) · exp(−r² / (2σ_k²))

    Each component k acts at its own 3D length scale σ_k:
      k=1 (compartments) → large σ → long-range
      k~5 (TADs)         → medium σ
      k>>5 (loops)       → small σ → short tether

    The cutoff is set to 3·σ_1 (longest-range component); shorter
    components decay naturally within that envelope.

    Parameters
    ----------
    A           : (N_beads, K) ndarray — from svd_decompose()
    sigmas      : (K,) per-component length scales [nm], descending
                  Use compute_sigma_eigenvalue() or compute_sigma_autocorr()
    k_scale     : global energy scale [kJ/mol]
    force_group : OpenMM force-group index

    Returns
    -------
    force : mm.CustomNonbondedForce
    """
    N_beads, K = A.shape
    if len(sigmas) != K:
        raise ValueError(f"sigmas length ({len(sigmas)}) must match K={K}")

    cutoff = 3.0 * float(sigmas[0])          # dominated by longest-range σ

    log.info("build_svd_force_multiscale: N=%d, K=%d, cutoff=%.3f nm",
             N_beads, K, cutoff)
    log.info("  σ: [%s] nm",
             ", ".join(f"{s:.3f}" for s in sigmas))

    # one Gaussian term per component, each with its own global sigma
    terms      = " + ".join(
        f"h{k+1}1*h{k+1}2*exp(-r^2/(2*hic_sigma_{k+1}^2))"
        for k in range(K)
    )
    expression = f"-hic_k * ({terms})"

    force = mm.CustomNonbondedForce(expression)
    force.addGlobalParameter("hic_k", k_scale)
    for k in range(K):
        force.addGlobalParameter(f"hic_sigma_{k+1}", float(sigmas[k]))
    for k in range(K):
        force.addPerParticleParameter(f"h{k+1}")
    for i in range(N_beads):
        force.addParticle(list(A[i]))

    # NoCutoff: Gaussian terms decay naturally; avoids cutoff-method mismatch on OpenCL/CUDA.
    force.setNonbondedMethod(mm.CustomNonbondedForce.NoCutoff)
    force.setForceGroup(force_group)

    log.info("  multi-scale SVD force ready  (%d particles, %d components)",
             N_beads, K)
    return force


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
    # The force formula is F_ij = k_scale·α/r·(c_ij − P(r)).
    # At large r (P≈0): F ≈ k_scale·α·c_ij/r  (attractive).
    # Typical force per bond at a 1 nm inter-bead separation:
    f_typical = k_scale * alpha * float(c_vals.mean()) / 1.0   # kJ/mol/nm
    bonds_per_bead = 2.0 * n_bonds / N_beads
    log.info(
        "  calibration: k_scale=%.2f kJ/mol, α=%.1f, "
        "<c_ij>=%.3f → F/bond@1nm≈%.1f kJ/mol/nm, "
        "%.1f bonds/bead → total~%.0f kJ/mol/nm per bead",
        k_scale, alpha, float(c_vals.mean()),
        f_typical, bonds_per_bead, f_typical * bonds_per_bead,
    )
    # kT ≈ 2.49 kJ/mol at 300 K; warn if total force per bead < ~10 kT/nm
    if f_typical * bonds_per_bead < 25.0:
        log.warning(
            "  Hi-C cross-entropy force may be too weak to drive folding "
            "(total force/bead < 10 kT/nm at 1 nm).  Consider increasing "
            "HIC_K_SCALE (current %.2f kJ/mol); 5–20 kJ/mol is typical.",
            k_scale,
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
# Layer 4 — Public entry point
# ═════════════════════════════════════════════════════════════════════════════

def build_hic_force(
    H_raw            : np.ndarray,
    N_beads          : int,
    r_comp           : float,
    mode             : Literal["svd", "svd_multiscale", "crossentropy"] = "svd",
    K                : int   = 10,
    threshold        : float = 0.01,
    alpha            : float = 3.0,
    k_scale          : float = 1.0,
    force_group      : int   = 1,
    already_balanced : bool  = False,
    sigma_beta       : float = 0.5,
    r_bead           : float = None,
) -> mm.Force:
    """
    Full pipeline: raw Hi-C array → OpenMM Force object.

    Pipeline
    --------
    H_raw
      │  symmetrize_and_clean()
      │  resize_matrix()              ← to N_beads × N_beads
      │  diagonal_normalize()         ← skip if already_balanced=True
      │
      ├─ mode='svd'
      │    compute_oe_matrix()
      │    svd_decompose()            → A  (N × K)
      │    build_svd_force(A)         → CustomNonbondedForce, single σ
      │
      ├─ mode='svd_multiscale'
      │    compute_oe_matrix()
      │    svd_decompose()            → lam, vecs, A
      │    compute_sigma_eigenvalue() → sigmas (K,)   [fast, default]
      │    compute_sigma_autocorr()   → sigmas (K,)   [if r_bead given]
      │    build_svd_force_multiscale → CustomNonbondedForce, per-component σ
      │
      └─ mode='crossentropy'
           build_crossentropy_force(C) → sparse CustomBondForce

    Parameters
    ----------
    H_raw            : (M, M) ndarray — raw Hi-C counts (any resolution).
                       Resampled to N_beads × N_beads if M ≠ N_beads.
    N_beads          : number of simulation beads.
    r_comp           : compartment / contact length scale [nm].
                       All SVD modes: σ_max = r_comp (single) or dominant σ.
                       CrossEntropy: r_c of the sigmoid.
    mode             : 'svd'           → single-σ CustomNonbondedForce
                       'svd_multiscale'→ per-component σ CustomNonbondedForce
                       'crossentropy'  → sparse CustomBondForce
    K                : (SVD modes) rank of eigendecomposition.
    threshold        : (crossentropy) minimum normalised c_ij for a bond.
    alpha            : (crossentropy) sigmoid steepness (2–4 typical).
    k_scale          : global energy scale [kJ/mol].
    force_group      : OpenMM force-group index.
    already_balanced : True → skip diagonal_normalize().
    sigma_beta       : (svd_multiscale) exponent for eigenvalue-based σ scaling.
                       Ignored if r_bead is provided (autocorr method used instead).
    r_bead           : (svd_multiscale, optional) bead radius [nm].
                       When given, uses autocorrelation-based σ (more principled).
                       When None, uses eigenvalue-based σ (faster).

    Returns
    -------
    force : mm.CustomNonbondedForce  (mode='svd' or 'svd_multiscale')
         or mm.CustomBondForce       (mode='crossentropy')

    Examples
    --------
    >>> # single-scale SVD
    >>> force = build_hic_force(hic_array, N_beads=500, r_comp=6.0,
    ...                         mode='svd', K=10, k_scale=1.0)
    >>> system.addForce(force)

    >>> # multi-scale SVD, eigenvalue method
    >>> force = build_hic_force(hic_array, N_beads=500, r_comp=6.0,
    ...                         mode='svd_multiscale', K=15, sigma_beta=0.5)
    >>> system.addForce(force)

    >>> # multi-scale SVD, autocorrelation method
    >>> force = build_hic_force(hic_array, N_beads=500, r_comp=6.0,
    ...                         mode='svd_multiscale', K=15, r_bead=1.0)
    >>> system.addForce(force)

    >>> # cross-entropy
    >>> force = build_hic_force(hic_array, N_beads=500, r_comp=2.0,
    ...                         mode='crossentropy', threshold=0.05)
    >>> system.addForce(force)
    """
    valid_modes = ("svd", "svd_multiscale", "crossentropy")
    if mode not in valid_modes:
        raise ValueError(
            f"mode must be one of {valid_modes}, got '{mode!r}'"
        )

    log_table(
        [
            ("Mode",          mode),
            ("N beads",       N_beads),
            ("r_comp",        f"{r_comp:.3f} nm"),
            ("K (SVD rank)",  K),
            ("k_scale",       f"{k_scale:.3f} kJ/mol"),
            ("Force group",   force_group),
            ("Input shape",   str(H_raw.shape)),
        ],
        title="Hi-C Force — build",
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

    # ── Layers 2+3: decompose + build ───────────────────────────────────────
    if mode == "svd":
        H_oe         = compute_oe_matrix(H)
        lam, vecs, A = svd_decompose(H_oe, K=K)
        force        = build_svd_force(A, r_comp,
                                       k_scale=k_scale,
                                       force_group=force_group)

    elif mode == "svd_multiscale":
        H_oe         = compute_oe_matrix(H)
        lam, vecs, A = svd_decompose(H_oe, K=K)

        if r_bead is not None:
            log.info("  sigma method: autocorrelation (r_bead=%.3f nm)", r_bead)
            sigmas = compute_sigma_autocorr(vecs, r_bead=r_bead)
        else:
            log.info("  sigma method: eigenvalue (beta=%.2f)", sigma_beta)
            sigmas = compute_sigma_eigenvalue(lam, sigma_max=r_comp,
                                              beta=sigma_beta)

        # clip: [one bead diameter, r_comp]
        lo     = r_bead if r_bead is not None else r_comp / 10.0
        sigmas = np.clip(sigmas, lo, r_comp)

        force = build_svd_force_multiscale(A, sigmas,
                                           k_scale=k_scale,
                                           force_group=force_group)

    else:  # mode == "crossentropy"
        force = build_crossentropy_force(H, r_comp,
                                         threshold=threshold,
                                         alpha=alpha,
                                         k_scale=k_scale,
                                         force_group=force_group)

    log.info("build_hic_force: done (%s force, group %d)", mode, force_group)
    log.info("═" * 60)
    return force
