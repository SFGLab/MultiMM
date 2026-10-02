"""
hic_force.py — Hi-C contact-probability force for MultiMM / OpenMM

Three layers: pre-processing (symmetrize_and_clean, diagonal_normalize,
resize_matrix; pure numpy), the OpenMM force builder
(build_crossentropy_force -> sparse CustomBondForce), and the public
entry point (build_hic_force: cleans -> resizes -> normalises -> builds).

Method: models contact *probability* as a monotone-decreasing function of
3-D distance, P_ij(r), and drives it up toward the observed contact
frequency c_ij with a one-sided, attraction-only loss:

    L = -Σ_{i<j} c_ij · log P_ij(r_ij)

F_i = -∇_i L. Since P_ij(r) strictly decreases with r, this gradient is
always attractive (pulls i, j closer), never repulsive — there is no term
that pushes a pair apart for being "too close" relative to its own c_ij.
Force magnitude is directly proportional to c_ij throughout: weak-evidence
pairs (small c_ij) pull weakly, strong ones (c_ij near 1) pull strongly,
and it naturally saturates (vanishes) as r -> 0, where P_ij(r) flattens
out — nothing else. (An earlier version used the full symmetric binary
cross-entropy, c_ij·log P + (1-c_ij)·log(1-P), which is self-regulating in
*both* directions: it also repels a pair whenever the current P_ij(r)
exceeds c_ij, i.e. whenever the structure happens to place that pair
closer than its own weak target implies. That is correct for a generative
P(r) model, but wrong as a *restraint*: any pair with a small-but-nonzero
c_ij that ends up nearby for any other reason (chain connectivity,
confinement, a neighbor's strong pull) gets actively pushed apart by the
Hi-C force itself, compounding with excluded volume into pervasive
unwanted "self-avoidance". The one-sided loss above removes that failure
mode entirely — keeping pairs apart is left to the excluded-volume/
backbone/confinement forces, as it should be.) P_ij(r) is pluggable
(`kernel=`):

  * "gaussian"    P = exp(-r² / (2σ²))                         (default)
  * "power_law"   P = 1 / (1 + (r / r_c)^α)                    (sigmoid; also "sigmoid")
  * "exponential" P = exp(-r / r_c)
  * "erfc"        P = ½ erfc((r - r_c) / (√2 σ_s))             (soft step)
  * "rouse"       P = erfc(r / √(2 s b²))                      (s = separation in beads,
                                                                  b = Kuhn length)

Each pair's loss contribution is also scaled by its own observed contact
strength, w_ij = c_ij, so weak-evidence pairs exert proportionally small
force instead of a hard constraint — directly, with no extra tunable
exponent or inclusion cutoff: every pair the data gives a nonzero c_ij for
(with or without the OE-enrichment normalisation — see `oe_normalize`)
gets a bond, weighted exactly as the data says.

Usage (model.py):
    from hic_force import build_hic_force
    force = build_hic_force(args.hic_matrix, N_beads=self.N,
                             rc=self.hic_rc, kernel="gaussian")
    self.system.addForce(force)
"""

from __future__ import annotations

import logging
import os
from typing import Optional

import numpy as np
from scipy.ndimage import zoom, gaussian_filter, median_filter
from scipy.special import erfc as _erfc, erfcinv as _erfcinv

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
    rc      : float,
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
        globals_ = {"hic_rc": rc, "hic_alpha": alpha}
        needs_sep = False

        def p_func(r, sep=None):
            return 1.0 / (1.0 + (r / rc) ** alpha)

        ref_r = rc

    elif kernel == "gaussian":
        p_expr = "P = exp(-(r^2) / (2 * hic_sigma^2));"
        globals_ = {"hic_sigma": sigma}
        needs_sep = False

        def p_func(r, sep=None):
            return np.exp(-(r ** 2) / (2.0 * sigma ** 2))

        ref_r = sigma

    elif kernel == "exponential":
        p_expr = "P = exp(-r / hic_rc);"
        globals_ = {"hic_rc": rc}
        needs_sep = False

        def p_func(r, sep=None):
            return np.exp(-r / rc)

        ref_r = rc

    elif kernel == "erfc":
        p_expr = "P = 0.5 * erfc((r - hic_rc) / (sqrt(2) * hic_sigmas));"
        globals_ = {"hic_rc": rc, "hic_sigmas": sigma_s}
        needs_sep = False

        def p_func(r, sep=None):
            return 0.5 * _erfc((r - rc) / (np.sqrt(2.0) * sigma_s))

        ref_r = rc

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
    rc      : float,
    alpha       : float           = 3.0,
    sigma       : Optional[float] = None,
    sigma_s     : Optional[float] = None,
    kuhn_length : Optional[float] = None,
    r_min       : Optional[float] = None,
):
    """
    Public helper so other modules (notably `validation.py`) can compute
    the exact same distance→contact-probability kernel P(r) used to build
    the Hi-C force, instead of a generic/unrelated proxy such as 1/(d+ε).

    Resolves None/naive-default parameters exactly like
    `build_crossentropy_force` does — including the dynamic-range
    calibration against `r_min` (see `_calibrate_kernel_dynamics`) — so a
    caller that passes the same (kernel, rc, alpha, r_min) the force was
    built with always gets back the exact same resolved kernel, never a
    mismatched one. Pass `r_min` (the excluded-volume floor distance) here
    whenever the force itself was built with one, so validation scores
    against what was actually optimized.

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

    sigma, alpha, sigma_s, kuhn_length, _ = _calibrate_kernel_dynamics(
        kernel, rc, alpha, sigma, sigma_s, kuhn_length, r_min=r_min,
    )

    _, _, needs_sep, p_func, _ = _kernel_spec(kernel, rc, alpha, sigma, sigma_s, kuhn_length)
    return p_func, needs_sep


# ── Dynamic-range calibration ─────────────────────────────────────────────────
# "rc" is auto-calibrated from the structure's OWN pairwise-distance
# distribution (auto_contact_scale: typically tens of nm), while the smallest
# distance two beads can actually approach is set by excluded volume — a
# couple of orders of magnitude smaller (e.g. ~0.1 nm by default).
#
# A first version of this calibration narrowed EVERY kernel's width to make
# P(r) visibly discriminate between r_min and rc. That was wrong for three of
# the five kernels, and caused a real regression (weak/unreasonable folding,
# simulated contact maps not resembling the experimental one):
#
#   The one-sided loss is U(r) = -k_scale * c_ij * log(P(r) + eps), and what
#   actually matters for the force is d(-log P)/dr, NOT how close P(r) looks
#   to 1 or 0 in absolute terms. For "power_law"/"sigmoid",
#   -log P = log(1 + (r/rc)^alpha) is locally FLAT near r=0 for alpha > 1
#   (dU/dr ~ alpha * r^(alpha-1) -> 0 as r -> 0): a high alpha genuinely
#   starves the force right where pairs are closest, which is the real
#   "saturation kills the force" failure this calibration exists to fix.
#
#   "gaussian", "erfc" and "rouse" do NOT have that failure mode: their
#   -log P is locally QUADRATIC in r (a harmonic spring), so the force is
#   ~ r / width^2 at every distance and never vanishes except exactly at
#   r = 0 — regardless of how saturated the raw P(r) value looks near
#   r_min. Narrowing their width to chase a lower P(rc) target doesn't fix a
#   vanishing force (there isn't one); it just raises 1/width^2, which
#   multiplies the force at EVERY distance, including the large distances
#   most pairs actually start at before folding. Measured concretely for a
#   typical rc=20 nm, r_min=0.1 nm: narrowing gaussian's sigma this way
#   multiplied its force by ~6x at every r (not just near r_min); narrowing
#   rouse's kuhn_length multiplied it by ~8-22x depending on genomic
#   separation. Both are silently equivalent to cranking HIC_K_SCALE several
#   times past its own "very high, risks freezing MD" warning threshold —
#   which is what produced structures that don't resemble the Hi-C target.
#
# So: only "power_law"/"sigmoid" (steepness alpha) and "erfc" (its sigma_s
# calibration solves for a WIDER, gentler step here, not a narrower one —
# safe) are actually recalibrated below. "gaussian" and "rouse" keep their
# original naive defaults (sigma=rc, kuhn_length=rc) unconditionally; P(r)
# at r_min/rc is still measured and logged for them so a poorly-scaled rc
# is still visible, it just isn't auto-"fixed" by stiffening the spring.
_DYNAMIC_RANGE_P_HIGH  = 0.95   # target P at r_min — "clearly in contact"
_DYNAMIC_RANGE_P_LOW   = 0.05   # target P at rc     — "clearly not in contact"
_DYNAMIC_RANGE_SEP_REF = 10.0   # reference genomic separation (beads), used
                                 # only to calibrate the rouse kernel's
                                 # kuhn_length — NOT used for per-pair P(r)
                                 # evaluation, which always uses the true
                                 # per-pair separation at runtime. Fixed (not
                                 # derived from N_beads/bond set) so model.py
                                 # and validation.py — which call into this
                                 # kernel machinery from different contexts —
                                 # always agree on the same calibrated value.


def _calibrate_kernel_dynamics(
    kernel      : str,
    rc      : float,
    alpha       : float,
    sigma       : Optional[float],
    sigma_s     : Optional[float],
    kuhn_length : Optional[float],
    r_min       : Optional[float] = None,
    sep_ref     : float           = _DYNAMIC_RANGE_SEP_REF,
    p_high      : float           = _DYNAMIC_RANGE_P_HIGH,
    p_low       : float           = _DYNAMIC_RANGE_P_LOW,
):
    """Resolve each kernel's None parameter, recalibrating steepness/width
    ONLY for the kernels where doing so actually fixes a vanishing force
    near r_min without also inflating it at large r (see module note
    above: "power_law"/"sigmoid" alpha, and "erfc" sigma_s). "gaussian" and
    "rouse" always get their plain naive default (sigma=rc / kuhn_length=rc)
    — recalibrating those would only stiffen an already-nonvanishing spring
    force at every distance, which is what caused the regression this
    revision fixes.

    Only touches a parameter that was left at its naive default
    (alpha exactly at the library default of 3.0, or sigma_s=None for
    "erfc") — anything the user set explicitly is resolved exactly as
    before and only measured, never overridden. Deterministic given
    (kernel, rc, alpha, r_min, sep_ref), so model.py (building the force)
    and validation.py (via `get_kernel_p_func`, scoring against it) always
    agree on the same resolved values without needing to pass them to each
    other explicitly.

    Returns (sigma, alpha, sigma_s, kuhn_length, diagnostics). diagnostics
    is None when r_min is not given (nothing to measure against); otherwise
    a dict with p_at_rmin/p_at_rc for logging, computed for every kernel
    regardless of whether it was actually recalibrated.
    """
    kernel = kernel.lower()
    naive_sigma_s = sigma_s is None
    naive_alpha   = (alpha == 3.0)   # library default — see docstring

    # Naive fallbacks first — unchanged behaviour whenever calibration can't
    # run (r_min unavailable) or doesn't apply to this kernel/parameter.
    # "gaussian"/"rouse" ALWAYS stay here (see module note): recalibrating
    # them stiffens the force everywhere, not just near r_min.
    sigma_r       = rc if sigma is None else sigma
    sigma_s_r     = 0.3 * rc if sigma_s is None else sigma_s
    kuhn_length_r = rc if kuhn_length is None else kuhn_length
    alpha_r       = alpha

    if r_min is None or r_min <= 0 or rc <= 0:
        return sigma_r, alpha_r, sigma_s_r, kuhn_length_r, None

    eps = 1e-9
    p_high = min(max(p_high, eps), 1.0 - eps)
    p_low  = min(max(p_low,  eps), 1.0 - eps)

    try:
        if kernel in ("power_law", "sigmoid") and naive_alpha:
            # P(rc) = 0.5 always (sigmoid midpoint, fixed by construction);
            # calibrate steepness instead so P(r_min) = p_high:
            #   alpha = ln((1-p_high)/p_high) / ln(r_min/rc)
            # Safe: -log(P) grows only ~ alpha*log(r/rc) at large r, so a
            # higher alpha never makes the far-field force explode.
            ratio = r_min / rc
            if 0.0 < ratio < 1.0:
                alpha_r = float(np.clip(
                    np.log((1.0 - p_high) / p_high) / np.log(ratio), 0.5, 50.0,
                ))

        elif kernel == "erfc" and naive_sigma_s:
            # P(r_min) = p_high => sigma_s = (r_min - rc) / (sqrt(2)*erfcinv(2*p_high))
            # In the typical regime (r_min << rc) this WIDENS sigma_s
            # relative to the naive 0.3*rc (which already overshoots
            # p_high, e.g. P(r_min)~0.9997 at the naive default) — a
            # gentler, not stiffer, step — so it's safe to apply.
            x = _erfcinv(2.0 * p_high)
            if abs(x) > eps:
                sigma_s_r = abs((r_min - rc) / (np.sqrt(2.0) * x))

        # "gaussian": -log(P) = r^2/(2 sigma^2) is a harmonic spring, force
        # ~ r/sigma^2 at every r, never vanishing near r_min regardless of
        # how saturated P(r_min) looks — nothing to fix, and narrowing
        # sigma would only multiply the force everywhere (measured ~6x for
        # a typical rc/r_min ratio). Left at its naive default always.
        #
        # "rouse": same harmonic-spring argument applies to erfc's far
        # tail, and shrinking kuhn_length this way was measured to inflate
        # the force ~8-22x depending on genomic separation. Left at its
        # naive default always.
        #
        # "exponential" has one parameter (rc itself doubles as the decay
        # constant) — no separate steepness to tune without redefining what
        # rc means for this kernel, so it's left alone; the diagnostics
        # below will still show whether its achievable range is poor.
    except (ValueError, FloatingPointError, ZeroDivisionError, OverflowError):
        alpha_r = alpha
        sigma_s_r = 0.3 * rc if sigma_s is None else sigma_s

    _, _, _, p_func, _ = _kernel_spec(kernel, rc, alpha_r, sigma_r, sigma_s_r, kuhn_length_r)
    diagnostics = {
        "p_at_rmin": float(p_func(r_min, sep_ref)),
        "p_at_rc":   float(p_func(rc, sep_ref)),
    }
    return sigma_r, alpha_r, sigma_s_r, kuhn_length_r, diagnostics


def auto_contact_scale(
    coords: np.ndarray,
    percentile: float = 10.0,
    n_samples: int = 20000,
    seed: int = 0,
) -> float:
    """Estimate a contact length scale from a structure's own pairwise-distance
    distribution, instead of reusing the Hi-C force's microscopic
    ``rc``/``sigma`` (tuned for bead-overlap, so it underflows to ~0 for
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
    c_ij = 0 exactly, which makes their force weight w_ij = c_ij exactly 0
    too: the cross-entropy term vanishes regardless of current
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


_DENOISE_MEDIAN_SIZE     = 3     # median-filter window, beads — removes pixels
                                 # with no support from their neighbours
_DENOISE_SIGMA           = 1.5   # Gaussian std, beads — stricter than validation's
                                 # display-only smoothing (sigma=1.0), since this
                                 # matrix is the force's actual optimization target.
_DENOISE_CLIP_PERCENTILE = 99.0  # safety-net cap for anything the median filter
                                 # doesn't catch, before the final smoothing pass


def denoise_contact_matrix(
    H: np.ndarray,
    sigma: float = _DENOISE_SIGMA,
    clip_percentile: float = _DENOISE_CLIP_PERCENTILE,
    median_size: int = _DENOISE_MEDIAN_SIZE,
) -> np.ndarray:
    """Stricter automatic denoising of a contact matrix — used both for the
    final c_ij matrix that actually drives the Hi-C force (applied *after*
    OE normalisation, if enabled, not on the raw contact counts — see
    `preprocess_hic_matrix`) and, for consistency, by validation.py when
    comparing simulated/RW contact maps against the experimental one.

    Three passes, in order:

    1. Median filter (``median_size``-bead window). This is the step that
       actually tells noise apart from real signal: a genuine compartment
       or TAD interaction shows up as a patch spanning *several* beads in
       both directions, so it survives a small local median unchanged. An
       isolated noisy pixel with no support from its neighbours — a single
       bead pair with an outlier count and nothing similar around it — does
       not survive: the median of its neighbourhood is close to background,
       so the spike is replaced by that background value instead of merely
       being blurred (which is all a plain Gaussian would do — it spreads
       an isolated spike's mass into its neighbours rather than removing
       it, and can even manufacture a small fake "patch" out of a single
       noisy pixel). This is the key fix for sparse, single-pixel noise
       that would otherwise still reach the Hi-C force as if it were a real
       interaction.
    2. A percentile safety-net cap (catches anything that still stands out,
       e.g. right at the matrix edge where the median window is truncated).
    3. A light Gaussian pass for the final smoothing/visual polish, then
       symmetrize and re-zero the diagonal.

    No new user-facing parameters: this always runs as part of
    ``preprocess_hic_matrix``.
    """
    H = median_filter(H, size=median_size, mode="nearest")
    nz = H[H > 0]
    if nz.size:
        cap = float(np.percentile(nz, clip_percentile))
        if cap > 0:
            H = np.minimum(H, cap)
    H_smoothed = gaussian_filter(H, sigma=sigma, mode="nearest")
    H_out = 0.5 * (H_smoothed + H_smoothed.T)
    np.fill_diagonal(H_out, 0.0)
    np.clip(H_out, 0.0, None, out=H_out)
    return H_out


def preprocess_hic_matrix(
    H_raw: np.ndarray,
    N_beads: int,
    already_balanced: bool = False,
    oe_normalize: bool = False,
    save_path: Optional[str] = None,
    chrom: Optional[str] = None,
) -> np.ndarray:
    """Run the same clean -> resize -> balance -> (optional) OE-floor ->
    denoise pipeline `build_hic_force` uses, without building the OpenMM
    force — so diagnostics can see exactly the c_ij values the force
    actually targets (see `_oe_normalize_matrix` for what `oe_normalize`
    does, and `denoise_contact_matrix` for the denoising step).

    save_path/chrom: if save_path is given, a before/after denoising
    comparison figure is saved to ``<save_path>/plots/hic_preprocessing.png``
    (both panels shown OE-normalized when ``oe_normalize=True``, since that
    is the actual quantity being denoised in that case).
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

    # ── denoise — applied AFTER OE normalisation (if enabled), since that is
    # the matrix that actually becomes the force's c_ij target, not the raw
    # counts. OE normalisation divides by a (possibly tiny) per-diagonal
    # mean, which can amplify rather than suppress shot noise, so smoothing
    # has to happen on the post-OE matrix to be effective. ────────────────
    H_before_denoise = H.copy()
    H = denoise_contact_matrix(H)
    n_suppressed = int(np.count_nonzero((H_before_denoise > 0) & (H <= 0)))
    log.info(
        "  Denoised c_ij matrix (%s OE normalisation): median filter (%dx%d bead) "
        "removed %d isolated/unsupported pixel(s) entirely, then outlier-capped at "
        "%.0fth pct + Gaussian sigma=%.1f bead — suppresses per-pixel noise before "
        "it drives the Hi-C force, while contiguous compartment/TAD-scale patches "
        "survive unchanged.",
        "after" if oe_normalize else "without",
        _DENOISE_MEDIAN_SIZE, _DENOISE_MEDIAN_SIZE, n_suppressed,
        _DENOISE_CLIP_PERCENTILE, _DENOISE_SIGMA,
    )

    if save_path is not None:
        try:
            from .plots import plot_hic_preprocessing
            plots_dir = os.path.join(save_path, "plots")
            plot_hic_preprocessing(
                H_before_denoise, H, plots_dir, chrom=chrom,
                sigma=_DENOISE_SIGMA, oe_normalized=oe_normalize,
            )
        except Exception as exc:  # pragma: no cover - plotting must never break the force build
            log.warning("Could not save Hi-C before/after denoising plot: %s", exc)

    return H


def build_crossentropy_force(
    C            : np.ndarray,
    rc       : float,
    kernel       : str             = "gaussian",
    alpha        : float           = 3.0,
    sigma        : Optional[float] = None,
    sigma_s      : Optional[float] = None,
    kuhn_length  : Optional[float] = None,
    k_scale      : float           = 1.0,
    force_group  : int             = 1,
    return_controller: bool        = False,
    noise_seed   : int             = 0,
    r_min        : Optional[float] = None,
):
    """
    Construct a sparse CustomBondForce from a one-sided, attraction-only
    loss between observed contact frequencies and a model contact
    probability P(r), with a pluggable distance→probability kernel.

    Physics
    -------
    Per-pair energy, weighted by the pair's own observed contact strength:

        U_ij(r) = -k_scale · c_ij · log(P + ε)

    where P = P(r) comes from the chosen `kernel` (see module docstring
    for the five options) and the force is F_i = -∇_i U summed over all
    bonds, which OpenMM derives automatically from the symbolic
    expression. Since P(r) strictly decreases with r, this gradient is
    always attractive (pulls i, j closer) and never repulsive, for any
    c_ij > 0 and any current distance — there is no term that pushes a
    pair apart. (This deliberately drops the `(1-c_ij)·log(1-P+ε)` term a
    full binary cross-entropy loss would have: that term actively repels
    a pair whenever its current P(r) exceeds c_ij, which — once every
    nonzero-c_ij pair gets a bond, not just the strongly-supported ones —
    meant huge numbers of weak-but-included pairs would get pushed apart
    just for happening to end up close together, compounding with
    excluded volume into unwanted pervasive "self-avoidance". Keeping
    pairs apart is excluded volume's job, not this force's.)

    The weight is exactly c_ij — no separate tunable exponent — so force
    scales continuously and directly with whatever the data says (the raw
    normalised contact frequency, or its OE-enrichment version if
    `oe_normalize` was applied upstream in `build_hic_force`): small c_ij
    means a weak pull, large c_ij means a strong one, nothing else. Every
    pair with a nonzero C[i,j] gets a bond (pairs that are exactly 0 — e.g.
    floored-out by OE normalisation — contribute literally zero force
    either way, via c_ij itself, so they're skipped purely to avoid
    building inert bonds, not as a tunable cutoff).

    Parameters
    ----------
    C            : (N_beads, N_beads) ndarray — normalised contact matrix
                   with values in [0, 1].  Typically output of
                   diagonal_normalize() divided by its maximum.
    rc       : characteristic contact distance [nm], used as r_c for
                   "power_law"/"sigmoid"/"exponential"/"erfc".  Ignored
                   by "gaussian" (use `sigma`) and "rouse" (use
                   `kuhn_length`), except as their fallback default.
    kernel       : one of VALID_KERNELS. Default "gaussian".
    alpha        : "power_law" sigmoid steepness (typical: 2–4).
    sigma        : "gaussian" width [nm]. Defaults to rc.
    sigma_s      : "erfc" softening width [nm]. Defaults to 0.3 * rc.
    kuhn_length  : "rouse" Kuhn (statistical segment) length [nm].
                   Defaults to rc.
    k_scale      : global energy scale [kJ/mol].
    force_group  : OpenMM force-group index.
    return_controller : if True, also return a `HiCNoiseController` that
                   can periodically redraw each bond's c_ij with fresh
                   noise around its original value (see that class).
    noise_seed   : RNG seed for the returned controller.
    r_min        : excluded-volume floor distance [nm] — the closest two
                   beads can realistically get (e.g. the bead's own EV
                   radius). When given, alpha ("power_law"/"sigmoid") or
                   sigma_s ("erfc"), if left at its default, is recalibrated
                   so P(r) actually discriminates near r_min instead of
                   starving the force there. sigma ("gaussian") and
                   kuhn_length ("rouse") are NOT recalibrated this way even
                   when r_min is given — for those two kernels the force
                   never vanishes near r_min regardless of P(r)'s absolute
                   value, so narrowing them would only make the force
                   stronger at every distance, not fix anything (see
                   `_calibrate_kernel_dynamics` module note). None (default)
                   keeps the old naive defaults unchanged for every kernel.

    Returns
    -------
    force : mm.CustomBondForce, or (force, HiCNoiseController) if
            return_controller=True.
    """
    kernel = kernel.lower()
    if kernel not in VALID_KERNELS:
        raise ValueError(f"Unknown Hi-C kernel '{kernel}'. Valid options: {VALID_KERNELS}")

    sigma, alpha, sigma_s, kuhn_length, dyn_diag = _calibrate_kernel_dynamics(
        kernel, rc, alpha, sigma, sigma_s, kuhn_length, r_min=r_min,
    )

    N_beads = C.shape[0]
    log.info(
        "build_crossentropy_force: N=%d, kernel=%s, rc=%.3f nm",
        N_beads, kernel, rc,
    )
    if dyn_diag is not None:
        _dyn_note = {
            "exponential": " — exponential has no separate steepness parameter to "
                           "calibrate; consider a different kernel if this range looks too narrow",
            "gaussian": " — measured only, not auto-tightened: gaussian's force "
                        "never vanishes near r_min regardless of this P value, so "
                        "narrowing sigma would only stiffen the force at every "
                        "distance (see module note on _calibrate_kernel_dynamics); "
                        "adjust HIC_GAUSSIAN_SIGMA by hand if this range looks wrong",
            "rouse": " — measured only, not auto-tightened: rouse's force never "
                     "vanishes near r_min regardless of this P value, so shrinking "
                     "kuhn_length would only stiffen the force at every distance "
                     "(see module note on _calibrate_kernel_dynamics); adjust "
                     "HIC_ROUSE_KUHN_LENGTH by hand if this range looks wrong",
        }.get(kernel, "")
        log.info(
            "  dynamic-range calibration (r_min=%.4f nm): P(r_min)=%.3f, P(rc)=%.3f "
            "(targets: ~%.2f / ~%.2f)%s",
            r_min, dyn_diag["p_at_rmin"], dyn_diag["p_at_rc"],
            _DYNAMIC_RANGE_P_HIGH, _DYNAMIC_RANGE_P_LOW, _dyn_note,
        )
        if dyn_diag["p_at_rmin"] < 0.5 and kernel not in ("gaussian", "rouse"):
            log.warning(
                "  Even the closest achievable distance (r_min=%.4f nm) only reaches "
                "P=%.3f for kernel=%s — rc (%.3f nm) may be too small relative to "
                "r_min, so the Hi-C force will barely act on most pairs.",
                r_min, dyn_diag["p_at_rmin"], kernel, rc,
            )

    # normalise to [0, 1]
    C_max = C.max()
    if C_max > 0:
        C = C / C_max
    else:
        log.warning("  contact matrix is all-zero — no bonds will be added")

    # upper-triangle pairs (exclude diagonal by construction). Pairs that are
    # exactly 0 are skipped — not a tunable cutoff, just avoiding building a
    # bond that would contribute literally zero force (weight = c_ij = 0)
    # either way; every pair the data gives a nonzero value for gets a bond.
    rows, cols = np.triu_indices(N_beads, k=1)
    mask       = C[rows, cols] > 0
    rows, cols = rows[mask], cols[mask]
    c_vals     = C[rows, cols]
    sep_vals   = (cols - rows).astype(np.float64)   # genomic separation, beads

    n_bonds = len(rows)
    log.info("  pairs with nonzero contact strength: %d  (%.2f%% of upper triangle)",
             n_bonds, 100.0 * n_bonds / (N_beads * (N_beads - 1) / 2))

    if n_bonds == 0:
        log.warning("  no bonds added — the contact matrix is all-zero "
                    "after normalisation; check matrix normalisation / "
                    "OE enrichment settings")

    p_expr, kernel_globals, needs_sep, p_func, ref_r = _kernel_spec(
        kernel, rc, alpha, sigma, sigma_s, kuhn_length,
    )

    # ── calibration diagnostics (approximate — numeric, kernel-agnostic) ─────
    # One-sided loss U = -k*c*log(P): |dU/dr| = k*c*|dP/dr|/P (no P*(1-P)
    # term — that belonged to the dropped two-sided cross-entropy form).
    w_mean = float(c_vals.mean()) if n_bonds else 0.0   # weight = c_ij directly
    dr     = max(ref_r * 1e-3, 1e-6)
    sep_ref = float(np.median(sep_vals)) if n_bonds else 1.0
    P_lo   = p_func(ref_r - dr, sep_ref)
    P_hi   = p_func(ref_r + dr, sep_ref)
    P_ref  = p_func(ref_r, sep_ref)
    dP_dr  = (P_hi - P_lo) / (2.0 * dr)
    denom  = max(P_ref, 1e-6)
    f_typical = k_scale * w_mean * abs(dP_dr) / denom      # kJ/mol/nm (approx)
    bonds_per_bead = 2.0 * n_bonds / N_beads
    log.info(
        "  calibration (kernel=%s): k_scale=%.2f kJ/mol, <c_ij>=%.3f → "
        "F/bond@r≈ref≈%.1f kJ/mol/nm, "
        "%.1f bonds/bead → total~%.0f kJ/mol/nm per bead",
        kernel, k_scale, w_mean,
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
            "5–20 kJ/mol is typical; the per-pair weight (= c_ij, a weak "
            "contact already pulls little) does not independently soften "
            "this.",
            k_scale,
        )

    # OpenMM expression — one-sided, attraction-only: no (1-c_ij)*log(1-P)
    # term, so there is no branch that can ever push a pair apart. The
    # per-pair weight is exactly c_ij (no separate exponent), substituted
    # directly rather than declared as its own auxiliary term.
    expression = (
        "-hic_k * c_ij * log(P + hic_eps);"
        f"{p_expr}"
    )

    force = mm.CustomBondForce(expression)
    force.addGlobalParameter("hic_k",    k_scale)
    force.addGlobalParameter("hic_eps",  1e-8)
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
    rc           : float,
    kernel           : str             = "gaussian",
    alpha            : float           = 3.0,
    sigma            : Optional[float] = None,
    sigma_s          : Optional[float] = None,
    kuhn_length      : Optional[float] = None,
    k_scale          : float           = 1.0,
    force_group      : int             = 1,
    already_balanced : bool            = False,
    oe_normalize     : bool            = False,
    return_controller: bool            = False,
    noise_seed       : int             = 0,
    save_path        : Optional[str]   = None,
    chrom            : Optional[str]   = None,
    r_min            : Optional[float] = None,
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
      │  denoise                      ← always; see preprocess_hic_matrix
      └─ build_crossentropy_force()   → sparse CustomBondForce

    Parameters
    ----------
    H_raw            : (M, M) ndarray — raw Hi-C counts (any resolution).
                       Resampled to N_beads × N_beads if M ≠ N_beads.
    N_beads          : number of simulation beads.
    rc           : characteristic contact distance [nm]. Used
                       directly by "power_law"/"sigmoid"/"exponential"/
                       "erfc" as r_c, and as the fallback default for
                       `sigma` / `kuhn_length` when those are not given.
    kernel           : distance→probability model. One of VALID_KERNELS.
                       Default "gaussian".
    alpha            : "power_law" sigmoid steepness (2–4 typical).
    sigma            : "gaussian" width [nm] (defaults to rc).
    sigma_s          : "erfc" softening width [nm] (defaults to 0.3*rc).
    kuhn_length      : "rouse" Kuhn length [nm] (defaults to rc).
    k_scale          : global energy scale [kJ/mol].
                       Recommended range: 5–20 kJ/mol.
                       Values > 30 kJ/mol risk freezing MD sampling.
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
    save_path        : if given, a before/after denoising comparison figure
                       is saved to ``<save_path>/plots/hic_preprocessing.png``
                       (both panels OE-normalized when oe_normalize=True).
    chrom            : optional chromosome label for the plot title.
    r_min            : excluded-volume floor distance [nm]. When given,
                       recalibrates alpha ("power_law"/"sigmoid") or
                       sigma_s ("erfc"), if left at their default, so P(r)
                       actually discriminates near r_min instead of
                       starving the force there. Does NOT touch sigma
                       ("gaussian") or kuhn_length ("rouse") — for those two
                       the force never vanishes near r_min in the first
                       place, so narrowing them would only stiffen the
                       force at every distance — see
                       `build_crossentropy_force`/`_calibrate_kernel_dynamics`.

    Returns
    -------
    force : mm.CustomBondForce, or (force, HiCNoiseController) if
            return_controller=True.

    Examples
    --------
    >>> force = build_hic_force(hic_array, N_beads=500, rc=0.15,
    ...                         kernel="gaussian", k_scale=10.0)
    >>> system.addForce(force)
    """
    log_table(
        [
            ("N beads",       N_beads),
            ("Kernel",        kernel),
            ("rc",        f"{rc:.3f} nm"),
            ("alpha",         alpha if kernel in ("power_law", "sigmoid") else "n/a"),
            ("k_scale",       f"{k_scale:.3f} kJ/mol"),
            ("Force group",   force_group),
            ("OE normalize",  oe_normalize),
            ("Input shape",   str(H_raw.shape)),
            ("r_min (EV floor)", f"{r_min:.4f} nm" if r_min is not None else "not given — naive kernel defaults"),
        ],
        title="Hi-C Force — build (contact-probability cross-entropy)",
        log_fn=log.info,
    )

    # ── Layer 1: pre-processing ──────────────────────────────────────────────
    H = preprocess_hic_matrix(
        H_raw, N_beads, already_balanced=already_balanced, oe_normalize=oe_normalize,
        save_path=save_path, chrom=chrom,
    )

    # ── Layer 2: build force ─────────────────────────────────────────────────
    force = build_crossentropy_force(
        H, rc,
        kernel=kernel,
        alpha=alpha,
        sigma=sigma,
        sigma_s=sigma_s,
        kuhn_length=kuhn_length,
        k_scale=k_scale,
        force_group=force_group,
        return_controller=return_controller,
        noise_seed=noise_seed,
        r_min=r_min,
    )

    log.info("build_hic_force: done (kernel=%s, group %d)", kernel, force_group)
    log.info("═" * 60)
    return force
