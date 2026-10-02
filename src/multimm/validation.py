import logging

from .logger import log_table, ProgressLogger

import numpy as np
import scipy.ndimage as ndimage
from scipy.stats import pearsonr, spearmanr, mannwhitneyu

from .utils import get_coordinates_cif
from .read_hic import pool_to_n_beads as _pool_to_n_beads
from .hic_force import get_kernel_p_func, auto_contact_scale, preprocess_hic_matrix, denoise_contact_matrix

logger = logging.getLogger(__name__)


# Progress bar: see logger.ProgressLogger (shared across the project, not
# re-defined per module).


# =============================================================================
# Hi-C model validation — experiment vs. model comparisons
# =============================================================================

# ── Low-level metric helpers ──────────────────────────────────────────────────

def _pool_matrix(m: np.ndarray, N: int) -> np.ndarray:
    """Re-sample *m* to an N×N matrix using fractional-overlap pooling.

    Delegates to :func:`read_hic.pool_to_n_beads` which uses a fractional
    overlap weight matrix (``np.linspace`` bin edges) so that every genomic bin
    contributes proportionally to the output — no bins are silently truncated
    when ``m.shape[0] % N != 0``.
    """
    if m.shape[0] == N:
        return m.astype(float)
    return _pool_to_n_beads(m, N)


def _smooth_matrix(m: np.ndarray, sigma: float = 1.0) -> np.ndarray:
    """Denoise a contact matrix, applied identically to simulated/
    experimental/RW maps before any metric is computed — so every metric
    below compares against the *denoised* heatmap, the same one the Hi-C
    force actually optimizes against (see `hic_force.preprocess_hic_matrix`
    / `hic_force.denoise_contact_matrix`), never the raw/noisy one.
    Suppresses the heavy-tailed, single-close-approach noise in the
    1/(d+ε)-style contact proxies (and per-pixel shot noise in the
    experimental matrix) without blurring out genuine TAD/compartment-scale
    structure; outlier pixels are capped before smoothing so a single noise
    spike isn't blurred outward into its neighbours.

    m: (N, N) contact matrix. sigma: Gaussian std (beads); 0/None disables.
    """
    if not sigma:
        return m
    return denoise_contact_matrix(m, sigma=sigma)


def diagonal_decay_profile(mat: np.ndarray, max_diag: int | None = None) -> np.ndarray:
    """Per-diagonal mean contact profile (diagonal decay).

    Parameters
    ----------
    mat : ndarray, shape (N, N)
    max_diag : int, optional
        Only compute diagonals 1..max_diag.  Defaults to N-1.

    Returns
    -------
    profile : ndarray, shape (max_diag,)
        profile[k] = mean of the k-th super-diagonal.
    """
    N = mat.shape[0]
    if max_diag is None:
        max_diag = N - 1
    profile = np.array([np.mean(np.diag(mat, k)) for k in range(1, max_diag + 1)])
    return profile


def insulation_score(mat: np.ndarray, window: int = 10) -> np.ndarray:
    """Sliding-window insulation score.

    score[i] = mean of the *window* × *window* sub-matrix centred on the
    diagonal at position i.  Low scores mark TAD boundaries.

    Parameters
    ----------
    mat : ndarray, shape (N, N)
    window : int
        Half-width of the sliding square (full width = 2*window).

    Returns
    -------
    score : ndarray, shape (N,)
    """
    N = mat.shape[0]
    score = np.full(N, np.nan)
    for i in range(window, N - window):
        block = mat[i - window : i, i : i + window]
        score[i] = np.mean(block)
    return score


def oe_matrix(mat: np.ndarray) -> np.ndarray:
    """Observed / Expected contact matrix (divide each diagonal by its mean)."""
    N = mat.shape[0]
    oe = np.zeros_like(mat, dtype=float)
    for k in range(N):
        diag = np.diag(mat, k)
        mean_k = diag.mean()
        if mean_k > 0:
            oe_diag = diag / mean_k
        else:
            oe_diag = diag
        idx = np.arange(N - k)
        oe[idx, idx + k] = oe_diag
        oe[idx + k, idx] = oe_diag
    return oe


def pc1_of_oe(mat: np.ndarray) -> np.ndarray:
    """First principal component (PC1) of the O/E-normalised contact matrix.

    The sign convention follows Hi-C practice: PC1 is flipped so that the
    sign correlates with gene density (positive ↔ A compartment), but since
    we compare two PC1 vectors the relative sign is still arbitrary; callers
    should correlate |pc1_sim| with |pc1_exp|, or use the absolute Pearson r.

    Parameters
    ----------
    mat : ndarray, shape (N, N)

    Returns
    -------
    pc1 : ndarray, shape (N,)
    """
    from sklearn.decomposition import PCA
    oe = oe_matrix(mat)
    # mean-centre rows
    oe_c = oe - oe.mean(axis=1, keepdims=True)
    pca = PCA(n_components=1)
    pc1 = pca.fit_transform(oe_c)[:, 0]
    return pc1


def inverse_contact_matrix(dist_map: np.ndarray, eps: float = 1e-3) -> np.ndarray:
    """Convert a pairwise distance matrix to a contact proxy via 1/(d + ε).

    Kept for backward compatibility / standalone use (e.g. notebooks). The
    validation pipeline itself no longer uses this generic proxy — see
    :func:`_coords_to_kernel_contact_f32`, which instead evaluates the exact
    same P(r) kernel used to build the Hi-C cross-entropy force, so that
    validation is numerically consistent with what the force actually
    optimises against.

    Parameters
    ----------
    dist_map : ndarray, shape (N, N)
    eps : float
        Regularisation to avoid division by zero.

    Returns
    -------
    contact : ndarray, shape (N, N)
    """
    return 1.0 / (dist_map + eps)


# ── Correlation helpers ───────────────────────────────────────────────────────

def _pearson(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    """Pearson r between two 1-D arrays; returns (r, p).  NaN-safe."""
    mask = np.isfinite(a) & np.isfinite(b)
    if mask.sum() < 3:
        return float("nan"), float("nan")
    r, p = pearsonr(a[mask], b[mask])
    return float(r), float(p)



def _spearman(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    """Spearman r between two 1-D arrays; returns (r, p).  NaN-safe."""
    mask = np.isfinite(a) & np.isfinite(b)
    if mask.sum() < 3:
        return float("nan"), float("nan")
    r, p = spearmanr(a[mask], b[mask])
    return float(r), float(p)


def _upper_tri(mat: np.ndarray) -> np.ndarray:
    """Return flattened upper triangle (k=1) of a square matrix."""
    idx = np.triu_indices_from(mat, k=1)
    return mat[idx]


def _check_contact_distance_monotonicity(
    coords: np.ndarray,
    contact_matrix: np.ndarray,
    log,
    n_samples: int = 8000,
    seed: int = 0,
    label: str = "structure",
) -> float:
    """Empirically verify that higher contact-proxy strength corresponds to
    a shorter 3-D distance for this structure/contact matrix — a sanity
    check on the pipeline, not the kernel (which is monotonic by
    construction). Logs the Spearman correlation between sampled contact
    values and their 3-D distances; strongly negative confirms the rule."""
    N = len(coords)
    if N < 3:
        return float("nan")
    rng = np.random.default_rng(seed)
    n_samples = min(n_samples, N * (N - 1) // 2)
    i_idx = rng.integers(0, N, size=n_samples * 2)
    j_idx = rng.integers(0, N, size=n_samples * 2)
    valid = i_idx != j_idx
    i_idx, j_idx = i_idx[valid][:n_samples], j_idx[valid][:n_samples]

    dists = np.linalg.norm(coords[i_idx] - coords[j_idx], axis=1)
    contacts = np.asarray(contact_matrix)[i_idx, j_idx]

    rho, p = _spearman(contacts, dists)
    ok = np.isfinite(rho) and rho < 0
    log.info(
        "  Sanity check (%s) — contact↔distance monotonicity: "
        "Spearman ρ=%.3f (p=%.1e)  %s",
        label, rho, p,
        "✓ higher contact ⇒ closer distance, as expected" if ok
        else "⚠ unexpected sign — inspect kernel/contact pipeline",
    )
    return rho


# ── Random-walk baseline ─────────────────────────────────────────────────────

def _rw_baseline_contact(
    N: int,
    p_func,
    needs_sep: bool = False,
    sep_matrix: np.ndarray | None = None,
    n_rw: int = 20,
    step_nm: float = 0.1,
    seed: int = 0,
    confine_radius_nm: float | None = None,
    log=None,
) -> np.ndarray:
    """Generate a null-model (random-walk ensemble) average contact map.

    Builds *n_rw* 3-D random walks of length *N* (bond length *step_nm*),
    optionally confined inside a sphere of radius *confine_radius_nm* via
    elastic reflection (recommended: set to the simulation's nuclear radius
    — unconfined walks are unrealistically large and inflate the RW
    baseline), computes each one's contact-probability proxy via the same
    P(r) kernel as the Hi-C force (``p_func``), and returns the float64
    N×N average.

    p_func: callable(r)->P, or callable(r, sep)->P if needs_sep. sep_matrix:
    (N, N) genomic-separation matrix, required when needs_sep. n_rw: number
    of realisations (default 20). step_nm: bond length (match
    POL_HARMONIC_BOND_R0). seed: RNG seed.
    """
    _log = log or logger
    _log.info("Generating random-walk null baseline over %d realisations …", n_rw)
    rng = np.random.default_rng(seed)
    acc: np.ndarray | None = None
    _progress = ProgressLogger(n_rw, _log, "RW realisations")
    for _rw_i in range(n_rw):
        # Independent 3-D Gaussian random walk with fixed bond length
        steps = rng.normal(0.0, 1.0, size=(N - 1, 3))
        steps *= step_nm / np.linalg.norm(steps, axis=1, keepdims=True)

        if confine_radius_nm is not None:
            R = float(confine_radius_nm)
            # Walk bead-by-bead with sphere-wall reflection
            pos = np.zeros(3)
            coords_list = [pos.copy()]
            for s in steps:
                candidate = pos + s
                dist = np.linalg.norm(candidate)
                if dist > R:
                    # Reflect: mirror the overshoot back along radial direction
                    radial_unit = candidate / (dist + 1e-30)
                    overshoot   = dist - R
                    candidate   = candidate - 2.0 * overshoot * radial_unit
                    # Safety: clamp (rare numerical edge case)
                    d2 = np.linalg.norm(candidate)
                    if d2 > R:
                        candidate = candidate * (R / (d2 + 1e-30))
                pos = candidate
                coords_list.append(pos.copy())
            coords = np.array(coords_list, dtype=np.float32)
        else:
            coords = np.vstack(
                [np.zeros((1, 3)), np.cumsum(steps, axis=0)]
            )

        frame = _coords_to_kernel_contact_f32(coords, p_func, sep_matrix if needs_sep else None)
        if acc is None:
            acc = frame.astype(np.float64)
        else:
            acc += frame.astype(np.float64)
        _progress.update(_rw_i)
    return acc / n_rw


# ── Fast ensemble accumulator ─────────────────────────────────────────────────

def _coords_to_inv_contact_f32(coords: np.ndarray, eps: float = 1e-3) -> np.ndarray:
    """Compute 1/(||r_i - r_j|| + ε) in float32 without materialising a
    separate distance matrix.

    Uses the Gram-matrix identity:
        D²_ij = ||r_i||² + ||r_j||² - 2 · r_i · r_j

    All operations run in float32 and are done in-place so peak memory is
    exactly one N×N float32 array plus the small coordinate buffer.

    Parameters
    ----------
    coords : (N, 3) float array (any dtype — cast internally)
    eps : float

    Returns
    -------
    inv : (N, N) float32
    """
    X  = np.asarray(coords, dtype=np.float32)
    sq = np.einsum("ij,ij->i", X, X)          # (N,) — squared norms
    # D² = sq_i + sq_j - 2 X X^T  (written in-place into the output buffer)
    inv = np.dot(X, X.T)                       # (N, N) — Gram matrix
    inv *= -2.0
    inv += sq[:, None]
    inv += sq[None, :]
    np.maximum(inv, 0.0, out=inv)              # numerical safety (may go tiny-negative)
    np.sqrt(inv, out=inv)                      # now D
    inv += eps                                 # D + ε
    np.reciprocal(inv, out=inv)               # 1 / (D + ε)
    return inv


def _coords_to_kernel_contact_f32(
    coords: np.ndarray,
    p_func,
    sep_matrix: np.ndarray | None = None,
    dtype=np.float32,
) -> np.ndarray:
    """Compute the contact-probability proxy P(r_ij) using the *same*
    distance→probability kernel that defines c_ij in the Hi-C cross-entropy
    force (see :func:`hic_force.get_kernel_p_func`), instead of a generic
    1/(d+ε) proxy.

    This makes validation numerically consistent with what the force
    actually targets: if the force was built with, say, the "rouse" kernel,
    validation now also measures contacts through P_rouse(r, sep) rather
    than through an unrelated inverse-distance heuristic.

    Uses the same Gram-matrix trick as :func:`_coords_to_inv_contact_f32` to
    get the pairwise distance matrix D without ever materialising a separate
    distance buffer, then applies ``p_func`` elementwise.

    Parameters
    ----------
    coords     : (N, 3) float array (any dtype — cast internally)
    p_func     : callable(r) -> P, or callable(r, sep) -> P when the kernel
                 needs genomic separation (e.g. "rouse")
    sep_matrix : (N, N) genomic-separation matrix (sep[i,j] = |i-j|), only
                 required when the kernel needs it (``needs_sep=True``)
    dtype      : output dtype (default float32, matching the distance calc)

    Returns
    -------
    P : (N, N) ndarray, contact-probability proxy
    """
    X  = np.asarray(coords, dtype=np.float32)
    sq = np.einsum("ij,ij->i", X, X)          # (N,) — squared norms
    D  = np.dot(X, X.T)                       # (N, N) — Gram matrix
    D *= -2.0
    D += sq[:, None]
    D += sq[None, :]
    np.maximum(D, 0.0, out=D)                 # numerical safety
    np.sqrt(D, out=D)                         # now D = pairwise distances

    if sep_matrix is not None:
        P = p_func(D, sep_matrix)
    else:
        P = p_func(D)
    return np.asarray(P, dtype=dtype)


def _accumulate_ensemble(
    cif_paths: list,
    p_func,
    needs_sep: bool = False,
    sep_matrix: np.ndarray | None = None,
    log=None,
) -> np.ndarray:
    """Streaming accumulation of the kernel-contact ensemble average.

    Converts each frame's distances to contact probability via the *same*
    P(r) kernel used to build the Hi-C force (``p_func``), then averages
    across frames — consistent with how the random-walk baseline and the
    force itself interpret distances.

    Memory cost: one N×N accumulator  +  one float32 N×N working buffer
    (both reused across frames — no extra allocations inside the loop).

    Returns
    -------
    avg : (N, N) float64  (upcast at the end for downstream precision)
    """
    _log = log or logger
    _log.info("Streaming kernel-contact accumulation over %d frames …", len(cif_paths))

    acc: np.ndarray | None = None
    n_total = len(cif_paths)
    _progress = ProgressLogger(n_total, _log, "Frames accumulated")

    for i, path in enumerate(cif_paths):
        coords = get_coordinates_cif(path)          # (N, 3) float64
        frame  = _coords_to_kernel_contact_f32(
            coords, p_func, sep_matrix if needs_sep else None,
        )
        if acc is None:
            acc = frame.astype(np.float64)           # first frame
        else:
            acc += frame                             # in-place accumulation
        _progress.update(i)

    if acc is None:
        raise ValueError("No frames were loaded.")

    acc /= len(cif_paths)                           # normalise in-place
    return acc


# ── Computer-vision similarity metrics ───────────────────────────────────────

def _ssim(a: np.ndarray, b: np.ndarray, win: int = 7) -> float:
    """Structural Similarity Index (SSIM) between two 2-D matrices.

    A sliding-window measure that decomposes similarity into luminance,
    contrast, and structure components.  Unlike Pearson r it is sensitive
    to *local* spatial patterns — a TAD shifted by a few bins lowers SSIM
    even if the global statistics are unchanged.

    Range: [-1, 1],  1 = identical.
    """
    try:
        from skimage.metrics import structural_similarity
        # Normalise using the *reference* matrix b's robust range (1st–99th
        # percentile).  Joint min/max lets outliers in either matrix squish
        # both images into a tiny common range, making a flat RW map look
        # nearly identical to a structured sim map → inflated RW SSIM.
        # By anchoring the range to b (experiment/reference) and clipping a
        # (simulation/RW), we preserve the actual dynamic range difference.
        lo   = np.percentile(b, 1)
        hi   = np.percentile(b, 99)
        span = hi - lo
        if span < 1e-12:
            span = 1e-12
        a_n = np.clip((a - lo) / span, 0.0, 1.0)
        b_n = np.clip((b - lo) / span, 0.0, 1.0)
        val, _ = structural_similarity(
            a_n, b_n,
            data_range=1.0,
            win_size=win | 1,     # must be odd
            full=True,
        )
        return float(val)
    except ImportError:
        logger.warning("scikit-image not available — SSIM skipped.")
        return float("nan")


def _gmsd(a: np.ndarray, b: np.ndarray) -> float:
    """Gradient Magnitude Similarity Deviation (GMSD).

    Measures how well *edges* (TAD boundaries, compartment borders) are
    preserved.  Uses a Prewitt-style finite-difference gradient.
    Lower GMSD = better edge alignment; 0 = identical gradients.

    Reference: Xue et al. (2014) IEEE Trans. Image Process.
    """
    def _grad_mag(m: np.ndarray) -> np.ndarray:
        gx = ndimage.prewitt(m, axis=0)
        gy = ndimage.prewitt(m, axis=1)
        return np.hypot(gx, gy)

    gm_a = _grad_mag(a.astype(float))
    gm_b = _grad_mag(b.astype(float))
    c    = 0.0026  # stability constant (same as original paper)
    gms  = (2.0 * gm_a * gm_b + c) / (gm_a ** 2 + gm_b ** 2 + c)
    return float(gms.std())   # deviation of the per-pixel GMS map


def _mutual_information(a: np.ndarray, b: np.ndarray, bins: int = 64) -> float:
    """Normalised Mutual Information between two matrices.

    Robust to monotonic rescalings — captures non-linear dependency that
    Pearson r misses.  Normalised to [0, 1] via NMI = 2·MI / (H_a + H_b).

    Higher is better;  0 = independent,  1 = perfectly predictable.
    """
    a_f = a.ravel().astype(float)
    b_f = b.ravel().astype(float)
    hist2d, _, _ = np.histogram2d(a_f, b_f, bins=bins)
    pxy  = hist2d / hist2d.sum()
    px   = pxy.sum(axis=1, keepdims=True)
    py   = pxy.sum(axis=0, keepdims=True)
    mask = pxy > 0
    mi   = float(np.sum(pxy[mask] * np.log(pxy[mask] / (px * py + 1e-300)[mask])))
    hx   = float(-np.sum(px[px > 0] * np.log(px[px > 0])))
    hy   = float(-np.sum(py[py > 0] * np.log(py[py > 0])))
    denom = hx + hy
    return mi / denom if denom > 0 else 0.0


def _hic_cv_metrics(
    sim_oe: np.ndarray,
    exp_oe: np.ndarray,
) -> dict:
    """Compute all computer-vision similarity metrics on the OE matrices.

    Returns a dict with keys: ssim, gmsd, nmi
    """
    return {
        "ssim": _ssim(sim_oe, exp_oe),
        "gmsd": _gmsd(sim_oe, exp_oe),
        "nmi":  _mutual_information(sim_oe, exp_oe),
    }


# ── Public validation API ─────────────────────────────────────────────────────

def _resolve_validation_kernel(
    coords: np.ndarray,
    kernel: str,
    rc: float,
    alpha: float,
    sigma: float | None,
    sigma_s: float | None,
    kuhn_length: float | None,
    auto_scale: bool,
    auto_scale_percentile: float,
    log,
    r_min: float | None = None,
):
    """Build the P(r) kernel used for validation, optionally auto-calibrating
    its length scale from the structure's own pairwise-distance distribution.

    ``rc`` is the Hi-C force's own microscopic bead-contact scale;
    reused directly here it underflows to ~0 beyond a few bead-spacings,
    giving a near-empty contact-proxy map. When ``auto_scale`` is True
    (default), the scale is instead calibrated from a percentile of the
    structure's real pairwise distances (:func:`hic_force.auto_contact_scale`),
    rescaling every length parameter by the same factor so only the
    absolute scale changes, not the kernel family/shape. Set
    ``auto_scale=False`` to use the literal force parameters instead.

    ``r_min``, when given, is passed straight through to
    :func:`hic_force.get_kernel_p_func` for its dynamic-range calibration
    (see ``hic_force._calibrate_kernel_dynamics``) — an absolute physical
    excluded-volume distance, so it is NOT rescaled by ``factor`` the way
    ``rc``/``sigma``/``sigma_s``/``kuhn_length`` are.
    """
    rc_eff, sigma_eff, sigma_s_eff, kuhn_length_eff = rc, sigma, sigma_s, kuhn_length
    if auto_scale:
        calibrated = auto_contact_scale(coords, percentile=auto_scale_percentile)
        factor = (calibrated / rc) if rc > 1e-12 else 1.0
        rc_eff = calibrated
        sigma_eff = (sigma * factor) if sigma is not None else None
        sigma_s_eff = (sigma_s * factor) if sigma_s is not None else None
        kuhn_length_eff = (kuhn_length * factor) if kuhn_length is not None else None
        log.info(
            "  Auto-calibrated validation contact scale: %.4f nm "
            "(was rc=%.4f nm; %.0fth percentile of sampled pairwise distances)",
            calibrated, rc, auto_scale_percentile,
        )
    return get_kernel_p_func(
        kernel, rc_eff, alpha=alpha, sigma=sigma_eff,
        sigma_s=sigma_s_eff, kuhn_length=kuhn_length_eff, r_min=r_min,
    )


def validate_hic_model(
    cif_path: str,
    hic_matrix: np.ndarray,
    insulation_window: int = 10,
    max_diag: int | None = None,
    save_path: str | None = None,
    n_rw: int = 20,
    rw_step_nm: float = 0.1,
    confine_radius_nm: float | None = None,
    smooth_sigma: float = 1.0,
    kernel: str = "gaussian",
    rc: float = 1.0,
    alpha: float = 3.0,
    sigma: float | None = None,
    sigma_s: float | None = None,
    kuhn_length: float | None = None,
    auto_scale: bool = True,
    auto_scale_percentile: float = 10.0,
    r_min: float | None = None,
    log=None,
) -> dict:
    """Validate a single simulated structure against experimental Hi-C.

    The simulated contact proxy is P(r) evaluated with the same
    kernel/parameters used to build the Hi-C force itself (see
    :func:`hic_force.get_kernel_p_func`), applied identically to the
    simulated structure and the random-walk null baseline. Three
    correlations are computed: diagonal decay (distance-decay law),
    insulation score (TAD boundaries), and PC1 of the O/E matrix (A/B
    compartmentalisation).

    insulation_window: half-width of the insulation sliding window (beads).
    max_diag: diagonals to include in the decay profile (default N//2).
    smooth_sigma: Gaussian smoothing (beads) applied identically to all
    matrices before any metric is computed; 0 disables it.
    kernel/rc/alpha/sigma/sigma_s/kuhn_length: same kernel family as the
    Hi-C force's own parameters. auto_scale (default True) recalibrates the
    absolute length scale from a percentile (auto_scale_percentile, default
    10.0) of the structure's own pairwise distances instead of reusing the
    force's microscopic rc directly, which would underflow to ~0 for
    most pairs — see :func:`hic_force.auto_contact_scale`.

    Caveat on "diagonal decay r": this correlates each structure's raw
    per-separation contact MEAN against experimental Hi-C's, which is
    dominated by the aggregate, bulk distance-decay trend — a property an
    untouched random walk already reproduces reasonably well purely from
    generic (ideal/self-avoiding) chain statistics, since real chromatin's
    own bulk decay law is close to that of an ideal polymer. None of the
    Hi-C force's kernels are designed to target this bulk trend (they add
    LOCAL deviations from it — compartments, TADs, loops), so a random-walk
    baseline tying or even beating the simulated structure on this one
    metric is not evidence the force isn't working. insulation_r, pc1_r and
    pearson_r/spearman_r are computed on the O/E-normalised matrix, which
    divides out that shared bulk trend first — these are the metrics that
    actually isolate the force's contribution, and are where a working Hi-C
    force should clearly separate from the random-walk baseline.

    Returns a dict with diag_decay_r/p, insulation_r/p, pc1_r/p,
    pearson_r/p, spearman_r/p, ssim, gmsd, nmi.
    """
    _log = log or logger

    # Load coordinates first: auto-calibration (below) needs the structure's
    # own pairwise-distance distribution to pick a sensible kernel scale.
    coords = get_coordinates_cif(cif_path)
    N      = len(coords)

    p_func, needs_sep = _resolve_validation_kernel(
        coords, kernel, rc, alpha, sigma, sigma_s, kuhn_length,
        auto_scale, auto_scale_percentile, _log, r_min=r_min,
    )

    # Use the fast Gram-matrix accumulator (float32, in-place, no temp arrays)
    sep_matrix  = np.abs(np.subtract.outer(np.arange(N), np.arange(N))) if needs_sep else None
    sim_contact = _coords_to_kernel_contact_f32(coords, p_func, sep_matrix).astype(np.float64, copy=False)

    # Sanity check: confirm "higher contact strength ⇒ closer 3-D distance"
    # actually holds for this structure's own contact map (the kernel
    # guarantees it analytically pair-by-pair; this verifies nothing broke
    # that guarantee downstream — see _check_contact_distance_monotonicity).
    _check_contact_distance_monotonicity(coords, sim_contact, _log, label="simulated structure")

    # Resize experimental matrix to match simulation resolution
    hic_r = _pool_matrix(hic_matrix, N)
    if max_diag is None:
        max_diag = N // 2

    # Smooth sim and exp identically before any derived metric — suppresses
    # single-bead noise from the heavy-tailed 1/(d+ε) proxy without erasing
    # TAD/compartment-scale structure. See _smooth_matrix.
    sim_contact = _smooth_matrix(sim_contact, smooth_sigma)
    hic_r       = _smooth_matrix(hic_r,       smooth_sigma)

    # ── 1. Diagonal decay ─────────────────────────────────────────────────────
    sim_decay = diagonal_decay_profile(sim_contact, max_diag)
    exp_decay = diagonal_decay_profile(hic_r,       max_diag)
    r_dd, p_dd = _pearson(sim_decay, exp_decay)

    # OE-normalize both matrices once — removes shared distance-decay baseline
    # so all subsequent metrics reflect structural features, not trivial decay.
    sim_oe = oe_matrix(sim_contact)
    exp_oe = oe_matrix(hic_r)

    # ── 2. Insulation score ───────────────────────────────────────────────────
    sim_ins = insulation_score(sim_oe, insulation_window)
    exp_ins = insulation_score(exp_oe, insulation_window)
    r_ins, p_ins = _pearson(sim_ins, exp_ins)

    # ── 3. PC1 (A/B compartments) ─────────────────────────────────────────────
    # Eigenvectors are defined only up to sign; use |r| (abs Pearson between
    # the raw PC1 vectors) instead of correlating |PC1|, which loses
    # compartment-identity information.
    sim_pc1 = pc1_of_oe(sim_contact)
    exp_pc1 = pc1_of_oe(hic_r)
    r_pc1_raw, p_pc1 = _pearson(sim_pc1, exp_pc1)
    r_pc1 = abs(r_pc1_raw)

    # ── 4. Direct matrix correlations (OE) ───────────────────────────────────
    sim_flat = _upper_tri(sim_oe)
    exp_flat = _upper_tri(exp_oe)
    r_pearson,  p_pearson  = _pearson(sim_flat, exp_flat)
    r_spearman, p_spearman = _spearman(sim_flat, exp_flat)

    # ── 5. Computer-vision structural similarity (OE) ─────────────────────────
    cv = _hic_cv_metrics(sim_oe, exp_oe)

    # ── 6. Random-walk null-model baseline ────────────────────────────────────
    # Generates n_rw independent random-walk chains with the same bond length
    # as the simulation, averages their contact maps, and computes the same
    # metrics — giving a physically grounded lower-bound for each score.
    try:
        rw_contact = _rw_baseline_contact(
            N, p_func, needs_sep=needs_sep, sep_matrix=sep_matrix,
            n_rw=n_rw, step_nm=rw_step_nm,
            confine_radius_nm=confine_radius_nm,
            log=_log,
        )
        rw_contact = _smooth_matrix(rw_contact, smooth_sigma)
        rw_oe      = oe_matrix(rw_contact)

        rw_decay            = diagonal_decay_profile(rw_contact, max_diag)
        r_dd_rw,  _         = _pearson(rw_decay, exp_decay)

        rw_ins              = insulation_score(rw_oe, insulation_window)
        r_ins_rw, _         = _pearson(rw_ins, exp_ins)

        rw_pc1              = pc1_of_oe(rw_contact)
        r_pc1_rw_raw, _     = _pearson(rw_pc1, exp_pc1)
        r_pc1_rw            = abs(r_pc1_rw_raw)

        rw_flat             = _upper_tri(rw_oe)
        r_pear_rw,  _       = _pearson(rw_flat, exp_flat)
        r_spear_rw, _       = _spearman(rw_flat, exp_flat)

        cv_rw               = _hic_cv_metrics(rw_oe, exp_oe)
        rw_ok = True
    except Exception as _rw_err:
        _log.warning("RW baseline failed: %s", _rw_err)
        rw_ok = False

    def _fmt(sim_val, rw_val, p=None):
        """Format  'sim  (RW: rw)'  or  'sim  (p=…)  (RW: rw)'."""
        rw_s = f"{rw_val:.4f}" if (rw_ok and not isinstance(rw_val, float | None) or rw_ok) else "n/a"
        if p is not None:
            return f"{sim_val:.4f}  (p={p:.2e})  [RW: {rw_s}]"
        return f"{sim_val:.4f}  [RW: {rw_s}]"

    if rw_ok:
        table_rows = [
            ("Metric",              "MultiMM  [RW baseline]"),
            ("Diagonal decay r",    _fmt(r_dd,       r_dd_rw,  p_dd)),
            ("Insulation score r",  _fmt(r_ins,      r_ins_rw, p_ins)),
            ("|PC1| r",             _fmt(r_pc1,      r_pc1_rw, p_pc1)),
            ("Pearson r (OE tri)",  _fmt(r_pearson,  r_pear_rw,  p_pearson)),
            ("Spearman r (OE tri)", _fmt(r_spearman, r_spear_rw, p_spearman)),
            ("SSIM",                f"{cv['ssim']:.4f}  (local pattern)  [RW: {cv_rw['ssim']:.4f}]"),
            ("GMSD ↓ better",       f"{cv['gmsd']:.4f}  (edge deviation)  [RW: {cv_rw['gmsd']:.4f}]"),
            ("NMI",                 f"{cv['nmi']:.4f}  (mutual info)  [RW: {cv_rw['nmi']:.4f}]"),
        ]
    else:
        table_rows = [
            ("Diagonal decay r",    f"{r_dd:.4f}  (p={p_dd:.2e})"),
            ("Insulation score r",  f"{r_ins:.4f}  (p={p_ins:.2e})"),
            ("|PC1| r",             f"{r_pc1:.4f}  (p={p_pc1:.2e})"),
            ("Pearson r (OE tri)",  f"{r_pearson:.4f}  (p={p_pearson:.2e})"),
            ("Spearman r (OE tri)", f"{r_spearman:.4f}  (p={p_spearman:.2e})"),
            ("SSIM",                f"{cv['ssim']:.4f}  (local pattern similarity)"),
            ("GMSD",                f"{cv['gmsd']:.4f}  (boundary edge deviation, ↓ better)"),
            ("NMI",                 f"{cv['nmi']:.4f}  (non-linear mutual information)"),
        ]

    log_table(table_rows, title="Hi-C Validation — single structure", log_fn=_log.info)
    if rw_ok and r_dd_rw >= r_dd and (r_ins > r_ins_rw or r_pc1 > r_pc1_rw or r_pearson > r_pear_rw):
        _log.info(
            "  Note: the random-walk baseline matched or beat the model on diagonal "
            "decay r (%.3f vs %.3f) while the model still leads on insulation/PC1/OE "
            "correlation above — this is expected, not a sign the Hi-C force isn't "
            "working.  Diagonal decay r reflects the bulk/aggregate distance-decay "
            "trend, which a random walk already reproduces well from generic polymer "
            "statistics; the Hi-C force instead targets the LOCAL deviations from "
            "that trend (compartments, TADs, loops) captured by the O/E-normalised "
            "metrics.  See validate_hic_model's docstring for details.",
            r_dd, r_dd_rw,
        )

    if save_path is not None:
        from .plots import plot_hic_comparison, plot_hic_validation_curves
        import os
        plots_dir = os.path.join(save_path, "plots")
        os.makedirs(plots_dir, exist_ok=True)
        plot_hic_comparison(
            sim_contact, hic_r, plots_dir,
            name="hic_comparison_single",
            rw_matrix=rw_contact if rw_ok else None,
        )
        plot_hic_validation_curves(
            sim_decay, exp_decay, rw_decay if rw_ok else None,
            sim_ins, exp_ins, rw_ins if rw_ok else None,
            sim_pc1, exp_pc1, rw_pc1 if rw_ok else None,
            plots_dir, name="hic_validation_curves_single",
            r_dd=r_dd, r_dd_rw=r_dd_rw if rw_ok else None,
            r_ins=r_ins, r_ins_rw=r_ins_rw if rw_ok else None,
            r_pc1=r_pc1, r_pc1_rw=r_pc1_rw if rw_ok else None,
        )

    return {
        "diag_decay_r": r_dd,       "diag_decay_p": p_dd,
        "insulation_r": r_ins,      "insulation_p": p_ins,
        "pc1_r":        r_pc1,      "pc1_p":        p_pc1,
        "pearson_r":    r_pearson,  "pearson_p":    p_pearson,
        "spearman_r":   r_spearman, "spearman_p":   p_spearman,
        "ssim":         cv["ssim"],
        "gmsd":         cv["gmsd"],
        "nmi":          cv["nmi"],
    }


def validate_hic_ensemble(
    cif_paths: list,
    hic_matrix: np.ndarray,
    insulation_window: int = 10,
    max_diag: int | None = None,
    save_path: str | None = None,
    n_rw: int = 20,
    rw_step_nm: float = 0.1,
    confine_radius_nm: float | None = None,
    smooth_sigma: float = 1.0,
    kernel: str = "gaussian",
    rc: float = 1.0,
    alpha: float = 3.0,
    sigma: float | None = None,
    sigma_s: float | None = None,
    kuhn_length: float | None = None,
    auto_scale: bool = True,
    auto_scale_percentile: float = 10.0,
    r_min: float | None = None,
    log=None,
) -> dict:
    """Validate an ensemble of simulated structures against experimental Hi-C.

    The simulated contact proxy is the average of the per-frame
    contact-probability P(r) across all trajectory frames:
    C_sim = mean_k[ P(D_k) ], using the same kernel as the Hi-C force (see
    :func:`hic_force.get_kernel_p_func`). The random-walk null baseline is
    evaluated through the same kernel. Metrics are the same three as
    ``validate_hic_model`` (diagonal decay, insulation, PC1 correlations).

    cif_paths: per-frame CIF paths from the MD trajectory. insulation_window,
    max_diag, smooth_sigma, kernel/rc/alpha/sigma/sigma_s/kuhn_length,
    auto_scale, auto_scale_percentile: same meaning as in
    ``validate_hic_model`` (auto-calibration here uses the ensemble's first
    frame as the representative structure).

    Caveat on "diagonal decay r" — see ``validate_hic_model``'s docstring:
    this metric is dominated by the bulk/aggregate distance-decay trend,
    which an untouched random walk already reproduces well from generic
    polymer statistics alone, and which none of the Hi-C force's kernels
    specifically target (they add LOCAL deviations on top — compartments,
    TADs, loops). A random-walk baseline tying or beating the simulated
    ensemble here is not evidence the force isn't working; insulation_r,
    pc1_r and pearson_r/spearman_r (computed on the O/E-normalised matrix,
    which divides out that shared bulk trend) are what actually isolate the
    force's contribution.

    Returns a dict with diag_decay_r/p, insulation_r/p, pc1_r/p,
    pearson_r/p, spearman_r/p, ssim, gmsd, nmi.
    """
    _log = log or logger

    # Peek at bead count + coords from the first frame so the separation
    # matrix (when the kernel needs one) can be built once and reused across
    # all frames, and so auto-calibration (below) has a representative
    # pairwise-distance distribution to calibrate against.
    _first_coords = get_coordinates_cif(cif_paths[0])
    N             = len(_first_coords)

    p_func, needs_sep = _resolve_validation_kernel(
        _first_coords, kernel, rc, alpha, sigma, sigma_s, kuhn_length,
        auto_scale, auto_scale_percentile, _log, r_min=r_min,
    )

    sep_matrix    = np.abs(np.subtract.outer(np.arange(N), np.arange(N))) if needs_sep else None

    # ── Fast streaming accumulation (float32, in-place, no temp arrays) ───────
    # Uses the Gram-matrix identity D²_ij = ||r_i||² + ||r_j||² − 2 r_i·r_j
    # so no separate N×N distance buffer is ever materialised; P(r) is then
    # evaluated with the same kernel used to build the Hi-C force.
    inv_avg = _accumulate_ensemble(
        cif_paths, p_func, needs_sep=needs_sep, sep_matrix=sep_matrix, log=_log,
    )

    # Sanity check against a representative frame (the first one): confirms
    # "higher contact strength ⇒ closer 3-D distance" survived the
    # ensemble-average-of-per-frame-heatmaps step, not just the single-frame
    # kernel guarantee — see _check_contact_distance_monotonicity.
    _check_contact_distance_monotonicity(
        _first_coords, inv_avg, _log, label="ensemble (vs. frame 1)",
    )

    hic_r = _pool_matrix(hic_matrix, N)
    if max_diag is None:
        max_diag = N // 2

    # Smooth sim and exp identically before any derived metric — suppresses
    # single-bead noise from the heavy-tailed 1/(d+ε) proxy without erasing
    # TAD/compartment-scale structure. See _smooth_matrix.
    inv_avg = _smooth_matrix(inv_avg, smooth_sigma)
    hic_r   = _smooth_matrix(hic_r,   smooth_sigma)

    # ── 1. Diagonal decay ─────────────────────────────────────────────────────
    sim_decay = diagonal_decay_profile(inv_avg, max_diag)
    exp_decay = diagonal_decay_profile(hic_r,   max_diag)
    r_dd, p_dd = _pearson(sim_decay, exp_decay)

    # OE-normalize both matrices once — removes shared distance-decay baseline
    # so all subsequent metrics reflect structural features, not trivial decay.
    sim_oe = oe_matrix(inv_avg)
    exp_oe = oe_matrix(hic_r)

    # ── 2. Insulation score ───────────────────────────────────────────────────
    sim_ins = insulation_score(sim_oe, insulation_window)
    exp_ins = insulation_score(exp_oe, insulation_window)
    r_ins, p_ins = _pearson(sim_ins, exp_ins)

    # ── 3. PC1 (A/B compartments) ─────────────────────────────────────────────
    sim_pc1 = pc1_of_oe(inv_avg)
    exp_pc1 = pc1_of_oe(hic_r)
    r_pc1_raw, p_pc1 = _pearson(sim_pc1, exp_pc1)
    r_pc1 = abs(r_pc1_raw)

    # ── 4. Direct matrix correlations (OE upper triangle) ────────────────────
    sim_flat = _upper_tri(sim_oe)
    exp_flat = _upper_tri(exp_oe)
    r_pearson,  p_pearson  = _pearson(sim_flat, exp_flat)
    r_spearman, p_spearman = _spearman(sim_flat, exp_flat)

    # ── 5. Computer-vision structural similarity (OE) ─────────────────────────
    cv = _hic_cv_metrics(sim_oe, exp_oe)

    # ── 6. Random-walk null-model baseline ────────────────────────────────────
    try:
        rw_contact = _rw_baseline_contact(
            N, p_func, needs_sep=needs_sep, sep_matrix=sep_matrix,
            n_rw=n_rw, step_nm=rw_step_nm,
            confine_radius_nm=confine_radius_nm,
            log=_log,
        )
        rw_contact = _smooth_matrix(rw_contact, smooth_sigma)
        rw_oe      = oe_matrix(rw_contact)

        rw_decay            = diagonal_decay_profile(rw_contact, max_diag)
        r_dd_rw,  _         = _pearson(rw_decay, exp_decay)

        rw_ins              = insulation_score(rw_oe, insulation_window)
        r_ins_rw, _         = _pearson(rw_ins, exp_ins)

        rw_pc1              = pc1_of_oe(rw_contact)
        r_pc1_rw_raw, _     = _pearson(rw_pc1, exp_pc1)
        r_pc1_rw            = abs(r_pc1_rw_raw)

        rw_flat             = _upper_tri(rw_oe)
        r_pear_rw,  _       = _pearson(rw_flat, exp_flat)
        r_spear_rw, _       = _spearman(rw_flat, exp_flat)

        cv_rw               = _hic_cv_metrics(rw_oe, exp_oe)
        rw_ok = True
    except Exception as _rw_err:
        _log.warning("RW baseline failed: %s", _rw_err)
        rw_ok = False

    def _fmt(sim_val, rw_val, p=None):
        rw_s = f"{rw_val:.4f}" if rw_ok else "n/a"
        if p is not None:
            return f"{sim_val:.4f}  (p={p:.2e})  [RW: {rw_s}]"
        return f"{sim_val:.4f}  [RW: {rw_s}]"

    if rw_ok:
        table_rows = [
            ("Frames averaged",     str(len(cif_paths))),
            ("Metric",              "MultiMM  [RW baseline]"),
            ("Diagonal decay r",    _fmt(r_dd,       r_dd_rw,   p_dd)),
            ("Insulation score r",  _fmt(r_ins,      r_ins_rw,  p_ins)),
            ("|PC1| r",             _fmt(r_pc1,      r_pc1_rw,  p_pc1)),
            ("Pearson r (OE tri)",  _fmt(r_pearson,  r_pear_rw,  p_pearson)),
            ("Spearman r (OE tri)", _fmt(r_spearman, r_spear_rw, p_spearman)),
            ("SSIM",                f"{cv['ssim']:.4f}  (local pattern)  [RW: {cv_rw['ssim']:.4f}]"),
            ("GMSD ↓ better",       f"{cv['gmsd']:.4f}  (edge deviation)  [RW: {cv_rw['gmsd']:.4f}]"),
            ("NMI",                 f"{cv['nmi']:.4f}  (mutual info)  [RW: {cv_rw['nmi']:.4f}]"),
        ]
    else:
        table_rows = [
            ("Frames averaged",     str(len(cif_paths))),
            ("Diagonal decay r",    f"{r_dd:.4f}  (p={p_dd:.2e})"),
            ("Insulation score r",  f"{r_ins:.4f}  (p={p_ins:.2e})"),
            ("|PC1| r",             f"{r_pc1:.4f}  (p={p_pc1:.2e})"),
            ("Pearson r (OE tri)",  f"{r_pearson:.4f}  (p={p_pearson:.2e})"),
            ("Spearman r (OE tri)", f"{r_spearman:.4f}  (p={p_spearman:.2e})"),
            ("SSIM",                f"{cv['ssim']:.4f}  (local pattern similarity)"),
            ("GMSD",                f"{cv['gmsd']:.4f}  (boundary edge deviation, ↓ better)"),
            ("NMI",                 f"{cv['nmi']:.4f}  (non-linear mutual information)"),
        ]

    log_table(table_rows, title="Hi-C Validation — ensemble", log_fn=_log.info)
    if rw_ok and r_dd_rw >= r_dd and (r_ins > r_ins_rw or r_pc1 > r_pc1_rw or r_pearson > r_pear_rw):
        _log.info(
            "  Note: the random-walk baseline matched or beat the model on diagonal "
            "decay r (%.3f vs %.3f) while the model still leads on insulation/PC1/OE "
            "correlation above — this is expected, not a sign the Hi-C force isn't "
            "working.  Diagonal decay r reflects the bulk/aggregate distance-decay "
            "trend, which a random walk already reproduces well from generic polymer "
            "statistics; the Hi-C force instead targets the LOCAL deviations from "
            "that trend (compartments, TADs, loops) captured by the O/E-normalised "
            "metrics.  See validate_hic_ensemble's docstring for details.",
            r_dd, r_dd_rw,
        )

    if save_path is not None:
        from .plots import plot_hic_comparison, plot_hic_validation_curves
        import os
        plots_dir = os.path.join(save_path, "plots")
        os.makedirs(plots_dir, exist_ok=True)
        plot_hic_comparison(
            inv_avg, hic_r, plots_dir,
            name="hic_comparison_ensemble",
            rw_matrix=rw_contact if rw_ok else None,
        )
        plot_hic_validation_curves(
            sim_decay, exp_decay, rw_decay if rw_ok else None,
            sim_ins, exp_ins, rw_ins if rw_ok else None,
            sim_pc1, exp_pc1, rw_pc1 if rw_ok else None,
            plots_dir, name="hic_validation_curves_ensemble",
            r_dd=r_dd, r_dd_rw=r_dd_rw if rw_ok else None,
            r_ins=r_ins, r_ins_rw=r_ins_rw if rw_ok else None,
            r_pc1=r_pc1, r_pc1_rw=r_pc1_rw if rw_ok else None,
        )

    return {
        "diag_decay_r": r_dd,       "diag_decay_p": p_dd,
        "insulation_r": r_ins,      "insulation_p": p_ins,
        "pc1_r":        r_pc1,      "pc1_p":        p_pc1,
        "pearson_r":    r_pearson,  "pearson_p":    p_pearson,
        "spearman_r":   r_spearman, "spearman_p":   p_spearman,
        "ssim":         cv["ssim"],
        "gmsd":         cv["gmsd"],
        "nmi":          cv["nmi"],
    }


# ── Input-vs-output diagnostics: loops (bedpe) and compartments (bed) ───────

def validate_loops(
    cif_path: str,
    ms,
    ns,
    save_path: str | None = None,
    name: str = "loops",
    n_background: int = 10,
    seed: int = 0,
    log=None,
) -> dict:
    """Diagnostic check: did the input loop anchors (from a .bedpe) actually
    end up close together in the output 3D structure?

    For every input loop (m, n) the 3D distance between beads m and n is
    compared against a size-matched background: for each loop, `n_background`
    random bead pairs with the *same genomic separation* |n-m| are sampled,
    so the comparison isolates "are loop anchors closer than generic pairs
    at the same genomic distance" rather than just "are nearby beads close"
    (which would be true of any polymer).

    Returns a dict with the median loop / background distance, the fold
    enrichment, and a Mann-Whitney U test p-value; also writes one
    consolidated figure when `save_path` is given.
    """
    _log = log or logger

    coords = get_coordinates_cif(cif_path)
    N = len(coords)

    ms = np.asarray(ms, dtype=int)
    ns = np.asarray(ns, dtype=int)
    valid = (ms >= 0) & (ns >= 0) & (ms < N) & (ns < N) & (ms != ns)
    ms, ns = ms[valid], ns[valid]

    if len(ms) == 0:
        _log.warning("validate_loops: no valid loop anchors in range, skipping.")
        return {}

    loop_dist = np.linalg.norm(coords[ms] - coords[ns], axis=1)
    loop_sep = np.abs(ns - ms)

    rng = np.random.default_rng(seed)
    bg_dist = np.empty(len(ms) * n_background)
    bg_sep = np.empty_like(bg_dist)
    k = 0
    for sep in loop_sep:
        if sep >= N:
            continue
        i0 = rng.integers(0, N - sep, size=n_background)
        j0 = i0 + sep
        bg_dist[k:k + n_background] = np.linalg.norm(coords[i0] - coords[j0], axis=1)
        bg_sep[k:k + n_background] = sep
        k += n_background
    bg_dist = bg_dist[:k]
    bg_sep = bg_sep[:k]

    median_loop = float(np.median(loop_dist))
    median_bg = float(np.median(bg_dist))
    fold_closer = median_bg / median_loop if median_loop > 0 else float("nan")

    try:
        _, p_value = mannwhitneyu(loop_dist, bg_dist, alternative="less")
        p_value = float(p_value)
    except Exception:
        p_value = float("nan")

    _log.info(
        f"Loop validation ({len(ms)} loops): median loop distance={median_loop:.4f}, "
        f"median background={median_bg:.4f}, fold-closer={fold_closer:.2f}x, "
        f"p(loops closer)={p_value:.2e}"
    )

    if save_path is not None:
        from .plots import plot_loop_validation
        import os
        plots_dir = os.path.join(save_path, "plots")
        os.makedirs(plots_dir, exist_ok=True)
        plot_loop_validation(
            loop_dist, bg_dist, loop_sep, bg_sep, plots_dir,
            name=name, fold_closer=fold_closer, p_value=p_value,
        )

    return {
        "n_loops": int(len(ms)),
        "median_loop_dist": median_loop,
        "median_bg_dist": median_bg,
        "fold_closer": fold_closer,
        "p_value": p_value,
    }


def validate_compartments(
    cif_path: str,
    Cs,
    save_path: str | None = None,
    name: str = "compartments",
    eps: float = 1e-3,
    log=None,
) -> dict:
    """Diagnostic check: does the structure's own derived A/B signal (PC1 of
    its OE-normalised contact matrix) agree with the input compartment track
    (from a .bed, via COMPARTMENT_PATH)?

    This is a self-consistency check between what was *fed in* as a
    structural restraint (or as context) and what the resulting 3D structure
    *actually looks like* once you call compartments on it the same way you
    would on real Hi-C data.
    """
    _log = log or logger

    coords = get_coordinates_cif(cif_path)
    N = len(coords)

    Cs = np.asarray(Cs)[:N]
    if len(Cs) < N:
        _log.warning("validate_compartments: compartment track shorter than structure, skipping.")
        return {}

    contact = _coords_to_inv_contact_f32(coords, eps=eps).astype(np.float64, copy=False)
    pc1 = pc1_of_oe(contact)

    valid = Cs != 0
    if valid.sum() < 10:
        _log.warning("validate_compartments: not enough labelled beads, skipping.")
        return {}

    r, p = _pearson(pc1[valid], Cs[valid].astype(float))

    # PC1 sign is arbitrary; align it to the input track for display and for
    # a sign-based "did we call the right compartment" accuracy metric.
    pc1_aligned = pc1 if (np.isnan(r) or r >= 0) else -pc1
    r_abs = abs(r) if not np.isnan(r) else float("nan")

    is_a = Cs[valid] > 0
    is_b = Cs[valid] < 0
    sign_call = pc1_aligned[valid] > 0
    accuracy = float(np.mean(sign_call == is_a))
    sens_a = float(np.mean(sign_call[is_a])) if is_a.sum() else float("nan")
    sens_b = float(np.mean(~sign_call[is_b])) if is_b.sum() else float("nan")

    _log.info(
        f"Compartment validation ({int(valid.sum())} labelled beads): "
        f"|r|={r_abs:.3f}, sign-call accuracy={accuracy:.1%} "
        f"(A sensitivity={sens_a:.1%}, B sensitivity={sens_b:.1%})"
    )

    if save_path is not None:
        from .plots import plot_compartment_validation
        import os
        plots_dir = os.path.join(save_path, "plots")
        os.makedirs(plots_dir, exist_ok=True)
        plot_compartment_validation(
            pc1_aligned[valid], Cs[valid], plots_dir,
            name=name, r=r_abs, accuracy=accuracy, sens_a=sens_a, sens_b=sens_b,
        )

    return {
        "n_beads": int(valid.sum()),
        "pc1_r": r_abs,
        "pc1_p": p,
        "sign_accuracy": accuracy,
        "a_sensitivity": sens_a,
        "b_sensitivity": sens_b,
    }


# =============================================================================
# Distance vs. experimental Hi-C strength — does the force actually do what
# it is supposed to: pull high-strength (enriched) pairs to small distance,
# and leave low-strength (background/depleted) pairs alone?
# =============================================================================

def validate_distance_vs_strength(
    cif_path: str,
    hic_matrix: np.ndarray,
    save_path: str | None = None,
    name: str = "distance_vs_strength",
    oe_normalize: bool = True,
    n_pairs: int = 20000,
    min_sep: int = 2,
    low_pct: float = 25.0,
    high_pct: float = 75.0,
    seed: int = 0,
    log=None,
) -> dict:
    """Check that pairs with higher experimental Hi-C strength end up at a
    smaller 3-D distance in the output structure — the basic physical claim
    the Hi-C force is built to enforce (see hic_force._oe_normalize_matrix).

    ``strength`` is the same c_ij the force targets: the experimental
    matrix run through :func:`hic_force.preprocess_hic_matrix` with the
    same ``oe_normalize`` setting as the force (``HIC_FORCE_OE``). With
    ``oe_normalize=True`` (default), strength = enrichment above the
    distance-decay background, floored at 0 — so "low strength" mostly
    means "at or below background" (the displayed heatmap's "towards -1"),
    not just "a smaller positive number".

    Two complementary checks on a random sample of bead pairs:
      * Classification — split pairs into a low-strength group
        (<= ``low_pct`` percentile) and a high-strength group
        (>= ``high_pct`` percentile), then test whether the high-strength
        group's distances are stochastically smaller (one-sided
        Mann-Whitney U). The U-statistic, rescaled, is an AUC: the
        probability that a random low-strength pair is farther apart than
        a random high-strength pair (0.5 = no better than chance, 1.0 =
        perfectly separated).
      * Regression — Pearson/Spearman correlation between strength and
        distance across all sampled pairs (expected negative: higher
        strength, shorter distance).

    Returns a dict with auc, mannwhitney_p, pearson_r/p, spearman_r/p,
    n_pairs, and the group sizes; also saves a two-panel diagnostic plot
    under ``save_path/plots`` when ``save_path`` is given.
    """
    _log = log or logger

    coords = get_coordinates_cif(cif_path)
    N = len(coords)
    H = preprocess_hic_matrix(hic_matrix, N, already_balanced=False, oe_normalize=oe_normalize)

    rng = np.random.default_rng(seed)
    rows, cols = np.triu_indices(N, k=max(1, min_sep))
    if rows.size == 0:
        _log.warning("validate_distance_vs_strength: no eligible pairs (N too small?)")
        return {}
    if rows.size > n_pairs:
        sel = rng.choice(rows.size, size=n_pairs, replace=False)
        rows, cols = rows[sel], cols[sel]

    strength = H[rows, cols]
    dist = np.linalg.norm(coords[rows] - coords[cols], axis=1)
    ok = np.isfinite(strength) & np.isfinite(dist)
    strength, dist = strength[ok], dist[ok]
    n_used = int(ok.sum())

    lo_thr = float(np.percentile(strength, low_pct))
    hi_thr = float(np.percentile(strength, high_pct))
    low_mask = strength <= lo_thr
    high_mask = strength >= hi_thr
    low_dist, high_dist = dist[low_mask], dist[high_mask]

    if low_dist.size >= 3 and high_dist.size >= 3:
        u_stat, p_mw = mannwhitneyu(low_dist, high_dist, alternative="greater")
        auc = float(u_stat / (low_dist.size * high_dist.size))
    else:
        auc, p_mw = float("nan"), float("nan")

    r_pearson, p_pearson = _pearson(strength, dist)
    rho_spearman, p_spearman = _spearman(strength, dist)

    _log.info(
        "Distance-vs-strength validation (%d pairs): AUC=%.3f (p=%.1e), "
        "Pearson r=%.3f (p=%.1e), Spearman ρ=%.3f (p=%.1e)",
        n_used, auc, p_mw, r_pearson, p_pearson, rho_spearman, p_spearman,
    )

    results = {
        "n_pairs": n_used,
        "low_group_n": int(low_dist.size),
        "high_group_n": int(high_dist.size),
        "auc": auc,
        "mannwhitney_p": p_mw,
        "pearson_r": r_pearson,
        "pearson_p": p_pearson,
        "spearman_r": rho_spearman,
        "spearman_p": p_spearman,
    }

    if save_path is not None:
        from .plots import plot_distance_vs_strength
        import os
        plots_dir = os.path.join(save_path, "plots")
        os.makedirs(plots_dir, exist_ok=True)
        plot_distance_vs_strength(
            strength, dist, low_mask, high_mask, results, plots_dir, name=name,
        )

    return results
