#!/usr/bin/env python3
"""
read_hic.py  —  Hi-C contact matrix loader for MultiMM / OpenMM
================================================================

Reads Hi-C data in three formats:

    .hic    — Juicer / Aiden-lab format (requires hicstraw)
    .cool   — single-resolution Cooler format (requires cooler)
    .mcool  — multi-resolution Cooler format (requires cooler)

Given a desired number of simulation beads N_beads, the script
automatically selects the best available resolution so that the
loaded matrix is as close to N_beads × N_beads as possible.
If the native bin count does not match exactly, a weighted
average-pooling step aggregates / interpolates to produce exactly
an N_beads × N_beads matrix.

Public API
----------
    list_resolutions(path)                      → dict
    get_chromosome_sizes(path)                  → dict
    choose_resolution(available, length, N)     → int
    read_hic_matrix(path, chrom, N_beads, ...)  → (np.ndarray, dict)
    plot_hic_matrix(H, ...)                     → matplotlib Figure

Install dependencies
--------------------
    pip install hic-straw cooler numpy scipy matplotlib

Example (bottom of this file)
------------------------------
    Run directly:  python read_hic.py
"""

from __future__ import annotations

import logging
import os
import pathlib
import warnings
from typing import Dict, List, Optional, Tuple

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter, median_filter

try:
    from logger import setup_logger
    setup_logger()
except ImportError:
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] %(levelname)-8s %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )

log = logging.getLogger(__name__)

from .utils import oe_matrix  # noqa: E402 — canonical OE-normalisation step


# ═════════════════════════════════════════════════════════════════════════════
# Layer 0 — Format detection & name normalisation
# ═════════════════════════════════════════════════════════════════════════════

def _detect_format(path: str) -> str:
    """Return 'hic', 'cool', or 'mcool' based on file extension.

    Raises
    ------
    FileNotFoundError  if the file does not exist
    ValueError         if the extension is not recognised
    """
    p = pathlib.Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Hi-C file not found: {path}")

    ext = p.suffix.lower()
    if ext == ".hic":
        return "hic"
    if ext == ".cool":
        return "cool"
    if ext == ".mcool":
        return "mcool"
    raise ValueError(
        f"Unrecognised extension '{ext}'.  "
        "Supported: .hic, .cool, .mcool"
    )


def _normalise_chrom(chrom: str, available: List[str]) -> str:
    """
    Match *chrom* to the actual name stored in the file.

    Handles 'chr1' ↔ '1', 'chrX' ↔ 'X', etc.  Case-insensitive fallback.

    Raises
    ------
    ValueError  if no match is found
    """
    if chrom in available:
        return chrom

    # try adding / stripping 'chr' prefix
    alt = chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"
    if alt in available:
        log.debug("chromosome '%s' normalised → '%s'", chrom, alt)
        return alt

    # case-insensitive search as last resort
    lower_map = {c.lower(): c for c in available}
    for candidate in (chrom.lower(), alt.lower()):
        if candidate in lower_map:
            found = lower_map[candidate]
            log.debug("chromosome '%s' normalised → '%s' (case fold)", chrom, found)
            return found

    raise ValueError(
        f"Chromosome '{chrom}' not found in file.\n"
        f"Available chromosomes: {available}"
    )


# ═════════════════════════════════════════════════════════════════════════════
# Layer 1 — File inspection  (no matrix loaded)
# ═════════════════════════════════════════════════════════════════════════════

def list_resolutions(path: str) -> Dict[str, object]:
    """
    Return the available resolutions (bp) for the given Hi-C file.

    Parameters
    ----------
    path : str  —  path to .hic, .cool or .mcool file

    Returns
    -------
    info : dict with keys
        'format'      : 'hic' | 'cool' | 'mcool'
        'resolutions' : sorted list of int (bp)
    """
    fmt = _detect_format(path)
    log.info("list_resolutions: %s  [%s]", path, fmt)

    if fmt == "hic":
        import hicstraw
        hic  = hicstraw.HiCFile(str(path))
        ress = sorted(int(r) for r in hic.getResolutions())

    elif fmt == "cool":
        import cooler
        c    = cooler.Cooler(str(path))
        ress = [int(c.binsize)]

    else:  # mcool
        import cooler
        keys = cooler.fileops.list_coolers(str(path))
        ress = sorted(
            int(k.split("/")[-1])
            for k in keys
            if k not in ("/", "")
            and k.split("/")[-1].isdigit()
        )
        if not ress:
            raise RuntimeError(
                f"No resolution groups found in {path}.  "
                "Expected paths like '/resolutions/5000'."
            )

    if not ress:
        raise RuntimeError(f"No resolutions found in {path}.")

    log.info("  available: %s bp", ress)
    return {"format": fmt, "resolutions": ress}


def get_chromosome_sizes(path: str) -> Dict[str, int]:
    """
    Return chromosome name → size in bp for the given file.

    Parameters
    ----------
    path : str  —  path to .hic, .cool or .mcool file

    Returns
    -------
    sizes : dict {chrom_name: length_bp}
    """
    fmt = _detect_format(path)

    if fmt == "hic":
        import hicstraw
        hic   = hicstraw.HiCFile(str(path))
        sizes = {
            c.name: int(c.length)
            for c in hic.getChromosomes()
            if c.name.lower() not in ("all", "assembly", "mt", "")
            and c.length > 0
        }

    elif fmt == "cool":
        import cooler
        c     = cooler.Cooler(str(path))
        sizes = {n: int(s) for n, s in zip(c.chromnames, c.chromsizes)}

    else:  # mcool
        import cooler
        keys = cooler.fileops.list_coolers(str(path))
        uri  = f"{path}::{keys[0]}"
        c    = cooler.Cooler(uri)
        sizes = {n: int(s) for n, s in zip(c.chromnames, c.chromsizes)}

    if not sizes:
        raise RuntimeError(f"No chromosomes found in {path}.")

    log.info("get_chromosome_sizes: %d chromosomes", len(sizes))
    return sizes


# ═════════════════════════════════════════════════════════════════════════════
# Layer 2 — Resolution selection
# ═════════════════════════════════════════════════════════════════════════════

def choose_resolution(
    available        : List[int],
    region_length_bp : int,
    N_beads          : int,
) -> int:
    """
    Pick the available resolution that makes the matrix closest to
    N_beads × N_beads (before average-pooling).

    Strategy: minimise |R - target| where target = region_length_bp / N_beads.
    Slightly prefer a finer resolution (more raw bins) so that the pooling
    step aggregates rather than interpolates, which is lossless.

    Parameters
    ----------
    available        : list of int — resolutions in the file (bp)
    region_length_bp : length of the genomic region (bp)
    N_beads          : desired number of simulation beads

    Returns
    -------
    resolution : int (bp)

    Raises
    ------
    ValueError  if region_length_bp <= 0 or N_beads <= 0
    """
    if region_length_bp <= 0:
        raise ValueError(f"region_length_bp must be > 0, got {region_length_bp}")
    if N_beads <= 0:
        raise ValueError(f"N_beads must be > 0, got {N_beads}")
    if not available:
        raise ValueError("available resolution list is empty")

    target  = region_length_bp / N_beads
    # prefer slightly finer (smaller R) so raw n_bins ≥ N_beads when possible
    best    = min(available, key=lambda r: (abs(r - target), r))
    n_bins  = region_length_bp // best

    if n_bins == 0:
        # all resolutions coarser than the region — pick finest
        best   = min(available)
        n_bins = max(region_length_bp // best, 1)
        log.warning(
            "Region (%d bp) is smaller than all available resolutions. "
            "Using finest resolution %d bp → %d bin(s).",
            region_length_bp, best, n_bins,
        )

    log.info(
        "choose_resolution: target %.0f bp → best %d bp → %d raw bins "
        "(requested N_beads=%d)",
        target, best, n_bins, N_beads,
    )
    return best


# ═════════════════════════════════════════════════════════════════════════════
# Layer 3 — Format-specific loaders  (return dense np.ndarray)
# ═════════════════════════════════════════════════════════════════════════════

def _load_from_hic(
    path          : str,
    chrom         : str,
    resolution    : int,
    start         : int,
    end           : int,
    normalization : str = "KR",
) -> np.ndarray:
    """
    Load a dense contact matrix from a .hic file via hicstraw.

    Falls back to 'NONE' normalisation if the requested one is unavailable.

    Parameters
    ----------
    path          : .hic file path
    chrom         : chromosome name as stored in the file
    resolution    : bin size (bp)
    start / end   : genomic coordinates (bp, 0-based, half-open)
    normalization : 'KR', 'SCALE', 'VC', 'VC_SQRT', or 'NONE'

    Returns
    -------
    H : (n_bins, n_bins) dense float64 ndarray, symmetric, diagonal = 0
    """
    import hicstraw

    log.info(
        "_load_from_hic: %s  %d–%d  res=%d  norm=%s",
        chrom, start, end, resolution, normalization,
    )

    hic       = hicstraw.HiCFile(str(path))
    n_bins    = max((end - start + resolution - 1) // resolution, 1)
    H         = np.zeros((n_bins, n_bins), dtype=np.float64)

    norms_to_try = [normalization] if normalization != "NONE" else []
    norms_to_try.append("NONE")

    records = None
    for norm in norms_to_try:
        try:
            mzd     = hic.getMatrixZoomData(
                chrom, chrom, "observed", norm, "BP", resolution
            )
            records = mzd.getRecords(start, end, start, end)
            if norm != normalization:
                log.warning(
                    "  normalisation '%s' unavailable — using '%s'",
                    normalization, norm,
                )
            break
        except Exception as exc:
            log.warning("  norm '%s' failed (%s)", norm, exc)

    if records is None:
        log.error("  all normalisations failed — returning zero matrix")
        return H

    n_loaded = 0
    for rec in records:
        val = float(rec.counts)
        if not np.isfinite(val) or val < 0:
            continue
        i = (rec.binX - start) // resolution
        j = (rec.binY - start) // resolution
        if 0 <= i < n_bins and 0 <= j < n_bins:
            H[i, j] += val
            if i != j:
                H[j, i] += val
            n_loaded += 1

    # symmetrize: take average where both i,j were filled independently
    H = 0.5 * (H + H.T)
    np.fill_diagonal(H, 0.0)

    log.info("  %d records → (%d, %d) matrix", n_loaded, n_bins, n_bins)
    return H


def _load_from_cool(
    path       : str,
    chrom      : str,
    resolution : int,
    start      : int,
    end        : int,
    fmt        : str = "cool",
) -> np.ndarray:
    """
    Load a dense contact matrix from a .cool / .mcool file via cooler.

    Tries KR-balance weights first; falls back to raw counts if unavailable
    or if balancing produces NaN-only rows.

    Parameters
    ----------
    path       : file path
    chrom      : chromosome name as stored in the file
    resolution : bin size — must be an available resolution in the file
    start/end  : genomic coordinates (bp, 0-based, half-open)
    fmt        : 'cool' or 'mcool'

    Returns
    -------
    H : (n_bins, n_bins) dense float64 ndarray, symmetric, diagonal = 0
    """
    import cooler

    uri = path if fmt == "cool" else f"{path}::resolutions/{resolution}"

    log.info("_load_from_cool: %s  %s  %d–%d", uri, chrom, start, end)

    c      = cooler.Cooler(uri)
    region = f"{chrom}:{start}-{end}"

    # ── try balanced, fall back to raw ───────────────────────────────────────
    H = None
    try:
        H_bal = np.array(c.matrix(balance=True).fetch(region), dtype=np.float64)
        nan_frac = np.isnan(H_bal).mean()
        if nan_frac > 0.9:
            log.warning(
                "  KR balance gave %.0f%% NaN — falling back to raw counts",
                100 * nan_frac,
            )
        else:
            H = H_bal
            if nan_frac > 0:
                log.info("  KR balance applied (%.1f%% NaN → 0)", 100 * nan_frac)
    except Exception as exc:
        log.warning("  KR balance unavailable (%s) — using raw counts", exc)

    if H is None:
        H = np.array(c.matrix(balance=False).fetch(region), dtype=np.float64)

    # NaN → 0, symmetrize, zero diagonal
    H = np.nan_to_num(H, nan=0.0, posinf=0.0, neginf=0.0)
    H = 0.5 * (H + H.T)
    np.fill_diagonal(H, 0.0)

    log.info("  loaded (%d, %d) matrix", *H.shape)
    return H


# ═════════════════════════════════════════════════════════════════════════════
# Layer 4 — Post-processing
# ═════════════════════════════════════════════════════════════════════════════

def handle_missing_bins(H: np.ndarray, max_gap: int = 3) -> np.ndarray:
    """
    Replace NaN / Inf and interpolate small consecutive empty-bin gaps.

    Empty bins (zero marginal) arise from unmappable or low-coverage regions.
    Gaps of ≤ max_gap consecutive bins are linearly interpolated from their
    non-empty neighbours; larger gaps are left as zeros.

    Parameters
    ----------
    H       : (N, N) contact matrix (may contain NaN / Inf)
    max_gap : maximum consecutive empty bins to fill by interpolation

    Returns
    -------
    H_clean : (N, N) float64 ndarray, no NaN / Inf
    """
    H = np.array(H, dtype=np.float64)

    # 1. replace non-finite
    bad = ~np.isfinite(H)
    if bad.any():
        log.warning(
            "  handle_missing_bins: %d non-finite entries → 0", int(bad.sum())
        )
        H[bad] = 0.0

    N         = H.shape[0]
    marginal  = H.sum(axis=1)
    empty     = marginal == 0.0
    n_empty   = int(empty.sum())

    if n_empty == 0:
        return H

    log.info(
        "  %d / %d bins empty (%.1f%%), interpolating gaps ≤ %d",
        n_empty, N, 100.0 * n_empty / N, max_gap,
    )

    # 2. interpolate small gaps
    filled = 0
    i = 0
    while i < N:
        if empty[i]:
            j = i
            while j < N and empty[j]:
                j += 1
            gap = j - i
            if gap <= max_gap:
                lo  = max(i - 1, 0)
                hi  = min(j, N - 1)
                if lo == i:      # gap starts at left edge
                    H[i:j, :] = H[hi, :]
                    H[:, i:j] = H[:, hi:hi+1]
                elif hi == j - 1:  # gap ends at right edge
                    H[i:j, :] = H[lo, :]
                    H[:, i:j] = H[:, lo:lo+1]
                else:
                    frac = np.linspace(0.0, 1.0, gap + 2)[1:-1]
                    for k, f in enumerate(frac):
                        row = (1.0 - f) * H[lo, :] + f * H[hi, :]
                        H[i + k, :] = row
                        H[:, i + k] = row
                filled += gap
            i = j
        else:
            i += 1

    # re-enforce symmetry after interpolation
    H = 0.5 * (H + H.T)
    np.fill_diagonal(H, 0.0)

    if filled:
        log.info("  interpolated %d empty bins", filled)

    return H


def _pool_weights(N_in: int, N_out: int) -> np.ndarray:
    """
    Build a (N_out, N_in) weighted resampling matrix.

    W[i, j]  =  fraction of input bin j that overlaps output bin i.
    Row sums equal exactly 1 by normalisation.

    Works for both:
      downsampling  (N_out < N_in) → weighted average pool
      upsampling    (N_out > N_in) → linear interpolation across boundaries

    The pooled matrix is then  H_out = W @ H_in @ W.T
    a single numpy matrix multiply — fast even for N ~ 3000.
    """
    edges = np.linspace(0.0, float(N_in), N_out + 1)
    W     = np.zeros((N_out, N_in), dtype=np.float64)

    for i in range(N_out):
        lo, hi = edges[i], edges[i + 1]
        j0 = int(np.floor(lo))
        j1 = int(np.ceil(hi))
        for j in range(j0, min(j1, N_in)):
            overlap  = min(float(j + 1), hi) - max(float(j), lo)
            W[i, j]  = max(overlap, 0.0)
        row_sum = W[i].sum()
        if row_sum > 0.0:
            W[i] /= row_sum

    return W


def pool_to_n_beads(H: np.ndarray, N_beads: int) -> np.ndarray:
    """
    Resize a contact matrix to **exactly** N_beads × N_beads using a
    single unified weighted-overlap resampling kernel.

    Both directions use the same ``_pool_weights`` matrix:

    * **Downsampling** (N_in > N_beads) — weighted average pool:
      each output bin receives a weighted average of the overlapping input
      bins.  Fractional-bin boundaries are handled correctly, so the result
      is exact for any integer or non-integer ratio.

    * **Upsampling** (N_in < N_beads) — linear interpolation:
      output bins smaller than one input bin receive a weighted combination
      of the (at most two) surrounding input bins.  A warning is logged
      because upsampling does not add real information; consider loading a
      finer resolution instead.

    The pooled matrix is computed as ``W @ H @ W.T`` — two dense matrix
    multiplications, fast even for N ~ 3000.

    Parameters
    ----------
    H       : (N_in, N_in) contact matrix (float64, symmetric, diagonal=0)
    N_beads : target output size

    Returns
    -------
    H_out : (N_beads, N_beads) float64 ndarray, symmetric, diagonal=0

    Guarantees
    ----------
    * Output shape is **always exactly** (N_beads, N_beads).
    * Output is symmetric:  H_out == H_out.T  (up to floating-point).
    * Diagonal is zeroed.
    """
    H    = np.asarray(H, dtype=np.float64)
    N_in = H.shape[0]

    if N_in == N_beads:
        out = H.copy()
        np.fill_diagonal(out, 0.0)
        return out

    if N_in < N_beads:
        log.warning(
            "  pool_to_n_beads: N_in=%d < N_beads=%d — upsampling via "
            "linear interpolation.  Consider a finer file resolution for "
            "better accuracy.",
            N_in, N_beads,
        )
    else:
        log.info(
            "  pool_to_n_beads: %d → %d (weighted average pool, ratio=%.3f)",
            N_in, N_beads, N_in / N_beads,
        )

    W     = _pool_weights(N_in, N_beads)   # (N_beads, N_in)
    H_out = W @ H @ W.T                    # (N_beads, N_beads)

    # enforce exact shape (guard against floating-point edge cases)
    assert H_out.shape == (N_beads, N_beads), (
        f"pool_to_n_beads: unexpected shape {H_out.shape} "
        f"(expected ({N_beads}, {N_beads}))"
    )

    H_out = 0.5 * (H_out + H_out.T)
    np.fill_diagonal(H_out, 0.0)
    return H_out


def symmetrize_and_clean(H: np.ndarray) -> np.ndarray:
    """Symmetrise H, zero NaN/Inf/negatives and the diagonal."""
    H = np.array(H, dtype=np.float64)
    if H.ndim != 2 or H.shape[0] != H.shape[1]:
        raise ValueError(f"Hi-C matrix must be square 2-D, got shape {H.shape}")

    H = 0.5 * (H + H.T)
    n_bad = int(np.count_nonzero(~np.isfinite(H)))
    if n_bad:
        log.warning("symmetrize_and_clean: %d non-finite entries -> 0", n_bad)
    np.nan_to_num(H, copy=False, nan=0.0, posinf=0.0, neginf=0.0)
    np.maximum(H, 0.0, out=H)
    np.fill_diagonal(H, 0.0)
    return H


def diagonal_normalize(H: np.ndarray, n_iter: int = 50, tol: float = 1e-6) -> np.ndarray:
    """Knight-Ruiz iterative row/column balancing: H <- D^-1 H D^-1 until
    marginals are ~1, or `tol` is reached (early stop), up to `n_iter`.
    """
    H = H.copy()
    it = 0
    for it in range(n_iter):
        row_sums = H.sum(axis=1)
        scale    = np.where(row_sums > 0, 1.0 / np.sqrt(row_sums), 0.0)
        H        = scale[:, None] * H * scale[None, :]
        np.nan_to_num(H, copy=False, nan=0.0, posinf=0.0, neginf=0.0)
        if np.abs(H.sum(axis=1) - 1.0).max() <= tol:
            break

    H = np.maximum(H, 0.0)
    np.fill_diagonal(H, 0.0)
    log.debug("diagonal_normalize: %d/%d iterations used", it + 1, n_iter)
    return H


def oe_enrichment_matrix(H: np.ndarray) -> np.ndarray:
    """OE enrichment above background: oe_matrix(H) - 1, floored at 0. Used
    as the Hi-C force's c_ij target when HIC_FORCE_OE=True — pairs at/below
    the distance-decay baseline get c_ij=0 ("loose"); only enriched pairs
    (OE > 1) attract.
    """
    return np.maximum(oe_matrix(H) - 1.0, 0.0)


_DENOISE_MEDIAN_SIZE     = 3     # median-filter window, beads
_DENOISE_SIGMA           = 1.5   # Gaussian smoothing std, beads
_DENOISE_CLIP_PERCENTILE = 99.0  # outlier safety-net cap


def denoise_contact_matrix(
    H: np.ndarray,
    sigma: float = _DENOISE_SIGMA,
    clip_percentile: float = _DENOISE_CLIP_PERCENTILE,
    median_size: int = _DENOISE_MEDIAN_SIZE,
) -> np.ndarray:
    """Denoise a contact matrix: median filter (drops isolated, unsupported
    pixels while real TAD/compartment patches survive) -> percentile cap ->
    light Gaussian smoothing -> re-symmetrize. Used both as the Hi-C force's
    final c_ij matrix and by validation.py for consistent comparisons.
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
    """Full Hi-C matrix prep pipeline: clean -> resize -> balance ->
    (optional) OE-floor -> denoise. This is what `hic_force.build_hic_force`
    runs before building the force, exposed here so diagnostics can see the
    exact c_ij values the force targets.

    save_path/chrom: if given, saves a before/after denoising plot to
    ``<save_path>/plots/hic_preprocessing.png``.
    """
    H = symmetrize_and_clean(H_raw)
    if H.shape[0] != N_beads:
        H = pool_to_n_beads(H, N_beads)
    if not already_balanced:
        H = diagonal_normalize(H)
    if oe_normalize:
        H = oe_enrichment_matrix(H)
        log.info("preprocess_hic_matrix: OE normalisation applied")

    # denoise AFTER OE (if enabled) — that's the matrix the force actually
    # targets, and OE division can amplify shot noise
    H_before_denoise = H.copy() if save_path is not None else None
    H = denoise_contact_matrix(H)
    n_suppressed = (
        int(np.count_nonzero((H_before_denoise > 0) & (H <= 0)))
        if H_before_denoise is not None else None
    )
    log.info(
        "preprocess_hic_matrix: denoised (%s OE) — median %dx%d, cap %.0fpct, "
        "gaussian sigma=%.1f%s",
        "after" if oe_normalize else "without",
        _DENOISE_MEDIAN_SIZE, _DENOISE_MEDIAN_SIZE,
        _DENOISE_CLIP_PERCENTILE, _DENOISE_SIGMA,
        f", suppressed {n_suppressed} px" if n_suppressed is not None else "",
    )

    if save_path is not None:
        try:
            from .plots import plot_hic_preprocessing
            plots_dir = os.path.join(save_path, "plots")
            plot_hic_preprocessing(
                H_before_denoise, H, plots_dir, chrom=chrom,
                sigma=_DENOISE_SIGMA, oe_normalized=oe_normalize,
            )
        except Exception as exc:  # pragma: no cover
            log.warning("Could not save Hi-C preprocessing plot: %s", exc)

    return H


# ═════════════════════════════════════════════════════════════════════════════
# Layer 5 — Public entry point
# ═════════════════════════════════════════════════════════════════════════════

def read_hic_matrix(
    path          : str,
    chrom         : str,
    N_beads       : int,
    region        : Optional[Tuple[int, int]] = None,
    normalization : str  = "KR",
    max_gap       : int  = 3,
    resize        : bool = True,
) -> Tuple[np.ndarray, dict]:
    """
    Load a Hi-C contact matrix and return an N_beads × N_beads array.

    Pipeline
    --------
    1. Detect format and list available resolutions
    2. Determine genomic region (full chromosome or user sub-region)
    3. choose_resolution()    → closest bin size to region_length / N_beads
    4. Load dense matrix      (_load_from_hic or _load_from_cool)
    5. handle_missing_bins()  — NaN / Inf / empty-bin interpolation
    6. pool_to_n_beads()      — weighted average pooling to N_beads × N_beads

    Note: denoising runs later, in ``preprocess_hic_matrix`` — this returns
    the raw (loaded/pooled, not yet denoised) matrix.

    Parameters
    ----------
    path          : path to .hic, .cool, or .mcool file
    chrom         : chromosome name, e.g. 'chr1' or '1'
    N_beads       : desired number of simulation beads
    region        : (start_bp, end_bp), 0-based half-open.  None → full chrom.
    normalization : for .hic files: 'KR', 'SCALE', 'VC', 'VC_SQRT', 'NONE'
                    (cooler uses its own stored balance weights)
    max_gap       : max consecutive empty bins to interpolate (default 3)
    resize        : if True (default), pool/resize to exactly N_beads × N_beads.
                    If False, return the native-resolution matrix.

    Returns
    -------
    H    : (N_beads, N_beads) float64 ndarray  (or native size if resize=False)
    meta : dict
        'chrom'      : chromosome name as stored in the file
        'start'      : region start (bp)
        'end'        : region end (bp)
        'resolution' : chosen resolution (bp)
        'n_bins_raw' : number of bins before resize
        'format'     : 'hic' | 'cool' | 'mcool'

    Examples
    --------
    # Full chromosome, 300 beads
    H, meta = read_hic_matrix('file.hic', 'chr21', N_beads=300)

    # Sub-region
    H, meta = read_hic_matrix('file.hic', 'chr1', N_beads=300,
                               region=(0, 50_000_000))

    # Pass directly to build_hic_force
    from hic_force import build_hic_force
    force = build_hic_force(H, N_beads=300, rc=6.0)
    """
    log.info("═" * 60)
    log.info("read_hic_matrix: %s  chrom=%s  N_beads=%d", path, chrom, N_beads)

    if N_beads <= 0:
        raise ValueError(f"N_beads must be > 0, got {N_beads}")

    path   = os.path.expanduser(str(path))
    fmt    = _detect_format(path)
    info   = list_resolutions(path)
    chrsz  = get_chromosome_sizes(path)
    chrom  = _normalise_chrom(chrom, list(chrsz.keys()))

    chrom_len = chrsz[chrom]
    if region is None:
        start, end = 0, chrom_len
    else:
        start, end = int(region[0]), int(region[1])
        if start < 0 or end > chrom_len or start >= end:
            raise ValueError(
                f"Invalid region ({start}, {end}) for {chrom} "
                f"(length {chrom_len} bp)"
            )

    log.info(
        "  region: %s:%d-%d  (%.2f Mb)", chrom, start, end,
        (end - start) / 1e6,
    )

    resolution = choose_resolution(info["resolutions"], end - start, N_beads)

    # ── load ──────────────────────────────────────────────────────────────────
    if fmt == "hic":
        H = _load_from_hic(path, chrom, resolution, start, end, normalization)
    else:
        H = _load_from_cool(path, chrom, resolution, start, end, fmt)

    H = np.asarray(H, dtype=np.float64)

    # ── clean ─────────────────────────────────────────────────────────────────
    np.fill_diagonal(H, 0.0)
    H = handle_missing_bins(H, max_gap=max_gap)
    n_bins_raw = H.shape[0]

    # ── resize to N_beads × N_beads via average pooling ───────────────────────
    if resize and n_bins_raw != N_beads:
        H = pool_to_n_beads(H, N_beads)

    meta = {
        "chrom"      : chrom,
        "start"      : start,
        "end"        : end,
        "resolution" : resolution,
        "n_bins_raw" : n_bins_raw,
        "format"     : fmt,
    }

    # ── final shape guarantee ─────────────────────────────────────────────────
    if resize:
        assert H.shape == (N_beads, N_beads), (
            f"read_hic_matrix: final shape {H.shape} ≠ ({N_beads}, {N_beads})"
        )

    log.info(
        "read_hic_matrix done → shape=(%d, %d)  non-zero=%.1f%%  "
        "range=[%.3e, %.3e]",
        *H.shape,
        100.0 * (H > 0).sum() / H.size,
        float(H.min()), float(H.max()),
    )
    log.info("═" * 60)

    # always print final shape explicitly so the caller can verify at a glance
    print(f"[read_hic_matrix] Final matrix shape: {H.shape}")

    return H, meta


# ═════════════════════════════════════════════════════════════════════════════
# Layer 6 — Visualisation
# ═════════════════════════════════════════════════════════════════════════════

def plot_hic_matrix(
    H         : np.ndarray,
    meta      : Optional[dict] = None,
    log_scale : bool = True,
    show_oe   : bool = True,
    vmax_pct  : float = 98.0,
    figsize   : Optional[Tuple[float, float]] = None,
    save_path : Optional[str] = None,
) -> plt.Figure:
    """
    Visualise a Hi-C contact matrix.

    Two-panel figure:
      Left   : balanced counts (log₁₀ scale optional, red colourmap)
      Right  : log₂ O/E matrix (seismic diverging colourmap)

    Parameters
    ----------
    H         : (N, N) contact matrix from read_hic_matrix
    meta      : metadata dict (used for title and axis tick labels)
    log_scale : if True, display raw counts on log₁₀ scale
    show_oe   : if True, add O/E panel
    vmax_pct  : percentile for colour saturation (default 98)
    figsize   : figure size; auto-computed if None
    save_path : if given, save to this path (PNG/PDF/SVG)

    Returns
    -------
    fig : matplotlib Figure
    """
    if H.size == 0:
        raise ValueError("H is empty")

    n_panels = 2 if show_oe else 1
    if figsize is None:
        figsize = (6 * n_panels + 1, 5.5)

    fig, axes = plt.subplots(
        1, n_panels, figsize=figsize, constrained_layout=True
    )
    if n_panels == 1:
        axes = [axes]

    # ── title ─────────────────────────────────────────────────────────────────
    if meta is not None:
        title = (
            f"{meta.get('chrom', '')}:"
            f"{meta.get('start', 0) / 1e6:.1f}–"
            f"{meta.get('end', 0) / 1e6:.1f} Mb  "
            f"({meta.get('resolution', 0) / 1e3:.0f} kb res, "
            f"{H.shape[0]} bins)"
        )
    else:
        title = f"Hi-C contact matrix  ({H.shape[0]} bins)"

    # ── Panel 1: balanced / raw counts ────────────────────────────────────────
    ax   = axes[0]
    data = H.copy().astype(np.float64)
    data[data <= 0] = np.nan

    if log_scale:
        with np.errstate(divide="ignore", invalid="ignore"):
            disp = np.log10(data)
        disp[~np.isfinite(disp)] = np.nan
        cbar_label = "log₁₀(contact freq.)"
    else:
        disp       = data
        cbar_label = "contact frequency"

    finite_disp = disp[np.isfinite(disp)]
    vmax = np.percentile(finite_disp, vmax_pct) if finite_disp.size else 1.0
    vmin = np.percentile(finite_disp, 2.0)       if finite_disp.size else 0.0

    im0 = ax.imshow(
        disp, cmap="Reds", vmin=vmin, vmax=vmax,
        origin="upper", interpolation="nearest", aspect="equal",
    )
    plt.colorbar(im0, ax=ax, label=cbar_label, fraction=0.046, pad=0.04)
    ax.set_title("Balanced counts" + (" (log₁₀)" if log_scale else ""))
    _add_axis_labels(ax, meta, H.shape[0])

    # ── Panel 2: log₂ O/E ─────────────────────────────────────────────────────
    if show_oe:
        ax2  = axes[1]
        H_oe = _compute_oe_for_plot(H)

        finite_oe = H_oe[np.isfinite(H_oe)]
        abs_max   = (
            np.percentile(np.abs(finite_oe), vmax_pct)
            if finite_oe.size else 1.0
        )
        abs_max = max(abs_max, 1e-6)

        im1 = ax2.imshow(
            H_oe, cmap="seismic", vmin=-abs_max, vmax=abs_max,
            origin="upper", interpolation="nearest", aspect="equal",
        )
        plt.colorbar(im1, ax=ax2, label="log₂(O/E)", fraction=0.046, pad=0.04)
        ax2.set_title("log₂ O/E")
        _add_axis_labels(ax2, meta, H.shape[0])

    fig.suptitle(title, fontsize=10, y=1.01)

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        log.info("plot_hic_matrix: saved → %s", save_path)

    return fig


# ── helpers ───────────────────────────────────────────────────────────────────

def _compute_oe_for_plot(H: np.ndarray, pseudocount: float = 1e-6) -> np.ndarray:
    """Compute log₂ O/E for visualisation (diagonal-based expected)."""
    N  = H.shape[0]
    Hp = H.astype(np.float64) + pseudocount

    # expected[k] = mean of k-th diagonal
    expected = np.ones(N, dtype=np.float64)
    for k in range(1, N):
        diag = np.diagonal(Hp, offset=k)
        if diag.size:
            expected[k] = diag.mean()

    sep = np.abs(np.arange(N)[:, None] - np.arange(N)[None, :])
    E   = expected[sep]

    OE = np.log2(Hp / E)
    np.fill_diagonal(OE, np.nan)
    return OE


def _add_axis_labels(
    ax      : plt.Axes,
    meta    : Optional[dict],
    N       : int,
    n_ticks : int = 5,
) -> None:
    """Add genomic-coordinate tick labels to a Hi-C heatmap axis."""
    if meta is None:
        ax.set_xlabel("bin index")
        ax.set_ylabel("bin index")
        return

    start   = meta.get("start", 0)
    end     = meta.get("end",   N)
    pos     = np.linspace(0, N - 1, n_ticks, dtype=int)
    mb      = np.linspace(start / 1e6, end / 1e6, n_ticks)
    labels  = [f"{v:.1f}" for v in mb]

    ax.set_xticks(pos)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(pos)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("position (Mb)")
    ax.set_ylabel("position (Mb)")


# ═════════════════════════════════════════════════════════════════════════════
# Example — run with GM12878 ENCODE dataset
# ═════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    import sys

    # ─────────────────────────────────────────────────────────────────────────
    # Configuration
    # ─────────────────────────────────────────────────────────────────────────
    # Primary file (.hic format — Juicer / Aiden-lab)
    HIC_FILE  = os.path.expanduser(
        "~/Data/ENCODE/Hi-C/GM12878/ENCSR742SAT/ENCFF522ZQR.hic"
    )
    # The format (.hic / .cool / .mcool) is auto-detected from the extension.
    # Change HIC_FILE to any supported format — the rest of the pipeline
    # (resolution selection, pooling, cleaning, visualisation) is unchanged.

    CHROM    = "chr6"        # short chromosome → fast to load
    N_BEADS  = 1000            # desired output matrix size
    REGION   = None           # None = full chromosome
                              # e.g. (10_000_000, 40_000_000) for a sub-region

    # ─────────────────────────────────────────────────────────────────────────
    # Sanity check
    # ─────────────────────────────────────────────────────────────────────────
    if not os.path.isfile(HIC_FILE):
        log.error(
            "File not found: %s\n"
            "  Edit HIC_FILE at the bottom of read_hic.py to point to your data.",
            HIC_FILE,
        )
        sys.exit(1)

    detected_fmt = _detect_format(HIC_FILE)
    print(f"\n══ Hi-C file: {HIC_FILE}")
    print(f"   Format detected: {detected_fmt.upper()}")

    # ─────────────────────────────────────────────────────────────────────────
    # Inspect available resolutions and chromosomes
    # ─────────────────────────────────────────────────────────────────────────
    print("\n── Available resolutions ──────────────────────────────────────")
    info = list_resolutions(HIC_FILE)
    for r in info["resolutions"]:
        print(f"  {r:>12,} bp   ({r / 1000:.0f} kb)")

    print("\n── Chromosome sizes ───────────────────────────────────────────")
    chrsz = get_chromosome_sizes(HIC_FILE)
    for name, size in sorted(chrsz.items(), key=lambda kv: (len(kv[0]), kv[0])):
        print(f"  {name:<8}  {size:>12,} bp  ({size / 1e6:.1f} Mb)")

    # ─────────────────────────────────────────────────────────────────────────
    # Load the contact matrix
    # ─────────────────────────────────────────────────────────────────────────
    print(f"\n── Loading {CHROM}  (N_beads={N_BEADS}) ────────────────────")

    H, meta = read_hic_matrix(
        HIC_FILE,
        chrom         = CHROM,
        N_beads       = N_BEADS,
        region        = REGION,
        normalization = "KR",   # .hic only; cooler uses stored balance weights
        max_gap       = 3,      # interpolate gaps ≤ 3 consecutive empty bins
        resize        = True,   # pool/interpolate to exactly N_BEADS × N_BEADS
    )

    # ─────────────────────────────────────────────────────────────────────────
    # Matrix summary — final shape is ALWAYS printed here
    # ─────────────────────────────────────────────────────────────────────────
    print()
    print("┌─────────────────────────────────────────────────────────┐")
    print(f"│  Final matrix shape  :  {H.shape}                          ")
    print(f"│  Format              :  {meta['format'].upper()}                            ")
    print(f"│  Resolution used     :  {meta['resolution']:,} bp  "
          f"({meta['resolution'] / 1000:.0f} kb)")
    print(f"│  Native bin count    :  {meta['n_bins_raw']}                          ")
    print(f"│  Pooling direction   :  "
          + ("downsampling" if meta['n_bins_raw'] >= N_BEADS else "upsampling")
          + f"  ({meta['n_bins_raw']} → {N_BEADS})")
    print(f"│  Non-zero entries    :  {(H > 0).sum():,} / {H.size:,}  "
          f"({100 * (H > 0).mean():.1f}%)")
    print(f"│  Value range         :  {H.min():.3e} – {H.max():.3e}")
    print("└─────────────────────────────────────────────────────────┘")

    # explicit assertion — will raise if something went wrong
    assert H.shape == (N_BEADS, N_BEADS), (
        f"ERROR: expected ({N_BEADS}, {N_BEADS}), got {H.shape}"
    )
    print(f"\n✓  Shape assertion passed: H.shape == ({N_BEADS}, {N_BEADS})")

    # ─────────────────────────────────────────────────────────────────────────
    # Visualise
    # ─────────────────────────────────────────────────────────────────────────
    print("\n── Plotting ───────────────────────────────────────────────────")
    out_png = f"hic_{CHROM}_{N_BEADS}beads.png"
    fig = plot_hic_matrix(
        H,
        meta      = meta,
        log_scale = True,
        show_oe   = True,
        vmax_pct  = 98.0,
        save_path = out_png,
    )
    plt.show()

    # ─────────────────────────────────────────────────────────────────────────
    # MultiMM integration snippet
    # ─────────────────────────────────────────────────────────────────────────
    print("\n── Ready for MultiMM ──────────────────────────────────────────")
    print("  from hic_force import build_hic_force")
    print(f"  # H.shape == {H.shape}  ← guaranteed N_beads × N_beads")
    print(f"  force = build_hic_force(H, N_beads={N_BEADS}, rc=6.0,")
    print("                          mode='svd', K=10)")
    print("  system.addForce(force)")

    print("""
── Supported formats (auto-detected from file extension) ───────
  .hic   → hicstraw  (KR / SCALE / VC / VC_SQRT / NONE norms)
  .cool  → cooler    (uses stored KR balance weights)
  .mcool → cooler    (picks closest resolution sub-group)

  To switch format, just change HIC_FILE above — resolution
  selection, average pooling, cleaning and visualisation all
  work identically regardless of format.
────────────────────────────────────────────────────────────────""")
