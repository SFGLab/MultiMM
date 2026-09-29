"""
quality_tests.py — post-simulation quality diagnostics for MultiMM
====================================================================

Each check is fast (pure-NumPy, no OpenMM calls needed), returns a
``QualityResult`` named-tuple, and emits a verdict + suggestion to the
logger.  ``run_quality_tests`` collects all results into a summary table
(printed via ``log_table``) and writes a CSV to
``<save_path>/metadata/quality_tests.csv``.

Checks implemented
------------------
1.  Energy stability          — detects divergence / NaN in the MD history
1b. MD mobility (RMSD)        — structural mobility during MD vs minimised structure
2.  Bond-distance distribution — adjacent-bead distances vs harmonic r0
3.  Angle distribution         — polymer bend angles vs equilibrium angle
4.  Excluded-volume overlaps   — counts bead pairs closer than EV threshold
5.  Compartment clustering     — A vs B 3-D spatial separation (if enabled)
5b. Compartment PC1            — |PC1| of model contact map vs input compartment signal
6.  Chromosome separation      — inter-chromosome centroid distance (if GW)
7.  Loop distance compliance   — loop anchor pairs within expected range (if loops)
7b. Loop signal enrichment     — loop anchors systematically closer than matched controls
8.  Container confinement      — fraction of beads inside spherical container
9.  B-lamina proximity         — B beads closer to shell than A beads
"""

from __future__ import annotations

import csv
import logging
import math
from typing import Any, Dict, List, NamedTuple, Optional, Sequence, Tuple

import numpy as np
from scipy.spatial import KDTree

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Public result type
# ---------------------------------------------------------------------------

class QualityResult(NamedTuple):
    name: str          # short test label
    status: str        # "PASS", "WARN", "FAIL", "SKIP"
    value: str         # human-readable measured value
    suggestion: str    # empty string when PASS


# ---------------------------------------------------------------------------
# 1. Energy stability
# ---------------------------------------------------------------------------

def check_energy_stability(md_history: Dict[str, List[float]]) -> QualityResult:
    """Inspect the MD potential-energy trace for divergence or NaN."""
    pot = md_history.get("potential", [])
    if len(pot) < 2:
        return QualityResult(
            "Energy stability", "SKIP",
            "no MD history", ""
        )

    arr = np.asarray(pot, dtype=float)
    has_nan = bool(np.any(~np.isfinite(arr)))
    start, end = arr[0], arr[-1]

    # Relative change over the trajectory
    rel_change = abs(end - start) / (abs(start) + 1e-12)
    # Divergence: energy grew more than 10× its initial magnitude
    diverged = bool(abs(end) > 10.0 * abs(start) + 1e4) or has_nan

    value = (
        f"E_start={start:.2e}  E_end={end:.2e}  "
        f"Δrel={rel_change:.1%}{'  NaN/Inf detected!' if has_nan else ''}"
    )

    if has_nan or diverged:
        return QualityResult(
            "Energy stability", "FAIL", value,
            "Energy diverged. Try: lower SIM_INTEGRATOR_STEP (e.g. 0.5 fs), "
            "lower SIM_TEMPERATURE, increase energy minimization tolerance, "
            "or reduce HIC_K_SCALE / EV_EPSILON."
        )
    if rel_change > 0.5:
        return QualityResult(
            "Energy stability", "WARN", value,
            "Energy changed >50% during MD. Consider longer equilibration or "
            "a smaller integrator step."
        )
    return QualityResult("Energy stability", "PASS", value, "")


# ---------------------------------------------------------------------------
# 2. Bond-distance distribution
# ---------------------------------------------------------------------------

def check_bond_distances(
    coords: np.ndarray,
    r0_nm: float,
    chr_ends: Optional[np.ndarray] = None,
    tol: float = 0.5,          # fraction of r0 tolerated as deviation
) -> QualityResult:
    """Check that adjacent-bead distances are close to the harmonic r0."""
    if chr_ends is None:
        chr_ends = np.array([0, len(coords)])

    dists: List[float] = []
    for s, e in zip(chr_ends[:-1], chr_ends[1:]):
        seg = coords[s:e]
        if len(seg) < 2:
            continue
        d = np.linalg.norm(np.diff(seg, axis=0), axis=1)
        dists.extend(d.tolist())

    if not dists:
        return QualityResult("Bond distances", "SKIP", "no bonds found", "")

    arr = np.asarray(dists)
    mean_d, std_d = float(arr.mean()), float(arr.std())
    frac_bad = float(np.mean(np.abs(arr - r0_nm) > tol * r0_nm))

    value = f"mean={mean_d:.4f} nm  std={std_d:.4f} nm  r0={r0_nm:.4f} nm  bad={frac_bad:.1%}"

    if frac_bad > 0.20:
        return QualityResult(
            "Bond distances", "FAIL", value,
            f"{frac_bad:.0%} of bonds deviate >{tol*100:.0f}% from r0. "
            "Try increasing POL_HARMONIC_BOND_K or running more minimization steps."
        )
    if frac_bad > 0.05:
        return QualityResult(
            "Bond distances", "WARN", value,
            f"{frac_bad:.0%} of bonds are far from r0 — consider a slightly stronger "
            "POL_HARMONIC_BOND_K or checking for force-field conflicts."
        )
    return QualityResult("Bond distances", "PASS", value, "")


# ---------------------------------------------------------------------------
# 3. Angle distribution
# ---------------------------------------------------------------------------

def check_angle_distribution(
    coords: np.ndarray,
    angle0_rad: float = math.pi,   # equilibrium angle (π = straight)
    chr_ends: Optional[np.ndarray] = None,
    warn_frac: float = 0.30,       # fraction of angles allowed to deviate >π/4
) -> QualityResult:
    """Check polymer stiffness: angles should cluster near the equilibrium."""
    if chr_ends is None:
        chr_ends = np.array([0, len(coords)])

    angles: List[float] = []
    for s, e in zip(chr_ends[:-1], chr_ends[1:]):
        seg = coords[s:e]
        if len(seg) < 3:
            continue
        v1 = seg[1:-1] - seg[:-2]
        v2 = seg[2:]   - seg[1:-1]
        n1 = np.linalg.norm(v1, axis=1, keepdims=True)
        n2 = np.linalg.norm(v2, axis=1, keepdims=True)
        mask = (n1.ravel() > 1e-12) & (n2.ravel() > 1e-12)
        cos_a = np.clip(
            np.sum(v1[mask] / n1[mask] * v2[mask] / n2[mask], axis=1),
            -1.0, 1.0
        )
        angles.extend(np.arccos(cos_a).tolist())

    if not angles:
        return QualityResult("Angle distribution", "SKIP", "no angles found", "")

    arr = np.asarray(angles)
    mean_a = float(arr.mean())
    frac_bent = float(np.mean(arr < angle0_rad - math.pi / 4))

    value = (
        f"mean={math.degrees(mean_a):.1f}°  "
        f"eq={math.degrees(angle0_rad):.1f}°  "
        f"bent(<135°)={frac_bent:.1%}"
    )

    if frac_bent > warn_frac:
        return QualityResult(
            "Angle distribution", "WARN", value,
            f"{frac_bent:.0%} of triplets are more bent than 135°. "
            "The polymer may be too flexible. Consider increasing "
            "POL_HARMONIC_ANGLE_CONSTANT_K."
        )
    return QualityResult("Angle distribution", "PASS", value, "")


# ---------------------------------------------------------------------------
# 4. Excluded-volume overlaps
# ---------------------------------------------------------------------------

def check_ev_overlaps(
    coords: np.ndarray,
    r0_nm: float,
    ev_fraction: float = 0.8,    # beads closer than ev_fraction * r0 are "overlapping"
    max_bad_frac: float = 0.02,
) -> QualityResult:
    """Count bead pairs closer than the EV hard-core threshold."""
    threshold = ev_fraction * r0_nm
    tree = KDTree(coords)
    pairs = tree.query_pairs(threshold, output_type="ndarray")
    # exclude bonded neighbours (|i-j| == 1) which are normally close
    if len(pairs):
        non_bonded = pairs[np.abs(pairs[:, 0] - pairs[:, 1]) > 1]
        n_bad = len(non_bonded)
    else:
        n_bad = 0

    n_beads = len(coords)
    n_pairs = n_beads * (n_beads - 1) // 2
    frac_bad = n_bad / max(n_pairs, 1)
    value = f"{n_bad} overlapping pairs  ({frac_bad:.3%} of all pairs)  threshold={threshold:.3f} nm"

    if frac_bad > max_bad_frac:
        return QualityResult(
            "EV overlaps", "WARN", value,
            f"{n_bad} bead pairs closer than {threshold:.3f} nm (EV threshold). "
            "Consider increasing EV_EPSILON or EV_POWER, or reducing competing "
            "attractive forces (HIC_K_SCALE, COB/loop strengths)."
        )
    return QualityResult("EV overlaps", "PASS", value, "")


# ---------------------------------------------------------------------------
# 5. Compartment clustering (A vs B spatial separation)
# ---------------------------------------------------------------------------

def check_compartment_clustering(
    coords: np.ndarray,
    compartments: np.ndarray,   # signed: positive = A, negative = B (or 1/-1/0)
    min_separation_ratio: float = 1.1,   # centroid_sep / mean_within must exceed this
) -> QualityResult:
    """Verify A and B compartments are spatially separated in 3D."""
    N = min(len(coords), len(compartments))
    C = compartments[:N]
    V = coords[:N]

    a_mask = C > 0
    b_mask = C < 0

    if a_mask.sum() < 5 or b_mask.sum() < 5:
        return QualityResult(
            "Compartment clustering", "SKIP",
            "too few A or B beads", ""
        )

    cA = V[a_mask].mean(axis=0)
    cB = V[b_mask].mean(axis=0)
    sep = float(np.linalg.norm(cA - cB))

    # within-compartment spread (mean distance to centroid)
    spread_A = float(np.linalg.norm(V[a_mask] - cA, axis=1).mean())
    spread_B = float(np.linalg.norm(V[b_mask] - cB, axis=1).mean())
    mean_spread = (spread_A + spread_B) / 2.0

    ratio = sep / (mean_spread + 1e-12)
    value = (
        f"A–B centroid sep={sep:.4f} nm  "
        f"mean within-spread={mean_spread:.4f} nm  ratio={ratio:.2f}"
    )

    if ratio < min_separation_ratio:
        return QualityResult(
            "Compartment clustering", "WARN", value,
            f"A and B compartments are not well separated (ratio={ratio:.2f} < {min_separation_ratio}). "
            "Try increasing COB_EA/COB_EB energy strengths, or check that the "
            "compartment BED file is correctly loaded."
        )
    return QualityResult("Compartment clustering", "PASS", value, "")


# ---------------------------------------------------------------------------
# 6. Chromosome separation
# ---------------------------------------------------------------------------

def check_chromosome_separation(
    coords: np.ndarray,
    chr_ends: np.ndarray,
    min_sep_nm: float = 0.05,   # minimum centroid–centroid distance expected
) -> QualityResult:
    """Verify that chromosomes occupy distinct spatial domains."""
    n_chrs = len(chr_ends) - 1
    if n_chrs < 2:
        return QualityResult(
            "Chromosome separation", "SKIP",
            "single chromosome or region", ""
        )

    centroids = np.array([
        coords[chr_ends[i]:chr_ends[i + 1]].mean(axis=0)
        for i in range(n_chrs)
    ])

    dists = []
    for i in range(n_chrs):
        for j in range(i + 1, n_chrs):
            dists.append(float(np.linalg.norm(centroids[i] - centroids[j])))

    arr = np.asarray(dists)
    min_d, mean_d = float(arr.min()), float(arr.mean())
    value = f"{n_chrs} chromosomes  min centroid sep={min_d:.4f} nm  mean={mean_d:.4f} nm"

    if min_d < min_sep_nm:
        return QualityResult(
            "Chromosome separation", "WARN", value,
            f"Some chromosomes have centroids <{min_sep_nm} nm apart. "
            "Consider enabling CHB_USE_CHROMOSOMAL_BLOCKS or increasing "
            "CHB_DE / CHB_KC."
        )
    return QualityResult("Chromosome separation", "PASS", value, "")


# ---------------------------------------------------------------------------
# 7. Loop distance compliance
# ---------------------------------------------------------------------------

def check_loop_distances(
    coords: np.ndarray,
    ms: np.ndarray,     # loop anchor indices (left)
    ns: np.ndarray,     # loop anchor indices (right)
    r0_loop_nm: float,  # LE_HARMONIC_BOND_R0 in nm
    tol: float = 3.0,   # allow up to tol * r0 (loops are soft)
) -> QualityResult:
    """Check that loop anchor pairs are within the expected distance."""
    if len(ms) == 0:
        return QualityResult("Loop distances", "SKIP", "no loops defined", "")

    ms_arr = np.asarray(ms, dtype=int)
    ns_arr = np.asarray(ns, dtype=int)

    valid = (ms_arr < len(coords)) & (ns_arr < len(coords))
    if valid.sum() == 0:
        return QualityResult("Loop distances", "SKIP", "loop indices out of range", "")

    d = np.linalg.norm(coords[ms_arr[valid]] - coords[ns_arr[valid]], axis=1)
    threshold = tol * r0_loop_nm
    frac_bad = float(np.mean(d > threshold))
    mean_d = float(d.mean())

    value = (
        f"{valid.sum()} loops  mean dist={mean_d:.4f} nm  "
        f"r0={r0_loop_nm:.4f} nm  >{tol:.0f}×r0: {frac_bad:.1%}"
    )

    if frac_bad > 0.30:
        return QualityResult(
            "Loop distances", "WARN", value,
            f"{frac_bad:.0%} of loops have anchors farther than {tol:.0f}×r0. "
            "Consider increasing LE_HARMONIC_BOND_K or reducing competing forces "
            "that stretch the polymer."
        )
    return QualityResult("Loop distances", "PASS", value, "")


# ---------------------------------------------------------------------------
# 8. Container confinement
# ---------------------------------------------------------------------------

def check_container_confinement(
    coords: np.ndarray,
    radius_nm: float,
    center: Optional[np.ndarray] = None,
) -> QualityResult:
    """Verify that beads are inside the spherical container."""
    if center is None:
        center = coords.mean(axis=0)

    r = np.linalg.norm(coords - center, axis=1)
    frac_outside = float(np.mean(r > radius_nm))
    max_r = float(r.max())
    value = f"R_container={radius_nm:.4f} nm  max_bead_r={max_r:.4f} nm  outside={frac_outside:.1%}"

    if frac_outside > 0.05:
        return QualityResult(
            "Container confinement", "WARN", value,
            f"{frac_outside:.1%} of beads are outside the container radius. "
            "Consider increasing SC_SCALE or reducing other forces that push "
            "beads outward."
        )
    return QualityResult("Container confinement", "PASS", value, "")


# ---------------------------------------------------------------------------
# 9. B-lamina proximity
# ---------------------------------------------------------------------------

def check_blamina_proximity(
    coords: np.ndarray,
    compartments: np.ndarray,
    nucleus_radius_nm: float,
    center: Optional[np.ndarray] = None,
    shell_frac: float = 0.75,   # B should be in outer shell_frac of radius
) -> QualityResult:
    """B-compartment beads should be closer to the nuclear envelope than A."""
    N = min(len(coords), len(compartments))
    C = compartments[:N]
    V = coords[:N]

    if center is None:
        center = V.mean(axis=0)

    r = np.linalg.norm(V - center, axis=1)
    a_mask = C > 0
    b_mask = C < 0

    if a_mask.sum() < 5 or b_mask.sum() < 5:
        return QualityResult(
            "B-lamina proximity", "SKIP",
            "too few A or B beads", ""
        )

    mean_r_A = float(r[a_mask].mean())
    mean_r_B = float(r[b_mask].mean())

    # B should be farther from center (closer to lamina)
    b_near_lamina_frac = float(np.mean(r[b_mask] > shell_frac * nucleus_radius_nm))
    a_near_lamina_frac = float(np.mean(r[a_mask] > shell_frac * nucleus_radius_nm))
    value = (
        f"mean r_A={mean_r_A:.4f} nm  mean r_B={mean_r_B:.4f} nm  "
        f"B in outer shell={b_near_lamina_frac:.1%}  A in outer shell={a_near_lamina_frac:.1%}"
    )

    if mean_r_B <= mean_r_A:
        return QualityResult(
            "B-lamina proximity", "WARN", value,
            "B compartment is NOT closer to the lamina than A. "
            "Consider enabling IBL_USE_B_LAMINA_INTERACTION or increasing "
            "IBL_SCALE."
        )
    if b_near_lamina_frac < 0.30:
        return QualityResult(
            "B-lamina proximity", "WARN", value,
            f"Only {b_near_lamina_frac:.0%} of B beads are in the outer "
            f"{(1-shell_frac)*100:.0f}% shell. Increase IBL_SCALE for stronger lamina attraction."
        )
    return QualityResult("B-lamina proximity", "PASS", value, "")


# ---------------------------------------------------------------------------
# 10. MD structural mobility (RMSD vs minimised)
# ---------------------------------------------------------------------------

def check_md_mobility(
    md_history: Dict[str, List[float]],
    frozen_threshold_nm: float = 0.005,   # mean RMSD below this → effectively frozen
    low_threshold_nm:    float = 0.03,    # mean RMSD below this → suspiciously low
) -> QualityResult:
    """Verify that the structure actually moves during MD (not frozen).

    Reads the per-sampling-step COM-removed RMSD relative to the energy-minimised
    structure that is accumulated in ``md_history["rmsd"]``.  A frozen structure
    (mean RMSD near zero) suggests the integrator step or temperature is too low,
    or that competing forces are overwhelming thermal fluctuations.

    Parameters
    ----------
    md_history         : dict populated by ``MultiMM.run_md``
    frozen_threshold_nm: RMSD below this → FAIL  (structure is frozen)
    low_threshold_nm   : RMSD below this → WARN  (very limited mobility)
    """
    rmsd_trace = md_history.get("rmsd", [])

    if len(rmsd_trace) < 2:
        return QualityResult(
            "MD mobility (RMSD)", "SKIP",
            "MD not run or RMSD not tracked", ""
        )

    arr      = np.asarray(rmsd_trace, dtype=float)
    mean_r   = float(arr.mean())
    max_r    = float(arr.max())
    final_r  = float(arr[-1])

    value = (
        f"mean RMSD={mean_r:.4f} nm  max={max_r:.4f} nm  "
        f"final={final_r:.4f} nm  over {len(arr)} steps"
    )

    if mean_r < frozen_threshold_nm:
        return QualityResult(
            "MD mobility (RMSD)", "FAIL", value,
            f"Structure is effectively frozen (mean RMSD={mean_r:.4f} nm < "
            f"{frozen_threshold_nm} nm).  Try: increase SIM_TEMPERATURE (e.g. 500–1000 K), "
            "increase SIM_INTEGRATOR_STEP, or reduce HIC_K_SCALE / LE_HARMONIC_BOND_K."
        )
    if mean_r < low_threshold_nm:
        return QualityResult(
            "MD mobility (RMSD)", "WARN", value,
            f"Very low structural mobility (mean RMSD={mean_r:.4f} nm < "
            f"{low_threshold_nm} nm).  The MD run may not be sampling conformational "
            "space effectively.  Consider higher SIM_TEMPERATURE or longer SIM_N_STEPS."
        )
    return QualityResult("MD mobility (RMSD)", "PASS", value, "")


# ---------------------------------------------------------------------------
# 11. Loop signal enrichment  (loop anchors vs matched non-loop pairs)
# ---------------------------------------------------------------------------

def check_loop_signal_enrichment(
    coords: np.ndarray,
    ms: np.ndarray,
    ns: np.ndarray,
    ds: Optional[np.ndarray] = None,
    le_fixed_distances: bool = False,
    n_control: int = 5000,
    seed: int = 0,
) -> QualityResult:
    """Check that loop anchors are systematically closer than matched controls.

    Two modes:

    * ``le_fixed_distances=True`` and ``ds`` provided: correlate actual 3D loop
      distances with the expected LE target distances (``ds``).  A high Pearson /
      Spearman r means the force field is pulling loop pairs to the right
      equilibrium distances.

    * Otherwise: compare mean 3D distance at loop anchors against a matched
      control set sampled at similar genomic separations.  Properly restrained
      loop pairs should be *closer* than random same-separation pairs.
    """
    ms_arr = np.asarray(ms, dtype=int)
    ns_arr = np.asarray(ns, dtype=int)
    n_loops = len(ms_arr)
    if n_loops == 0:
        return QualityResult("Loop enrichment", "SKIP", "no loop data", "")

    loop_dists = np.linalg.norm(coords[ms_arr] - coords[ns_arr], axis=1)

    if le_fixed_distances and ds is not None and len(ds) == n_loops:
        # -- distance correlation mode --
        from scipy.stats import pearsonr, spearmanr
        ds_arr = np.asarray(ds, dtype=float)
        valid = np.isfinite(loop_dists) & np.isfinite(ds_arr) & (ds_arr > 0)
        if valid.sum() < 3:
            return QualityResult("Loop enrichment", "SKIP", "too few valid loop pairs", "")
        r_p, p_p = pearsonr(loop_dists[valid], ds_arr[valid])
        r_s, _ = spearmanr(loop_dists[valid], ds_arr[valid])
        mean_err = float(np.mean(np.abs(loop_dists[valid] - ds_arr[valid])))
        value = (
            f"Pearson r={r_p:.3f} (p={p_p:.2e})  Spearman ρ={r_s:.3f}  "
            f"mean |Δd|={mean_err:.4f} nm  n={int(valid.sum())} loops"
        )
        if r_p < 0.25:
            return QualityResult(
                "Loop enrichment", "WARN", value,
                "Weak correlation between model loop distances and expected LE distances. "
                "Consider increasing LE_HARMONIC_BOND_K or running more minimisation steps."
            )
        return QualityResult("Loop enrichment", "PASS", value, "")

    else:
        # -- enrichment mode: loop vs matched random pairs --
        N = len(coords)
        genomic_seps = np.abs(ms_arr - ns_arr)
        rng = np.random.default_rng(seed)
        ctrl_i = rng.integers(0, N, size=n_control)
        ctrl_sep = rng.integers(
            int(genomic_seps.min()), int(genomic_seps.max()) + 1, size=n_control
        )
        ctrl_j = np.clip(ctrl_i + ctrl_sep, 0, N - 1)
        ctrl_dists = np.linalg.norm(coords[ctrl_i] - coords[ctrl_j], axis=1)

        mean_loop = float(loop_dists.mean())
        mean_ctrl = float(ctrl_dists.mean())
        enrichment = (mean_ctrl - mean_loop) / (mean_ctrl + 1e-12)

        value = (
            f"mean loop dist={mean_loop:.4f} nm  "
            f"mean control dist={mean_ctrl:.4f} nm  "
            f"enrichment (closer)={enrichment*100:.1f}%  n={n_loops} loops"
        )
        if enrichment < 0:
            return QualityResult(
                "Loop enrichment", "FAIL", value,
                "Loop anchors are FARTHER apart than genomically-matched random pairs. "
                "Loop forces may not be effective. Check LOOPS_PATH and LE_HARMONIC_BOND_K."
            )
        if enrichment < 0.05:
            return QualityResult(
                "Loop enrichment", "WARN", value,
                f"Loop anchors only {enrichment*100:.1f}% closer than random pairs — "
                "marginal enrichment.  Consider increasing LE_HARMONIC_BOND_K."
            )
        return QualityResult("Loop enrichment", "PASS", value, "")


# ---------------------------------------------------------------------------
# 12. Compartment PC1 correlation with input signal
# ---------------------------------------------------------------------------

def check_compartment_pc1(
    coords: np.ndarray,
    compartments: np.ndarray,
    n_sample: int = 2000,
) -> QualityResult:
    """Correlate the leading PC of the model contact map with input compartment labels.

    Reuses ``_coords_to_inv_contact_f32`` and ``pc1_of_oe`` already implemented
    in ``validation.py``, so there is no duplication of matrix algebra.

    The sign of PC1 is arbitrary; we always compare absolute values.
    When compartment data are discrete (≤ 4 unique values, e.g. A/B labels),
    we measure sign-agreement of signed PC1 with signed compartment labels.
    For continuous compartment signals we report Spearman ρ(|PC1|, |data|).
    """
    N = len(coords)
    n_comp = min(N, len(compartments))
    if n_comp < 10:
        return QualityResult("Compartment PC1", "SKIP", "not enough beads", "")

    # Subsample for speed on large structures
    if n_comp > n_sample:
        rng = np.random.default_rng(42)
        idx = rng.choice(n_comp, size=n_sample, replace=False)
        idx.sort()
        sub_coords = coords[idx]
        sub_comp   = compartments[idx]
    else:
        sub_coords = coords[:n_comp]
        sub_comp   = compartments[:n_comp]

    # Build inverse-contact matrix via fast Gram-matrix trick (from validation.py)
    try:
        from .validation import _coords_to_inv_contact_f32, pc1_of_oe
        inv_c = _coords_to_inv_contact_f32(sub_coords).astype(np.float64)
        pc1_model = pc1_of_oe(inv_c)                  # sign-arbitrary, shape (n,)
    except Exception as _e:
        return QualityResult("Compartment PC1", "SKIP", f"matrix build failed: {_e}", "")

    comp_signal = sub_comp.astype(float)
    n_unique = len(np.unique(comp_signal[np.isfinite(comp_signal)]))
    is_discrete = n_unique <= 4

    if is_discrete:
        # Sign-agreement (pick orientation that best matches data)
        sign_model = np.sign(pc1_model)
        sign_comp  = np.sign(comp_signal)
        agr_fwd = float(np.mean(sign_model == sign_comp))
        agr_rev = float(np.mean(-sign_model == sign_comp))
        agreement = max(agr_fwd, agr_rev)
        value = (
            f"discrete ({n_unique} labels)  "
            f"sign agreement={agreement:.1%}  n={len(comp_signal)}"
        )
        if agreement < 0.55:
            return QualityResult(
                "Compartment PC1", "WARN", value,
                f"PC1 sign agreement ({agreement:.0%}) barely above chance. "
                "Compartment forces may be insufficient — increase COB_EA/COB_EB."
            )
    else:
        # Continuous: Spearman ρ between |PC1| and |data|
        from scipy.stats import spearmanr
        abs_model = np.abs(pc1_model)
        abs_comp  = np.abs(comp_signal)
        valid = np.isfinite(abs_comp) & np.isfinite(abs_model)
        if valid.sum() < 5:
            return QualityResult("Compartment PC1", "SKIP", "too few finite values", "")
        rho, p = spearmanr(abs_model[valid], abs_comp[valid])
        value = (
            f"continuous  Spearman ρ(|PC1|, |data|)={rho:.3f}  p={p:.2e}  n={int(valid.sum())}"
        )
        if rho < 0.2:
            return QualityResult(
                "Compartment PC1", "WARN", value,
                f"Weak compartment PC1 correlation (ρ={rho:.2f}). "
                "Increase COB_EA/COB_EB or verify COMPARTMENT_PATH is correctly loaded."
            )

    return QualityResult("Compartment PC1", "PASS", value, "")


# ---------------------------------------------------------------------------
# Master runner
# ---------------------------------------------------------------------------

def run_quality_tests(
    coords: np.ndarray,
    args: Any,
    md_history: Optional[Dict[str, List[float]]] = None,
    compartments: Optional[np.ndarray] = None,
    chr_ends: Optional[np.ndarray] = None,
    ms: Optional[np.ndarray] = None,
    ns: Optional[np.ndarray] = None,
    ds: Optional[np.ndarray] = None,
    nucleus_radius_nm: Optional[float] = None,
    save_path: str = "results/",
) -> List[QualityResult]:
    """
    Run all applicable quality tests and report results.

    Parameters
    ----------
    coords          : (N, 3) array of bead positions in nm
    args            : SimulationConfig object (or any object with the same fields)
    md_history      : dict with keys 'step', 'potential', 'kinetic', 'total', 'temperature'
    compartments    : (N,) signed array (positive=A, negative=B) or None
    chr_ends        : chromosome boundary indices array, or None
    ms, ns          : loop anchor index arrays, or None
    ds              : loop equilibrium distances (nm) per pair, or None; used by
                      check_loop_signal_enrichment when LE_FIXED_DISTANCES=True
    nucleus_radius_nm : nuclear envelope radius in nm (used for container / lamina checks)
    save_path       : directory to write the CSV report
    """
    from .logger import log_table   # imported here to avoid circular imports at module level

    results: List[QualityResult] = []

    # ── helper: extract nm value from OpenMM Quantity or plain float ─────────
    def _nm(q: Any) -> float:
        try:
            from openmm.unit import nanometers
            if hasattr(q, "value_in_unit"):
                return float(q.value_in_unit(nanometers))
        except Exception:
            pass
        return float(q)

    # =========================================================================
    # BASIC TESTS — always run; these must always appear in the summary table.
    # Each check is wrapped individually so one crash cannot suppress the others.
    # =========================================================================

    # ── 1. Energy stability ───────────────────────────────────────────────────
    try:
        if md_history and args.SIM_RUN_MD:
            results.append(check_energy_stability(md_history))
        else:
            results.append(QualityResult("Energy stability", "SKIP", "MD not run", ""))
    except Exception as _e:
        results.append(QualityResult("Energy stability", "SKIP", f"check failed: {_e}", ""))

    # ── 1b. MD mobility (RMSD vs minimised) ──────────────────────────────────
    try:
        if md_history and args.SIM_RUN_MD:
            results.append(check_md_mobility(md_history))
        else:
            results.append(QualityResult("MD mobility (RMSD)", "SKIP", "MD not run", ""))
    except Exception as _e:
        results.append(QualityResult("MD mobility (RMSD)", "SKIP", f"check failed: {_e}", ""))

    # ── 2. Bond distances ─────────────────────────────────────────────────────
    try:
        if args.POL_USE_HARMONIC_BOND:
            r0 = _nm(args.POL_HARMONIC_BOND_R0)
            results.append(check_bond_distances(coords, r0, chr_ends=chr_ends))
        else:
            results.append(QualityResult("Bond distances", "SKIP", "harmonic bond disabled", ""))
    except Exception as _e:
        results.append(QualityResult("Bond distances", "SKIP", f"check failed: {_e}", ""))

    # ── 3. Angle distribution ─────────────────────────────────────────────────
    try:
        if args.POL_USE_HARMONIC_ANGLE:
            try:
                from openmm.unit import radians as _rad_unit
                if hasattr(args.POL_HARMONIC_ANGLE_R0, "value_in_unit"):
                    angle0 = float(args.POL_HARMONIC_ANGLE_R0.value_in_unit(_rad_unit))
                else:
                    angle0 = float(args.POL_HARMONIC_ANGLE_R0)
            except Exception:
                angle0 = math.pi
            results.append(check_angle_distribution(coords, angle0, chr_ends=chr_ends))
        else:
            results.append(QualityResult("Angle distribution", "SKIP", "harmonic angle disabled", ""))
    except Exception as _e:
        results.append(QualityResult("Angle distribution", "SKIP", f"check failed: {_e}", ""))

    # ── 4. EV overlaps ───────────────────────────────────────────────────────
    try:
        if args.EV_USE_EXCLUDED_VOLUME:
            r0 = _nm(args.POL_HARMONIC_BOND_R0)
            results.append(check_ev_overlaps(coords, r0))
        else:
            results.append(QualityResult("EV overlaps", "SKIP", "EV disabled", ""))
    except Exception as _e:
        results.append(QualityResult("EV overlaps", "SKIP", f"check failed: {_e}", ""))

    # =========================================================================
    # OPTIONAL TESTS — run when the corresponding forces / data are enabled.
    # A single try/except covers the remaining group; if any fails the summary
    # still shows the four basic tests above.
    # =========================================================================
    try:
        # ── 5. Compartment clustering ─────────────────────────────────────────
        comp_enabled = (
            args.COB_USE_COMPARTMENT_BLOCKS
            or args.SCB_USE_SUBCOMPARTMENT_BLOCKS
            or args.IBL_USE_B_LAMINA_INTERACTION
        )
        if comp_enabled and compartments is not None and len(compartments) > 0:
            results.append(check_compartment_clustering(coords, compartments))
        else:
            results.append(QualityResult("Compartment clustering", "SKIP",
                                         "no compartment force or no compartment data", ""))

        # ── 5b. Compartment PC1 correlation ──────────────────────────────────
        if compartments is not None and len(compartments) > 0:
            results.append(check_compartment_pc1(coords, compartments))
        else:
            results.append(QualityResult("Compartment PC1", "SKIP",
                                         "no compartment data", ""))

        # ── 6. Chromosome separation ──────────────────────────────────────────
        chb_enabled = args.CHB_USE_CHROMOSOMAL_BLOCKS or args.CF_USE_CENTRAL_FORCE
        if chb_enabled and chr_ends is not None and len(chr_ends) > 2:
            results.append(check_chromosome_separation(coords, chr_ends))
        else:
            results.append(QualityResult("Chromosome separation", "SKIP",
                                         "single region or chromosome force disabled", ""))

        # ── 7. Loop distance compliance ───────────────────────────────────────
        loops_enabled = args.LE_USE_HARMONIC_BOND and ms is not None and ns is not None
        if loops_enabled:
            r0_loop = _nm(args.LE_HARMONIC_BOND_R0)
            results.append(check_loop_distances(coords, ms, ns, r0_loop))
        else:
            results.append(QualityResult("Loop distances", "SKIP",
                                         "loop force disabled or no loop data", ""))

        # ── 7b. Loop signal enrichment ────────────────────────────────────────
        if loops_enabled:
            le_fixed = getattr(args, "LE_FIXED_DISTANCES", False)
            results.append(
                check_loop_signal_enrichment(
                    coords, ms, ns,
                    ds=ds,
                    le_fixed_distances=le_fixed,
                )
            )
        else:
            results.append(QualityResult("Loop enrichment", "SKIP",
                                         "loop force disabled or no loop data", ""))

        # ── 8. Container confinement ──────────────────────────────────────────
        if args.SC_USE_SPHERICAL_CONTAINER and nucleus_radius_nm is not None:
            results.append(check_container_confinement(coords, nucleus_radius_nm))
        else:
            results.append(QualityResult("Container confinement", "SKIP",
                                         "container disabled or radius unknown", ""))

        # ── 9. B-lamina proximity ─────────────────────────────────────────────
        blamina_enabled = args.IBL_USE_B_LAMINA_INTERACTION
        if (
            blamina_enabled
            and compartments is not None
            and len(compartments) > 0
            and nucleus_radius_nm is not None
        ):
            results.append(check_blamina_proximity(coords, compartments, nucleus_radius_nm))
        else:
            results.append(QualityResult("B-lamina proximity", "SKIP",
                                         "B-lamina disabled or missing data", ""))
    except Exception as _opt_exc:
        log.warning("Optional quality checks raised an exception: %s", _opt_exc)

    # ── Summary table (logger) ────────────────────────────────────────────────
    STATUS_ICON = {"PASS": "✅", "WARN": "⚠", "FAIL": "✗", "SKIP": "–"}
    table_rows = [
        (f"{STATUS_ICON.get(r.status, '?')} {r.status:<4}  {r.name}", r.value)
        for r in results
    ]
    log_table(table_rows, title="Quality tests — summary", log_fn=log.info)

    # Emit suggestions for non-passing tests
    for r in results:
        if r.suggestion:
            log.warning("  [%s] %s → %s", r.status, r.name, r.suggestion)

    # ── CSV output ────────────────────────────────────────────────────────────
    import os
    csv_path = os.path.join(save_path, "metadata", "quality_tests.csv")
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)
    try:
        with open(csv_path, "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(["test", "status", "value", "suggestion"])
            for r in results:
                writer.writerow([r.name, r.status, r.value, r.suggestion])
        log.info("Quality tests CSV written to %s", csv_path)
    except OSError as exc:
        log.warning("Could not write quality_tests.csv: %s", exc)

    return results
