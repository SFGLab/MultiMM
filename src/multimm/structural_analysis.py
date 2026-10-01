"""
Standalone structural-analysis module for MultiMM.

This isolates all the "post-simulation polymer analysis" logic (radius of
gyration, end-to-end distance, shape descriptors, bond/angle statistics,
distance-vs-genomic-separation scaling, local compaction, etc.) from the
general-purpose plotting helpers in `plots.py`, so it can be read, tested
and extended on its own.

Public entry point: `analyze_structure(V, save_path, name="structure")`.
It writes:
  - a plain-text report (`<save_path>/analysis/<name>_report.txt`)
  - ONE consolidated, high-quality multi-panel overview figure
    (`<save_path>/analysis/plots/<name>_overview.png`)
instead of the many small, inconsistently-binned plots this used to
produce.
"""

import os
import logging

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial import ConvexHull, distance
from scipy.stats import gaussian_kde, gamma as _gamma_dist

logger = logging.getLogger(__name__)

# Boltzmann constant in MD-conventional kJ/(mol*K) — matches the OpenMM unit
# system (energy in kJ/mol when mass is in amu and velocity in nm/ps), so no
# unit conversion is needed when combining with OpenMM-extracted velocities
# and masses.
_KB_KJ_PER_MOL_K = 0.0083144626181532


def _smart_bin_count(data, lo=15, hi=40):
    """
    Choose a sane number of histogram bins via the Freedman-Diaconis rule,
    clamped to [lo, hi] so distributions stay readable instead of showing
    hundreds of noisy, chaotic-looking bars.
    """
    data = np.asarray(data)
    data = data[np.isfinite(data)]

    if data.size < 2:
        return lo

    iqr = np.subtract(*np.percentile(data, [75, 25]))

    if iqr <= 0:
        return lo

    bin_width = 2 * iqr / (data.size ** (1 / 3))

    if bin_width <= 0:
        return lo

    n_bins = int(np.ceil((data.max() - data.min()) / bin_width))

    return int(np.clip(n_bins, lo, hi))


def _hist_with_kde(ax, data, color, xlabel, title, bin_range=None):
    """Smoothed, smart-binned histogram + KDE overlay for one panel."""
    n_bins = _smart_bin_count(data)

    ax.hist(
        data, bins=n_bins, range=bin_range, density=True,
        color=color, alpha=0.55, edgecolor="white", linewidth=0.6,
    )

    try:
        kde = gaussian_kde(data)
        xs = np.linspace(np.min(data), np.max(data), 200)
        ax.plot(xs, kde(xs), color=color, linewidth=2.0)
    except Exception:
        pass

    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Density")


def analyze_structure(V, save_path, name="structure"):
    """
    Advanced structural analysis for polymer-like 3D structures.

    Outputs:
    - detailed text report
    - a single consolidated, high-quality multi-panel figure (PNG only)
      instead of many small chaotic plots
    """

    sns.set_style("whitegrid")
    plt.rcParams.update({"font.size": 11})

    palette = sns.color_palette("mako", 6)

    # ------------------------------------------------------------
    # safety
    # ------------------------------------------------------------
    V = np.asarray(V, dtype=np.float64)
    V = V[np.isfinite(V).all(axis=1)]

    N = len(V)

    base = os.path.join(save_path, "analysis")
    os.makedirs(base, exist_ok=True)

    # center
    R_cm = np.mean(V, axis=0)
    Vc = V - R_cm

    # radius of gyration
    Rg = np.sqrt(np.mean(np.sum(Vc**2, axis=1)))

    # end-to-end
    Ree = np.linalg.norm(V[-1] - V[0])

    # pairwise distances
    dmat = distance.cdist(V, V)
    mean_dist = np.mean(dmat)

    # convex hull
    try:
        hull = ConvexHull(V)
        volume = hull.volume
    except Exception:
        volume = np.nan

    density = N / volume if volume > 0 else np.nan

    # gyration tensor
    G = np.dot(Vc.T, Vc) / N
    eigvals = np.sort(np.linalg.eigvalsh(G))

    l1, l2, l3 = eigvals

    asphericity = l3 - 0.5 * (l1 + l2)
    acylindricity = l2 - l1

    # bond lengths
    bonds = np.linalg.norm(np.diff(V, axis=0), axis=1)

    # angles (stiffness)
    v1 = V[1:-1] - V[:-2]
    v2 = V[2:] - V[1:-1]

    cos_angles = np.sum(v1 * v2, axis=1) / (
        np.linalg.norm(v1, axis=1) * np.linalg.norm(v2, axis=1) + 1e-8
    )

    angles = np.arccos(np.clip(cos_angles, -1, 1))

    # distance vs genomic separation
    separations = []
    spatial_dists = []

    for s in range(1, min(500, N // 2)):
        idx = np.arange(N - s)
        d = np.linalg.norm(V[idx + s] - V[idx], axis=1)

        separations.append(s)
        spatial_dists.append(np.mean(d))

    separations = np.array(separations)
    spatial_dists = np.array(spatial_dists)

    # local compaction (sliding window Rg)
    window = max(10, N // 100)
    local_rg = []

    for i in range(N - window):
        chunk = V[i:i + window]
        cm = np.mean(chunk, axis=0)
        local_rg.append(np.sqrt(np.mean(np.sum((chunk - cm) ** 2, axis=1))))

    local_rg = np.array(local_rg)

    # REPORT
    report_path = os.path.join(base, f"{name}_report.txt")

    with open(report_path, "w") as f:

        f.write("===== STRUCTURE ANALYSIS =====\n\n")

        f.write(f"N beads: {N}\n\n")

        f.write("---- Global ----\n")
        f.write(f"Rg: {Rg:.4f}\n")
        f.write(f"Ree: {Ree:.4f}\n")
        f.write(f"Mean distance: {mean_dist:.4f}\n\n")

        f.write("---- Volume ----\n")
        f.write(f"Volume: {volume:.4f}\n")
        f.write(f"Density: {density:.6f}\n\n")

        f.write("---- Shape ----\n")
        f.write(f"Eigenvalues: {eigvals}\n")
        f.write(f"Asphericity: {asphericity:.6f}\n")
        f.write(f"Acylindricity: {acylindricity:.6f}\n\n")

        f.write("---- Local properties ----\n")
        f.write(f"Mean bond length: {np.mean(bonds):.4f}\n")
        f.write(f"Mean angle (rad): {np.mean(angles):.4f}\n\n")

        f.write("Interpretation:\n")
        f.write("Rg ~ size of polymer\n")
        f.write("Distance vs separation → scaling law\n")
        f.write("Angles → stiffness\n")
        f.write("Local Rg → domain compaction\n")

    # ------------------------------------------------------------
    # ONE consolidated, high-quality figure (2 x 3 panels, PNG only)
    # ------------------------------------------------------------
    os.makedirs(base + "/plots", exist_ok=True)

    r = np.linalg.norm(Vc, axis=1)
    angles_deg = np.degrees(angles)

    fig, axes = plt.subplots(2, 3, figsize=(16, 9))

    # (1) bond lengths
    _hist_with_kde(
        axes[0, 0], bonds, palette[0],
        xlabel="Bond length", title="Bond Length Distribution",
    )

    # (2) angles (degrees — easier to read than radians)
    _hist_with_kde(
        axes[0, 1], angles_deg, palette[1],
        xlabel="Angle (°)", title="Angle Distribution", bin_range=(0, 180),
    )

    # (3) radial distribution from center of mass
    _hist_with_kde(
        axes[0, 2], r, palette[2],
        xlabel="Distance from COM", title="Radial Distribution",
    )

    # (4) scaling law, log-log, with fitted exponent (polymer physics)
    ax = axes[1, 0]
    valid = (separations > 0) & (spatial_dists > 0)

    if valid.sum() > 2:
        log_s = np.log10(separations[valid])
        log_r = np.log10(spatial_dists[valid])
        slope, intercept = np.polyfit(log_s, log_r, 1)
        fit_line = 10 ** (intercept + slope * log_s)

        ax.loglog(separations[valid], spatial_dists[valid], "o", ms=3,
                  color=palette[3], alpha=0.6, label="Data")
        ax.loglog(separations[valid], fit_line, "-", color="black",
                  linewidth=1.8, label=f"Fit: R(s) ~ s^{slope:.2f}")
        ax.legend(fontsize=9, frameon=False)
    else:
        ax.loglog(separations, spatial_dists, "-", color=palette[3])

    ax.set_title("Distance-Separation Scaling (log-log)", fontsize=11, fontweight="bold")
    ax.set_xlabel("Genomic separation s (beads)")
    ax.set_ylabel("Mean spatial distance R(s)")

    # (5) local compaction (sliding-window Rg along the chain)
    ax = axes[1, 1]
    ax.plot(local_rg, color=palette[4], linewidth=1.2)
    ax.axhline(np.mean(local_rg), color="black", linestyle="--", linewidth=1,
               label=f"Mean = {np.mean(local_rg):.3f}")
    ax.fill_between(np.arange(len(local_rg)), local_rg, np.mean(local_rg),
                     color=palette[4], alpha=0.15)
    ax.legend(fontsize=9, frameon=False)
    ax.set_title("Local Compaction (Sliding-Window Rg)", fontsize=11, fontweight="bold")
    ax.set_xlabel("Bead index")
    ax.set_ylabel("Local Rg")

    # (6) summary scorecard panel
    ax = axes[1, 2]
    ax.axis("off")
    ax.set_title("Summary", fontsize=11, fontweight="bold", loc="left")

    summary_lines = [
        ("N beads", f"{N}"),
        ("Radius of gyration (Rg)", f"{Rg:.4f}"),
        ("End-to-end distance (Ree)", f"{Ree:.4f}"),
        ("Mean pairwise distance", f"{mean_dist:.4f}"),
        ("Convex hull volume", f"{volume:.4f}"),
        ("Density (N/volume)", f"{density:.6f}"),
        ("Asphericity", f"{asphericity:.6f}"),
        ("Acylindricity", f"{acylindricity:.6f}"),
        ("Mean bond length", f"{np.mean(bonds):.4f}"),
        ("Mean angle", f"{np.mean(angles_deg):.1f}°"),
    ]

    y = 0.95
    for label, value in summary_lines:
        ax.text(0.0, y, label, transform=ax.transAxes, fontsize=10, ha="left", va="top")
        ax.text(1.0, y, value, transform=ax.transAxes, fontsize=10, ha="right", va="top",
                fontweight="bold", color=palette[5])
        y -= 0.095

    fig.suptitle(f"Structural Analysis — {name}", fontsize=15, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(base + f"/plots/{name}_overview.png", dpi=250)
    plt.close(fig)

    # ------------------------------------------------------------
    return {
        "Rg": Rg,
        "Ree": Ree,
        "volume": volume,
        "density": density,
        "asphericity": asphericity,
        "acylindricity": acylindricity,
    }


def _hist_with_fit(ax, data, pdf_fn, color, xlabel, title, fit_label=None,
                    ref_pdf_fn=None, ref_label=None, ref_color="#898781"):
    """Smart-binned histogram + an overlaid *theoretical* PDF curve (rather
    than an empirical KDE — see :func:`_hist_with_kde` for that case).

    Used for the velocity/energy diagnostics, where the question is not
    "what does the empirical density look like" but "does the empirical
    density match the Maxwell-Boltzmann prediction".  An optional second
    reference curve (``ref_pdf_fn``) overlays the theoretical prediction at
    the simulation's *target* temperature, so a mismatch between the
    estimated and target curves is immediately visible.
    """
    n_bins = _smart_bin_count(data)
    ax.hist(
        data, bins=n_bins, density=True,
        color=color, alpha=0.55, edgecolor="white", linewidth=0.6,
        label="Data",
    )
    xs = np.linspace(np.min(data), np.max(data), 300)
    try:
        ax.plot(xs, pdf_fn(xs), color=color, linewidth=2.2,
                label=fit_label or "Boltzmann fit")
    except Exception:
        pass
    if ref_pdf_fn is not None:
        try:
            ax.plot(xs, ref_pdf_fn(xs), color=ref_color, linewidth=1.6, linestyle="--",
                    label=ref_label or "Target T")
        except Exception:
            pass
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Density")
    ax.legend(fontsize=8, frameon=False)


def analyze_dynamics(
    velocities,
    masses,
    save_path,
    name="dynamics",
    target_temperature=None,
    rmsf=None,
    log=None,
):
    """
    Diagnostic analysis of bead *velocities* / kinetic fluctuations — a
    companion to :func:`analyze_structure` that answers "how much are the
    beads actually moving", rather than "what does the structure look
    like".

    Three independent theoretical-distribution checks are performed against
    the Maxwell-Boltzmann prediction, using a mass-weighted variable
    (``p = sqrt(m) * v``) so a single pair of curves is valid even when
    bead masses are heterogeneous:

      1. Mass-weighted speed  ``s = sqrt(m) * |v|``  →  Maxwell speed pdf
         ``f(s) = sqrt(2/pi) s^2/a^3 exp(-s^2/2a^2)``,  a = sqrt(kB T).
      2. Mass-weighted velocity components (vx, vy, vz pooled)  →  Normal
         pdf  N(0, kB T)  (one Cartesian component of a Boltzmann gas).
      3. Per-bead kinetic energy ``KE = 0.5 m v^2``  →  Gamma(3/2, kB T)
         (the 3-DOF energy distribution; mass-invariant already, no
         weighting needed).

    ``kB*T`` is estimated directly from the data via equipartition
    (``kB*T_est = (2/3) * mean(KE)``) and drawn as the solid fit; when
    ``target_temperature`` is supplied, the theoretical curve at the
    simulation's *set-point* temperature is overlaid as a dashed reference,
    so a mismatch between actual and target thermal energy is immediately
    visible (e.g. insufficient equilibration, a thermostat that hasn't
    caught up, or a frozen subset of beads).

    ``rmsf``, when supplied (per-bead root-mean-square fluctuation in nm,
    accumulated across the MD trajectory), is plotted per bead index to
    give a *global* picture of which regions of the chain move the most
    over the whole run — complementary to the single-frame speed snapshot.

    Parameters
    ----------
    velocities : (N, 3) ndarray, nm/ps (OpenMM convention)
    masses     : (N,) ndarray, amu (OpenMM convention) — matched index-wise
                 to ``velocities``.  Zero-mass entries (virtual sites) are
                 excluded from the velocity/energy statistics.
    save_path  : str — same root as :func:`analyze_structure`
    name       : str — output file stem
    target_temperature : float, optional — set-point T in Kelvin
    rmsf       : (N,) ndarray, optional — per-bead RMSF in nm across frames
    log        : logging.Logger, optional
    """
    _log = log or logger

    sns.set_style("whitegrid")
    plt.rcParams.update({"font.size": 11})
    palette = sns.color_palette("mako", 6)

    velocities = np.asarray(velocities, dtype=np.float64)
    masses = np.asarray(masses, dtype=np.float64)

    if velocities.ndim != 2 or velocities.shape[1] != 3 or len(masses) != len(velocities):
        _log.warning("analyze_dynamics: velocities/masses shape mismatch, skipping.")
        return {}

    finite = np.isfinite(velocities).all(axis=1) & np.isfinite(masses) & (masses > 0)
    if finite.sum() < 10:
        _log.warning("analyze_dynamics: not enough valid (mass>0, finite) beads, skipping.")
        return {}

    v = velocities[finite]
    m = masses[finite]
    N = len(v)

    speed = np.linalg.norm(v, axis=1)                       # (N,) nm/ps
    ke = 0.5 * m * speed ** 2                                # (N,) kJ/mol (per-bead KE)

    sqrt_m = np.sqrt(m)
    speed_scaled = sqrt_m * speed                            # ~ Maxwell(a=sqrt(kT))
    components_scaled = (v * sqrt_m[:, None]).ravel()        # ~ N(0, kT), 3N samples

    mean_ke = float(np.mean(ke))
    kT_est = (2.0 / 3.0) * mean_ke                           # equipartition: <KE> = 1.5 kT
    T_est = kT_est / _KB_KJ_PER_MOL_K if kT_est > 0 else float("nan")

    kT_target = None
    if target_temperature is not None and np.isfinite(target_temperature):
        kT_target = _KB_KJ_PER_MOL_K * float(target_temperature)

    def _maxwell_speed_pdf(a2):
        def f(s):
            return np.sqrt(2.0 / np.pi) * s ** 2 / (a2 ** 1.5) * np.exp(-(s ** 2) / (2.0 * a2))
        return f

    def _gaussian_pdf(var):
        def f(x):
            return 1.0 / np.sqrt(2.0 * np.pi * var) * np.exp(-(x ** 2) / (2.0 * var))
        return f

    def _ke_pdf(scale):
        return lambda x: _gamma_dist.pdf(x, a=1.5, scale=scale)

    base = os.path.join(save_path, "analysis")
    os.makedirs(base + "/plots", exist_ok=True)

    fig, axes = plt.subplots(2, 3, figsize=(16, 9))

    # (1) mass-weighted speed distribution vs Maxwell-Boltzmann
    _hist_with_fit(
        axes[0, 0], speed_scaled,
        pdf_fn=_maxwell_speed_pdf(kT_est) if kT_est > 0 else (lambda x: np.zeros_like(x)),
        ref_pdf_fn=_maxwell_speed_pdf(kT_target) if kT_target else None,
        color=palette[0],
        xlabel=r"$\sqrt{m}\,|v|$  (mass-weighted speed)",
        title="Speed Distribution vs Maxwell-Boltzmann",
        fit_label=f"Fit (T≈{T_est:.1f} K)" if np.isfinite(T_est) else "Fit",
        ref_label=f"Target T = {target_temperature:.0f} K" if target_temperature else None,
    )

    # (2) mass-weighted velocity components vs Gaussian (1-D Boltzmann)
    _hist_with_fit(
        axes[0, 1], components_scaled,
        pdf_fn=_gaussian_pdf(kT_est) if kT_est > 0 else (lambda x: np.zeros_like(x)),
        ref_pdf_fn=_gaussian_pdf(kT_target) if kT_target else None,
        color=palette[1],
        xlabel=r"$\sqrt{m}\,v_{x,y,z}$  (mass-weighted component)",
        title="Velocity Components vs Boltzmann (Gaussian)",
        fit_label=f"Fit (T≈{T_est:.1f} K)" if np.isfinite(T_est) else "Fit",
        ref_label=f"Target T = {target_temperature:.0f} K" if target_temperature else None,
    )

    # (3) per-bead kinetic energy vs the 3-DOF Boltzmann energy distribution
    _hist_with_fit(
        axes[0, 2], ke,
        pdf_fn=_ke_pdf(kT_est) if kT_est > 0 else (lambda x: np.zeros_like(x)),
        ref_pdf_fn=_ke_pdf(kT_target) if kT_target else None,
        color=palette[2],
        xlabel="Kinetic energy per bead (kJ/mol)",
        title="Kinetic Energy vs Boltzmann (Γ(3/2, kT))",
        fit_label=f"Fit (T≈{T_est:.1f} K)" if np.isfinite(T_est) else "Fit",
        ref_label=f"Target T = {target_temperature:.0f} K" if target_temperature else None,
    )

    # (4) per-bead fluctuation across the whole trajectory (RMSF) — "how
    # much does each part of the chain move, over the whole run"
    ax = axes[1, 0]
    has_rmsf = rmsf is not None and np.isfinite(np.asarray(rmsf)).any()
    if has_rmsf:
        rmsf = np.asarray(rmsf, dtype=np.float64)
        idx = np.arange(len(rmsf))
        ax.plot(idx, rmsf, color=palette[3], linewidth=1.2)
        mean_rmsf = float(np.nanmean(rmsf))
        ax.axhline(mean_rmsf, color="black", linestyle="--", linewidth=1,
                   label=f"Mean = {mean_rmsf:.3f} nm")
        ax.fill_between(idx, rmsf, mean_rmsf, color=palette[3], alpha=0.15)
        ax.legend(fontsize=9, frameon=False)
        ax.set_title("Per-Bead Fluctuation Across Trajectory (RMSF)", fontsize=11, fontweight="bold")
        ax.set_xlabel("Bead index")
        ax.set_ylabel("RMSF (nm)")
    else:
        ax.axis("off")
        ax.text(0.5, 0.5, "RMSF not available\n(no trajectory accumulated)",
                transform=ax.transAxes, ha="center", va="center", fontsize=10, color="#898781")

    # (5) instantaneous speed by bead index — single-frame snapshot of "who
    # is moving fast right now", complementary to the trajectory-wide RMSF
    ax = axes[1, 1]
    full_speed = np.full(len(velocities), np.nan)
    full_speed[finite] = speed
    idx_all = np.arange(len(velocities))
    ax.plot(idx_all, full_speed, color=palette[4], linewidth=0.8, alpha=0.85)
    ax.axhline(np.nanmean(full_speed), color="black", linestyle="--", linewidth=1,
               label=f"Mean = {np.nanmean(full_speed):.3f} nm/ps")
    ax.legend(fontsize=9, frameon=False)
    ax.set_title("Instantaneous Speed by Bead Index", fontsize=11, fontweight="bold")
    ax.set_xlabel("Bead index")
    ax.set_ylabel("Speed |v| (nm/ps)")

    # (6) summary scorecard
    ax = axes[1, 2]
    ax.axis("off")
    ax.set_title("Summary", fontsize=11, fontweight="bold", loc="left")

    summary_lines = [
        ("N beads (m>0)", f"{N}"),
        ("Estimated T (equipartition)", f"{T_est:.1f} K" if np.isfinite(T_est) else "—"),
    ]
    if target_temperature is not None:
        delta = T_est - target_temperature if np.isfinite(T_est) else float("nan")
        summary_lines.append(("Target T", f"{target_temperature:.1f} K"))
        summary_lines.append(("ΔT (est. − target)", f"{delta:+.1f} K" if np.isfinite(delta) else "—"))
    summary_lines += [
        ("Mean speed |v|", f"{np.mean(speed):.4f} nm/ps"),
        ("Mean kinetic energy", f"{mean_ke:.4f} kJ/mol"),
    ]
    if has_rmsf:
        summary_lines.append(("Mean RMSF", f"{np.nanmean(rmsf):.4f} nm"))
        summary_lines.append(("Max RMSF (most mobile bead)", f"{np.nanmax(rmsf):.4f} nm"))

    y = 0.92
    for label, value in summary_lines:
        ax.text(0.0, y, label, transform=ax.transAxes, fontsize=10, ha="left", va="top")
        ax.text(1.0, y, value, transform=ax.transAxes, fontsize=10, ha="right", va="top",
                fontweight="bold", color=palette[5])
        y -= 0.12

    fig.suptitle(f"Bead Dynamics Diagnostics — {name}", fontsize=15, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(base + f"/plots/{name}_dynamics.png", dpi=250)
    plt.close(fig)

    _log.info(
        f"Dynamics diagnostics: T_est≈{T_est:.1f} K"
        + (f" (target {target_temperature:.1f} K)" if target_temperature is not None else "")
        + f", mean speed={np.mean(speed):.4f} nm/ps, mean KE={mean_ke:.4f} kJ/mol"
    )

    return {
        "T_estimated": T_est,
        "target_temperature": float(target_temperature) if target_temperature is not None else None,
        "mean_speed": float(np.mean(speed)),
        "mean_kinetic_energy": mean_ke,
        "mean_rmsf": float(np.nanmean(rmsf)) if has_rmsf else None,
    }
