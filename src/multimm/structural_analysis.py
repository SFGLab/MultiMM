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
from scipy.stats import gaussian_kde

logger = logging.getLogger(__name__)


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
