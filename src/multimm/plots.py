import logging
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyvista as pv
import seaborn as sns
from matplotlib.pyplot import figure
from scipy.spatial import ConvexHull, distance
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from scipy.stats import gaussian_kde, rankdata
from matplotlib.lines import Line2D
import matplotlib.colors as mcolors
from mpl_toolkits.mplot3d import Axes3D
from .utils import get_coordinates_cif, oe_matrix
from .hic_force import get_boltzmann_p_func, auto_contact_scale

logger = logging.getLogger(__name__)

# Single colormap used for every Hi-C-style heatmap in this module, paired
# with _hic_display_transform below, so all heatmaps share one normalised,
# consistently-coloured display pipeline.
#
# Custom blue <-> red diverging map with a visibly-inked warm-gray midpoint,
# replacing matplotlib's built-in "coolwarm". coolwarm fades to near-white
# across a broad band around its center — fine for a scalar field where most
# values sit near an extreme, but wrong here: in a Hi-C O/E matrix MOST pairs
# sit near "expected" (log2 ratio ~ 0), i.e. exactly in that broad pale band,
# so the whole heatmap read as washed out/whitish. This map keeps the
# near-zero region a visible muted gray instead of letting it fade into the
# plot background, while still ramping to fully saturated blue/red at the
# ±1 (clipped log2 O/E) extremes.
_HIC_CMAP = mcolors.LinearSegmentedColormap.from_list(
    "multimm_hic_div",
    [
        (0.00, "#0d366b"),   # darkest blue  — strongly depleted
        (0.22, "#2a78d6"),   # blue
        (0.50, "#d7d4cc"),   # muted warm gray — visible "nothing", not white
        (0.78, "#e34948"),   # red
        (1.00, "#7a1414"),   # darkest red   — strongly enriched
    ],
    N=256,
)


pv.set_jupyter_backend("server")
color_dict = {-2: "#bf0020", -1: "#e36a24", 1: "#20c8e6", 2: "#181385", 0: "#ffffff"}
comp_dict = {-2: "B2", -1: "B1", 1: "A2", 2: "A1", 0: "no compartment"}

def _render_chain_image(points, colors=None, r=0.2, cmap="coolwarm", zoom=1.0):
    """Render a polymer chain via `viz_structure()` to a temp screenshot for
    embedding in a matplotlib panel. `zoom` > 1.0 tightens the camera for a
    small subplot. Returns None (never raises) if offscreen rendering fails,
    so callers can fall back to a plain scatter."""
    import tempfile

    try:
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            tmp_path = tmp.name
        try:
            viz_structure(points, colors=colors, r=r, cmap=cmap, save_path=tmp_path, zoom=zoom)
            img = plt.imread(tmp_path)
        finally:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)
        return img
    except Exception as exc:
        logger.warning(f"pyvista chain rendering unavailable, falling back to scatter: {exc}")
        return None


def plot_projection(struct_3D, Cs, save_path, name="structure"):
    """Chromatin structural analysis centered on COM, in one 2x3 figure:
    (1) 3D structure, (2) PCA projection, (3) radial distribution from COM,
    (4) density landscape in PCA space, (5) radial distance by
    subcompartment or PCA explained variance, (6) summary scorecard.
    Subcompartment state is signed (B-type negative, A-type positive) and
    shown with one diverging blue<->red colormap throughout.
    """

    sns.set_style("whitegrid")
    plt.rcParams.update({"font.size": 11})

    diverging_cmap = plt.get_cmap("RdBu_r")
    sequential_cmap = "Blues"
    accent = "#2a78d6"  # single-series accent, consistent with the sequential hue

    # preprocessing (STRICT ALIGNMENT GUARANTEE)
    X = np.asarray(struct_3D, dtype=np.float64)
    has_comps = Cs is not None

    if has_comps:
        Cs = np.asarray(Cs)

        N = min(len(X), len(Cs))
        X = X[:N]
        Cs = Cs[:N]

        mask = np.isfinite(X).all(axis=1)

        X = X[mask]
        Cs = Cs[mask]

        # remove undefined compartments
        valid = Cs != 0
        X = X[valid]
        Cs = Cs[valid]

        has_comps = len(X) > 0
    else:
        mask = np.isfinite(X).all(axis=1)
        X = X[mask]

    if len(X) < 4:
        logger.warning(f"plot_projection: not enough valid points for '{name}', skipping.")
        return

    # CENTER OF MASS SHIFT
    com = X.mean(axis=0)
    Xc = X - com

    # PCA (COM-centered)
    pca = PCA(n_components=min(3, Xc.shape[1]))
    X_pca = pca.fit_transform(Xc)
    r = np.linalg.norm(Xc, axis=1)
    X0 = Xc - Xc.mean(axis=0)
    G = (X0.T @ X0) / len(X0)
    eigvals = np.linalg.eigvalsh(G)
    anisotropy_scalar = np.sqrt(eigvals.max() / (eigvals.min() + 1e-12))
    explained_var = pca.explained_variance_ratio_

    df = pd.DataFrame({
        "x": Xc[:, 0], "y": Xc[:, 1], "z": Xc[:, 2],
        "pc1": X_pca[:, 0], "pc2": X_pca[:, 1],
        "r_com": r,
    })

    if has_comps:
        df["subcomp"] = Cs
        unique_sub = np.sort(df.subcomp.unique())
        abs_max = np.max(np.abs(unique_sub)) if len(unique_sub) > 0 else 1.0
        sub_norm = mcolors.Normalize(vmin=-abs_max, vmax=abs_max)
        point_colors = diverging_cmap(sub_norm(df.subcomp.values))
    else:
        r_norm = mcolors.Normalize(vmin=r.min(), vmax=r.max())
        point_colors = plt.get_cmap(sequential_cmap)(r_norm(r))

    base = os.path.join(save_path, "plots")
    os.makedirs(base, exist_ok=True)

    fig = plt.figure(figsize=(17, 10))
    gs = fig.add_gridspec(2, 3)

    # ---- (1) 3D structure, COM-centered — rendered as a CONNECTED chain
    # (adjacent beads joined by a tube) via the existing viz_structure()/
    # polyline_from_points() pyvista pipeline, then embedded as an image;
    # falls back to a plain scatter if the off-screen renderer is
    # unavailable for any reason.
    ax1 = fig.add_subplot(gs[0, 0])
    chain_img = _render_chain_image(
        Xc,
        colors=(df.subcomp.values if has_comps else None),
        r=0.2,
        cmap="coolwarm",
        zoom=1.5,  # tighten in beyond PyVista's default auto-fit framing so the
                   # chain is clearly visible rather than small in a sea of margin
    )
    if chain_img is not None:
        ax1.imshow(chain_img)
        ax1.axis("off")
    else:
        ax1.remove()
        ax1 = fig.add_subplot(gs[0, 0], projection="3d")
        ax1.plot(Xc[:, 0], Xc[:, 1], Xc[:, 2], color="#c3c2b7", linewidth=0.6, alpha=0.6, zorder=1)
        ax1.scatter(Xc[:, 0], Xc[:, 1], Xc[:, 2], c=point_colors, s=4, alpha=0.8, linewidths=0, zorder=2)
        ax1.set_xlabel("x (nm)")
        ax1.set_ylabel("y (nm)")
        ax1.set_zlabel("z (nm)")
    ax1.set_title("3D Structure (COM-centered, chain)", fontsize=11, fontweight="bold")

    # ---- (2) PCA projection ----
    ax2 = fig.add_subplot(gs[0, 1])
    sc = ax2.scatter(df.pc1, df.pc2, c=point_colors, s=8, alpha=0.75, linewidths=0)
    ax2.set_title("PCA Projection (PC1 vs PC2)", fontsize=11, fontweight="bold")
    ax2.set_xlabel(f"PC1 ({explained_var[0]*100:.0f}% var)")
    ax2.set_ylabel(f"PC2 ({explained_var[1]*100:.0f}% var)" if len(explained_var) > 1 else "PC2")

    if has_comps:
        cbar = fig.colorbar(plt.cm.ScalarMappable(norm=sub_norm, cmap=diverging_cmap), ax=ax2, fraction=0.046, pad=0.04)
        cbar.set_label("Subcompartment (B ← 0 → A)")
    else:
        cbar = fig.colorbar(plt.cm.ScalarMappable(norm=r_norm, cmap=sequential_cmap), ax=ax2, fraction=0.046, pad=0.04)
        cbar.set_label("Distance from COM")

    # ---- (3) radial distribution from COM ----
    ax3 = fig.add_subplot(gs[0, 2])
    if has_comps:
        for scv in unique_sub:
            if scv == 0:
                continue
            sub = df[df.subcomp == scv]
            if len(sub) < 5:
                continue
            color = diverging_cmap(sub_norm(scv))
            sns.kdeplot(sub.r_com, ax=ax3, fill=True, alpha=0.25, color=color, linewidth=1.8,
                        label=comp_dict.get(scv, str(scv)))
        ax3.legend(fontsize=8, frameon=False, title="Subcompartment")
    else:
        _hist_with_kde(ax3, r, accent, xlabel="Distance from COM", title="")
    ax3.set_title("Radial Distribution from COM", fontsize=11, fontweight="bold")
    ax3.set_xlabel("Distance from COM")
    ax3.set_ylabel("Density")

    # ---- (4) free-energy landscape in PCA space ----
    # F(PC1, PC2) = -ln[ P(PC1, PC2) / P_max ]  — the conformational analogue
    # of a potential-of-mean-force surface (in units of kT, the standard MD
    # convention): the most *populated* region of conformational space sits
    # at F=0, and less-visited regions rise above it. A colorbar makes the
    # unit and direction explicit instead of leaving a bare, unlabeled
    # density plot that could be mistaken for a probability or a generic
    # heatmap.
    ax4 = fig.add_subplot(gs[1, 0])
    x = df.pc1.values
    y = df.pc2.values
    Xg, Yg = np.mgrid[x.min():x.max():150j, y.min():y.max():150j]
    pos = np.vstack([Xg.ravel(), Yg.ravel()])

    try:
        kde_all = gaussian_kde(np.vstack([x, y]))
        Z_all = kde_all(pos).reshape(Xg.shape)
        Z_all = np.clip(Z_all, 1e-300, None)
        F = -np.log(Z_all / Z_all.max())
        # Clip the free-energy scale at a generous but finite ceiling so a
        # handful of near-empty bins (F -> inf) don't wash out the contrast
        # across the populated basin.
        f_ceiling = np.percentile(F, 99.0)
        F = np.clip(F, 0.0, f_ceiling)

        cf = ax4.contourf(Xg, Yg, F, levels=20, cmap="viridis")
        ax4.contour(Xg, Yg, F, levels=8, colors="white", linewidths=0.4, alpha=0.4)
        cbar4 = fig.colorbar(cf, ax=ax4, fraction=0.046, pad=0.04)
        cbar4.set_label("Relative free energy  F = -ln(P/P$_{max}$)  [kT]")
    except Exception:
        ax4.text(0.5, 0.5, "Free-energy landscape unavailable\n(insufficient spread)",
                  transform=ax4.transAxes, ha="center", va="center", fontsize=9, color="#898781")

    if has_comps:
        # Overlay per-subcompartment occupancy contours (outline only, no
        # fill) on top of the shared free-energy background, so this panel
        # answers both "where is the stable basin" (background) and "which
        # subcompartment lives there" (overlay) at once.
        for scv in unique_sub:
            if scv == 0:
                continue
            sub = df[df.subcomp == scv]
            if len(sub) < 10:
                continue
            try:
                kde = gaussian_kde([sub.pc1, sub.pc2])
                Zc = kde(pos).reshape(Xg.shape)
            except Exception:
                continue
            color = diverging_cmap(sub_norm(scv))
            ax4.contour(Xg, Yg, Zc, levels=4, colors=[color], linewidths=1.4, alpha=0.95)

        legend_elements = [
            Line2D([0], [0], color=diverging_cmap(sub_norm(v)), lw=2,
                   label=comp_dict.get(v, str(v)))
            for v in unique_sub if v != 0
        ]
        ax4.legend(handles=legend_elements, frameon=False, fontsize=8, title="Subcompartment (contour overlay)")
        ax4.set_title("Free-Energy Landscape + Subcompartment Occupancy", fontsize=11, fontweight="bold")
    else:
        ax4.set_title("Free-Energy Landscape (PCA space)", fontsize=11, fontweight="bold")
    ax4.set_xlabel("PC1")
    ax4.set_ylabel("PC2")

    # ---- (5) radial by subcompartment, or PCA explained variance ----
    ax5 = fig.add_subplot(gs[1, 1])
    if has_comps:
        sns.violinplot(
            data=df, x="subcomp", y="r_com",
            palette=[diverging_cmap(sub_norm(v)) for v in unique_sub],
            inner=None, cut=0, ax=ax5,
        )
        sns.stripplot(data=df, x="subcomp", y="r_com", color="black", alpha=0.2, size=1.3, ax=ax5)
        ax5.set_xticks(range(len(unique_sub)))
        ax5.set_xticklabels([comp_dict.get(v, str(v)) for v in unique_sub])
        ax5.set_title("Radial Distance by Subcompartment", fontsize=11, fontweight="bold")
        ax5.set_xlabel("Subcompartment state")
        ax5.set_ylabel("Distance from COM")
    else:
        n_comp = len(explained_var)
        ax5.bar(range(1, n_comp + 1), explained_var * 100, color=accent, alpha=0.85,
                edgecolor="white")
        ax5.set_xticks(range(1, n_comp + 1))
        ax5.set_xticklabels([f"PC{i}" for i in range(1, n_comp + 1)])
        ax5.set_title("PCA Explained Variance", fontsize=11, fontweight="bold")
        ax5.set_ylabel("Variance explained (%)")

    # ---- (6) summary scorecard ----
    ax6 = fig.add_subplot(gs[1, 2])
    ax6.axis("off")
    ax6.set_title("Summary", fontsize=11, fontweight="bold", loc="left")

    summary_lines = [
        ("N beads", f"{len(X)}"),
        ("Anisotropy (λmax/λmin)^0.5", f"{anisotropy_scalar:.3f}"),
        ("PC1 explained var.", f"{explained_var[0]*100:.1f}%"),
        ("PC2 explained var.", f"{explained_var[1]*100:.1f}%" if len(explained_var) > 1 else "—"),
        ("Mean dist. from COM", f"{np.mean(r):.4f}"),
    ]
    if has_comps:
        n_a = int(np.sum(df.subcomp > 0))
        n_b = int(np.sum(df.subcomp < 0))
        summary_lines.append(("A-type beads", f"{n_a}"))
        summary_lines.append(("B-type beads", f"{n_b}"))

    y = 0.92
    for label, value in summary_lines:
        ax6.text(0.0, y, label, transform=ax6.transAxes, fontsize=10, ha="left", va="top")
        ax6.text(1.0, y, value, transform=ax6.transAxes, fontsize=10, ha="right", va="top",
                 fontweight="bold", color=accent)
        y -= 0.13

    fig.suptitle(f"Structural Projection — {name}", fontsize=15, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(os.path.join(base, f"{name}_projection.png"), dpi=250)
    plt.close(fig)


def plot_compartment_aggregation(same_d, diff_d, purity, chance_purity, save_dir,
                                  name="compartment_aggregation", p_value=None, effect=None, k_eff=10):
    """One consolidated figure checking whether beads sharing a compartment
    sign (A-A / B-B — from EITHER an input .bed track or Hi-C-derived PC1)
    actually end up closer together in 3D than different-compartment (A-B)
    beads. This is distinct from plot_compartment_validation (which checks
    1D PC1-vs-input-track agreement) — here we measure real 3D spatial
    segregation. Pure plotting function: the pairwise-distance/purity
    arrays and statistics are computed by validation.validate_compartment_
    aggregation, matching the plot_loop_validation / plot_compartment_
    validation convention (stats in validation.py, drawing here).
    """
    os.makedirs(save_dir, exist_ok=True)
    sns.set_style("whitegrid")

    same_color, diff_color = "#2a78d6", "#c0392b"
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))

    ax1 = axes[0]
    parts = ax1.violinplot([same_d, diff_d], showmedians=True, widths=0.8)
    for pc, color in zip(parts["bodies"], [same_color, diff_color]):
        pc.set_facecolor(color)
        pc.set_alpha(0.55)
        pc.set_edgecolor(color)
    for key in ("cbars", "cmins", "cmaxes", "cmedians"):
        parts[key].set_color("#4a4a45")
    ax1.set_xticks([1, 2])
    ax1.set_xticklabels(["Same compartment\n(A-A / B-B)", "Different compartment\n(A-B)"])
    ax1.set_ylabel("Pairwise 3D distance")
    ax1.set_title("Spatial Separation by Compartment", fontsize=11, fontweight="bold")
    has_p = p_value is not None and np.isfinite(p_value)
    sig_label = "n.s." if not has_p else ("p < 0.001" if p_value < 1e-3 else f"p = {p_value:.3g}")
    eff_label = f"\nEffect size r = {effect:.3f}" if (effect is not None and np.isfinite(effect)) else ""
    ax1.text(0.5, 0.98, f"Mann-Whitney U (same < diff): {sig_label}{eff_label}",
             transform=ax1.transAxes, ha="center", va="top", fontsize=9,
             bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#c3c2b7", alpha=0.9))

    ax2 = axes[1]
    sns.histplot(purity, bins=np.linspace(0, 1, 21), color=same_color, alpha=0.75, ax=ax2, stat="density")
    ax2.axvline(chance_purity, color=diff_color, linestyle="--", linewidth=1.6,
                label=f"Chance level ({chance_purity:.2f})")
    ax2.axvline(float(np.mean(purity)), color=same_color, linestyle="-", linewidth=1.8,
                label=f"Observed mean ({np.mean(purity):.2f})")
    ax2.set_xlim(0, 1)
    ax2.set_xlabel(f"Fraction of {k_eff} nearest 3D neighbors\nsharing the bead's compartment sign")
    ax2.set_ylabel("Density")
    ax2.set_title("Local Compartment Purity", fontsize=11, fontweight="bold")
    ax2.legend(frameon=False, fontsize=9)

    fig.suptitle("Compartment Aggregation Validation (3D spatial segregation)", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])

    out_path = os.path.join(save_dir, f"{name}.png")
    fig.savefig(out_path, dpi=250, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Saved compartment aggregation plot → {out_path}")


def _save_plotter(plotter, save_path):
    """
    Save PyVista scene in multiple formats.
    """
    os.makedirs(os.path.dirname(save_path), exist_ok=True)

    # PyVista supports direct image export
    plotter.show(screenshot=save_path + ".png")

    # optional additional formats via export (vtk scene)
    plotter.export_vtkjs(save_path + ".vtkjs")


def polyline_from_points(points):
    poly = pv.PolyData()
    poly.points = points

    the_cell = np.arange(0, len(points), dtype=np.int_)
    the_cell = np.insert(the_cell, 0, len(points))
    poly.lines = the_cell

    return poly


def viz_structure(V, colors=None, r=0.1, cmap="coolwarm", save_path=None, zoom=1.0, legend_labels=None):
    """
    Visualize structure V and optionally save it to a file.

    `zoom` > 1.0 tightens the camera in on the structure beyond PyVista's
    default auto-fit framing (which leaves generous margins around the
    bounding box) — useful for embedding the render into a small subplot
    panel where the chain should visibly fill the frame rather than sit as
    a small shape surrounded by empty space. 1.0 (default) keeps the
    original auto-fit framing unchanged.

    `legend_labels`: optional (neg_label, pos_label) pair (e.g. ("B",
    "A")) — when given together with `save_path`, a small frameless
    legend is composited onto the saved screenshot so signed/compartment
    colouring is actually readable (viz_structure otherwise never shows a
    scalar bar/key).
    """

    logger.info(
        f"Visualizing structure: N={len(V)}, "
        f"colored={colors is not None}, "
        f"save_path={save_path}"
    )

    polyline = polyline_from_points(V)
    polyline["scalars"] = np.arange(polyline.n_points)

    if colors is not None and len(colors) > 0:

        colors = np.array(colors[: len(V)])

        logger.info("Color mapping enabled (signed scheme: neg/zero/pos)")

        # ------------------------------------------------------------
        # NEW: signed piecewise normalization
        # ------------------------------------------------------------

        color_values = np.zeros(len(colors), dtype=float)

        neg = colors < 0
        pos = colors > 0
        zero = colors == 0

        # normalize negatives -> [0, 1]. Guard against a degenerate span
        # (all negatives share one magnitude, e.g. binary +-1 Hi-C-derived
        # compartments): dividing by ~0 would silently give 0 for every
        # element, which is indistinguishable from the "no compartment"
        # (zero) case below — default to 0.0 instead (already the correct,
        # fully-saturated "B" end).
        if np.any(neg):
            nmin, nmax = colors[neg].min(), colors[neg].max()
            nspan = nmax - nmin
            color_values[neg] = (colors[neg] - nmin) / nspan if nspan > 1e-12 else 0.0

        # normalize positives -> [0, 1]. Same degenerate-span guard, but
        # defaulting to 1.0 (fully-saturated "A" end) — 0.0 here would
        # collide with the zero/unassigned scalar of 0.5 after the 0.5+0.5*
        # mapping below just as badly as on the negative side.
        if np.any(pos):
            pmin, pmax = colors[pos].min(), colors[pos].max()
            pspan = pmax - pmin
            color_values[pos] = (colors[pos] - pmin) / pspan if pspan > 1e-12 else 1.0

        # store sign mask separately (IMPORTANT for colormap)
        polyline["colors_raw"] = colors
        polyline["colors_norm"] = color_values

        # encode sign explicitly in scalars:
        # - negative -> [0, 0.5]
        # - zero     -> exactly 0.5
        # - positive -> [0.5, 1]
        scalar = np.zeros(len(colors), dtype=float)

        if np.any(neg):
            scalar[neg] = 0.5 * color_values[neg]

        if np.any(pos):
            scalar[pos] = 0.5 + 0.5 * color_values[pos]

        scalar[zero] = 0.5

        polyline["colors"] = scalar
        polymer = polyline.tube(radius=r)

        cmap = "coolwarm"  # keep diverging map for rendering

    else:
        logger.info("No coloring applied (uniform rendering)")
        polymer = polyline.tube(radius=r)

    plotter = pv.Plotter(off_screen=True if save_path else False)

    plotter.add_mesh(
        polymer,
        smooth_shading=True,
        cmap=cmap,
        scalars="colors" if colors is not None else None,
        show_scalar_bar=False,
    )

    if zoom != 1.0:
        # reset_camera() first performs the same auto-fit framing `.show()`
        # would otherwise apply on first render; doing it explicitly here
        # lets .zoom() tighten in *beyond* that default fit, rather than
        # zooming from an unset/identity camera that then gets overwritten
        # by .show()'s own auto-fit.
        plotter.reset_camera()
        plotter.camera.zoom(zoom)

    if save_path:
        logger.info(f"Saving visualization to: {save_path}")
        plotter.show(screenshot=save_path)
        if colors is not None and legend_labels is not None:
            _composite_signed_legend(save_path, cmap, legend_labels)
    else:
        logger.info("Displaying visualization interactively")
        plotter.show()

    plotter.close()

    logger.info("Visualization finished")


def _composite_signed_legend(save_path, cmap_name, legend_labels):
    """Overlay a small frameless legend (colored dot proxies at the
    colormap's negative/positive extremes) onto an already-saved
    screenshot. `viz_structure` renders via an off-screen PyVista plotter
    with show_scalar_bar=False, so this is the only way a signed/
    compartment-colored render gets a legible key."""
    neg_label, pos_label = legend_labels[0], legend_labels[1]
    try:
        cmap = plt.get_cmap(cmap_name)
        img = plt.imread(save_path)
        h, w = img.shape[0], img.shape[1]
        fig = plt.figure(figsize=(w / 150, h / 150), dpi=150)
        ax = fig.add_axes([0, 0, 1, 1])
        ax.imshow(img)
        ax.axis("off")
        handles = [
            Line2D([0], [0], marker="o", linestyle="none", markerfacecolor=cmap(1.0),
                   markeredgecolor="none", markersize=12, label=pos_label),
            Line2D([0], [0], marker="o", linestyle="none", markerfacecolor=cmap(0.0),
                   markeredgecolor="none", markersize=12, label=neg_label),
        ]
        ax.legend(handles=handles, loc="upper right", frameon=False, fontsize=13, labelcolor="#2b2b28")
        fig.savefig(save_path, dpi=150)
        plt.close(fig)
    except Exception as exc:
        logger.warning(f"Could not composite legend onto {save_path}: {exc}")

def save_chimera_cmd(start, end, total_residues, cmd_filename="coloring.cmd"):
    """
    Create a Chimera .cmd file:
    - Color residues outside the given region blue.
    - Color residues inside the region red.
    """

    logger.info(
        f"Writing Chimera cmd: {cmd_filename} | "
        f"region={start}-{end} | total_residues={total_residues}"
    )

    with open(cmd_filename, "w") as f:

        # Color all residues blue first (except highlighted region)
        if start > 1:
            logger.info(f"Coloring blue: 1-{start-1}")
            f.write(f"color blue :1-{start-1}\n")

        if end < total_residues:
            logger.info(f"Coloring blue: {end+1}-{total_residues}")
            f.write(f"color blue :{end+1}-{total_residues}\n")

        # Highlight region
        logger.info(f"Coloring red region: {start}-{end}")
        f.write(f"color red :{start}-{end}\n")

        f.write("focus\n")

    logger.info("Chimera cmd file written successfully")

def viz_gene_structure(V, start, end, r=0.1, cmap="coolwarm", save_path=None):
    """Visualize structure V, highlight a continuous region in red, rest in
    blue."""
    polyline = polyline_from_points(V)
    polyline["scalars"] = np.arange(polyline.n_points)

    # Create colors: 0 for blue, 1 for red
    colors = np.zeros(len(V))
    colors[start : end + 1] = 1  # Mark the highlighted region

    polyline["colors"] = colors

    # Create tube
    polymer = polyline.tube(radius=r)

    # Create plotter
    plotter = pv.Plotter(off_screen=True if save_path else False)
    plotter.add_mesh(
        polymer,
        smooth_shading=True,
        scalars="colors",
        cmap=["blue", "red"],  # Explicit color map
        show_scalar_bar=False,
        clim=[0, 1],  # Force colors 0 and 1
    )

    if save_path:
        plotter.show(screenshot=save_path)
    else:
        plotter.show()


def viz_chroms(sim_path, r=0.1, comps=True):
    logger.info(f"Chromosome visualization started: {sim_path}")

    cif_path = sim_path + "model/MultiMM_minimized.cif"
    chrom_idxs_path = sim_path + "metadata/chrom_idxs.npy"
    chrom_comps_path = sim_path + "metadata/compartments.npy"
    chrom_ends_path = sim_path + "metadata/chrom_lengths.npy"

    chrom_idxs = np.load(chrom_idxs_path)
    chrom_ends = np.load(chrom_ends_path)

    logger.info(f"Loaded chrom_idxs: {len(chrom_idxs)}, chrom_ends: {len(chrom_ends)}")

    if comps:
        comps_array = np.load(chrom_comps_path)
        logger.info(f"Loaded compartments array: shape={comps_array.shape}")

    V = get_coordinates_cif(cif_path)
    N = len(V)

    logger.info(f"Structure loaded: N_beads={N}")

    chroms = np.zeros(N)

    for i in range(len(chrom_ends) - 1):
        start, end = chrom_ends[i], chrom_ends[i + 1]
        chroms[start:end] = chrom_idxs[i]

    logger.info(f"Chromosome assignment completed over {len(chrom_ends)-1} segments")

    viz_structure(
        V,
        chroms[: len(V)],
        cmap="gist_ncar",
        r=r,
        save_path=sim_path + "plots/minimized_structure_chromosomes.png",
    )

    logger.info("Chromosome-colored structure saved")

    if comps:
        viz_structure(
            V,
            comps_array[: len(V)],
            cmap="coolwarm",
            r=r,
            save_path=sim_path + "plots/minimized_structure_compartments.png",
            legend_labels=("B (dense)", "A (sparse)"),
        )
        logger.info("Compartment-colored structure saved")

    logger.info("Chromosome visualization finished successfully")

# oe_matrix() now lives in utils.py (the single canonical OE-normalisation
# implementation, shared with validation.py and hic_force.py) — imported above.


def _hic_log2_oe(oe_matrix: np.ndarray):
    """log2(O/E) of an OE-normalised matrix, floored so zero/background
    entries don't produce -inf, plus a mask of which entries were
    genuinely OBSERVED (oe_matrix > 0) rather than floored background.

    The floor itself is necessarily somewhat arbitrary (the 0.1th
    percentile of positive entries), so every floored/zero pixel ends up
    at exactly the same large-magnitude log2 value. That's fine as a
    floor, but it must never be allowed into a *percentile* computation
    used to pick a display clip (see ``_hic_display_transform`` /
    ``_joint_hic_clip``) — enough floored background pixels at that one
    extreme value would drag the percentile there too, producing a clip
    that reflects the floor, not the data's real dynamic range, and
    washing out every genuinely-observed, moderate feature.
    """
    finite_pos = oe_matrix[np.isfinite(oe_matrix) & (oe_matrix > 0)]
    floor = float(np.percentile(finite_pos, 0.1)) if finite_pos.size else 1e-6
    floor = max(floor, 1e-10)
    log2_m = np.log2(np.clip(oe_matrix, floor, None))
    observed = np.isfinite(oe_matrix) & (oe_matrix > 0)
    return log2_m, observed


def _joint_hic_clip(oe_matrices: "list[np.ndarray]", pct: float = 99.0) -> float:
    """Shared colour-scale clip for ``_hic_display_transform``, computed
    from the POOLED, genuinely-observed log2(O/E) values (see
    ``_hic_log2_oe``) across every panel being compared — not just one of
    them.

    A clip sized to only one panel (e.g. only the experimental matrix)
    breaks as soon as another panel's real dynamic range differs: too
    narrow a clip for a panel with stronger contrast collapses it into
    flat, oversaturated colour blocks, while that same narrow clip makes a
    panel with real-but-moderate structure (often the experimental one,
    once denoised) look almost uniformly gray. Pooling means every panel's
    actual signal is represented in the one shared scale, so all of them
    stay legible at once.
    """
    pooled = []
    for oe in oe_matrices:
        log2_m, observed = _hic_log2_oe(oe)
        vals = log2_m[observed & np.isfinite(log2_m)]
        if vals.size:
            pooled.append(vals)
    if not pooled:
        return 1.0
    clip = float(np.percentile(np.abs(np.concatenate(pooled)), pct))
    return clip if clip > 1e-10 else 1.0


def _hic_display_transform(oe_matrix: np.ndarray, clip: "float | None" = None, pct: float = 99.0, gamma: float = 0.6):
    """Shared display normalisation for every Hi-C-style heatmap: log2
    fold-enrichment, clipped to [-1, 1] via a shared ``clip`` (so panels
    being compared map the same log2(O/E) value to the same colour), then
    a sign-preserving power stretch (`sign(x)*|x|**gamma`) for contrast.
    Keeps real magnitude information, unlike a rank transform (which maps
    every matrix to a uniform distribution and erases how much contrast is
    actually there — making noise and real signal look equally "detailed").

    If ``clip`` is not given, it is picked from THIS matrix's own
    genuinely-observed entries only (see ``_hic_log2_oe``) — floored
    background pixels are excluded so they can't drag the clip to an
    arbitrary, non-representative value. For a multi-panel comparison,
    prefer ``_joint_hic_clip`` over several matrices and pass its result in
    here as ``clip`` for each one, so they share one colour scale sized to
    all of them, not just the first.

    Returns (display, clip); pass the returned ``clip`` back in for other
    matrices in the same comparison so they share one colour scale.
    """
    log2_m, observed = _hic_log2_oe(oe_matrix)

    if clip is None:
        src = log2_m[observed & np.isfinite(log2_m)]
        if not src.size:
            src = log2_m[np.isfinite(log2_m)]
        clip = float(np.percentile(np.abs(src), pct)) if src.size else 1.0
        clip = clip if clip > 1e-10 else 1.0

    normalized = np.clip(log2_m / clip, -1.0, 1.0)
    display = np.sign(normalized) * np.abs(normalized) ** gamma
    return display, clip


def get_heatmap(
    cif_file,
    viz=False,
    save=False,
    save_path=None,
    vmax=None,
    vmin=None,
    log_scale=True,
    oe_normalize=True,
    reorder_by_diagonal=False,
    name="structure",
    rc=None,
    alpha=4.0,
    auto_scale=True,
    auto_scale_percentile=10.0,
    kernel="power_law",
):
    """Compute and visualize a contact/interaction heatmap from a 3D structure.

    Uses the same distance->contact-probability inversion as the Hi-C
    Boltzmann-PMF force (:func:`hic_force.get_boltzmann_p_func`), so this
    map and the ensemble-averaged proxy in :func:`plot_hic_comparison` are
    computed by identical methodology and directly comparable.

    rc/alpha/kernel mirror the Hi-C force's own parameters (``kernel`` should
    match HIC_BOLTZMANN_KERNEL for a like-for-like comparison). ``rc=None`` is fine
    when ``auto_scale=True`` (default): the length scale is instead
    recalibrated from a percentile of this structure's own pairwise-distance
    distribution (:func:`hic_force.auto_contact_scale`,
    ``auto_scale_percentile``, default 10.0) rather than the force's
    deliberately microscopic default, which would underflow to ~0 contact
    probability for most pairs here.

    ``oe_normalize`` applies Observed/Expected normalisation (recommended
    for comparison with experimental Hi-C). ``log_scale``, ``vmin``,
    ``vmax`` are accepted for backward compatibility but only affect the
    returned matrix, not the displayed image (always shown via
    ``_hic_display_transform``).
    """

    # ------------------------------------------------------------
    # Output dir
    # ------------------------------------------------------------
    base_dir = save_path
    os.makedirs(base_dir, exist_ok=True)

    def _save_local(fig, name):
        fig.savefig(name + ".png", dpi=300, bbox_inches="tight")

    # ------------------------------------------------------------
    # Load structure
    # ------------------------------------------------------------
    V = get_coordinates_cif(cif_file)
    logger.info(f"Loaded structure: shape={V.shape}, file={cif_file}")

    # ------------------------------------------------------------
    # Distance → contact proxy, via the SAME distance->contact inversion used
    # to build the Hi-C Boltzmann-PMF force, with its length scale
    # auto-calibrated from this structure's own pairwise-distance
    # distribution (see docstring).
    # ------------------------------------------------------------
    D = distance.cdist(V, V, metric="euclidean")

    rc_eff = rc if rc is not None else 1.0
    if auto_scale:
        calibrated = auto_contact_scale(V, percentile=auto_scale_percentile)
        rc_eff = calibrated
        logger.info(
            f"Auto-calibrated contact scale: {calibrated:.4f} nm "
            f"({auto_scale_percentile:.0f}th percentile of pairwise distances)"
        )

    p_func = get_boltzmann_p_func(rc_eff, alpha=alpha, kernel=kernel)
    mat = p_func(D)

    logger.info(
        f"Raw contact matrix: min={mat.min():.3e}, max={mat.max():.3e}, "
        f"mean={mat.mean():.3e}"
    )

    if log_scale:
        mat = np.log1p(mat)
        logger.info("Applied log1p transform to contact matrix")

    # ------------------------------------------------------------
    # OE normalisation — removes distance-decay baseline so the
    # simulated map is directly comparable to experimental Hi-C
    # ------------------------------------------------------------
    if oe_normalize:
        mat = oe_matrix(mat)
        logger.info("Applied OE (Observed/Expected) diagonal normalisation")

    # ------------------------------------------------------------
    # Optional reordering
    # ------------------------------------------------------------
    if reorder_by_diagonal:
        order = np.argsort(np.linalg.norm(V - np.mean(V, axis=0), axis=1))
        mat = mat[np.ix_(order, order)]
        logger.info("Reordered matrix by distance-from-centroid sorting")

    # ------------------------------------------------------------
    # Visualization — signed log2(O/E), clipped + gamma-stretched to
    # [-1, 1] (see _hic_display_transform) with the single shared, soft
    # diverging colormap used for every Hi-C-style heatmap in this module,
    # so this map and plot_hic_comparison's panels are visually and
    # methodologically consistent, not just individually "nice-looking".
    # ------------------------------------------------------------
    if viz:
        sns.set_style("white")
        fig, ax = plt.subplots(figsize=(7, 6))

        disp, _ = _hic_display_transform(mat)
        im = ax.imshow(disp, cmap=_HIC_CMAP, vmin=-1.0, vmax=1.0, interpolation="nearest")

        ax.set_title("Structure-derived Contact Map (OE-normalised)", fontsize=12, fontweight="bold")
        ax.set_xlabel("Bead index")
        ax.set_ylabel("Bead index")

        cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("log2(O/E), gamma-stretched  (enriched ↑ / depleted ↓)")

        ax.set_aspect("equal")
        ax.tick_params(length=0)

        if save and save_path is not None:
            logger.info(f"Saving heatmap to: {save_path}/{name}_contact_map.png")
            _save_local(fig, save_path + f"/{name}_contact_map")

        plt.close(fig)

    logger.info("Heatmap computation finished")
    return mat

def plot_md_thermo(history, save_path, target_temperature=None):
    """
    Plot energy + temperature (+ RMSD, when available) evolution from MD.

    Two stacked panels:
      - top: potential / kinetic / total energy (left axis) and temperature
        (right, twin axis) with ONE combined legend covering every curve.
      - bottom: RMSD relative to the minimised structure, if the history
        contains it (one line is dropped cleanly otherwise).

    `target_temperature` (optional, in kelvin) draws a thin reference line
    at the simulation's set-point temperature so deviations are easy to spot.
    """

    logger.info("Creating MD thermodynamics plot...")

    sns.set_style("whitegrid")

    # explicit numeric cast — plotting raw (possibly mixed-type) history
    # values directly can make matplotlib mistake them for categorical data
    steps = np.asarray(history["step"], dtype=float)
    rmsd = np.asarray(history.get("rmsd", []), dtype=float)
    has_rmsd = len(rmsd) == len(steps) and len(rmsd) > 0

    palette = sns.color_palette("deep")

    if has_rmsd:
        fig, (ax1, ax3) = plt.subplots(
            2, 1, figsize=(9, 7), sharex=True,
            gridspec_kw={"height_ratios": [2.2, 1], "hspace": 0.08},
        )
    else:
        fig, ax1 = plt.subplots(figsize=(9, 5))
        ax3 = None

    # ---- top panel: energies (left axis) + temperature (right axis) ----
    potential = np.asarray(history["potential"], dtype=float)
    kinetic = np.asarray(history["kinetic"], dtype=float)
    total = np.asarray(history["total"], dtype=float)
    temperature = np.asarray(history["temperature"], dtype=float)

    l1, = ax1.plot(steps, potential, color=palette[0], linewidth=1.6,
                    label="Potential energy")
    l2, = ax1.plot(steps, kinetic, color=palette[1], linewidth=1.6,
                    label="Kinetic energy")
    l3, = ax1.plot(steps, total, color=palette[2], linewidth=2.0,
                    label="Total energy")

    ax1.set_ylabel("Energy (kJ/mol)")

    ax2 = ax1.twinx()
    l4, = ax2.plot(steps, temperature, color=palette[3], linestyle="--",
                   linewidth=1.6, label="Temperature")
    ax2.set_ylabel("Temperature (K)")
    ax2.grid(False)

    handles = [l1, l2, l3, l4]

    if target_temperature is not None:
        l5 = ax2.axhline(target_temperature, color="black", linestyle=":", linewidth=1.2,
                          label=f"Target T = {target_temperature:.0f} K")
        handles.append(l5)

    if ax3 is not None:
        # ---- bottom panel: RMSD vs the minimised structure ----
        ax3.plot(steps, rmsd, color=palette[4], linewidth=1.6, label="RMSD vs minimised")
        ax3.fill_between(steps, rmsd, color=palette[4], alpha=0.15)
        ax3.set_xlabel("Step")
        ax3.set_ylabel("RMSD (nm)")
        handles.append(ax3.get_lines()[0])
        ax3.grid(True, alpha=0.4)
    else:
        ax1.set_xlabel("Step")

    ax1.grid(True, alpha=0.4)

    # figure-level title + ONE combined legend covering every curve on
    # both panels and both y-axes, placed above the title so nothing overlaps
    fig.suptitle("MultiMM MD Thermodynamics", fontsize=14, fontweight="bold", y=0.99)
    fig.legend(
        handles=handles, labels=[h.get_label() for h in handles],
        loc="upper center", bbox_to_anchor=(0.5, 0.95),
        ncol=min(len(handles), 3), frameon=True, fontsize=9,
    )

    # subplots_adjust (not tight_layout) because ax2 is a twinx() axes,
    # which tight_layout can't account for and warns about.
    fig.subplots_adjust(top=0.84 if ax3 is not None else 0.88, bottom=0.08, left=0.1, right=0.9)

    out = os.path.join(save_path, "plots/md_thermodynamics.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)

    logger.info(f"MD thermodynamics plot saved to: {out}")


# Canonical term order (matches the "Forcefield — active terms" table in
# model.py's add_forcefield) — fixes each term's color regardless of which
# subset is actually enabled in a given run, so a term's color never changes
# run-to-run just because other terms were toggled off.
_ENERGY_TERM_ORDER = [
    "Excluded volume",
    "Harmonic bonds",
    "Harmonic angles",
    "Loop extrusion",
    "Compartment blocks",
    "Subcompartment blocks",
    "Chromosomal blocks",
    "Spherical container",
    "B-lamina interaction",
    "Central force",
    "Hi-C guided force",
]


def plot_energy_components(history, save_path):
    """Plot each active force term's potential energy vs. time, one line per
    term (EV, bonds, angles, loop extrusion, ... — whatever was enabled for
    this run; see model.py's _register_force/run_md).

    Single axis (all terms share the same kJ/mol units), one fixed color per
    term name (see _ENERGY_TERM_ORDER) so colors stay stable across runs with
    different enabled forces, and a legend since there are always >= 2 terms
    in practice.
    """
    comps = history.get("energy_components", {})
    comps = {name: vals for name, vals in comps.items() if len(vals) == len(history["step"])}
    if not comps:
        logger.info("No energy-component history to plot — skipping.")
        return

    logger.info("Creating energy-components plot...")

    sns.set_style("whitegrid")
    # explicit numeric cast — see plot_md_thermo for why
    steps = np.asarray(history["step"], dtype=float)
    comps = {name: np.asarray(vals, dtype=float) for name, vals in comps.items()}

    # fixed color per canonical term name; any unexpected/extra name still
    # gets a stable color by falling back to its position in the palette.
    # husl (not tab20) so adjacent terms are hue-separated even when only
    # two neighbouring terms end up active in a given run — tab20 pairs
    # adjacent indices as dark/light of the SAME hue, which can look nearly
    # identical at a glance.
    full_palette = sns.color_palette("husl", len(_ENERGY_TERM_ORDER))
    color_of = dict(zip(_ENERGY_TERM_ORDER, full_palette))

    # plot in canonical order first, then any leftover names, so the legend
    # order is always the same regardless of dict insertion order
    ordered_names = [n for n in _ENERGY_TERM_ORDER if n in comps]
    ordered_names += [n for n in comps if n not in _ENERGY_TERM_ORDER]

    fig, ax = plt.subplots(figsize=(9, 5))
    for name in ordered_names:
        color = color_of.get(name, (0.5, 0.5, 0.5))
        ax.plot(steps, comps[name], linewidth=1.6, label=name, color=color)

    ax.set_xlabel("Step")
    ax.set_ylabel("Potential energy (kJ/mol)")
    ax.grid(True, alpha=0.4)
    ax.set_title("MultiMM — Energy Components", fontsize=14, fontweight="bold")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12),
              ncol=min(len(ordered_names), 4), frameon=True, fontsize=9)

    out = os.path.join(save_path, "plots/energy_components.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)

    logger.info(f"Energy-components plot saved to: {out}")

from .structural_analysis import analyze_structure, analyze_dynamics, _hist_with_kde  # noqa: F401  (re-exported for backward compatibility)

# ── Hi-C comparison heatmap ───────────────────────────────────────────────────

def plot_hic_comparison(
    sim_matrix: "np.ndarray",
    exp_matrix: "np.ndarray",
    save_dir: str,
    name: str = "hic_comparison",
    low_percentile: float = 1.0,
    high_percentile: float = 99.0,
    rw_matrix: "np.ndarray | None" = None,
    shared_scale: bool = False,
) -> None:
    """Save a side-by-side heatmap figure comparing simulated vs. experimental
    Hi-C (plus an optional random-walk null-model panel if *rw_matrix* is given).

    Each matrix is OE-normalised then run through ``_hic_display_transform``
    (signed log2(O/E), gamma-stretched). By default (``shared_scale=False``)
    every panel is clipped to its OWN percentile, so each one uses its full
    colour range and stays legible regardless of how its absolute dynamic
    range compares to the others — the simulated contact-proxy map and the
    random-walk null are both far sparser than a denoised experimental Hi-C
    matrix, so a clip sized to (or pooled with) the experimental panel can
    crush them to a near-uniform, washed-out white rather than showing their
    real, smaller-but-genuine structure.

    Set ``shared_scale=True`` to instead use one clip pooled across every
    displayed panel's ``high_percentile`` (see ``_joint_hic_clip``), which
    keeps relative-magnitude comparisons literal at the cost of exactly this
    washing-out risk for whichever panel has the narrower dynamic range.

    sim_matrix/exp_matrix: (N, N) contact matrices at matching resolution.
    low_percentile/high_percentile: colour-clipping percentiles.
    rw_matrix: optional (N, N) random-walk ensemble contact map.
    """
    import os
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors

    os.makedirs(save_dir, exist_ok=True)

    # sim_matrix is the ensemble-averaged contact proxy: each frame is first
    # converted to its own contact-probability heatmap via the Hi-C force's
    # distance->contact inversion, then averaged across frames ("average of
    # heatmaps", not "heatmap of the average structure") — see _accumulate_ensemble /
    # _rw_baseline_contact in validation.py. OE normalisation is applied to
    # the raw (linear) matrices; log2(O/E) below is the display transform.
    sim_oe = oe_matrix(sim_matrix)
    exp_oe = oe_matrix(exp_matrix)
    rw_oe  = oe_matrix(rw_matrix) if rw_matrix is not None else None

    if shared_scale:
        # Pooled clip across every panel (see _joint_hic_clip) — literal
        # magnitude comparison, but a panel with a narrower natural dynamic
        # range gets crushed toward white (see docstring above).
        panels_for_clip = [exp_oe, sim_oe] + ([rw_oe] if rw_oe is not None else [])
        clip = _joint_hic_clip(panels_for_clip, pct=high_percentile)
        exp_disp, _ = _hic_display_transform(exp_oe, clip=clip)
        sim_disp, _ = _hic_display_transform(sim_oe, clip=clip)
        rw_disp, _  = (_hic_display_transform(rw_oe, clip=clip) if rw_oe is not None else (None, None))
        scale_msg = "shared colour scale pooled across"
    else:
        # Independent clip per panel (default) — every panel uses its own
        # genuinely-observed dynamic range, so each stays visible on its own.
        exp_disp, _ = _hic_display_transform(exp_oe, pct=high_percentile)
        sim_disp, _ = _hic_display_transform(sim_oe, pct=high_percentile)
        rw_disp, _  = (_hic_display_transform(rw_oe, pct=high_percentile) if rw_oe is not None else (None, None))
        scale_msg = "independent colour scale per panel across"

    disp_matrices = [exp_disp, sim_disp]
    titles        = ["Experimental Hi-C", "Simulated (ensemble-averaged contact proxy)"]

    if rw_oe is not None:
        disp_matrices.append(rw_disp)
        titles.append("Random Walk (null model)")
        logger.info("Applied OE normalisation to sim, exp, and RW matrices for comparison "
                     "(%s all three)", scale_msg)
    else:
        logger.info("Applied OE normalisation to sim and exp matrices for comparison "
                     "(%s both)", scale_msg)

    n_panels  = len(disp_matrices)
    fig_width = 6 * n_panels   # 12 for 2 panels, 18 for 3

    fig, axes = plt.subplots(
        1, n_panels,
        figsize=(fig_width, 5),
        dpi=150,
        constrained_layout=True,
    )

    if n_panels == 1:
        axes = [axes]

    for ax, mat, title in zip(axes, disp_matrices, titles):
        im = ax.imshow(mat, cmap=_HIC_CMAP, vmin=-1.0, vmax=1.0, origin="upper", aspect="auto")
        ax.set_title(title, fontsize=13, fontweight="bold")
        ax.set_xlabel("Genomic bin", fontsize=11)
        ax.set_ylabel("Genomic bin", fontsize=11)
        ax.tick_params(labelsize=9)
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        scale_label = "shared scale" if shared_scale else "own scale"
        cbar.set_label(f"log2(O/E), {scale_label}, gamma-stretched (enriched ↑ / depleted ↓)", fontsize=9)

    out_path = os.path.join(save_dir, f"{name}.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Saved Hi-C comparison heatmap → {out_path}")


# ── Hi-C denoising preview (before / after Gaussian smoothing) ──────────────

def plot_hic_preprocessing(
    H_before: "np.ndarray",
    H_after: "np.ndarray",
    save_dir: str,
    chrom: "str | None" = None,
    sigma: float = 1.0,
    oe_normalized: bool = False,
    name: str = "hic_preprocessing",
) -> None:
    """Side-by-side heatmap of the Hi-C matrix before vs. after the
    denoising step in ``read_hic.preprocess_hic_matrix`` — the actual c_ij
    target the force (and validation) optimizes against.

    If ``oe_normalized=True``, both panels are already OE-normalized (so
    "before" means "OE-normalized, not yet denoised", not "raw") — keeps
    both panels on the same scale. Colour scale is taken from the denoised
    (after) matrix's 99th percentile, so noise pixels clip to the top
    colour in the "before" panel instead of washing out real structure.

    H_before/H_after: (N, N), same shape/scale. chrom: optional label.
    sigma: Gaussian std (beads) used, shown in the title. oe_normalized:
    whether both panels are OE enrichment rather than raw counts.
    """
    import os
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(save_dir, exist_ok=True)

    if np.any(H_after > 0):
        vmax = float(np.percentile(H_after[H_after > 0], 99.0))
    elif np.any(H_before > 0):
        # Denoising left nothing (degenerate/empty case) — fall back to the
        # raw matrix's own scale rather than leaving vmax at 0.
        vmax = float(np.percentile(H_before[H_before > 0], 99.0))
    else:
        vmax = 1.0
    vmax = max(vmax, 1e-12)

    chrom_label = f" ({chrom})" if chrom else ""
    quantity = "OE-normalized Hi-C" if oe_normalized else "Raw Hi-C"
    cbar_label = "OE enrichment (above background)" if oe_normalized else "contact strength"

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), dpi=150, constrained_layout=True)

    for ax, mat, title in zip(
        axes,
        [H_before, H_after],
        [f"{quantity}{chrom_label}, before denoising",
         f"{quantity}{chrom_label}, denoised  (median+Gaussian σ={sigma:g})"],
    ):
        im = ax.imshow(mat, cmap="Reds", vmin=0.0, vmax=vmax, origin="upper", aspect="equal",
                        interpolation="nearest")
        ax.set_title(title, fontsize=12, fontweight="bold")
        ax.set_xlabel("Genomic bin", fontsize=10)
        ax.set_ylabel("Genomic bin", fontsize=10)
        ax.tick_params(labelsize=8)
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, extend="max")
        cbar.set_label(cbar_label, fontsize=9)

    out_path = os.path.join(save_dir, f"{name}.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Saved Hi-C before/after denoising preview → {out_path}")


# ── Hi-C validation curves (decay / insulation / PC1) ────────────────────────

def plot_hic_validation_curves(
    sim_decay, exp_decay, rw_decay,
    sim_ins, exp_ins, rw_ins,
    sim_pc1, exp_pc1, rw_pc1,
    save_dir, name="hic_validation_curves",
    r_dd=None, r_dd_rw=None,
    r_ins=None, r_ins_rw=None,
    r_pc1=None, r_pc1_rw=None,
):
    """One figure with the three Hi-C validation curves side by side:
    diagonal decay, insulation score, and PC1 (A/B compartment), each
    comparing MultiMM against experimental Hi-C and an optional RW baseline.

    Decay is plotted log-log normalised to 1.0 at the shortest separation
    (display-only, not fed into r_dd). Insulation and PC1 arrive from
    validation.py already smoothed and range-normalised — [0, 1] for
    insulation, [-1, 1] (sign-preserving) for PC1 — the exact same signal
    the reported r_ins/r_pc1 was computed on, so what's plotted always
    matches what's scored. PC1 also arrives already sign-aligned to the
    real A/B convention (validation.py: utils.align_pc1_sign, against each
    matrix's own contact density), so no further re-alignment happens here
    — the plotted curves match the signed r_pc1 reported.
    """
    os.makedirs(save_dir, exist_ok=True)
    sns.set_style("whitegrid")

    ink_exp = "#0b0b0b"       # experimental = ground truth, neutral ink
    color_sim = "#2a78d6"     # MultiMM = categorical slot 1 (blue)
    color_rw = "#898781"      # RW null baseline = muted gray, dashed

    has_rw = rw_decay is not None

    def _log_bin_average(x, y, n_bins=60):
        """Average *y* over log-spaced bins of *x* (smooths the noisy tail
        of diagonal-decay curves at large genomic separation)."""
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        mask = (x > 0) & np.isfinite(y)
        x, y = x[mask], y[mask]
        if len(x) < 4:
            return x, y
        edges = np.logspace(np.log10(x.min()), np.log10(x.max()), n_bins + 1)
        idx = np.digitize(x, edges)
        xs, ys = [], []
        for b in range(1, len(edges)):
            m = idx == b
            if m.sum() > 0:
                xs.append(x[m].mean())
                ys.append(y[m].mean())
        return np.array(xs), np.array(ys)

    fig, axes = plt.subplots(1, 3, figsize=(17, 5))

    # ---- (1) diagonal decay, log-log — normalised to 1.0 at shortest separation ----
    # Isolates decay rate from each curve's unrelated absolute scale.
    ax = axes[0]
    s = np.arange(1, len(exp_decay) + 1)

    def _plot_decay(s, y, color, label, linestyle="-", linewidth=2.0, raw_alpha=0.18):
        xs_raw, ys_raw = s, np.asarray(y, dtype=float)
        ref = ys_raw[np.isfinite(ys_raw) & (ys_raw > 0)]
        ref_val = ref[0] if ref.size else 1.0
        ax.loglog(xs_raw, ys_raw / ref_val, color=color, linewidth=1.0, alpha=raw_alpha)
        xs, ys = _log_bin_average(s, y)
        ys = ys / ref_val
        ax.loglog(xs, ys, color=color, linestyle=linestyle, linewidth=linewidth, label=label)

    _plot_decay(s, exp_decay, ink_exp, "Experimental")
    _plot_decay(s, sim_decay, color_sim,
                f"MultiMM (r={r_dd:.2f})" if r_dd is not None else "MultiMM")
    if has_rw:
        _plot_decay(s, rw_decay, color_rw,
                    f"Random walk (r={r_dd_rw:.2f})" if r_dd_rw is not None else "Random walk",
                    linestyle="--", linewidth=1.6)
    ax.axhline(1.0, color="#c3c2b7", linewidth=0.8, linestyle=":", zorder=0)
    ax.set_title("Diagonal Decay (log-log, normalised)", fontsize=12, fontweight="bold")
    ax.set_xlabel("Genomic separation s (beads)")
    ax.set_ylabel("Mean contact (normalised to s=1)")
    ax.legend(fontsize=9, frameon=False)

    # ---- (2) insulation score — already smoothed + [0, 1] normalised ----
    # (see validation.py: same signal the reported r_ins was computed on)
    ax = axes[1]
    idx = np.arange(len(exp_ins))

    exp_ins_s = np.asarray(exp_ins, dtype=float)
    sim_ins_s = np.asarray(sim_ins, dtype=float)

    mask = np.isfinite(sim_ins_s) & np.isfinite(exp_ins_s)
    if has_rw:
        rw_ins_s = np.asarray(rw_ins, dtype=float)
        mask &= np.isfinite(rw_ins_s)

    ax.plot(idx[mask], exp_ins_s[mask], color=ink_exp, linewidth=1.8, label="Experimental")
    ax.plot(idx[mask], sim_ins_s[mask], color=color_sim, linewidth=1.6, alpha=0.9,
            label=f"MultiMM (r={r_ins:.2f})" if r_ins is not None else "MultiMM")
    ax.fill_between(idx[mask], exp_ins_s[mask], sim_ins_s[mask],
                     color=color_sim, alpha=0.08)
    if has_rw:
        ax.plot(idx[mask], rw_ins_s[mask], color=color_rw, linestyle="--", linewidth=1.3,
                label=f"Random walk (r={r_ins_rw:.2f})" if r_ins_rw is not None else "Random walk")
    ax.set_title("Insulation Score (TAD boundaries, normalised)", fontsize=12, fontweight="bold")
    ax.set_xlabel("Bead index")
    ax.set_ylabel("Insulation score (normalised 0–1)")
    ax.legend(fontsize=9, frameon=False)

    # ---- (3) PC1 — A/B compartment signal (already smoothed, [-1,1] ----
    # normalised, and sign-aligned to the real A/B convention; see
    # validation.py / utils.align_pc1_sign). Plotted as-is — no re-alignment
    # here, so the curves match the signed r_pc1 reported.
    ax = axes[2]

    exp_pc1_s = np.asarray(exp_pc1, dtype=float)
    sim_pc1_s = np.asarray(sim_pc1, dtype=float)

    idx = np.arange(len(exp_pc1_s))
    ax.axhline(0, color="#c3c2b7", linewidth=1.0)
    ax.plot(idx, exp_pc1_s, color=ink_exp, linewidth=1.8, label="Experimental")
    ax.plot(idx, sim_pc1_s, color=color_sim, linewidth=1.6, alpha=0.9,
            label=f"MultiMM (r={r_pc1:.2f})" if r_pc1 is not None else "MultiMM")
    ax.fill_between(idx, exp_pc1_s, sim_pc1_s, color=color_sim, alpha=0.08)
    if has_rw:
        rw_pc1_s = np.asarray(rw_pc1, dtype=float)
        ax.plot(idx, rw_pc1_s, color=color_rw, linestyle="--", linewidth=1.3,
                label=f"Random walk (r={r_pc1_rw:.2f})" if r_pc1_rw is not None else "Random walk")
    ax.set_ylim(-1.05, 1.05)
    ax.set_title("PC1 — A/B Compartment Signal (normalised)", fontsize=12, fontweight="bold")
    ax.set_xlabel("Bead index")
    ax.set_ylabel("PC1 (normalised to [-1, 1], sign-aligned)")
    ax.legend(fontsize=9, frameon=False)

    fig.suptitle("Hi-C Validation Curves", fontsize=15, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.94])

    out_path = os.path.join(save_dir, f"{name}.png")
    fig.savefig(out_path, dpi=250, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Saved Hi-C validation curves → {out_path}")


# ── Input-vs-output diagnostics: loops & compartments ────────────────────────

def plot_loop_validation(loop_dist, bg_dist, loop_sep, bg_sep, save_dir,
                          name="loop_validation", fold_closer=None, p_value=None):
    """One consolidated figure checking whether input loop anchors (.bedpe)
    ended up close together in 3D, relative to a genomic-separation-matched
    background of non-loop pairs.
    """
    os.makedirs(save_dir, exist_ok=True)
    sns.set_style("whitegrid")

    color_loop = "#2a78d6"
    color_bg = "#898781"

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # ---- (1) distance distributions ----
    ax = axes[0]
    df = pd.DataFrame({
        "distance": np.concatenate([loop_dist, bg_dist]),
        "group": ["Loop anchors"] * len(loop_dist) + ["Background (matched)"] * len(bg_dist),
    })
    sns.violinplot(data=df, x="group", y="distance", ax=ax, cut=0, inner=None,
                    palette=[color_loop, color_bg])
    sns.boxplot(data=df, x="group", y="distance", ax=ax, width=0.12,
                showcaps=True, boxprops={"facecolor": "white", "alpha": 0.7},
                whiskerprops={"linewidth": 1.2}, showfliers=False)
    ax.set_xlabel("")
    ax.set_ylabel("3D distance")
    title = "Loop Anchors vs Separation-Matched Background"
    if fold_closer is not None and np.isfinite(fold_closer):
        subtitle = f"{fold_closer:.2f}× closer"
        if p_value is not None and np.isfinite(p_value):
            subtitle += f"  (p={p_value:.1e})"
        title += f"\n{subtitle}"
    ax.set_title(title, fontsize=11, fontweight="bold")

    # ---- (2) distance vs genomic separation, loops over the expected scaling ----
    ax = axes[1]
    order = np.argsort(bg_sep)
    bg_sep_sorted, bg_dist_sorted = bg_sep[order], bg_dist[order]
    # binned running mean = the "expected" scaling curve
    bins = np.linspace(bg_sep_sorted.min(), bg_sep_sorted.max(), 25)
    bin_idx = np.digitize(bg_sep_sorted, bins)
    bin_centers, bin_means = [], []
    for b in range(1, len(bins)):
        m = bin_idx == b
        if m.sum() > 0:
            bin_centers.append(bg_sep_sorted[m].mean())
            bin_means.append(bg_dist_sorted[m].mean())

    ax.scatter(bg_sep, bg_dist, s=6, alpha=0.15, color=color_bg, label="Background pairs")
    ax.plot(bin_centers, bin_means, color=color_bg, linewidth=2.0, label="Expected (background mean)")
    ax.scatter(loop_sep, loop_dist, s=16, alpha=0.85, color=color_loop, label="Loop anchors",
               edgecolors="white", linewidths=0.3)
    ax.set_xlabel("Genomic separation (beads)")
    ax.set_ylabel("3D distance")
    ax.set_title("Loop Anchors vs Expected Polymer Scaling", fontsize=11, fontweight="bold")
    ax.legend(fontsize=9, frameon=False)

    fig.suptitle("Loop Validation (input .bedpe vs output structure)", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])

    out_path = os.path.join(save_dir, f"{name}.png")
    fig.savefig(out_path, dpi=250, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Saved loop validation plot → {out_path}")


def plot_compartment_validation(pc1_aligned, Cs, save_dir, name="compartment_validation",
                                 r=None, accuracy=None, sens_a=None, sens_b=None):
    """One consolidated figure checking whether the structure's own derived
    PC1 signal agrees with the input compartment track (.bed).
    """
    os.makedirs(save_dir, exist_ok=True)
    sns.set_style("whitegrid")

    diverging_cmap = plt.get_cmap("RdBu_r")
    unique_sub = np.sort(np.unique(Cs))
    abs_max = np.max(np.abs(unique_sub)) if len(unique_sub) else 1.0
    norm = mcolors.Normalize(vmin=-abs_max, vmax=abs_max)
    labels = comp_dict  # {-2: "B2", -1: "B1", 1: "A2", 2: "A1", 0: "no compartment"}

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # ---- (1) PC1 by input compartment label ----
    ax = axes[0]
    df = pd.DataFrame({"pc1": pc1_aligned, "subcomp": Cs})
    sns.violinplot(
        data=df, x="subcomp", y="pc1", order=unique_sub, ax=ax, cut=0, inner=None,
        palette=[diverging_cmap(norm(v)) for v in unique_sub],
    )
    sns.stripplot(data=df, x="subcomp", y="pc1", order=unique_sub, ax=ax,
                  color="black", alpha=0.2, size=1.5)
    ax.axhline(0, color="#c3c2b7", linewidth=1.0)
    ax.set_xticks(range(len(unique_sub)))
    ax.set_xticklabels([labels.get(v, str(v)) for v in unique_sub])
    ax.set_xlabel("Input compartment label")
    ax.set_ylabel("Structure-derived PC1 (sign-aligned)")
    title = "Derived PC1 by Input Compartment"
    if r is not None and np.isfinite(r):
        title += f"\n|r| = {r:.2f}"
    ax.set_title(title, fontsize=11, fontweight="bold")

    # ---- (2) agreement scorecard ----
    ax = axes[1]
    metrics = [
        ("Overall sign\naccuracy", accuracy),
        ("A-compartment\nsensitivity", sens_a),
        ("B-compartment\nsensitivity", sens_b),
    ]
    metrics = [(lbl, val) for lbl, val in metrics if val is not None and np.isfinite(val)]
    if metrics:
        labels_, values_ = zip(*metrics)
        bars = ax.bar(labels_, [v * 100 for v in values_], color="#2a78d6", alpha=0.85,
                       edgecolor="white", width=0.55)
        for b, v in zip(bars, values_):
            ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 1.5, f"{v*100:.0f}%",
                    ha="center", fontsize=10, fontweight="bold")
        ax.axhline(50, color="#c3c2b7", linestyle="--", linewidth=1.0, label="Chance level")
        ax.set_ylim(0, 105)
        ax.legend(fontsize=9, frameon=False, loc="lower right")
    ax.set_ylabel("Agreement (%)")
    ax.set_title("Input ↔ Output Agreement", fontsize=11, fontweight="bold")

    fig.suptitle("Compartment Validation (input .bed vs output structure)", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])

    out_path = os.path.join(save_dir, f"{name}.png")
    fig.savefig(out_path, dpi=250, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Saved compartment validation plot → {out_path}")


def plot_distance_vs_strength(strength, dist, low_mask, high_mask, results,
                               save_dir, name="distance_vs_strength", seed=0):
    """Two-panel diagnostic for validate_distance_vs_strength: does higher
    experimental Hi-C strength actually give a smaller 3-D distance?

    Left: classification — distance distributions for the low- vs
    high-strength groups (violin + box), with each group's median distance
    labelled directly and a bracket showing the fold-difference between
    them, so the separation reads at a glance rather than only from the
    AUC number. Right: regression — a (subsampled) scatter of strength vs
    distance, plus a quantile-binned median trend line (with an IQR band)
    so the strength→distance relationship is visible through the scatter's
    skewed density, on a symlog x-axis.

    Strength here is OE-enrichment above background, floored at 0: most
    sampled pairs sit at or very near exactly 0, with a long thin tail of
    genuinely enriched pairs. On a linear axis that piles almost every
    point into a hairline band at x=0 ("everything collapses in a line").
    A symlog scale (linear near 0, logarithmic beyond a small threshold)
    keeps the floored mass readable while spreading the enriched tail out,
    and the binned trend line makes the regression legible despite the
    uneven point density.
    """
    os.makedirs(save_dir, exist_ok=True)
    sns.set_style("whitegrid")

    strength = np.asarray(strength, dtype=float)
    dist = np.asarray(dist, dtype=float)
    low_mask = np.asarray(low_mask, dtype=bool)
    high_mask = np.asarray(high_mask, dtype=bool)

    color_low = "#898781"   # background/depleted — muted gray
    color_high = "#c0392b"  # enriched — warm red (attractive)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

    # ---- (1) classification: distance distribution, low- vs high-strength ----
    ax = axes[0]
    low_dist, high_dist = dist[low_mask], dist[high_mask]
    df = pd.DataFrame({
        "distance": np.concatenate([low_dist, high_dist]),
        "group": (["Low strength\n(background/depleted)"] * int(low_mask.sum())
                  + ["High strength\n(enriched)"] * int(high_mask.sum())),
    })
    if len(df):
        sns.violinplot(data=df, x="group", y="distance", ax=ax, cut=0, inner=None,
                        palette=[color_low, color_high], width=0.8)
        sns.boxplot(data=df, x="group", y="distance", ax=ax, width=0.10,
                    showcaps=True, boxprops={"facecolor": "white", "alpha": 0.85},
                    whiskerprops={"linewidth": 1.2}, showfliers=False)

    if low_dist.size and high_dist.size:
        med_low, med_high = float(np.median(low_dist)), float(np.median(high_dist))
        top_low, top_high = float(np.nanmax(low_dist)), float(np.nanmax(high_dist))
        y_top = max(top_low, top_high)
        y_span = y_top - float(np.nanmin(df["distance"])) if len(df) else 1.0
        label_y_low = top_low + 0.04 * max(y_span, 1e-6)
        label_y_high = top_high + 0.04 * max(y_span, 1e-6)
        bracket_y = y_top + 0.14 * max(y_span, 1e-6)
        text_y = bracket_y + 0.05 * max(y_span, 1e-6)
        box_kw = dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="none", alpha=0.85)
        # per-group median label, placed just above that group's own violin top
        ax.text(0, label_y_low, f"median={med_low:.2f}", ha="center", va="bottom",
                fontsize=9, fontweight="bold", color=color_low, bbox=box_kw)
        ax.text(1, label_y_high, f"median={med_high:.2f}", ha="center", va="bottom",
                fontsize=9, fontweight="bold", color=color_high, bbox=box_kw)
        # bracket + fold-difference between the two medians, above both labels
        ax.plot([0, 0, 1, 1], [bracket_y, bracket_y + 0.01 * y_span,
                                bracket_y + 0.01 * y_span, bracket_y],
                color="#333333", linewidth=1.0)
        if med_high > 1e-12:
            fold = med_low / med_high
            ax.text(0.5, text_y, f"{fold:.1f}× closer", ha="center", va="bottom",
                    fontsize=10, fontweight="bold", color="#333333")
        ax.set_ylim(top=text_y + 0.15 * max(y_span, 1e-6))

    ax.set_xlabel("")
    ax.set_ylabel("3D distance")
    auc = results.get("auc", float("nan"))
    p_mw = results.get("mannwhitney_p", float("nan"))
    title = "Does high strength → smaller distance?"
    if np.isfinite(auc):
        title += f"\nAUC={auc:.2f}" + (f"  (p={p_mw:.1e})" if np.isfinite(p_mw) else "")
    ax.set_title(title, fontsize=11, fontweight="bold")

    # ---- (2) regression: strength vs distance, plain scatter ----
    ax = axes[1]
    finite = np.isfinite(strength) & np.isfinite(dist)
    s_plot, d_plot = strength[finite], dist[finite]
    max_points = 4000
    if s_plot.size > max_points:
        rng = np.random.default_rng(seed)
        sel = rng.choice(s_plot.size, size=max_points, replace=False)
        s_plot, d_plot = s_plot[sel], d_plot[sel]
    ax.scatter(s_plot, d_plot, s=10, alpha=0.25, color="#2a78d6",
               edgecolors="none", zorder=2, label="pairs (subsampled)")

    # symlog x-axis: linear in [-linthresh, linthresh] (so the floored-at-0
    # mass isn't discarded or infinitely compressed), logarithmic beyond it
    # (so the enriched tail actually spreads out instead of hugging x=0).
    pos = s_plot[s_plot > 0]
    if pos.size >= 5:
        linthresh = max(float(np.percentile(pos, 5)), 1e-6)
    else:
        spread = float(np.nanmax(s_plot) - np.nanmin(s_plot)) if s_plot.size else 1.0
        linthresh = max(spread * 1e-3, 1e-6)
    ax.set_xscale("symlog", linthresh=linthresh, linscale=1.5)

    # Quantile-binned median trend (equal-count bins, so the degenerate
    # pile-up near 0 doesn't just become one giant uninformative bin): makes
    # the strength→distance relationship legible despite the skewed density.
    n_bins = int(np.clip(s_plot.size // 150, 6, 14))
    if s_plot.size >= 3 * n_bins:
        order = np.argsort(s_plot)
        s_sorted, d_sorted = s_plot[order], d_plot[order]
        edges = np.linspace(0, s_sorted.size, n_bins + 1).astype(int)
        bin_x, bin_med, bin_lo, bin_hi = [], [], [], []
        for i in range(n_bins):
            sl = slice(edges[i], edges[i + 1])
            if sl.stop - sl.start < 2:
                continue
            bin_x.append(float(np.median(s_sorted[sl])))
            bin_med.append(float(np.median(d_sorted[sl])))
            bin_lo.append(float(np.percentile(d_sorted[sl], 25)))
            bin_hi.append(float(np.percentile(d_sorted[sl], 75)))
        if bin_x:
            ax.fill_between(bin_x, bin_lo, bin_hi, color="#d35400", alpha=0.18,
                             zorder=3, linewidth=0, label="binned IQR")
            ax.plot(bin_x, bin_med, color="#d35400", linewidth=2.2, marker="o",
                    markersize=4.5, zorder=4, label="binned median")

    ax.set_xlabel("Experimental strength (c_ij target, symlog scale)")
    ax.set_ylabel("3D distance")
    ax.legend(loc="best", fontsize=8, framealpha=0.85)
    r_p = results.get("pearson_r", float("nan"))
    r_s = results.get("spearman_r", float("nan"))
    subtitle = []
    if np.isfinite(r_p):
        subtitle.append(f"Pearson r={r_p:.2f}")
    if np.isfinite(r_s):
        subtitle.append(f"Spearman ρ={r_s:.2f}")
    ax.set_title("Strength vs Distance"
                 + (f"\n{', '.join(subtitle)}" if subtitle else ""),
                 fontsize=11, fontweight="bold")

    fig.suptitle("Distance vs Experimental Hi-C Strength", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])

    out_path = os.path.join(save_dir, f"{name}.png")
    fig.savefig(out_path, dpi=250, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Saved distance-vs-strength validation plot → {out_path}")
