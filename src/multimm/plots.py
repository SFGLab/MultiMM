import logging
import os

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyvista as pv
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.pyplot import figure
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import ConvexHull, distance
from scipy.stats import gaussian_kde
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

from .utils import get_coordinates_cif

logger = logging.getLogger(__name__)


pv.set_jupyter_backend("server")
color_dict = {-2: "#bf0020", -1: "#e36a24", 1: "#20c8e6", 2: "#181385", 0: "#ffffff"}
comp_dict = {-2: "B2", -1: "B1", 1: "A2", 2: "A1", 0: "no compartment"}


def plot_projection(struct_3D, Cs, save_path):
    """Chromatin structural analysis centered on COM.

    Includes:
    - PCA embedding
    - 3D structure
    - radial compaction (COM-based)
    - anisotropy via gyration tensor
    - density landscapes
    - subcompartment-dependent structure
    """
    sns.set_style("whitegrid")
    plt.rcParams.update({"figure.dpi": 600, "font.size": 11})

    # preprocessing (STRICT ALIGNMENT GUARANTEE)
    X = np.asarray(struct_3D, dtype=np.float64)
    Cs = np.asarray(Cs)

    N = min(len(X), len(Cs))
    X = X[:N]
    Cs = Cs[:N]

    mask = np.isfinite(X).all(axis=1)

    X = X[mask]
    Cs = Cs[mask]

    # remove invalid compartments early (IMPORTANT)
    valid = Cs != 0
    X = X[valid]
    Cs = Cs[valid]

    # CENTER OF MASS SHIFT (IMPORTANT CHANGE)
    com = X.mean(axis=0)
    Xc = X - com  # everything now COM-based

    # PCA (COM-centered)
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(Xc)
    r = np.linalg.norm(Xc, axis=1)
    X0 = Xc - Xc.mean(axis=0)
    G = (X0.T @ X0) / len(X0)
    eigvals = np.linalg.eigvalsh(G)
    anisotropy_scalar = np.sqrt(eigvals.max() / (eigvals.min() + 1e-12))

    df = pd.DataFrame(
        {
            "x": Xc[:, 0],
            "y": Xc[:, 1],
            "z": Xc[:, 2],
            "pc1": X_pca[:, 0],
            "pc2": X_pca[:, 1],
            "r_com": r,
            "anisotropy": anisotropy_scalar,
            "subcomp": Cs,
        }
    )

    #df = df[df["subcomp"] != 0]

    # output
    base = os.path.join(save_path, "plots")
    os.makedirs(base, exist_ok=True)

    def save(fig, name):
        fig.savefig(os.path.join(base, name + ".png"), dpi=600)
        fig.savefig(os.path.join(base, name + ".pdf"), dpi=600)
        plt.close(fig)

    # 1. PCA projection
    fig, ax = plt.subplots(figsize=(7, 6))
    sc = ax.scatter(df.pc1, df.pc2, c=df.subcomp, s=10, cmap="Spectral", alpha=0.7)

    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label("Subcompartment state")

    ax.set_title("Chromatin PCA (COM-centered configuration)")
    ax.set_xlabel("PC1 (collective mode)")
    ax.set_ylabel("PC2 (collective mode)")

    save(fig, "pca_projection")

    # 2. 3D structure (COM-centered)
    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot(111, projection="3d")

    sc = ax.scatter(Xc[:, 0], Xc[:, 1], Xc[:, 2], c=df.subcomp, cmap="Spectral", s=4, alpha=0.7)

    cbar = fig.colorbar(sc, ax=ax, shrink=0.6)
    cbar.set_label("Subcompartment state")
    ax.set_title("3D Chromatin Structure (center-of-mass frame)")
    ax.set_xlabel("X - COM")
    ax.set_ylabel("Y - COM")
    ax.set_zlabel("Z - COM")
    save(fig, "structure_3D_com")

    # 3. Radial compaction (COM-based)
    fig, ax = plt.subplots(figsize=(7, 4))

    sns.kdeplot(data=df, x="r_com", hue="subcomp", fill=True, alpha=0.5, palette="Spectral", ax=ax)

    ax.set_title("Radial Compaction from Center of Mass")
    ax.set_xlabel("Distance from COM")
    ax.set_ylabel("Density")
    save(fig, "radial_com")

    # 4. PCA density landscape
    fig, ax = plt.subplots(figsize=(7, 6))

    sns.kdeplot(x=df.pc1, y=df.pc2, cmap="mako", fill=True, levels=60, ax=ax)

    ax.set_title("Free-energy-like landscape (PCA space)")
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    save(fig, "pca_density")

    # 5. radial vs subcompartment (IMPROVED: distribution + raw structure)
    fig, ax = plt.subplots(figsize=(7, 4))
    unique_sub = np.sort(df.subcomp.unique())
    abs_max = np.max(np.abs(unique_sub)) if len(unique_sub) > 0 else 1.0
    norm = mcolors.Normalize(vmin=-abs_max, vmax=abs_max)
    cmap = plt.get_cmap("coolwarm")
    sns.violinplot(
        data=df, x="subcomp", y="r_com", palette=[cmap(norm(v)) for v in unique_sub], inner=None, cut=0, ax=ax
    )
    sns.stripplot(data=df, x="subcomp", y="r_com", color="black", alpha=0.25, size=1.5, ax=ax)
    ax.set_title("Radial Distribution by Subcompartment (COM frame)")
    ax.set_xlabel("Subcompartment state")
    ax.set_ylabel("Distance from COM")
    save(fig, "radial_by_subcomp")

    # 7. axis correlations (structure signature, COM-centered, density-based)
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    pairs = [
        (df.x, df.y, "X–Y plane"),
        (df.x, df.z, "X–Z plane"),
        (df.y, df.z, "Y–Z plane"),
    ]
    for ax, (a, b, title) in zip(axes, pairs):
        # replace scatter with 2D density (structure, not noise)
        sns.kdeplot(x=a, y=b, ax=ax, fill=True, cmap="mako", levels=40, thresh=0.05, bw_adjust=0.6)
        # light contour overlay (adds geometry readability)
        sns.kdeplot(x=a, y=b, ax=ax, color="white", levels=6, linewidths=0.6, alpha=0.6, bw_adjust=0.6)
        ax.set_title(title)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.grid(True, alpha=0.2)
    # shared labels (cleaner than repeating)
    axes[0].set_ylabel("Coordinate axis (nm)")
    axes[1].set_xlabel("Coordinate axis (nm)")
    fig.suptitle("Coordinate Correlations in COM frame (density representation)", y=1.02)
    save(fig, "axis_correlations")

    # 8. PCA KDE per subcompartment (signed colors, white background)
    fig, ax = plt.subplots(figsize=(7, 6))
    # enforce clean white background (important for KDE readability)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    x = df.pc1.values
    y = df.pc2.values
    Xg, Yg = np.mgrid[x.min() : x.max() : 200j, y.min() : y.max() : 200j]
    pos = np.vstack([Xg.ravel(), Yg.ravel()])
    unique_sub = np.sort(df.subcomp.unique())

    # ------------------------------------------------------------
    # SIGN-BASED colormap (this is the key fix)
    # negative → blue, positive → red
    # ------------------------------------------------------------
    abs_max = np.max(np.abs(unique_sub)) if len(unique_sub) > 0 else 1.0
    norm = mcolors.Normalize(vmin=-abs_max, vmax=abs_max)
    cmap = plt.get_cmap("coolwarm")

    for scv in unique_sub:

        sub = df[df.subcomp == scv]
        if len(sub) < 10:
            continue

        kde = gaussian_kde([sub.pc1, sub.pc2])
        Z = kde(pos).reshape(Xg.shape)

        color = cmap(norm(scv))

        ax.contourf(Xg, Yg, Z, levels=3, alpha=0.10, colors=[color])

        ax.contour(Xg, Yg, Z, levels=5, colors=[color], linewidths=1.2, alpha=0.9)

    # ------------------------------------------------------------
    # legend (sign-based meaning preserved)
    # ------------------------------------------------------------
    legend_elements = [Line2D([0], [0], color=cmap(norm(v)), lw=2, label=f"subcomp {v}") for v in unique_sub if v != 0]
    ax.legend(handles=legend_elements, frameon=True, fontsize=9)
    ax.set_title("Subcompartment density in PCA space")
    ax.set_xlabel("PC1 (collective chromatin mode)")
    ax.set_ylabel("PC2 (collective chromatin mode)")
    save(fig, "pca_kde_subcomp")


def _save_plotter(plotter, save_path):
    """Save PyVista scene in multiple formats."""
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


def viz_structure(V, colors=None, r=0.1, cmap="coolwarm", save_path=None):
    """Visualize structure V and optionally save it to a file."""
    logger.info(f"Visualizing structure: N={len(V)}, " f"colored={colors is not None}, " f"save_path={save_path}")

    polyline = polyline_from_points(V)
    polyline["scalars"] = np.arange(polyline.n_points)

    if colors is not None and len(colors) > 0:

        colors = np.array(colors[: len(V)])

        logger.info("Color mapping enabled (signed scheme: neg/zero/pos)")

        logger.info(f"Color mapping enabled: min={colors_min}, max={colors_max}, diff={diff}")
        # ------------------------------------------------------------
        # NEW: signed piecewise normalization
        # ------------------------------------------------------------

        color_values = np.zeros(len(colors), dtype=float)

        neg = colors < 0
        pos = colors > 0
        zero = colors == 0

        # normalize negatives -> [0, 1]
        if np.any(neg):
            nmin, nmax = colors[neg].min(), colors[neg].max()
            color_values[neg] = (colors[neg] - nmin) / (nmax - nmin + 1e-12)

        # normalize positives -> [0, 1]
        if np.any(pos):
            pmin, pmax = colors[pos].min(), colors[pos].max()
            color_values[pos] = (colors[pos] - pmin) / (pmax - pmin + 1e-12)

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

    if save_path:
        logger.info(f"Saving visualization to: {save_path}")
        plotter.show(screenshot=save_path)
    else:
        logger.info("Displaying visualization interactively")
        plotter.show()

    plotter.close()

    logger.info("Visualization finished")


def save_chimera_cmd(start, end, total_residues, cmd_filename="coloring.cmd"):
    """
    Create a Chimera .cmd file:
    - Color residues outside the given region blue.
    - Color residues inside the region red.
    """
    logger.info(f"Writing Chimera cmd: {cmd_filename} | " f"region={start}-{end} | total_residues={total_residues}")

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
        )
        logger.info("Compartment-colored structure saved")

    logger.info("Chromosome visualization finished successfully")


def get_heatmap(
    cif_file,
    viz=False,
    save=False,
    save_path=None,
    vmax=None,
    vmin=None,
    log_scale=True,
    reorder_by_diagonal=False,
    name="structure",
):
    """Compute and visualize contact/interaction heatmap from 3D structure."""
    # ------------------------------------------------------------
    # Output dir (UNCHANGED)
    # ------------------------------------------------------------
    base_dir = save_path
    os.makedirs(base_dir, exist_ok=True)

    def _save_local(fig, name):
        path = name
        fig.savefig(path + ".png", dpi=300)
        fig.savefig(path + ".pdf")
        fig.savefig(path + ".svg")

    # ------------------------------------------------------------
    # Load structure
    # ------------------------------------------------------------
    V = get_coordinates_cif(cif_file)
    logger.info(f"Loaded structure: shape={V.shape}, file={cif_file}")

    # ------------------------------------------------------------
    # Distance → contact
    # ------------------------------------------------------------
    mat = distance.cdist(V, V, metric="euclidean")
    mat = 1.0 / (mat + 1) ** (2 / 3)

    logger.info(f"Raw contact matrix: min={mat.min():.3e}, max={mat.max():.3e}, " f"mean={mat.mean():.3e}")

    if log_scale:
        mat = np.log1p(mat)
        logger.info("Applied log1p transform to contact matrix")

    # ------------------------------------------------------------
    # Optional reordering
    # ------------------------------------------------------------
    if reorder_by_diagonal:
        order = np.argsort(np.linalg.norm(V - np.mean(V, axis=0), axis=1))
        mat = mat[np.ix_(order, order)]
        logger.info("Reordered matrix by distance-from-centroid sorting")

    # ------------------------------------------------------------
    # Visualization
    # ------------------------------------------------------------
    if viz:

        from matplotlib.colors import PowerNorm

        fig, ax = plt.subplots(figsize=(8, 7), dpi=600)

        im = ax.imshow(mat, cmap="coolwarm", interpolation="nearest", norm=PowerNorm(gamma=0.4))

        ax.set_title("Structure-derived Contact Map", fontsize=12)
        ax.set_xlabel("Bead index")
        ax.set_ylabel("Bead index")

        cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("Contact strength (Inverse Distance)")

        ax.set_aspect("equal")
        ax.tick_params(length=0)

        # --------------------------------------------------------
        # Save (UNCHANGED)
        # --------------------------------------------------------
        if save and save_path is not None:
            logger.info(f"Saving heatmap to: {save_path}/{name}_contact_map.*")
            _save_local(fig, save_path + f"/{name}_contact_map")

        plt.close(fig)

    logger.info("Heatmap computation finished")
    return mat


def plot_md_thermo(history, save_path):
    """Plot energy + temperature evolution from MD."""
    logger.info("Creating MD thermodynamics plot...")

    steps = history["step"]

    fig, ax1 = plt.subplots(figsize=(6, 4))

    ax1.plot(steps, history["potential"], label="Potential energy")
    ax1.plot(steps, history["kinetic"], label="Kinetic energy")
    ax1.plot(steps, history["total"], label="Total energy")

    ax1.set_xlabel("Step")
    ax1.set_ylabel("Energy")
    ax1.legend()
    ax1.grid(True)

    ax2 = ax1.twinx()
    ax2.plot(steps, history["temperature"], "k--", label="Temperature")
    ax2.set_ylabel("Temperature")

    plt.title("MultiMM MD Thermodynamics")

    out = os.path.join(save_path, "plots/md_thermodynamics.png")
    plt.savefig(out, dpi=300)
    plt.close()

    logger.info(f"MD thermodynamics plot saved to: {out}")


def analyze_structure(V, save_path, name="structure"):
    """Advanced structural analysis for polymer-like 3D structures.

    Outputs:
    - detailed report (text)
    - multiple physically meaningful plots
    """
    # ------------------------------------------------------------
    # safety
    # ------------------------------------------------------------
    V = np.asarray(V, dtype=np.float64)
    V = V[np.isfinite(V).all(axis=1)]

    N = len(V)

    base = os.path.join(save_path, "analysis")
    os.makedirs(base, exist_ok=True)

    def _save_local(fig, fname):
        path = fname
        fig.savefig(path + ".png", dpi=300)
        fig.savefig(path + ".pdf")
        fig.savefig(path + ".svg")
        plt.close(fig)

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
    except:
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

    cos_angles = np.sum(v1 * v2, axis=1) / (np.linalg.norm(v1, axis=1) * np.linalg.norm(v2, axis=1) + 1e-8)

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
        chunk = V[i : i + window]
        cm = np.mean(chunk, axis=0)
        local_rg.append(np.sqrt(np.mean(np.sum((chunk - cm) ** 2, axis=1))))

    local_rg = np.array(local_rg)

    # REPORT
    report_path = os.path.join(base, f"{name}_report.txt")
    os.makedirs(base, exist_ok=True)

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

    # PLOTS
    os.makedirs(base + "/plots", exist_ok=True)

    # 1. bond lengths
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(bonds, bins=80)
    ax.set_title("Bond Length Distribution")
    ax.set_xlabel("Bond length")
    _save_local(fig, base + f"/plots/{name}_bonds")

    # ------------------------------------------------------------

    # 2. angles
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(angles, bins=80)
    ax.set_title("Angle Distribution")
    ax.set_xlabel("Angle (rad)")
    _save_local(fig, base + f"/plots/{name}_angles")

    # ------------------------------------------------------------

    # 3. distance vs genomic separation (VERY important)
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.plot(separations, spatial_dists)
    ax.set_title("Distance vs Genomic Separation")
    ax.set_xlabel("Genomic separation (beads)")
    ax.set_ylabel("Mean spatial distance")
    _save_local(fig, base + f"/plots/{name}_scaling")

    # ------------------------------------------------------------

    # 4. log-log scaling (polymer physics!!)
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.loglog(separations, spatial_dists)
    ax.set_title("Scaling (log-log)")
    ax.set_xlabel("s")
    ax.set_ylabel("R(s)")
    _save_local(fig, base + f"/plots/{name}_scaling_loglog")

    # ------------------------------------------------------------

    # 5. local compaction
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(local_rg)
    ax.set_title("Local Compaction (Sliding Rg)")
    ax.set_xlabel("Bead index")
    ax.set_ylabel("Local Rg")
    _save_local(fig, base + f"/plots/{name}_local_compaction")

    # ------------------------------------------------------------

    # 6. radial distribution
    r = np.linalg.norm(Vc, axis=1)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(r, bins=60)
    ax.set_title("Radial Distribution")
    ax.set_xlabel("Distance from COM")
    _save_local(fig, base + f"/plots/{name}_radial")

    # ------------------------------------------------------------
    return {
        "Rg": Rg,
        "Ree": Ree,
        "volume": volume,
        "density": density,
        "asphericity": asphericity,
        "acylindricity": acylindricity,
    }
