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
from mpl_toolkits.mplot3d import Axes3D
from .utils import get_coordinates_cif

logger = logging.getLogger(__name__)


pv.set_jupyter_backend("server")
color_dict = {-2: "#bf0020", -1: "#e36a24", 1: "#20c8e6", 2: "#181385", 0: "#ffffff"}
comp_dict = {-2: "B2", -1: "B1", 1: "A2", 2: "A1", 0: "no compartment"}

def plot_projection(struct_3D, Cs, save_path):
    """
    Enhanced structural analysis:
    - PCA projection
    - 3D structure
    - radial + component distributions
    - anisotropy proxy
    - PCA density map

    Removed: t-SNE (as requested)
    """

    sns.set_style("whitegrid")
    plt.rcParams.update({
        "figure.dpi": 600,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "font.size": 11
    })

    # ------------------------------------------------------------
    # Safety / preprocessing
    # ------------------------------------------------------------
    struct_3D = np.asarray(struct_3D, dtype=np.float64)
    Cs = np.asarray(Cs)

    N = struct_3D.shape[0]
    Cs = Cs[:N]

    mask_valid = np.isfinite(struct_3D).all(axis=1)
    struct_3D = struct_3D[mask_valid]
    Cs = Cs[mask_valid]

    # ------------------------------------------------------------
    # PCA
    # ------------------------------------------------------------
    pca = PCA(n_components=2)
    struct_2d = pca.fit_transform(struct_3D)

    # ------------------------------------------------------------
    # Geometry features (IMPORTANT for intuition)
    # ------------------------------------------------------------
    r = np.linalg.norm(struct_3D, axis=1)

    # anisotropy proxy (per-point spread along axes)
    anisotropy = np.abs(struct_3D[:, 0]) + np.abs(struct_3D[:, 1]) + np.abs(struct_3D[:, 2])

    # ------------------------------------------------------------
    # Dataframe
    # ------------------------------------------------------------
    df = pd.DataFrame({
        "x": struct_3D[:, 0],
        "y": struct_3D[:, 1],
        "z": struct_3D[:, 2],
        "x_pca": struct_2d[:, 0],
        "y_pca": struct_2d[:, 1],
        "distance": r,
        "anisotropy": anisotropy,
        "subcomp": Cs
    })

    df = df[df["subcomp"] != 0.0]

    # ------------------------------------------------------------
    # Output dir
    # ------------------------------------------------------------
    base_dir = os.path.join(save_path, "plots")
    os.makedirs(base_dir, exist_ok=True)

    def _save_local(fig, name):
        path = os.path.join(base_dir, name)
        fig.savefig(path + ".png", dpi=600)
        fig.savefig(path + ".pdf", dpi=600)
        fig.savefig(path + ".svg", dpi=600)

    # ============================================================
    # 1. PCA projection (with proper colorbar)
    # ============================================================
    fig, ax = plt.subplots(figsize=(7, 6))

    sc = ax.scatter(
        df["x_pca"],
        df["y_pca"],
        c=df["subcomp"],
        cmap="Spectral",
        s=12,
        alpha=0.6
    )

    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label("Subcompartment index")

    ax.set_title("PCA Projection")
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")

    _save_local(fig, "PCA_projection")
    plt.close(fig)

    # ============================================================
    # 2. 3D structure (main physical object)
    # ============================================================

    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot(111, projection="3d")

    sc = ax.scatter(
        df["x"],
        df["y"],
        df["z"],
        c=df["subcomp"],
        cmap="Spectral",
        s=4,
        alpha=0.7
    )

    cbar = fig.colorbar(sc, ax=ax, shrink=0.6)
    cbar.set_label("Subcompartment")

    ax.set_title("3D Chromatin Structure")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")

    _save_local(fig, "structure_3D")
    plt.close(fig)

    # ============================================================
    # 3. Radial distribution (compaction signal)
    # ============================================================
    fig, ax = plt.subplots(figsize=(7, 4))

    sns.kdeplot(
        data=df,
        x="distance",
        hue="subcomp",
        fill=True,
        palette="Spectral",
        alpha=0.5,
        ax=ax
    )

    ax.set_title("Radial Compaction Profile")
    ax.set_xlabel("Distance from origin")

    _save_local(fig, "radial_distribution")
    plt.close(fig)

    # ============================================================
    # 4. PCA density landscape (free energy-like intuition)
    # ============================================================
    fig, ax = plt.subplots(figsize=(7, 6))

    sc = sns.kdeplot(
        data=df,
        x="x_pca",
        y="y_pca",
        cmap="mako",
        fill=True,
        levels=50,
        thresh=0.05,
        ax=ax
    )

    ax.set_title("PCA Density Landscape")

    _save_local(fig, "pca_density")
    plt.close(fig)

    # ============================================================
    # 5. Component-wise spatial spread (NEW: very informative)
    # ============================================================
    fig, ax = plt.subplots(figsize=(7, 4))

    sns.boxplot(
        data=df,
        x="subcomp",
        y="distance",
        palette="Spectral",
        ax=ax
    )

    ax.set_title("Radial Distance by Subcompartment")
    ax.set_xlabel("Subcompartment")
    ax.set_ylabel("Radius")

    _save_local(fig, "radial_by_subcomp")
    plt.close(fig)

    # ============================================================
    # 6. Anisotropy distribution (NEW: shape diagnostics)
    # ============================================================
    fig, ax = plt.subplots(figsize=(7, 4))

    sns.kdeplot(
        data=df,
        x="anisotropy",
        hue="subcomp",
        fill=True,
        palette="Spectral",
        alpha=0.5,
        ax=ax
    )

    ax.set_title("Structural Anisotropy Distribution")

    _save_local(fig, "anisotropy_distribution")
    plt.close(fig)

    # ============================================================
    # 7. Axis projections (NEW: structure elongation intuition)
    # ============================================================
    fig, ax = plt.subplots(figsize=(7, 4))

    ax.scatter(df["x"], df["y"], s=3, alpha=0.4, label="XY")
    ax.scatter(df["x"], df["z"], s=3, alpha=0.4, label="XZ")
    ax.scatter(df["y"], df["z"], s=3, alpha=0.4, label="YZ")

    ax.set_title("Pairwise Coordinate Projections")
    ax.legend()

    _save_local(fig, "axis_projections")
    plt.close(fig)

    # ============================================================
    # NEW: PCA contour KDE per subcompartment (clean scientific plot)
    # ============================================================
    fig, ax = plt.subplots(figsize=(7, 6))

    subcomps = np.sort(df["subcomp"].unique())

    # consistent colormap
    cmap = plt.get_cmap("Spectral", len(subcomps))

    # grid for KDE evaluation
    x = df["x_pca"].values
    y = df["y_pca"].values

    xmin, xmax = x.min(), x.max()
    ymin, ymax = y.min(), y.max()

    X, Y = np.mgrid[
        xmin:xmax:200j,
        ymin:ymax:200j
    ]

    positions = np.vstack([X.ravel(), Y.ravel()])

    for i, sc_val in enumerate(subcomps):

        sub = df[df["subcomp"] == sc_val]

        if len(sub) < 10:
            continue  # avoid unstable KDE

        values = np.vstack([sub["x_pca"], sub["y_pca"]])

        kde = gaussian_kde(values, bw_method=0.2)
        Z = np.reshape(kde(positions).T, X.shape)

        # contour lines (clean + publication style)
        ax.contour(
            X, Y, Z,
            levels=5,
            colors=[cmap(i)],
            linewidths=1.2,
            alpha=0.9
        )

        # optional: faint fill for intuition
        ax.contourf(
            X, Y, Z,
            levels=3,
            alpha=0.08,
            colors=[cmap(i)]
        )

    ax.set_title("Subcompartment KDE Contours (PCA space)")
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    
    # clean frame
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    _save_local(fig, "pca_kde_contours")
    plt.close(fig)

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


def viz_structure(V, colors=None, r=0.1, cmap="coolwarm", save_path=None):
    """
    Visualize structure V and optionally save it to a file.
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

        colors_min = np.min(colors)
        colors_max = np.max(colors)
        diff = colors_max - colors_min

        logger.info(
            f"Color mapping enabled: min={colors_min}, max={colors_max}, diff={diff}"
        )

        if diff > 0:
            color_values = (colors - colors_min) / diff
        else:
            color_values = np.zeros_like(colors, dtype=float)
            logger.warning("Flat color array detected (all values identical)")

        polyline["colors"] = color_values
        polymer = polyline.tube(radius=r)

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
    name="structure"
):
    """
    Compute and visualize contact/interaction heatmap from 3D structure.
    """

    # ------------------------------------------------------------
    # Output dir (UNCHANGED)
    # ------------------------------------------------------------
    base_dir = save_path
    os.makedirs(base_dir, exist_ok=True)

    def _save_local(fig, name):
        path = name
        fig.savefig(path + ".png", dpi=300, bbox_inches="tight")
        fig.savefig(path + ".pdf", bbox_inches="tight")
        fig.savefig(path + ".svg", bbox_inches="tight")

    # ------------------------------------------------------------
    # Load structure
    # ------------------------------------------------------------
    V = get_coordinates_cif(cif_file)
    logger.info(f"Loaded structure: shape={V.shape}, file={cif_file}")

    # ------------------------------------------------------------
    # Distance → contact
    # ------------------------------------------------------------
    mat = distance.cdist(V, V, metric="euclidean")
    mat = 1.0 / (mat + 1)**(2/3)

    logger.info(
        f"Raw contact matrix: min={mat.min():.3e}, max={mat.max():.3e}, "
        f"mean={mat.mean():.3e}"
    )

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

        im = ax.imshow(
            mat,
            cmap="coolwarm",
            interpolation="nearest",
            norm=PowerNorm(gamma=0.4)
        )

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

def analyze_structure(V, save_path, name="structure"):
    """
    Advanced structural analysis for polymer-like 3D structures.

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
        local_rg.append(np.sqrt(np.mean(np.sum((chunk - cm)**2, axis=1))))

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
    os.makedirs(base+'/plots', exist_ok=True)

    # 1. bond lengths
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(bonds, bins=80)
    ax.set_title("Bond Length Distribution")
    ax.set_xlabel("Bond length")
    _save_local(fig, base+f"/plots/{name}_bonds")

    # ------------------------------------------------------------

    # 2. angles
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(angles, bins=80)
    ax.set_title("Angle Distribution")
    ax.set_xlabel("Angle (rad)")
    _save_local(fig, base+f"/plots/{name}_angles")

    # ------------------------------------------------------------

    # 3. distance vs genomic separation (VERY important)
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.plot(separations, spatial_dists)
    ax.set_title("Distance vs Genomic Separation")
    ax.set_xlabel("Genomic separation (beads)")
    ax.set_ylabel("Mean spatial distance")
    _save_local(fig, base+f"/plots/{name}_scaling")

    # ------------------------------------------------------------

    # 4. log-log scaling (polymer physics!!)
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.loglog(separations, spatial_dists)
    ax.set_title("Scaling (log-log)")
    ax.set_xlabel("s")
    ax.set_ylabel("R(s)")
    _save_local(fig, base+f"/plots/{name}_scaling_loglog")

    # ------------------------------------------------------------

    # 5. local compaction
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(local_rg)
    ax.set_title("Local Compaction (Sliding Rg)")
    ax.set_xlabel("Bead index")
    ax.set_ylabel("Local Rg")
    _save_local(fig, base+f"/plots/{name}_local_compaction")

    # ------------------------------------------------------------

    # 6. radial distribution
    r = np.linalg.norm(Vc, axis=1)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(r, bins=60)
    ax.set_title("Radial Distribution")
    ax.set_xlabel("Distance from COM")
    _save_local(fig, base+f"/plots/{name}_radial")

    # ------------------------------------------------------------
    return {
        "Rg": Rg,
        "Ree": Ree,
        "volume": volume,
        "density": density,
        "asphericity": asphericity,
        "acylindricity": acylindricity,
    }