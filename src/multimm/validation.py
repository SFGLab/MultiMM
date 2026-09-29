import logging

from .logger import log_table

import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as ndimage
import seaborn as sns
from scipy.spatial import KDTree
from scipy.stats import pearsonr, spearmanr
from sklearn.decomposition import PCA
from tqdm import tqdm

from .utils import (
    get_coordinates_cif,
    min_max_normalize,
    standarize,
    structure_to_heatmap,
    mean_downsample,
    rescale_matrix,
    remove_zero_rows_and_columns,
    remove_diagonals,
    compute_compartments,
)
from .read_hic import pool_to_n_beads as _pool_to_n_beads

logger = logging.getLogger(__name__)


# Metrics for comparisons
def calculate_correlation(matrix1, matrix2):
    # Flatten matrices and calculate Pearson correlation
    flat1 = matrix1.flatten()
    flat2 = matrix2.flatten()
    correlation, p_val = pearsonr(flat1, flat2)
    return correlation, p_val


def rv_coefficient(matrix1, matrix2):
    """Computes the RV coefficient between two matrices.

    Parameters:
    matrix1 (ndarray): First input matrix.
    matrix2 (ndarray): Second input matrix, must have the same shape as matrix1.

    Returns:
    float: The RV coefficient, a measure of similarity between the matrices.
    """
    if matrix1.shape != matrix2.shape:
        raise ValueError("Matrices must have the same dimensions.")

    # Compute inner products
    matrix1_flat = matrix1.flatten()
    matrix2_flat = matrix2.flatten()

    # Calculate the dot products
    dot_product_12 = np.dot(matrix1_flat, matrix2_flat)
    dot_product_11 = np.dot(matrix1_flat, matrix1_flat)
    dot_product_22 = np.dot(matrix2_flat, matrix2_flat)

    # Compute the RV coefficient
    rv = dot_product_12 / np.sqrt(dot_product_11 * dot_product_22)

    return rv


def mantel_test(matrix1, matrix2, permutations=1000):
    """Computes the Mantel test correlation between two distance matrices.

    Parameters:
    matrix1 (ndarray): First distance matrix.
    matrix2 (ndarray): Second distance matrix, must have the same shape as matrix1.
    permutations (int): Number of permutations for significance testing.

    Returns:
    tuple: Mantel correlation coefficient and p-value.
    """
    if matrix1.shape != matrix2.shape:
        raise ValueError("Matrices must have the same dimensions.")

    # Flatten upper triangular parts of the matrices to avoid redundancy
    triu_indices = np.triu_indices_from(matrix1, k=1)
    flat_matrix1 = matrix1[triu_indices]
    flat_matrix2 = matrix2[triu_indices]

    # Compute the Pearson correlation for the original matrices
    mantel_corr, _ = pearsonr(flat_matrix1, flat_matrix2)

    # Permutation testing
    permuted_corrs = []
    for _ in range(permutations):
        np.random.shuffle(flat_matrix2)
        permuted_corr, _ = pearsonr(flat_matrix1, flat_matrix2)
        permuted_corrs.append(permuted_corr)

    # Calculate p-value: proportion of permuted correlations greater than or equal to the observed
    permuted_corrs = np.array(permuted_corrs)
    p_value = np.sum(permuted_corrs >= mantel_corr) / permutations

    return mantel_corr, p_value


# Metrics with moving windows
def fast_pearson_correlation(m1, m2):
    """Compute Pearson correlation efficiently between two flattened arrays."""
    m1_flat = m1.flatten()
    m2_flat = m2.flatten()

    # Compute covariance and standard deviations
    cov = np.mean((m1_flat - m1_flat.mean()) * (m2_flat - m2_flat.mean()))
    std_m1 = np.std(m1_flat)
    std_m2 = np.std(m2_flat)

    return cov / (std_m1 * std_m2)


def compute_pearson_correlation(m1, m2, window_size):
    """Compute the average Pearson correlation for non-overlapping windows of a
    given size."""
    N = m1.shape[0]
    correlations = []

    # Iterate over non-overlapping windows
    for i in range(0, N, window_size):
        for j in range(0, N, window_size):
            # Make sure the window fits completely within the matrix
            if i + window_size <= N and j + window_size <= N:
                sub_m1 = m1[i : i + window_size, j : j + window_size]
                sub_m2 = m2[i : i + window_size, j : j + window_size]

                # Compute Pearson correlation for the submatrices
                correlation = fast_pearson_correlation(sub_m1, sub_m2)
                correlations.append(correlation)

    return np.mean(correlations)


def correlation_vs_window_size(m1, m2):
    """Compute the average Pearson correlation for selected window sizes and
    plot the result."""
    N = m1.shape[0]

    # Select window sizes that divide N evenly
    window_sizes = np.arange(100, 3 * N // 4, 100)
    avg_correlations = []

    for window_size in tqdm(window_sizes):
        avg_correlation = compute_pearson_correlation(m1, m2, window_size)
        avg_correlations.append(avg_correlation)

    # Plotting the results
    plt.plot(window_sizes, avg_correlations, "k-o")
    plt.xlabel("Window Size")
    plt.ylabel("Average Pearson Correlation")
    plt.grid(True)
    plt.show()
    return window_sizes, avg_correlations


# Random Walks
def random_walk_3d(N):
    """Generate a random walk with N points in 3D.

    Parameters:
    N (int): Number of points in the random walk

    Returns:
    np.ndarray: Array of shape (N, 3) with coordinates of the walk
    """
    if N < 1:
        raise ValueError("N must be at least 1")

    # Directions for movement in 3D: x, y, z
    directions = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [-1, 0, 0], [0, -1, 0], [0, 0, -1]])

    # Initialize the walk
    walk = np.zeros((N, 3), dtype=int)

    current_position = np.array([0, 0, 0])

    for i in range(1, N):
        direction = directions[np.random.choice(len(directions))]
        new_position = current_position + direction
        walk[i] = new_position
        current_position = new_position

    return walk


def generate_self_avoiding_walk(N, step_size=1.0, max_backtracks=50):
    """Generates a self-avoiding random walk in 3D with N points, where
    consecutive points have a constant distance (step_size) and beads avoid
    each other, using a backtracking approach.

    Parameters:
    N (int): Number of points in the structure.
    step_size (float): Constant distance between consecutive points.
    max_backtracks (int): Maximum number of steps to backtrack if stuck.

    Returns:
    walk (ndarray): The generated 3D structure with shape (N, 3).
    """
    # Initialize the walk with the first point at the origin
    walk = np.zeros((N, 3))

    # Directions for movement: 6 possible unit steps in 3D
    directions = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]])

    current_index = 1
    backtracks = 0

    while current_index < N:
        # Generate a list of possible next positions
        np.random.shuffle(directions)  # Shuffle directions to add randomness
        placed = False

        for direction in directions:
            # Calculate the next point position
            next_point = walk[current_index - 1] + direction * step_size

            # Check for self-avoidance: point should not overlap with any existing points
            if not np.any(np.all(np.isclose(walk[:current_index], next_point, atol=1e-6), axis=1)):
                # Place the next point if it's valid
                walk[current_index] = next_point
                current_index += 1
                placed = True
                backtracks = 0  # Reset backtracks counter after a successful placement
                break

        # If no valid placement was found, backtrack
        if not placed:
            current_index -= 1
            backtracks += 1

            # If the number of backtracks exceeds the limit, stop to avoid infinite loops
            if backtracks > max_backtracks:
                logger.info(f"Exceeded maximum backtracks. Unable to complete the walk with {N} points.")
                return walk[:current_index]  # Return the walk up to the last valid point

    return walk


def pca_downsample(V, n):
    pca = PCA(n_components=3)
    V_reduced = pca.fit_transform(V)
    indices = np.linspace(0, V_reduced.shape[0] - 1, n, dtype=int)
    V_downgraded = V_reduced[indices, :]
    return V_downgraded


# Heatmap Comparison
def find_local_maxima(heatmap, min_distance=1):
    # Find local maxima
    maxima = ndimage.maximum_filter(heatmap, size=min_distance) == heatmap
    maxima_positions = np.transpose(np.nonzero(maxima))  # (N, 2) array of positions
    maxima_intensities = heatmap[maxima_positions[:, 0], maxima_positions[:, 1]]

    return maxima_positions, maxima_intensities


def compare_maxima_positions(pos1, pos2, distance_threshold=1):
    # Create KDTree for the second set of positions
    tree = KDTree(pos2)

    matched_pairs = []  # To store pairs of matched indices

    for idx, point in enumerate(pos1):
        # Find the nearest point in pos2 within the distance threshold
        distances, indices = tree.query(point, distance_upper_bound=distance_threshold)
        if distances != np.inf:  # Valid match within threshold
            matched_pairs.append((idx, indices))

    return matched_pairs


def analyze_heatmaps(heatmap1, heatmap2, min_distance=1, distance_threshold=1):
    # Step 1: Find local maxima for both heatmaps
    pos1, intensities1 = find_local_maxima(heatmap1, min_distance=min_distance)
    pos2, intensities2 = find_local_maxima(heatmap2, min_distance=min_distance)

    # Step 2: Compare positions of maxima and find matches
    matched_pairs = compare_maxima_positions(pos1, pos2, distance_threshold=distance_threshold)

    # Extract matched intensities
    matched_intensities1 = np.array([intensities1[i1] for i1, _ in matched_pairs])
    matched_intensities2 = np.array([intensities2[i2] for _, i2 in matched_pairs])

    # Step 3: Calculate percentage of common maxima
    percentage_common_maxima = (len(matched_pairs) / len(pos1)) * 100

    # Step 4: Compute correlation of intensities (Pearson correlation)
    if len(matched_intensities1) > 1 and len(matched_intensities2) > 1:
        correlation, _ = pearsonr(matched_intensities1, matched_intensities2)
    else:
        correlation = np.nan  # In case of insufficient data

    return percentage_common_maxima, correlation


def compare_matrices(m, mr, exp_m, viz=True):
    # Remove empty rows
    exp_m, rs, cs = remove_zero_rows_and_columns(exp_m)
    m = np.delete(m, rs, axis=0)
    m = np.delete(m, cs, axis=1)
    mr = np.delete(mr, rs, axis=0)
    mr = np.delete(mr, cs, axis=1)

    # Normalize
    m, mr, exp_m = standarize(m), standarize(mr), standarize(exp_m)
    m, mr, exp_m = min_max_normalize(m), min_max_normalize(mr), min_max_normalize(exp_m)

    # Compute Compartments
    eignvec_sim, _ = compute_compartments(m)
    eignvec_rw, _ = compute_compartments(mr)
    eignvec_exp1, eignvec_exp2 = compute_compartments(exp_m)

    # Compute correlation with the experimental eigenvector
    corr_sim1, pvalue_sim1 = pearsonr(eignvec_sim, eignvec_exp1)
    corr_rw1, pvalue_rw1 = pearsonr(eignvec_rw, eignvec_exp1)
    corr_sim2, pvalue_sim2 = pearsonr(eignvec_sim, eignvec_exp2)
    corr_rw2, pvalue_rw2 = pearsonr(eignvec_rw, eignvec_exp2)

    if viz:
        logger.info("Correlation of simulation with the first eigenvector:", corr_sim1)
        logger.info("Correlation of random walk with the first eigenvector:", corr_rw1)
        logger.info("Correlation of simulation with the second eigenvector:", corr_sim2)
        logger.info("Correlation of random walk with the second eigenvector:", corr_rw2)

    return np.abs(corr_sim1), np.abs(corr_rw1), np.abs(corr_sim2), np.abs(corr_rw2)


def pipeline_single_ensemble(V, Vr, exp_m, viz=True):
    # Downgrade structures to have same number of points as the experimental heatmaps
    V_reduced = mean_downsample(V, len(exp_m))
    Vr_reduced = mean_downsample(Vr, len(exp_m))

    # Compute simulated heatmaps
    m = structure_to_heatmap(V_reduced)
    mr = structure_to_heatmap(Vr_reduced)

    # Compute the statistics
    corr_sim1, corr_rw1, corr_sim2, corr_rw2 = compare_matrices(m, mr, exp_m, viz)


def ensemble_pipeline_boxplot(ensemble_path, exp_path, N_chroms=22, N_ens=20, viz=True):
    # Lists to hold correlation values
    Cs_sim, Cs_rw = list(), list()

    for i in range(N_chroms):
        # Load experimental data
        exp_m = np.nan_to_num(np.load(exp_path + f"/chrom{i+1}_primary+replicate_ice_norm_50kb.npy"))
        exp_m = remove_diagonals(exp_m, 5)
        L = len(exp_m)

        # Calculate heatmaps from experimental structures
        logger.info(f"Chromosome {i+1}:")
        logger.info("Calculating heatmaps from experimental structures...")
        corrs_sim, corrs_rw = [], []
        for j in tqdm(range(N_ens)):
            V = get_coordinates_cif(ensemble_path + f"/ens_{j+1}/chromosomes/MultiMM_minimized_chr{i+1}.cif")
            V_reduced = mean_downsample(V, L)
            m = structure_to_heatmap(V_reduced)

            Vr = random_walk_3d(len(V))
            Vr_reduced = mean_downsample(Vr, L)
            mr = structure_to_heatmap(Vr_reduced)

            # Calculate correlations
            corr_sim1, corr_rw1, _, _ = compare_matrices(m, mr, exp_m, False)
            corrs_sim.append(np.abs(corr_sim1))
            corrs_rw.append(np.abs(corr_rw1))

        Cs_sim.append(corrs_sim)
        Cs_rw.append(corrs_rw)
        logger.info(f"Chromosome {i+1} done!\n")

    if viz:
        # Box plot for each chromosome
        plt.figure(figsize=(20, 5), dpi=200)

        # Prepare data for box plots
        data_sim = [Cs_sim[i] for i in range(N_chroms)]
        data_rw = [Cs_rw[i] for i in range(N_chroms)]

        box_sim = plt.boxplot(
            data_sim,
            positions=np.arange(N_chroms) - 0.2,
            widths=0.4,
            patch_artist=True,
            boxprops=dict(facecolor="blue", color="blue"),
            medianprops=dict(color="black"),
        )
        box_rw = plt.boxplot(
            data_rw,
            positions=np.arange(N_chroms) + 0.2,
            widths=0.4,
            patch_artist=True,
            boxprops=dict(facecolor="red", color="red"),
            medianprops=dict(color="black"),
        )

        plt.xticks(np.arange(N_chroms), [f"chr{i + 1}" for i in range(N_chroms)])
        plt.xlabel("Chromosomes", fontsize=16)
        plt.ylabel("Correlation with 1st Eigenvector", fontsize=14)
        plt.legend(
            [box_sim["boxes"][0], box_rw["boxes"][0]],
            ["Simulation", "Random Walk"],
            loc="upper right",
        )
        plt.savefig("heatmap_correlation_boxplots.pdf", format="pdf", dpi=200)
        plt.savefig("heatmap_correlation_boxplots.svg", format="svg", dpi=200)
        plt.show()


def ensemble_pipeline_bars(ensemble_path, exp_path, N_chroms=22, N_ens=20, viz=True):
    # Run loop for each chromosome
    Cs_sim1, Cs_sim2 = list(), list()
    Cs_rw1, Cs_rw2 = list(), list()

    for i in range(N_chroms):
        # Experimental path
        exp_m = np.nan_to_num(np.load(exp_path + f"/chrom{i+1}_primary+replicate_ice_norm_50kb.npy"))
        L = len(exp_m)

        # Average the heatmaps of each ensemble
        avg_m = 0
        logger.info(f"Chromosome {i+1}:")
        logger.info("Calculating heatmaps from experimental structures...")
        for j in tqdm(range(N_ens)):
            V = get_coordinates_cif(ensemble_path + f"/ens_{j+1}/chromosomes/MultiMM_minimized_chr{i+1}.cif")
            V_reduced = mean_downsample(V, L)
            m = structure_to_heatmap(V_reduced)
            avg_m += m
        avg_m /= N_ens

        # Average the heatmaps of each random walk
        avg_mr = 0
        logger.info("Calculating heatmaps from random structures...")
        for j in tqdm(range(N_ens)):
            Vr = random_walk_3d(len(V))
            Vr_reduced = mean_downsample(Vr, L)
            mr = structure_to_heatmap(Vr_reduced)
            avg_mr += mr
        avg_mr /= N_ens

        avg_m = remove_diagonals(avg_m, 1)
        avg_mr = remove_diagonals(avg_mr, 1)
        exp_m = remove_diagonals(exp_m, 1)

        # Compare them
        corr_sim1, corr_rw1, corr_sim2, corr_rw2 = compare_matrices(avg_m, avg_mr, exp_m, False)
        Cs_sim1.append(corr_sim1)
        Cs_sim2.append(corr_sim2)
        Cs_rw1.append(corr_rw1)
        Cs_rw2.append(corr_rw2)

        if viz:
            logger.info("Correlation of simulation with the first eigenvector:", corr_sim1)
            logger.info("Correlation of random walk with the first eigenvector:", corr_rw1)
            logger.info("Correlation of simulation with the second eigenvector:", corr_sim2)
            logger.info("Correlation of random walk with the second eigenvector:", corr_rw2)
        logger.info(f"Chromosome {i+1} done!\n")

    if viz:
        chroms = ["chr" + str(i) for i in range(1, N_chroms + 1)]
        X_axis = np.arange(N_chroms)
        plt.figure(figsize=(20, 5), dpi=200)
        plt.bar(X_axis - 0.2, Cs_sim1, 0.4, label="Simulation", color="blue")
        plt.bar(X_axis + 0.2, Cs_rw1, 0.4, label="Random Walk", color="red")
        plt.xticks(X_axis, chroms)
        plt.xlabel("Chromosomes", fontsize=16)
        plt.legend()
        plt.ylabel("Correlation with First Eigenvector", fontsize=14)
        plt.savefig("corr_1st_eigenvec.pdf", format="pdf", dpi=200)
        plt.savefig("corr_1st_eigenvec.svg", format="svg", dpi=200)

        plt.show()

        plt.figure(figsize=(20, 5), dpi=200)
        plt.bar(X_axis - 0.2, Cs_sim2, 0.4, label="Simulation", color="blue")
        plt.bar(X_axis + 0.2, Cs_rw2, 0.4, label="Random Walk", color="red")

        plt.xticks(X_axis, chroms)
        plt.xlabel("Chromosomes", fontsize=16)
        plt.legend()
        plt.ylabel("Correlation with Second Eigenvector", fontsize=14)
        plt.savefig("corr_2st_eigenvec.pdf", format="pdf", dpi=200)
        plt.savefig("corr_2st_eigenvec.svg", format="svg", dpi=200)
        plt.show()


def regions_pipeline(regions_dir, chroms, starts, ends, N_ens=1000):
    # Experimental path
    corrs_sim, corrs_rw = [], []
    pvals_sim, pvals_rw = [], []
    ps_sim, ints_sim = [], []
    ps_rw, ints_rw = [], []

    # Average the heatmaps of each ensemble
    logger.info("Calculating heatmaps from experimental structures...")
    for i in tqdm(range(N_ens)):
        try:
            exp_m = np.nan_to_num(
                np.load(
                    f"/home/skorsak/Data/Rao/regs/primary+replicate_ice_norm_25kb_chrom{chroms[i]}_{starts[i]}_{ends[i]}.npy"
                )
            )
            L = len(exp_m)
        except Exception:
            logger.info(f"Problem with chromosome {chroms[i]}_{starts[i]}_{ends[i]} experimental data")
            continue

        try:
            V = get_coordinates_cif(
                regions_dir + f"/regens_chrom{chroms[i]}_region_{starts[i]}_{ends[i]}/MultiMM_minimized.cif"
            )
        except Exception:
            logger.info(f"Problem with chromosome {chroms[i]}_{starts[i]}_{ends[i]} simulated data")
            continue
        V_reduced = mean_downsample(V, L)
        m = structure_to_heatmap(V_reduced)

        Vr = random_walk_3d(len(V))
        Vr_reduced = mean_downsample(Vr, L)
        mr = structure_to_heatmap(Vr_reduced)

        exp_m, rs, cs = remove_zero_rows_and_columns(exp_m)
        m = np.delete(m, rs, axis=0)
        m = np.delete(m, cs, axis=1)
        mr = np.delete(mr, rs, axis=0)
        mr = np.delete(mr, cs, axis=1)

        m = remove_diagonals(m, 1)
        mr = remove_diagonals(mr, 1)
        exp_m = remove_diagonals(exp_m, 1)

        m = (m - np.mean(m)) / np.std(m)
        mr = (mr - np.mean(mr)) / np.std(mr)
        exp_m = (exp_m - np.mean(exp_m)) / np.std(exp_m)
        m = (m - np.min(m)) / (np.max(m) - np.min(m))
        mr = (mr - np.min(mr)) / (np.max(mr) - np.min(mr))
        exp_m = (exp_m - np.min(exp_m)) / (np.max(exp_m) - np.min(exp_m))

        p_sim, int_sim = analyze_heatmaps(
            remove_diagonals(m, 4),
            remove_diagonals(exp_m, 4),
            min_distance=5,
            distance_threshold=5,
        )
        p_rw, int_rw = analyze_heatmaps(
            remove_diagonals(mr, 4),
            remove_diagonals(exp_m, 4),
            min_distance=5,
            distance_threshold=5,
        )

        corr_V, pval_sim = calculate_correlation(m, exp_m)
        corr_Vr, pval_rw = calculate_correlation(mr, exp_m)
        corrs_sim.append(corr_V)
        corrs_rw.append(corr_Vr)
        pvals_sim.append(pval_sim)
        pvals_rw.append(pval_rw)
        ps_sim.append(p_sim)
        ps_rw.append(p_rw)
        ints_sim.append(int_sim)
        ints_rw.append(int_rw)

    # Create the violin plot
    data = [corrs_sim, corrs_rw]
    plt.figure(figsize=(6, 9))
    sns.violinplot(data=data)
    plt.xticks([0, 1], ["Simulation", "Random Walk"], fontsize=16)
    plt.ylabel("Correlation with Experimental Data", fontsize=16)
    plt.savefig("violin.svg", format="svg", dpi=200)
    plt.savefig("violin.pdf", format="pdf", dpi=200)
    plt.show()

    # Create the violin plot
    data = [np.array(ps_sim) / 100, np.array(ps_rw) / 100]
    plt.figure(figsize=(6, 9))
    sns.violinplot(data=data)
    plt.xticks([0, 1], ["Simulation", "Random Walk"], fontsize=16)
    plt.ylabel("Percentage of Common Loops", fontsize=16)
    plt.savefig("violin_ps.pdf", format="pdf", dpi=200)
    plt.show()

    # Create the violin plot
    data = [ints_sim, ints_rw]
    plt.figure(figsize=(6, 9))
    sns.violinplot(data=data)
    plt.xticks([0, 1], ["Simulation", "Random Walk"], fontsize=16)
    plt.ylabel("Peak Intensity Correlation", fontsize=16)
    plt.savefig("violin_ints.pdf", format="pdf", dpi=200)
    plt.show()


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


# ── Random-walk baseline ─────────────────────────────────────────────────────

def _rw_baseline_contact(
    N: int,
    n_rw: int = 20,
    step_nm: float = 0.1,
    eps: float = 1e-3,
    seed: int = 0,
    confine_radius_nm: float | None = None,
) -> np.ndarray:
    """Generate a null-model (random-walk ensemble) average contact map.

    Builds *n_rw* simple 3-D random walks of length *N* (bond length
    *step_nm* nm) optionally confined inside a sphere of radius
    *confine_radius_nm*, computes the 1/(d+ε) contact proxy for each,
    and returns the average as a float64 N×N matrix.

    Confinement is implemented as elastic reflection: at each step, if the
    bead would land outside the sphere, it is reflected back along the
    radial direction.  This gives a realistic confined-polymer null model
    comparable to the nuclear boundary used in the actual MD simulation.

    Without confinement an N=1000 chain with step_nm=0.1 has end-to-end
    ≈ 3.16 nm, far larger than a typical nucleus (~1 nm at simulation
    scale), so the unconfined walk produces artificially high diagonal
    correlations simply because beads at similar genomic position are
    close in 3-D by construction (not by folding) — inflating the RW
    baseline.

    Parameters
    ----------
    N                  : number of beads
    n_rw               : number of independent realisations to average
                         (default 20 — enough for stable mean)
    step_nm            : bond length in nm (match POL_HARMONIC_BOND_R0)
    eps                : regularisation in 1/(d+ε)
    seed               : numpy random seed for reproducibility
    confine_radius_nm  : sphere radius in nm.  None → unconfined (legacy).
                         Recommended: set to the simulation nuclear radius
                         (radius2 attribute of the Model object).
    """
    rng = np.random.default_rng(seed)
    acc: np.ndarray | None = None
    for _ in range(n_rw):
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

        frame = _coords_to_inv_contact_f32(coords, eps=eps)
        if acc is None:
            acc = frame.astype(np.float64)
        else:
            acc += frame.astype(np.float64)
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


def _accumulate_ensemble(
    cif_paths: list,
    eps: float = 1e-3,
    log=None,
) -> np.ndarray:
    """Streaming accumulation of the inverse-contact ensemble average.

    Memory cost: one float32 N×N accumulator  +  one float32 N×N working buffer
    (both reused across frames — no extra allocations inside the loop).

    Returns
    -------
    inv_avg : (N, N) float64  (upcast at the end for downstream precision)
    """
    _log = log or logger
    _log.info("Streaming inverse-contact accumulation over %d frames …", len(cif_paths))

    acc: np.ndarray | None = None

    for path in tqdm(cif_paths, desc="Frames", leave=False):
        coords = get_coordinates_cif(path)          # (N, 3) float64
        frame  = _coords_to_inv_contact_f32(coords, eps=eps)   # (N, N) float32, in-place
        if acc is None:
            acc = frame                              # first frame — reuse buffer
        else:
            acc += frame                             # in-place accumulation

    if acc is None:
        raise ValueError("No frames were loaded.")

    acc /= len(cif_paths)                           # normalise in-place
    return acc.astype(np.float64, copy=False)       # upcast once


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

def validate_hic_model(
    cif_path: str,
    hic_matrix: np.ndarray,
    insulation_window: int = 10,
    max_diag: int | None = None,
    save_path: str | None = None,
    n_rw: int = 20,
    rw_step_nm: float = 0.1,
    confine_radius_nm: float | None = None,
    log=None,
) -> dict:
    """Validate a single simulated structure against experimental Hi-C.

    Three metrics are computed between the simulated contact proxy
    (1/(distance + ε)) and the experimental Hi-C matrix:

    1. **Diagonal decay correlation** — Pearson r between per-diagonal mean
       contact profiles.  Captures whether the overall distance-decay law of
       contacts is reproduced.
    2. **Insulation score correlation** — Pearson r between the sliding-window
       insulation score vectors.  Measures TAD boundary agreement.
    3. **PC1 correlation** — Pearson r between the first principal component of
       the O/E-normalised contact matrices.  Reflects A/B compartmentalisation.

    Parameters
    ----------
    cif_path : str
        Path to a MultiMM output .cif file.
    hic_matrix : ndarray, shape (M, M)
        Experimental Hi-C contact matrix (raw counts or normalised).
    insulation_window : int
        Half-width of the insulation-score sliding window (in beads).
    max_diag : int or None
        How many diagonals to include in the decay profile.  Defaults to N//2.
    log : logging.Logger, optional

    Returns
    -------
    dict with keys:
        'diag_decay_r', 'diag_decay_p',
        'insulation_r', 'insulation_p',
        'pc1_r',        'pc1_p',
        'pearson_r',    'pearson_p',
        'spearman_r',   'spearman_p',
        'ssim',         'gmsd',         'nmi'
    """
    _log = log or logger

    # Use the fast Gram-matrix accumulator (float32, in-place, no temp arrays)
    coords      = get_coordinates_cif(cif_path)
    sim_contact = _coords_to_inv_contact_f32(coords).astype(np.float64, copy=False)
    N           = sim_contact.shape[0]

    # Resize experimental matrix to match simulation resolution
    hic_r = _pool_matrix(hic_matrix, N)
    if max_diag is None:
        max_diag = N // 2

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
            N, n_rw=n_rw, step_nm=rw_step_nm,
            confine_radius_nm=confine_radius_nm,
        )
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

    if save_path is not None:
        from .plots import plot_hic_comparison
        import os
        plots_dir = os.path.join(save_path, "plots")
        os.makedirs(plots_dir, exist_ok=True)
        plot_hic_comparison(
            sim_contact, hic_r, plots_dir,
            name="hic_comparison_single",
            rw_matrix=rw_contact if rw_ok else None,
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
    eps: float = 1e-3,
    save_path: str | None = None,
    n_rw: int = 20,
    rw_step_nm: float = 0.1,
    confine_radius_nm: float | None = None,
    log=None,
) -> dict:
    """Validate an ensemble of simulated structures against experimental Hi-C.

    The simulated contact proxy is the **inverse-average** of distance maps
    across all trajectory frames:

        C_sim = mean_k[ 1 / (D_k + ε) ]

    This is the standard ensemble-to-contact conversion for polymer simulations.
    Three metrics are then computed identically to ``validate_hic_model``.

    Parameters
    ----------
    cif_paths : list of str
        Paths to per-frame CIF files from the MD trajectory.
    hic_matrix : ndarray, shape (M, M)
        Experimental Hi-C contact matrix.
    insulation_window : int
        Half-width of the insulation-score sliding window.
    max_diag : int or None
        Number of diagonals in the decay profile.
    eps : float
        Regularisation for 1/(d + ε).
    log : logging.Logger, optional

    Returns
    -------
    dict with keys:
        'diag_decay_r', 'diag_decay_p',
        'insulation_r', 'insulation_p',
        'pc1_r',        'pc1_p',
        'pearson_r',    'pearson_p',
        'spearman_r',   'spearman_p',
        'ssim',         'gmsd',         'nmi'
    """
    _log = log or logger

    # ── Fast streaming accumulation (float32, in-place, no temp arrays) ───────
    # Uses the Gram-matrix identity D²_ij = ||r_i||² + ||r_j||² − 2 r_i·r_j
    # so no separate N×N distance buffer is ever materialised.
    # Peak RAM = 2 × N² × 4 bytes (accumulator + one working frame).
    inv_avg = _accumulate_ensemble(cif_paths, eps=eps, log=_log)

    N     = inv_avg.shape[0]
    hic_r = _pool_matrix(hic_matrix, N)
    if max_diag is None:
        max_diag = N // 2

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
            N, n_rw=n_rw, step_nm=rw_step_nm, eps=eps,
            confine_radius_nm=confine_radius_nm,
        )
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

    if save_path is not None:
        from .plots import plot_hic_comparison
        import os
        plots_dir = os.path.join(save_path, "plots")
        os.makedirs(plots_dir, exist_ok=True)
        plot_hic_comparison(
            inv_avg, hic_r, plots_dir,
            name="hic_comparison_ensemble",
            rw_matrix=rw_contact if rw_ok else None,
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
