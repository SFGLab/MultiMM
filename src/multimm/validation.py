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
    """Down-sample or up-sample *m* to an N×N matrix by mean-pooling."""
    M = m.shape[0]
    if M == N:
        return m.astype(float)
    if M > N:
        # mean-pool: reshape into (N, block, N, block)
        block = M // N
        trimmed = m[: N * block, : N * block]
        return trimmed.reshape(N, block, N, block).mean(axis=(1, 3))
    # upsample (simple repeat — rare in practice)
    factor = N // M
    return np.repeat(np.repeat(m, factor, axis=0), factor, axis=1).astype(float)


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


# ── Public validation API ─────────────────────────────────────────────────────

def validate_hic_model(
    cif_path: str,
    hic_matrix: np.ndarray,
    insulation_window: int = 10,
    max_diag: int | None = None,
    save_path: str | None = None,
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
        'pc1_r',        'pc1_p'
    """
    from .utils import get_coordinates_cif, model_distance_heatmap
    _log = log or logger

    coords   = get_coordinates_cif(cif_path)
    dist_map = model_distance_heatmap(coords)
    N        = dist_map.shape[0]

    # Resize experimental matrix to match simulation resolution
    hic_r = _pool_matrix(hic_matrix, N)
    if max_diag is None:
        max_diag = N // 2

    # Contact proxy: 1/(d + ε)
    sim_contact = inverse_contact_matrix(dist_map)

    # ── 1. Diagonal decay ─────────────────────────────────────────────────────
    sim_decay = diagonal_decay_profile(sim_contact, max_diag)
    exp_decay = diagonal_decay_profile(hic_r,       max_diag)
    r_dd, p_dd = _pearson(sim_decay, exp_decay)

    # ── 2. Insulation score ───────────────────────────────────────────────────
    sim_ins = insulation_score(sim_contact, insulation_window)
    exp_ins = insulation_score(hic_r,       insulation_window)
    r_ins, p_ins = _pearson(sim_ins, exp_ins)

    # ── 3. PC1 (A/B compartments) ─────────────────────────────────────────────
    sim_pc1 = pc1_of_oe(sim_contact)
    exp_pc1 = pc1_of_oe(hic_r)
    r_pc1, p_pc1 = _pearson(np.abs(sim_pc1), np.abs(exp_pc1))

    log_table(
        [
            ("Diagonal decay r",    f"{r_dd:.4f}  (p={p_dd:.2e})"),
            ("Insulation score r",  f"{r_ins:.4f}  (p={p_ins:.2e})"),
            ("|PC1| r",             f"{r_pc1:.4f}  (p={p_pc1:.2e})"),
        ],
        title="Hi-C Validation — single structure",
        log_fn=_log.info,
    )

    # ── 4. Direct matrix similarity ───────────────────────────────────────────
    sim_flat = _upper_tri(sim_contact)
    exp_flat = _upper_tri(hic_r)
    r_pearson, p_pearson = _pearson(sim_flat, exp_flat)
    r_spearman, p_spearman = _spearman(sim_flat, exp_flat)

    log_table(
        [
            ("Pearson r (matrix)",   f"{r_pearson:.4f}  (p={p_pearson:.2e})"),
            ("Spearman r (matrix)",  f"{r_spearman:.4f}  (p={p_spearman:.2e})"),
        ],
        title="Hi-C Validation — direct similarity",
        log_fn=_log.info,
    )

    if save_path is not None:
        from .plots import plot_hic_comparison
        import os
        plots_dir = os.path.join(save_path, "plots")
        os.makedirs(plots_dir, exist_ok=True)
        plot_hic_comparison(sim_contact, hic_r, plots_dir, name="hic_comparison_single")

    return {
        "diag_decay_r": r_dd,      "diag_decay_p": p_dd,
        "insulation_r": r_ins,     "insulation_p": p_ins,
        "pc1_r":        r_pc1,     "pc1_p":        p_pc1,
        "pearson_r":    r_pearson, "pearson_p":    p_pearson,
        "spearman_r":   r_spearman,"spearman_p":   p_spearman,
    }


def validate_hic_ensemble(
    cif_paths: list,
    hic_matrix: np.ndarray,
    insulation_window: int = 10,
    max_diag: int | None = None,
    eps: float = 1e-3,
    save_path: str | None = None,
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
        'pc1_r',        'pc1_p'
    """
    from .utils import get_coordinates_cif, model_distance_heatmap
    _log = log or logger

    _log.info(f"Computing inverse-average contact map from {len(cif_paths)} frames …")
    inv_avg = None
    for path in cif_paths:
        coords   = get_coordinates_cif(path)
        dist_map = model_distance_heatmap(coords)
        inv      = inverse_contact_matrix(dist_map, eps=eps)
        if inv_avg is None:
            inv_avg = inv
        else:
            inv_avg += inv
    inv_avg /= len(cif_paths)

    N     = inv_avg.shape[0]
    hic_r = _pool_matrix(hic_matrix, N)
    if max_diag is None:
        max_diag = N // 2

    # ── 1. Diagonal decay ─────────────────────────────────────────────────────
    sim_decay = diagonal_decay_profile(inv_avg, max_diag)
    exp_decay = diagonal_decay_profile(hic_r,   max_diag)
    r_dd, p_dd = _pearson(sim_decay, exp_decay)

    # ── 2. Insulation score ───────────────────────────────────────────────────
    sim_ins = insulation_score(inv_avg, insulation_window)
    exp_ins = insulation_score(hic_r,   insulation_window)
    r_ins, p_ins = _pearson(sim_ins, exp_ins)

    # ── 3. PC1 (A/B compartments) ─────────────────────────────────────────────
    sim_pc1 = pc1_of_oe(inv_avg)
    exp_pc1 = pc1_of_oe(hic_r)
    r_pc1, p_pc1 = _pearson(np.abs(sim_pc1), np.abs(exp_pc1))

    log_table(
        [
            ("Frames averaged",     str(len(cif_paths))),
            ("Diagonal decay r",    f"{r_dd:.4f}  (p={p_dd:.2e})"),
            ("Insulation score r",  f"{r_ins:.4f}  (p={p_ins:.2e})"),
            ("|PC1| r",             f"{r_pc1:.4f}  (p={p_pc1:.2e})"),
        ],
        title="Hi-C Validation — ensemble",
        log_fn=_log.info,
    )

    # ── 4. Direct matrix similarity ───────────────────────────────────────────
    sim_flat = _upper_tri(inv_avg)
    exp_flat = _upper_tri(hic_r)
    r_pearson, p_pearson = _pearson(sim_flat, exp_flat)
    r_spearman, p_spearman = _spearman(sim_flat, exp_flat)

    log_table(
        [
            ("Pearson r (matrix)",   f"{r_pearson:.4f}  (p={p_pearson:.2e})"),
            ("Spearman r (matrix)",  f"{r_spearman:.4f}  (p={p_spearman:.2e})"),
        ],
        title="Hi-C Validation — direct similarity",
        log_fn=_log.info,
    )

    if save_path is not None:
        from .plots import plot_hic_comparison
        import os
        plots_dir = os.path.join(save_path, "plots")
        os.makedirs(plots_dir, exist_ok=True)
        plot_hic_comparison(inv_avg, hic_r, plots_dir, name="hic_comparison_ensemble")

    return {
        "diag_decay_r": r_dd,      "diag_decay_p": p_dd,
        "insulation_r": r_ins,     "insulation_p": p_ins,
        "pc1_r":        r_pc1,     "pc1_p":        p_pc1,
        "pearson_r":    r_pearson, "pearson_p":    p_pearson,
        "spearman_r":   r_spearman,"spearman_p":   p_spearman,
    }
