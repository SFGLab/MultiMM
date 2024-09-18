from utils import *
from tqdm import tqdm
from scipy.stats import pearsonr, spearmanr
from scipy.spatial import distance
from scipy.spatial.distance import pdist, squareform
from scipy.sparse import csr_matrix
from scipy.stats import pearsonr
import numpy as np
from sklearn.neighbors import NearestNeighbors
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from scipy.sparse.linalg import eigsh

# Metrics for comparisons
def calculate_correlation(matrix1, matrix2):
    # Flatten matrices and calculate Pearson correlation
    flat1 = matrix1.flatten()
    flat2 = matrix2.flatten()
    correlation, p_val = pearsonr(flat1, flat2)
    return correlation, p_val

def rv_coefficient(matrix1, matrix2):
    """
    Computes the RV coefficient between two matrices.
    
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
    """
    Computes the Mantel test correlation between two distance matrices.
    
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
    """
    Compute Pearson correlation efficiently between two flattened arrays.
    """
    m1_flat = m1.flatten()
    m2_flat = m2.flatten()
    
    # Compute covariance and standard deviations
    cov = np.mean((m1_flat - m1_flat.mean()) * (m2_flat - m2_flat.mean()))
    std_m1 = np.std(m1_flat)
    std_m2 = np.std(m2_flat)
    
    return cov / (std_m1 * std_m2)

def compute_pearson_correlation(m1, m2, window_size):
    """
    Compute the average Pearson correlation for non-overlapping windows of a given size.
    """
    N = m1.shape[0]
    correlations = []
    
    # Iterate over non-overlapping windows
    for i in range(0, N, window_size):
        for j in range(0, N, window_size):
            # Make sure the window fits completely within the matrix
            if i + window_size <= N and j + window_size <= N:
                sub_m1 = m1[i:i+window_size, j:j+window_size]
                sub_m2 = m2[i:i+window_size, j:j+window_size]
                
                # Compute Pearson correlation for the submatrices
                correlation = fast_pearson_correlation(sub_m1, sub_m2)
                correlations.append(correlation)
    
    return np.mean(correlations)

def correlation_vs_window_size(m1, m2):
    """
    Compute the average Pearson correlation for selected window sizes and plot the result.
    """
    N = m1.shape[0]
    
    # Select window sizes that divide N evenly
    window_sizes = np.arange(100,3*N//4,100)
    avg_correlations = []

    for window_size in tqdm(window_sizes):
        avg_correlation = compute_pearson_correlation(m1, m2, window_size)
        avg_correlations.append(avg_correlation)
    
    # Plotting the results
    plt.plot(window_sizes, avg_correlations,'k-o')
    plt.xlabel('Window Size')
    plt.ylabel('Average Pearson Correlation')
    plt.grid(True)
    plt.show()
    return window_sizes, avg_correlations

# Random Walks
def random_walk_3d(N):
    """
    Generate a random walk with N points in 3D.
    
    Parameters:
    N (int): Number of points in the random walk
    
    Returns:
    np.ndarray: Array of shape (N, 3) with coordinates of the walk
    """
    if N < 1:
        raise ValueError("N must be at least 1")

    # Directions for movement in 3D: x, y, z
    directions = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], 
                           [-1, 0, 0], [0, -1, 0], [0, 0, -1]])
    
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
    """
    Generates a self-avoiding random walk in 3D with N points, where consecutive points
    have a constant distance (step_size) and beads avoid each other, using a backtracking approach.
    
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
    directions = np.array([[1, 0, 0], [-1, 0, 0], 
                           [0, 1, 0], [0, -1, 0], 
                           [0, 0, 1], [0, 0, -1]])
    
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
                print(f"Exceeded maximum backtracks. Unable to complete the walk with {N} points.")
                return walk[:current_index]  # Return the walk up to the last valid point

    return walk

def structure_to_heatmap(V, use_sparse=False):
    # Calculate pairwise Euclidean distances
    dist_matrix = distance.cdist(V, V, 'euclidean')
    
    # Invert distances, avoiding division by zero
    inv_dist_matrix = 1.0 / (dist_matrix+1)**3/2
    
    return inv_dist_matrix

def rescale_matrix(matrix, target_size):
    # Rescale or coarse-grain the matrix to the target size (nxn)
    N = matrix.shape[0]
    indices = np.linspace(0, N - 1, target_size, dtype=int)
    rescaled_matrix = matrix[np.ix_(indices, indices)]
    return rescaled_matrix

# Preprocessing
def mean_downsample(V, target_size):
    """
    Downsamples a 3D structure from Nx3 to nx3 by computing moving averages.
    
    Parameters:
    V (ndarray): The input array of shape (N, 3).
    target_size (int): The desired number of rows in the downsampled array (n).
    
    Returns:
    V_downsampled (ndarray): The downsampled array of shape (n, 3).
    """
    N, dims = V.shape
    assert dims == 3, "Input array must have shape (N, 3)"
    assert target_size < N, "Target size must be less than the original size"
    
    # Calculate the window size for downsampling
    window_size = N / target_size
    
    # Initialize the downsampled structure
    V_downsampled = np.zeros((target_size, dims))
    
    # Compute moving averages for each downsampled point
    for i in range(target_size):
        start_idx = int(i * window_size)
        end_idx = int(min(start_idx + window_size, N))
        V_downsampled[i] = np.mean(V[start_idx:end_idx], axis=0)
    
    return V_downsampled

def pca_downsample(V, n):
    pca = PCA(n_components=3)
    V_reduced = pca.fit_transform(V)
    indices = np.linspace(0, V_reduced.shape[0] - 1, n, dtype=int)
    V_downgraded = V_reduced[indices, :]
    return V_downgraded

def remove_zero_rows_and_columns(matrix):
    # Convert input to numpy array if it isn't already
    matrix = np.array(matrix)
    
    # Find rows where all elements are zero
    zero_rows = np.where(~matrix.any(axis=1))[0]
    
    # Find columns where all elements are zero
    zero_columns = np.where(~matrix.any(axis=0))[0]
    
    # Remove the zero rows and columns
    matrix = np.delete(matrix, zero_rows, axis=0)
    matrix = np.delete(matrix, zero_columns, axis=1)
    
    return matrix, zero_rows, zero_columns

def remove_diagonals(matrix, n_diag):
    """
    Removes the specified number of diagonals from a square matrix.
    
    Parameters:
    - matrix (np.array): The input square matrix.
    - n_diag (int): Number of diagonals to remove.
    
    Returns:
    - modified_matrix (np.array): Matrix with specified diagonals removed.
    """
    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError("The input matrix must be square.")
    mean = np.mean(matrix)
    modified_matrix = matrix.copy()
    
    # Create a mask for diagonals to be removed
    n = matrix.shape[0]
    
    for d in range(n_diag + 1):
        # Mask for main diagonal and `d` diagonals above and below
        mask = np.zeros_like(matrix, dtype=bool)
        
        # Main diagonal and diagonals above and below
        for i in range(n):
            if i + d < n:
                mask[i, i + d] = True  # Diagonals above main
                mask[i + d, i] = True  # Symmetric diagonals below main
            if i - d >= 0:
                mask[i, i - d] = True  # Diagonals below main
                mask[i - d, i] = True  # Symmetric diagonals above main
        
        # Remove diagonals by setting them to zero
        modified_matrix[mask] = mean
    
    return modified_matrix

def minMax(v,Max=1,Min=0):
    return (Max-Min)*(v-np.min(v))/(np.max(v)-np.min(v))+Min

def standarize(v):
    return (v-np.mean(v))/np.std(v)

# Compute eigenvectors
def compute_compartments(matrix):
    """
    Computes the correlation matrix from a Hi-C contact matrix,
    performs eigenvalue decomposition, and returns the eigenvector
    corresponding to the highest eigenvalue.
    
    Parameters:
    - matrix: 2D numpy array, the Hi-C contact matrix (should be square and symmetric)
    
    Returns:
    - eigenvector: 1D numpy array, the eigenvector corresponding to the highest eigenvalue
    """
    # Ensure matrix is symmetric
    assert matrix.shape[0] == matrix.shape[1], "Matrix must be square"
    assert np.allclose(matrix, matrix.T), "Matrix is not symmetric"
    
    # Compute the correlation matrix
    # Normalizing the matrix (remove diagonal for correlation computation)
    matrix = np.nan_to_num(matrix)  # Replace NaNs with 0
    np.fill_diagonal(matrix, 0)     # Remove the diagonal
    
    # Compute the correlation matrix
    correlation_matrix = np.corrcoef(matrix, rowvar=False)
    # correlation_matrix = matrix
    
    # Perform eigenvalue decomposition
    eigvals, eigvecs = eigsh(correlation_matrix)
    
    # Get indices of the top two eigenvalues
    indices = np.argsort(eigvals)[::-1]  # Sort in descending order
    top_two_indices = indices[:2]
    
    # Get the top two eigenvectors
    eigenvector1 = eigvecs[:, top_two_indices[0]]
    eigenvector2 = eigvecs[:, top_two_indices[1]]
    
    return eigenvector1, eigenvector2

def compare_matrices(m,mr,exp_m, viz=True):
    # Remove empty rows
    exp_m, rs, cs = remove_zero_rows_and_columns(exp_m)
    m = np.delete(m, rs, axis=0)
    m = np.delete(m, cs, axis=1)
    mr = np.delete(mr, rs, axis=0)
    mr = np.delete(mr, cs, axis=1)

    # Normalize
    m, mr, exp_m = standarize(m), standarize(mr), standarize(exp_m)
    m, mr, exp_m = minMax(m), minMax(mr), minMax(exp_m)

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
        print('Correlation of simulation with the first eigenvector:',corr_sim1)
        print('Correlation of random walk with the first eigenvector:',corr_rw1)
        print('Correlation of simulation with the second eigenvector:',corr_sim2)
        print('Correlation of random walk with the second eigenvector:',corr_rw2)

    return corr_sim1, corr_rw1, corr_sim2, corr_rw2

def pipeline_single_ensemble(V,Vr,exp_m,viz=True):
    # Downgrade structures to have same number of points as the experimental heatmaps
    V_reduced = mean_downsample(V, len(exp_m))
    Vr_reduced = mean_downsample(Vr, len(exp_m))

    # Compute simulated heatmaps
    m = structure_to_heatmap(V_reduced)
    mr = structure_to_heatmap(Vr_reduced)

    # Compute the statistics
    corr_sim1, corr_rw1, corr_sim2, corr_rw2 = compare_matrices(m,mr,exp_m,viz)

def ensemble_pipeline(ensemble_path,exp_m,N_chroms=22,N_ens=100,viz=True):
    # Run loop for each chromosome
    Cs_sim1, Cs_sim2 = list(), list()
    Cs_rw1, Cs_rw2 = list(), list()
    avg_ms, avg_mrs = list(), list()
    L = len(exp_m)
    for i in range(N_chroms):
        # Average the heatmaps of each ensemble
        avg_m = 0
        print(f'Chromosome {i+1}:')
        print('Calculating heatmaps from experimental structures...')
        for j in tqdm(range(N_ens)):
            V = get_coordinates_cif(ensemble_path+f'/ens_{j+1}/chromosomes/MultiMM_minimized_chr{i+1}.cif')
            V_reduced = mean_downsample(V, L)
            m = structure_to_heatmap(V_reduced)
            avg_m += m
        avg_m /= N_ens
        avg_ms.append(avg_m)

        # Average the heatmaps of each random walk
        avg_mr = 0
        print('Calculating heatmaps from random structures...')
        for j in tqdm(range(N_ens)):
            Vr = random_walk_3d(len(V))
            Vr_reduced = mean_downsample(Vr, L)
            mr = structure_to_heatmap(Vr_reduced)
            avg_mr += mr
        avg_mr /= N_ens 
        avg_mrs.append(avg_mr)

        # Compare them
        corr_sim1, corr_rw1, corr_sim2, corr_rw2 = compare_matrices(avg_m,avg_mr,exp_m,False)
        Cs_sim1.append(corr_sim1)
        Cs_sim2.append(corr_sim2)
        Cs_rw1.append(corr_rw1)
        Cs_rw2.append(corr_rw2)

        if viz:
            print('Correlation of simulation with the first eigenvector:',corr_sim1)
            print('Correlation of random walk with the first eigenvector:',corr_rw1)
            print('Correlation of simulation with the second eigenvector:',corr_sim2)
            print('Correlation of random walk with the second eigenvector:',corr_rw2)
        print(f'Chromosome {i+1} done!\n')

    if viz:
        X_axis = np.arange(N_chroms) 
        plt.figure(figsize=(16, 6),dpi=200)
        plt.bar(X_axis - 0.2, Cs_sim1, 0.4, label = 'Simulation', color='blue') 
        plt.bar(X_axis + 0.2, Cs_rw1, 0.4, label = 'Random Walk', color= 'red') 

        plt.xticks(X_axis, chroms)
        plt.xlabel("Chromosomes") 
        plt.ylabel("Correlation") 
        plt.legend()
        plt.show()
    
    return avg_ms, avg_mrs

# Import structures
path = '/home/skorsak/Data/simulation_results/gw_ens'
exp_m = np.nan_to_num(np.load('/home/skorsak/Data/Rao/chrom5_primary+replicate_ice_norm_25kb.npy'))

avg_ms, avg_mrs = ensemble_pipeline(path,exp_m)