from utils import *
from tqdm import tqdm
from scipy.stats import pearsonr, spearmanr
from scipy.spatial import distance
from scipy.spatial.distance import pdist, squareform
from scipy.sparse import csr_matrix
from scipy.stats import pearsonr
import scipy.ndimage as ndimage
from scipy.spatial.distance import cdist
from scipy.spatial import KDTree
import numpy as np
from sklearn.neighbors import NearestNeighbors
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from scipy.sparse.linalg import eigsh
import seaborn as sns

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

    return np.abs(corr_sim1), np.abs(corr_rw1), np.abs(corr_sim2), np.abs(corr_rw2)

def pipeline_single_ensemble(V,Vr,exp_m,viz=True):
    # Downgrade structures to have same number of points as the experimental heatmaps
    V_reduced = mean_downsample(V, len(exp_m))
    Vr_reduced = mean_downsample(Vr, len(exp_m))

    # Compute simulated heatmaps
    m = structure_to_heatmap(V_reduced)
    mr = structure_to_heatmap(Vr_reduced)

    # Compute the statistics
    corr_sim1, corr_rw1, corr_sim2, corr_rw2 = compare_matrices(m,mr,exp_m,viz)

def ensemble_pipeline_boxplot(ensemble_path, exp_path, N_chroms=22, N_ens=20, viz=True):
    # Lists to hold correlation values
    Cs_sim, Cs_rw = list(), list()
    
    for i in range(N_chroms):
        # Load experimental data
        exp_m = np.nan_to_num(np.load(exp_path + f'/chrom{i+1}_primary+replicate_ice_norm_50kb.npy'))
        exp_m = remove_diagonals(exp_m, 5)
        L = len(exp_m)
        
        # Calculate heatmaps from experimental structures
        print(f'Chromosome {i+1}:')
        print('Calculating heatmaps from experimental structures...')
        corrs_sim, corrs_rw = [], []
        for j in tqdm(range(N_ens)):
            V = get_coordinates_cif(ensemble_path + f'/ens_{j+1}/chromosomes/MultiMM_minimized_chr{i+1}.cif')
            V_reduced = mean_downsample(V, L)
            m = structure_to_heatmap(V_reduced)

            Vr = random_walk_3d(len(V))
            Vr_reduced = mean_downsample(Vr, L)
            mr = structure_to_heatmap(Vr_reduced)
        
            # Calculate correlations
            corr_sim1, corr_rw1 ,_ ,_ = compare_matrices(m, mr, exp_m, False)
            corrs_sim.append(np.abs(corr_sim1))
            corrs_rw.append(np.abs(corr_rw1))

        Cs_sim.append(corrs_sim)
        Cs_rw.append(corrs_rw)
        print(f'Chromosome {i+1} done!\n')

    if viz:
        # Box plot for each chromosome
        plt.figure(figsize=(20, 5), dpi=200)
        
        # Prepare data for box plots
        data_sim = [Cs_sim[i] for i in range(N_chroms)]
        data_rw = [Cs_rw[i] for i in range(N_chroms)]
        
        box_sim = plt.boxplot(data_sim, positions=np.arange(N_chroms) - 0.2, widths=0.4, patch_artist=True,
                              boxprops=dict(facecolor='blue', color='blue'), medianprops=dict(color='black'))
        box_rw = plt.boxplot(data_rw, positions=np.arange(N_chroms) + 0.2, widths=0.4, patch_artist=True,
                             boxprops=dict(facecolor='red', color='red'), medianprops=dict(color='black'))

        plt.xticks(np.arange(N_chroms), [f'chr{i + 1}' for i in range(N_chroms)])
        plt.xlabel("Chromosomes", fontsize=16)
        plt.ylabel("Correlation with 1st Eigenvector", fontsize=14)
        plt.legend([box_sim["boxes"][0], box_rw["boxes"][0]], ['Simulation', 'Random Walk'], loc='upper right')
        plt.savefig('heatmap_correlation_boxplots.pdf', format='pdf', dpi=200)
        plt.savefig('heatmap_correlation_boxplots.svg', format='svg', dpi=200)
        plt.show()

def ensemble_pipeline_bars(ensemble_path,exp_path,N_chroms=22,N_ens=20,viz=True):
    # Run loop for each chromosome
    Cs_sim1, Cs_sim2 = list(), list()
    Cs_rw1, Cs_rw2 = list(), list()
    
    for i in range(N_chroms):
        # Experimental path
        exp_m = np.nan_to_num(np.load(exp_path+f'/chrom{i+1}_primary+replicate_ice_norm_50kb.npy'))
        L = len(exp_m)

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

        # Average the heatmaps of each random walk
        avg_mr = 0
        print('Calculating heatmaps from random structures...')
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
        chroms = ['chr'+str(i) for i in range(1,N_chroms+1)]
        X_axis = np.arange(N_chroms) 
        plt.figure(figsize=(20, 5),dpi=200)
        plt.bar(X_axis - 0.2, Cs_sim1, 0.4, label = 'Simulation', color='blue') 
        plt.bar(X_axis + 0.2, Cs_rw1, 0.4, label = 'Random Walk', color= 'red')
        plt.xticks(X_axis, chroms)
        plt.xlabel("Chromosomes",fontsize=16) 
        plt.legend()
        plt.ylabel("Correlation with First Eigenvector",fontsize=14) 
        plt.savefig('corr_1st_eigenvec.pdf',format='pdf',dpi=200)
        plt.savefig('corr_1st_eigenvec.svg',format='svg',dpi=200)
        
        plt.show()

        plt.figure(figsize=(20, 5),dpi=200)
        plt.bar(X_axis - 0.2, Cs_sim2, 0.4, label = 'Simulation', color='blue') 
        plt.bar(X_axis + 0.2, Cs_rw2, 0.4, label = 'Random Walk', color= 'red') 

        plt.xticks(X_axis, chroms)
        plt.xlabel("Chromosomes",fontsize=16)
        plt.legend()
        plt.ylabel("Correlation with Second Eigenvector",fontsize=14)
        plt.savefig('corr_2st_eigenvec.pdf',format='pdf',dpi=200)
        plt.savefig('corr_2st_eigenvec.svg',format='svg',dpi=200)
        plt.show()

def regions_pipeline(regions_dir,chroms,starts,ends,N_ens=1000):
    # Experimental path
    corrs_sim, corrs_rw = [], []
    pvals_sim, pvals_rw = [], []
    ps_sim, ints_sim = [], []
    ps_rw, ints_rw = [], []

    # Average the heatmaps of each ensemble
    avg_m = 0
    print('Calculating heatmaps from experimental structures...')
    for i in tqdm(range(N_ens)):
        try:
            exp_m = np.nan_to_num(np.load(f'/home/skorsak/Data/Rao/regs/primary+replicate_ice_norm_25kb_chrom{chroms[i]}_{starts[i]}_{ends[i]}.npy'))
            L = len(exp_m)
        except:
            print(f"Problem with chromosome {chroms[i]}_{starts[i]}_{ends[i]} experimental data")
            continue
        
        try:
            V = get_coordinates_cif(regions_dir+f'/regens_chrom{chroms[i]}_region_{starts[i]}_{ends[i]}/MultiMM_minimized.cif')
        except:
            print(f"Problem with chromosome {chroms[i]}_{starts[i]}_{ends[i]} simulated data")
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

        m = (m-np.mean(m))/np.std(m)
        mr = (mr-np.mean(mr))/np.std(mr)
        exp_m = (exp_m-np.mean(exp_m))/np.std(exp_m)
        m = (m-np.min(m))/(np.max(m)-np.min(m))
        mr = (mr-np.min(mr))/(np.max(mr)-np.min(mr))
        exp_m = (exp_m-np.min(exp_m))/(np.max(exp_m)-np.min(exp_m))

        p_sim, int_sim = analyze_heatmaps(remove_diagonals(m,4),remove_diagonals(exp_m,4), min_distance=5, distance_threshold=5)
        p_rw, int_rw = analyze_heatmaps(remove_diagonals(mr,4), remove_diagonals(exp_m,4), min_distance=5, distance_threshold=5)

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
    plt.xticks([0, 1], ['Simulation', 'Random Walk'],fontsize=16)
    plt.ylabel('Correlation with Experimental Data',fontsize=16)
    plt.savefig('violin.svg',format='svg',dpi=200)
    plt.savefig('violin.pdf',format='pdf',dpi=200)
    plt.show()

    # Create the violin plot
    data = [np.array(ps_sim)/100, np.array(ps_rw)/100]
    plt.figure(figsize=(6, 9))
    sns.violinplot(data=data)
    plt.xticks([0, 1], ['Simulation', 'Random Walk'],fontsize=16)
    plt.ylabel('Percentage of Common Loops',fontsize=16)
    plt.savefig('violin_ps.pdf',format='pdf',dpi=200)
    plt.show()

    # Create the violin plot
    data = [ints_sim, ints_rw]
    plt.figure(figsize=(6, 9))
    sns.violinplot(data=data)
    plt.xticks([0, 1], ['Simulation', 'Random Walk'],fontsize=16)
    plt.ylabel('Peak Intensity Correlation',fontsize=16)
    plt.savefig('violin_ints.pdf',format='pdf',dpi=200)
    plt.show()

# Import structures
sim_path = '/home/skorsak/Data/simulation_results/gw_ens_5M'
exp_path = '/home/skorsak/Data/Rao'

ensemble_pipeline_bars(sim_path,exp_path)

# # Load data
# regs_path = '/home/skorsak/Data/simulation_results/ensembles_of_regions'
# starts = np.load('/home/skorsak/Data/simulation_results/starts.npy')
# ends = np.load('/home/skorsak/Data/simulation_results/ends.npy')
# chroms = np.load('/home/skorsak/Data/simulation_results/chroms.npy')

# regions_pipeline(regs_path,chroms,starts,ends)
