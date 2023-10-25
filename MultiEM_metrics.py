from MultiEM_utils import *
from scipy.spatial import ConvexHull
import numpy as np
from sklearn.metrics import silhouette_score, davies_bouldin_score
from sklearn.cluster import KMeans
from scipy.stats import entropy
from ellipsoid_measure import calc_ellipsoid_ratio2

def gyration_radius(V):
    '''
    Import the structure coordinates V and export radius of gyration.
    We assume the same mass for all atoms.
    
    V (numpy array): Dimension MXN, where M number of V and N
                     is the always 3 as the dimension of our space.
    Rg (float): Radius of gyration. 
    '''
    r_mean, N = np.average(V,axis=0), len(V)
    Rg_sq = 1/N*np.sum(np.linalg.norm(V-r_mean,axis=0)**2)
    Rg = np.sqrt(Rg_sq)
    return Rg

def get_stats(V):
    '''
    Imports structure and returns basic statistics of it.
    '''
    res = stats.describe(V)
    return np.concatenate([
        [
            res.minmax[0],
            res.minmax[1],
            res.mean,
            res.variance,
            res.skewness,
            res.kurtosis
        ],
        np.percentile(V, q=[10, 25, 50, 75, 90, 95])
    ])

# Define a function to calculate the box-counting fractal dimension
def fractal_dimension_3d(V, scale=30, viz=False):
    '''
    Import structure and return fractal dimension.

    ---- Made by ChatGPT. --- ;)
    '''
    min_coords = np.min(V, axis=0)
    max_coords = np.max(V, axis=0)
    dims = []
    scales = []

    for current_scale in range(1, scale + 1):
        # Divide the space into boxes of the current scale
        box_size = (max_coords - min_coords) / current_scale
        n_boxes = np.ceil((max_coords - min_coords) / box_size).astype(int)+1

        # Create an empty grid to count V in each box
        grid = np.zeros(n_boxes, dtype=int)

        for point in V:
            # Calculate the box in which the point falls
            box_x = int((point[0] - min_coords[0]) / box_size[0])
            box_y = int((point[1] - min_coords[1]) / box_size[1])
            box_z = int((point[2] - min_coords[2]) / box_size[2])

            # Update the grid
            grid[box_x, box_y, box_z] += 1

        # Count non-empty boxes
        N = np.sum(grid > 0)

        # Store the results for this scale
        dims.append(np.log(N))
        scales.append(np.log(1.0 / current_scale))

    # Fit a linear regression to determine the fractal dimension
    coeffs = np.polyfit(scales, dims, 1)

    if viz:
        # Visualize the box counting and the linear regression
        plt.scatter(scales, dims, label='Box Counting Data')
        plt.plot(scales, np.polyval(coeffs, scales), color='red', label='Linear Fit')
        plt.xlabel('Log(1 / Scale)')
        plt.ylabel('Log(Number of Boxes)')
        plt.legend()
        plt.show()
    return -coeffs[0]

def sphericity(V):
    # Calculate the convex hull of the structure
    hull = ConvexHull(V)

    # Calculate the volume of the structure
    volume = hull.volume

    # Calculate the surface area of the structure
    surface_area = hull.area

    # Calculate sphericity using the formula
    sphericity = (np.pi**(1/3) * (6 * volume)**(2/3)) / surface_area

    return sphericity

def anisotropy(V):
    '''
    Gives anisotropy of a provided structure.

    If the value is close to 0 is isotropic, if it is close to 1
    it is anisotropic.
    '''
    # Calculate the center of mass (centroid) of the structure
    centroid = np.mean(V, axis=0)

    # Translate the structure to the centroid
    centered_V = V - centroid

    # Calculate the moment of inertia tensor
    moment_of_inertia = np.dot(centered_V.T, centered_V)

    # Calculate the eigenvalues of the moment of inertia tensor
    eigenvalues = np.linalg.eigvals(moment_of_inertia)

    # Calculate the anisotropy
    anisotropy = 1.5 * (max(eigenvalues) - min(eigenvalues)) / np.sum(eigenvalues)

    return anisotropy

def density(V, stick_radius=1): # check if stick radius is correct
    # Calculate the convex hull of the structure
    hull = ConvexHull(V)

    # Calculate the volume of the structure
    volume = hull.volume

    # Calculate the number of points
    num_points = len(V)

    # Calculate the volume of each stick and the total volume
    stick_volume = np.pi * stick_radius**2
    total_stick_volume = (num_points-1) * stick_volume

    # Calculate the density metric (number of points per unit volume)
    density_metric = num_points / volume

    return density_metric

def aspect_ratio(V):
    '''
    A value of 1 indicates that the structure is isotropic (uniform in all directions), 
    while values greater than 1 indicate anisotropy, with the extent of elongation along
    the longest principal component relative to the shortest one.
    '''
    # Calculate the mean of the data points
    mean_point = np.mean(V, axis=0)

    # Center the data by subtracting the mean
    centered_data = V - mean_point

    # Calculate the covariance matrix
    covariance_matrix = np.cov(centered_data, rowvar=False)

    # Perform PCA to find the principal components
    eigenvalues, _ = np.linalg.eigh(covariance_matrix)

    # Calculate the aspect ratio
    longest_component = np.max(eigenvalues)
    shortest_component = np.min(eigenvalues)
    aspect_ratio = longest_component / shortest_component

    return aspect_ratio

def compute_clustering_metrics_without_labels(V, n_clusters):
    # Perform K-Means clustering on the data
    kmeans = KMeans(n_clusters=n_clusters, random_state=0)
    cluster_labels = kmeans.fit_predict(V)

    # Calculate the Silhouette Score
    silhouette = silhouette_score(V, cluster_labels, metric='euclidean')

    # Calculate the Davies-Bouldin Index
    davies_bouldin = davies_bouldin_score(V, cluster_labels)

    return silhouette, davies_bouldin

def calculate_entropy(V):
    # Calculate the number of points in the structure
    num_points = len(V)

    # Calculate the distribution of points (you can use a 3D histogram)
    x, y, z = V[:, 0], V[:, 1], V[:, 2]
    hist, edges = np.histogramdd((x, y, z), bins=10)  # You can adjust the number of bins

    # Normalize the histogram to obtain a probability distribution
    prob_distribution = hist / num_points

    # Calculate the Shannon entropy
    entropy_value = entropy(prob_distribution.flatten())

    return entropy_value

def calculate_hurst_exponent(V, max_lag=10):
    """
    Calculate the Hurst exponent for a 3D spatial structure using a basic R/S analysis.
    
    Parameters:
    - V: Numpy array containing the 3D spatial structure data.
    - max_lag: Maximum lag for R/S analysis.

    Returns:
    - Hurst exponent estimate.
    """

    N = len(V)

    hurst_sum = 0

    for lag in range(1, max_lag + 1):
        lag_sum = 0

        for i in range(N - lag):
            lag_sum += np.linalg.norm(V[i + lag] - V[i])

        hurst_sum += lag_sum / (N - lag)

    hurst_sum /= max_lag

    return hurst_sum

def save_metrics(V,path_name=''):
    Rg = gyration_radius(V)
    D = fractal_dimension_3d(V)
    asp = aspect_ratio(V)
    sph = sphericity(V)
    ani = anisotropy(V)
    entr = calculate_entropy(V)
    ell = calc_ellipsoid_ratio2(V)

    f = open(path_name+'metrics_info.txt', "w")
    f.write(f'Gyration radius: {Rg}.\n')
    f.write(f'Fractal Dimension: {D}.\n')
    f.write(f'Aspect ratio: {asp}.\n')
    f.write(f'Sphericity: {sph}.\n')
    f.write(f'Anisotropy: {ani}.\n')
    f.write(f'Entropy: {entr}.\n')
    f.write(f'Ellipsoid ratio: {ell}.\n')
    f.close()

