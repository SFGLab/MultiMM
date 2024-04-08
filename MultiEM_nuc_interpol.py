import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from MultiEM_utils import *
from MultiEM_init_tools import *

def interpolate_structure_with_nucleosomes(V, bw_array, N_nuc=3):
    """
    Interpolate the 3D structure V with nucleosomes.
    """
    # Calculate mean of bw_array
    mean_signal = np.mean(bw_array)
    
    # Normalize bw_array
    norm_bw_array = bw_array / mean_signal
    
    # Initialize interpolated structure
    interpolated_structure = []
    
    # Interpolate each segment with nucleosomes
    for i in range(len(V) - 1):
        # Get segment endpoints
        start_point = V[i]
        end_point = V[i + 1]
        
        # Calculate segment length
        segment_length = np.linalg.norm(end_point - start_point)
        
        # Calculate number of nucleosomes in segment
        num_nucleosomes = int(np.round(norm_bw_array[i] * N_nuc))
        
        # Generate helices for nucleosomes
        helices = generate_nucleosome_helices(start_point, end_point, num_nucleosomes)
        
        # Append helices to interpolated structure
        interpolated_structure.extend(helices)
    
    return np.concatenate(np.array(interpolated_structure))

def generate_nucleosome_helices(start_point, end_point, num_nucleosomes, turns=1.6):
    """
    Generate helices for nucleosomes in a segment.
    """
    # Calculate segment vector
    segment_vector = end_point - start_point
    
    # Calculate helix parameters
    helix_radius = np.linalg.norm(segment_vector) / (2 * np.pi * turns)
    helix_axis = segment_vector / np.linalg.norm(segment_vector)
    helix_height = np.linalg.norm(segment_vector) / num_nucleosomes
    
    # Initialize helices
    helices = []
    
    # Generate helices for each nucleosome
    for i in range(num_nucleosomes):
        # Calculate helix center
        helix_center = start_point + (i + 0.5) * helix_height * helix_axis
        
        # Calculate helix angle
        helix_angle = np.pi / 6 if i % 2 == 0 else -np.pi / 6  # Zig-zag configuration
        
        # Generate helix coordinates
        theta = np.linspace(0, 2 * np.pi, 100)
        x = helix_center[0] + helix_radius * np.cos(theta)
        y = helix_center[1] + helix_radius * np.sin(theta)
        z = helix_center[2] + helix_height * theta / (2 * np.pi) + helix_angle
        
        helix = np.vstack([x, y, z]).T
        helices.append(helix)
    
    return helices

# Example data
V = get_coordinates_cif('/home/skorsak/Templates/MultiEM-main/test/MultiEM_minimized.cif')
N = len(V)
bw_array = np.random.rand(N)  # Mock signal array

# Interpolate structure with nucleosomes
interpolated_structure = interpolate_structure_with_nucleosomes(V, bw_array)
points = np.arange(0,len(interpolated_structure))

# Plot the result
fig = plt.figure(figsize=(30,30))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(interpolated_structure[:, 0], interpolated_structure[:, 1], interpolated_structure[:, 2], s=0.5, c=points, cmap='rainbow',edgecolors='k')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Interpolated Structure with Nucleosomes')
plt.savefig('/home/skorsak/Templates/MultiEM-main/structure.pdf',format='pdf',dpi=300)
plt.show()

# write_mmcif_chrom(interpolated_structure,'/home/skorsak/Templates/MultiEM_main/with_nucs.cif')