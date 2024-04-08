import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import plotly.graph_objects as go
from MultiEM_utils import *
from MultiEM_init_tools import *
from tqdm import tqdm

def interpolate_structure_with_nucleosomes(V, bw_array, N_nuc=5):
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
    sign = 1 
    print('Building nucleosome structure...')
    for i in tqdm(range(len(V) - 1)):
        # Get segment endpoints
        start_point = V[i]
        end_point = V[i + 1]
        
        # Calculate segment length
        segment_length = np.linalg.norm(end_point - start_point)
        
        # Calculate number of nucleosomes in segment
        num_nucleosomes = int(np.round(norm_bw_array[i] * N_nuc))
        
        # Generate helices for nucleosomes
        helices, sign = generate_nucleosome_helices(start_point, end_point, num_nucleosomes, 1.6, sign)
        
        # Append helices to interpolated structure
        interpolated_structure.extend(helices)
    
    return np.concatenate(np.array(interpolated_structure))

def rotate_with_matrix(V,theta,axis='x'):
    # Define rotation
    if axis=='x':
        R = np.array([[1, 0, 0],
                      [0, np.cos(theta), -np.sin(theta)],
                      [0, np.sin(theta), np.cos(theta)]])
    elif axis=='y':
        R = np.array([[np.cos(theta), 0, np.sin(theta)],
                [0, 1, 0],
                [-np.sin(theta), 0, np.cos(theta)]])
    elif axis=='z':
        R = np.array([[np.cos(theta), -np.sin(theta), 0],
                    [np.sin(theta), np.cos(theta), 0],
                    [0, 0, 1]])

    # Apply the rotation to the structure V
    V_rotated = np.dot(V, R)
    return V_rotated

def make_helix(r,theta,z0):
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    z = z0 * theta[::-1] / (2 * np.pi)
    return np.vstack([x, y, z]).T

def rotate_structure(V,phi):
    x, y ,z = V[:,0], V[:,1], V[:,2]
    x = x
    y = y*np.cos(phi)-z*np.sin(phi)
    z = y*np.sin(phi)+z*np.cos(phi)
    return np.vstack([x, y, z]).T

def generate_nucleosome_helices(start_point, end_point, num_nucleosomes, turns=1.6, sign=1):
    """
    Generate helices for nucleosomes in a segment.
    """
    # Calculate segment vector
    segment_vector = end_point - start_point

    # Calculate helix parameters
    helix_radius = np.linalg.norm(segment_vector) / (2 * np.pi) /5
    helix_axis = segment_vector / np.linalg.norm(segment_vector)
    helix_height = 0.04

    # Initialize helices
    theta = np.linspace(0, turns* 2 * np.pi, 20)
    helices = list()

    # Generate helices for each nucleosome
    for i in range(num_nucleosomes):
        # Calculate helix center
        phi = np.pi/2
        helix_center = start_point + (i + 0.5) * helix_height * helix_axis

        # Calculate helix angle
        helix = make_helix(helix_radius,theta,helix_height)
        helix = helix+helix_axis

        # Generate helix coordinates        
        helices.append(helix)
    if num_nucleosomes%2==0 or num_nucleosomes==1: sign*=-1
    return helices, sign

# Example data
V = get_coordinates_cif('/mnt/raid/codes/mine/MultiEM-main/Trios_ensembles/CHS_d/MultiEM_minimized.cif')
N = len(V)
bw_array = np.random.rand(N)  # Mock signal array

# Interpolate structure with nucleosomes
iV = interpolate_structure_with_nucleosomes(V, bw_array)
points = np.arange(0,len(iV))

fig = go.Figure(data=go.Scatter3d(
    x=iV[:5000,0], y=iV[:5000,1], z=iV[:5000,2],
    marker=dict(
        size=2,
        color=points[:5000],
        colorscale='rainbow',
    ),
    line=dict(
        color='darkblue',
        width=1
    )
))

fig.update_layout(
    width=800,
    height=700,
    autosize=False,
    scene=dict(
        camera=dict(
            up=dict(
                x=0,
                y=0,
                z=1
            ),
            eye=dict(
                x=0,
                y=1.0707,
                z=1,
            )
        ),
        aspectratio = dict( x=1, y=1, z=0.7 ),
        aspectmode = 'manual'
    ),
)

fig.show()

# # Plot the result
# fig = plt.figure(figsize=(10,10),dpi=200)
# ax = fig.add_subplot(111, projection='3d')
# ax.plot(interpolated_structure[:1000, 0], interpolated_structure[:1000, 1], interpolated_structure[:1000,2])

# ax.set_xlabel('X')
# ax.set_ylabel('Y')
# ax.set_zlabel('Z')
# ax.set_title('Interpolated Structure with Nucleosomes')
# plt.savefig('/mnt/raid/codes/mine/MultiEM-main/structure.pdf',format='pdf',dpi=300)
# plt.show()

# write_mmcif_chrom(interpolated_structure,'/mnt/raid/codes/mine/MultiEM-main/with_nucs.cif')