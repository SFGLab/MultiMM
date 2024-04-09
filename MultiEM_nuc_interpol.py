import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import plotly.graph_objects as go
from MultiEM_utils import *
from MultiEM_init_tools import *
from tqdm import tqdm
import torch

def amplify(structure, scale=10):
    return [(s[0]*scale, s[1]*scale, s[2]*scale) for s in structure]

def get_orig_helix(n=40, r=3.25, a=0.2):
    points = []
    for t in np.linspace(0, 4*np.pi, n):
        points.append((r*np.cos(t)-r, r*np.sin(t), a*t))
    return points

def orthonormalize(vectors):
    """
        Orthonormalizes the vectors using gram schmidt procedure.

        Parameters:
            vectors: torch tensor, size (dimension, n_vectors)
                    they must be linearly independant
        Returns:
            orthonormalized_vectors: torch tensor, size (dimension, n_vectors)
    """
    assert (vectors.size(1) <= vectors.size(0)), 'number of vectors must be smaller or equal to the dimension'
    orthonormalized_vectors = torch.zeros_like(vectors)
    orthonormalized_vectors[:, 0] = vectors[:, 0] / torch.norm(vectors[:, 0], p=2)

    for i in range(1, orthonormalized_vectors.size(1)):
        vector = vectors[:, i]
        V = orthonormalized_vectors[:, :i]
        PV_vector= torch.mv(V, torch.mv(V.t(), vector))
        orthonormalized_vectors[:, i] = (vector - PV_vector) / torch.norm(vector - PV_vector, p=2)

    return orthonormalized_vectors

def move_structure_to(struct, p1, p2, x0=np.array([None])):
    if p1[0]==p2[0] and p1[1]==p2[1] and p1[2]==p2[2]:
        raise(Exception("Starting point and the ending points must be different!"))
    v_x = [1,0,0]
    v_y = [0,1,0]
    v_z = [0,0,1]
    w_z = [p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]]
    if w_z[0] != 0:
        w_x = [-w_z[1]/w_z[0], 1, 0]
        w_y = [-w_z[2]/w_z[0], 0, 1]
    elif w_z[1] != 0:
        w_x = [1, -w_z[0]/w_z[1], 0]
        w_y = [0, -w_z[2]/w_z[1], 1]
    else:
        w_x = [1, 0, -w_z[0]/w_z[2]]
        w_y = [0, 1, -w_z[1]/w_z[2]]
    A = torch.transpose(torch.tensor([w_z, w_x, w_y]),0,1)
    A_norm = np.array(orthonormalize(A))
    w_x = list(A_norm[:,1])
    w_y = list(A_norm[:,2])
    w_z = list(A_norm[:,0])
    if x0.all()==None: x0=p1
    new_helix = []
    for p in struct:
        new_helix.append((x0[0]+p[0]*w_x[0]+p[1]*w_y[0]+p[2]*w_z[0],
                         x0[1]+p[0]*w_x[1]+p[1]*w_y[1]+p[2]*w_z[1],
                         x0[2]+p[0]*w_x[2]+p[1]*w_y[2]+p[2]*w_z[2]))
    return np.array(new_helix)

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
    print('Done you have the whole structure with nucleosomes')
    
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

def make_helix(r,theta,z0,sign):
    x = r * (np.cos(theta)-1) if sign==1 else r * (-np.cos(theta)+1) # to avoid DNA crossing
    y = r * np.sin(theta)
    z = z0 * theta / (2 * np.pi)
    return np.vstack([x, y, z]).T

def generate_nucleosome_helices(start_point, end_point, num_nucleosomes, turns=1.6, sign=1):
    """
    Generate helices for nucleosomes in a segment.
    """
    # Calculate segment vector
    segment_vector = end_point - start_point

    # Calculate helix parameters
    helix_radius = np.linalg.norm(segment_vector)/np.pi
    helix_axis = segment_vector / np.linalg.norm(segment_vector)
    helix_height = 0.2

    # Initialize helices
    theta = np.linspace(0, turns* 2 * np.pi, 10)
    helices = list()

    helix_points = [start_point]
    for i in range(num_nucleosomes):
        helix_points.append(start_point +  (i + 1.0/num_nucleosomes) * helix_axis* helix_height)

    # Generate helices for each nucleosome
    for i in range(num_nucleosomes):
        # Calculate helix center
        phi = np.pi/2
        
        # Calculate helix angle
        zigzag_angle = sign*np.pi/3
        helix = make_helix(helix_radius,theta,helix_height,sign)
        helix = move_structure_to(helix,helix_points[i],helix_points[i+1])+[0,0,zigzag_angle]
        sign*=-1

        # Generate helix coordinates        
        helices.append(helix)

    return helices, sign

# Example data
V = get_coordinates_cif('/mnt/raid/codes/mine/MultiEM-main/Trios_ensembles/CHS_d/MultiEM_minimized.cif')
N = len(V)
bw_array = np.random.rand(N)  # Mock signal array

# Interpolate structure with nucleosomes
iV = interpolate_structure_with_nucleosomes(V, bw_array)
points = np.arange(0,len(iV))

fig = go.Figure(data=go.Scatter3d(
    x=iV[:1000,0], y=iV[:1000,1], z=iV[:1000,2],
    marker=dict(
        size=2,
        color=points[:1000],
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