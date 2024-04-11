#########################################################
######## SEBASTIAN KORSAK, WARSAW 2024 ##################
######## SPECIAL THANKS TO KRZYSTOF BANECKI #############
#########################################################

import numpy as np
from MultiEM_utils import *
from MultiEM_init_tools import *
from tqdm import tqdm
import pyvista as pv
import torch

def import_bw(bw_path,N_beads,coords=None,chrom=None,viz=False,binary=False,path='',sign=1,norm=False):
    '''
    Imports .BigWig data and outputs compartments.

    It assumes that higher signal coresponds to B compartment.

    In case that you would like to switch the sign then add flag sign=-1.
    '''
    # Open file
    bw = pyBigWig.open(bw_path)
    chroms_set = np.fromiter(bw.chroms().keys(),dtype='S20')
    n_chroms=24 if b'chrY' in chroms_set else 23
    print('Number of chromosomes:',n_chroms)

    # Compute the total length of chromosomes
    if chrom==None:
        chrom_length = 0
        lengths = list()
        for i in range(n_chroms):
            chrom_length += bw.chroms(chrs[i])
            lengths.append(bw.chroms(chrs[i]))
        lengths = np.array(lengths)
        resolution = chrom_length//(2*N_beads)
        polymer_lengths = lengths//resolution
        np.save(path+'chrom_lengths.npy',polymer_lengths)

    # Import the downgraded signal
    print('Importing bw signal...')
    if chrom==None:
        genomewide_signal = list()
        for i in tqdm(range(n_chroms)):
            signal = bw.values(chrs[i],0,-1, numpy=True)
            signal = np.nan_to_num(signal, copy=True, nan=0.0, posinf=0.0, neginf=0.0)
            genomewide_signal.append(compute_averages(signal, polymer_lengths[i]))
        genomewide_signal = np.concatenate(genomewide_signal)
    else:
        genomewide_signal = bw.values(chrom,coords[0],coords[1], numpy=True)
    bw.close()

    genomewide_signal = compute_averages(genomewide_signal,N_beads)
    if norm: genomewide_signal = (genomewide_signal-np.mean(genomewide_signal))/np.std(genomewide_signal)
    
    # Transform signal to binary or adjuct it to have zero mean
    if binary:
        genomewide_signal[genomewide_signal>0] = -1
        genomewide_signal[genomewide_signal<=0] = 1
        
        # Subtitute zeros with random spin states
        mask = genomewide_signal==0
        n_zeros = np.count_nonzero(mask)
        nums = np.array(rd.choices([-1,1],k=n_zeros))
        genomewide_signal[mask] = nums
    
    print('Done!\n')

    # Plotting
    if viz:
        figure(figsize=(25, 5), dpi=100)
        xax = np.arange(len(genomewide_signal))
        plt.fill_between(xax,genomewide_signal, where=(genomewide_signal>np.mean(genomewide_signal)),alpha=0.50,color='purple')
        plt.fill_between(xax,genomewide_signal, where=(genomewide_signal<np.mean(genomewide_signal)),alpha=0.50,color='orange')
        lines = np.cumsum(polymer_lengths//2)
        for i in range(n_chroms):
            plt.axvline(x=lines[i], color='b')
        plt.xlabel('Genomic Distance',fontsize=16)
        plt.ylabel('BW signal renormalized',fontsize=16)
        plt.ylim((np.mean(genomewide_signal)-np.std(genomewide_signal),np.mean(genomewide_signal)+np.std(genomewide_signal)))
        plt.grid()
        plt.show()

    np.save(path+'signal.npy',genomewide_signal)
    
    return sign*genomewide_signal

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

def interpolate_structure_with_nucleosomes(V, bw_array, max_Nnuc=3):
    """
    Interpolate the 3D structure V with nucleosomes.
    """
    # Calculate mean of bw_array
    mean_signal = np.mean(bw_array)
    std_signal = np.std(bw_array)
    
    # Normalize bw_array
    norm_bw_array = -(bw_array - mean_signal)/std_signal
    
    # Initialize interpolated structure
    interpolated_structure = []
    
    # Interpolate each segment with nucleosomes
    sign, phi = 1, 0 
    print('Building nucleosome structure...')
    for i in tqdm(range(len(V) - 1)):
        # Get segment endpoints
        start_point = V[i]
        end_point = V[i + 1]
        
        # Calculate segment length
        segment_length = np.linalg.norm(end_point - start_point)
        
        # Calculate number of nucleosomes in segment
        num_nucleosomes = int(np.round(norm_bw_array[i] * max_Nnuc))
        
        # Generate helices for nucleosomes
        helices, sign, phi = generate_nucleosome_helices(start_point, end_point, num_nucleosomes, phi, 1.6, sign)
        
        # Append helices to interpolated structure
        interpolated_structure.extend(helices)
    print('Done! You have the whole structure with nucleosomes. ;)')
    
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

def generate_nucleosome_helices(start_point, end_point, num_nucleosomes, phi ,turns=1.6, sign=1):
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
        # Calculate helix angle
        if sign==1: phi+=np.pi/3
        zigzag_displacement = sign*np.pi/6
        zz_add = np.array([zigzag_displacement*np.sin(phi),zigzag_displacement*np.cos(phi),0])
        helix = make_helix(helix_radius,theta,helix_height,sign)+zz_add
        helix = move_structure_to(helix,helix_points[i],helix_points[i+1])
        sign*=-1

        # Generate helix coordinates        
        helices.append(helix)

    return helices, sign, phi

def polyline_from_points(points):
    poly = pv.PolyData()
    poly.points = points
    the_cell = np.arange(0, len(points), dtype=np.int_)
    the_cell = np.insert(the_cell, 0, len(points))
    poly.lines = the_cell
    return poly

def viz_structure(V, colors=None):
    polyline = polyline_from_points(V)
    polyline["scalars"] = np.arange(polyline.n_points)

    if colors is not None:
        cmap = "coolwarm"
        color_values = (colors - np.min(colors)) / (np.max(colors) - np.min(colors))  # Normalize colors
        polyline["colors"] = color_values  # Set colors as point scalars
        polymer = polyline.tube(radius=0.3)
        polymer.plot(smooth_shading=True, cmap=cmap, scalars="colors")
    else:
        polymer = polyline.tube(radius=0.1)
        polymer.plot(smooth_shading=True)

def viz_chroms(cif_path,chrom_ends_path):
    chrom_ends = np.load(chrom_ends_path)
    V = get_coordinates_cif(cif_path)
    N = len(V)
    chroms = np.zeros(N)
    for chrom, i in enumerate(range(len(chrom_ends)-1)):
        start, end = chrom_ends[i], chrom_ends[i+1]
        chroms[start:end] = chrom
    viz_structure(V,chroms[:len(V)])
    
def main():
    # Example data
    V = get_coordinates_cif('/home/skorsak/Templates/MultiEM-main/test/MultiEM_minimized.cif')
    N = len(V)
    bw_array = import_bw('/home/skorsak/Documents/data/encode/ATAC-Seq/ENCSR637XSC_GM12878/ENCFF667MDI_pval.bigWig',N)  # Mock signal array

    # Interpolate structure with nucleosomes
    iV = interpolate_structure_with_nucleosomes(V, bw_array)
    points = np.arange(0,len(iV))
    print('Final Length of Nucleosome Interpolated Structure:',len(iV))

    viz_structure(iV)