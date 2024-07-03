#########################################################
######## SEBASTIAN KORSAK, WARSAW 2024 ##################
######## SPECIAL THANKS TO KRZYSTOF BANECKI #############
#########################################################

import numpy as np
from MultiEM_utils import *
from MultiEM_init_tools import *
from MultiEM_plots import *
from tqdm import tqdm
import torch

class NucleosomeInterpolation:
    def __init__(self,V,bw,max_nucs_per_bead=4,zig_zag_displacement=np.pi/6,points_per_nuc=20,phi_norm=np.pi):
        self.V, self.bw = V, bw
        self.max_nucs_per_bead, self.nuc_points = max_nucs_per_bead, points_per_nuc
        self.zigzag_d, self.phi_norm = zig_zag_displacement, phi_norm

    def amplify(self,structure, scale=10):
        return [(s[0]*scale, s[1]*scale, s[2]*scale) for s in structure]
    
    def make_helix(self,r,theta,z0):
        x = r * (np.cos(theta)-1)*self.sign
        y = r * np.sin(theta)
        z = z0 * theta / (2 * np.pi)
        return np.vstack([x, y, z]).T
    
    def min_max_scale(self,array):
        # Find the minimum and maximum values of the array
        min_val = np.min(array)
        max_val = np.max(array)
        
        # Scale the array to the interval [-1, 1]
        scaled_array = (array - min_val) / (max_val - min_val)
        return scaled_array

    def orthonormalize(self,vectors):
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

    def move_structure_to(self, struct, p1, p2, x0=np.array([None])):
        if p1[0]==p2[0] and p1[1]==p2[1] and p1[2]==p2[2]:
            raise(Exception("Starting point and the ending points must be different!"))
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
        A_norm = np.array(self.orthonormalize(A))
        w_x = list(A_norm[:,1])
        w_y = list(A_norm[:,2])
        w_z = list(A_norm[:,0])
        if x0.all()==None: x0=p1
        new_helix = []
        for p in struct:
            new_helix.append((x0[0]+p[0]*w_x[0]+p[1]*w_y[0]+p[2]*w_z[0],
                            x0[1]+p[0]*w_x[1]+p[1]*w_y[1]+p[2]*w_z[1],
                            x0[2]+p[0]*w_x[2]+p[1]*w_y[2]+p[2]*w_z[2]))
        return new_helix

    def interpolate_structure_with_nucleosomes(self):
        """
        Interpolate the 3D structure V with nucleosomes.
        """    
        # Normalize bw_array
        self.bw[self.bw>np.mean(self.bw)+6*np.std(self.bw)] = np.mean(self.bw)+6*np.std(self.bw)
        norm_bw_array = self.min_max_scale(self.bw)
        
        # Initialize interpolated structure
        interpolated_structure = []
        
        # Interpolate each segment with nucleosomes
        self.sign, self.phi = 1, 0
        print('Building nucleosome structure...')
        for i in tqdm(range(len(self.V) - 1)):
            # Get segment endpoints
            start_point = self.V[i]
            end_point = self.V[i + 1]
            
            # Calculate segment length
            segment_length = np.linalg.norm(end_point - start_point)
            
            # Calculate number of nucleosomes in segment
            num_nucleosomes = int(np.round(norm_bw_array[i] * self.max_nucs_per_bead))
            
            # Generate helices that represent nucleosomes
            if num_nucleosomes>0:
                helices = self.single_bead_nucgenerator(start_point, end_point, num_nucleosomes)
                interpolated_structure.extend(helices)
            else:
                helices = [[(tuple((self.V[i]+self.V[i+1])/2))]]
                interpolated_structure.extend(helices)
        
        # Flatten the nested list to get a list of coordinate tuples
        flattened_coords = [coord for sublist in interpolated_structure for coord in sublist]
        print('Done! You have the whole structure with nucleosomes. ;)')
        return np.array(flattened_coords)

    def single_bead_nucgenerator(self, start_point, end_point, num_nucleosomes ,turns=1.6):
        """
        Generate helices for nucleosomes in a segment.
        """
        # Calculate segment vector
        segment_vector = end_point - start_point

        # Calculate helix parameters
        helix_radius = np.linalg.norm(segment_vector)/np.pi
        helix_axis = segment_vector / np.linalg.norm(segment_vector)
        helix_height = helix_radius/0.965

        # Initialize helices
        theta = np.linspace(0, turns* 2 * np.pi, self.nuc_points)
        helices = list()

        helix_points = [start_point]
        for i in range(num_nucleosomes):
            helix_points.append(start_point +  (i + 1.0/num_nucleosomes) * helix_axis* helix_height)

        # Generate helices for each nucleosome
        for i in range(num_nucleosomes):
            if self.sign==1: 
                self.phi+=(self.phi_norm/num_nucleosomes)
                self.phi = self.phi%(np.pi/2)
            zigzag_displacement = -self.sign*self.zigzag_d
            add_vec = [zigzag_displacement*np.cos(self.phi),zigzag_displacement*np.sin(self.phi),0]
            h = self.make_helix(helix_radius, theta, helix_height)+add_vec
            helix = self.move_structure_to(h,helix_points[i],helix_points[i+1])
            self.sign*=-1

            # Generate helix coordinates     
            helices.append(helix)

        return helices
    
def main():
    # Example data
    V = get_coordinates_cif('/home/skorsak/Data/simulation_results/GM12878_GW/MultiEM_minimized.cif')
    print('Initial granularity of structure =',len(V))
    bw = import_bw('/home/skorsak/Data/encode/ATAC-Seq/ENCSR637XSC_GM12878/ENCFF667MDI_pval.bigWig',len(V))  # Mock self.signal array

    # Interpolate structure with nucleosomes
    nuc_interpol = NucleosomeInterpolation(V[:1000],bw[:1000])
    iV = nuc_interpol.interpolate_structure_with_nucleosomes()
    print('Final Length of Nucleosome Interpolated Structure:',len(iV))
    viz_structure(iV,r=0.1)