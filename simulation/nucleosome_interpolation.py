#########################################################
######## SEBASTIAN KORSAK, KRZYSZTOF BANECKI ############
#################### WARSAW 2024  #######################
#########################################################

import numpy as np
from .utils import *
from .initial_structure_tools import *
from .plots import *
from tqdm import tqdm


def makeUnit(x):
    """Normalize vector to norm 1"""
    return x / np.linalg.norm(x)


def get_perpendicular(vec):
    """Get random perpendicular vector to vec"""
    if vec[0] != 0 or vec[1] != 0:
        return np.array([vec[1], -vec[0], 0])
    else:
        return np.array([vec[2], 0, -vec[0]])


def get_perp_component(x, v):
    """Component of x orthogonal to v. Result is perpendicular to v"""
    return x - np.dot(x, v) / np.dot(v, v) * v


class NucleosomeInterpolation:
    def __init__(self, V, bw, nuc_radius=0.1, points_per_nuc=20, phi_norm=np.pi/5):
        self.V, self.bw = V, bw
        self.max_nucs_per_bead, self.nuc_points = int(np.ceil(1/(2*nuc_radius))), points_per_nuc
        self.nuc_r, self.phi_norm = nuc_radius, phi_norm
    
    def make_helix(self, r, theta, z0):
        x = r * (-np.cos(theta)+1)
        y = r * np.sin(theta)
        z = z0 * theta / theta[-1]
        return np.vstack([x, y, z]).T
    
    def min_max_scale(self,array):
        # Find the minimum and maximum values of the array
        min_val = np.min(array)
        max_val = np.max(array)
        
        # Scale the array to the interval [-1, 1]
        scaled_array = (array - min_val) / (max_val - min_val)
        return scaled_array

    def move_structure_to(self, struct, p0, p1, p2):
        """
        The structure will be placed assuming new X-axis would be defined by vector p2-p1.
        New Y-axis would be defined by part of p0-p1 vector orthogonal to p2-p1.
        """
        if p1[0]==p2[0] and p1[1]==p2[1] and p1[2]==p2[2]:
            raise(Exception("Starting point and the ending point must be different!"))
        if p1[0]==p0[0] and p1[1]==p0[1] and p1[2]==p0[2]:
            raise(Exception("Starting point and the reference point must be different!"))
        
        w_x = makeUnit(p2 - p1)
        w_y = makeUnit(get_perp_component(p1 - p0, w_x))
        w_z = makeUnit(np.cross(w_x, w_y))

        new_helix = [p1 + p[0]*w_x + p[1]*w_y + p[2]*w_z for p in struct]
        return new_helix

    def interpolate_structure_with_nucleosomes(self, mode="random"):
        """
        Interpolate the 3D structure V with nucleosomes.
        """    
        # Normalize bw_array
        bw_signal = np.log(self.bw + 1e-6)
        if not np.all(bw_signal == bw_signal[0]):
            bw_signal = self.min_max_scale(bw_signal)
        elif self.bw[0] == 0:
            bw_signal = np.zeros(bw_signal.shape)
        else:
            bw_signal = np.ones(bw_signal.shape)
        
        # Initialize interpolated structure
        interpolated_structure = []
        
        # Interpolate each segment with nucleosomes
        print('Building nucleosome structure...')
        prev_zigzag = None
        for i in tqdm(range(len(self.V) - 1)):
            # Get segment endpoints
            start_point = self.V[i]
            end_point = self.V[i + 1]
            
            # Calculate number of nucleosomes in segment
            num_nucleosomes = int(np.round(bw_signal[i] * self.max_nucs_per_bead))
            
            # Generate helices that represent nucleosomes
            interpolated_structure.append([start_point])
            if num_nucleosomes > 0:
                helices, prev_zigzag = self.single_bead_nucgenerator(start_point, end_point, num_nucleosomes, 
                                                                     prev_zigzag_vec=prev_zigzag)
                interpolated_structure.extend(helices)
            else:
                prev_zigzag = None
        interpolated_structure.append([self.V[-1]])

        # Flatten the nested list to get a list of coordinate tuples
        flattened_coords = [coord for sublist in interpolated_structure for coord in sublist]
        print('Done! You have the whole structure with nucleosomes. ;)')
        return np.array(flattened_coords)

    def single_bead_nucgenerator(self, start_point, end_point, num_nucleosomes, prev_zigzag_vec=None, turns=1.65, mode="random"):
        """
        Generate helices for nucleosomes in a segment.
        """
        # Calculate segment vector
        segment_vector = end_point - start_point
        segment_vector_norm = makeUnit(segment_vector)

        # Calculate nucleosome and linker parameters
        linker_len = self.nuc_r * 3.45
        nuc_height = self.nuc_r * 1

        # Initialize helices
        theta = np.linspace(0, turns * 2 * np.pi, self.nuc_points)
        nucleosome = self.make_helix(self.nuc_r, theta, nuc_height)
        helices = list()
        
        # Generate helices for each nucleosome
        if prev_zigzag_vec is None:
            zigzag_vec1 = makeUnit(get_perpendicular(segment_vector))
        else:
            zigzag_vec1 = makeUnit(get_perp_component(prev_zigzag_vec, segment_vector))
            if all(v==0 for v in zigzag_vec1):
                zigzag_vec1 = makeUnit(get_perpendicular(segment_vector))
        zigzag_vec2 = makeUnit(np.cross(zigzag_vec1, segment_vector))
        phi = 0
        for i in range(num_nucleosomes):
            helix_point = start_point + (i+1)/(num_nucleosomes+1)*segment_vector
            zigzag_vec = linker_len/2*(np.cos(phi)*zigzag_vec1 + np.sin(phi)*zigzag_vec2)
            if mode == "random": zigzag_vec *= np.random.uniform(0.5, 1.5)
            p1 = helix_point + zigzag_vec - nuc_height/2*segment_vector_norm
            p2 = helix_point + zigzag_vec + nuc_height/2*segment_vector_norm
            helix = self.move_structure_to(nucleosome, helix_point, p1, p2)
            helices.append(helix)
            if mode == "random":
                phi += np.pi + np.random.uniform(self.phi_norm, 2*self.phi_norm)*(np.random.randint(2)*2-1)
            else:
                phi += np.pi if i%2 == 0 else np.pi + self.phi_norm

        # zigzag vec for another segment
        zigzag_vec = np.cos(phi)*zigzag_vec1 + np.sin(phi)*zigzag_vec2

        return helices, zigzag_vec
    
def main():
    # Example data
    V = get_coordinates_cif('/home/skorsak/Data/simulation_results/GM12878_GW/MultiMM_minimized.cif')
    print('Initial granularity of structure =', len(V))
    bw = import_bw('/home/skorsak/Data/encode/ATAC-Seq/ENCSR637XSC_GM12878/ENCFF667MDI_pval.bigWig', len(V))  # Mock self.signal array

    # Interpolate structure with nucleosomes
    nuc_interpol = NucleosomeInterpolation(V, bw)
    iV = nuc_interpol.interpolate_structure_with_nucleosomes(mode="random")
    print('Final Length of Nucleosome Interpolated Structure:', len(iV))
    viz_structure(iV, r=0.1)
