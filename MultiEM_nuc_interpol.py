#########################################################
######## SEBASTIAN KORSAK, WARSAW 2024 ##################
######## SPECIAL THANKS TO KRZYSTOF BANECKI #############
#########################################################

import numpy as np
from MultiEM_utils import *
from MultiEM_init_tools import *
from MultiEM_plots import *
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


class NucleosomeInterpolation:
    def __init__(self,V,bw,max_nucs_per_bead=4,zig_zag_displacement=0.2,points_per_nuc=20,phi_norm=np.pi/5):
        self.V, self.bw = V, bw
        self.max_nucs_per_bead, self.nuc_points = max_nucs_per_bead, points_per_nuc
        self.zigzag_d, self.phi_norm = zig_zag_displacement, phi_norm

    def amplify(self,structure, scale=10):
        return [(s[0]*scale, s[1]*scale, s[2]*scale) for s in structure]
    
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

    def move_structure_to(self, struct, p0, p1, p2, x0=np.array([None])):
        """
        The structure will be placed assuming new X-axis would be defined by vector p2-p1.
        New Y-axis would be defined by part of p0-p1 vector orthogonal to p2-p1.
        """
        if p1[0]==p2[0] and p1[1]==p2[1] and p1[2]==p2[2]:
            raise(Exception("Starting point and the ending points must be different!"))
        w_x = p2 - p1
        v_01 = p1 - p0
        w_y = v_01 - np.dot(w_x, v_01)/np.linalg.norm(w_x)**2 * w_x
        w_z = np.cross(w_x, w_y)

        w_x = makeUnit(w_x) if np.linalg.norm(w_x)>0 else w_x
        w_y = makeUnit(w_y) if np.linalg.norm(w_y)>0 else w_y
        w_z = makeUnit(w_z) if np.linalg.norm(w_z)>0 else w_z

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
            interpolated_structure.append([start_point])
            if num_nucleosomes>0:
                helices = self.single_bead_nucgenerator(start_point, end_point, num_nucleosomes)
                interpolated_structure.extend(helices)
            else:
                helices = [[(tuple((self.V[i]+self.V[i+1])/2))]]
                interpolated_structure.extend(helices)
        interpolated_structure.append([self.V[-1]])

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
        helix_radius = 0.1 # np.linalg.norm(segment_vector)/np.pi
        helix_axis = makeUnit(segment_vector)
        helix_height = helix_radius/0.965

        # Initialize helices
        theta = np.linspace(0, turns* 2 * np.pi, self.nuc_points)
        helices = list()
        
        # Generate helices for each nucleosome
        zigzag_vec1 = makeUnit(get_perpendicular(segment_vector))
        zigzag_vec2 = makeUnit(np.cross(zigzag_vec1, segment_vector))
        phi = 0
        for i in range(num_nucleosomes):
            helix_point = start_point + (i+1)/(num_nucleosomes+1)*segment_vector
            zigzag_vec = self.zigzag_d*(np.cos(phi)*zigzag_vec1 + np.sin(phi)*zigzag_vec2)
            p1 = helix_point + zigzag_vec - helix_radius/2*helix_axis
            p2 = helix_point + zigzag_vec + helix_radius/2*helix_axis
            h = self.make_helix(helix_radius, theta, helix_height)
            helix = self.move_structure_to(h, helix_point, p1, p2)
            
            helices.append(helix)
            phi += np.pi if i%2 == 0 else np.pi + self.phi_norm

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