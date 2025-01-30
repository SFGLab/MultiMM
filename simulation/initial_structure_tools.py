import numpy as np
import pandas as pd
from scipy import interpolate
import matplotlib.pyplot as plt
from hilbertcurve.hilbertcurve import HilbertCurve
import warnings
from .utils import *
warnings.simplefilter(action='ignore', category=FutureWarning)

subcomp_dict={-2:'B1',-1:'B2',0:'O',1:'A1',2:'A2'}

def encode_pure(digits, value):
    """encodes value using the given digits"""
    assert value >= 0
    if value == 0:
        return digits[0]
    n = len(digits)
    result = []
    while value != 0:
        rest = value // n
        result.append(digits[value - rest * n])
        value = rest
    result.reverse()
    return "".join(result)


def decode_pure(digits_values, s):
    """decodes the string s using the digit, value associations for each character"""
    result = 0
    n = len(digits_values)
    for c in s:
        result *= n
        result += digits_values[c]
    return result


def hy36encode(width, value):
    """encodes value as base-10/upper-case base-36/lower-case base-36 hybrid"""
    i = value
    if i >= 1 - 10 ** (width - 1):
        if i < 10 ** width:
            return ("%%%dd" % width) % i
        i -= 10 ** width
        if i < 26 * 36 ** (width - 1):
            i += 10 * 36 ** (width - 1)
            return encode_pure(digits_upper, i)
        i -= 26 * 36 ** (width - 1)
        if i < 26 * 36 ** (width - 1):
            i += 10 * 36 ** (width - 1)
            return encode_pure(digits_lower, i)
    raise ValueError("value out of range.")


def hy36decode(width, s):
    """decodes base-10/upper-case base-36/lower-case base-36 hybrid"""
    if len(s) == width:
        f = s[0]
        if f == "-" or f == " " or f.isdigit():
            try:
                return int(s)
            except ValueError:
                pass
            if s == " " * width:
                return 0
        elif f in digits_upper_values:
            try:
                return decode_pure(
                    digits_values=digits_upper_values, s=s) - 10 * 36 ** (width - 1) + 10 ** width
            except KeyError:
                pass
        elif f in digits_lower_values:
            try:
                return decode_pure(
                    digits_values=digits_lower_values, s=s) + 16 * 36 ** (width - 1) + 10 ** width
            except KeyError:
                pass
    raise ValueError("invalid number literal.")

def find_element_indexes(arr, elem):
    indexes = np.where(arr == elem)[0]
    ranges = []
    start_index = indexes[0]
    
    for i in range(1, len(indexes)):
        if indexes[i] != indexes[i-1] + 1:
            if start_index == indexes[i-1]:
                ranges.append(str(start_index))
            else:
                ranges.append(f"{start_index}-{indexes[i-1]}")
            start_index = indexes[i]

    # Handle the last range
    if start_index == indexes[-1]:
        ranges.append(str(start_index))
    else:
        ranges.append(f"{start_index}-{indexes[-1]}")

    return ', '.join(ranges)

def write_cmm(comps,name='MultiMM_compartment_colors.cmd'):
    comp_old = 2
    counter, start = 0, 0
    comp_dict = {-2:'#bf0020', -1:'#ba5062', 1:'#4e4c87',2:'#181385',0:'#fafcfc'}
    spins = np.unique(comps)
    lines = []

    for s in spins:
        lines.append('color '+comp_dict[int(s)]+' :')
    
    for i, s in enumerate(spins):
        positions = find_element_indexes(comps, s)
        lines[i] = lines[i]+positions
    
    content=''
    for i in range(len(lines)):
        content+=lines[i]+'\n'
    
    with open(name,'w') as fp:
        fp.write(content)
    fp.close()

def read_compartments(file,ch,reg,res,binary=True):
    file = pd.read_csv('/mnt/raid/data/compartments/RAO_GM12878_subcomp_hg38.bed',header=None,sep='\t')
    file = file[(file[1]>reg[0])&(file[2]<reg[1])&(file[0]==ch)].reset_index(drop=True)

    coords, comp_list = np.zeros((len(file),2)), list()
    for i in range(len(file)):
        coords[i,0] = int(file[1][i]//res)
        coords[i,1] = int(file[2][i]//res)
        if binary:
            comp_list.append(file[3][i][0])
        else:
            comp_list.append(file[3][i])
    return coords, comp_list

def generate_hilbert_curve(n_points,p=8,n=3,displacement_sigma=0.1,scale=6,viz=False,add_noise=False):
    hilbert_curve = HilbertCurve(p, n)

    distances = list(range(n_points))
    points = np.array(hilbert_curve.points_from_distances(distances))
    if add_noise:
        displacement = np.random.normal(loc=0.0, scale=displacement_sigma, size=n_points*3).reshape(n_points,3)
        V_interpol = V_interpol + displacement
    
    return points

def polymer_circle(n: int, z_stretch: float = 1.0, radius: float = 5.0) -> np.ndarray:
    points = []
    angle_increment = 360 / float(n)
    radius = 1 / (2 * np.sin(np.radians(angle_increment) / 2.)) if radius==None else radius
    z_stretch = z_stretch / n
    z = 0
    for i in range(n):
        x = radius * np.cos(angle_increment * i * np.pi / 180)
        y = radius * np.sin(angle_increment * i * np.pi / 180)
        if z_stretch != 0:
            z += z_stretch
        points.append((x, y, z))
    points = np.array(points)
    return points

def helix_structure(N_beads, radius=1, pitch=2):
    theta = np.linspace(0, 4 * np.pi, N_beads)  # 2 full turns
    x = radius * np.cos(theta)
    y = radius * np.sin(theta)
    z = np.linspace(0, pitch * N_beads, N_beads)
    V = np.column_stack((x, y, z))
    return V

def spiral_structure(N_beads, initial_radius=1, pitch=1, growth_factor=0.05):
    theta = np.linspace(0, 4 * np.pi, N_beads)
    radius = initial_radius + growth_factor * np.arange(N_beads)
    x = radius * np.cos(theta)
    y = radius * np.sin(theta)
    z = np.linspace(0, pitch * N_beads, N_beads)
    V = np.column_stack((x, y, z))
    return V

def sphere_surface_structure(N_beads, radius=1):
    phi = np.random.uniform(0, 2 * np.pi, N_beads)
    costheta = np.random.uniform(-1, 1, N_beads)
    u = np.random.uniform(0, 1, N_beads)
    
    theta = np.arccos(costheta)
    r = radius * u ** (1/3)

    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    
    V = np.column_stack((x, y, z))
    return V

def confined_random_walk(N_beads, box_size=5):
    V = np.zeros((N_beads, 3))
    for i in range(1, N_beads):
        step = np.random.choice([-1, 1], size=3)  # Random step in x, y, z
        V[i] = V[i-1] + step
        # Keep the points within a confined box (folding in)
        V[i] = np.clip(V[i], -box_size, box_size)
    return V

def trefoil_knot_structure(N_beads, scale=5):
    t = np.linspace(0, 2 * np.pi, N_beads)
    x = scale * (np.sin(t) + 2 * np.sin(2 * t))
    y = scale * (np.cos(t) - 2 * np.cos(2 * t))
    z = -scale * np.sin(3 * t)
    
    V = np.column_stack((x, y, z))
    return V

def random_walk_structure(N_beads, step_size=1):
    # Initialize the structure array
    V = np.zeros((N_beads, 3))
    
    # Loop over each bead, starting from the second one
    for i in range(1, N_beads):
        # Randomly pick a direction for each step
        step_direction = np.random.normal(size=3)
        step_direction /= np.linalg.norm(step_direction)  # Normalize to make it unit length
        
        # Move the current bead from the last bead by a fixed step size
        V[i] = V[i-1] + step_size * step_direction
    
    return V

def compute_init_struct(N_beads,mode='rw',scale=5):
    match mode:
        case 'rw':
            return random_walk_structure(N_beads)
        case 'confined_rw':
            return confined_random_walk(N_beads)
        case 'knot':
            return trefoil_knot_structure(N_beads)
        case 'self_avoiding_rw':
            return self_avoiding_random_walk(N_beads)
        case 'circle':
            return polymer_circle(N_beads,50,5)
        case 'helix':
            return helix_structure(N_beads)
        case 'spiral':
            return spiral_structure(N_beads)
        case 'sphere':
            return sphere_surface_structure(N_beads)
        case 'hilbert':
            return generate_hilbert_curve(N_beads,scale=scale)
        case _:
            return IndentationError('Invalid option for initial structure.')

def build_init_mmcif(n_dna,chrom_ends,psf=True,path='',curve='hilbert',scale=5):
    # Define the initial coordinates of histones and the structure of DNA
    dna_points = compute_init_struct(n_dna,mode=curve,scale=scale)
    # Write the positions in .mmcif file
    atoms = ''
    
    ## DNA beads
    for i in range(n_dna):
        chain_idx = np.searchsorted(chrom_ends,i)
        if i in chrom_ends: chain_idx+=1
        [atom_type,res_name, atom_name, cl] = ['HETATM','ALB', 'CB', chr(65+chain_idx)] if (i in chrom_ends) | (i in chrom_ends-1) else ['ATOM','ALA', 'CA', chr(65+chain_idx)]
        atoms += ('{0:} {1:} {2:} {3:} {4:} {5:} {6:} {7:} {8:} '
                  '{9:} {10:.3f} {11:.3f} {12:.3f}\n'.format(atom_type, i+1, 'D', atom_name,\
                                                             '.', res_name, cl, chain_idx, i+1, '?',\
                                                             dna_points[i,0], dna_points[i,1], dna_points[i,2]))
    
    # Write connections
    connects = ''
    
    for i in range(n_dna-1):
        chain_idx = np.searchsorted(chrom_ends,i)
        if i in chrom_ends: chain_idx+=1
        if i not in chrom_ends-1:
            [atom_type1,res_name1, atom_name1, cl1] = ['HETATM','ALB', 'CB', chr(65+chain_idx)] if (i in chrom_ends) else ['ATOM','ALA', 'CA', chr(65+chain_idx)]
            [atom_type2,res_name2, atom_name2, cl2] = ['HETATM','ALB', 'CB', chr(65+chain_idx)] if (i+1 in chrom_ends-1) else ['ATOM','ALA', 'CA', chr(65+chain_idx)]
            connects += f'D{i+1} covale {res_name1} {cl1} {i+1} {atom_name1} {res_name2} {cl2} {i+2} {atom_name2}\n'

    # Save files
    mmcif_file_name = path+'MultiMM_init.cif'
    atomhead = mmcif_atomhead()
    conhead = mmcif_connecthead()
    mmcif_file_content = atomhead+atoms+'\n'+conhead+connects

    with open(mmcif_file_name, 'w') as f:
        f.write(mmcif_file_content)

    if psf:
        generate_psf(n_dna,path+'MultiMM.psf')

    print("File {} saved...".format(mmcif_file_name))

def write_mmcif(coords,chrom_ends,path):
    # Write the positions in .mmcif file
    atoms = ''
    
    ## DNA beads
    for i in range(len(coords)):
        chain_idx = np.searchsorted(chrom_ends,i)
        if i in chrom_ends: chain_idx+=1
        [res_name, atom_name, cl] = ['ALB', 'CB', chr(65+chain_idx)] if (i in chrom_ends) | (i in chrom_ends-1) else ['ALA', 'CA', chr(65+chain_idx)]
        atoms += ('{0:} {1:} {2:} {3:} {4:} {5:} {6:} {7:} {8:} '
                  '{9:} {10:.3f} {11:.3f} {12:.3f}\n'.format('ATOM', i+1, 'D', atom_name,\
                                                             '.', res_name, cl, chain_idx, i+1, '?',\
                                                             coords[i,0], coords[i,1], coords[i,2]))
    
    # Write connections
    connects = ''
    
    for i in range(len(coords)-1):
        chain_idx = np.searchsorted(chrom_ends,i)
        if i in chrom_ends: chain_idx+=1
        if i not in chrom_ends-1:
            [res_name1, atom_name1, cl1] = ['ALB', 'CB', chr(65+chain_idx)] if (i in chrom_ends) else ['ALA', 'CA', chr(65+chain_idx)]
            [res_name2, atom_name2, cl2] = ['ALB', 'CB', chr(65+chain_idx)] if (i+1 in chrom_ends-1) else ['ALA', 'CA', chr(65+chain_idx)]
            connects += f'D{i+1} covale {res_name1} {cl1} {i+1} {atom_name1} {res_name2} {cl2} {i+2} {atom_name2}\n'

    # Save files
    ## .pdb
    atomhead = mmcif_atomhead()
    conhead = mmcif_connecthead()
    mmcif_file_content = atomhead+atoms#+conhead+connects
    
    f = open(path, "w")
    f.write(mmcif_file_content)
    f.close()

def write_mmcif_chrom(coords,path):
    # Write the positions in .mmcif file
    atoms = ''
    
    ## DNA beads
    for i in range(len(coords)):
        [res_name, atom_name, cl] = ['ALA', 'CA', 'A'] if (i!=0 and i!=len(coords)-1) else ['ALB', 'CA', 'A']
        atoms += ('{0:} {1:} {2:} {3:} {4:} {5:} {6:} {7:} {8:} '
                  '{9:} {10:.3f} {11:.3f} {12:.3f}\n'.format('ATOM', i+1, 'D', atom_name,\
                                                             '.', res_name, cl, 1, i+1, '?',\
                                                             coords[i,0], coords[i,1], coords[i,2]))
    
    # Write connections
    connects = ''
    
    for i in range(len(coords)-1):
        [res_name1, atom_name1, cl1] = ['ALA', 'CA', 'A'] if (i!=0 and i!=(len(coords)-1)) else ['ALB', 'CA', 'A']
        [res_name2, atom_name2, cl2] = ['ALA', 'CA', 'A'] if ((i+1)!=0 and (i+1)!=(len(coords)-1)) else ['ALB', 'CA', 'A']
        connects += f'D{i+1} covale {res_name1} {cl1} {i+1} {atom_name1} {res_name2} {cl2} {i+2} {atom_name2}\n'

    # Save files
    ## .pdb
    atomhead = mmcif_atomhead()
    conhead = mmcif_connecthead()
    mmcif_file_content = atomhead+atoms+conhead+connects
    
    f = open(path, "w")
    f.write(mmcif_file_content)
    f.close()

def generate_psf(n: int, file_name='stochastic_LE.psf', title="No title provided"):
    """
    Saves PSF file. Useful for trajectories in DCD file format.
    :param n: number of points
    :param file_name: PSF file name
    :param title: Human readable string. Required in PSF file.
    :return: List with string records of PSF file.
    """
    assert len(title) < 40, "provided title in psf file is too long."
    # noinspection PyListCreation
    lines = ['PSF CMAP\n']
    lines.append('\n')
    lines.append('      1 !NTITLE\n')
    lines.append('REMARKS {}\n'.format(title))
    lines.append('\n')
    lines.append('{:>8} !NATOM\n'.format(n))
    for k in range(1, n + 1):
        lines.append('{:>8} BEAD {:<5} ALA  CA   A      0.000000        1.00 0           0\n'.format(k, k))
    lines.append('\n')
    lines.append('{:>8} !NBOND: bonds\n'.format(n - 1))
    for i in range(1, n):
        lines.append('{:>8}{:>8}\n'.format(i, i + 1))
    with open(file_name, 'w') as f:
        f.writelines(lines)

def mmcif_atomhead():
    head = """data_MultiMM
# 
_entry.id MultiMM
# 
_audit_conform.dict_name       mmcif_pdbx.dic 
_audit_conform.dict_version    5.296 
_audit_conform.dict_location   http://mmcif.pdb.org/dictionaries/ascii/mmcif_pdbx.dic 
# ----------- ATOMS ----------------
loop_
_atom_site.group_PDB 
_atom_site.id 
_atom_site.type_symbol 
_atom_site.label_atom_id 
_atom_site.label_alt_id 
_atom_site.label_comp_id 
_atom_site.label_asym_id 
_atom_site.label_entity_id 
_atom_site.label_seq_id 
_atom_site.pdbx_PDB_ins_code 
_atom_site.Cartn_x 
_atom_site.Cartn_y 
_atom_site.Cartn_z
"""
    return head

def mmcif_connecthead():
    content = """#
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_label_atom_id
"""
    return content

def get_coordinates_cif(file):
    '''
    It returns the corrdinate matrix V (N,3) of a .pdb file.
    The main problem of this function is that coordiantes are not always in 
    the same column position of a .pdb file. Do changes appropriatelly,
    in case that the data aren't stored correctly. 
    
    Input:
    file (str): the path of the .cif file.
    
    Output:
    V (np.array): the matrix of coordinates
    '''
    V = list()
    
    with open(file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("ATOM"):
                columns = line.split()
                x = eval(columns[10])
                y = eval(columns[11])
                z = eval(columns[12])
                V.append([x, y, z])
    
    return np.array(V)

def get_df_cif(file):
    '''
    It returns the corrdinate matrix V (N,3) of a .pdb file.
    The main problem of this function is that coordiantes are not always in 
    the same column position of a .pdb file. Do changes appropriatelly,
    in case that the data aren't stored correctly. 
    
    Input:
    file (str): the path of the .cif file.
    
    Output:
    V (np.array): the matrix of coordinates
    '''
    df = pd.DataFrame(columns=['x', 'y', 'z', 'residue','atom','res_idxs'])
    xs, ys, zs, residues, atoms, res_idxs = list(), list(), list(), list(), list(), list()
    with open(file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            columns = line.split()
            if line.startswith("ATOM"):
                xs.append(eval(columns[10]))
                ys.append(eval(columns[11]))
                zs.append(eval(columns[12]))
                residues.append(columns[5])
                atoms.append(columns[3])
                res_idxs.append(int(columns[8]))
    df['x'] = xs
    df['y'] = ys
    df['z'] = zs
    df['residue'] = residues
    df['atom'] = atoms
    df['res_idxs'] = res_idxs
    return df

def random_walk(n: int,uni_lim=1, R1=0, R2=np.inf) -> np.ndarray:
    r = (R1+R2)/2
    points = [[r,0,0]]
    count = 0
    while count<n-1:
        v = np.random.uniform(-uni_lim, uni_lim, 3)
        vec = [points[-1][0] + v[0], points[-1][1] + v[1], points[-1][2] + v[2]]
        if np.linalg.norm(vec)>=R1 and np.linalg.norm(vec)<=R2:
            points.append(vec)
            count+=1
    return np.array(points)

def self_avoiding_random_walk(n: int, step: float = 1.0, bead_radius: float = 0.5, epsilon: float = 0.001, two_dimensions=False) -> np.ndarray:
    potential_new_step = [0, 0, 0]
    while True:
        points = [np.array([0, 0, 0])]
        for _ in tqdm(range(n - 1)):
            step_is_ok = False
            trials = 0
            while not step_is_ok and trials < 1000:
                potential_new_step = points[-1] + step * random_versor()
                if two_dimensions:
                    potential_new_step[2] = 0
                for j in points:
                    d = dist(j, potential_new_step)
                    if d < 2 * bead_radius - epsilon:
                        trials += 1
                        break
                else:
                    step_is_ok = True
            points.append(potential_new_step)
        points = np.array(points)
        return points

def plot_structure(struct):
    # syntax for 3-D projection
    ax = plt.axes(projection ='3d')
    
    # defining all 3 axis
    z = struct[:,2]
    x = struct[:,0]
    y = struct[:,1]
    
    # plotting
    ax.plot3D(x, y, z, 'green')
    ax.set_title('3D line plot')
    plt.show()

def load_from_file(path,num_pts):
    V = get_coordinates_cif(path)
    x_sample, y_sample, z_sample = V[:,0], V[:,1], V[:,2]
    tck, u = interpolate.splprep([x_sample,y_sample,z_sample], s=2)
    x_knots, y_knots, z_knots = interpolate.splev(tck[0], tck)
    u_fine = np.linspace(0,1,num_pts)
    x_fine, y_fine, z_fine = interpolate.splev(u_fine, tck)
    points = 20*np.vstack((x_fine, y_fine, z_fine)).T
    return points
