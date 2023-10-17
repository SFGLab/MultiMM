from structuregenerator.generator import self_avoiding_random_walk, line, polymer_circle, baseball
import numpy as np
import pandas as pd
from points_io.hybrid36 import hy36encode
from mpl_toolkits import mplot3d
from scipy import interpolate
import matplotlib.pyplot as plt
# from hilbert import decode
import plotly.express as px
from hilbertcurve.hilbertcurve import HilbertCurve

def write_cmm(comps,name='MultiEM_compartment_colors.cmd'):
    comp_old = 2
    counter, start = 0, 0
    comp_dict = {-1:'blue', 1:'red'}
    content = ''

    for i, comp in enumerate(comps):
        if comp_old==comp:
            counter+=1
        elif i!=0:
            content+=f'color {comp_dict[comp_old]} :{start}-{start+counter+1}\n'
            counter, start = 0, i
        comp_old=comp

    content+=f'color {comp_dict[comp]} :{start}-{start+counter+1}\n'
    with open(name, 'w') as f:
        f.write(content)

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

def generate_hilbert_curve(n_points,p=10,n=3,viz=False):
    p=7; n=3
    hilbert_curve = HilbertCurve(p, n)
    distances = list(range(n_points))
    points = hilbert_curve.points_from_distances(distances)
    if viz:
        fig = plt.figure()
        ax = plt.axes(projection='3d')

        z = np.array(points)[:,2]
        x = np.array(points)[:,0]
        y = np.array(points)[:,1]

        ax.plot3D (x, y, z, 'green')
        ax.set_title('Hilbert Curve')
        plt.show()
    return np.array(points)

def build_init_mmcif(n_dna,psf=True,path=''):
    # Define the initial coordinates of histones and the structure of DNA
    dna_points = generate_hilbert_curve(n_dna)
    
    # Write the positions in .mmcif file
    atoms = ''
    
    ## DNA beads
    for i in range(n_dna):
        [atom_type,res_name, atom_name, cl] = ['ATOM','ALA', 'CA', 'A'] if (i!=0 and i!=(n_dna-1)) else ['HETATM','ALB', 'CB', 'A']
        atoms += ('{0:} {1:} {2:} {3:} {4:} {5:} {6:} {7:} {8:} '
                  '{9:} {10:.3f} {11:.3f} {12:.3f}\n'.format(atom_type, i+1, 'D', atom_name,\
                                                             '.', res_name, cl, 1, i+1, '?',\
                                                             dna_points[i,0], dna_points[i,1], dna_points[i,2]))
    
    # Write connections
    connects = ''
    
    for i in range(n_dna-1):
        [atom_type1,res_name1, atom_name1, cl1] = ['ATOM','ALA', 'CA', 'A'] if (i!=0 and i!=(n_dna-1)) else ['HETATM','ALB', 'CB', 'A']
        [atom_type2,res_name2, atom_name2, cl2] = ['ATOM','ALA', 'CA', 'A'] if ((i+1)!=0 and (i+1)!=(n_dna-1)) else ['HETATM','ALB', 'CB', 'A']
        connects += f'D{i+1} covale {res_name1} {cl1} {i+1} {atom_name1} {res_name2} {cl2} {i+2} {atom_name2}\n'

    # Save files
    ## .pdb
    mmcif_file_name = path+'MultiEM_init.cif'
    atomhead = mmcif_atomhead()
    conhead = mmcif_connecthead()
    mmcif_file_content = atomhead+atoms+'\n'+conhead+connects

    with open(mmcif_file_name, 'w') as f:
        f.write(mmcif_file_content)

    if psf:
        generate_psf(n_dna,path+'MultiEM.psf')

    print("File {} saved...".format(mmcif_file_name))

def write_mmcif(coords,path):
    # Write the positions in .mmcif file
    atoms = ''
    
    ## DNA beads
    for i in range(len(coords)):
        [res_name, atom_name, cl] = ['ALA', 'CA', 'A'] if (i!=0 and i!=len(coords)-1) else ['ALB', 'CB', 'A']
        atoms += ('{0:} {1:} {2:} {3:} {4:} {5:} {6:} {7:} {8:} '
                  '{9:} {10:.3f} {11:.3f} {12:.3f}\n'.format('ATOM', i+1, 'D', atom_name,\
                                                             '.', res_name, cl, 1, i+1, '?',\
                                                             coords[i,0], coords[i,1], coords[i,2]))
    
    # Write connections
    connects = ''
    
    for i in range(len(coords)-1):
        [res_name1, atom_name1, cl1] = ['ALA', 'CA', 'A'] if (i!=0 and i!=(len(coords)-1)) else ['ALB', 'CB', 'A']
        [res_name2, atom_name2, cl2] = ['ALA', 'CA', 'A'] if ((i+1)!=0 and (i+1)!=(len(coords)-1)) else ['ALB', 'CB', 'A']
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
    head = """data_nucsim
# 
_entry.id nucsim
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

def plot_structure(path,with_nucleosomes=False):
    df = get_df_cif(file=path)
    fig = px.line_3d(df[df['residue']=='ALA'], x="x", y="y", z="z")
    fig.update_layout(
        autosize=True,
        width=1000,
        height=900,
        margin=dict(r=0, l=0, b=0, t=0))
    fig.show()

# def hilbert_curve(num_pts,num_bits=2,scale=50):
#     # The maximum Hilbert integer.
#     num_dims=3
#     max_h = 2**(num_bits*num_dims)
    
#     # Generate a sequence of Hilbert integers.
#     hilberts = np.arange(max_h)

#     # Compute the 2-dimensional locations.
#     locs = scale*decode(hilberts, num_dims, num_bits)
    
#     # Spline interpolation
#     x_sample, y_sample, z_sample = locs[:,0], locs[:,1], locs[:,2]
#     tck, u = interpolate.splprep([x_sample,y_sample,z_sample], s=2)
#     x_knots, y_knots, z_knots = interpolate.splev(tck[0], tck)
#     u_fine = np.linspace(0,1,num_pts)
#     x_fine, y_fine, z_fine = interpolate.splev(u_fine, tck)
#     points = np.vstack((x_fine, y_fine, z_fine)).T
#     return points

def load_from_file(path,num_pts):
    V = get_coordinates_cif(path)
    x_sample, y_sample, z_sample = V[:,0], V[:,1], V[:,2]
    tck, u = interpolate.splprep([x_sample,y_sample,z_sample], s=2)
    x_knots, y_knots, z_knots = interpolate.splev(tck[0], tck)
    u_fine = np.linspace(0,1,num_pts)
    x_fine, y_fine, z_fine = interpolate.splev(u_fine, tck)
    points = 20*np.vstack((x_fine, y_fine, z_fine)).T
    return points