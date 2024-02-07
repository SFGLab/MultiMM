from structuregenerator.generator import self_avoiding_random_walk, line, polymer_circle
import numpy as np
from points_io.hybrid36 import hy36encode
from nucs_helixes_banecki import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from scipy import interpolate
from tqdm import tqdm

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

def interpolate_structure(n_dna,path):
    print('Interpolating...')
    V = 40*get_coordinates_cif(path)
    x_sim, y_sim, z_sim = V[:,0], V[:,1], V[:,2]
    tck, u = interpolate.splprep(x=[x_sim,y_sim,z_sim], s=2)
    # x_knots, y_knots, z_knots = interpolate.splev(tck[0], tck)
    u_fine = np.linspace(0,1,n_dna)
    x_fine, y_fine, z_fine = interpolate.splev(u_fine, tck)
    V_interpol = np.vstack((x_fine,y_fine,z_fine)).T
    return V_interpol

def plot_3d_struct(V):
    x, y, z = V[:,0], V[:,1], V[:,2]
    ax = plt.figure().add_subplot(projection='3d')
    ax.plot(x, y, z, c='red', alpha=0.9)
    ax.legend()
    plt.show()

def build_init_mmcif_nucs(n_nucs,n_dna,entry,exit,mode='circle',psf=False,path=None):
    # Define the initial coordinates of histones and the structure of DNA
    hists = {0:'HA',1:'HB',2:'HC',3:'HD'}
    if mode=='circle':
        dna_points = polymer_circle(n_dna,10)
    elif mode=='path':
        dna_points = interpolate_structure(n_dna,path)
    else:
        raise Exception("You can choose between 'circle' or 'path' for mode value.")
    dna_points, nucs_pos = add_helixes(dna_points,entry,exit)
    p_x, p_y, p_z = nucs_pos[:,0], nucs_pos[:,1], nucs_pos[:,2]
    
    # Write the positions in .mmcif file
    atoms = ''
    hist_list = []
    n_hist = 4
    ## Histone beads
    print('Creating histones...')
    for k in tqdm(range(n_nucs)):
        # Some parameters
        histones = np.array(polymer_circle(n_hist))
        histone_points = move_structure_to(histones,dna_points[entry[k]],dna_points[exit[k]],nucs_pos[k])

        # Add 4 histone atoms
        for i in range(n_hist):
            x_coord, y_coord, z_coord = histone_points[i]
            hist_list.append((x_coord,y_coord,z_coord))
            atoms += ('{0:} {1:} {2:} {3:} {4:} {5:} {6:} {7:} {8:} '
                      '{9:} {10:.3f} {11:.3f} {12:.3f}\n'.format('HETATM',4*k+i+1,'H',hists[i],\
                                                            '.','HIS','A',1,k+1,'?',\
                                                            x_coord,y_coord,z_coord))

    print('Creating DNA...')
    ## DNA beads
    for i in tqdm(range(n_dna)):
        if i==0 or i==(n_dna-1):
            res_name, atom_name = 'ALB', 'CB'
        else:
            res_name, atom_name = 'ALA', 'CA'
        atoms += ('{0:} {1:} {2:} {3:} {4:} {5:} {6:} {7:} {8:} '
                '{9:} {10:.3f} {11:.3f} {12:.3f}\n'.format('ATOM', n_nucs*n_hist+i+1, 'C', atom_name,\
                                                            '.', res_name, 'A', 1, n_nucs+1+i, '?',\
                                                            dna_points[i][0], dna_points[i][1], dna_points[i][2]))
    print('Writting connections...')
    # Write connections
    connects = ''

    ## For histone atoms
    count = 0
    
    for k in tqdm(range(n_nucs)):
        for i in range(0,4):
            for j in range(i+1,4):
                count+=1
                connects += f'H{count} covale HXX A {k+1} {hists[i]} HXX A {k+1} {hists[j]}\n'

    connects += f'D1 covale ALB A {n_nucs+1} CB ALA A {n_nucs+2} CA\n'
    count = 1
    for k in tqdm(range(1,n_dna-2)):
        count+=1
        connects += f'D{count} covale ALA A {n_nucs+k+1} CA ALA A {n_nucs+k+2} CA\n'
    count+=1
    connects += f'D{count} covale ALA A {n_dna+n_nucs-1} CA ALB A {n_dna+n_nucs} CB\n'

    # Save files
    ## .pdb
    mmcif_file_name = 'dna_histones.cif'
    atomhead = mmcif_atomhead()
    conhead = mmcif_connecthead()
    mmcif_file_content = atomhead+atoms+conhead+connects

    with open(mmcif_file_name, 'w') as f:
        f.write(mmcif_file_content)

    print("File {} saved...".format(mmcif_file_name))

def psf_maker_nucs(n_nucs,n_dna,file_name,title="No title provided"):
    assert len(title) < 40, "provided title in psf file is too long."

    # parameters
    hists = {0:'HIA',1:'HIB',2:'HIC',3:'HID'}
    n_his = 4
    
    # noinspection PyListCreation
    lines = ['PSF CMAP\n']
    lines.append('\n')
    lines.append('      1 !NTITLE\n')
    lines.append('REMARKS {}\n'.format(title))
    lines.append('\n')
    lines.append('{:>8} !NATOM\n'.format(n_his*n_nucs+n_dna))

    # atoms
    ## histone atoms
    # Write the positions in .mmcif file
    atoms = ''
    n_hist = 4
    ## Histone beads
    for k in range(n_nucs):
        # Some parameters
        histone_points = polymer_circle(n_hist)
        
        # Add 4 histone atoms
        for i in range(n_hist):
            lines.append('{:>8} HIS  {:<5}HXX  H     {:<3}   0.000000        1.00 0           0\n'.format(hy36encode(4,4*k+i+1), hy36encode(4,k+1), hists[i]))

    ## DNA atoms
    lines.append('{:>8} DNA  {:<5}ALB  {:<2}   C      0.000000        1.00 0           0\n'.format(hy36encode(4,n_nucs*n_his+1), hy36encode(4,n_nucs+1), 'CB'))
    for k in range(1,n_dna-1):
        lines.append('{:>8} DNA  {:<5}ALA  {:<2}   C      0.000000        1.00 0           0\n'.format(hy36encode(4,n_nucs*n_his+k+1), hy36encode(4,n_nucs+k+1), 'CA'))
    lines.append('{:>8} DNA  {:<5}ALB  {:<2}   C      0.000000        1.00 0           0\n'.format(hy36encode(4,n_nucs*n_his+n_dna), hy36encode(4,n_nucs+n_dna), 'CB'))
    lines.append('\n')
    lines.append('{:>8} !NBOND: bonds\n'.format(n_dna + 6*n_nucs - 2))
    
    # Write connections
    ## For histone atoms
    for k in range(n_nucs):
        for i in range(0,4):
            for j in range(i+1,4):
                lines.append('{:>8}{:>8}\n'.format(hy36encode(4,4*k+i+1), hy36encode(4,4*k+j+1)))
    ## for dna atoms            
    for k in range(1,n_dna-1):
        lines.append('{:>8}{:>8}\n'.format(hy36encode(4,n_nucs*n_his+k+1), hy36encode(4,n_nucs*n_his+k+2)))
    
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