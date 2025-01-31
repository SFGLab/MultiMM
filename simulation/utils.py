#########################################################################
########### CREATOR: SEBASTIAN KORSAK, WARSAW 2024 ######################
#########################################################################

from matplotlib.pyplot import figure
from matplotlib.colors import to_hex
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.spatial import distance
from tqdm import tqdm
import pyBigWig
import random as rd
from itertools import groupby
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

pd.options.mode.chained_assignment = None

chrs = {0:'chr1',1:'chr2',2:'chr3',3:'chr4',4:'chr5',5:'chr6',6:'chr7',7:'chr8',
          8:'chr9',9:'chr10',10:'chr11',11:'chr12',12:'chr13',13:'chr14',14:'chr15',
          15:'chr16',16:'chr17',17:'chr18',18:'chr19',19:'chr20',20:'chr21',21:'chr22',
          22:'chrX',23:'chrY'}

chrom_lengths_array = np.array([0,248387328,242696752,201105948,193574945,
                 182045439,172126628,160567428,146259331,
                 150617247,134758134,135127769,133324548,
                 113566686,101161492,99753195,96330374,
                 84276897,80542538,61707364,66210255,
                 45090682,51324926,154259566,62460029])

chrom_sizes = {'chr1':248387328,'chr2':242696752,'chr3':201105948,'chr4':193574945,
               'chr5':182045439,'chr6':172126628,'chr7':160567428,'chr8':146259331,
               'chr9':150617247,'chr10':134758134,'chr11':135127769,'chr12':133324548,
               'chr13':113566686,'chr14':101161492,'chr15':99753195,'chr16':96330374,
               'chr17':84276897,'chr18':80542538,'chr19':61707364,'chr20':66210255,
               'chr21':45090682,'chr22':51324926,'chrX':154259566,'chrY':62460029}

def min_max_normalize(matrix, Min=0, Max=1):
    # Calculate the minimum and maximum values of the matrix
    matrix = np.nan_to_num(matrix)
    min_val = np.min(matrix)
    max_val = np.max(matrix)

    # Normalize the matrix using the min-max formula
    normalized_matrix = Min + (Max - Min) * ((matrix - min_val) / (max_val - min_val))

    return normalized_matrix

chrom_strength = 1-min_max_normalize(chrom_lengths_array[1:])

def get_coordinates_mm(mm_vec):
    '''
    It returns the corrdinate matrix V (N,3) of a .pdb file.
    The main problem of this function is that coordiantes are not always in 
    the same column position of a .pdb file. Do changes appropriatelly,
    in case that the data aren't stored correctly. 
    
    Input:
    file (Openmm Qunatity): an OpenMM vector of the form
    Quantity(value=[Vec3(x=0.16963918507099152, y=0.9815883636474609, z=-1.4776774644851685), 
    Vec3(x=0.1548253297805786, y=0.9109517931938171, z=-1.4084612131118774), 
    Vec3(x=0.14006929099559784, y=0.8403329849243164, z=-1.3392155170440674), 
    Vec3(x=0.12535107135772705, y=0.7697405219078064, z=-1.269935131072998),
    ...,
    unit=nanometer)
    
    Output:
    V (np.array): the matrix of coordinates
    '''
    V = list()

    for i in range(len(mm_vec)):
        x, y ,z = mm_vec[i][0]._value, mm_vec[i][1]._value, mm_vec[i][2]._value
        V.append([x, y, z])
    
    return np.array(V)

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

def get_heatmap(mm_vec,viz=False,save=False,path=''):
    '''
    It returns the corrdinate matrix V (N,3) of a .pdb file.
    The main problem of this function is that coordiantes are not always in 
    the same column position of a .pdb file. Do changes appropriatelly,
    in case that the data aren't stored correctly.
    
    Input:
    file (Openmm Qunatity): an OpenMM vector of the form 
    Quantity(value=[Vec3(x=0.16963918507099152, y=0.9815883636474609, z=-1.4776774644851685),
    Vec3(x=0.1548253297805786, y=0.9109517931938171, z=-1.4084612131118774),
    Vec3(x=0.14006929099559784, y=0.8403329849243164, z=-1.3392155170440674),
    Vec3(x=0.12535107135772705, y=0.7697405219078064, z=-1.269935131072998),
    ...,
    unit=nanometer)
    
    Output:
    H (np.array): a heatmap of the 3D structure.
    '''
    V = get_coordinates_mm(mm_vec)
    mat = distance.cdist(V, V, 'euclidean') # this is the way \--/
    mat = 1/(mat+1)
    
    if save: np.save('sim_heat.npy',mat)
    
    if viz:
        figure(figsize=(25, 20))
        plt.imshow(mat,cmap="coolwarm",vmax=np.average(mat)+3*np.std(mat),vmin=np.average(mat)-3*np.std(mat))
        plt.savefig(path+'heatmap.svg',format='svg',dpi=500)
        plt.savefig(path+'heatmap.pdf',format='pdf',dpi=500)
        plt.close()
    return mat

def compute_averages(arr1, N2):
    # Calculate the window size
    window_size = len(arr1) // N2

    # Reshape the array into a 2D array with the specified window size
    reshaped_arr = arr1[:N2 * window_size].reshape(N2, -1)

    # Calculate the average along the specified axis (axis=1)
    averaged_arr = np.mean(reshaped_arr, axis=1)

    return averaged_arr

def import_bed(bed_file,N_beads,coords=None,chrom=None,save_path='',shuffle=False,seed=0,n_chroms = 22):
    # Load compartment dataset
    np.random.seed(seed)
    comps_df = pd.read_csv(bed_file,header=None,sep='\t')

    # Find maximum coordinate of each chromosome
    print('Cleaning and transforming subcompartments dataframe...')
    if chrom!=None:
        comps_df = comps_df[(comps_df[0]==chrom)&(comps_df[1]>coords[0]) & (comps_df[2]<coords[1])].reset_index(drop=True)
    chrom_idxs = np.arange(n_chroms).astype(int)
    if shuffle: np.random.shuffle(chrom_idxs)
    chrom_ends = np.cumsum(np.insert(chrom_lengths_array[1:][chrom_idxs], 0, 0)) if chrom==None else np.array([0,chrom_sizes[chrom]])
    
    # Sum bigger chromosomes with the maximum values of previous chromosomes
    if chrom==None:
        for count, i in enumerate(chrom_idxs):
            comps_df[1][comps_df[0]==chrs[i]]=comps_df[1][comps_df[0]==chrs[i]]+chrom_ends[count]
            comps_df[2][comps_df[0]==chrs[i]]=comps_df[2][comps_df[0]==chrs[i]]+chrom_ends[count]

    # Convert genomic coordinates to simulation beads
    resolution = chrom_ends[-1]//N_beads if chrom==None else (coords[1]-coords[0])//N_beads
    chrom_ends = np.array(chrom_ends)//resolution
    chrom_ends[-1] = N_beads
    np.save(save_path+'chrom_lengths.npy',chrom_ends)
    if chrom!=None:
        comps_df[1], comps_df[2] = comps_df[1]-coords[0], comps_df[2]-coords[0]
    comps_df[1], comps_df[2] = comps_df[1]//resolution, comps_df[2]//resolution
    
    # Convert compartemnts to vector
    print('Building subcompartments_array...')
    comps_array = np.zeros(N_beads)
    for i in tqdm(range(len(comps_df))):
        if comps_df[3][i].startswith('A.1') or comps_df[3][i].startswith('A1'):
            val = 2
        elif comps_df[3][i].startswith('A.2') or comps_df[3][i].startswith('A2') or comps_df[3][i].startswith('A'):
            val = 1
        elif comps_df[3][i].startswith('B.2') or comps_df[3][i].startswith('B2'):
            val = -2
        elif comps_df[3][i].startswith('B.1') or comps_df[3][i].startswith('B1') or comps_df[3][i].startswith('B'):
            val = -1

        comps_array[comps_df[1][i]:comps_df[2][i]] = val
    np.save(save_path+'compartments.npy',comps_array)
    np.save(save_path+'chrom_idxs.npy',chrom_idxs)
    print('Done')
    return comps_array.astype(int), chrom_ends.astype(int), chrom_idxs.astype(int)

def align_comps(comps,ms,chrom_ends):
    for i in range(len(chrom_ends)-1):
        start, end = chrom_ends[i], chrom_ends[i+1]
        mms = ms[(start<ms)&(ms<end)]
        comps_with_loops = comps[mms]
        Aloops = np.count_nonzero(comps_with_loops>0)
        Bloops = np.count_nonzero(comps_with_loops<0)
        if Aloops>Bloops: comps[start:end] = -comps[start:end]
    return comps

def integers_to_hex_colors(start, end):
    # Generate a range of integers
    integers = np.arange(start, end + 1)

    # Map each integer to a rainbow color and convert to hex format
    rgb_colors = plt.cm.rainbow(integers / max(integers))
    hex_colors = [to_hex(color) for color in rgb_colors]

    return hex_colors

def write_chrom_colors(chrom_ends,chrom_idxs,name='MultiMM_chromosome_colors.cmd'):    
    colors = integers_to_hex_colors(0, len(chrom_ends)+1)
    
    content = ''
    for i in range(len(chrom_ends)-1):
        content+=f'color {colors[chrom_idxs[i]]} :.{chr(64+1+i)}\n'

    with open(name, 'w') as f:
        f.write(content)

def min_max_trans(x):
    return (x-x.min())/(x.max()-x.min())

def process_file(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            columns = line.strip().split('\t')
            # Repeat the second and fourth columns, and add a new column with 1
            new_line = f"{columns[0]}\t{columns[1]}\t{columns[1]}\t{columns[2]}\t{columns[3]}\t{columns[3]}\t1\n"
            outfile.write(new_line)

def downsample_arrays(ms, ns, cs, ds, down_prob):
    assert len(ms) == len(ns) == len(cs) == len(ds), "Arrays must have the same length"
    
    # Generate a mask of indices to keep
    mask = np.random.rand(len(ms)) < down_prob
    indices = np.where(mask)[0]

    # Apply the same mask to all arrays
    ms_downsampled = ms[indices]
    ns_downsampled = ns[indices]
    cs_downsampled = cs[indices]
    ds_downsampled = ds[indices]

    return ms_downsampled, ns_downsampled, cs_downsampled, ds_downsampled

def import_mns_from_bedpe(bedpe_file, N_beads, coords=None, chrom=None, threshold=0, min_loop_dist=2, path='', down_prob=1.0, shuffle=False, seed=0, n_chroms=22):
    # Import loops
    np.random.seed(seed)
    loops = pd.read_csv(bedpe_file,header=None,sep='\t')
    chrom_idxs = np.arange(n_chroms).astype(int)
    if shuffle: np.random.shuffle(chrom_idxs)
    if chrom!=None:
        loops = loops[(loops[0]==chrom)&(loops[1]>coords[0])&(loops[2]<coords[1])&(loops[4]>coords[0])&(loops[5]<coords[1])].reset_index(drop=True)
    chrom_ends = np.cumsum(np.insert(chrom_lengths_array[1:][chrom_idxs], 0, 0)) if chrom==None else np.array([0,chrom_sizes[chrom]])
    print('Cleaning and transforming loops dataframe...')
    
    # Sum bigger chromosomes with the maximum values of previous chromosomes
    if chrom==None:
        for count, i in enumerate(chrom_idxs):
            loops[1][loops[0]==chrs[i]]=loops[1][loops[0]==chrs[i]]+chrom_ends[count]
            loops[2][loops[0]==chrs[i]]=loops[2][loops[0]==chrs[i]]+chrom_ends[count]
            loops[4][loops[3]==chrs[i]]=loops[4][loops[3]==chrs[i]]+chrom_ends[count]
            loops[5][loops[3]==chrs[i]]=loops[5][loops[3]==chrs[i]]+chrom_ends[count]
    
    # Convert genomic coordinates to simulation beads
    resolution = int(np.max(loops[5].values))//N_beads if chrom==None else (coords[1]-coords[0])//N_beads
    chrom_ends = np.array(chrom_ends)//resolution
    chrom_ends[-1] = N_beads
    np.save(path+'chrom_lengths.npy',chrom_ends)
    if chrom!=None:
        loops[1], loops[2], loops[4], loops[5] = loops[1]-coords[0], loops[2]-coords[0], loops[4]-coords[0], loops[5]-coords[0]
    loops[1], loops[2], loops[4], loops[5] = loops[1]//resolution, loops[2]//resolution, loops[4]//resolution, loops[5]//resolution
    loops['ms'] = (loops[1].values+loops[2].values)//2
    loops['ns'] = (loops[4].values+loops[5].values)//2
    loops['Total Count'] = loops.groupby(['ms', 'ns'])[6].transform('mean')
    counts = loops['Total Count'].values
    
    # Filter the ones above the threshold
    print('Importing loops...')
    mns, cs = np.vstack((loops['ms'].values[counts>threshold], loops['ns'].values[counts>threshold])), counts[counts>threshold]
    mns, idxs = np.unique(mns,axis=1,return_index=True)
    cs = cs[idxs]
    ms, ns = mns[0,:], mns[1,:]
    ms[ms>=N_beads],ns[ns>=N_beads]=N_beads-1, N_beads-1
    ms,ns,cs = ms[ns>ms+min_loop_dist], ns[ns>ms+min_loop_dist], cs[ns>ms+min_loop_dist]
    ds = 0.1+0.1*min_max_trans(1/cs**2/3) if not np.all(cs==cs[0]) else np.ones(len(ms)) 

    # Perform some data cleaning
    mask = (ns-ms)!=0
    ms = ms[mask]
    ns = ns[mask]
    ds= ds[mask]
    cs = cs[mask]

    if down_prob<1.0:
        ms, ns, cs, ds = downsample_arrays(ms, ns, cs, ds, down_prob)

    avg_ls = np.average(ns-ms)
    print('Average loop size:',avg_ls)

    N_loops = len(ms)
    np.save(path+'chrom_idxs.npy',chrom_idxs)
    np.save(path+'ms.npy',ms)
    np.save(path+'ns.npy',ns)
    np.save(path+'ds.npy',ds)
    print('Done! Number of loops is ',N_loops)
    return ms.astype(int), ns.astype(int), ds, chrom_ends.astype(int), chrom_idxs.astype(int)

def generate_arrays(N_loops, N, l=6):
    # Generate array ms with random integers between 0 and N (exclusive)
    ms = np.random.randint(0, N, size=N_loops)

    # Generate array ns by adding a random integer from an exponential distribution with average l
    ns = ms + np.round(np.random.exponential(l, size=N_loops)).astype(int)
    ns = np.maximum(ns, 3)
    ns = np.minimum(ns, N-1)

    # Define array ks with weights randomly distributed in the scale of 50 to 3000
    ks = np.random.uniform(50, 3000, N_loops)

    return ms, ns, ks

def shuffle_blocks(array):
    # Identify the unique blocks of repeating elements
    unique_blocks = [list(g) for k, g in groupby(array)]
    
    # Shuffle the unique blocks
    np.random.shuffle(unique_blocks)
    
    # Concatenate the shuffled blocks to create the final shuffled array
    shuffled_array = [elem for block in unique_blocks for elem in block]

    return shuffled_array

def import_bw(bw_path,N_beads,coords=None,chrom=None,viz=False,binary=False,path='',norm=False,shuffle=False,seed=0,n_chroms=22):
    '''
    Imports .BigWig data and outputs compartments.

    It assumes that higher signal coresponds to B compartment.

    In case that you would like to switch the sign then add flag sign=-1.
    '''
    # Open file
    np.random.seed(seed)
    bw = pyBigWig.open(bw_path)
    chrom_idxs = np.arange(n_chroms).astype(int)
    if shuffle: np.random.shuffle(chrom_idxs)
    print('Number of chromosomes:',n_chroms)

    # Compute the total length of chromosomes
    if chrom==None:
        chrom_length = 0
        lengths = list()
        for i in range(n_chroms):
            chrom_length += bw.chroms(chrs[chrom_idxs[i]])
            lengths.append(bw.chroms(chrs[chrom_idxs[i]]))
        lengths = np.array(lengths)
        resolution = chrom_length//(2*N_beads)
        polymer_lengths = lengths//resolution
        np.save(path+'chrom_lengths.npy',polymer_lengths)

    # Import the downgraded signal
    print('Importing bw signal...')
    if chrom==None:
        genomewide_signal = list()
        for i in tqdm(range(n_chroms)):
            signal = bw.values(chrs[chrom_idxs[i]],0,-1, numpy=True)
            signal = np.nan_to_num(signal, copy=True, nan=0.0, posinf=0.0, neginf=0.0)
            genomewide_signal.append(compute_averages(signal, polymer_lengths[i]))
        genomewide_signal = np.concatenate(genomewide_signal)
    else:
        genomewide_signal = bw.values(chrom,coords[0],coords[1], numpy=True)
        genomewide_signal = np.nan_to_num(genomewide_signal, copy=True, nan=0.0, posinf=0.0, neginf=0.0)
    bw.close()

    genomewide_signal = compute_averages(genomewide_signal,N_beads)
    if norm: genomewide_signal = (genomewide_signal-np.mean(genomewide_signal)+3*np.std(genomewide_signal))/np.std(genomewide_signal)
    
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
    
    return genomewide_signal
