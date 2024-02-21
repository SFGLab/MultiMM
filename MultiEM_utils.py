#########################################################################
########### CREATOR: SEBASTIAN KORSAK, WARSAW 2024 ######################
#########################################################################

from matplotlib.pyplot import figure
from matplotlib.colors import to_hex
import matplotlib.pyplot as plt
import os
import pandas as pd
import numpy as np
from scipy.sparse import coo_matrix
from scipy.spatial import distance
from tqdm import tqdm
import pyBigWig
from scipy import stats
import random as rd
from itertools import groupby
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

pd.options.mode.chained_assignment = None

chrs = {0:'chr1',1:'chr2',2:'chr3',3:'chr4',4:'chr5',5:'chr6',6:'chr7',7:'chr8',
          8:'chr9',9:'chr10',10:'chr11',11:'chr12',12:'chr13',13:'chr14',14:'chr15',
          15:'chr16',16:'chr17',17:'chr18',18:'chr19',19:'chr20',20:'chr21',21:'chr22',
          22:'chrX',23:'chrY'}

chrom_lengths_array = [0,248387328,242696752,201105948,193574945,
                 182045439,172126628,160567428,146259331,
                 150617247,134758134,135127769,133324548,
                 113566686,101161492,99753195,96330374,
                 84276897,80542538,61707364,66210255,
                 45090682,51324926,154259566,62460029]

chrom_sizes = {'chr1':248387328,'chr2':242696752,'chr3':201105948,'chr4':193574945,
               'chr5':182045439,'chr6':172126628,'chr7':160567428,'chr8':146259331,
               'chr9':150617247,'chr10':134758134,'chr11':135127769,'chr12':133324548,
               'chr13':113566686,'chr14':101161492,'chr15':99753195,'chr16':96330374,
               'chr17':84276897,'chr18':80542538,'chr19':61707364,'chr20':66210255,
               'chr21':45090682,'chr22':51324926,'chrX':154259566,'chrY':62460029}

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

def import_bw(bw_path,N_beads,n_chroms,viz=False,binary=False,path=''):
    genomewide_signal = list()
    bw = pyBigWig.open(bw_path)

    # Compute the total length of chromosomes
    chrom_length = 0
    lengths = list()
    for i in range(n_chroms):
        chrom_length += bw.chroms(chrs[i])
        lengths.append(bw.chroms(chrs[i]))
    step = chrom_length//N_beads-5*n_chroms

    # Just import signal
    mean=0
    print('Computing mean of bw signal...')
    for i in tqdm(range(n_chroms)):
        mean += bw.stats(chrs[i])[0]
    mean = mean/n_chroms

    # Import the downgraded signal
    polymer_lengths = [0]
    count=0
    print('Importing bw signal...')
    for i in tqdm(range(n_chroms)):
        signal = bw.values(chrs[i],0,-1)
        new_signal =  list()
        for i in range(step,len(signal)+1,step):
            new_signal.append(np.average(signal[(i-step):i]))
            count+=1
        polymer_lengths.append(count)
        genomewide_signal.append(new_signal)
    genomewide_signal=np.concatenate(genomewide_signal)
    print('Length of signal:',len(genomewide_signal))

    # Write chromosomes file
    np.save(path+'chrom_lengths.npy',polymer_lengths)
    
    # Transform signal to binary or adjuct it to have zero mean
    if binary:
        genomewide_signal[genomewide_signal<mean]=-1
        genomewide_signal[genomewide_signal>mean]=1
        
        # Subtitute zeros with random spin states
        mask = genomewide_signal==0
        n_zeros = np.count_nonzero(mask)
        nums = np.array(rd.choices([-1,1],k=n_zeros))
        genomewide_signal[mask] = nums
    else:
        genomewide_signal = genomewide_signal-mean
        max_val = np.max([np.max(genomewide_signal),-np.min(genomewide_signal)])
        genomewide_signal = genomewide_signal/max_val
    print('Done!\n')

    if viz:
        figure(figsize=(25, 5), dpi=600)
        xax = np.arange(len(genomewide_signal))
        plt.fill_between(xax,genomewide_signal, where=(new_signal>0),alpha=0.50,color='purple')
        plt.fill_between(xax,genomewide_signal, where=(new_signal<0),alpha=0.50,color='orange')
        plt.xlabel('Genomic Distance (x250kb)',fontsize=16)
        plt.ylabel('BW signal renormalized',fontsize=16)
        plt.grid()
        plt.close()

    np.save(path+'comps.npy',genomewide_signal[:N_beads])
    
    return genomewide_signal[:N_beads], polymer_lengths

def import_compartments_from_Calder(bed_file,N_beads,coords=None,chrom=None,save_path=''):
    # Load compartment dataset
    comps_df = pd.read_csv(bed_file,header=None,sep='\t')

    # Find maximum coordinate of each chromosome
    print('Cleaning and transforming subcompartments dataframe...')
    chrom_ends = np.cumsum(chrom_lengths_array) if chrom==None else np.array([0,chrom_sizes[chrom]])
    if chrom!=None:
        comps_df = comps_df[(comps_df[0]==chrom)&(comps_df[1]>coords[0]) & (comps_df[2]<coords[1])].reset_index(drop=True)
    n_chroms = len(np.unique(comps_df[0].values))
    
    # Sum bigger chromosomes with the maximum values of previous chromosomes
    if chrom==None:
        for i in tqdm(range(n_chroms)):
            comps_df[1][comps_df[0]==chrs[i]]=comps_df[1][comps_df[0]==chrs[i]]+chrom_ends[i]
            comps_df[2][comps_df[0]==chrs[i]]=comps_df[2][comps_df[0]==chrs[i]]+chrom_ends[i]

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
        if comps_df[3][i].startswith('A.1'):
            val = -2
        elif comps_df[3][i].startswith('A.2'):
            val = -1
        elif comps_df[3][i].startswith('B.1'):
            val = 1
        elif comps_df[3][i].startswith('B.2'):
            val = 2
        comps_array[comps_df[1][i]:comps_df[2][i]] = val
    np.save(save_path+'compartments.npy',comps_array)
    print('Done')
    return comps_array.astype(int), chrom_ends.astype(int)

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

def write_chrom_colors(chrom_ends,name='MultiEM_chromosome_colors.cmd'):    
    colors = integers_to_hex_colors(0, len(chrom_ends)-1)
    
    content = ''
    for i in range(len(chrom_ends)-1):
        content+=f'color {colors[i]} :.{chr(64+1+i)}\n'

    with open(name, 'w') as f:
        f.write(content)

def min_max_trans(x):
    return (x-x.min())/(x.max()-x.min())

def import_mns_from_bedpe(bedpe_file,N_beads,coords=None,chrom=None,threshold=3,viz=False,min_loop_dist=0,path=''):
    # Import loops
    loops = pd.read_csv(bedpe_file,header=None,sep='\t')
    n_chroms = len(np.unique(loops[0].values))
    chroms = list(chrs[i] for i in range(n_chroms)) if chrom==None else [chrom]
    if chrom!=None:
        loops = loops[(loops[0]==chrom)&(loops[1]>coords[0])&(loops[2]<coords[1])&(loops[4]>coords[0])&(loops[5]<coords[1])].reset_index(drop=True)
    chrom_ends = np.cumsum(chrom_lengths_array) if chrom==None else np.array([0,chrom_sizes[chrom]])
    
    print('Cleaning and transforming loops dataframe...')
    # Sum bigger chromosomes with the maximum values of previous chromosomes
    if chrom==None:
        for i in tqdm(range(n_chroms)):
            loops[1][loops[0]==chrs[i]]=loops[1][loops[0]==chrs[i]]+chrom_ends[i]
            loops[2][loops[0]==chrs[i]]=loops[2][loops[0]==chrs[i]]+chrom_ends[i]
            loops[4][loops[3]==chrs[i]]=loops[4][loops[3]==chrs[i]]+chrom_ends[i]
            loops[5][loops[3]==chrs[i]]=loops[5][loops[3]==chrs[i]]+chrom_ends[i]
    
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
    # zs = np.abs(np.log10(np.abs((cs-np.mean(cs)))/np.std(cs)))
    # zs[zs>np.mean(zs)+np.std(zs)] = np.mean(zs)+np.std(zs)
    # ks = 1000+599000*min_max_trans(zs)
    ds = 0.05+0.15*min_max_trans(1/cs**2/3)

    # Perform some data cleaning
    mask = (ns-ms)!=0
    ms = ms[mask]
    ns = ns[mask]
    avg_ls = np.average(ns-ms)
    print('Average loop size:',avg_ls)
    ds= ds[mask]
    cs = cs[mask]
    N_loops = len(ms)
    np.save(path+'ms.npy',ms)
    np.save(path+'ns.npy',ns)
    np.save(path+'ds.npy',ds)
    print('Done! Number of loops is ',N_loops)
    return ms, ns, ds, chrom_ends

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