#########################################################################
########### CREATOR: SEBASTIAN KORSAK, WARSAW 2022 ######################
#########################################################################

import hicstraw as straw # Install it with:  python3 -m pip install hic-straw
from matplotlib.pyplot import figure
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
import networkx as nx

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

chrom_colors = ['#ab0215','#f50a0a','#f5540a','#f5b20a','#e9f50a','#b6f50a','#3df50a',
              '#0af59b','#0af5f1','#0abef5','#0a70f5','#0a2df5','#770af5',
              '#a60af5','#ce0af5','#f50ae1','#e37dcd','#d494c6','#d9abcf',
              '#d9bdd3','#e6dce4','#f7f7f7','#1c4e4f','#4e4f1c']

def make_folder(N_beads,chrom,region):
    folder_name = f'stochastic_model_Nbeads_{N_beads}_chr_{chrom}_region_{region[0]}_{region[1]}'
    try:
        os.mkdir(folder_name)
        os.mkdir(folder_name+'/plots')
        os.mkdir(folder_name+'/other')
        os.mkdir(folder_name+'/pdbs')
        os.mkdir(folder_name+'/heatmaps')
    except OSError as error:
        print(f'Directory with name "{folder_name}" already exists! No problem lets continue!')
    return folder_name

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

def get_hic(primary,chrom,resolution,th1=10,th2=40,normalization="NONE",name=None,viz=False,N_beads=None,colormap='bright_red',binary=False):
    try:
        os.mkdir('heatmaps')
        os.mkdir('heatmaps/heat_plots')
        os.mkdir('heatmaps/arrays')
    except OSError as error:
        print(f'Directories already exist.')

    # assumes KR normalization and BP resolutions
    result = straw.straw(normalization, primary, chrom, chrom, "BP", resolution)
    binX = result[0]
    binY = result[1]
    counts = result[2]
    N = len(binX)
    row_indices, col_indices, data = list(), list(), list()
    for i in tqdm(range(N)):
        row_indices.append(binX[i])
        col_indices.append(binY[i])
        if binary:
            data.append(1) if counts[i]>th1 and counts[i]<th2 else data.append(0)
        else:
            data.append(counts[i])
        if binX[i] != binY[i]:
            row_indices.append(binY[i])
            col_indices.append(binX[i])
            if binary:
                data.append(1) if counts[i]>th1 and counts[i]<th2 else data.append(0)
            else:
                data.append(counts[i])
    row_indices = np.asarray(row_indices) / resolution
    col_indices = np.asarray(col_indices) / resolution
    max_size = int(max(np.max(row_indices), np.max(col_indices))) + 1
    matrix = coo_matrix((data, (row_indices.astype(int), col_indices.astype(int))),
                         shape=(max_size, max_size)).toarray()
    
    matrix = matrix.T + matrix  # make matrix symmetric
    matrix = (matrix-matrix.min())/(matrix.max()-matrix.min())
    print('full matrix shape:',matrix.shape)
    matrix[np.isnan(matrix)] = 0
    matrix[np.isinf(matrix)] = 0
    np.fill_diagonal(matrix, 0)

    if viz:
        counts = np.array(counts)
        plt.hist(counts[counts<50],bins=20)
        plt.grid()
        plt.show()
        
        figure(figsize=(15, 15))
        plt.imshow(matrix,cmap='gnuplot2_r')
        plt.title('Vizualized Experimental Heatmap',fontsize=20)
        plt.colorbar()
        plt.show()

    return matrix

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

    np.save(path+'genomewide_signal.npy',genomewide_signal[:N_beads])
    
    return genomewide_signal[:N_beads], polymer_lengths

def import_compartments_from_bed(bed_file,N_beads,n_chroms,path):
    comps_df = pd.read_csv(bed_file,header=None,sep='\t')

    # Find maximum coordinate of each chromosome
    print('Cleaning and transforming subcompartments dataframe...')
    chrom_ends = np.cumsum(chrom_lengths_array)
    # s, chrom_ends = 0, [0]
    # for i in tqdm(range(n_chroms)):
    #     max_chr = np.max(comps_df[2][comps_df[0]==chrs[i]])
    #     s+=max_chr
    #     chrom_ends.append(s)
    # chrom_ends=np.array(chrom_ends)
    
    # Sum bigger chromosomes with the maximum values of previous chromosomes
    for i in tqdm(range(n_chroms)):
        comps_df[1][comps_df[0]==chrs[i]]=comps_df[1][comps_df[0]==chrs[i]]+chrom_ends[i]
        comps_df[2][comps_df[0]==chrs[i]]=comps_df[2][comps_df[0]==chrs[i]]+chrom_ends[i]

    # Convert genomic coordinates to simulation beads
    resolution = int(np.max(comps_df[2].values))//N_beads
    chrom_ends = np.array(chrom_ends)//resolution
    chrom_ends[-1] = N_beads
    np.save(path+'chrom_lengths.npy',chrom_ends)
    comps_df[1], comps_df[2] = comps_df[1]//resolution, comps_df[2]//resolution

    # Convert compartemnts to vector
    print('Building subcompartments_array...')
    comps_array = np.zeros(N_beads)
    for i in tqdm(range(len(comps_df))):
        if comps_df[3][i].startswith('A.1'):
            val = 2
        elif comps_df[3][i].startswith('A.2'):
            val = 1
        elif comps_df[3][i].startswith('B.1'):
            val = -1
        elif comps_df[3][i].startswith('B.2'):
            val = -2
        comps_array[comps_df[1][i]:comps_df[2][i]] = val
    np.save(path+'genomewide_signal.npy',comps_array)
    print('Done')
    return comps_array, chrom_ends

def write_chrom_colors(chrom_ends,name='MultiEM_chromosome_colors.cmd'):    
    content = ''
    for i in range(len(chrom_ends)-1):
        content+=f'color {chrom_colors[i]} :{chrom_ends[i]}-{chrom_ends[i+1]}\n'

    with open(name, 'w') as f:
        f.write(content)

def min_max_trans(x):
    return (x-x.min())/(x.max()-x.min())

def import_mns_from_bedpe(bedpe_file,N_beads,n_chroms,threshold=3,viz=False,mode='kd',min_loop_dist=3,path=''):
    ds, ks, cs  = None, None, None
    # Import loops
    chroms = list(chrs[i] for i in range(n_chroms))
    loops = pd.read_csv(bedpe_file,header=None,sep='\t')
    loops = loops[loops[0].isin(chroms) & loops[3].isin(chroms)].reset_index(drop=True)
    
    print('Cleaning and transforming loops dataframe...')
    # Find maximum coordinate of each chromosome
    s, chrom_ends = 0, [0]
    for i in tqdm(range(n_chroms)):
        max1 = np.max(loops[2][loops[0]==chrs[i]])
        max2 = np.max(loops[5][loops[3]==chrs[i]])
        s+=np.max([max1,max2])
        chrom_ends.append(int(s))
    
    # Sum bigger chromosomes with the maximum values of previous chromosomes
    for i in tqdm(range(n_chroms)):
        loops[1][loops[0]==chrs[i]]=loops[1][loops[0]==chrs[i]]+chrom_ends[i]
        loops[2][loops[0]==chrs[i]]=loops[2][loops[0]==chrs[i]]+chrom_ends[i]
        loops[4][loops[3]==chrs[i]]=loops[4][loops[3]==chrs[i]]+chrom_ends[i]
        loops[5][loops[3]==chrs[i]]=loops[5][loops[3]==chrs[i]]+chrom_ends[i]
    
    # Convert genomic coordinates to simulation beads
    resolution = int(np.max(loops[5].values))//N_beads
    chrom_ends = np.array(chrom_ends)//resolution
    chrom_ends[-1] = N_beads
    np.save(path+'chrom_lengths.npy',chrom_ends)
    loops[1], loops[2], loops[4], loops[5] = loops[1]//resolution, loops[2]//resolution, loops[4]//resolution, loops[5]//resolution
    loops['ms'] = (loops[1].values+loops[2].values)//2
    loops['ns'] = (loops[4].values+loops[5].values)//2
    loops['Total Count'] = loops.groupby(['ms', 'ns'])[6].transform('sum')
    counts = loops['Total Count'].values
    
    # Filter the ones above the threshold
    print('Importing loops...')
    mns, cs = np.vstack((loops['ms'].values[counts>threshold], loops['ns'].values[counts>threshold])), counts[counts>threshold]
    mns, idxs = np.unique(mns,axis=1,return_index=True)
    cs = cs[idxs]
    ms, ns = mns[0,:], mns[1,:]
    ms[ms>=N_beads],ns[ns>=N_beads]=N_beads-1, N_beads-1
    ms,ns,cs = ms[ns>ms+min_loop_dist], ns[ns>ms+min_loop_dist], cs[ns>ms+min_loop_dist]
    if mode=='k':
        zs = np.abs(np.log10(np.abs((cs-np.mean(cs)))/np.std(cs)))
        zs[zs>np.mean(zs)+np.std(zs)] = np.mean(zs)+np.std(zs)
        ds, ks = None, 50+2950*min_max_trans(zs)
    elif mode=='d':
        ks=None
        ds = 1/cs**(1/3)
        ds, ks = 0.1+0.3*min_max_trans(ds), None
    elif mode=='kd':
        zs = np.abs(np.log10(np.abs((cs-np.mean(cs)))/np.std(cs)))
        zs[zs>np.mean(zs)+np.std(zs)] = np.mean(zs)+np.std(zs)
        ds = 1/cs**(1/3)
        ds, ks = 0.1+0.3*min_max_trans(ds), 50+2950*min_max_trans(zs)
    else:
        raise InterruptedError("The mode of loop generator should be either 'k' so as to generate Hook constats, or 'd' for equillibrium distances, or 'kd' for both.")
    
    if viz:
        print('min k:',np.min(ks))
        print('max k:',np.max(ks))
        plt.hist(cs,bins=30)
        plt.xlabel('counts')
        plt.show()
        plt.hist(ks,bins=30)
        plt.xlabel('strengths')
        plt.show()

    # Perform some data cleaning
    mask = (ns-ms)!=0
    ms = ms[mask]
    ns = ns[mask]
    if mode=='kd' or mode=='k': ks= ks[mask]
    if mode=='kd' or mode=='d':ds = ds[mask]
    cs = cs[mask]
    N_loops = len(ms)
    np.save(path+'ms.npy',ms)
    np.save(path+'ns.npy',ns)
    np.save(path+'ks.npy',ks)
    np.save(path+'ds.npy',ds)
    print('Done! Number of loops is ',N_loops)
    return ms, ns, ds, ks, cs, chrom_ends

def import_mns_from_txt(txt_file,N_beads,n_chroms,threshold=1,min_loop_dist=5,path='',mode='kd',viz=False,graph_mode=True):
    '''
    This function is importer of (signle-cell).
    '''
    ds, ks, cs  = None, None, None
    # Import loops
    chroms = list(chrs[i] for i in range(n_chroms))
    loops = pd.read_csv(txt_file,header=None,sep='\t')
    loops = loops[loops[0].isin(chroms) & loops[3].isin(chroms)].reset_index(drop=True)
    
    print('Cleaning and transforming loops dataframe...')
    # Find maximum coordinate of each chromosome
    s, chrom_ends = 0, [0]
    for i in tqdm(range(n_chroms)):
        max1 = np.max(loops[2][loops[0]==chrs[i]])
        max2 = np.max(loops[5][loops[3]==chrs[i]])
        s+=np.max([max1,max2])
        chrom_ends.append(s)
    
    # Sum bigger chromosomes with the maximum values of previous chromosomes
    for i in tqdm(range(n_chroms)):
        loops[1][loops[0]==chrs[i]]=loops[1][loops[0]==chrs[i]]+chrom_ends[i]
        loops[2][loops[0]==chrs[i]]=loops[2][loops[0]==chrs[i]]+chrom_ends[i]
        loops[4][loops[3]==chrs[i]]=loops[4][loops[3]==chrs[i]]+chrom_ends[i]
        loops[5][loops[3]==chrs[i]]=loops[5][loops[3]==chrs[i]]+chrom_ends[i]
    
    # Convert genomic coordinates to simulation beads
    resolution = np.max(loops[5].values)//N_beads
    chrom_ends = np.array(chrom_ends)//resolution
    chrom_ends[-1] = N_beads
    np.save(path+'chrom_lengths.npy',chrom_ends)
    loops[1], loops[2], loops[4], loops[5] = loops[1]//resolution, loops[2]//resolution, loops[4]//resolution, loops[5]//resolution
    loops['ms'] = (loops[1].values+loops[2].values)//2
    loops['ns'] = (loops[4].values+loops[5].values)//2
    loops['count'] = 1
    loops['Total Count'] = loops.groupby(['ms', 'ns'])['count'].transform('sum')
    counts = loops['Total Count'].values
    
    # Filter the ones above the threshold
    print('Importing loops...')
    mns, cs = np.vstack((loops['ms'].values[counts>threshold], loops['ns'].values[counts>threshold])), counts[counts>threshold]
    mns, idxs = np.unique(mns,axis=1,return_index=True)
    cs = cs[idxs]
    ms, ns = mns[0,:], mns[1,:]
    ms[ms>=N_beads],ns[ns>=N_beads]=N_beads-1, N_beads-1
    ms,ns,cs = ms[ns>ms+min_loop_dist], ns[ns>ms+min_loop_dist], cs[ns>ms+min_loop_dist]
    if mode=='k':
        zs = np.abs(np.log10(np.abs((cs-np.mean(cs)))/np.std(cs)))
        zs[zs>np.mean(zs)+np.std(zs)] = np.mean(zs)+np.std(zs)
        ds, ks = None, 50+2950*min_max_trans(zs)
    elif mode=='d':
        ks=None
        ds = 1/cs**(1/3)
        ds, ks = 0.1+0.3*min_max_trans(ds), None
    elif mode=='kd':
        zs = np.abs(np.log10(np.abs((cs-np.mean(cs)))/np.std(cs)))
        zs[zs>np.mean(zs)+np.std(zs)] = np.mean(zs)+np.std(zs)
        ds = 1/cs**(1/3)
        ds, ks = 0.1+0.3*min_max_trans(ds), 50+2950*min_max_trans(zs)
    else:
        raise InterruptedError("The mode of loop generator should be either 'k' so as to generate Hook constats, or 'd' for equillibrium distances, or 'kd' for both.")
    
    if viz:
        print('min k:',np.min(ks))
        print('max k:',np.max(ks))
        plt.hist(cs,bins=30)
        plt.xlabel('counts')
        plt.show()
        plt.hist(ks,bins=30)
        plt.xlabel('strengths')
        plt.show()

    # Perform some data cleaning
    mask = (ns-ms)!=0
    ms = ms[mask]
    ns = ns[mask]
    if mode=='kd' or mode=='k': ks= ks[mask]
    if mode=='kd' or mode=='d':ds = ds[mask]
    cs = cs[mask]
    N_loops = len(ms)
    np.save(path+'ms.npy',ms)
    np.save(path+'ns.npy',ns)
    np.save(path+'ks.npy',ks)
    np.save(path+'ds.npy',ds)
    print('Done! Number of loops is ',N_loops)

    if graph_mode:
        # Build the graph
        H = nx.Graph()
        H.add_nodes_from(np.arange(N_beads)) # Add nodes

        # Avoid the ones in diagonal
        mask = (ns-ms)!=0
        ms, ns = ms[mask], ns[mask]

        # Connect adjacent nodes
        adj_w = np.full(N_beads,0.2)
        x, y = np.arange(N_beads-1), np.arange(1,N_beads)

        # Add loops with weight proportional to 1/c
        H.add_weighted_edges_from(zip(x,y,adj_w))
        H.add_weighted_edges_from(zip(ms,ns,0.1/cs+0.1))
        
        # Define the matrix with shortest paths
        L = np.zeros((N_beads,N_beads))
        ds = list()
        for m, n in zip(ms,ns):
            d = nx.shortest_path_length(H,source=m,target=n,weight='weight')
            L[m,n], L[n,m] = d, d
            ds.append(d)
        ds = np.array(ds)
    return ms, ns, ds, ks, cs, chrom_ends

