import straw # Install it with:  python3 -m pip install hic-straw
import numpy as np
import os
import tqdm
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm
from scipy.sparse import coo_matrix
from matplotlib.pyplot import figure
from numpy import linalg as LA
import plotly.express as px

def extract_data_along_hic_diagonal(primary,chrom,resolution,region,normalization="NONE",name=None,N_beads=None,colormap='bright_red'):
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
        data.append(counts[i])
        if binX[i] != binY[i]:
            row_indices.append(binY[i])
            col_indices.append(binX[i])
            data.append(counts[i])
    row_indices = np.asarray(row_indices) / resolution
    col_indices = np.asarray(col_indices) / resolution
    max_size = int(max(np.max(row_indices), np.max(col_indices))) + 1
    matrix = coo_matrix((data, (row_indices.astype(int), col_indices.astype(int))),
                            shape=(max_size, max_size)).toarray()
    matrix[np.isnan(matrix)] = 0
    matrix[np.isinf(matrix)] = 0
    matrix_full = matrix.T + matrix  # make matrix symmetric
    np.fill_diagonal(matrix_full, np.diag(matrix))  # prevent doubling of diagonal from prior step

    matrix_region = matrix_full[region[0]//resolution:region[1]//resolution,region[0]//resolution:region[1]//resolution]
    
    print('full matrix shape:',matrix_full.shape)

    return matrix_full, matrix_region

def vix_matrix(M,err=1):
    figure(figsize=(14, 13))
    plt.imshow(M,cmap="coolwarm",vmax=np.average(M)+err*np.std(M),vmin=np.average(M)-err*np.std(M))
    plt.xlabel('Genomic Distance',fontsize=16)
    plt.ylabel('Genomic Distance',fontsize=16)
    plt.colorbar()
    plt.show()

res = 100000
chrom = '1'
reg = [221989754,222662002]

M, m = extract_data_along_hic_diagonal(primary='/mnt/raid/data/encode/HiC/ENCSR968KAY_HiC/ENCFF256UOW.hic',\
                                        chrom=chrom,\
                                        resolution=res,\
                                        region=reg,\
                                        name = 'HiC',\
                                        normalization="NONE")
vix_matrix(M)
print('Length of matrix of region:',len(m))

D = 100/(M+1)**1/3

N = len(D)
d2 = np.sum(D**2,axis=0)/N-np.sum(D**2)/(2*N**2)

G = np.zeros((N,N))
for i in tqdm(range(N)):
    for j in range(N):
        G[i,j]=0.5*(d2[i]+d2[j]-D[i,j]**2)

eigenvalues, eigenvectors = LA.eig(G)
eigenvalues, eigenvectors = np.sort(eigenvalues)[::-1], eigenvectors[np.argsort(eigenvalues)[::-1]]
S = np.zeros((N,3))
for i in range(3):
    S[:,i] = np.sqrt(eigenvalues[i])*eigenvectors[i]/np.sum(eigenvectors[i])

df = pd.DataFrame()
x,y,z = S[:,0], S[:,1], S[:,2]
df['x'], df['y'], df['z'] = x,y,z
figure(figsize=(25, 25))
fig = px.line_3d(df, x='x', y='y', z='z')
fig.show()