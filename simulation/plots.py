import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import pandas as pd
from .utils import *
from sklearn.decomposition import PCA
import pyvista as pv

pv.set_jupyter_backend("server")
color_dict = {-2:'#bf0020', -1:'#e36a24', 1:'#20c8e6',2:'#181385',0:'#ffffff'}
comp_dict = {-2:'B2', -1:'B1', 1:'A2',2:'A1',0:'no compartment'}

def plot_projection(struct_3D,Cs,save_path):
    # Dimensionality Reduction
    pca = PCA(n_components=2)
    struct_2d = pca.fit_transform(struct_3D)
    colors, comps = [], []
    for c in Cs[:len(struct_3D)]: 
        colors.append(color_dict[c])
        comps.append(comp_dict[c])
    
    # Calculate Distances
    dists = list()
    for vec in struct_3D: dists.append(np.linalg.norm(vec))
    dists = np.array(dists)

    # Make dataframe
    df = pd.DataFrame()
    df['x'], df['y'], df['z'] = struct_3D[:,0], struct_3D[:,1], struct_3D[:,2]
    df['x_PCA'],df['y_PCA'] = struct_2d[:,0], struct_2d[:,1]
    df['distance'] = dists
    df['subcomp'] = Cs
    df['subcomp_text'] = comps
    df.drop(df[df['subcomp']==0.0].index,inplace=True)

    # Plot Distribution
    figure(figsize=(8, 8), dpi=100)
    sns.scatterplot(data=df,x='x_PCA',y='y_PCA',hue='subcomp', palette="coolwarm", alpha=0.5)
    plt.xlabel('First Pricipal Component')
    plt.ylabel('Second Pricipal Component')
    plt.title('Scatter Plot of PCA 2D Projection')
    plt.savefig(save_path+'plots/PCA.svg',format='svg',dpi=100)
    plt.close()

    # Plot more stuff
    figure(figsize=(8, 8), dpi=100)
    sns.kdeplot(data=df, x='x', y='y', palette="coolwarm", hue='subcomp', bw_adjust=.5)
    plt.title('Subcompartment 2D Density Plot')
    plt.savefig(save_path+'plots/density_subcomp.svg',format='svg',dpi=100)
    plt.close()
    
    figure(figsize=(8, 5), dpi=100)
    sns.kdeplot(data=df, x='distance', hue='subcomp', fill=True, palette='coolwarm')
    plt.title('Subcompartment Density Plot')
    plt.savefig(save_path+'plots/kde_subcomp.svg',format='svg',dpi=100)
    plt.close()

    figure(figsize=(8, 5), dpi=100)
    sns.kdeplot(data=df,x='distance',fill=True)
    plt.title('Density Plot')
    plt.savefig(save_path+'plots/kde.svg',format='svg',dpi=100)
    plt.close()
    
    figure(figsize=(10, 8), dpi=100)
    sns.kdeplot(data=df, x='x_PCA', y='y_PCA', cmap="gnuplot2", shade=True,cbar=True)
    plt.title('2D Density Plot')
    plt.savefig(save_path+'plots/density.svg',format='svg',dpi=100)
    plt.close()

def polyline_from_points(points):
    poly = pv.PolyData()
    poly.points = points
    the_cell = np.arange(0, len(points), dtype=np.int_)
    the_cell = np.insert(the_cell, 0, len(points))
    poly.lines = the_cell
    return poly

def viz_structure(V, colors=None, r=0.1, cmap='coolwarm'):
    polyline = polyline_from_points(V)
    polyline["scalars"] = np.arange(polyline.n_points)

    if colors is not None:
        colors = colors[:len(V)]
        color_values = (colors - np.min(colors)) / (np.max(colors) - np.min(colors))  # Normalize colors
        polyline["colors"] = color_values  # Set colors as point scalars
        polymer = polyline.tube(radius=r)
        polymer.plot(smooth_shading=True, cmap=cmap, scalars="colors", show_scalar_bar=False)
    else:
        polymer = polyline.tube(radius=r)
        polymer.plot(smooth_shading=True, show_scalar_bar=False)

def viz_chroms(sim_path,r=0.1,comps=True):
    cif_path = sim_path + '/MultiMM_minimized.cif'
    chrom_idxs_path = sim_path + '/chrom_idxs.npy'
    chrom_comps_path = sim_path + '/compartments.npy'
    chrom_ends_path = sim_path + '/chrom_lengths.npy'
    chrom_idxs = np.load(chrom_idxs_path)
    if comps: comps = np.load(chrom_comps_path)
    chrom_ends = np.load(chrom_ends_path)
    V = get_coordinates_cif(cif_path)
    N = len(V)
    chroms = np.zeros(N)
    for i in range(len(chrom_ends)-1):
        start, end = chrom_ends[i], chrom_ends[i+1]
        chroms[start:end] = chrom_idxs[i]
    viz_structure(V,chroms[:len(V)],cmap='gist_ncar',r=r)
    if comps: viz_structure(V,comps[:len(V)],cmap='coolwarm',r=r)

def get_heatmap(cif_file,viz=False,th=1,save=False,save_path=None,vmax=1,vmin=0):
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
    V = get_coordinates_cif(cif_file)
    print('Matrix shape:',V.shape)
    mat = distance.cdist(V, V, 'euclidean') # this is the way \--/
    mat = 1/(mat+1)

    if viz:
        figure(figsize=(15, 12),dpi=500)
        plt.imshow(mat,cmap="Reds",vmax=vmax,vmin=vmin)
        if save: plt.savefig(save_path,format='svg',dpi=500)
        plt.colorbar()
        plt.show()
        if save: np.save(save_path.replace("svg", "npy"),mat)
    return mat
