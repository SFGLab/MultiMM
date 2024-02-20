import seaborn as sns 
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import pandas as pd
from MultiEM_utils import *
from sklearn.decomposition import PCA

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

    # Plot Distribution
    figure(figsize=(8, 8), dpi=200)
    sns.jointplot(data=df,x='x_PCA',y='y_PCA',hue='subcomp', palette="seismic", kind='hist',alpha=0.5)
    plt.xlabel('First Pricipal Component')
    plt.ylabel('Second Pricipal Component')
    plt.savefig(save_path+'PCA.pdf',format='pdf',dpi=200)
    plt.show()

    # Plot more stuff
    figure(figsize=(8, 8), dpi=200)
    sns.kdeplot(data=df, x='x', y='y', cmap="Reds", shade=True)
    plt.savefig(save_path+'density.pdf',format='pdf',dpi=200)
    plt.show()

    figure(figsize=(8, 8), dpi=200)
    sns.kdeplot(data=df, x='x', y='y', palette="seismic", hue='subcomp', bw_adjust=.5)
    plt.savefig(save_path+'density_subcomp.pdf',format='pdf',dpi=200)
    plt.show()
    
    figure(figsize=(8, 5), dpi=200)
    sns.kdeplot(data=df, x='distance', hue='subcomp', fill=True, palette='seismic')
    plt.savefig(save_path+'kde_subcomp.pdf',format='pdf',dpi=200)
    plt.show()

    figure(figsize=(8, 5), dpi=200)
    sns.kdeplot(data=df,x='distance',fill=True)
    plt.savefig(save_path+'kde.pdf',format='pdf',dpi=200)
    plt.show()