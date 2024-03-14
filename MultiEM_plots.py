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
    df.drop(df[df['subcomp']==0.0].index,inplace=True)

    # Plot Distribution
    figure(figsize=(8, 8), dpi=250)
    sns.scatterplot(data=df,x='x_PCA',y='y_PCA',hue='subcomp', palette="coolwarm_r", alpha=0.5)
    plt.xlabel('First Pricipal Component')
    plt.ylabel('Second Pricipal Component')
    plt.title('Scatter Plot of PCA 2D Projection')
    plt.savefig(save_path+'plots/PCA.pdf',format='pdf',dpi=200)
    plt.close()

    # Plot more stuff
    figure(figsize=(8, 8), dpi=250)
    sns.kdeplot(data=df, x='x', y='y', palette="coolwarm_r", hue='subcomp', bw_adjust=.5)
    plt.title('Subcompartment 2D Density Plot')
    plt.savefig(save_path+'plots/density_subcomp.pdf',format='pdf',dpi=200)
    plt.close()
    
    figure(figsize=(8, 5), dpi=200)
    sns.kdeplot(data=df, x='distance', hue='subcomp', fill=True, palette='coolwarm_r')
    plt.title('Subcompartment Density Plot')
    plt.savefig(save_path+'plots/kde_subcomp.pdf',format='pdf',dpi=200)
    plt.close()

    figure(figsize=(8, 5), dpi=200)
    sns.kdeplot(data=df,x='distance',fill=True)
    plt.title('Density Plot')
    plt.savefig(save_path+'plots/kde.pdf',format='pdf',dpi=200)
    plt.close()
    
    figure(figsize=(10, 8), dpi=200)
    sns.kdeplot(data=df, x='x_PCA', y='y_PCA', cmap="gnuplot2", shade=True,cbar=True)
    plt.title('2D Density Plot')
    plt.savefig(save_path+'plots/density.pdf',format='pdf',dpi=200)
    plt.close()
