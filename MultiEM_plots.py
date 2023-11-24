import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import seaborn as sns
import pandas as pd
from MultiEM_utils import *
from ellipsoid_measure import *
from scipy.interpolate import splrep, BSpline

def plot_chrom_metrics_from_txt(path_names):
    # Explore the details of the file so as to initialize the data
    f = open(path_names[0]+f'/chromosomes_info/{chrs[0]}_metrics_info.txt', "r")
    metric_names = list()
    lines = f.readlines()
    for line in lines:
        name, value = line.split(': ')
        metric_names.append(name)
    num_metrics = len(metric_names)
    num_chroms = len(chrs)-1
    num_phases = len(path_names)
    f.close()
    data = np.zeros((num_metrics,num_chroms,num_phases))

    # Build the data
    for phase_count, p in enumerate(path_names):
        for i in range(num_chroms):
            f = open(p+f'/chromosomes_info/{chrs[i]}_metrics_info.txt', "r")
            lines = f.readlines()
            for metric_count, line in enumerate(lines):
                _, value = line.split(': ')
                value = eval(value[:-2])
                data[metric_count,i,phase_count] = value
            f.close()

    # Now plot the data
    for i_metric in range(num_metrics):
        xaxis_list, colors, metric_values_per_chrom = [] ,[], []
        for i_chrom in range(num_chroms):
            metric_values_per_chrom+=list(data[i_metric,i_chrom,:])
            for i_phase in range(num_phases):
                xaxis_list+=[2*(i_phase+0.1)+1.7*i_chrom/num_chroms]
                colors+=[chrom_colors[i_chrom]]
        figure(figsize=(10, 7))
        plt.scatter(xaxis_list,metric_values_per_chrom,c=colors,marker='o')
        plt.tick_params(
                axis='x',          # changes apply to the x-axis
                which='both',      # both major and minor ticks are affected
                bottom=False,      # ticks along the bottom edge are off
                top=False,         # ticks along the top edge are off
                labelbottom=False) # labels along the bottom edge are off
        for i_phase in range(num_phases+1):
            plt.axvline(x = 2*i_phase, color = 'b')
        plt.title(metric_names[i_metric],fontsize=18)
        plt.xlabel('Chromosomes',fontsize=16)
        plt.show()

def add_ell_volume(folder_name,N_G1=25, N_S=41, N_G2M=35, n_chroms=23):
    # G1 phase
    for i in range(N_G1):
        for j in range(n_chroms):
            f = open(folder_name+f'/metacell{i}_k562_G1/chromosomes_info/{chrs[j]}_metrics_info.txt', "a")
            V = get_coordinates_cif(folder_name+f'/metacell{i}_k562_G1/chromosomes/MultiEM_minimized_{chrs[j]}.cif')
            vol = calc_ellipsoid_ratio2(V)
            f.write(f"\nEllipsoid volume: {vol}.")
            f.close()

    # S phase
    for i in range(N_S):
        for j in range(n_chroms):
            f = open(folder_name+f'/metacell{i}_k562_S/chromosomes_info/{chrs[j]}_metrics_info.txt', "a")
            V = get_coordinates_cif(folder_name+f'/metacell{i}_k562_S/chromosomes/MultiEM_minimized_{chrs[j]}.cif')
            vol = calc_ellipsoid_ratio2(V)
            f.write(f"\nEllipsoid volume: {vol}.")
            f.close()

    # G2M phase
    for i in range(N_G2M):
        for j in range(n_chroms):
            f = open(folder_name+f'/metacell{i}_k562_G2M/chromosomes_info/{chrs[j]}_metrics_info.txt', "a")
            V = get_coordinates_cif(folder_name+f'/metacell{i}_k562_G2M/chromosomes/MultiEM_minimized_{chrs[j]}.cif')
            vol = calc_ellipsoid_ratio2(V)
            f.write(f"\nEllipsoid volume: {vol}.")
            f.close()

def plot_metacell_metrics_from_txt(folder_name,N_G1=25, N_S=41, N_G2M=35, n_chroms=23):
    # Explore the details of the file so as to initialize the data
    G1_Rgs, G1_sph_ratios = np.zeros((N_G1,n_chroms)), np.zeros((N_G1,n_chroms)) # num of metacells, num of chromosomes
    S_Rgs, S_sph_ratios = np.zeros((N_S,n_chroms)), np.zeros((N_S,n_chroms))
    G2M_Rgs, G2M_sph_ratios = np.zeros((N_G2M,n_chroms)), np.zeros((N_G2M,n_chroms))

    # G1 phase
    for i in range(N_G1):
        for j in range(n_chroms):
            f = open(folder_name+f'/metacell{i}_k562_G1/chromosomes_info/{chrs[j]}_metrics_info.txt', "r")
            lines = f.readlines()
            for line in lines:
                try:
                    name, value = line.split(': ')
                    if name=='Gyration radius': G1_Rgs[i,j] = eval(value[:-2])
                    if name=='Ellipsoid volume': G1_sph_ratios[i,j] = eval(value[:-2])
                except:
                    pass
            f.close()

    # S phase
    for i in range(N_S):
        for j in range(n_chroms):
            f = open(folder_name+f'/metacell{i}_k562_S/chromosomes_info/{chrs[j]}_metrics_info.txt', "r")
            lines = f.readlines()
            for line in lines:
                try:
                    name, value = line.split(': ')
                    if name=='Gyration radius': S_Rgs[i,j] = eval(value[:-2])
                    if name=='Ellipsoid volume': S_sph_ratios[i,j] = eval(value[:-2])
                except:
                    pass
            f.close()

    # G2M phase
    for i in range(N_G2M):
        for j in range(n_chroms):
            f = open(folder_name+f'/metacell{i}_k562_G2M/chromosomes_info/{chrs[j]}_metrics_info.txt', "r")
            lines = f.readlines()
            for line in lines:
                try:
                    name, value = line.split(': ')
                    if name=='Gyration radius': G2M_Rgs[i,j] = eval(value[:-2])
                    if name=='Ellipsoid volume': G2M_sph_ratios[i,j] = eval(value[:-2])
                except:
                    pass
            f.close()

    # plot
    ## Gyration radius
    fig, axes = plt.subplots(5, 5)
    fig.set_figheight(15)
    fig.set_figwidth(25)
    fig.suptitle('Gyration Radius',fontsize=35)
    fig.tight_layout()
    for i in range(n_chroms):
        Rgs = np.hstack((G1_Rgs[:,i],S_Rgs[:,i],G2M_Rgs[:,i]))
        x = np.arange(len(Rgs))
        regression_model = np.poly1d(np.polyfit(x, Rgs, 5))
        reg_inv = np.linspace(0,len(Rgs),10000)
        axes[i//5,i%5].plot(Rgs,'o')
        axes[i//5,i%5].plot(reg_inv,regression_model(reg_inv),'r-',lw=2)
        axes[i//5,i%5].set_xlabel(f'Chrom {chrs[i]}')
        axes[i//5,i%5].axvline(x = N_G1+0.5, color = 'g')
        axes[i//5,i%5].axvline(x = N_G1+N_S+0.5, color = 'g')
    plt.show()

    ## Ellipsoid ratio
    fig, axes = plt.subplots(5, 5)
    fig.set_figheight(15)
    fig.set_figwidth(25)
    fig.suptitle('Ellipsoid Volume',fontsize=35)
    fig.tight_layout()
    for i in range(n_chroms):
        ERs = np.hstack((G1_sph_ratios[:,i],S_sph_ratios[:,i],G2M_sph_ratios[:,i]))
        mask = ERs<np.mean(ERs)+0.5*np.std(ERs)
        ERs = ERs[mask]
        x = np.arange(len(ERs))
        regression_model = np.poly1d(np.polyfit(x, ERs, 5))
        reg_inv = np.linspace(0,len(ERs),10000)
        axes[i//5,i%5].plot(ERs,'o')
        axes[i//5,i%5].plot(reg_inv,regression_model(reg_inv),'r-',lw=2)
        axes[i//5,i%5].set_xlabel(f'Chrom {chrs[i]}')
        axes[i//5,i%5].axvline(x = N_G1+0.5, color = 'g')
        axes[i//5,i%5].axvline(x = N_G1+N_S+0.5, color = 'g')
    plt.show()

# # Example
# path_names = ['G1_txt',
#               'S_txt',
#               'G2M_txt']
# plot_chrom_metrics_from_txt(path_names) 

plot_metacell_metrics_from_txt('metacells')