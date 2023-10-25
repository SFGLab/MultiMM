import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import seaborn as sns
import pandas as pd
from MultiEM_utils import *

def plot_chrom_metrics_from_txt(path_names):
    # Explore the details of the file so as to initialize the data
    f = open(path_names[0]+f'/chromosomes_info/{chrs[0]}_metrics_info.txt', "r")
    metric_names = list()
    lines = f.readlines()
    for line in lines:
        name, value = line.split(': ')
        metric_names.append(name)
    num_metrics = len(metric_names)
    num_chroms = len(chrs)
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

# Example
path_names = ['k562_G1_combined_all_k_big_structure',
              'k562_S_combined_all_k_big_structure',
              'k562_G2M_combined_all_k_big_structure']
plot_chrom_metrics_from_txt(path_names) 