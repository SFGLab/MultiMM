import pandas as pd
import numpy as np
from tqdm import tqdm

def NucHunter_loader(tsv_file,region,chrom,start_point=10):
    nuc_csv = pd.read_csv(tsv_file,sep='\t')
    nuc_pos = nuc_csv[nuc_csv['chromosome']==chrom]['position'].values
    set1 = np.where(nuc_pos>=region[0])[0]
    set2 = np.where(nuc_pos<=region[1])[0]
    idxs_in_region = np.intersect1d(set1,set2)
    nuc_pos_region = nuc_pos[idxs_in_region]
    min_pos = np.min(nuc_pos)
    average_dna_length = np.average(nuc_csv['fuzziness'])
    return nuc_pos_region//2-np.min(nuc_pos_region//2)+10

def drop_bedpe_columns(bedpe_file,out_name):
    bedpe_df = pd.read_csv(bedpe_file,sep='\t',header=None)

    content = ''
    print('Writing content...')
    for i in tqdm(range(len(bedpe_df))):
        content += f'{bedpe_df[0][i]}\t{bedpe_df[1][i]}\t{bedpe_df[2][i]}\n'
    print('Content Done!!\n\nWritting the file...')

    with open(out_name, 'w') as f:
        f.write(content)
    print('File done! :D')

def danpos_to_bedpe(file,output,chrom,region):
    nuc_csv = pd.read_csv(file,sep='\t')
    nuc_csv = nuc_csv[nuc_csv['chr']==chrom]
    nuc_csv = nuc_csv[nuc_csv['start']>=region[0]]
    nuc_csv = nuc_csv[nuc_csv['end']<=region[1]]
    starts = nuc_csv['start'].values
    ends = nuc_csv['end'].values

    content = ''
    for i in tqdm(range(len(nuc_csv))):
        start, end = starts[i], ends[i]
        content += f'{chrom}\t{start}\t{end}\n'

    with open(output, 'w') as f:
        f.write(content)

def bedpe_to_array(bedpe_file):
    nuc_csv = pd.read_csv(bedpe_file,sep='\t',header=None)
    start_pos = nuc_csv[1].values
    end_pos = nuc_csv[2].values
    print('Average nucleosome size',np.average(end_pos-start_pos))
    print('Average distance between nucleosomes:',np.average(start_pos[1:]-end_pos[:-1]))
    print('Minimum distance between nucleosomes:',np.min(start_pos[1:]-end_pos[:-1]))
    nuc_pos = 40*start_pos//145
    nuc_pos = nuc_pos-np.min(nuc_pos)+1
    return nuc_pos

def puffin_to_array(nucs_path):
    nuc_df = pd.read_csv(nucs_path,sep='\t')
    nuc_df = nuc_df.sort_values(by=['# Location']).reset_index(drop=True)
    starts = nuc_df['# Location'].values
    ends = nuc_df['# Location'].values+146
    starts = 40*starts//146
    ends = 40*ends//146
    return starts.astype(int), ends.astype(int)

def read_loops(file,chrom,region):
    loop_csv = pd.read_csv(file,sep='\t',header=None)
    loop_csv = loop_csv[loop_csv[0]==chrom]
    loop_csv = loop_csv[loop_csv[1]>=region[0]]
    loop_csv = loop_csv[loop_csv[5]<=region[1]]
    return IndentationError