import os, glob, re
import configparser
import random
import numpy as np

chrom_sizes = {'chr1':248387328,'chr2':242696752,'chr3':201105948,'chr4':193574945,
               'chr5':182045439,'chr6':172126628,'chr7':160567428,'chr8':146259331,
               'chr9':150617247,'chr10':134758134,'chr11':135127769,'chr12':133324548,
               'chr13':113566686,'chr14':101161492,'chr15':99753195,'chr16':96330374,
               'chr17':84276897,'chr18':80542538,'chr19':61707364,'chr20':66210255,
               'chr21':45090682,'chr22':51324926,'chrX':154259566,'chrY':62460029}

def random_region_generator(num, size):
    chroms, regions = list(), list()
    for i in range(num):
        chrom = 'chr'+str(random.choice([1,2,3,4,5,6,7,8]))
        rd = int(np.random.uniform(0,248387328-size-2,1))
        region = [rd,rd+size]
        chroms.append(chrom)
        regions.append(region)
    np.save('ens_regions.npy',regions)
    np.save('ens_chroms.npy',chroms)
    return chroms, regions

def modify_config(N_beads,config_path,chrom,region,loop_path,compartment_path,out_path):
    ## Change config.ini file
    config = configparser.RawConfigParser()
    config.read(config_path)
    config.set('Main', 'INITIAL_STRUCTURE_TYPE', 'circle')
    config.set('Main', 'CHROM', chrom)
    config.set('Main', 'LOC_START', region[0])
    config.set('Main', 'LOC_END', region[1])
    config.set('Main', 'LOOPS_PATH',loop_path)
    config.set('Main', 'COMPARTMENT_PATH',compartment_path)
    config.set('Main', 'OUT_PATH',out_path)
    with open('config.ini', 'w') as configfile:
        config.write(configfile)

def run_MultiEM(N_beads,chroms,regions):
    '''
    A function that runs MultiEM.
    '''
    families = ['YRB_f','YRB_m','YRB_d','PUR_f','PUR_m','PUR_d','CHS_f','CHS_m','CHS_d']
    idxs = ['GM19239','GM19238','GM19240','HG00731','HG00732','HG00733','HG00512','HG00513','HG00514']
    N_fam = len(families)
    N = len(regions)

    for i in range(N):
        for j in range(N_fam):
            loop_path = f'/mnt/raid/data/Trios/loops_pet1+_CTCFmotif_final/{idxs[j]}_ctcfmotif_1mb_pet1.bedpe'
            compartment_path = f'/mnt/raid/data/Trios/calder_HiChIP_subcomp/{families[j]}.bed'
            out_path = f'{families[j]}_{chroms[i]}_{regions[i][0]}_{regions[i][1]}'
            modify_config(N_beads,'/mnt/raid/codes/mine/MultiEM-main/config.ini',chroms[i],regions[i],loop_path,compartment_path,out_path)
            cmd = 'python MultiEM.py -c config.ini'
            os.system(cmd)

def main():
    chroms, regions = random_region_generator(10000, 1000000)
    run_MultiEM(2000,chroms,regions)

if __name__=='__main__':
    main()

