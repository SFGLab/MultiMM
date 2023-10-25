from MultiEM import *
from MultiEM_utils

N_beads = 50000
resolution = np.sum(chrom_lengths_array)//N_beads
chrom_ends = np.cumsum(chrom_lengths_array//resolution)
chrom_ends[-1] = N_beads

for i in range(35):
	# Run simulation
	path = f'metacell{i}_k562_G1/'
	lp = f'/mnt/raid/data/single_cell/k562_metacells/k562.PET.txt.G1_{i}.rmblacklist.txt'
	md = MultiEM(N_beads=N_beads,loop_path=lp,path=path)
	md.run_pipeline(run_MD=False,build_init_struct=True,Temperature=310,
	                init_struct_path=None,plots=False)
