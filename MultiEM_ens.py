from MultiEM import *

chrom_lengths = [0,248387328,242696752,201105948,193574945,
                 182045439,172126628,160567428,146259331,
                 150617247,134758134,135127769,133324548,
                 113566686,101161492,99753195,96330374,
                 84276897,80542538,61707364,66210255,
                 45090682,51324926,154259566,62460029]

N_beads, n_chrom = 50000, 24
resolution = np.sum(chrom_lengths)//N_beads
chrom_ends = np.cumsum(chrom_lengths//resolution)
chrom_ends[-1] = N_beads

for i in range(35):
	path = f'metacell{i}_k562_G1/'
	try:
	    os.mkdir(path)
	    os.mkdir(path+'chromosomes')
	except OSError as error:
	    print("Folder 'chromosomes' already exists!")
	write_chrom_colors(chrom_ends,name=path+'MultiEM_chromosome_colors.cmd')

	# Load from files
	ms, ns = import_mns_from_txt(txt_file=f'/mnt/raid/data/single_cell/k562_metacells/k562.PET.txt.G1_{i}.rmblacklist.txt',
	                             N_beads=N_beads,n_chroms=n_chrom,path=path)

	# Run simulation
	md = MultiEM(Cs=None,chrom_ends=chrom_ends,ms=ms,ns=ns,ks=None,ds=None,path=path)
	md.run_pipeline(run_MD=False,build_init_struct=True,Temperature=310,
	                init_struct_path=None,plots=False)
