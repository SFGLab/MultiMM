# MultiMM: An OpenMM-based software for whole genome 3D structure reconstruction

MultiMM is an OpenMM model designed for modeling the 3D structure of the whole *human* genome. Its distinguishing feature is that it is multiscale, meaning it aims to model different levels of chromatin organization, from smaller scales (nucleosomes) to the level of chromosomal territories. The algorithm is both fast and accurate. A key feature enabling its speed is GPU parallelization via OpenMM, along with smart assumptions that assist the optimizer in finding the global minimum. One such fundamental assumption is the use of a Hilbert curve as the initial structure. This helps MultiMM to converge faster because the initial structure is already highly compacted.

![GW_em](https://github.com/user-attachments/assets/3d019616-2c0f-4bfc-a792-d87fcdc0d96a)

After running the MultiMM model, users obtain a genome-wide structure. Users can color different chromosomes or compartments, or they can visualize individual chromosomes separately. MultiMM is simple to use, as it is based on modifying a single configuration file. Here we explain its usage.

![MultiMM_scales](https://github.com/user-attachments/assets/a0d14ddc-41bf-4d14-8a8c-aa418fe575b5)

## Key Features

- OpenMM based.
- User-friendly software available via PyPI.
- Efficient simulation of chromatin interactions using 3D conformations.
- Can simulate the scales of nucleosomes, TADs, compartments, chromosomal territories and interactions with lamina. Scalable simulations across different force fields and resolution levels.
- Compatible with modern GPU-accelerated simulation libraries. CPU acceleration can also be done.
- Possibility of creation of ensembles of 3D structures of chromatin.

## About Operating Systems

MultiMM has been tested primarily on Linux-based operating systems, with successful tests in Ubuntu, Debian, and Red Hat-based distributions. It is also possible to run it on macOS, though without CUDA support, which is helpful for accelerating computations. We do not recommend running MultiMM on Windows systems.

## Installation
MultiMM can be easily installed with pip:

```bash
pip install MultiMM
```

PyPI software: https://pypi.org/project/MultiMM/.

## Input Data

MultiMM relies on three types of datasets:

- Loop interactions in `bedpe` format (mandatory).
- Compartmentalization data in `bed` format (optional).
- ATAC-Seq p-value data in `.BigWig` format (optional).

### Loops (bedpe file)

For **loop interactions**, users must provide a file describing interactions between anchor 1 and anchor 2, along with their strength.

The file must be in .bedpe format, contain 7 columns, and must not include a header.

```text
chr10	100225000	100230000	chr10	100420000	100425000	95
chr10	100225000	100230000	chr10	101005000	101010000	56
chr10	101190000	101195000	chr10	101370000	101375000	152
chr10	101190000	101200000	chr10	101470000	101480000	181
chr10	101600000	101605000	chr10	101805000	101810000	152
```
The file can include interactions from all chromosomes; MultiMM will automatically handle them.

For single-cell data, users should prepare the file by setting the second and third columns identical (as well as the fifth and sixth columns) and set the strength value to 1.


### Compartments

For **(sub)compartment interactions**, the file should be in the format produced by the CALDER software: https://github.com/CSOgroup/CALDER2. Users do not need to run CALDER specifically, but the file format must match. The file should contain at least the first four columns with chromosome, regions, and the subcompartment label. Example:

```text
chr1	700001	900000	A.1.2.2.2.2.2.2	0.875	.	700001	900000	#FF4848
chr1	900001	1400000	A.1.1.1.1.2.1.1.1.1.1	1	.	900001	1400000	#FF0000
chr1	1400001	1850000	A.1.1.1.1.2.1.2.2.2.1	1	.	1400001	1850000	#FF0000
chr1	1850001	2100000	B.1.1.2.2.1.2.1	0.5	.	1850001	2100000	#DADAFF
```

### Nucleosomes

For ATAC-Seq data, users must provide a BigWig file containing p-values. The pyBigWig library is required to read BigWig files. Note: pyBigWig is not compatible with Windows systems.

### Definition of a region based on a gene (optional)

To model a genomic region around a specific gene, you can load a `.tsv` file containing gene information. Specify either the gene name or the gene ID.

The `.tsv` file should be formatted as follows:

```text
gene_id	gene_name	chromosome	start	end
ENSG00000160072	ATAD3B	chr1	1471765	1497848
ENSG00000279928	DDX11L17	chr1	182696	184174
ENSG00000228037		chr1	2581560	2584533
ENSG00000142611	PRDM16	chr1	3069168	3438621
ENSG00000284616		chr1	5301928	5307394
ENSG00000157911	PEX10	chr1	2403964	2413797
ENSG00000269896		chr1	2350414	2352820
ENSG00000228463		chr1	257864	359681
ENSG00000260972		chr1	5492978	5494674
ENSG00000224340		chr1	10054445	10054781
```

In case that you specify the gene, MultiMM will output a visualization with the gene as well. For example, in the folowing region we can see the polymer structure that is modelled (around the gene) and the gene with red color.

The MultiMM model targets the gene file automatically, so you do not have to provide it. However, you can optionally change it.

![minimized_structure_gene_coloring](https://github.com/user-attachments/assets/15666286-0162-4fdd-a875-6b50b65049fb)


> **Note:** Currently, MultiMM is designed to work with human genome data. While it may be possible to run the code on other organisms with additional debugging and modifications, full support for other species is planned for future versions. MultiMM can process various types of datasets. It is capable of calling loops from a range of experiments, including Hi-C, scHi-C, ChIA-PET, and Hi-ChIP. However, we cannot guarantee that the default parameters are optimal for every dataset. Therefore, users are encouraged to test the software carefully and verify the convergence of the algorithm with their own data.Before adjusting any parameters, please read the method paper thoroughly to understand the role and impact of each force.

## Usage
All the model's parameters are specified in a `config.ini` file. This file should have the following format:

```ini
[Main]

; Platform selection
PLATFORM = OpenCL

; Input data
FORCEFIELD_PATH = forcefields/ff.xml
LOOPS_PATH = /home/skorsak/Data/Rao/GSE63525_GM12878_primary+replicate_HiCCUPS_looplist_hg19.bedpe
COMPARTMENT_PATH = /home/skorsak/Data/Rao/subcompartments_primary_replicate/sub_compartments/all_sub_compartments.bed
ATACSEQ_PATH = /home/skorsak/Data/encode/ATAC-Seq/ENCSR637XSC_GM12878/ENCFF667MDI_pval.bigWig
OUT_PATH = application_note

; Simulation Parameters
N_BEADS = 50000
SHUFFLE_CHROMS = True
NUC_DO_INTERPOLATION = True

; Enable forcefield for genome-wide simulation
SC_USE_SPHERICAL_CONTAINER = True
CHB_USE_CHROMOSOMAL_BLOCKS = True
SCB_USE_SUBCOMPARTMENT_BLOCKS = True
IBL_USE_B_LAMINA_INTERACTION = True
CF_USE_CENTRAL_FORCE = True

; Simulation Parameters
SIM_RUN_MD = True
SIM_SAMPLING_STEP = 50
SIM_N_STEPS = 1000
TRJ_FRAMES = 100
```

After specifying the parameters and forces, users can run the following command in the terminal:

```bash
MultiMM -c config.ini
```

The software will output a folder with the resulting structure and plots showing compartment distribution.

Example data can be found on Google Drive: https://drive.google.com/drive/folders/1nFAPE4pCaHpeL5nw6nq0VvfUFoc24aXm?usp=sharing. Note that this data is publicly available from Rao et al. The subcompartment predictions were made using CALDER, and the ATAC-Seq data is from ENCODE.

In the `examples` folder, we provide example configuration files for different modeling scenarios.

## Helper argument based on modelling levels

In the newer version of MultiMM, we have included the `MODELLING_LEVEL` argument. This is a magic argument that helps users who are new to molecular modelling to specify the parameters of the model depending on the resolution they would like to model.

Therefore, we have the following modelling levels:

- `GENE`: The user should provide the gene of interest and a `.bedpe` file path, and MultiMM will model the gene (with the default gene_window). In this modelling level compartment forces are neglected.
- `REGION`: The user has to povide the chromosome and the coordinates of interest. They can provide compartment interactions as well (optionally). MultiMM models only the specified region in the genome.
- `CHROMOSOME`: If the user would like to model a whole chromosome. Then the user provides the chromosome number, and the model internally specifies the start and ending coordinates. The user can import compartments as well.
- `GW`: For genome-wide simulation. In this case, MultiMM models all chromosomes, and thus users do not have to provide any chromosome or coordinates. This is the most computationally expensive option and the computation can take from minutes to hours depending on the hardware.

This argument specifies the number of simulation beads appropriatelly as well. No matter what number of beads the user will give, if `MODELLING_LEVEL` is specified, `N_BEADS` will change to the default one. Namely, for gene it is 1000 beads, for a region it is 5000 beads, for chromosome 20000 beads and for the whole genome 200000 beads.

This is a convenient argument for new users. Nevertheness, we suggest more advanced users to not use this arguent, if they need more freedom in the definition of the parameters of their choice.

## Vizualization

For vizualization purposes, if you would like to import the whole genome structure you may use the command,

```python
import simulation.plots as splt

splt.viz_chroms(sim_path)
```

For `sim_path` you should add the output folder directory path (add `comps=Fase` in case that you do not need compartment coloring). Otherwise, in case that the user would like to model a particular region, without using compartment or chromosome coloring (pretty much any cif structure), they can type,

```python
import simulation.plots as splt
import simulation.utils as suts

V = suts.get_coordinates_cif(cif_path)
splt.viz_structure(V)
``` 

We would like once again to thank people who developed `pyvista` library and allow us to have fast and good vizualizations of large chromatin structure.

## Simulation Arguments

MultiMM has numerous configurable parameters. Below is a description of each argument and its default values. The defaults have been tested for genome-wide simulation but can be modified in the configuration file if needed. Units are typically assumed based on OpenMM conventions, though explicit unit specification is not required.

Here we can see the long table of the simulation arguments. Somtimes MultiMM might not work for some choices of arguments. For example:

* If you would like to model lamina interaction having disabled compartment interactions.
* If you do not provide appropriate data i.e. for compartmentalization but you have enabled (sub)compartent-specific forcefield.

Therefore, it is advisable to read the paper and understand well the meaning of each force before you start running simulations. MultiMM is a research model, not a market product and thus it requires a level of expertise (despite the easiness of usage) to underatand and run it.


| Argument Name                | Type         | Value       | Units         | Description |
|------------------------------|--------------|-------------|---------------|-------------|
| PLATFORM                     | str          | CPU        | CPU          | name of the platform. Available choices: CPU, OpenCL, CUDA. |
| CPU_THREADS                  | int          | None        | None          | Number of CPU threads (in case that CPU is chosen as platform). |
| DEVICE                       | str          | None        | None          | device index for CUDA or OpenCL (count from 0) |
| MODELLING_LEVEL               | str          | None        | None          | Helping function to specify parameters of simulation. Choose 'GENE', 'REGION', 'CHROM' or 'GW' depending on the resolution of interest. |
| INITIAL_STRUCTURE_PATH       | str          | None        | None          | Path to CIF file. |
| BUILD_INITIAL_STRUCTURE      | bool         | True        | None          | To build a new initial structure. |
| INITIAL_STRUCTURE_TYPE       | str          | hilbert     | None          | you can choose between: hilbert, circle, rw, confined_rw, self_avoiding_rw, helix, spiral, sphere, knot. |
| GENERATE_ENSEMBLE            | bool         | False       | False         | True if you need to create an ensemble of structures. |
| N_ENSEMBLE                   | int          | None        | None          | Number of samples that you would like to create. |
| DOWNSAMPLING_PROB            | float        | 1.0        | 1.0          | Probability of downsampling (from 0 to 1). There is no downsampling by default (p=1). |
| FORCEFIELD_PATH              | str          | None        | None          | Path to XML file with forcefield. |
| N_BEADS                      | int          | 50000       | None          | Number of Simulation Beads. |
| COMPARTMENT_PATH             | str          | None        | None          | It should be  a `.bed` file with subcompartments from Calder. |
| LOOPS_PATH                   | str          | None        | None          | A `.bedpe` file path with loops. It is required. |
| ATACSEQ_PATH                 | str          | None        | None          | A `.bw` or `.BigWig` file path with atacseq data for nucleosome interpolation. It is not required. |
| OUT_PATH                     | str          | results     | None          | Output folder name. |
| GENE_TSV      | str  | ''    | default_path   | A .tsv with genes and their locations in the genome. This targets in the interan directory by default. |
| GENE_NAME     | str  | ''    | None   | The name of the gene of interest. |
| GENE_ID       | str  | ''    | None   | The id of the gene of interest. |
| GENE_WINDOW   | int  | 100000| bp     | The window around the area of the gene of interest. |
| LOC_START                    | int          | None        | None          | Starting region coordinate (*in case that you do not want to model the whole genome*). |
| LOC_END                      | int          | None        | None          | Ending region coordinate (*in case that you do not want to model the whole genome*). |
| CHROM                        | str          | None        | None          | Chromosome that corresponds the the modelling region of interest (*in case that you do not want to model the whole genome*). |
| SHUFFLE_CHROMS               | bool         | False       | None          | Shuffle the chromosomes. |
| SHUFFLING_SEED               | int          | 0           | None          | Shuffling random seed. |
| SAVE_PLOTS                   | bool         | True        | None          | Save plots. |
| POL_USE_HARMONIC_BOND        | bool         | True        | None          | Use harmonic bond interaction for consecutive beads ($i,i\pm 1$). |
| POL_HARMONIC_BOND_R0         | float        | 0.1         | nm          | Harmonic bond distance equilibrium constant. |
| POL_HARMONIC_BOND_K          | float        | 300000.0    | kJ/mol/nm^2   | harmonic bond force constant |
| POL_USE_HARMONIC_ANGLE       | bool         | True        | None          | Use harmonic angle interaction for consecutuve beads ($i,i\pm 1, i\pm 2$). |
| POL_HARMONIC_ANGLE_R0        | float        | pi          | None (radians)   | Equilibrium angle of harmonic angle force. |
| POL_HARMONIC_ANGLE_CONSTANT_K| float        | 100.0       | kJ/mol/radian^2 | Harmonic angle force constant. |
| LE_USE_HARMONIC_BOND         | bool         | True        | None          | Use harmonic bond interaction for long range loops. |
| LE_FIXED_DISTANCES           | bool         | False       | None          | For fixed distances between loops. False if you want to correlate with the heatmap strength. |
| LE_HARMONIC_BOND_R0          | float        | 0.1         | None          | Harmonic bond distance equilibrium constant for long-range loops. |
| LE_HARMONIC_BOND_K           | float        | 30000.0     | kJ/mol/nm^2   | Harmonic bond force constant for long-range loops. |
| EV_USE_EXCLUDED_VOLUME       | bool         | True        | None          | Use excluded volume $V_{ev}=\epsilon\left(\frac{\sigma}{r+r_{\text{small}}}\right)^{\alpha}$. |
| EV_EPSILON                   | float        | 100.0       | kJ/mol          | Epsilon parameter - strength of excluded volume. |
| EV_R_SMALL                   | float        | 0.05        | nm          | Add something small in denominator to make it not exploding all the time. |
| EV_POWER                     | float        | 3.0         | None          | Power $\alpha$ in the exponent of EV potential. |
| SC_USE_SPHERICAL_CONTAINER   | bool         | False       | None          | Use Spherical container. |
| SC_RADIUS1                   | float        | None        | nm          | Inner spherical container radius. |
| SC_RADIUS2                   | float        | None        | nm          | Outer spherical container radius. |
| SC_SCALE                     | float        | 1000.0      | kJ/mol/nm^2   | Spherical container scaling factor |
| CHB_USE_CHROMOSOMAL_BLOCKS   | bool         | False       | None          | Use Chromosomal Blocks. |
| CHB_KC                       | float        | 0.3         | nm^(-4)       | Block copolymer width parameter. |
| CHB_DE                       | float        | 1e-5        | kJ/mol        | Energy factor for block copolymer chromosomal model. |
| COB_USE_COMPARTMENT_BLOCKS   | bool         | False       | kJ/mol        | Use Compartment Blocks. |
| COB_DISTANCE                 | float        | None        | nm        | Block copolymer equilibrium distance for chromosomal blocks. |
| COB_EA                       | float        | 1.0         | kJ/mol         | Energy strength for A compartment. |
| COB_EB                       | float        | 2.0         | kJ/mol         | Energy strength for B compartment. |
| SCB_USE_SUBCOMPARTMENT_BLOCKS| bool         | False       | None          | Use Subcompartment Blocks. |
| SCB_DISTANCE                 | float        | None        | nm          | Block copolymer equilibrium distance for chromosomal blocks. |
| SCB_EA1                      | float        | 1.0         | kJ/mol        | Energy strength for A1 compartment. |
| SCB_EA2                      | float        | 1.33        | kJ/mol        | Energy strength for A2 compartment. |
| SCB_EB1                      | float        | 1.66        | kJ/mol        | Energy strength for B1 compartment. |
| SCB_EB2                      | float        | 2.0         | kJ/mol        | Energy strength for B2 compartment. |
| IBL_USE_B_LAMINA_INTERACTION | bool         | False       | None          | Interactions of B compartment with lamina. |
| IBL_SCALE                    | float        | 400.0       | kJ/mol        | Scaling factor for B comoartment interaction with lamina. |
| CF_USE_CENTRAL_FORCE         | bool         | False       | None          | Attraction of smaller chromosomes to the nucleolus (in the center). |
| CF_STRENGTH                  | float        | 10.0       | kJ/mol        | Strength of Attraction of Smaller Chromosomes |
| NUC_DO_INTERPOLATION         | bool         | False       | None          | Attraction of smaller chromosomes. |
| NUC_RADIUS                   | float        | 0.1         | None          | The radius of the single nucleosome helix. |
| POINTS_PER_NUC               | int          | 20          | None          | The number of points that consist a nucleosome helix. |
| PHI_NORM                     | float        | pi/5        | None          | Zig zag angle. |
| SIM_RUN_MD                   | bool         | False       | None          | Do you want to run MD simulation? |
| SIM_N_STEPS                  | int          | 10000        | None          | Number of steps in MD simulation |
| SIM_ERROR_TOLERANCE          | float        | 0.01        | None          | Error tolerance for variable MD simulation |
| SIM_AMD_ALPHA                | float        | 100.0       | kJ/mol          | Alpha of AMD simulation. |
| SIM_AMD_E                    | float        | 1000.0      | kJ/mol          | E (energy) of AMD simulation. |
| SIM_SAMPLING_STEP            | int          | 100         | None          | It determines in t how many steps we save a structure. |
| SIM_INTEGRATOR_TYPE          | str          | langevin    | None          | Possible choices: variable_nagevin, langevin, variable_verlet, verlet, amd, brownian. |
| SIM_INTEGRATOR_STEP          | Quantity     | 1      | fsec  | The step of integrator. |
| SIM_FRICTION_COEFF           | float        | 0.5         | 1/psec          | Friction coefficient (Used only with langevin and brownian integrators). |
| SIM_SET_INITIAL_VELOCITIES   | bool         | False       | None          | Sets initial velocities based on Boltzmann distribution. |
| SIM_TEMPERATURE              | Quantity     | 310  | kelvin        | Simulation temperature |
| TRJ_FRAMES                   | int          | 2000        | None          | Number of trajectory frames to save. |

## Output Directory
The output directory is organized in the following folders,

```
config_auto.ini
├── md_frames
│   ├── frame_1_100.cif
├── metadata
│   ├── chimera_gene_coloring.cmd
│   ├── chrom_idxs.npy
│   ├── chrom_lengths.npy
│   ├── ds.npy
│   ├── ms.npy
│   ├── MultiMM_annealing.dcd
│   ├── MultiMM_init.cif
│   ├── MultiMM.psf
│   ├── ns.npy
│   └── parameters.txt
├── model
│   ├── MultiMM_afterMD.cif
│   └── MultiMM_minimized.cif
├── plots
│   ├── initial_structure_gene_coloring.png
│   ├── initial_structure.png
│   ├── minimized_structure_gene_coloring.png
│   ├── minimized_structure.png
│   ├── structure_afterMD_gene_coloring.png
│   └── structure_afterMD.png
```

In `md_frames`, the frames of md dynamics can be found in case that md simulation is enabled. 

In `metadata` you can find the initial structure and other produced numpy arrays. For example, `ms`, `ns`, are the left and right locations of loops in the region of interest. `ds` is the loop strength converted to distance. `psf` and `dcd` files are for the visualization of the trajectory in UCSF chimera software: https://www.cgl.ucsf.edu/chimera/, and `chimera_gene_coloring.cmd` is genetrated to give the coloring with the red region to be the gene of interest.

In the `model` is the resulted minimized structure and the structure after the MD simulation. If it is genomewide simulation, it would output the structures of each chromosome in a folder `chromosomes`.

In `plots` directory you can find plots of the structures. Note that the initial structure plot is only to see the initial stucture used in the simulation. Initial structure does not have direct biological meaning.

## Citation and Contribution

The software is freely distributed under the GNU license and is available for use in research, in accordance with the open-source license of MultiMM. If you use the software for research or other purposes, please cite the following paper:

- Korsak, Sevastianos, Krzysztof Banecki, and Dariusz Plewczynski. "Multiscale molecular modeling of chromatin with MultiMM: From nucleosomes to the whole genome." Computational and Structural Biotechnology Journal 23 (2024): 3537–3548.

If you would like to contribute to the development of this model or suggest improvements, we encourage you to contact the authors. Additionally, if you encounter any issues while running the software, your feedback is highly appreciated, and we are happy to assist you.
