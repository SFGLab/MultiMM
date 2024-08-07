# MultiEM: An OpenMM-based software for the whole genome 3D structure reconstruction
MultiMM is an OpenMM model for the modelling of the whole genome structure. What makes it different from other models is that it is multiscale, which means that it aims to model different scales of chromatin starting from smaller scales (nucleosomes) to the scale of chromosomal territories. The algorithm is fast and accurate. Key point for the fast modelling is the GPU parallelization of the OpenMM and smart assumptions that assist the optimizer to find the global minimum. The most fundamental of these assumptions, is the assumption of Hilbert curve as an initial structure. This helps MultiMM to minimize faster because the initial structure is already highly compacted.

![GW_em](https://github.com/user-attachments/assets/3d019616-2c0f-4bfc-a792-d87fcdc0d96a)

Having run MultiMM model, the user obtains a genomewide structure. The user can color different chromosomes, or different compartments, or they can even vizualize seprately chromosomes. The usage of MultiEM is very simple since it is based on the modification of a single configuration file. Here we exaplain its usage.

![MultiEM_scales](https://github.com/user-attachments/assets/a0d14ddc-41bf-4d14-8a8c-aa418fe575b5)


## About Operating System

MultiMM is tested mainly in Linux-based operating systems. It has been tested succesfully in Ubuntu, Debian and Red-Hat based distributions. It is possible also to run it on MacOS, but without CUDA support, which is very helpful for the acceleration of the computations. We suggest the user to not run MultiMM on Windows computers.

## Installation
Create a python 3.10 environment and type,

```
pip install -r requirements.txt
```

## Input Data

MultiMM model relies on three types of datasets:

* Loop interactions in `bedpe` format (mandatory).
* Compartmentalization data in `bed` format (optional).
* ATAC-Seq p-value data in `.BigWig ` format (optional).

For **loop interactions**, the user needs to provide a fie with interactions between anchor 1 and anchor 2, and their strength. Therefore, the file should be in `.bedpe` format, it should have 7 columns, it should not contain a header, and it should look like that,

```
chr10	100225000	100230000	chr10	100420000	100425000	95
chr10	100225000	100230000	chr10	101005000	101010000	56
chr10	101190000	101195000	chr10	101370000	101375000	152
chr10	101190000	101200000	chr10	101470000	101480000	181
chr10	101600000	101605000	chr10	101805000	101810000	152
```

The file may contain all chromosomes, and MultiMM can model them automatically. In case of single cell data, the user can probide a file with the second and third column identical (the same for fifth and sixth column), and have everywhere strength 1.

For **(sub)compartment interactions**, the file should look like the one produced from CALDER software: https://github.com/CSOgroup/CALDER2. It does not have to be called by CALDER but it should be in the same format. Specifically, the user needs to provide at least the first four columns of the file, containing chromosome, regions, and the subcompartments label. Therefore, it should look like this,

```
chr1	700001	900000	A.1.2.2.2.2.2	0.875	.	700001	900000	#FF4848
chr1	900001	1400000	A.1.1.1.1.2.1.1.1.1.1	1	.	900001	1400000	#FF0000
chr1	1400001	1850000	A.1.1.1.1.2.1.2.2.2.1	1	.	1400001	1850000	#FF0000
chr1	1850001	2100000	B.1.1.2.2.1.2.1	0.5	.	1850001	2100000	#DADAFF
```

For ATAC-Seq data the user should provide a file with p-value in BigWig format. It is needed to have the library pyBigWig which does not work in Windows operating systems.

**Attention!** For now MultiMM works only for human genome data. The code probaby can run for other organisms as well with a little debugging and modifications.

## Usage
All the parameters of the model are saved within a `config.ini` file. This file should have the following form,

```
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

; Enable forcefield for GW simulation
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

Having specified the paramameters and the forces that you would like to use, you can run on terminal the following command,

```
python run.py -c config.ini
```

The sotware will return you a folder with the resulting structure, and some plots that show how compartments are distributed.

Example data can be found in the Google Drive: https://drive.google.com/drive/folders/1nFAPE4pCaHpeL5nw6nq0VvfUFoc24aXm?usp=sharing. Please note that these data are not produced in our laboratory, and they are publicly available data from Rao et al, with predicted subcompartments with CALDER. The ATAC-Seq data are from ENCODE. 

In `examples` folder we provide example configuration files that can be used in different modelling cases.

## The long table of simulation arguments
There is a big amount of parameters defined in the aguments of MultiMM. Here we provide the default values and a description of each one of them. The default values are tested to work with the genome-wide simulation, and they can be changes in the configuration file, if they do not satisfy user's need. In general it is not needed to specify the units, but it is important to know what kind of OpenMM units are assumed.

Note that despite the fact that (sub)compartment forcefields are by default disabled, they should be enabled in the configuration file, in case that compartmentalization data would be provided. Similarly, for 

| Argument Name                | Type         | Value       | Units         | Description |
|------------------------------|--------------|-------------|---------------|-------------|
| PLATFORM                     | str          | None        | None          | name of the platform. Available choices: CPU, OpenCL, CUDA. |
| DEVICE                       | str          | None        | None          | device index for CUDA or OpenCL (count from 0) |
| INITIAL_STRUCTURE_PATH       | str          | None        | None          | Path to CIF file. |
| BUILD_INITIAL_STRUCTURE      | bool         | True        | None          | To build a new initial structure. |
| INITIAL_STRUCTURE_TYPE       | str          | hilbert     | None          | you can choose between: hilbert, circle. |
| FORCEFIELD_PATH              | str          | None        | None          | Path to XML file with forcefield. |
| N_BEADS                      | int          | 50000       | None          | Number of Simulation Beads. |
| COMPARTMENT_PATH             | str          | None        | None          | It should be  a `.bed` file with subcompartments from Calder. |
| LOOPS_PATH                   | str          | None        | None          | A `.bedpe` file path with loops. It is required. |
| ATACSEQ_PATH                 | str          | None        | None          | A `.bw` or `.BigWig` file path with atacseq data for nucleosome interpolation. It is not required. |
| OUT_PATH                     | str          | results     | None          | Output folder name. |
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
| CHB_DE                       | float        | 1e-3        | kJ/mol        | Energy factor for block copolymer chromosomal model. |
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
| IBL_SCALE                    | float        | 100.0       | kJ/mol        | Scaling factor for B comoartment interaction with lamina. |
| CF_USE_CENTRAL_FORCE         | bool         | False       | None          | Attraction of smaller chromosomes to the nucleolus (in the center). |
| CF_STRENGTH                  | float        | 100.0       | kJ/mol        | Strength of Interaction |
| NUC_DO_INTERPOLATION         | bool         | False       | None          | Attraction of smaller chromosomes. |
| MAX_NUCS_PER_BEAD            | int          | 4           | None          | Maximum amount of nucleosomes per single bead for nucleosome interpolation. |
| NUC_RADIUS                   | float        | 0.1         | None          | The radius of the single nucleosome helix. |
| POINTS_PER_NUC               | int          | 20          | None          | The number of points that consist a nucleosome helix. |
| PHI_NORM                     | float        | pi/5        | None          | Zig zag angle. |
| SIM_RUN_MD                   | bool         | False       | None          | Do you want to run MD simulation? |
| SIM_N_STEPS                  | int          | None        | None          | Number of steps in MD simulation |
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

## Copyrights

The software is freely distributed and everybody can use it how they want, improve it, or apply it to their research interests. In case that the software would be used for research, we would like that you will cite our paper:

- Korsak, Sevastianos, Krzysztof Banecki, and Dariusz Plewczynski. "Multiscale Molecular Modelling of Chromatin with MultiMM: From Nucleosomes to the Whole Genome." bioRxiv (2024): 2024-07.

Please, communicate with the authors in case that you would like to contribute in this model, and you would like to improve it.
