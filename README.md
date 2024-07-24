# MultiEM: An OpenMM-based software for the whole genome 3D structure reconstruction
MultiMM is an OpenMM model for the modelling of the whole genome structure. What makes it different from other models is that it is multiscale, which means that it aims to model different scales of chromatin starting from smaller scales (nucleosomes) to the scale of chromosomal territories. The algorithm is fast and accurate. Key point for the fast modelling is the GPU parallelization of the OpenMM and smart assumptions that assist the optimizer to find the global minimum. The most fundamental of these assumptions, is the assumption of Hilbert curve as an initial structure. This helps MultiMM to minimize faster because the initial structure is already highly compacted.

![GW_em_lowres](https://github.com/user-attachments/assets/13a432e7-8412-420b-bdfb-c3eef1410e9b)

Having run MultiMM model, the user obtains a genomewide structure. The user can color different chromosomes, or different compartments, or they can even vizualize seprately chromosomes. The usage of MultiEM is very simple since it is based on a single configuration file. In this documentation we show examples of its usage.

## Required datasets

MultiMM model relies on three types of datasets:

* Loop interactions in `bedpe` format (mandatory).
* Compartmentalization data in `bed` format (optional).
* ATAC-Seq p-value data in `.BigWig ` format (optional).

For calling subcompartments, we suggest our users to use Calder software: https://github.com/CSOgroup/CALDER2.

## About Operating System

MultiMM is tested mainly in Linux-based operating systems. It has been tested succesfully in Ubuntu, Debian and Red-Hat based distributions. It is possible also to run it on MacOS, but without CUDA support, which is very helpful for the acceleration of the computations. We suggest the user to not run MultiMM on Windows computers.

## Installation
Create a python 3.10 environment and type,

```
pip install -r requirements.txt
```

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
NUC_DO_INTERPOLATION = True
```

Having specified the paramameters and the forces that you would like to use, you can run on terminal the following command,

```
python run.py -c config.ini
```

The sotware will return you a folder with the resulting structure, and some plots that show how compartments are distributed.

## The long table of simulation arguments
There is a big amount of parameters defined in the aguments of MultiMM. Here we provide the default values and a description of each one of them. The default values are tested to work with the genome-wide simulation, and they can be changes in the configuration file, if they do not satisfy user's need.

| Argument Name                | Type         | Value       | Units         | Description |
|------------------------------|--------------|-------------|---------------|-------------|
| PLATFORM                     | str          | None        | None          | name of the platform. Available choices: {' '.join(available_platforms)} |
| DEVICE                       | str          | None        | None          | device index for CUDA or OpenCL (count from 0) |
| INITIAL_STRUCTURE_PATH       | str          | None        | None          | Path to CIF file. |
| BUILD_INITIAL_STRUCTURE      | bool         | True        | None          | To build a new initial structure. |
| INITIAL_STRUCTURE_TYPE       | str          | hilbert     | None          | you can choose between: hilbert, circle. |
| FORCEFIELD_PATH              | str          | None        | None          | Path to XML file with forcefield. |
| N_BEADS                      | int          | 50000       | None          | Number of Simulation Beads. |
| COMPARTMENT_PATH             | str          | None        | None          | It can be either .bed file with subcompartments from Calder or .BigWig signal. |
| LOOPS_PATH                   | str          | None        | None          | A .bedpe file path with loops. It is required. |
| ATACSEQ_PATH                 | str          | None        | None          | A .bw or .BigWig file path with atacseq data. It is not required. |
| OUT_PATH                     | str          | results     | None          | Output folder name. |
| LOC_START                    | int          | None        | None          | Starting region coordinate. |
| LOC_END                      | int          | None        | None          | Ending region coordinate. |
| CHROM                        | str          | None        | None          | Chromosome that corresponds the the modelling region of interest (in case that you do not want to model the whole genome). |
| SHUFFLE_CHROMS               | bool         | False       | None          | Shuffle the chromosomes. |
| SHUFFLING_SEED               | int          | 0           | None          | Shuffling random seed. |
| SAVE_PLOTS                   | bool         | True        | None          | Save plots. |
| POL_USE_HARMONIC_BOND        | bool         | True        | None          | Use harmonic bond interaction. |
| POL_HARMONIC_BOND_R0         | float        | 0.1         | None          | harmonic bond distance equilibrium constant |
| POL_HARMONIC_BOND_K          | float        | 300000.0    | kJ/mol/nm^2   | harmonic bond force constant |
| POL_USE_HARMONIC_ANGLE       | bool         | True        | None          | Use harmonic angle interaction. |
| POL_HARMONIC_ANGLE_R0        | float        | pi          | None          | harmonic angle distance equilibrium constant |
| POL_HARMONIC_ANGLE_CONSTANT_K| float        | 100.0       | kJ/mol/radian^2 | harmonic angle force constant |
| LE_USE_HARMONIC_BOND         | bool         | True        | None          | Use harmonic bond interaction for long range loops. |
| LE_FIXED_DISTANCES           | bool         | False       | None          | For fixed distances between loops. False if you want to correlate with the hatmap strength. |
| LE_HARMONIC_BOND_R0          | float        | 0.1         | None          | harmonic bond distance equilibrium constant |
| LE_HARMONIC_BOND_K           | float        | 30000.0     | kJ/mol/nm^2   | harmonic bond force constant |
| EV_USE_EXCLUDED_VOLUME       | bool         | True        | None          | Use excluded volume. |
| EV_EPSILON                   | float        | 100.0       | None          | Epsilon parameter. |
| EV_R_SMALL                   | float        | 0.05        | None          | Add something small in denominator to make it not exploding all the time. |
| EV_POWER                     | float        | 3.0         | None          | Power in the exponent of EV potential. |
| SC_USE_SPHERICAL_CONTAINER   | bool         | False       | None          | Use Spherical container |
| SC_RADIUS1                   | float        | None        | None          | Spherical container radius, |
| SC_RADIUS2                   | float        | None        | None          | Spherical container radius, |
| SC_SCALE                     | float        | 1000.0      | None          | Spherical container scaling factor |
| CHB_USE_CHROMOSOMAL_BLOCKS   | bool         | False       | None          | Use Chromosomal Blocks. |
| CHB_KC                       | float        | 0.3         | None          | Block copolymer width parameter. |
| CHB_DE                       | float        | 1e-3        | None          | Energy factor for block copolymer chromosomal model. |
| COB_USE_COMPARTMENT_BLOCKS   | bool         | False       | None          | Use Compartment Blocks. |
| COB_DISTANCE                 | float        | None        | None          | Block copolymer equilibrium distance for chromosomal blocks. |
| COB_EA                       | float        | 1.0         | None          | Energy strength for A compartment. |
| COB_EB                       | float        | 2.0         | None          | Energy strength for B compartment. |
| SCB_USE_SUBCOMPARTMENT_BLOCKS| bool         | False       | None          | Use Subcompartment Blocks. |
| SCB_DISTANCE                 | float        | None        | None          | Block copolymer equilibrium distance for chromosomal blocks. |
| SCB_EA1                      | float        | 1.0         | None          | Energy strength for A1 compartment. |
| SCB_EA2                      | float        | 1.33        | None          | Energy strength for A2 compartment. |
| SCB_EB1                      | float        | 1.66        | None          | Energy strength for B1 compartment. |
| SCB_EB2                      | float        | 2.0         | None          | Energy strength for B2 compartment. |
| IBL_USE_B_LAMINA_INTERACTION | bool         | False       | None          | Interactions of B compartment with lamina. |
| IBL_SCALE                    | float        | 100.0       | None          | Scaling factor for B comoartment interaction with lamina. |
| CF_USE_CENTRAL_FORCE         | bool         | False       | None          | Attraction of smaller chromosomes. |
| CF_STRENGTH                  | float        | 100.0       | None          | Strength of Interaction |
| NUC_DO_INTERPOLATION         | bool         | False       | None          | Attraction of smaller chromosomes. |
| MAX_NUCS_PER_BEAD            | int          | 4           | None          | Maximum amount of nucleosomes per single bead. |
| NUC_RADIUS                   | float        | 0.1         | None          | The radius of the single nucleosome helix. |
| POINTS_PER_NUC               | int          | 20          | None          | The number of points that consist a nucleosome helix. |
| PHI_NORM                     | float        | pi/5        | None          | Zig zag angle. |
| SIM_RUN_MD                   | bool         | False       | None          | Do you want to run MD simulation? |
| SIM_N_STEPS                  | int          | None        | None          | Number of steps in MD simulation |
| SIM_ERROR_TOLERANCE          | float        | 0.01        | None          | Error tolerance for variable MD simulation |
| SIM_AMD_ALPHA                | float        | 100.0       | None          | Alpha of AMD simulation. |
| SIM_AMD_E                    | float        | 1000.0      | None          | E (energy) of AMD simulation. |
| SIM_SAMPLING_STEP            | int          | 100         | None          | It determines in t how many steps we save a structure. |
| SIM_INTEGRATOR_TYPE          | str          | langevin    | None          | Alternative: langevin, verlet |
| SIM_INTEGRATOR_STEP          | Quantity     | 1 femtosecond | femtosecond  | The step of integrator. |
| SIM_FRICTION_COEFF           | float        | 0.5         | None          | Friction coefficient (Used only with langevin integrator) |
| SIM_SET_INITIAL_VELOCITIES   | bool         | False       | None          | Sets initial velocities based on Boltzmann distribution |
| SIM_TEMPERATURE              | Quantity     | 310 kelvin  | kelvin        | Simulation temperature |
| TRJ_FRAMES                   | int          | 2000        | None          | Number of trajectory frames to save. |


## Copyrights
Please cite our article if you would like to base your research in this software.
