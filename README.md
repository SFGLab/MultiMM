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
PLATFORM = CUDA

; Input data
FORCEFIELD_PATH = forcefields/ff.xml
LOOPS_PATH = /mnt/raid/data/Trios/ChiA-PiPE_Loops/loops_pet3+/HG00733_PUR_C_CTCF_1mb_pet3.bedpe
COMPARTMENT_PATH = /mnt/raid/data/Trios/calder_HiChIP_subcomp/PUR_d.bed
OUT_PATH = test

; Simulation Parameters
N_BEADS = 50000

; Enable forcefield for GW simulation
SC_USE_SPHERICAL_CONTAINER = True
CHB_USE_CHROMOSOMAL_BLOCKS = True
SCB_USE_SUBCOMPARTMENT_BLOCKS = True
IAL_USE_A_LAMINA_INTERACTION = False
IBL_USE_B_LAMINA_INTERACTION = True
CF_USE_CENTRAL_FORCE = True
```

Having specified the paramameters that you would like to use, you can run on terminal the following command,

```
python run.py -c config.ini
```

The sotware will return you a folder with the resulting structure, and some plots that show how compartments are distributed.

## Copyrights
Please cite our article if you would like to base your research in this software.
