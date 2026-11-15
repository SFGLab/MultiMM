# All arguments

*This page is generated from `src/multimm/config.py` by `docs_tools/gen_arguments.py`, so names, types and defaults always match the code. Descriptions are shortened here; `MultiMM --help` and the model pages give the full text.*

Types: `quantity` is a number with a unit written as text, for example `0.1 nanometer` or `310 kelvin` ([syntax](../model/scales.md)); `bool` accepts `True/False`, `yes/no`, `1/0`. "*none*" means the argument is empty by default. Arguments marked **unused** are accepted but not read anywhere in the current code.

## Run and hardware

Where and how the computation runs.

| Argument | Type | Default | Description |
|---|---|---|---|
| `PLATFORM` | str | `CPU` | Name of the platform. |
| `CPU_THREADS` | int | *none* | The number of CPU threads. |
| `DEVICE` | str | *none* | Device index for CUDA or OpenCL (count from 0). **(unused)** |

## Modelling setup

What is modelled and where results go.

| Argument | Type | Default | Description |
|---|---|---|---|
| `MODELLING_LEVEL` | str | *none* | Choose 'GENE' or 'REGION' for gene or TAD level, 'CHROM' for chromosome level, and 'GW' for genome level. |
| `N_BEADS` | int | `50000` | Number of Simulation Beads. |
| `CHROM` | str | *none* | Chromosome that corresponds the modelling region of interest. |
| `LOC_START` | int | *none* | Starting region coordinate. |
| `LOC_END` | int | *none* | Ending region coordinate. |
| `OUT_PATH` | str | `results` | Output folder name. |
| `SAVE_PLOTS` | bool | `True` | Save plots. |
| `FORCEFIELD_PATH` | str | *bundled file* | Path to XML file with forcefield. |

## Input data

Optional experimental inputs ([details](input.md)).

| Argument | Type | Default | Description |
|---|---|---|---|
| `LOOPS_PATH` | str | *none* | A .bedpe file path with loops. |
| `COMPARTMENT_PATH` | str | *none* | A .bed file with (sub)compartments from Calder (or the same format), or a .bw/.bigwig signal track from which … |
| `HIC_PATH` | str | *none* | Path to a Hi-C contact file (.hic, .cool, or .mcool). |
| `ATACSEQ_PATH` | str | *none* | A .bw or .BigWig file path with atacseq data. |
| `DOWNSAMPLING_PROB` | float | `1.0` | Probability of downsampling contacts (from 0 to 1). |
| `COMPARTMENT_FLIP_PROB` | float | `0.0` | Probability of flipping compartment identity per bead (A↔B). |
| `COMPARTMENT_NOISE_STD` | float | `0.0` | Standard deviation of Gaussian noise applied to compartment field before discretization. |

## Initial structure

Starting conformation ([model](../model/initial-structure.md)).

| Argument | Type | Default | Description |
|---|---|---|---|
| `BUILD_INITIAL_STRUCTURE` | bool | `True` | To build a new initial structure. |
| `INITIAL_STRUCTURE_TYPE` | enum | `hilbert` | You can choose between. |
| `INITIAL_STRUCTURE_PATH` | str | *none* | Path to CIF file. |

## Ensembles and randomness

Many structures, shuffling and seeds.

| Argument | Type | Default | Description |
|---|---|---|---|
| `GENERATE_ENSEMBLE` | bool | `False` | Generate an ensemble of structures instead of one. |
| `N_ENSEMBLE` | int | *none* | Number of samples of structures that you would like to calculate. |
| `SHUFFLE_CHROMS` | bool | `False` | Shuffle the chromosomes. |
| `SHUFFLING_SEED` | int | `0` | Shuffling random seed. |

## Backbone

Bonds, angles and excluded volume ([model](../model/backbone.md)).

| Argument | Type | Default | Description |
|---|---|---|---|
| `POL_USE_HARMONIC_BOND` | bool | `True` | Use harmonic bond interaction. |
| `POL_HARMONIC_BOND_R0` | quantity | `0.1 nanometer` | Harmonic bond distance equilibrium constant. |
| `POL_HARMONIC_BOND_K` | quantity | `300000.0 kilojoules_per_mole/nanometer**2` | Harmonic bond force constant (fixed unit. |
| `POL_USE_HARMONIC_ANGLE` | bool | `True` | Use harmonic angle interaction. |
| `POL_HARMONIC_ANGLE_R0` | quantity | `3.141592653589793 radian` | Harmonic angle distance equilibrium constant. |
| `POL_HARMONIC_ANGLE_CONSTANT_K` | quantity | `100.0 kilojoules_per_mole/radian**2` | Harmonic angle force constant. |
| `EV_USE_EXCLUDED_VOLUME` | bool | `True` | Use excluded volume. |
| `EV_FORCE_TYPE` | str | `powerlaw` | Excluded volume functional form. Options: powerlaw, gaussian_core. |
| `EV_EPSILON` | float | `100.0` | Epsilon parameter. |
| `EV_R_SMALL` | float | `0.05` | Add something small in denominator to make it not exploding all the time. |
| `EV_POWER` | float | `6.0` | Power in the exponent of EV potential. |

## Loops

Loop bonds from a `.bedpe` file ([model](../model/loops.md)).

| Argument | Type | Default | Description |
|---|---|---|---|
| `LE_USE_HARMONIC_BOND` | bool | `True` | Use harmonic bond interaction for long range loops. |
| `LE_LOOP_FORCE_TYPE` | str | `harmonic` | Loop extrusion bond functional form. Options: harmonic, fene_soft, gaussian_tether. |
| `LE_FIXED_DISTANCES` | bool | `False` | For fixed distances between loops. |
| `LE_HARMONIC_BOND_R0` | quantity | `0.1 nanometer` | Harmonic bond distance equilibrium constant. |
| `LE_HARMONIC_BOND_K` | quantity | `30000.0 kilojoules_per_mole/nanometer**2` | Harmonic bond force constant (fixed unit. |

## Compartments and chromosomes

A/B, subcompartment and chromosome-territory forces ([model](../model/compartments.md)).

| Argument | Type | Default | Description |
|---|---|---|---|
| `COB_USE_COMPARTMENT_BLOCKS` | bool | `False` | Use Compartment Blocks. |
| `COB_FORCE_TYPE` | str | `gaussian` | Compartment block interaction functional form. Options: gaussian, yukawa, theta. |
| `COB_DISTANCE` | quantity | *none* | Block copolymer equilibrium distance for chromosomal blocks. **(unused)** |
| `COB_EA` | float | `1.0` | Energy strength for A compartment. |
| `COB_EB` | float | `2.0` | Energy strength for B compartment. |
| `COB_CUTOFF` | float | `5.0` | Cutoff of the compartment and subcompartment forces, in bond lengths b0 (5 = 0.5 nm). |
| `SCB_USE_SUBCOMPARTMENT_BLOCKS` | bool | `False` | Use Subcompartment Blocks. |
| `SCB_FORCE_TYPE` | str | `gaussian` | Subcompartment block interaction functional form. Options: gaussian, yukawa, theta. |
| `SCB_DISTANCE` | quantity | *none* | Block copolymer equilibrium distance for chromosomal blocks. **(unused)** |
| `SCB_EA1` | float | `1.0` | Energy strength for A1 compartment. |
| `SCB_EA2` | float | `1.33` | Energy strength for A2 compartment. |
| `SCB_EB1` | float | `1.66` | Energy strength for B1 compartment. |
| `SCB_EB2` | float | `2.0` | Energy strength for B2 compartment. |
| `CHB_USE_CHROMOSOMAL_BLOCKS` | bool | `False` | Use Chromosomal Blocks. |
| `CHB_FORCE_TYPE` | str | `polynomial` | Chromosome self-attraction kernel controlling global compaction into globules. Options: polynomial, gaussian, saturating. |
| `CHB_KC` | float | `0.3` | Block copolymer width parameter. |
| `CHB_DE` | float | `1e-04` | Energy factor for block copolymer chromosomal model. |

## Hi-C force

Restraints from a Hi-C map ([model](../model/hic-force.md)).

| Argument | Type | Default | Description |
|---|---|---|---|
| `HIC_USE_FORCE` | bool | `False` | Apply Hi-C contact-guided force to the simulation. |
| `HIC_NORMALIZATION` | str | `KR` | Hi-C matrix normalisation method. Options: KR, VC, VC_SQRT, NONE. |
| `HIC_K_SCALE` | float | `40.0` | Global energy scale for the Hi-C force [kJ/mol]. |
| `HIC_RC` | float | *none* | Explicit contact-radius scale [nm] for the Hi-C force's own target-distance mapping. |
| `HIC_AUTO_SCALE` | bool | `True` | If True (default, recommended), recalibrate the Hi-C force's contact-radius scale (rc) from the median (50th … |
| `HIC_BOLTZMANN_KERNEL` | str | `exponential` | P(r) shape assumed for the Hi-C Boltzmann-PMF's equilibrium pair-distance distribution, i.e. the c_ij -> … Options: exponential, power_law, sigmoid. |
| `HIC_BOLTZMANN_ALPHA` | float | `4.0` | Hi-C scaling-law exponent converting contact strength to a target 3-D distance via classic Boltzmann-inversion. |
| `HIC_FORCE_OE` | bool | `False` | If False (default), target raw (KR-balanced) contact frequency directly. |
| `HIC_NOISE_INTENSITY` | float | `0.0` | Std-dev of Gaussian noise (same [0, 1] scale as c_ij) added to every Hi-C bond's contact strength, redrawn … |
| `HIC_MAX_GAP` | int | `10` | Longest run of empty Hi-C bins (in bins) filled by interpolation. |
| `HIC_INSULATION_WINDOW` | int | `10` | Half-width (beads) of the sliding square used by Hi-C validation's insulation-score metric. |
| `HIC_BLOCK_COPOLYMER` | bool | `False` | Opt-in. |
| `HIC_BLOCK_COPOLYMER_MIN_BP` | float | `5_000_000` | Minimum modelled region size (bp) for HIC_BLOCK_COPOLYMER to stay enabled. |

## Nuclear organisation

Container, lamina and nucleolus ([model](../model/lamina.md)).

| Argument | Type | Default | Description |
|---|---|---|---|
| `SC_USE_SPHERICAL_CONTAINER` | bool | `False` | Use Spherical container. |
| `SC_SCALE` | float | `1000.0` | Spherical container scaling factor. |
| `SC_RADIUS1` | quantity | *none* | Spherical container radius,. **(unused)** |
| `SC_RADIUS2` | quantity | *none* | Spherical container radius,. **(unused)** |
| `IBL_USE_B_LAMINA_INTERACTION` | bool | `False` | Interactions of B compartment with lamina. |
| `BLAMINA_FORCE_TYPE` | str | `sin` | B-lamina interaction functional form. Options: sin, gaussian_shell, harmonic_shell, logistic_shell. |
| `IBL_SCALE` | float | `400.0` | Scaling factor for B compartment interaction with lamina. |
| `CF_USE_CENTRAL_FORCE` | bool | `False` | Attraction of smaller chromosomes. |
| `CENTRAL_FORCE_TYPE` | str | `harmonic` | Central nucleolar attraction functional form controlling radial bias toward nucleus center. Options: harmonic, gaussian, logistic. |
| `CF_STRENGTH` | float | `20.0` | Strength of Interaction. |

## Gene mode

Modelling around a gene ([advanced](../model/advanced.md)).

| Argument | Type | Default | Description |
|---|---|---|---|
| `GENE_TSV` | str | *bundled file* | A .tsv with genes and their locations in the genome. |
| `GENE_NAME` | str | *none* | The name of the gene of interest. |
| `GENE_ID` | str | *none* | The id of the gene of interest. |
| `GENE_WINDOW` | int | `100000` | The window around of the area around the gene of interest. |

## Nucleosome interpolation

Post-processing with nucleosomes ([advanced](../model/advanced.md)).

| Argument | Type | Default | Description |
|---|---|---|---|
| `NUC_DO_INTERPOLATION` | bool | `False` | Attraction of smaller chromosomes. |
| `MAX_NUCS_PER_BEAD` | int | `4` | Maximum amount of nucleosomes per single bead. |
| `NUC_RADIUS` | float | `0.1` | The radius of the single nucleosome helix. |
| `POINTS_PER_NUC` | int | `20` | The number of points that consist a nucleosome helix. |
| `PHI_NORM` | float | `0.6283185307179586` | Zig zag angle. |

## Molecular dynamics

Optional MD after minimisation.

| Argument | Type | Default | Description |
|---|---|---|---|
| `SIM_RUN_MD` | bool | `False` | Do you want to run MD simulation? |
| `SIM_N_STEPS` | int | `10000` | Number of steps in MD simulation. |
| `TRJ_FRAMES` | int | None | *none* | Number of CIF frames to save during MD. |
| `SIM_SAMPLING_STEP` | int | `100` | Save a structure every this many MD steps. |
| `SIM_INTEGRATOR_TYPE` | str | `langevin` | Alternative. |
| `SIM_INTEGRATOR_STEP` | quantity | `1 femtosecond` | The step of integrator. |
| `SIM_FRICTION_COEFF` | float | `0.5` | Friction coefficient. |
| `SIM_TEMPERATURE` | quantity | `310 kelvin` | Simulation temperature. |
| `SIM_SET_INITIAL_VELOCITIES` | bool | `True` | Initialize the MD simulation with a random initial velocity field, drawn from the Maxwell-Boltzmann … |
| `SIM_ERROR_TOLERANCE` | float | `0.01` | Error tolerance for variable MD simulation. |
| `SIM_AMD_ALPHA` | float | `100.0` | Alpha of AMD simulation. |
| `SIM_AMD_E` | float | `1000.0` | E (energy) of AMD simulation. |
