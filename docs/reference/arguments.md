# All arguments

*This page is generated from `src/multimm/config.py` by `docs_tools/gen_arguments.py`, so names, types, defaults and descriptions always match the code.*

Types: `quantity` is a number with a unit written as text, for example `0.1 nanometer` or `310 kelvin` ([syntax](../model/scales.md)); `bool` accepts `True/False`, `yes/no`, `1/0`. "*none*" means the argument is empty by default. Arguments marked **unused** are accepted but not read anywhere in the current code.

## Run and hardware

Where and how the computation runs.

| Argument | Type | Default | Description |
|---|---|---|---|
| `PLATFORM` | str | `CPU` | name of the platform. Available choices: Reference CPU |
| `CPU_THREADS` | int | *none* | The number of CPU threads (in case you would like to specify them). |
| `DEVICE` | str | *none* | device index for CUDA or OpenCL (count from 0) **(unused)** |

## Modelling setup

What is modelled and where results go.

| Argument | Type | Default | Description |
|---|---|---|---|
| `MODELLING_LEVEL` | str | *none* | Choose 'GENE' or 'REGION' for gene or TAD level, 'CHROM' for chromosome leve, and 'GW' for genome level. It will setup some parameters for you and print you helpful comments. |
| `N_BEADS` | int | `50000` | Number of Simulation Beads. |
| `CHROM` | str | *none* | Chromosome that corresponds the the modelling region of interest (in case that you do not want to model the whole genome). |
| `LOC_START` | int | *none* | Starting region coordinate. |
| `LOC_END` | int | *none* | Ending region coordinate. |
| `OUT_PATH` | str | `results` | Output folder name. |
| `SAVE_PLOTS` | bool | `True` | Save plots. |
| `FORCEFIELD_PATH` | str | *bundled file* | Path to XML file with forcefield. |

## Input data

Optional experimental inputs ([details](input.md)).

| Argument | Type | Default | Description |
|---|---|---|---|
| `LOOPS_PATH` | str | *none* | A .bedpe file path with loops. Optional — simulation runs without loops if not provided. |
| `COMPARTMENT_PATH` | str | *none* | A .bed file with (sub)compartments from Calder (or the same format), or a .bw/.bigwig signal track (e.g. an eigenvector/PC1 track) from which A/B compartment calls are derived. |
| `HIC_PATH` | str | *none* | Path to a Hi-C contact file (.hic, .cool, or .mcool). Required when HIC_USE_FORCE=True. |
| `ATACSEQ_PATH` | str | *none* | A .bw or .BigWig file path with atacseq data. It is not required. |
| `DOWNSAMPLING_PROB` | float | `1.0` | Probability of downsampling contacts (from 0 to 1). |
| `COMPARTMENT_FLIP_PROB` | float | `0.0` | Probability of flipping compartment identity per bead (A↔B). Applied after BED parsing. |
| `COMPARTMENT_NOISE_STD` | float | `0.0` | Standard deviation of Gaussian noise applied to compartment field before discretization. |

## Initial structure

Starting conformation ([model](../model/initial-structure.md)).

| Argument | Type | Default | Description |
|---|---|---|---|
| `BUILD_INITIAL_STRUCTURE` | bool | `True` | To build a new initial structure. |
| `INITIAL_STRUCTURE_TYPE` | enum | `hilbert` | you can choose between: hilbert, circle, rw, confined_rw, knot, self_avoiding_rw, spiral, sphere. |
| `INITIAL_STRUCTURE_PATH` | str | *none* | Path to CIF file. |

## Ensembles and randomness

Many structures, shuffling and seeds.

| Argument | Type | Default | Description |
|---|---|---|---|
| `GENERATE_ENSEMBLE` | bool | `False` | Default value: false. True in case that you would like to have an ensemble of structures instead of one. Better to disable it for large simulations that require long computational time. Moreover it is better to start random walk initial structure in case of true value. |
| `N_ENSEMBLE` | int | *none* | Number of samples of structures that you would like to calculate. |
| `SHUFFLE_CHROMS` | bool | `False` | Shuffle the chromosomes. |
| `SHUFFLING_SEED` | int | `0` | Shuffling random seed. |

## Backbone

Bonds, angles and excluded volume ([model](../model/backbone.md)).

| Argument | Type | Default | Description |
|---|---|---|---|
| `POL_USE_HARMONIC_BOND` | bool | `True` | Use harmonic bond interaction. |
| `POL_HARMONIC_BOND_R0` | quantity | `0.1 nanometer` | harmonic bond distance equilibrium constant |
| `POL_HARMONIC_BOND_K` | quantity | `300000.0 kilojoules_per_mole/nanometer**2` | harmonic bond force constant (fixed unit: kJ/mol/nm^2) |
| `POL_USE_HARMONIC_ANGLE` | bool | `True` | Use harmonic angle interaction. |
| `POL_HARMONIC_ANGLE_R0` | quantity | `3.141592653589793 radian` | harmonic angle distance equilibrium constant |
| `POL_HARMONIC_ANGLE_CONSTANT_K` | quantity | `100.0 kilojoules_per_mole/radian**2` | harmonic angle force constant (fixed unit: kJ/mol/radian^2) |
| `EV_USE_EXCLUDED_VOLUME` | bool | `True` | Use excluded volume. |
| `EV_FORCE_TYPE` | str | `powerlaw` | Excluded volume functional form. Options: powerlaw (default), gaussian_core. |
| `EV_EPSILON` | float | `100.0` | Epsilon parameter. |
| `EV_R_SMALL` | float | `0.05` | Add something small in denominator to make it not exploding all the time. |
| `EV_POWER` | float | `6.0` | Power in the exponent of EV potential. |

## Loops

Loop bonds from a `.bedpe` file ([model](../model/loops.md)).

| Argument | Type | Default | Description |
|---|---|---|---|
| `LE_USE_HARMONIC_BOND` | bool | `True` | Use harmonic bond interaction for long range loops. |
| `LE_LOOP_FORCE_TYPE` | str | `harmonic` | Loop extrusion bond functional form. Options: harmonic (default), fene_soft, gaussian_tether. |
| `LE_FIXED_DISTANCES` | bool | `False` | For fixed distances between loops. False if you want to correlate with the hatmap strength. |
| `LE_HARMONIC_BOND_R0` | quantity | `0.1 nanometer` | harmonic bond distance equilibrium constant |
| `LE_HARMONIC_BOND_K` | quantity | `30000.0 kilojoules_per_mole/nanometer**2` | harmonic bond force constant (fixed unit: kJ/mol/nm^2) |

## Compartments and chromosomes

A/B, subcompartment and chromosome-territory forces ([model](../model/compartments.md)).

| Argument | Type | Default | Description |
|---|---|---|---|
| `COB_USE_COMPARTMENT_BLOCKS` | bool | `False` | Use Compartment Blocks. |
| `COB_FORCE_TYPE` | str | `gaussian` | Compartment block interaction functional form. Options: gaussian (default), yukawa, theta. |
| `COB_DISTANCE` | quantity | *none* | Block copolymer equilibrium distance for chromosomal blocks. **(unused)** |
| `COB_EA` | float | `1.0` | Energy strength for A compartment. |
| `COB_EB` | float | `2.0` | Energy strength for B compartment. |
| `SCB_USE_SUBCOMPARTMENT_BLOCKS` | bool | `False` | Use Subcompartment Blocks. |
| `SCB_FORCE_TYPE` | str | `gaussian` | Subcompartment block interaction functional form. Options: gaussian (default), yukawa, theta. |
| `SCB_DISTANCE` | quantity | *none* | Block copolymer equilibrium distance for chromosomal blocks. **(unused)** |
| `SCB_EA1` | float | `1.0` | Energy strength for A1 compartment. |
| `SCB_EA2` | float | `1.33` | Energy strength for A2 compartment. |
| `SCB_EB1` | float | `1.66` | Energy strength for B1 compartment. |
| `SCB_EB2` | float | `2.0` | Energy strength for B2 compartment. |
| `CHB_USE_CHROMOSOMAL_BLOCKS` | bool | `False` | Use Chromosomal Blocks. |
| `CHB_FORCE_TYPE` | str | `polynomial` | Chromosome self-attraction kernel controlling global compaction into globules. Options: polynomial (default, handcrafted potential), gaussian (soft collapse kernel), saturating (soft-core bounded attraction). |
| `CHB_KC` | float | `0.3` | Block copolymer width parameter. |
| `CHB_DE` | float | `1e-04` | Energy factor for block copolymer chromosomal model. |

## Hi-C force

Restraints from a Hi-C map ([model](../model/hic-force.md)).

| Argument | Type | Default | Description |
|---|---|---|---|
| `HIC_USE_FORCE` | bool | `False` | Apply Hi-C contact-guided force to the simulation. |
| `HIC_NORMALIZATION` | str | `KR` | Hi-C matrix normalisation method. Options: KR (default), VC, VC_SQRT, NONE. |
| `HIC_K_SCALE` | float | `40.0` | Global energy scale for the Hi-C force [kJ/mol]. Recommended range: 20–80 (default 40); higher values (up to ~160) consistently improved insulation-score validation in testing with no diagonal-decay cost, for N_BEADS >= 300 -- below that, gains are less consistent, so lower values may suit small systems better. Above 200 triggers a runtime warning. Each pair's force is already weighted by its own observed contact strength c_ij, so weak-evidence pairs stay proportionally soft without a separate knob. |
| `HIC_RC` | float | *none* | Explicit contact-radius scale [nm] for the Hi-C force's own target-distance mapping — independent of r_comp (the compartment/subcompartment block-copolymer force's interaction range, set in set_radiuses() from nucleus geometry). The two forces are physically unrelated and are no longer tied to the same value. Default None: the Hi-C force picks its own scale instead — auto-calibrated from the initial structure's pairwise distances when HIC_AUTO_SCALE is True (recommended, default), or a nucleus-scale fallback otherwise. Set this explicitly to pin the Hi-C force's contact radius to a fixed value and skip auto-calibration entirely, regardless of HIC_AUTO_SCALE. |
| `HIC_AUTO_SCALE` | bool | `True` | If True (default, recommended), recalibrate the Hi-C force's contact-radius scale (rc) from the median (50th percentile) of the actual initial structure's pairwise-distance distribution (hic_force.auto_contact_scale) instead of a fixed default, which can leave the force with no gradient at realistic bead separations and barely change the structure from its initial state. The 50th-percentile choice is fixed, not a separate tunable knob: it's simply the right scale for most contacted pairs — not just the closest — to start within the force's reach. |
| `HIC_BOLTZMANN_KERNEL` | str | `exponential` | P(r) shape assumed for the Hi-C Boltzmann-PMF's equilibrium pair-distance distribution, i.e. the c_ij -> r_target inversion (see hic_force.get_boltzmann_p_func / _boltzmann_r_target). Options: exponential (default — classic Boltzmann distribution, c ~ exp(-r/lambda), lambda = rc/alpha), power_law (Hi-C scaling law, c ~ r^-alpha), sigmoid (bounded logistic contact probability, steepness set by alpha/rc). All three are exact functional inverses of their own P(r), reuse HIC_BOLTZMANN_ALPHA as their steepness knob, and are restrained with the same harmonic well — only the strength<->distance mapping's shape changes. |
| `HIC_BOLTZMANN_ALPHA` | float | `4.0` | Hi-C scaling-law exponent converting contact strength to a target 3-D distance via classic Boltzmann-inversion: r_target = r_min * c_ij^(-1/alpha), clipped to [r_min, rc] (see hic_force.build_boltzmann_force), then restrained there with a harmonic well weighted by c_ij — a genuine two-sided restraint (pulls if farther than its target, pushes if closer). Typical literature range: 3-4; higher values make the strength→distance mapping steeper (weak contacts hit the rc cap sooner), for every HIC_BOLTZMANN_KERNEL. |
| `HIC_FORCE_OE` | bool | `False` | If False (default), target raw (KR-balanced) contact frequency directly. If True, target enrichment-above-background (OE ratio minus 1, floored at 0) instead: pairs at/below the expected distance-decay baseline get zero force weight ('loose'); only enriched pairs attract. OE scores higher on PC1/OE-Pearson metrics in testing, but raw frequency gave better overall results in practice (diagonal decay + insulation), hence the default. |
| `HIC_NOISE_INTENSITY` | float | `0.0` | Std-dev of Gaussian noise (same [0, 1] scale as c_ij) added to every Hi-C bond's contact strength, redrawn once per saved MD frame around its ORIGINAL data-derived value (never drifting cumulatively). This lets different contacts take turns pulling strongest from frame to frame, nudging the structure to explore nearby configurations instead of settling into one exact attractor, while every draw stays anchored to the real data. 0 (default) disables it; try 0.05-0.2 for mild exploration. Only applied during MD (SIM_RUN_MD=True). |
| `HIC_MAX_GAP` | int | `10` | Maximum gap fraction (in %) tolerated when interpolating missing Hi-C bins. |
| `HIC_INSULATION_WINDOW` | int | `10` | Half-width (beads) of the sliding square used by Hi-C validation's insulation-score metric (validation.insulation_score) — the full window is 2x this value. Should roughly match the bead-scale size of a real TAD/domain for your resolution and N_BEADS; too small or too large relative to the actual domain size weakens the insulation_r correlation reported in hic_validation.npy even when the force itself is working well, since the metric is then measuring boundaries at the wrong genomic scale. |
| `HIC_BLOCK_COPOLYMER` | bool | `False` | Opt-in: derive A/B compartments straight from the Hi-C matrix's own PC1 (sign-aligned to local contact density: dense -> B, sparse -> A) and feed them into the same block-copolymer force normally built from a .bed file — the Boltzmann force alone reproduces diagonal decay/insulation well but tends to miss PC1 correlation, so this adds it back with no .bed file needed. Suggested use: enable it when you have HIC_USE_FORCE=True, no compartment .bed file, and a region large enough to actually contain compartments (see HIC_BLOCK_COPOLYMER_MIN_BP) — e.g. a whole chromosome or genome-wide run; leave it off for TAD/region-scale runs, where it has no effect anyway. Only takes effect when HIC_USE_FORCE=True and no COMPARTMENT_PATH is given; auto-disables (with a warning) when the modelled region is below HIC_BLOCK_COPOLYMER_MIN_BP. Raises an error if a .bed-based compartment force (COB_USE_COMPARTMENT_BLOCKS / SCB_USE_SUBCOMPARTMENT_BLOCKS) is also enabled — pick one source of compartments, not both. |
| `HIC_BLOCK_COPOLYMER_MIN_BP` | float | `5_000_000` | Minimum modelled region size (bp) for HIC_BLOCK_COPOLYMER to stay enabled. A/B compartments need several alternating domains to even be visible; below this a region is TAD-scale, not compartment-scale, so deriving compartments from it doesn't make sense. |

## Nuclear organisation

Container, lamina and nucleolus ([model](../model/lamina.md)).

| Argument | Type | Default | Description |
|---|---|---|---|
| `SC_USE_SPHERICAL_CONTAINER` | bool | `False` | Use Spherical container |
| `SC_SCALE` | float | `1000.0` | Spherical container scaling factor |
| `SC_RADIUS1` | quantity | *none* | Spherical container radius, **(unused)** |
| `SC_RADIUS2` | quantity | *none* | Spherical container radius, **(unused)** |
| `IBL_USE_B_LAMINA_INTERACTION` | bool | `False` | Interactions of B compartment with lamina. |
| `BLAMINA_FORCE_TYPE` | str | `sin` | B-lamina interaction functional form. Options: sin (default), gaussian_shell, harmonic_shell, logistic_shell. |
| `IBL_SCALE` | float | `400.0` | Scaling factor for B comoartment interaction with lamina. |
| `CF_USE_CENTRAL_FORCE` | bool | `False` | Attraction of smaller chromosomes. |
| `CENTRAL_FORCE_TYPE` | str | `harmonic` | Central nucleolar attraction functional form controlling radial bias toward nucleus center. Encodes chromosome-size dependent positioning. Options: harmonic (default, quadratic confinement around R1), gaussian (soft nucleolar enrichment field), logistic (soft-core radial partitioning with smooth boundary). |
| `CF_STRENGTH` | float | `20.0` | Strength of Interaction |

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
| `SIM_N_STEPS` | int | `10000` | Number of steps in MD simulation |
| `TRJ_FRAMES` | int | None | *none* | Number of CIF frames to save during MD. When set, SIM_SAMPLING_STEP is overridden to SIM_N_STEPS // TRJ_FRAMES so that exactly TRJ_FRAMES structures are written. When None, SIM_SAMPLING_STEP is used as-is. |
| `SIM_SAMPLING_STEP` | int | `100` | It determines in t how many steps we save a structure. |
| `SIM_INTEGRATOR_TYPE` | str | `langevin` | Alternative: langevin, verlet |
| `SIM_INTEGRATOR_STEP` | quantity | `1 femtosecond` | The step of integrator. |
| `SIM_FRICTION_COEFF` | float | `0.5` | Friction coefficient (Used only with langevin integrator) |
| `SIM_TEMPERATURE` | quantity | `310 kelvin` | Simulation temperature |
| `SIM_SET_INITIAL_VELOCITIES` | bool | `True` | Initialize the MD simulation with a random initial velocity field, drawn from the Maxwell-Boltzmann distribution at SIM_TEMPERATURE (OpenMM's setVelocitiesToTemperature), seeded by SHUFFLING_SEED. True by default so every run starts from a physically realistic, randomized velocity field instead of the all-zero velocities OpenMM uses otherwise. |
| `SIM_ERROR_TOLERANCE` | float | `0.01` | Error tolerance for variable MD simulation |
| `SIM_AMD_ALPHA` | float | `100.0` | Alpha of AMD simulation. |
| `SIM_AMD_E` | float | `1000.0` | E (energy) of AMD simulation. |
