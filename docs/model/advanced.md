# Advanced options

This page collects options that change *what* is modelled or *how* the result is sampled, rather than the energy terms themselves.

## Modelling levels

`MODELLING_LEVEL` is a shortcut that sets many parameters at once for a typical scale. The preset is applied after the configuration is read and **overwrites** the listed parameters (the changes are printed at the start of the run). Everything else keeps your value.

| Level | Typical use | Beads | Container | Compartment force | B-lamina | MD |
|---|---|---|---|---|---|---|
| `GENE` | one gene and its surroundings | 1000 | off | off | off | on, 10 000 steps |
| `REGION` (or `LOC`) | a TAD-scale or Mb-scale region | 5000 | off | on if `COMPARTMENT_PATH` is given | off | on, 10 000 steps |
| `CHROM` (or `CHROMOSOME`) | a whole chromosome | 20 000 | off | on if `COMPARTMENT_PATH` is given | off | on, 10 000 steps |
| `GW` (or `GENOME`) | the whole genome | 200 000 | **on** | on if `COMPARTMENT_PATH` is given | on if `COMPARTMENT_PATH` is given | **off** |

All presets also switch off the chromosome blocks, the subcompartment force and the central force. `CHROM` additionally sets `LOC_START = 1` and `LOC_END` to the length of `CHROM`. Leave `MODELLING_LEVEL` empty to configure everything yourself.

!!! warning "Presets override your settings"
    A preset overwrites, for example, `N_BEADS`. If you want 20 000 beads on a region, do not use a level; set the parameters directly.

## Gene modelling

With `MODELLING_LEVEL = GENE`, the modelled region is built around a gene instead of from `CHROM`, `LOC_START` and `LOC_END`. Give the gene by name (`GENE_NAME`) or by identifier (`GENE_ID`). **Several genes** can be given as a comma-separated list (for example `GENE_NAME = MYC, PVT1`); the region then covers all of them plus the window (they must lie on the same chromosome). If a name has several annotations, the longest one is used. its position is looked up in an annotation table (`GENE_TSV`; the default is a bundled hg38 annotation). The region is the gene body extended by `GENE_WINDOW` base pairs on each side (default 100 000), cut to `N_BEADS` = 1000 beads.

**Every annotated gene in the modelled window** (also for `MODELLING_LEVEL = REGION`) is mapped onto the beads and reported, not only the requested ones:

| Output | Content |
|---|---|
| `metadata/genes_in_region.tsv` | one row per gene: coordinates, first and last bead, whether it is a target, its radius of gyration `rg` (how compact the gene is in 3D) and `dist_to_target` (distance of its centre to the target genes) |
| `plots/genes_3d.png` | the structure with every gene in its own colour (targets red, the rest of the chain grey) and gene-name labels |
| `plots/genes_track.png` | left: gene track along the beads, with overlapping genes stacked; right: matrix of 3D distances between gene centres, showing which genes are spatial neighbours |
| `metadata/chimera_genes_coloring.cmd` | ChimeraX script with the same colours |

The target gene alone is also highlighted in `plots/structure_afterMD_gene_coloring.png` and `metadata/chimera_gene_coloring.cmd`. Gene mode is meant to be used with local data: a loop file and/or a Hi-C map for the region.

```ini
[Main]
MODELLING_LEVEL = GENE
GENE_NAME       = MYC
GENE_WINDOW     = 200000
LOOPS_PATH      = loops.bedpe
OUT_PATH        = results/MYC
```

## Ensembles and randomness

One run gives one structure, but chromatin in a population of cells is heterogeneous. An **ensemble** is a set of structures that differ in the random ingredients of the model.

With `GENERATE_ENSEMBLE = True` and `N_ENSEMBLE = n`, MultiMM runs the whole pipeline $n$ times. Run $i$ is stored in `run_<i>/` inside `OUT_PATH` (numbers are zero-padded) and uses `SHUFFLING_SEED = i`. The seed controls every random choice: the random-walk initial structures, the loops removed by `DOWNSAMPLING_PROB`, the compartment flips and noise (`COMPARTMENT_FLIP_PROB`, `COMPARTMENT_NOISE_STD`) and the Hi-C noise (`HIC_NOISE_INTENSITY`). For the structures to differ, **at least one of these must be on**; a random-walk initial structure is the simplest source of variety. Ensembles multiply the cost by $n$, so they are best used for small systems.

`SHUFFLE_CHROMS` randomly permutes the order of the chromosomes in a genome-wide model, which probes how much the arrangement of territories depends on the starting order.

## Molecular dynamics

### Why minimisation, and why dynamics?

**Minimisation** answers the question: *what is the most comfortable shape for this chain?* Starting from the initial structure, OpenMM moves every bead a little in the direction that lowers the energy, again and again, until no small move helps (it uses an L-BFGS optimiser, a standard downhill method). The result is one low-energy structure, which is often all you need for a figure or for a large system.

A real chromosome is not frozen in one shape. It is constantly jiggled by thermal motion, and a population of cells contains many different shapes. **Molecular dynamics (MD)** mimics this: instead of walking only downhill, the beads move as if they were really at temperature $T$. The result is a *movie* of plausible shapes around the minimum, which lets you measure how much things fluctuate and gives many structures instead of one.

### How MD works: the integration scheme

MD follows Newton's law, "force equals mass times acceleration". The force on each bead is minus the slope of the energy, $\mathbf{F}_i=-\nabla_i E$. Computers cannot solve this continuously, so they advance time in tiny steps $\Delta t$ (the *time step*, here 1 femtosecond). A step has three moves, known as the **Verlet** scheme (shown here in its *velocity* form; OpenMM uses the equivalent leap-frog arrangement):

$$
\mathbf{v}\leftarrow\mathbf{v}+\tfrac{\Delta t}{2m}\mathbf{F}(\mathbf{r}),\qquad
\mathbf{r}\leftarrow\mathbf{r}+\Delta t\,\mathbf{v},\qquad
\mathbf{v}\leftarrow\mathbf{v}+\tfrac{\Delta t}{2m}\mathbf{F}(\mathbf{r}_\text{new}),
$$

that is: kick the velocities with the current force, let the beads drift, then kick again with the new force. Here $m$ is the bead mass. Plain Verlet conserves energy and so has no temperature control: the system would keep whatever temperature it started with.

To keep a fixed temperature, the default **Langevin** integrator adds two ingredients to each step: a *friction* that slows the beads down (coefficient $\gamma$, `SIM_FRICTION_COEFF`) and *random kicks* that stir them. Together they act as a thermostat, as if the chromatin sat in a thick fluid at temperature $T$:

$$
m\,\dot{\mathbf{v}}=-\nabla E-\gamma m\,\mathbf{v}+\sqrt{2\gamma m k_BT}\;\boldsymbol{\xi}(t),
$$

where $\boldsymbol{\xi}(t)$ is random noise. Friction removes energy, the noise adds it, and the balance gives exactly the temperature $T$. If the friction is large, the motion becomes pure diffusion, which is what the **Brownian** integrator does (it ignores velocities). Chromatin lives in a very viscous nucleus, so this is often a good description.

**Which one should I use?** Keep `langevin` (the default) unless you have a reason. Use `brownian` for slow, over-damped motion, `verlet` only for testing (it has no thermostat), and `amd` to explore more broadly by lowering the energy barriers.

By default MultiMM only **minimises** the energy (`SIM_RUN_MD = False`). The result is the nearest low-energy structure. With `SIM_RUN_MD = True` it continues with dynamics at temperature `SIM_TEMPERATURE` (310 K), which adds thermal fluctuations and allows many frames to be saved.

| Key | Default | Meaning |
|---|---|---|
| `SIM_INTEGRATOR_TYPE` | `langevin` | `langevin` (stochastic dynamics with friction, the default), `brownian`, `verlet` (no thermostat), `variable_verlet`, `variable_langevin` (adaptive time step), `amd` (accelerated MD) |
| `SIM_INTEGRATOR_STEP` | 1 femtosecond | time step |
| `SIM_FRICTION_COEFF` | 0.5 | friction in 1/ps (Langevin, Brownian) |
| `SIM_N_STEPS` | 10 000 | number of steps |
| `SIM_SAMPLING_STEP` | 100 | save a frame every this many steps |
| `TRJ_FRAMES` | none | if set, exactly this many frames are saved (it overrides `SIM_SAMPLING_STEP`) |
| `SIM_SET_INITIAL_VELOCITIES` | `True` | draw initial velocities from the temperature |
| `SIM_AMD_ALPHA`, `SIM_AMD_E` | 100.0, 1000.0 | AMD boost parameters |

Frames are written to `md_frames/`. The length scale of the model is not physical ([Units and scales](scales.md)), so the simulated time should be read as a relative, not an absolute, time.

!!! warning "Known issue: variable integrators"
    `variable_verlet` and `variable_langevin` read `SIM_ERROR_TOLERANCE` from the wrong object and currently fail with an error. Use `langevin` or `verlet`.

## Platforms

`PLATFORM` selects the OpenMM backend: `CPU` (default), `Reference` (slow, for testing), and, when available in your OpenMM installation, `CUDA` or `OpenCL` on a GPU. Large genome-wide systems need a GPU. `CPU_THREADS` sets the number of CPU threads. All custom nonbonded forces are given a cutoff so that they run on the GPU platforms.

## Nucleosome interpolation

The beads of MultiMM are coarse (kilobases). **Nucleosome interpolation** is an optional post-processing step that decorates the minimised structure with explicit nucleosomes, to produce an atomistic-looking fibre for visualisation. It does not change the modelled structure.

**Number of nucleosomes.** From an ATAC-seq signal track (`ATACSEQ_PATH`, a `.bigWig`), the log signal $\log(a_i+10^{-6})$ is rescaled to $\tilde a_i\in[0,1]$. The number of nucleosomes placed between bead $i$ and $i+1$ is proportional to this signal (a higher value gives more nucleosomes):

$$
n_i=\operatorname{round}\!\bigl(\tilde a_i\,n_\text{max}\bigr),\qquad n_\text{max}=\left\lceil \frac{1}{2\,r_\text{nuc}}\right\rceil ,
$$

where $r_\text{nuc}$ is the nucleosome radius (`NUC_RADIUS`, 0.1) and $n_\text{max}$ is the maximum number of nucleosomes per bead (5 for the default).

**Shape.** Each nucleosome is a short **helix** of radius $r_\text{nuc}$, height $r_\text{nuc}$ and **1.65 turns** (the DNA wraps about 1.65 times around the histone core), drawn with `POINTS_PER_NUC` points (20). Successive nucleosomes are placed at equal spacing along the segment between two beads, offset sideways in a **zigzag** by half a linker length ($3.45\,r_\text{nuc}/2$, randomly scaled between 0.5 and 1.5), and rotated by a constant angle `PHI_NORM` ($\pi/5$) between neighbours, which mimics the zigzag chromatin fibre.

!!! warning "Known issue: not functional in this release"
    The step cannot currently run. The modelling-level code always resets `NUC_DO_INTERPOLATION` to `False` and `ATACSEQ_PATH` to none, and, in addition, the call that builds the nucleosomes passes its arguments in a different order from the function signature (`MAX_NUCS_PER_BEAD` is also ignored: $n_\text{max}$ is always computed from `NUC_RADIUS`). The output file `model/MultiMM_minimized_with_nucs.cif` is therefore not produced. The description above is how the method is designed to work.

!!! info "In MultiMM"
    | Key | Default |
    |---|---|
    | `NUC_DO_INTERPOLATION` | `False` |
    | `ATACSEQ_PATH` | none |
    | `NUC_RADIUS` | 0.1 |
    | `POINTS_PER_NUC` | 20 |
    | `PHI_NORM` | 0.6283 ($\pi/5$) |
    | `MAX_NUCS_PER_BEAD` | 4 (ignored) |
