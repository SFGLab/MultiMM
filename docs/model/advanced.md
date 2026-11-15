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

### Why an ensemble?

A single run gives you **one** structure. But a Hi-C map is not the picture of one nucleus: it is an *average* over millions of cells, and in each cell the chromatin is folded a little differently (loops form and dissolve, domains wobble). One model, however good, is just one possible member of that crowd. An **ensemble** is a set of structures, each one a plausible conformation, that together describe the variety. Use it when you care about *how much* the structure varies (for example the spread of the distance between two loci), or when you want to compare with a population average such as Hi-C, which is best reproduced by averaging many structures rather than looking at one.

!!! tip "In plain words"
    Think of taking many photographs of a crowd instead of one: each photo is different, but together they show what the crowd looks like. An ensemble is that stack of photos for a piece of chromatin.

### How to switch it on

```ini
GENERATE_ENSEMBLE = True
N_ENSEMBLE        = 20       ; number of structures (required)
```

`N_ENSEMBLE` must be set whenever `GENERATE_ENSEMBLE` is `True`; MultiMM stops with a clear message otherwise.

### How it works

MultiMM simply runs the whole pipeline $n$ times (`N_ENSEMBLE` = $n$), one after another. Run number $i$ ($i=0,\dots,n-1$) is written to its own folder `run_<i>/` inside `OUT_PATH` (numbers are zero-padded, so `run_00 … run_19`), with the usual layout of [outputs](../reference/output.md), and uses the random seed `SHUFFLING_SEED = i`. A **seed** is the starting number of a pseudo-random generator: the same seed always gives the same "random" choices, so every ensemble is **reproducible**, while different seeds give different structures.

The seed feeds every random ingredient of the model:

| Random ingredient | Switched on by | Effect |
|---|---|---|
| Random-walk initial structure | `INITIAL_STRUCTURE_TYPE` = `rw`, `confined_rw`, `self_avoiding_rw`, ... | A different starting conformation for each run |
| Initial velocities of the beads | `SIM_RUN_MD` with `SIM_SET_INITIAL_VELOCITIES` | A different dynamical trajectory for each run |
| Loops randomly removed | `DOWNSAMPLING_PROB` < 1 | Each run keeps a different random subset of loops, mimicking that not every loop is present in every cell |
| Compartment flips | `COMPARTMENT_FLIP_PROB` > 0 | Some beads switch A ↔ B at random |
| Compartment noise | `COMPARTMENT_NOISE_STD` > 0 | Compartment boundaries shift slightly |
| Hi-C noise | `HIC_NOISE_INTENSITY` > 0 | The strength of each Hi-C restraint is redrawn at every saved frame |
| Chromosome order | `SHUFFLE_CHROMS` | The order of chromosomes in a genome-wide model is permuted |

!!! warning "The runs must be able to differ"
    With the default settings the initial structure is the deterministic Hilbert curve and there is no noise, so a run that only minimises the energy gives the **same** structure every time and an ensemble would be $n$ identical copies. To get variety, turn on **at least one** ingredient from the table: a random-walk start, molecular dynamics with random initial velocities, or one of the noise options.

### Cost and tips

An ensemble multiplies the computing time by $n$, so it is best for small and medium systems (a gene or a region); for a whole genome, think twice. Start with a small $n$ (5 to 10) to check that the structures indeed differ, then increase it. Because the structures are saved as ordinary files, you can average their contact maps afterwards or measure any distance across the runs.

`SHUFFLE_CHROMS` randomly permutes the order of the chromosomes in a genome-wide model, which probes how much the arrangement of territories depends on the starting order.

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
