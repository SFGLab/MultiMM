# MultiMM — Example configuration files

This folder contains ready-to-use configuration files covering the most common
modelling scenarios.  Copy the relevant file, update the file paths to point to
your data, and run:

```bash
MultiMM config.ini          # or point to any .ini file here
```

---

## Files at a glance

| File | Use-case | Key data needed |
|---|---|---|
| `config_gw.ini` | Full genome (all chromosomes) | loops + subcompartments |
| `config_specific_region.ini` | Single TAD / gene locus | loops (+ optional compartments) |
| `config_single_cell.ini` | Genome-wide single-cell | single-cell .bedpe loops |
| `config_hic.ini` | Hi-C contact-force (region) | .hic or .cool matrix |

---

## config_gw.ini — Genome-wide simulation

Models all chromosomes as one polymer chain inside a spherical nucleus.
All major force terms are available:

- **Loop extrusion** (`LE_USE_HARMONIC_BOND`) — enforces CTCF/cohesion loops
- **Subcompartment blocks** (`SCB_USE_SUBCOMPARTMENT_BLOCKS`) — A1/A2/B1/B2 phase separation
- **B-lamina interaction** (`IBL_USE_B_LAMINA_INTERACTION`) — anchors B-compartment to the nuclear periphery
- **Chromosomal blocks** (`CHB_USE_CHROMOSOMAL_BLOCKS`) — self-compaction of each chromosome into a globule
- **Central force** (`CF_USE_CENTRAL_FORCE`) — size-dependent nucleolar positioning
- **Spherical container** (`SC_USE_SPHERICAL_CONTAINER`) — hard nuclear boundary

Required: `LOOPS_PATH` (.bedpe) and `COMPARTMENT_PATH` (.bed, Calder format).

---

## config_specific_region.ini — TAD / sub-chromosomal region

Zooms into a defined genomic interval (`CHROM`, `LOC_START`, `LOC_END`).
Loop extrusion shapes local TAD and loop structure; compartment blocks can be
added if subcompartment data covering the region is available.
Nucleosome interpolation (`NUC_DO_INTERPOLATION`) adds fine-scale chromatin
fibre detail using an ATAC-seq signal.

Required: `LOOPS_PATH` (.bedpe) filtered to the target region.

---

## config_single_cell.ini — Single-cell genome-wide

Designed for sparse single-cell contact data (Dip-C, scHi-C, scChIA-PET).
Because signal is binary (contact / no contact), loop bond distances are fixed
(`LE_FIXED_DISTANCES = True`) rather than scaled by contact frequency.
Compartment forces are disabled (no subcompartment labels available per cell);
nuclear confinement and central force replace them.

Required: `LOOPS_PATH` (.bedpe) with single-cell contact pairs.

---

## config_hic.ini — Hi-C contact-force simulation

Uses a raw Hi-C contact matrix (`.hic` or `.cool`/`.mcool`) as a direct
structural restraint via a low-rank force decomposition.  The pipeline:

1. `read_hic.py` loads and normalises the matrix, auto-selects resolution,
   and resizes it to exactly `N_BEADS × N_BEADS` via weighted average pooling.
2. `hic_force.py` decomposes the O/E-normalised matrix by SVD and builds an
   OpenMM `CustomNonbondedForce` that attracts beads with correlated
   eigenvector components.

Three force modes are available:

| `HIC_FORCE_MODE` | Description |
|---|---|
| `svd` | Single-σ CustomNonbondedForce.  Fast, good default. |
| `svd_multiscale` | Per-component σ scaled by eigenvalue.  More accurate. |
| `crossentropy` | Sparse CustomBondForce.  Best for sharp loop peaks. |

Required: `HIC_PATH` pointing to a `.hic`, `.cool`, or `.mcool` file.

---

## All configuration parameters

A full list of available parameters with defaults and descriptions is printed by:

```bash
python -c "from multimm.config import SimulationConfig; help(SimulationConfig)"
```

or can be found in `src/multimm/config.py`.
