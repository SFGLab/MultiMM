![MultiMM](https://github.com/user-attachments/assets/3996e49c-7db4-4194-8af9-493bb82ab87a)

# MultiMM: An OpenMM-based software for whole-genome 3D structure reconstruction

MultiMM is an OpenMM model designed for modeling the 3D structure of the whole *human* genome. Its distinguishing feature is that it is multiscale, meaning it aims to model different levels of chromatin organization, from smaller scales (nucleosomes) to the level of chromosomal territories. The algorithm is both fast and accurate. A key feature enabling its speed is GPU parallelization via OpenMM, along with smart assumptions that assist the optimizer in finding the global minimum. One such fundamental assumption is the use of a Hilbert curve as the initial structure. This helps MultiMM converge faster because the initial structure is already highly compacted.

![GW_em](https://github.com/user-attachments/assets/3d019616-2c0f-4bfc-a792-d87fcdc0d96a)

After running MultiMM, users obtain a genome-wide 3D structure. Chromosomes, compartments, or individual genomic regions can be colored and visualized separately. MultiMM is simple to use: all parameters are controlled through a single configuration file.

![MultiMM_scales](https://github.com/user-attachments/assets/a0d14ddc-41bf-4d14-8a8c-aa418fe575b5)

The workflow is illustrated in the schematic above. The user provides chromatin loop calls from a 3C-type experiment and, optionally, compartment annotations. MultiMM imports an initial structure, preprocesses the input, applies a physically motivated force field, and produces a 3D structure. When ATAC-seq data are supplied, nucleosome interpolation is applied as a post-processing step.

![image](https://github.com/user-attachments/assets/4a446111-241d-427d-b568-c03e7a2c63c4)

---

## Key Features

- OpenMM-based simulation engine with GPU acceleration (CUDA / OpenCL) and CPU fallback.
- User-friendly installation via PyPI; all parameters set in a single `config.ini` file.
- Multiscale: nucleosome → TAD → compartment → chromosome territory → whole nucleus.
- Hi-C contact-guided force field: raw `.hic` / `.cool` / `.mcool` matrices used directly as structural restraints via a classic **Boltzmann-inversion** CustomBondForce — each pair's contact strength is converted into a target 3-D distance via the Hi-C scaling law and restrained there with a harmonic well weighted by its own observed contact strength, so weakly-supported pairs stay soft rather than acting as hard constraints.
- Ensemble generation: multiple independent structures from a single run.
- Nucleosome interpolation from ATAC-seq signal.
- **Comprehensive Hi-C validation** computed automatically: diagonal decay correlation, insulation score correlation, PC1 compartment correlation (sign-aligned to contact density), Pearson / Spearman OE-matrix correlation, SSIM, GMSD, NMI — each reported alongside a random-walk null-model baseline for direct comparison.
- **Post-simulation quality control suite** (10 checks): energy stability, bond distances, angle distribution, excluded-volume overlaps, compartment clustering, chromosome separation, loop-distance compliance, container confinement, B-lamina proximity, and MD structural mobility (RMSD vs. minimised structure).
- Structured logger with coloured, time-stamped output, section banners, and per-stage success messages.

---

## About Operating Systems

MultiMM has been tested primarily on Linux (Ubuntu, Debian, Red Hat). It can also run on macOS, though without CUDA support. Running on Windows is not recommended.

---

## Installation

```bash
pip install MultiMM
```

PyPI package: https://pypi.org/project/MultiMM/

---

## Model Overview

Chromatin is represented as a coarse-grained polymer. The total energy `E` decomposes into physically motivated terms:

```
E = E_backbone + E_loops + E_block + E_excluded + E_confinement + E_chromosomal + E_HiC
```

Each term encodes a distinct biological mechanism:

| Term | Mechanism |
|---|---|
| `E_backbone` | Polymer connectivity and stiffness |
| `E_loops` | Long-range loop extrusion or experimental contact constraints |
| `E_block` | Compartment and subcompartment phase separation |
| `E_excluded` | Steric repulsion between beads |
| `E_confinement` | Nuclear geometry: spherical container and lamina affinity |
| `E_chromosomal` | Chromosome territory formation and global compaction |
| `E_HiC` | Data-driven contact restraint from a raw Hi-C matrix |

---

### Polymer Backbone

The backbone encodes chain connectivity and local rigidity through two standard terms.

**Harmonic bond** (nearest-neighbor connectivity):

```
E_bond = sum_i  (k_b / 2) * (r_{i,i+1} - r0)^2
```

where `r_{i,i+1}` is the distance between consecutive beads, `r0` is the equilibrium bond length, and `k_b` is the bond stiffness.

**Harmonic angle** (chain stiffness / persistence length):

```
E_angle = sum_i  (k_theta / 2) * (theta_i - theta0)^2
```

where `theta_i` is the angle formed by three consecutive beads `(i, i+1, i+2)`, `theta0` is the preferred angle, and `k_theta` controls bending rigidity. Together these reproduce a discretized worm-like chain.

---

### Loop Interactions

Long-range loop constraints tether pairs of beads `(m, n)` identified from loop-calling experiments. Three functional forms are available.

**Harmonic** (default):

```
E_loops_harmonic = sum_{(m,n)}  (k / 2) * (r_mn - r0)^2
```

**Soft FENE-like** (bounded, avoids divergence at large extension):

```
E_loops_fene = sum_{(m,n)}  k * (r_mn - r0)^2 / (1 + alpha * (r_mn - r0)^2)
```

**Gaussian tether** (smooth, fully bounded):

```
E_loops_gaussian = sum_{(m,n)}  k * (1 - exp(-(r_mn - r0)^2 / sigma^2))
```

Loop bond strengths can be fixed (`LE_FIXED_DISTANCES = True`) or scaled by experimental contact frequency (`LE_FIXED_DISTANCES = False`). Providing `LOOPS_PATH` is optional; the simulation runs without loops if it is not supplied.

---

### Block-Copolymer Compartmentalization

State-dependent pairwise attractions drive compartment phase separation. Each bead carries a label `s_i` representing compartment identity.

**Compartment level (A/B):**

```
E_comp = -sum_{i<j}  eps(s_i, s_j) * exp(-r_ij^2 / (2 * r_c^2))
```

The coupling `eps(s_i, s_j)` is attractive for like compartments (A–A, B–B) and weak or repulsive otherwise, reproducing large-scale A/B segregation.

**Subcompartment level (A1/A2/B1/B2):**

```
E_sub = -sum_{i<j}  eps_ab * exp(-r_ij^2 / (2 * r_sc^2)),   s_i = alpha, s_j = beta
```

This promotes finer microphase separation inside A/B compartments.

**Chromosome territory (self-compaction):**

```
E_chrom = sum_{i<j}  delta(chi_i, chi_j) * V(r_ij)
```

where `chi_i` is the chromosome label and `V(r)` is a soft attractive potential acting only between beads on the same chromosome. The default polynomial form is:

```
V(r) = dE * (k_C * r^4 - r^3 + r^2)
```

**Alternative interaction kernels** (experimental):

| Kernel | Expression | Effect |
|---|---|---|
| Yukawa | `V(r) ~ -exp(-r/lambda) / r` | Screened, longer-range |
| Power-law | `V(r) ~ -1 / (r^alpha + eps)` | Scale-free attraction |
| Theta (contact) | `V(r) ~ -Theta(r_c - r)` | Binary hard-cutoff |
| Saturating | `V(r) ~ -1 / (1 + k_C * r^2)` | Bounded, prevents over-collapse |

---

### Nuclear Geometry and Lamina Interactions

**Spherical container** — soft penalty for excursions outside the nuclear shell:

```
E_container = C * sum_i [ max(0, r_i - R2)^2 + max(0, R1 - r_i)^2 ]
```

where `r_i` is the radial distance from the nuclear center. This confines chromatin between radii `R1` and `R2`.

**B-lamina interaction** — anchors B-compartment chromatin to the nuclear periphery:

```
E_lamina = -sum_i  B(s_i) * V(r_i)
```

where `B(s_i)` selects B-compartment beads. Available radial profiles:

| Mode | Expression | Description |
|---|---|---|
| `sin` (default) | `V(r) = sin^8(pi*(r-R1)/(R2-R1)) - 1` | Sharp peripheral preference |
| `gaussian_shell` | `V(r) ~ -(exp(-(r-R1)^2/(2*s^2)) + exp(-(r-R2)^2/(2*s^2)))` | Localized at both boundaries |
| `harmonic_shell` | `V(r) ~ (r - r0)^2,  r0 = (R1+R2)/2` | Pulls toward mid-shell |
| `logistic_shell` | Smooth sigmoidal walls at `R1` and `R2` | Smooth boundary transition |

---

### Hi-C Contact-Guided Force

When a raw Hi-C contact matrix is provided (`HIC_PATH`), it's used directly as a structural restraint via a sparse Boltzmann-inversion `CustomBondForce` — an alternative or complement to loop-extrusion forces, with no explicit loop calls required.

**Pipeline:** the matrix is loaded (`.hic` / `.cool` / `.mcool`) and resampled to `N_beads x N_beads` (`read_hic.py`), then Knight–Ruiz balanced to a normalised `c_ij ∈ [0, 1]`. With `HIC_FORCE_OE=True`, an Observed/Expected step (divide each diagonal by its mean, subtract 1, floor at 0) is applied on top, so `c_ij` targets relative enrichment over the distance-decay background instead of absolute contact frequency. A bond is built for every pair with nonzero `c_ij` — no separate sparsity cutoff.

> **Why `HIC_FORCE_OE=False` is the default:** OE scores higher on PC1/OE-Pearson metrics in isolated testing, but raw (KR-balanced) contact frequency gave better overall results in practice — stronger diagonal decay and insulation score — so it's the default. Set `HIC_FORCE_OE=True` if compartment/TAD-level enrichment metrics matter more to you than the bulk distance-decay trend.

The restraint is a classic **Boltzmann inversion**, `U(r) = -k_B T ln P(r)`, assuming the equilibrium distance distribution for a pair is Gaussian around a data-derived target. Each `c_ij` converts to a target distance via one of three interchangeable `P(r)` shapes, selected by `HIC_BOLTZMANN_KERNEL`:

| Kernel | `P(r)` |
|---|---|
| `exponential` (default) | `clip(exp(-(r-r_min)/lambda), 0, 1)`, `lambda = rc/alpha` — the classic Boltzmann distribution |
| `power_law` | `clip((r_min/r)^alpha, 0, 1)` — the standard Hi-C scaling law, `c ∝ r^(-α)` |
| `sigmoid` | `1 / (1 + exp(k*(r-r0)))`, `r0 = (r_min+rc)/2`, `k = alpha/rc` — bounded, logistic contact probability |

All three are exact functional inverses of their own `P(r)`, share `HIC_BOLTZMANN_ALPHA` as a steepness knob, and are clipped to `[r_min, rc]` before being restrained with the same **two-sided** harmonic well — `U_ij(r) = 1/2 * k_scale * c_ij * (r - r_target,ij)^2` — which pulls a pair farther than its target *and* pushes one that has overshot closer. `r_min` is the excluded-volume floor distance; `rc` is `HIC_RC` when set explicitly, otherwise auto-calibrated from the initial structure's own pairwise distances (`HIC_AUTO_SCALE=True`, the default).

Recommended `HIC_K_SCALE`: 20–80 kJ mol⁻¹ (default 40). Higher values (up to ~160) consistently improved insulation-score validation in testing with no diagonal-decay cost, for `N_BEADS` ≳ 300 — below that, gains were less consistent, so lower values may suit small systems better. `HIC_BOLTZMANN_ALPHA` (typical range 3–4) sets how steep the strength→distance mapping is for every kernel.

> **Spherical-shell artifact & `HIC_BOLTZMANN_TOL_FRAC`:** with `HIC_FORCE_OE=False`, real Hi-C's long power-law decay tail means almost every pair gets some nonzero `c_ij`, and the weakest/longest-range ones often end up with nearly the same `r_target` (frequently capped at `rc`). Restraining a large share of all pairs to one *exact* shared distance is only satisfiable in 3-D by spreading beads over a spherical shell (the same mechanism behind the Thomson problem), producing an unnaturally round, hollow-looking structure instead of a graded globule. `HIC_BOLTZMANN_TOL_FRAC` (float, default `0.0`) widens each pair's well into a flat-bottom band — zero force within `±tol_frac * r_target`, harmonic beyond it — so weak-evidence pairs get genuine slack instead of false precision. Try `0.15–0.3` if your structures look shell-like; `0` reproduces the original exact harmonic well.

> **Compartments (PC1):** this per-pair force optimises individual target distances, not the whole-matrix PC1 signal — on its own it reproduces diagonal decay and insulation well but tends to under-shoot PC1 correlation. `HIC_BLOCK_COPOLYMER` (below, opt-in) fixes this; the older alternative is the dedicated bed-based compartment force (`COB_USE_COMPARTMENT_BLOCKS`, driven by `COMPARTMENT_PATH`).

#### Compartments from the Hi-C matrix itself (`HIC_BLOCK_COPOLYMER`)

`HIC_BLOCK_COPOLYMER` (default `False`, opt-in) derives A/B compartments straight from the Hi-C matrix instead of requiring a separate `.bed` file: PC1 of the O/E matrix is computed, sign-aligned to each bead's own local contact density (dense → B, sparse → A — the real Hi-C convention), discretized to the same `+1`/`-1` labels `import_bed()` builds from `.bed` files, and fed into the *same* block-copolymer force (`add_compartment_blocks`, `COB_EA`/`COB_EB`) normally reserved for bed-based compartments.

**When to use it:** turn it on for a whole-chromosome or genome-wide run with `HIC_USE_FORCE=True` and no compartment `.bed` file on hand — it's the easiest way to get compartment-level structure without one. Leave it off for a TAD/region-scale run (it has no effect there anyway, see the size gate below) or whenever you already have a `.bed` compartment track, which is usually the more reliable source.

It only activates when `HIC_USE_FORCE=True` and no `COMPARTMENT_PATH` is given, and auto-disables itself (with a warning) when the modelled region is below `HIC_BLOCK_COPOLYMER_MIN_BP` (default 5 Mb) — a region that small is TAD-scale, not compartment-scale, so there's nothing for PC1 to resolve. Two guardrails: a warning if `HIC_FORCE_OE=True` is also set (the Boltzmann force already over-weights compartments under OE, so stacking the block-copolymer force on top risks double-counting them — `HIC_FORCE_OE=False` is recommended alongside `HIC_BLOCK_COPOLYMER`); and an error if a bed-based compartment force (`COB_USE_COMPARTMENT_BLOCKS` / `SCB_USE_SUBCOMPARTMENT_BLOCKS`) is enabled at the same time — pick one source of compartments, not both.

**Softer by default than the `.bed` path:** Hi-C-derived PC1 labels are coarser and noisier than a curated `.bed` compartment call, so `COB_EA`/`COB_EB` (tuned for `.bed` data) over-aggregate A/B segregation if applied at full strength here. `HIC_BLOCK_COPOLYMER_STRENGTH_SCALE` (default `0.4`) scales `COB_EA`/`COB_EB` down for this path only — a `.bed`-based compartment force always uses `COB_EA`/`COB_EB` unscaled. Set it to `1.0` for full strength, or lower still if structures still look over-aggregated.

**Validation** — after simulation, MultiMM reports seven metrics against the experimental Hi-C matrix, each alongside a random-walk null baseline, using the same `HIC_BOLTZMANN_KERNEL` `P(r)` the force itself was built with (`hic_force.get_boltzmann_p_func`):

| Metric | What it measures |
|---|---|
| Diagonal decay correlation | Whether contact frequency falls off with genomic distance at the right rate |
| Insulation score correlation | Agreement of TAD boundary positions (smoothed, normalized to [0, 1] before scoring) |
| PC1 compartment correlation | A/B compartment identity — PC1 sign-aligned to contact density (same convention as `HIC_BLOCK_COPOLYMER`) before scoring, smoothed and normalized to [-1, 1] |
| Pearson / Spearman OE correlation | Global agreement of the observed-over-expected matrices |
| SSIM | Structural similarity (local contrast, luminance, structure) |
| GMSD | Edge/boundary sharpness agreement |
| NMI | Normalised mutual information |

Results are saved to `metadata/hic_validation.npy` (keys: `diagonal_decay_r`, `insulation_r`, `pc1_r`, `pearson_oe_r`, `spearman_oe_r`, `ssim`, `gmsd`, `nmi`, plus `_rw_*` null-baseline variants and `_p` p-values where applicable).

---

### Nucleosome-Scale Interpolation

After coarse-grained optimization, nucleosome positions are interpolated using a beads-on-a-string zigzag model. Each nucleosome is represented as a helix with 1.65 DNA turns. The number of nucleosomes per bead is derived from normalized ATAC-seq signal, enforcing nucleosome-rich regions in low-accessibility chromatin.

---

## Internal Parameter Definitions

All geometric and interaction scales are derived from a single microscopic length scale — the polymer bond length `b0` (`POL_HARMONIC_BOND_R0`).

**Nuclear radius** (dense globule scaling):

```
R2 = b0 * N^(1/3)
```

which enforces constant monomer density `N / R2^3 ≈ const`.

**Inner compartment radius** (fixed volume fraction `f`):

```
R1 = R2 * f^(1/3)
```

**Compartment interaction length scale:**

```
r_c ~ O(b0) ≈ 1.5 * b0
```

ensuring interactions remain local relative to the polymer backbone.

**Loop equilibrium distances** — either globally fixed or derived from experimental loop lengths `d_i`:

```
r0_i ∈ { r0_global,  d_i }
```

**State variables** — each bead carries a compartment label `s_i` and chromosome label `chi_i`, which modulate interactions through selection rules:

```
E_ij ∝ delta(s_i, s_j)       (compartment-selective)
E_ij ∝ delta(chi_i, chi_j)   (chromosome-selective)
```

---

## Input Data

MultiMM accepts four types of input files:

| File | Format | Parameter | Required |
|---|---|---|---|
| Chromatin loops | `.bedpe` | `LOOPS_PATH` | Optional |
| Compartment labels | `.bed` (CALDER format) | `COMPARTMENT_PATH` | Optional |
| ATAC-seq signal | `.bw` / `.BigWig` | `ATACSEQ_PATH` | Optional |
| Hi-C contact matrix | `.hic` / `.cool` / `.mcool` | `HIC_PATH` | Required if `HIC_USE_FORCE=True` |

### Loops (`.bedpe`)

Seven-column file, no header, one interaction per row:

```
chr10  100225000  100230000  chr10  100420000  100425000  95
chr10  100225000  100230000  chr10  101005000  101010000  56
chr10  101190000  101195000  chr10  101370000  101375000  152
```

Columns 1–3: first anchor (chrom, start, end); columns 4–6: second anchor; column 7: contact strength. The file may contain all chromosomes; MultiMM selects the relevant subset automatically.

For single-cell data, set columns 2 and 3 to the same value (and columns 5 and 6 likewise) and use strength = 1.

### Compartments (`.bed`, CALDER format)

Produced by [CALDER2](https://github.com/CSOgroup/CALDER2). The file must contain at least four columns: chrom, start, end, label.

```
chr1  700001   900000   A.1.2.2.2.2.2.2
chr1  900001  1400000   A.1.1.1.1.2.1.1.1.1.1
chr1 1400001  1850000   A.1.1.1.1.2.1.2.2.2.1
chr1 1850001  2100000   B.1.1.2.2.1.2.1
```

### Hi-C contact matrix (`.hic` / `.cool` / `.mcool`)

Any standard Juicer `.hic` file or Cooler `.cool` / `.mcool` file is accepted. MultiMM automatically selects the best available resolution for the requested region and resamples to `N_beads x N_beads`. KR, VC, VC_SQRT, and NONE normalizations are supported.

### ATAC-seq signal (`.bw` / `.BigWig`)

A p-value BigWig track. Required only for nucleosome interpolation (`NUC_DO_INTERPOLATION = True`). The pyBigWig library is required; note that pyBigWig is not compatible with Windows.

### Gene-based region definition (optional)

To model a genomic window around a specific gene, provide a `.tsv` with gene annotations and set `GENE_NAME` or `GENE_ID`:

```
gene_id           gene_name  chromosome  start     end
ENSG00000160072   ATAD3B     chr1        1471765   1497848
ENSG00000142611   PRDM16     chr1        3069168   3438621
```

When a gene is specified, MultiMM generates visualizations with the gene highlighted in red.

![minimized_structure_gene_coloring](https://github.com/user-attachments/assets/15666286-0162-4fdd-a875-6b50b65049fb)

> **Note:** MultiMM is designed for human genome data. Other organisms may work with additional modifications, but full support is not guaranteed. MultiMM can process data from Hi-C, scHi-C, ChIA-PET, and Hi-ChIP experiments. Default parameters are optimized for population-average Hi-C data; users are encouraged to validate convergence on their own datasets. Read the method paper before adjusting force parameters.

---

## Usage

All parameters are specified in a `config.ini` file. The minimal genome-wide example:

```ini
[Main]

PLATFORM = OpenCL

; Input data
FORCEFIELD_PATH  = forcefields/ff.xml
LOOPS_PATH       = /path/to/loops.bedpe
COMPARTMENT_PATH = /path/to/calder_subcompartments.bed
ATACSEQ_PATH     = /path/to/atac.bw
OUT_PATH         = results

; Bead resolution
N_BEADS          = 50000
SHUFFLE_CHROMS   = True

; Force field — genome-wide
SC_USE_SPHERICAL_CONTAINER    = True
CHB_USE_CHROMOSOMAL_BLOCKS    = True
SCB_USE_SUBCOMPARTMENT_BLOCKS = True
IBL_USE_B_LAMINA_INTERACTION  = True
CF_USE_CENTRAL_FORCE          = True
NUC_DO_INTERPOLATION          = True

; MD annealing (optional)
SIM_RUN_MD        = True
SIM_N_STEPS       = 1000
SIM_SAMPLING_STEP = 50
TRJ_FRAMES        = 100
```

Run with:

```bash
MultiMM -c config.ini
```

Example data (GM12878, Rao et al.; CALDER subcompartments; ENCODE ATAC-seq) is available at:  
https://drive.google.com/drive/folders/1nFAPE4pCaHpeL5nw6nq0VvfUFoc24aXm?usp=sharing

Ready-to-use configuration files for common scenarios are in the `examples/` folder, kept in sync with the current `SimulationConfig` field set (see `src/multimm/config.py`).

> **Note:** every key in a `config.ini` file must match a field defined in `SimulationConfig` — a typo, a removed/renamed field (e.g. the old `HIC_ALPHA` / `HIC_THRESHOLD`), or anything else not recognised raises a clear `ValueError` at startup naming the offending file and every bad key, with a "did you mean ...?" suggestion (closest real field name) where one exists, instead of being silently ignored. A missing or mistyped `-c`/`--config_file` path is also caught explicitly (it used to fail silently and just run with class defaults), and a malformed INI file surfaces `configparser`'s own file/line error.

---

## Modelling Levels

The `MODELLING_LEVEL` parameter automatically configures resolution and forces for common use cases:

| Level | Scope | Default `N_BEADS` | Forces active |
|---|---|---|---|
| `GENE` | ±100 kb window around a gene | 1 000 | Backbone + loops |
| `REGION` | User-defined chromosomal interval | 5 000 | Backbone + loops (+ optional compartments) |
| `CHROM` | Full chromosome | 20 000 | Backbone + loops + compartments |
| `GW` | All chromosomes | 200 000 | Full force field |

When `MODELLING_LEVEL` is set, the specified `N_BEADS` value is overridden. Advanced users who need fine-grained control should leave this parameter unset.

---

## Visualization

**Genome-wide structure with chromosome coloring:**

```python
import multimm.plots as splt
splt.viz_chroms(sim_path)          # add comps=False to disable compartment coloring
```

**Any single region or CIF structure:**

```python
import multimm.plots as splt
import multimm.utils as suts

V = suts.get_coordinates_cif(cif_path)
splt.viz_structure(V)
```

Visualization is powered by [PyVista](https://pyvista.org/).

---

## Simulation Arguments

### Platform and Device

| Parameter | Type | Default | Description |
|---|---|---|---|
| `PLATFORM` | str | `CPU` | Compute platform: `CPU`, `OpenCL`, `CUDA`, `Reference` |
| `CPU_THREADS` | int | None | Number of CPU threads (CPU platform only) |
| `DEVICE` | str | `""` | Device index for CUDA / OpenCL (count from 0) |

### Input and Output

| Parameter | Type | Default | Description |
|---|---|---|---|
| `FORCEFIELD_PATH` | str | bundled | Path to OpenMM XML forcefield |
| `LOOPS_PATH` | str | None | `.bedpe` loop file (optional) |
| `COMPARTMENT_PATH` | str | None | `.bed` subcompartment file (CALDER format) |
| `ATACSEQ_PATH` | str | None | `.bw` / `.BigWig` ATAC-seq p-value track |
| `HIC_PATH` | str | None | `.hic` / `.cool` / `.mcool` contact matrix |
| `OUT_PATH` | str | `results` | Output directory |
| `INITIAL_STRUCTURE_PATH` | str | `""` | Path to existing `.cif` initial structure |
| `GENE_TSV` | str | bundled | Gene annotation `.tsv` |
| `GENE_NAME` | str | `""` | Gene name for region definition |
| `GENE_ID` | str | `""` | Ensembl gene ID for region definition |

### Region of Interest

| Parameter | Type | Default | Description |
|---|---|---|---|
| `CHROM` | str | None | Chromosome (e.g. `chr1`); leave blank for genome-wide |
| `LOC_START` | int | None | Start coordinate (bp) |
| `LOC_END` | int | None | End coordinate (bp) |
| `GENE_WINDOW` | int | 100 000 | Flanking window around gene (bp) |
| `MODELLING_LEVEL` | str | `""` | Auto-configure: `GENE`, `REGION`, `CHROM`, `GW` |

### Initial Structure

| Parameter | Type | Default | Description |
|---|---|---|---|
| `BUILD_INITIAL_STRUCTURE` | bool | `True` | Build a new initial structure |
| `INITIAL_STRUCTURE_TYPE` | str | `hilbert` | `hilbert`, `circle`, `rw`, `confined_rw`, `self_avoiding_rw`, `helix`, `spiral`, `sphere`, `knot` |

### Simulation

| Parameter | Type | Default | Units | Description |
|---|---|---|---|---|
| `N_BEADS` | int | 50 000 | — | Number of coarse-grained beads |
| `SHUFFLE_CHROMS` | bool | `False` | — | Randomize chromosome order |
| `SHUFFLING_SEED` | int | 0 | — | Random seed for chromosome shuffling |
| `SIM_RUN_MD` | bool | `False` | — | Run MD annealing after energy minimization |
| `SIM_N_STEPS` | int | 10 000 | — | Number of MD steps |
| `SIM_SAMPLING_STEP` | int | 100 | — | Steps between saved trajectory frames |
| `SIM_TEMPERATURE` | Quantity | 310 | K | Simulation temperature |
| `SIM_SET_INITIAL_VELOCITIES` | bool | `True` | — | Randomize the initial velocity field (Maxwell-Boltzmann at `SIM_TEMPERATURE`, seeded by `SHUFFLING_SEED`) instead of starting from all-zero velocities |
| `SIM_INTEGRATOR_TYPE` | str | `langevin` | — | `langevin`, `verlet`, `brownian` |
| `SIM_INTEGRATOR_STEP` | Quantity | 1 | fs | Integrator time step |
| `SIM_FRICTION_COEFF` | float | 0.5 | ps⁻¹ | Friction coefficient (Langevin / Brownian) |
| `TRJ_FRAMES` | int \| None | `None` | — | When set, overrides `SIM_SAMPLING_STEP` to `SIM_N_STEPS // TRJ_FRAMES` so exactly this many CIF frames are saved; `None` uses `SIM_SAMPLING_STEP` as-is |

### Polymer Backbone

| Parameter | Type | Default | Units | Description |
|---|---|---|---|---|
| `POL_USE_HARMONIC_BOND` | bool | `True` | — | Harmonic bond between consecutive beads |
| `POL_HARMONIC_BOND_R0` | Quantity | 0.1 | nm | Equilibrium bond length |
| `POL_HARMONIC_BOND_K` | Quantity | 300 000 | kJ mol⁻¹ nm⁻² | Bond stiffness |
| `POL_USE_HARMONIC_ANGLE` | bool | `True` | — | Harmonic angle (bending rigidity) |
| `POL_HARMONIC_ANGLE_R0` | Quantity | pi | rad | Equilibrium angle |
| `POL_HARMONIC_ANGLE_CONSTANT_K` | Quantity | 100 | kJ mol⁻¹ rad⁻² | Bending stiffness |

### Excluded Volume

| Parameter | Type | Default | Units | Description |
|---|---|---|---|---|
| `EV_USE_EXCLUDED_VOLUME` | bool | `True` | — | Steric repulsion |
| `EV_FORCE_TYPE` | str | `powerlaw` | — | `powerlaw`, `gaussian_core` |
| `EV_EPSILON` | float | 100.0 | kJ mol⁻¹ | Repulsion strength |
| `EV_R_SMALL` | float | 0.05 | nm | Regularization radius |
| `EV_POWER` | float | 6.0 | — | Exponent of power-law repulsion |

### Loop Extrusion

| Parameter | Type | Default | Units | Description |
|---|---|---|---|---|
| `LE_USE_HARMONIC_BOND` | bool | `True` | — | Enable loop bonds (requires `LOOPS_PATH`) |
| `LE_FIXED_DISTANCES` | bool | `False` | — | Fix loop distances; if `False`, scale by contact strength |
| `LE_HARMONIC_BOND_R0` | Quantity | 0.1 | nm | Loop equilibrium distance |
| `LE_HARMONIC_BOND_K` | Quantity | 30 000 | kJ mol⁻¹ nm⁻² | Loop bond stiffness |
| `LE_LOOP_FORCE_TYPE` | str | `harmonic` | — | `harmonic`, `fene_soft`, `gaussian_tether` |

### Hi-C Contact Force

| Parameter | Type | Default | Description |
|---|---|---|---|
| `HIC_USE_FORCE` | bool | `False` | Use Hi-C matrix as structural restraint |
| `HIC_PATH` | str | None | Path to `.hic`, `.cool`, or `.mcool` file |
| `HIC_NORMALIZATION` | str | `KR` | Matrix normalization: `KR`, `VC`, `VC_SQRT`, `NONE` |
| `HIC_K_SCALE` | float | 40.0 | Global energy scale (kJ mol⁻¹) — the harmonic well's stiffness. Recommended: 20–80 kJ mol⁻¹ (higher within that for `N_BEADS` ≳ 300). Above 200 kJ mol⁻¹ triggers a runtime warning. Per-pair weight is always `c_ij` itself (soft at low `c_ij`, firm at high `c_ij`) — not a separate tunable exponent. |
| `HIC_RC` | float | None | Explicit contact-radius scale `rc` [nm] — the upper clip bound for target distances, fully independent of `r_comp` (the compartment/subcompartment force's range). None (default) → auto-calibrated instead (see `HIC_AUTO_SCALE`). |
| `HIC_AUTO_SCALE` | bool | `True` | When `HIC_RC` is unset, recalibrate `rc` from the median of the initial structure's own pairwise distances, instead of a fixed nucleus-scale guess that can leave the force with no gradient. |
| `HIC_BOLTZMANN_ALPHA` | float | 4.0 | Hi-C scaling-law exponent converting contact strength to a target distance, shared by every `HIC_BOLTZMANN_KERNEL` as its steepness knob. Typical literature range: 3–4; higher values make the strength→distance mapping steeper. |
| `HIC_BOLTZMANN_KERNEL` | str | `exponential` | `P(r)` shape for the `c_ij -> r_target` inversion: `exponential` (classic Boltzmann distribution), `power_law` (Hi-C scaling law), or `sigmoid` (bounded logistic contact probability) — see the kernel table above. |
| `HIC_BOLTZMANN_TOL_FRAC` | float | 0.2 | Flat-bottom tolerance, as a fraction of each pair's own `r_target` (e.g. `0.2` → ±20% zero-force zone, harmonic beyond it). `0` (default) is the original exact two-sided well. Raise this (try `0.15–0.3`) if structures look like a spherical shell — see the note above. |
| `HIC_FORCE_OE` | bool | `False` | If `False` (default), target raw (KR-balanced) contact frequency directly. If `True`, apply Observed/Expected normalisation first, so the force optimises for relative enrichment over the distance-decay background instead — see the note above. |
| `HIC_MAX_GAP` | int | 10 | Maximum gap fraction (%) tolerated when interpolating missing bins |
| `HIC_INSULATION_WINDOW` | int | 10 | Half-width (beads) of the sliding window used by the insulation-score validation metric. Match it to your real TAD/domain size in beads — mismatched window size weakens `insulation_r` even when the force is working well. |
| `HIC_BLOCK_COPOLYMER` | bool | `False` | Opt-in: derive A/B compartments from the Hi-C matrix's own (density-aligned) PC1 and feed them into the block-copolymer force, instead of requiring a `.bed` file — see the dedicated section above. Suggested for whole-chromosome/genome-wide runs with no compartment `.bed` on hand. Only active with `HIC_USE_FORCE=True` and no `COMPARTMENT_PATH`; auto-disables below `HIC_BLOCK_COPOLYMER_MIN_BP`; errors if a bed-based compartment force is enabled too. |
| `HIC_BLOCK_COPOLYMER_MIN_BP` | float | 5,000,000 | Minimum modelled region size (bp) for `HIC_BLOCK_COPOLYMER` to stay enabled — below this the region is TAD-scale, not compartment-scale. |
| `HIC_BLOCK_COPOLYMER_STRENGTH_SCALE` | float | 0.4 | Scales `COB_EA`/`COB_EB` down for Hi-C-derived (`HIC_BLOCK_COPOLYMER`) compartments only — never affects the `.bed`-based path. `1.0` = full strength (same as `.bed`); lower softens A/B segregation. |

### Compartment and Subcompartment Forces

| Parameter | Type | Default | Units | Description |
|---|---|---|---|---|
| `COB_USE_COMPARTMENT_BLOCKS` | bool | `False` | — | A/B compartment phase separation |
| `COB_FORCE_TYPE` | str | `gaussian` | — | `gaussian`, `yukawa`, `powerlaw`, `theta` |
| `COB_EA` | float | 1.0 | kJ mol⁻¹ | A-compartment attraction strength |
| `COB_EB` | float | 2.0 | kJ mol⁻¹ | B-compartment attraction strength |
| `SCB_USE_SUBCOMPARTMENT_BLOCKS` | bool | `False` | — | A1/A2/B1/B2 subcompartment separation |
| `SCB_FORCE_TYPE` | str | `gaussian` | — | `gaussian`, `yukawa`, `powerlaw`, `theta` |
| `SCB_EA1` | float | 1.0 | kJ mol⁻¹ | A1 attraction strength |
| `SCB_EA2` | float | 1.33 | kJ mol⁻¹ | A2 attraction strength |
| `SCB_EB1` | float | 1.66 | kJ mol⁻¹ | B1 attraction strength |
| `SCB_EB2` | float | 2.0 | kJ mol⁻¹ | B2 attraction strength |

### Spherical Container

| Parameter | Type | Default | Units | Description |
|---|---|---|---|---|
| `SC_USE_SPHERICAL_CONTAINER` | bool | `False` | — | Enable spherical nuclear boundary |
| `SC_RADIUS1` | Quantity | auto | nm | Inner radius (nucleolus boundary); derived from `N` if unset |
| `SC_RADIUS2` | Quantity | auto | nm | Outer radius (nuclear boundary); derived from `N` if unset |
| `SC_SCALE` | float | 1 000 | kJ mol⁻¹ nm⁻² | Container wall stiffness |

### B-Lamina Interaction

| Parameter | Type | Default | Units | Description |
|---|---|---|---|---|
| `IBL_USE_B_LAMINA_INTERACTION` | bool | `False` | — | Anchor B-chromatin to nuclear periphery |
| `IBL_SCALE` | float | 400.0 | kJ mol⁻¹ | Lamina interaction strength |
| `BLAMINA_FORCE_TYPE` | str | `sin` | — | `sin`, `gaussian_shell`, `harmonic_shell`, `logistic_shell` |

### Chromosomal Blocks

| Parameter | Type | Default | Units | Description |
|---|---|---|---|---|
| `CHB_USE_CHROMOSOMAL_BLOCKS` | bool | `False` | — | Chromosome territory self-compaction |
| `CHB_KC` | float | 0.3 | nm⁻⁴ | Block copolymer width parameter |
| `CHB_DE` | float | 1×10⁻⁴ | kJ mol⁻¹ | Energy factor |
| `CHB_FORCE_TYPE` | str | `polynomial` | — | `polynomial`, `gaussian`, `saturating` |

### Central Force (Nucleolar Positioning)

| Parameter | Type | Default | Units | Description |
|---|---|---|---|---|
| `CF_USE_CENTRAL_FORCE` | bool | `False` | — | Size-dependent radial attraction to nucleus center |
| `CF_STRENGTH` | float | 20.0 | kJ mol⁻¹ | Attraction strength |
| `CENTRAL_FORCE_TYPE` | str | `harmonic` | — | `harmonic`, `gaussian`, `logistic` |

### Ensemble Generation

| Parameter | Type | Default | Description |
|---|---|---|---|
| `GENERATE_ENSEMBLE` | bool | `False` | Produce multiple independent structures |
| `N_ENSEMBLE` | int | None | Number of structures |
| `DOWNSAMPLING_PROB` | float | 1.0 | Fraction of loop contacts retained per structure |
| `COMPARTMENT_FLIP_PROB` | float | 0.0 | Probability of stochastic A↔B compartment flip per bead |
| `COMPARTMENT_NOISE_STD` | float | 0.0 | Gaussian noise on compartment field before discretization |

### Nucleosome Interpolation

| Parameter | Type | Default | Description |
|---|---|---|---|
| `NUC_DO_INTERPOLATION` | bool | `False` | Enable nucleosome interpolation (requires `ATACSEQ_PATH`) |
| `MAX_NUCS_PER_BEAD` | int | 4 | Maximum nucleosomes per coarse-grained bead |
| `NUC_RADIUS` | float | 0.1 | Nucleosome helix radius |
| `POINTS_PER_NUC` | int | 20 | Points per nucleosome helix |
| `PHI_NORM` | float | pi/5 | Zigzag angle |

---

## Output Directory Structure

```
OUT_PATH/
├── config_auto.ini              # copy of all parameters used
├── md_frames/
│   └── frame_1_100.cif          # trajectory frames (if SIM_RUN_MD = True)
├── metadata/
│   ├── chrom_idxs.npy
│   ├── chrom_lengths.npy
│   ├── ms.npy                   # left loop anchor bead indices
│   ├── ns.npy                   # right loop anchor bead indices
│   ├── ds.npy                   # loop equilibrium distances (nm)
│   ├── hic_validation.npy       # Hi-C validation metrics dict (if HIC_USE_FORCE = True)
│   ├── quality_tests.csv        # post-simulation quality check results
│   ├── parameters.txt           # human-readable parameter log
│   ├── MultiMM_init.cif         # initial structure
│   ├── MultiMM.psf              # UCSF Chimera topology
│   ├── MultiMM_annealing.dcd    # MD trajectory (Chimera/VMD format)
│   └── chimera_gene_coloring.cmd
├── model/
│   ├── MultiMM_minimized.cif    # energy-minimized structure
│   └── MultiMM_afterMD.cif      # structure after MD annealing
└── plots/
    ├── initial_structure.png
    ├── minimized_structure.png
    └── structure_afterMD.png
```

For genome-wide runs, an additional `chromosomes/` folder contains per-chromosome CIF files.

`hic_validation.npy` stores a dictionary with metrics keys `diagonal_decay_r`, `insulation_r`, `pc1_r`, `pearson_oe_r`, `spearman_oe_r`, `ssim`, `gmsd`, `nmi`; `_rw_*` variants hold the random-walk null-model baseline for each, and `_p` suffixes give p-values where available.

`quality_tests.csv` contains one row per quality check with columns `test`, `status` (PASS / WARN / FAIL / SKIP), `value`, and `suggestion`.

UCSF Chimera trajectory visualization: https://www.cgl.ucsf.edu/chimera/

---

## Recent Updates ❗

> Items marked *experimental* may change API or behaviour in future releases.

- **New `HIC_BLOCK_COPOLYMER` (opt-in, default `False`):** derives A/B compartments directly from the Hi-C matrix's own PC1 — sign-aligned to local contact density, discretized the same way `import_bed()` reads a `.bed` file — and feeds them into the existing block-copolymer force, fixing the Boltzmann force's main weak spot (PC1 correlation) without needing separate compartment data. Best suited to whole-chromosome/genome-wide runs with no compartment `.bed` on hand. Validation's PC1 correlation now uses this same density-based sign alignment too, so `pc1_r` is a real signed Pearson r, not `abs(r)`. See the dedicated section above for the region-size gate and the `HIC_FORCE_OE`/bed-compartment-force guardrails.
- **`HIC_FORCE_OE` defaults to `False` (raw contact frequency), `HIC_K_SCALE` default raised 20→40, plus a new `HIC_INSULATION_WINDOW` knob:** OE scored higher on isolated PC1/OE-Pearson metrics in testing, but raw frequency gave better overall results in practice (stronger diagonal decay + insulation) — see the note in the Hi-C Contact-Guided Force section. A parameter sweep across polymer sizes (100–1000 beads) also showed `HIC_K_SCALE` values of 40–80 consistently beat the old default of 20 on insulation score with no diagonal-decay cost, for `N_BEADS` ≳ 300. `HIC_INSULATION_WINDOW` (default 10 beads) lets the insulation-score window match your actual TAD/domain size.
- **MD temperature now always measured from kinetic energy:** the plotted/recorded `temperature` series is computed every frame via the equipartition theorem, `T = 2*KE / (dof * k_B)`, with `dof` accounting for particle count, constraints, and COM-motion removal — it no longer reads the integrator's thermostat set-point (`SIM_TEMPERATURE` still appears as a dashed reference line, unchanged).
- **New energy-components plot:** `plots/energy_components.png` shows each active force term's own potential energy over the MD trajectory (excluded volume, bonds, angles, loop extrusion, compartment/subcompartment/chromosomal blocks, spherical container, B-lamina, central force, Hi-C force — whichever are enabled), via dedicated OpenMM force groups, one fixed color per term.
- **Hi-C force simplified to Boltzmann-PMF only, with pluggable `P(r)` kernels:** a single classic Boltzmann-inversion restraint (`HIC_BOLTZMANN_ALPHA`, `HIC_K_SCALE`) replaces the old cross-entropy/SVD machinery. `HIC_BOLTZMANN_KERNEL` selects the equilibrium pair-distance distribution shape: `exponential` (default — the classic Boltzmann distribution), `power_law` (Hi-C scaling law), or `sigmoid` (bounded logistic contact probability); all three are exact functional inverses of their own `P(r)` and share `HIC_BOLTZMANN_ALPHA` as their steepness knob.
- **Hi-C validation suite:** seven metrics (diagonal decay, insulation, |PC1|, Pearson/Spearman OE, SSIM, GMSD, NMI) computed automatically against a random-walk baseline, saved to `metadata/hic_validation.npy`. PC1 and insulation score are heavily smoothed and range-normalized (`[-1, 1]` sign-preserving / `[0, 1]`) before correlating, so curves show the general trend (where the real minima/maxima are) rather than bead-to-bead noise, and the reported r always matches `hic_validation_curves_*.png`. Simulated/experimental/random-walk heatmaps are denoised identically — same function, same strength — as the Hi-C matrix fed into the force itself.
---

## Citation

If you use MultiMM in your research, please cite:

> Korsak, Sevastianos, Krzysztof Banecki, and Dariusz Plewczynski. "Multiscale molecular modeling of chromatin with MultiMM: From nucleosomes to the whole genome." *Computational and Structural Biotechnology Journal* 23 (2024): 3537–3548.

The software is freely distributed under the GNU General Public License v3.

For questions, bug reports, or contributions, please contact the authors or open an issue on GitHub.