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
- Hi-C contact-guided force field: raw `.hic` / `.cool` / `.mcool` matrices used directly as structural restraints via a **cross-entropy** CustomBondForce that models contact probability with a pluggable distance kernel (`gaussian` by default; `power_law`/`sigmoid`, `exponential`, `erfc`, and separation-aware `rouse` also available) and weights each pair's force by its own observed contact strength, so weakly-supported pairs stay soft rather than acting as hard constraints.
- Ensemble generation: multiple independent structures from a single run.
- Nucleosome interpolation from ATAC-seq signal.
- **Comprehensive Hi-C validation** computed automatically: diagonal decay correlation, insulation score correlation, |PC1| compartment correlation, Pearson / Spearman OE-matrix correlation, SSIM, GMSD, NMI — each reported alongside a random-walk null-model baseline for direct comparison.
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

When a raw Hi-C contact matrix is provided (`HIC_PATH`), it is used directly as a structural restraint via a sparse cross-entropy `CustomBondForce`. This is an alternative or complement to loop-extrusion forces and works without requiring explicit loop calls.

**Pipeline:**

1. The matrix is loaded from `.hic`, `.cool`, or `.mcool`, auto-selecting resolution and resampling to exactly `N_beads x N_beads` via weighted average pooling (`read_hic.py`).
2. Knight–Ruiz iterative balancing: `H <- D^{-1} H D^{-1}` until row marginals ≈ 1, giving the normalised contact matrix `c_ij ∈ [0, 1]`. If `HIC_FORCE_OE=True`, an additional Observed/Expected normalisation (divide each diagonal by its mean) is applied on top, so `c_ij` targets relative enrichment over the distance-decay background instead of absolute contact frequency (default `False`).
3. A sparse set of bonds is built for every pair with `c_ij ≥ HIC_THRESHOLD` (a sparsity cutoff only — O(M) bonds, M ≪ N²).

Instead of converting contacts into fixed target distances, the force models the contact *probability* directly as a distance kernel `P_ij(r)` and minimises the binary cross-entropy (negative log-likelihood) against the observed `c_ij`:

```
U_ij(r) = -k * w_ij * [c_ij * log(P + eps) + (1 - c_ij) * log(1 - P + eps)]
w_ij    = c_ij ^ HIC_WEIGHT_POWER
```

The force is `F_i = -∇_i U`, a self-regulating residual: it pushes pairs apart when the model overshoots the observed contact probability (`P > c_ij`) and pulls them together when it undershoots (`P < c_ij`), vanishing once the structure matches the data. The weight `w_ij` additionally scales *both* branches of the loss by the pair's own observed contact strength, so pairs barely above `HIC_THRESHOLD` exert a correspondingly tiny force (soft) while well-supported contacts (`c_ij` near 1) behave like a firm restraint (hard) — with `HIC_WEIGHT_POWER = 1.0` (default) the force is directly proportional to `c_ij`.

`P_ij(r)` is selected with `HIC_KERNEL`:

| Kernel | `P(r)` | Notes |
|---|---|---|
| `gaussian` (default) | `exp(-r² / 2σ²)` | Width `HIC_GAUSSIAN_SIGMA` (defaults to `r_comp`). Decays faster than power-law at large `r`; good for short-range loops. |
| `power_law` / `sigmoid` | `1 / (1 + (r/r_c)^α)` | Steepness `HIC_POWERLAW_ALPHA`. At `α=2` a Lorentzian; large `α` approaches a step. |
| `exponential` | `exp(-r/r_c)` | Constant log-derivative — force doesn't vanish at long range (Yukawa-like). |
| `erfc` | `0.5·erfc((r - r_c) / (√2·σ_s))` | Soft step centred on `r_c`, width `HIC_ERFC_SIGMA`; closest to a binary contact definition. |
| `rouse` | `erfc(r / √(2·s·b²))` | Separation-aware: `s = \|i-j\|` in beads, Kuhn length `HIC_ROUSE_KUHN_LENGTH`. Each pair's effective contact threshold follows the Gaussian-chain law, automatically tracking the expected diagonal decay. |

Each kernel's extra parameter is named `HIC_<KERNEL>_*` and only takes effect when that kernel is selected via `HIC_KERNEL` — e.g. setting `HIC_ERFC_SIGMA` has no effect unless `HIC_KERNEL=erfc`.

All kernels plug into the same weighted cross-entropy loss above; only `P(r)` changes. Recommended `HIC_K_SCALE`: 5–20 kJ mol⁻¹ (weak contacts are already softened independently via `HIC_WEIGHT_POWER`, so `HIC_K_SCALE` mainly sets how firm the *well-supported* contacts are).

> **Note on compartments (PC1):** this per-pair force optimises individual contact probabilities, not the mean diagonal decay profile or the spectral (PC1) structure of the full contact map — a structure where every `P_ij ≈ c_ij` can still have the wrong A/B compartment signal, since PC1 is a covariance property of the whole matrix rather than of individual pairs. If `|PC1| r` validation is weak, pair this force with the dedicated compartment force (`COB_USE_COMPARTMENT_BLOCKS`, driven by `COMPARTMENT_PATH`) rather than expecting the Hi-C force alone to recover compartment identity.

**Hi-C validation** — after simulation, MultiMM automatically computes seven metrics comparing the model contact map against the experimental Hi-C matrix, each reported alongside a random-walk null-model baseline:

| Metric | What it measures |
|---|---|
| Diagonal decay correlation | Whether contact frequency falls off with genomic distance at the right rate |
| Insulation score correlation | Agreement of TAD boundary positions |
| \|PC1\| compartment correlation | A/B compartment identity (Pearson r of \|PC1\| vectors) |
| Pearson / Spearman OE correlation | Global agreement of the observed-over-expected matrices |
| SSIM | Structural similarity (local contrast, luminance, structure) |
| GMSD | Edge/boundary sharpness agreement |
| NMI | Normalised mutual information |

Results are saved to `metadata/hic_validation.npy` (a dict; keys: `diagonal_decay_r`, `insulation_r`, `pc1_r`, `pearson_oe_r`, `spearman_oe_r`, `ssim`, `gmsd`, `nmi`, plus `_rw_*` variants for the null baseline and `_p` suffix for p-values where applicable).

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

Ready-to-use configuration files for common scenarios are in the `examples/` folder.

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
| `TRJ_FRAMES` | int | 2 000 | — | Total trajectory frames to save |

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
| `HIC_KERNEL` | str | `gaussian` | Distance→probability kernel: `gaussian`, `power_law`/`sigmoid`, `exponential`, `erfc`, `rouse` |
| `HIC_K_SCALE` | float | 5.0 | Global energy scale (kJ mol⁻¹). Recommended: 5–20 kJ mol⁻¹. Values > 30 kJ mol⁻¹ risk freezing MD thermal sampling (runtime warning). |
| `HIC_WEIGHT_POWER` | float | 1.0 | Exponent β in per-pair weight `c_ij^β`; 1.0 = force proportional to observed contact strength (soft at low `c_ij`, firm at high `c_ij`) |
| `HIC_POWERLAW_ALPHA` | float | 3.0 | Sigmoid steepness; only takes effect when `HIC_KERNEL=power_law`/`sigmoid` |
| `HIC_GAUSSIAN_SIGMA` | float | None (→`r_comp`) | Width σ; only takes effect when `HIC_KERNEL=gaussian` |
| `HIC_ERFC_SIGMA` | float | None (→`0.3·r_comp`) | Softening width; only takes effect when `HIC_KERNEL=erfc` |
| `HIC_ROUSE_KUHN_LENGTH` | float | None (→`r_comp`) | Kuhn length `b`; only takes effect when `HIC_KERNEL=rouse` |
| `HIC_THRESHOLD` | float | 0.01 | Minimum `c_ij` to build a bond at all (sparsity cutoff only; force strength still scales continuously with `c_ij` via `HIC_WEIGHT_POWER`) |
| `HIC_FORCE_OE` | bool | `False` | If `True`, apply Observed/Expected normalisation to the Hi-C matrix before it is passed as the `c_ij` target, so the force optimises for relative contact enrichment over the distance-decay background rather than absolute contact frequency. In general leave this `False`. |
| `HIC_MAX_GAP` | int | 10 | Maximum gap fraction (%) tolerated when interpolating missing bins |

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

- **Hi-C validation bugfix:** `get_coordinates_cif` only read `ATOM` lines, silently dropping MultiMM's two `HETATM` chain-terminus beads from every loaded structure and desyncing bead↔genomic-bin alignment in Hi-C validation (the cause of a random-walk baseline scoring as well as, or better than, the real simulation on PC1/Pearson/SSIM/NMI). Now reads `ATOM` and `HETATM`. Also added light Gaussian smoothing (`smooth_sigma`, default 1 bead, 0 to disable) to sim/exp/RW matrices before computing metrics.
- **Hi-C force — pluggable kernels, soft weighting:** `HIC_KERNEL` selects `gaussian` (new default), `power_law`/`sigmoid`, `exponential`, `erfc`, or `rouse`. Force is now weighted by `c_ij^HIC_WEIGHT_POWER` (default 1.0) so weak contacts stay soft instead of acting as hard constraints; `HIC_THRESHOLD` is sparsity-only now. Backward compatible via `HIC_KERNEL=power_law`, `HIC_WEIGHT_POWER=0`.
- **Hi-C force simplified to crossentropy-only:** SVD-based modes (`svd`, `svd_multiscale`) removed; a single, sparser, numerically-stabler cross-entropy `CustomBondForce` remains, tuned via `HIC_K_SCALE`.
- **SVD force sigma bugfix:** `sigma` was `r_comp/3` instead of a fraction of the nuclear radius, making old Hi-C-guided SVD runs indistinguishable from a random walk. Fixed; re-run pre-fix SVD results. Defaults also updated (`HIC_N_COMPONENTS` 5→15, `HIC_K_SCALE` 130→20).
- **SVD force anti-collapse controls:** added `HIC_SIGMA_SCALE` (default 0.20) and `HIC_MIN_SEPARATION` (default 5) to stop the dominant eigenvector collapsing the structure to a point.
- **Hi-C contact force** *(experimental)*: restrain directly from a raw contact matrix (KR-balanced, OE-normalised) instead of called loops.
- **Coloured structured logger:** time-stamped, colour-coded log lines, section banners, and summary tables.
- **Hi-C validation suite:** seven metrics (diagonal decay, insulation, |PC1|, Pearson/Spearman OE, SSIM, GMSD, NMI) computed automatically against a random-walk baseline, saved to `metadata/hic_validation.npy`.
- **Quality control helpers** *(experimental)*: up to 12 post-simulation checks (energy, bonds/angles, excluded volume, compartments, territories, loops, confinement) → `metadata/quality_tests.csv`.
- **Additional data modalities:** `.cool`/`.mcool` support alongside `.hic`, KR/VC/VC_SQRT/NONE normalisation, ensemble generation with compartment noise and loop downsampling.

---

## Citation

If you use MultiMM in your research, please cite:

> Korsak, Sevastianos, Krzysztof Banecki, and Dariusz Plewczynski. "Multiscale molecular modeling of chromatin with MultiMM: From nucleosomes to the whole genome." *Computational and Structural Biotechnology Journal* 23 (2024): 3537–3548.

The software is freely distributed under the GNU General Public License v3.

For questions, bug reports, or contributions, please contact the authors or open an issue on GitHub.
