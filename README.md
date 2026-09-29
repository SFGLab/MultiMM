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
- Hi-C contact-guided force field: raw `.hic` / `.cool` / `.mcool` matrices can be used directly as structural restraints via SVD or cross-entropy decomposition.
- Ensemble generation: multiple independent structures from a single run.
- Nucleosome interpolation from ATAC-seq signal.
- Validated against experimental Hi-C: distance-map correlation and compartment eigenvector correlation are computed automatically.

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

Chromatin is represented as a coarse-grained polymer. The total energy $E$ decomposes into physically motivated terms:

$$E = E_{\text{backbone}} + E_{\text{loops}} + E_{\text{block}} + E_{\text{excluded}} + E_{\text{confinement}} + E_{\text{chromosomal}} + E_{\text{Hi-C}}$$

Each term encodes a distinct biological mechanism:

| Term | Mechanism |
|---|---|
| $E_{\text{backbone}}$ | Polymer connectivity and stiffness |
| $E_{\text{loops}}$ | Long-range loop extrusion or experimental contact constraints |
| $E_{\text{block}}$ | Compartment and subcompartment phase separation |
| $E_{\text{excluded}}$ | Steric repulsion between beads |
| $E_{\text{confinement}}$ | Nuclear geometry: spherical container and lamina affinity |
| $E_{\text{chromosomal}}$ | Chromosome territory formation and global compaction |
| $E_{\text{Hi-C}}$ | Data-driven contact restraint from a raw Hi-C matrix |

---

### Polymer Backbone

The backbone encodes chain connectivity and local rigidity through two standard terms.

**Harmonic bond** (nearest-neighbor connectivity):

$$E_{\text{bond}} = \sum_{i} \frac{k_b}{2} \left( r_{i,i+1} - r_0 \right)^2$$

where $r_{i,i+1}$ is the distance between consecutive beads, $r_0$ is the equilibrium bond length, and $k_b$ is the bond stiffness.

**Harmonic angle** (chain stiffness / persistence length):

$$E_{\text{angle}} = \sum_{i} \frac{k_\theta}{2} \left( \theta_i - \theta_0 \right)^2$$

where $\theta_i$ is the angle formed by three consecutive beads $(i, i+1, i+2)$, $\theta_0$ is the preferred angle, and $k_\theta$ controls bending rigidity. Together these reproduce a discretized worm-like chain.

---

### Loop Interactions

Long-range loop constraints tether pairs of beads $(m, n)$ identified from loop-calling experiments. Three functional forms are available.

**Harmonic** (default):

$$E_{\text{loops}}^{\text{harmonic}} = \sum_{(m,n)} \frac{k}{2} \left( r_{mn} - r_0 \right)^2$$

**Soft FENE-like** (bounded, avoids divergence at large extension):

$$E_{\text{loops}}^{\text{fene}} = \sum_{(m,n)} \frac{k (r_{mn} - r_0)^2}{1 + \alpha (r_{mn} - r_0)^2}$$

**Gaussian tether** (smooth, fully bounded):

$$E_{\text{loops}}^{\text{gaussian}} = \sum_{(m,n)} k \left( 1 - e^{-(r_{mn} - r_0)^2 / \sigma^2} \right)$$

Loop bond strengths can be fixed (`LE_FIXED_DISTANCES = True`) or scaled by experimental contact frequency (`LE_FIXED_DISTANCES = False`). Providing `LOOPS_PATH` is optional; the simulation runs without loops if it is not supplied.

---

### Block-Copolymer Compartmentalization

State-dependent pairwise attractions drive compartment phase separation. Each bead carries a label $s_i$ representing compartment identity.

**Compartment level (A/B):**

$$E_{\text{comp}} = -\sum_{i<j} \epsilon(s_i, s_j) \exp\!\left( -\frac{r_{ij}^2}{2r_c^2} \right)$$

The coupling $\epsilon(s_i, s_j)$ is attractive for like compartments (A–A, B–B) and weak or repulsive otherwise, reproducing large-scale A/B segregation.

**Subcompartment level (A1/A2/B1/B2):**

$$E_{\text{sub}} = -\sum_{i<j} \epsilon_{\alpha\beta} \exp\!\left( -\frac{r_{ij}^2}{2r_{sc}^2} \right), \quad s_i = \alpha,\; s_j = \beta$$

This promotes finer microphase separation inside A/B compartments.

**Chromosome territory (self-compaction):**

$$E_{\text{chrom}} = \sum_{i<j} \delta_{\chi_i, \chi_j} \, V(r_{ij})$$

where $\chi_i$ is the chromosome label and $V(r)$ is a soft attractive potential acting only between beads on the same chromosome. The default polynomial form is:

$$V(r) = dE \left( k_C r^4 - r^3 + r^2 \right)$$

**Alternative interaction kernels** (experimental):

| Kernel | Expression | Effect |
|---|---|---|
| Yukawa | $V(r) \sim -e^{-r/\lambda}/r$ | Screened, longer-range |
| Power-law | $V(r) \sim -1/(r^\alpha + \varepsilon)$ | Scale-free attraction |
| Theta (contact) | $V(r) \sim -\Theta(r_c - r)$ | Binary hard-cutoff |
| Saturating | $V(r) \sim -1/(1 + k_C r^2)$ | Bounded, prevents over-collapse |

---

### Nuclear Geometry and Lamina Interactions

**Spherical container** — soft penalty for excursions outside the nuclear shell:

$$E_{\text{container}} = C \sum_i \left[ \max(0,\, r_i - R_2)^2 + \max(0,\, R_1 - r_i)^2 \right]$$

where $r_i$ is the radial distance from the nuclear center. This confines chromatin between radii $R_1$ and $R_2$.

**B-lamina interaction** — anchors B-compartment chromatin to the nuclear periphery:

$$E_{\text{lamina}} = -\sum_i B(s_i)\, V(r_i)$$

where $B(s_i)$ selects B-compartment beads. Available radial profiles:

| Mode | Expression | Description |
|---|---|---|
| `sin` (default) | $V(r) = \sin^8\!\left(\tfrac{\pi(r-R_1)}{R_2-R_1}\right) - 1$ | Sharp peripheral preference |
| `gaussian_shell` | $V(r) \sim -\left[e^{-(r-R_1)^2/2\sigma^2} + e^{-(r-R_2)^2/2\sigma^2}\right]$ | Localized at both boundaries |
| `harmonic_shell` | $V(r) \sim (r - r_0)^2,\quad r_0 = (R_1+R_2)/2$ | Pulls toward mid-shell |
| `logistic_shell` | Smooth sigmoidal walls at $R_1$ and $R_2$ | Smooth boundary transition |

---

### Hi-C Contact-Guided Force

When a raw Hi-C contact matrix is provided (`HIC_PATH`), it is used directly as a structural restraint through a low-rank decomposition. This is an alternative or complement to loop-extrusion forces and works without requiring explicit loop calls.

**Pipeline:**

1. The matrix is loaded from `.hic`, `.cool`, or `.mcool`, auto-selecting resolution and resampling to exactly $N_\text{beads} \times N_\text{beads}$ via weighted average pooling (`read_hic.py`).
2. Knight–Ruiz iterative balancing: $H \leftarrow D^{-1} H D^{-1}$ until row marginals ≈ 1.
3. Observed/Expected normalisation: $OE[i,j] = \log\!\left(H[i,j] / E[i,j]\right)$, where $E[i,j]$ is the genome-wide mean contact frequency at separation $|i-j|$.
4. Truncated eigendecomposition of the symmetric OE matrix: $OE \approx \sum_{k=1}^{K} \lambda_k \mathbf{v}_k \mathbf{v}_k^T$.

Per-particle scalars encoding sign and magnitude:

$$a_k(i) = \text{sign}(\lambda_k) \cdot \sqrt{|\lambda_k|} \cdot v_k(i)$$

Three force modes are available (`HIC_FORCE_MODE`):

**`svd`** — single-scale Gaussian CustomNonbondedForce:

$$U(i,j;\,r) = -k \left[\sum_{k=1}^K a_k(i)\, a_k(j)\right] \exp\!\left(-\frac{r^2}{2\sigma^2}\right)$$

where $\sigma = r_c / 3$ and the cutoff is $3\sigma$. Same-type beads accumulate a positive dot product and attract; opposite types repel.

**`svd_multiscale`** — each eigenvector acts at its own spatial scale $\sigma_k$:

$$U(i,j;\,r) = -k \sum_{k=1}^K a_k(i)\, a_k(j) \exp\!\left(-\frac{r^2}{2\sigma_k^2}\right)$$

The per-component length scales $\sigma_k$ can be assigned by eigenvalue magnitude ($\sigma_k = \sigma_{\max}(|\lambda_k|/|\lambda_1|)^\beta$) or by eigenvector genomic autocorrelation length.

**`crossentropy`** — sparse CustomBondForce for sharp contact peaks:

$$P(r) = \frac{1}{1 + (r/r_c)^\alpha}, \qquad U_{ij}(r) = -k\left[c_{ij}\log(P + \varepsilon) + (1-c_{ij})\log(1-P+\varepsilon)\right]$$

Only pairs with $c_{ij} \geq \text{threshold}$ receive a bond, keeping the bond list sparse.

**Hi-C validation** — after simulation, MultiMM automatically computes:
- Pearson $r$ between the model distance map and $\log(1 + \text{Hi-C})$ contact frequency (upper triangle, excluding near-diagonal bands).
- Pearson $r$ between $|\text{PC1}|$ of the model and $|\text{PC1}|$ of the experimental Hi-C matrix.

Results are saved to `metadata/hic_validation.npy`.

---

### Nucleosome-Scale Interpolation

After coarse-grained optimization, nucleosome positions are interpolated using a beads-on-a-string zigzag model. Each nucleosome is represented as a helix with 1.65 DNA turns. The number of nucleosomes per bead is derived from normalized ATAC-seq signal, enforcing nucleosome-rich regions in low-accessibility chromatin.

---

## Internal Parameter Definitions

All geometric and interaction scales are derived from a single microscopic length scale — the polymer bond length $b_0$ (`POL_HARMONIC_BOND_R0`).

**Nuclear radius** (dense globule scaling):

$$R_2 = b_0\, N^{1/3}$$

which enforces constant monomer density $N/R_2^3 \approx \text{const}$.

**Inner compartment radius** (fixed volume fraction $f$):

$$R_1 = R_2\, f^{1/3}$$

**Compartment interaction length scale:**

$$r_c \sim \mathcal{O}(b_0) \approx 1.5\, b_0$$

ensuring interactions remain local relative to the polymer backbone.

**Loop equilibrium distances** — either globally fixed or derived from experimental loop lengths $d_i$:

$$r_0^{(i)} \in \{r_0^{\text{global}},\; d_i\}$$

**State variables** — each bead carries a compartment label $s_i$ and chromosome label $\chi_i$, which modulate interactions through selection rules:

$$E_{ij} \propto \delta(s_i,\, s_j), \qquad E_{ij} \propto \delta(\chi_i,\, \chi_j)$$

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

Any standard Juicer `.hic` file or Cooler `.cool` / `.mcool` file is accepted. MultiMM automatically selects the best available resolution for the requested region and resamples to $N_\text{beads} \times N_\text{beads}$. KR, VC, VC_SQRT, and NONE normalizations are supported.

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
import simulation.plots as splt
splt.viz_chroms(sim_path)          # add comps=False to disable compartment coloring
```

**Any single region or CIF structure:**

```python
import simulation.plots as splt
import simulation.utils as suts

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
| `SIM_INTEGRATOR_TYPE` | str | `langevin` | — | `langevin`, `verlet`, `brownian` |
| `SIM_INTEGRATOR_STEP` | Quantity | 1 | fs | Integrator time step |
| `SIM_FRICTION_COEFF` | float | 0.5 | ps⁻¹ | Friction coefficient (Langevin / Brownian) |
| `SIM_SET_INITIAL_VELOCITIES` | bool | `False` | — | Initialize velocities from Boltzmann distribution |
| `TRJ_FRAMES` | int | 2 000 | — | Total trajectory frames to save |

### Polymer Backbone

| Parameter | Type | Default | Units | Description |
|---|---|---|---|---|
| `POL_USE_HARMONIC_BOND` | bool | `True` | — | Harmonic bond between consecutive beads |
| `POL_HARMONIC_BOND_R0` | Quantity | 0.1 | nm | Equilibrium bond length |
| `POL_HARMONIC_BOND_K` | Quantity | 300 000 | kJ mol⁻¹ nm⁻² | Bond stiffness |
| `POL_USE_HARMONIC_ANGLE` | bool | `True` | — | Harmonic angle (bending rigidity) |
| `POL_HARMONIC_ANGLE_R0` | Quantity | π | rad | Equilibrium angle |
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
| `HIC_FORCE_MODE` | str | `crossentropy` | `svd`, `svd_multiscale`, `crossentropy` |
| `HIC_NORMALIZATION` | str | `KR` | Matrix normalization: `KR`, `VC`, `VC_SQRT`, `NONE` |
| `HIC_N_COMPONENTS` | int | 10 | Number of SVD eigenvectors retained |
| `HIC_K_SCALE` | float | 1.0 | Global energy scale (kJ mol⁻¹) |
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
| `SC_RADIUS1` | Quantity | auto | nm | Inner radius (nucleolus boundary); derived from $N$ if unset |
| `SC_RADIUS2` | Quantity | auto | nm | Outer radius (nuclear boundary); derived from $N$ if unset |
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
| `PHI_NORM` | float | π/5 | Zigzag angle |

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
│   ├── ds.npy                   # loop strengths (as distances)
│   ├── hic_validation.npy       # Hi-C validation metrics (if HIC_USE_FORCE = True)
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

`hic_validation.npy` stores a dictionary with keys `heatmap_r`, `heatmap_p` (Pearson $r$ and $p$-value between model distance map and Hi-C) and `eigvec_r`, `eigvec_p` ($|$PC1$|$ correlation).

UCSF Chimera trajectory visualization: https://www.cgl.ucsf.edu/chimera/

---

## Citation

If you use MultiMM in your research, please cite:

> Korsak, Sevastianos, Krzysztof Banecki, and Dariusz Plewczynski. "Multiscale molecular modeling of chromatin with MultiMM: From nucleosomes to the whole genome." *Computational and Structural Biotechnology Journal* 23 (2024): 3537–3548.

The software is freely distributed under the GNU General Public License v3.

For questions, bug reports, or contributions, please contact the authors or open an issue on GitHub.
