# Units and scales

Every force in MultiMM needs a length or an energy to be meaningful ("how far apart are bonded beads?", "how big is the nucleus?", "how strong is an attraction compared with thermal motion?"). Rather than asking you to choose all of these independently, MultiMM derives almost every scale from **one microscopic length and the number of beads**. This page lists those definitions in one place.

!!! tip "In plain words"
    It is like the scale on a map: we choose one length (the distance between neighbouring beads) and then measure every other size, such as the radius of the nucleus, in terms of it. You only choose how many beads you want; MultiMM works out the rest.

## One length unit: the bond length

The basic length is the **bond length** $b_0$ (`POL_HARMONIC_BOND_R0`, default $0.1$ nm). It is the preferred distance between consecutive beads, and the other geometric quantities are written in terms of it.

!!! warning "Lengths are model units"
    A bead represents a stretch of chromatin of several kilobases, whose true physical size is tens of nanometres, but $b_0$ is only $0.1$ nm. The "nm" in the configuration is a label for the length unit of the simulation, chosen so that the numbers behave well in OpenMM. Compare structures with each other, and with Hi-C, in **relative** terms. To compare with microscopy, rescale lengths by the physical size of one bead (or by the radius of your nucleus).

## Genomic resolution

The number of beads $N$ (`N_BEADS`) sets the genomic length represented by one bead:

$$
\text{resolution}=\frac{\text{modelled genomic length}}{N}\ \text{bp per bead}.
$$

For a region the length is `LOC_END - LOC_START`; for a whole chromosome it is its length; for the whole genome it is the sum of the chromosome lengths. Examples: a 13 Mb region with $N=5000$ gives about 2.6 kb per bead; the autosomes with $N=50\,000$ give about 60 kb per bead.

All input data are mapped onto this grid. A genomic coordinate $x$ becomes the bead index $\lfloor (x-x_\text{start})/\text{resolution}\rfloor$. A loop anchor, a compartment interval and a Hi-C bin are all converted this way.

## The geometric scales

From $b_0$ and $N$ MultiMM computes the size of the nucleus. The idea is a **dense globule of constant density**: a volume that grows in proportion to the number of beads, so the radius grows as $N^{1/3}$.

| Quantity | Definition | Meaning |
|---|---|---|
| Nuclear radius $R_2$ | $R_2=b_0\,N^{1/3}$ | Outer wall of the nucleus (nuclear envelope, lamina) |
| Inner radius $R_1$ | $R_1=R_2\,f^{1/3}$ with $f=0.2$ | Surface of an excluded central sphere (the nucleolus) that holds 20% of the volume |
| Compartment range $r_c$ | $r_c=R_2/3$ | Width of the compartment and subcompartment attraction |
| Contact cutoff $r_\text{cut}$ | $r_\text{cut}=5\,b_0$ | Hard cutoff beyond which compartment forces are ignored |
| Excluded-volume size $\sigma$ | `LE_HARMONIC_BOND_R0` (default 0.1) | Length scale of bead repulsion |
| Loop distance $d$ | between $0.1$ and $0.2$ | Rest length of a loop bond (see [Loops](loops.md)) |
| Hi-C contact radius $r_\text{c}^\text{Hi-C}$ | median pairwise distance of the initial structure, or `HIC_RC` | Largest target distance of the Hi-C force |
| Hi-C minimum distance $r_\text{min}$ | `LE_HARMONIC_BOND_R0` | Smallest target distance of the Hi-C force |

The nucleus is **centred on the centroid of the initial structure**, which is the origin of the radial coordinate $r$ in all radial forces.

!!! info "In MultiMM"
    Only $b_0$ and $N$ are free. The nucleolus fraction $f=0.2$ is fixed. $R_1$ and $R_2$ always follow from $b_0$ and $N$; the keys `SC_RADIUS1` and `SC_RADIUS2` exist but are not used.

A worked example, for the default $b_0=0.1$ nm:

| $N$ | $R_2$ | $R_1$ | $r_c=R_2/3$ | $r_\text{cut}$ |
|---|---|---|---|---|
| 1 000 | 1.00 | 0.58 | 0.33 | 0.5 |
| 5 000 | 1.71 | 1.00 | 0.57 | 0.5 |
| 50 000 | 3.68 | 2.15 | 1.23 | 0.5 |
| 200 000 | 5.85 | 3.42 | 1.95 | 0.5 |

!!! note "Why the cutoff does not grow with $N$"
    $R_2$ and $r_c$ grow with the system, but the compartment cutoff stays at $5\,b_0$. Compartment forces are meant to be **local contact attractions**: without a cutoff, OpenMM would sum every pair in the system, and the total attraction would then grow with $N$ instead of with the local density, which makes large structures collapse. Inside the cutoff the Gaussian profile is shallow, so the force acts almost like a contact attraction. See [Compartments](compartments.md).

## Energy and temperature

Energies are in **kJ/mol** and the temperature is 310 K (`SIM_TEMPERATURE`). The thermal energy is

$$
k_BT = 0.008314\ \tfrac{\text{kJ}}{\text{mol K}}\times 310\ \text{K}\approx 2.58\ \text{kJ/mol}.
$$

It helps to read force strengths in units of $k_BT$:

| Force | Default | In units of $k_BT$ |
|---|---|---|
| Bond stiffness $k_b$ | $3\times10^{5}$ kJ mol⁻¹ nm⁻² | thermal bond stretch $\sqrt{k_BT/k_b}\approx0.003$ nm, about 3% of $b_0$: bonds are nearly rigid |
| Angle stiffness $k_\theta$ | 100 kJ mol⁻¹ rad⁻² | rms bending of about 0.23 rad per joint |
| Excluded volume $\varepsilon$ | 100 kJ/mol | repulsion of $\approx3.4\,k_BT$ at $r=b_0$, $\approx39\,k_BT$ at $r=b_0/2$ |
| Compartment attraction $E_A$, $E_B$ | 1 and 2 kJ/mol | 0.39 and 0.78 $k_BT$ per contact: weak, so compartments form collectively |
| Hi-C stiffness | 40 | see [Hi-C force](hic-force.md) |

Individual attractions are much weaker than $k_BT$ on purpose. A bead has many neighbours, so many weak contacts add up to a strong collective preference, which is the mechanism behind phase separation.

## Mass, time and friction

Each bead carries the theoretical mass of 10 kb of chromatin, $16\,427.889$ amu (set in the bundled force field file `ff.xml`). The integration step is `SIM_INTEGRATOR_STEP` (default 1 fs) and the Langevin friction is `SIM_FRICTION_COEFF` (default 0.5 ps⁻¹). Because lengths are in model units, simulated time should be read as an optimisation clock rather than as real time.

## Writing quantities in the config file

Parameters that carry a unit are written as a number followed by an OpenMM unit name:

```ini
POL_HARMONIC_BOND_R0  = 0.1 nanometer
POL_HARMONIC_BOND_K   = 300000.0 kilojoules_per_mole/nanometer**2
SIM_INTEGRATOR_STEP   = 1 femtosecond
SIM_TEMPERATURE       = 310 kelvin
```

Parameters shown as plain numbers (`EV_EPSILON`, `COB_EA`, `HIC_K_SCALE`, ...) are interpreted in the unit given in the [argument table](../reference/arguments.md).
