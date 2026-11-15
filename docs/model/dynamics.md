# Minimisation and molecular dynamics

The forces described in the previous pages define an **energy** $E$ for every possible shape of the chromatin chain. This page explains the next step: how MultiMM finds good shapes, and how it can then make the chain move like a real, jiggling polymer. Everything here is controlled by a handful of `SIM_*` settings, and the defaults work for most runs.

!!! tip "In plain words"
    Think of $E$ as a landscape of hills and valleys, where every position on the map is one possible folding of the chain. **Minimisation** is a hiker walking downhill to the bottom of the nearest valley: you get one comfortable shape. **Molecular dynamics (MD)** is the same hiker, but now warm and restless: the chain keeps jiggling around the valley floor, and the hiker takes a photo every now and then. You get a *movie* of plausible shapes.

## Why two steps?

Starting from the [initial structure](initial-structure.md), the chain is compact and tidy, but it does not yet agree with your data. Minimisation fixes this quickly. It is cheap, deterministic, and for large systems (a whole chromosome or genome) it is often all you need.

A real chromosome, however, is never frozen in a single shape. Thermal motion makes it fluctuate constantly, and the cells in a sample each have a slightly different structure. MD mimics this. It is most useful for small and medium systems (a gene, a region) where you want

- an **ensemble** of structures rather than one,
- a measure of **how flexible** each part is (the RMSF, see below),
- a check that the structure is **stable**, not just a lucky point of the landscape.

## Minimisation

MultiMM calls OpenMM's `minimizeEnergy`, which uses the **L-BFGS** method: at every step it moves all beads a little in the direction that lowers $E$ the fastest, and it uses the memory of earlier steps to choose smart step sizes. It stops when the largest force on any bead is below a tolerance (10 kJ/mol/nm, the OpenMM default). The result is written to `model/MultiMM_minimized.cif`.

Minimisation always finishes in a *local* minimum, the nearest valley, which is why a good starting structure matters.

## What is integration?

To make the chain move we need Newton's law: *force = mass × acceleration*. For bead $i$ with position $\mathbf{r}_i$ and velocity $\mathbf{v}_i$,

$$
m\,\frac{d\mathbf{v}_i}{dt}=\mathbf{F}_i=-\nabla_{\mathbf{r}_i}E,
\qquad
\frac{d\mathbf{r}_i}{dt}=\mathbf{v}_i .
$$

In words: the energy tells us the force (a bead is pushed in the direction in which $E$ decreases), the force changes the velocity, and the velocity changes the position. For more than two beads there is no formula for the answer, so the computer advances time in many tiny steps of length $\Delta t$ (the **time step**). **Integration** means repeating this update step after step, and an **integrator** is the recipe for one update.

A recipe has to be good in two ways. It must be *accurate*, so that errors do not pile up, and it must be *stable*, so that the energy does not slowly drift. The simplest recipe (move along the current velocity, then update the velocity, called Euler) is neither, which is why MD programs use better ones such as the **Verlet** scheme:

$$
\mathbf{v}\leftarrow\mathbf{v}+\tfrac{\Delta t}{2m}\mathbf{F}(\mathbf{r}),\qquad
\mathbf{r}\leftarrow\mathbf{r}+\Delta t\,\mathbf{v},\qquad
\mathbf{v}\leftarrow\mathbf{v}+\tfrac{\Delta t}{2m}\mathbf{F}(\mathbf{r}_\text{new}) .
$$

First kick the velocity with the current force, then let the bead drift, then kick again with the new force. (OpenMM uses the equivalent leap-frog arrangement, which gives the same trajectory.) The Verlet scheme conserves energy very well, but it has no temperature control: it keeps whatever energy the system started with.

### The time step

The time step must be small compared with the fastest motion in the system, which here is the vibration of the stiff bonds. A step that is too large makes the energy explode (you will see `NaN` in the log). MultiMM uses **1 femtosecond** (`SIM_INTEGRATOR_STEP`), which is very safe. Each bead has the mass of 10 kb of chromatin (16 428 atomic mass units).

!!! warning "Time is relative"
    The lengths of the model are not physical nanometres ([Units and scales](scales.md)), so the simulated time is not real time either. Use MD to compare *relative* fluctuations and to produce an ensemble, not to measure absolute speeds.

### Temperature

Real chromatin sits in a bath at about 310 K, which keeps it jiggling. The integrators below that have a **thermostat** (Langevin, Brownian) imitate this: they add a little friction and random kicks, and the balance between the two keeps the average kinetic energy at $\tfrac32k_BT$ per bead. The temperature is set by `SIM_TEMPERATURE` (default 310 K); it never needs to be changed unless you want to explore other conditions. Remember that $k_BT\approx2.58$ kJ/mol at 310 K is the yardstick for how strong an interaction is.

## The integrators

`SIM_INTEGRATOR_TYPE` chooses the recipe. The default is `langevin`.

| `SIM_INTEGRATOR_TYPE` | Thermostat | Settings it uses | Use it for |
|---|---|---|---|
| `langevin` (**default**) | yes | `SIM_TEMPERATURE`, `SIM_FRICTION_COEFF`, `SIM_INTEGRATOR_STEP` | almost everything |
| `brownian` | yes (overdamped) | `SIM_TEMPERATURE`, `SIM_FRICTION_COEFF`, `SIM_INTEGRATOR_STEP` | slow, viscous motion; no inertia |
| `verlet` | no | `SIM_INTEGRATOR_STEP` | tests; checking energy conservation |
| `variable_langevin` | yes | `SIM_TEMPERATURE`, `SIM_FRICTION_COEFF`, `SIM_ERROR_TOLERANCE` | automatic time step with a thermostat |
| `variable_verlet` | no | `SIM_ERROR_TOLERANCE` | automatic time step, no thermostat |
| `amd` | no | `SIM_INTEGRATOR_STEP`, `SIM_AMD_ALPHA`, `SIM_AMD_E` | exploring beyond the nearest valley |

### Langevin (default)

The Langevin integrator adds two ingredients to Newton's law. A **friction** $-\gamma m\mathbf{v}$ slows the beads down, as if they moved through a thick liquid, and a **random force** kicks them, as if molecules of the liquid were hitting them:

$$
m\,\dot{\mathbf{v}}=-\nabla E-\gamma m\,\mathbf{v}+\sqrt{2\gamma m k_BT}\;\boldsymbol{\xi}(t),
$$

where $\boldsymbol{\xi}(t)$ is random noise with zero average. Friction removes energy, the noise puts it back, and the balance is exactly the temperature $T$ (the *fluctuation-dissipation* relation: the stronger the friction, the stronger the kicks needed).

- **`SIM_FRICTION_COEFF`** is $\gamma$, in 1/ps (default 0.5). Its inverse $1/\gamma=2$ ps is the time it takes a bead to forget its velocity. A *small* $\gamma$ gives inertial, quickly moving beads with a weak thermostat; a *large* $\gamma$ gives slow, sticky motion. The static result (the average structure) does not depend on $\gamma$, only the speed of exploring does.
- **`SIM_TEMPERATURE`** is $T$.
- **`SIM_INTEGRATOR_STEP`** is $\Delta t$.

### Brownian

When friction is so large that velocities die out immediately, inertia does not matter and the beads simply *diffuse* while being pushed by the force. This is the **overdamped** limit, and is what the Brownian integrator computes:

$$
\mathbf{r}\leftarrow\mathbf{r}+\frac{\Delta t}{\gamma m}\,\mathbf{F}+\sqrt{\frac{2k_BT\,\Delta t}{\gamma m}}\;\boldsymbol{\xi}.
$$

The first term is a drift along the force, the second a random walk. Chromatin lives in a very viscous nucleus, so this is often a realistic picture. It uses the same settings as Langevin, but velocities play no role.

### Verlet

The plain Verlet scheme above, with **no thermostat**. It conserves the total energy: the temperature is whatever the starting velocities (and the energy released by the structure relaxing) give, and it can drift. It uses only `SIM_INTEGRATOR_STEP`. Choose it for tests and to check that the time step is small enough (the total energy in the plots should stay flat).

### Variable-step integrators

`variable_verlet` and `variable_langevin` choose the time step *by themselves*. After every step they estimate the integration error by comparing the result with a simpler method, and they shrink or grow $\Delta t$ to keep that error below `SIM_ERROR_TOLERANCE` (default 0.01). They have **no fixed time step**, so `SIM_INTEGRATOR_STEP` is ignored, and the simulated time per MD step varies. They are useful when you are unsure which step is safe. A smaller tolerance means a more careful, slower run.

### Accelerated MD (`amd`)

Chromatin can sit in a valley for a very long time because the hills around it are high. **Accelerated MD** lifts the valley floors, so the chain crosses the hills more easily. Whenever the energy $V$ is below a threshold $E_\text{thr}$, the force is rescaled as if the energy were

$$
V^*=V+\frac{(E_\text{thr}-V)^2}{\alpha+E_\text{thr}-V},
$$

and unchanged otherwise (OpenMM, following Hamelberg et al. 2007).

- **`SIM_AMD_E`** is the threshold $E_\text{thr}$ in kJ/mol (default 1000). Valleys deeper than this are lifted.
- **`SIM_AMD_ALPHA`** is $\alpha$ in kJ/mol (default 100). A *small* $\alpha$ flattens the landscape more strongly (a stronger boost); a *large* $\alpha$ is gentler.

The numbers must be adapted to the size of your system, because the total energy grows with the number of beads. AMD has no thermostat, and the structures it samples come from a modified landscape, so use it for exploration, not for reading off thermodynamic averages.

!!! tip "Which one should I use?"
    Keep `langevin` unless you have a reason. Use `brownian` if you want purely diffusive motion, `verlet` only for tests, and `amd` only to explore further than the nearest valley. If a run produces `NaN`, reduce `SIM_INTEGRATOR_STEP` or try a `variable_*` integrator.

## How long, and what is saved

| Key | Default | Meaning |
|---|---|---|
| `SIM_RUN_MD` | `False` | run MD after the minimisation (presets such as `GENE` switch it on) |
| `SIM_N_STEPS` | 10 000 | total number of time steps |
| `SIM_SAMPLING_STEP` | 100 | save a frame every this many steps (so 100 frames by default) |
| `TRJ_FRAMES` | none | if set, exactly this many frames are saved (it overrides `SIM_SAMPLING_STEP`) |
| `SIM_SET_INITIAL_VELOCITIES` | `True` | draw the starting velocities from the temperature, so that the system starts at the right $T$ |

The frames are written to `md_frames/`, together with a trajectory file `metadata/MultiMM_annealing.dcd`, and the last structure is `model/MultiMM_afterMD.cif`. With `SAVE_PLOTS` on, MultiMM also reports how the run went:

- the **energy of every force term** over time (`plots/energy_components.png`), which shows which force is doing the work,
- the **potential, kinetic and total energy** and the **temperature**, computed from the kinetic energy (it should hover around `SIM_TEMPERATURE`),
- the **RMSD** from the minimised structure (how far the chain wanders),
- the **RMSF** of every bead (how flexible each part is), and the distribution of velocities, which should look like a bell curve.

!!! tip "How do I know that the run is long enough?"
    Look at the energy components and the RMSD: when they stop drifting and just fluctuate around a level, the structure has relaxed. If they are still going up or down at the end, increase `SIM_N_STEPS`.
