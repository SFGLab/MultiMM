# Writing a config file

The config file is where you tell MultiMM what to do. It is just a text file with lines like `N_BEADS = 1000`, and you only write the settings you want to change. If you make a typo, MultiMM will stop and suggest the right name, so there is no harm in trying.

A config file is a plain text file in **INI format**: a list of `KEY = value` lines. You only need to write the keys you want to change; everything else keeps its default (see [Arguments](arguments.md) for all keys and defaults).

## Rules

- **One section.** Put everything under a header such as `[Main]`. The name of the section does not matter.
- **Keys are the argument names** in capitals, exactly as in [Arguments](arguments.md). An unknown key stops the run and the error suggests the closest valid name, so typos are caught early.
- **Comments** start with `;` or `#`, either on their own line or after a value.
- **Booleans** can be `True`/`False`, `yes`/`no` or `1`/`0`.
- **Numbers** are written plainly. Large integers may use underscores for readability, for example `160_000_000`.
- **Quantities** (a value with a unit) are written as a number, a space and a unit name: `0.1 nanometer`, `310 kelvin`, `1 femtosecond`, `300000.0 kilojoules_per_mole/nanometer**2`. A bare number is *not* accepted for these keys.
- **Empty or missing** means "not used". An empty value for an optional path, for example `HIC_PATH =`, is the same as not setting it.
- **Paths** can be absolute or relative to the folder where you run the command.
- **Chromosomes** are written as `chr1` or just `1`.

## A first example: a region with loops and compartments

```ini
[Main]
; --- what to model ---
CHROM     = chr1
LOC_START = 160_000_000
LOC_END   = 173_000_000
N_BEADS   = 5000                 ; resolution is (LOC_END - LOC_START) / N_BEADS

; --- data ---
LOOPS_PATH       = data/loops.bedpe
COMPARTMENT_PATH = data/compartments.bed

; --- forces: loops on, compartment force on ---
LE_USE_HARMONIC_BOND       = True
COB_USE_COMPARTMENT_BLOCKS = True

; --- output ---
OUT_PATH = results/chr1_region
```

## Driving the structure with a Hi-C map

Without loops the loop term must be switched off, and the Hi-C force switched on:

```ini
[Main]
CHROM     = chr1
LOC_START = 160_000_000
LOC_END   = 173_000_000
N_BEADS   = 5000

HIC_PATH      = data/map.hic
HIC_USE_FORCE = True
HIC_K_SCALE   = 40.0
HIC_BOLTZMANN_KERNEL = exponential

LE_USE_HARMONIC_BOND = False     ; no loop file in this run
OUT_PATH = results/chr1_hic
```

## A genome-wide model

```ini
[Main]
MODELLING_LEVEL = GW             ; preset: 200000 beads, container on, MD off
LOOPS_PATH       = data/loops.bedpe
COMPARTMENT_PATH = data/subcompartments.bed
SCB_USE_SUBCOMPARTMENT_BLOCKS = True
PLATFORM = CUDA
OUT_PATH = results/genome
```

Remember that a preset overwrites some parameters ([Advanced options](../model/advanced.md)).

## Changing physical parameters

Backbone, force-field and simulation parameters have sensible defaults and mostly do not need to be set. When you do, remember the units and the model scales ([Units and scales](../model/scales.md)):

```ini
POL_HARMONIC_BOND_R0 = 0.1 nanometer
POL_HARMONIC_ANGLE_CONSTANT_K = 200.0 kilojoules_per_mole/radian**2
EV_FORCE_TYPE = gaussian_core
LE_LOOP_FORCE_TYPE = fene_soft
SIM_TEMPERATURE = 310 kelvin
```

## Dependencies between keys

Some switches only make sense together with an input file or with each other. In particular:

| If you enable | you also need |
|---|---|
| `LE_USE_HARMONIC_BOND` (on by default) | `LOOPS_PATH` |
| `COB_USE_COMPARTMENT_BLOCKS`, `SCB_USE_SUBCOMPARTMENT_BLOCKS` | `COMPARTMENT_PATH` |
| `IBL_USE_B_LAMINA_INTERACTION` | `COMPARTMENT_PATH` |
| `HIC_USE_FORCE` | `HIC_PATH` |
| `HIC_BLOCK_COPOLYMER` | `HIC_USE_FORCE` and no compartment force |
| `MODELLING_LEVEL = GENE` | `GENE_NAME` or `GENE_ID` |
| `CHROM` | `LOC_START` and `LOC_END` (or the `CHROM` preset) |

!!! tip "Keep a record"
    Every run saves the complete, resolved configuration (defaults included) as `metadata/config_auto.ini`, so any result can be reproduced by passing that file to `-c`.
