# Running MultiMM

You run MultiMM by writing a short text file that says *what you want* (which data, how many beads, which forces) and giving it to the program. Nothing else is needed.

MultiMM is a command-line tool. After [installation](../intro/installation.md) the command `MultiMM` is available.

## The basic command

Describe the run in a [config file](config-file.md) and pass it with `-c`:

```bash
MultiMM -c config.ini
```

The result is written to the folder given by `OUT_PATH` ([output](output.md)). Running `MultiMM` with no config file uses only the defaults and will stop with an explanation, because some input (at least a loop file, see below) is needed.

## Overriding values from the command line

Every [argument](arguments.md) is also a command-line option, written in lower case with two dashes. The command line wins over the config file, which wins over the defaults:

```bash
MultiMM -c config.ini --n_beads 2000 --out_path results/test
MultiMM -c config.ini --sim_run_md True
```

The order is therefore: **defaults, then the config file, then the command line**. This is convenient for scans, for example looping over `--shuffling_seed` or `--n_beads`.

## What happens during a run

1. **Configuration.** The arguments are read and checked. A mistyped key in the config file stops the run with a "did you mean ...?" suggestion, and files that are named but do not exist are reported. If `MODELLING_LEVEL` is set, the preset overwrites some parameters and lists them ([Advanced options](../model/advanced.md)).
2. **Data loading.** Loops, compartments, the Hi-C map and the gene region are read and converted to bead coordinates.
3. **Initial structure.** A starting conformation is built, or read from `INITIAL_STRUCTURE_PATH` ([Initial structure](../model/initial-structure.md)).
4. **Force field.** The enabled terms are added; a summary table shows which ones are active.
5. **Minimisation.** The energy is minimised with OpenMM. This always happens and gives `MultiMM_minimized.cif`.
6. **Molecular dynamics** (only if `SIM_RUN_MD = True`). Frames are saved along the way.
7. **Validation and plots.** The structure is compared with the input data, and diagnostic plots are produced (`SAVE_PLOTS`).

## How much it tells you (`-v`)

By default MultiMM stays quiet: it prints the banner and the eight stages of the pipeline (input data, initial structure, forces, energy minimization, MD simulation, analysis, plotting, outputs), each with a progress bar. Add `-v` flags when you want to see more.

| Flag | What you see |
|---|---|
| *(none)* | stages and progress bars only; warnings and errors are always shown |
| `-v` | plus the main messages: what was loaded, which forces were added |
| `-vv` | plus tables of statistics for the input and output data |
| `-vvv` | everything, including debug messages and OpenMM's step-by-step energy report |

```bash
MultiMM -c config.ini -vv
```

Whatever the console level, the full log is written to `metadata/multimm.log` in your output folder, so nothing is lost. Set `NO_COLOR=1` to turn colours off.

## Checks before you start

MultiMM refuses to run, with a clear message, in these cases:

- `LE_USE_HARMONIC_BOND = True` (the default!) but no `LOOPS_PATH`. If you do not have loops, set `LE_USE_HARMONIC_BOND = False`.
- `HIC_USE_FORCE = True` without `HIC_PATH`.
- `COB_USE_COMPARTMENT_BLOCKS = True` without `COMPARTMENT_PATH`.
- A path that does not exist.

## Speed and size

The cost grows with the number of beads $N$ and, for the Hi-C force, with the number of contacted pairs. As rough guidance: a few thousand beads run comfortably on a laptop CPU; tens of thousands benefit from a GPU (`PLATFORM = CUDA` or `OpenCL`); the genome-wide preset with 200 000 beads should be run on a GPU and, by default, is only minimised, not simulated. Plots of the heat map are skipped for $N\ge50\,000$.

## Examples shipped with the code

The folder `examples/` of the repository contains ready-to-edit configurations: `config_gw.ini` (genome-wide) and `config_hic.ini` (Hi-C-driven region). Copy one, change the file paths and run it with `-c`.
