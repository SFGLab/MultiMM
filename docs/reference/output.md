# Output

Everything is written below the folder `OUT_PATH` (default `results`).

```text
OUT_PATH/
├── model/                 final structures
│   ├── MultiMM_minimized.cif
│   ├── MultiMM_afterMD.cif            (only if SIM_RUN_MD = True)
│   └── chromosomes/                   (genome-wide runs: one file per chromosome)
├── md_frames/             saved frames of the MD trajectory (CIF files)
├── plots/                 diagnostic figures (if SAVE_PLOTS = True)
│   └── chromosomes/
└── metadata/              inputs as seen by the model, validation and logs
```

## Structures

The 3D models are **mmCIF** files in which every bead is a particle, so any molecular viewer (ChimeraX, PyMOL, VMD) can open them. `MultiMM_minimized.cif` is the energy-minimised structure and always exists. `MultiMM_afterMD.cif` is the last frame of the dynamics. The frames in `md_frames/` follow the sampling set by `SIM_SAMPLING_STEP` or `TRJ_FRAMES` ([Advanced options](../model/advanced.md#molecular-dynamics)). In a genome-wide run each chromosome is a separate chain, and also saved as its own file in `chromosomes/`. Coordinates are in the model units ([Units and scales](../model/scales.md)).

## Metadata

| File | Content |
|---|---|
| `config_auto.ini` | the full resolved configuration (defaults, config file, command line and presets) |
| `parameters.txt` | the same parameters as plain text |
| `multimm.log`, `output.log` | the log of the run |
| `MultiMM_init.cif` | the generated initial structure |
| `ms.npy`, `ns.npy`, `ds.npy` | loop anchors (bead indices) and target distances used by the model |
| `compartments.npy` | per-bead compartment states $s_i$ ([Compartments](../model/compartments.md)) |
| `compartments_from_hic.npy` | compartments derived from Hi-C PC1, if `HIC_BLOCK_COPOLYMER` was used |
| `chrom_lengths.npy`, `chrom_idxs.npy` | chromosome boundaries (in beads) and their order |
| `loop_validation.npy` | how close loop anchors ended up compared with random pairs |
| `compartment_validation.npy`, `compartment_aggregation.npy` | agreement of the structure with the input compartments, and the 3D clustering of A and B beads |
| `hic_validation.npy` | metrics comparing the structure's contact map with the input Hi-C map |
| `MultiMM_compartment_colors.cmd`, `MultiMM_chromosome_colors.cmd`, `chimera_gene_coloring.cmd` | ChimeraX scripts that colour beads by compartment, chromosome or gene |

Which of these appear depends on the inputs you used. The `.npy` files are NumPy arrays: load them with `numpy.load`.

## Plots

With `SAVE_PLOTS = True` (default) the `plots/` folder holds figures of the structure coloured by compartment or gene, the energy of each force term over the run, and the validation plots (for example the contact map of the structure next to the input map). The exact set depends on which forces and inputs were used.

## Ensembles

With `GENERATE_ENSEMBLE = True` each run gets its own subfolder `run_0`, `run_1`, ..., each with the layout above ([Ensembles](../model/advanced.md#ensembles-and-randomness)).
