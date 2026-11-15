## MultiMM 2.1.0

### Highlight: a Hi-C force

MultiMM can now fold a structure directly from a Hi-C contact map. Set `HIC_PATH` to a `.hic`, `.cool` or `.mcool` file and `HIC_USE_FORCE = True`, and the model adds a restraint between every pair of beads that were seen in contact.

The idea is *Boltzmann inversion*. A pair that touches often should sit close, and a pair that touches rarely should sit far apart. Each contact frequency is converted into a target distance, and a spring pulls the pair toward it. The stronger the contact, the stiffer the spring. The result reproduces the diagonal decay of the map, TADs and loops, and you don't need to supply loop calls.

- **Three distance kernels** (`HIC_BOLTZMANN_KERNEL`): `exponential` (default), `power_law` and `sigmoid`, with the steepness set by `HIC_BOLTZMANN_ALPHA`.
- **Automatic scale** (`HIC_AUTO_SCALE`): the contact radius is calibrated from the initial structure, so you don't have to tune it by hand. You can also set it yourself with `HIC_RC`.
- **Normalisation and clean-up:** `KR`, `VC`, `VC_SQRT` or no normalisation, plus interpolation across short gaps of empty bins (`HIC_MAX_GAP`).
- **Raw frequency or enrichment:** by default the raw contact frequency is used; `HIC_FORCE_OE` uses enrichment above the distance-decay background instead.
- **Noise for ensembles** (`HIC_NOISE_INTENSITY`): the restraint strengths are redrawn at each saved frame, so repeated runs explore nearby conformations.
- **Compartments from Hi-C** (`HIC_BLOCK_COPOLYMER`): A/B compartments can be derived from the map itself, so no `.bed` file is needed.
- **Built-in validation:** after a run, the model's contact map is compared with the input map (correlation of contacts, P(s) decay, insulation score and others), with the metrics saved in `metadata/`.

### Other changes

- **Console output:** four verbosity levels. With no flag you see only the pipeline stages and progress bars; `-v` adds details, `-vv` adds statistics tables, and `-vvv` is fully verbose. There is also a new colour theme and banner. The full log is always saved to `metadata/multimm.log`.
- **Lamina and central forces:** the B-lamina profiles (`sin`, `gaussian_shell`, `harmonic_shell`, `logistic_shell`) were corrected so the wells sit at the nuclear walls. All three central-force forms now have their minimum at R1.
- **Gene mode:** regions around several genes can be annotated, and plots mark each gene on the structure.
- **New option:** `COB_CUTOFF` sets the range of the compartment and subcompartment forces (default 5 bond lengths, the same as before).
- **Bug fixes:**
  - the variable-step integrators used the wrong error tolerance;
  - forces sharing a parameter name caused errors when combined;
  - `GENERATE_ENSEMBLE` without `N_ENSEMBLE` now stops with a clear message;
  - input file extensions are accepted in any letter case.
- **Installation:** `hic-straw`, which needs a compiler, is now optional: `pip install "MultiMM[hic]"` for `.hic` files. `.cool` and `.mcool` need nothing extra.
- **Documentation:** a new documentation site covers the model, every force and all arguments. It also has pages on minimisation and dynamics and on running ensembles. The test suite is now fast and covers the signs and shapes of all force terms.
