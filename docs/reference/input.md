# Input data

Every input is optional in principle, but the force you want decides which ones you need ([dependencies](config-file.md#dependencies-between-keys)). Coordinates are in **base pairs** on the human reference genome (hg38 chromosome sizes are built in, chromosomes `chr1` to `chr22`).

## Loops (`LOOPS_PATH`)

A tab-separated `.bedpe` file **without a header**. Each row is one loop:

| Column | Content |
|---|---|
| 1 | chromosome of the first anchor |
| 2, 3 | start and end of the first anchor |
| 4 | chromosome of the second anchor (not used; loops are assumed intra-chromosomal) |
| 5, 6 | start and end of the second anchor |
| 7 | loop strength (a count or score; use 1 for single-cell data) |

```text
chr1	160100000	160110000	chr1	160500000	160510000	37.5
chr1	160640000	160650000	chr1	161200000	161210000	12.0
```

For a region run only loops lying completely inside `LOC_START`-`LOC_END` are used. Only the ranking of the strengths matters ([Loops](../model/loops.md)). Loops are mandatory unless `LE_USE_HARMONIC_BOND = False`.

## Compartments (`COMPARTMENT_PATH`)

Two formats are recognised by the file extension.

**CALDER-style `.bed`** (tab-separated, no header): chromosome, start, end and a label in column 4. Labels starting with `A.1`, `A.2`, `B.1`, `B.2` give the four subcompartments; plain `A` and `B` are also accepted. Rows that overlap the region are used, clipped at its edges.

```text
chr1	160000000	160250000	A.1.1
chr1	160250000	160400000	B.1.2
```

**Signal track `.bw` / `.bigWig`**, for example a Hi-C eigenvector or a compartment score. The signal is averaged per bead, centred and scaled to $[-1,1]$; *higher* values are called **B** and *lower* values **A**. Only A and B are distinguished (no subcompartments). Beads with no data stay unassigned.

Short stretches of unassigned beads between annotated ones are filled from their neighbours. The details of how labels are used are in [Compartments](../model/compartments.md).

## Hi-C map (`HIC_PATH`)

A `.hic`, `.cool` or `.mcool` file. The chromosome and region are taken from `CHROM`, `LOC_START` and `LOC_END`, and the map is balanced according to `HIC_NORMALIZATION` (`KR`, `VC`, `VC_SQRT` or `NONE`). The resolution of the file does not have to match `N_BEADS`; the matrix is resized. See [Hi-C force](../model/hic-force.md).

## ATAC-seq track (`ATACSEQ_PATH`)

A `.bw` / `.bigWig` track, used only by [nucleosome interpolation](../model/advanced.md#nucleosome-interpolation).

## Gene annotation (`GENE_TSV`)

A table of gene coordinates used to find `GENE_NAME` or `GENE_ID`. A hg38 annotation is bundled with the package and is the default.

## Initial structure (`INITIAL_STRUCTURE_PATH`)

An mmCIF (`.cif`) file with one bead per atom, used instead of the generated starting structure when `BUILD_INITIAL_STRUCTURE = False` ([Initial structure](../model/initial-structure.md)). The number of beads must equal `N_BEADS`.

## Force-field file (`FORCEFIELD_PATH`)

An OpenMM XML file that defines the bead type and its mass. The bundled file is used by default and rarely needs to be changed.
