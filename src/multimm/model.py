import logging
import os
import sys
import time

import numpy as np
import openmm as mm
from openmm.app import DCDReporter, ForceField, PDBxFile, Simulation, StateDataReporter
from openmm.unit import Quantity, nanometers, kelvin, nanometer, picosecond, dalton

from .initial_structure_tools import build_init_mmcif, write_cmm, write_mmcif_chrom
from .nucleosome_interpolation import NucleosomeInterpolation
from .utils import *
from .plots import *
from .read_hic import read_hic_matrix
from .hic_force import build_hic_force, auto_contact_scale
from .validation import (
    validate_hic_model, validate_hic_ensemble, validate_loops, validate_compartments,
    validate_compartment_aggregation, validate_distance_vs_strength,
)
from .logger import log_table, log_section, log_success
from .quality_tests import run_quality_tests

logger = logging.getLogger(__name__)

def _is_empty(val) -> bool:
    return val is None or str(val).strip() == "" or str(val).lower() == "none"


class MultiMM:
    def __init__(self, args):
        """
        Input data:
        ------------
        args: list of arguments imported from config.ini file.
        """
        self.md_history = {
            "step": [],
            "potential": [],
            "kinetic": [],
            "total": [],
            "temperature": [],
            "rmsd": [],          # per-sampling-step RMSD vs minimised structure (nm)
            "energy_components": {},   # per-term energy, filled in once forces are built
        }

        # Positions saved right after energy minimisation (used for MD mobility check)
        self.minimized_positions = None

        # name -> OpenMM force-group id, one per active force term — lets
        # run_md() query each term's energy separately (getState(groups={id}))
        # for the per-component energy-vs-time plot. Populated by
        # _register_force() as each add_*force() method runs.
        self.force_groups = {}
        self._next_group = 0

        # Import args
        self.args = args

        # Output folder
        self.ms, self.ns, self.ds, self.chr_ends, self.Cs = None, None, None, None, None

        # Make save directory
        self.save_path = args.OUT_PATH + "/"
        # Create main save directory and subdirectories if they don't exist
        os.makedirs(os.path.join(self.save_path, "md_frames"), exist_ok=True)
        os.makedirs(os.path.join(self.save_path, "plots"), exist_ok=True)
        if _is_empty(args.GENE_ID) and _is_empty(args.GENE_NAME) and args.LOC_START is None:
            os.makedirs(os.path.join(self.save_path, "plots", "chromosomes"), exist_ok=True)
        os.makedirs(os.path.join(self.save_path, "metadata"), exist_ok=True)
        os.makedirs(os.path.join(self.save_path, "model"), exist_ok=True)
        if _is_empty(args.GENE_ID) and _is_empty(args.GENE_NAME) and args.LOC_START is None:
            os.makedirs(os.path.join(self.save_path, "model", "chromosomes"), exist_ok=True)

        chrom = args.CHROM
        if _is_empty(chrom):
            chrom = None
        coords = [args.LOC_START, args.LOC_END] if (args.LOC_START is not None and args.LOC_END is not None) else None
        if chrom is not None and coords is None:
            if chrom in chrom_sizes:
                coords = [0, chrom_sizes[chrom]]

        if (args.GENE_TSV is not None) and (str(args.MODELLING_LEVEL).lower() == "gene"):
            if args.GENE_ID is not None and str(args.GENE_ID).lower() != "none" and str(args.GENE_ID).lower() != "":
                logger.info(f"Gene ID: {args.GENE_ID}")
                chrom, coords, gene_coords = get_gene_region(
                    gene_tsv=args.GENE_TSV,
                    gene_id=args.GENE_ID,
                    window_size=args.GENE_WINDOW,
                )
                self.gene_start, self.gene_end = ((gene_coords[0] - coords[0]) * self.args.N_BEADS) // (
                    coords[1] - coords[0]
                ), ((gene_coords[1] - coords[0]) * self.args.N_BEADS) // (coords[1] - coords[0])
                logger.info(
                    f"We model the region {coords[0]}-{coords[1]} of chrom {chrom} of the gene {args.GENE_ID}.\n"
                )
            elif (
                args.GENE_NAME is not None
                and str(args.GENE_NAME).lower() != "none"
                and str(args.GENE_NAME).lower() != ""
            ):
                logger.info(f"Gene name: {args.GENE_NAME}")
                chrom, coords, gene_coords = get_gene_region(
                    gene_tsv=args.GENE_TSV,
                    gene_name=args.GENE_NAME,
                    window_size=args.GENE_WINDOW,
                )
                self.gene_start, self.gene_end = ((gene_coords[0] - coords[0]) * self.args.N_BEADS) // (
                    coords[1] - coords[0]
                ), ((gene_coords[1] - coords[0]) * self.args.N_BEADS) // (coords[1] - coords[0])
                logger.info(
                    f"We model the region {coords[0]}-{coords[1]} of chrom {chrom} of the gene {args.GENE_NAME}.\n"
                )
            else:
                raise ValueError("You did not provide gene name or ID.")

        # Compartments
        # if args.EIGENVECTOR_TSV!=None and args.EIGENVECTOR_TSV.lower()!='none' and args.EIGENVECTOR_TSV.lower()!='':
        #     if args.EIGENVECTOR_TSV.lower().endswith('.tsv'):
        #         self.Cs, self.chr_ends = get_eigenvector(args.EIGENVECTOR_TSV, args.N_BEADS, chrom, coords)
        #     else:
        #         raise ValueError('Eigenvector should be in tsv format.')
        if args.COMPARTMENT_PATH:
            if args.COMPARTMENT_PATH.lower().endswith(".bed"):
                self.Cs, self.chr_ends, self.chrom_idxs = import_bed(
                    bed_file=args.COMPARTMENT_PATH,
                    N_beads=self.args.N_BEADS,
                    chrom=chrom,
                    coords=coords,
                    save_path=self.save_path,
                    shuffle=args.SHUFFLE_CHROMS,
                    seed=args.SHUFFLING_SEED,
                    flip_prob=args.COMPARTMENT_FLIP_PROB,
                    noise_strength=args.COMPARTMENT_NOISE_STD,
                )
            else:
                raise ValueError("Compartments file should be in .bed format.")

        # Loops (optional)
        if not _is_empty(args.LOOPS_PATH):
            if str(args.LOOPS_PATH).lower().endswith(".bedpe"):
                self.ms, self.ns, self.ds, self.chr_ends, self.chrom_idxs = import_mns_from_bedpe(
                    bedpe_file=args.LOOPS_PATH,
                    N_beads=self.args.N_BEADS,
                    coords=coords,
                    chrom=chrom,
                    path=self.save_path,
                    shuffle=args.SHUFFLE_CHROMS,
                    seed=args.SHUFFLING_SEED,
                    down_prob=args.DOWNSAMPLING_PROB,
                )
            else:
                raise ValueError("LOOPS_PATH must point to a .bedpe file.")
        else:
            if args.LE_USE_HARMONIC_BOND:
                raise ValueError(
                    "LE_USE_HARMONIC_BOND=True but no LOOPS_PATH provided. "
                    "Either supply a loops file or disable LE_USE_HARMONIC_BOND."
                )
            logger.info("No loops file provided — loop extrusion force will be skipped.")

        # Ensure chr_ends is always a numpy array (build_init_mmcif does arithmetic on it)
        if self.chr_ends is not None and not isinstance(self.chr_ends, np.ndarray):
            self.chr_ends = np.asarray(self.chr_ends, dtype=int)
            logger.debug("chr_ends converted to numpy array (was list/other sequence).")

        # Fallback chr_ends/chrom_idxs when neither loops nor compartments populated them
        if self.chr_ends is None:
            # Single-region or no-data run: treat the whole bead array as one segment
            self.chr_ends = np.array([0, self.args.N_BEADS])
            self.chrom_idxs = [0]
            logger.info(
                "chr_ends not set by loops or compartments; using single-segment fallback "
                f"[0, {self.args.N_BEADS}]."
            )

        # ── validation warnings ───────────────────────────────────────────────
        if args.HIC_USE_FORCE and _is_empty(args.HIC_PATH):
            raise ValueError(
                "HIC_USE_FORCE=True but no HIC_PATH provided. "
                "Either supply a Hi-C file or set HIC_USE_FORCE=False."
            )
        if args.COB_USE_COMPARTMENT_BLOCKS and _is_empty(args.COMPARTMENT_PATH):
            raise ValueError(
                "COB_USE_COMPARTMENT_BLOCKS=True but no COMPARTMENT_PATH provided."
            )
        if args.SCB_USE_SUBCOMPARTMENT_BLOCKS and _is_empty(args.COMPARTMENT_PATH):
            raise ValueError(
                "SCB_USE_SUBCOMPARTMENT_BLOCKS=True but no COMPARTMENT_PATH provided."
            )
        if args.IBL_USE_B_LAMINA_INTERACTION and _is_empty(args.COMPARTMENT_PATH):
            raise ValueError(
                "IBL_USE_B_LAMINA_INTERACTION=True but no COMPARTMENT_PATH provided."
            )
        if (
            args.HIC_USE_FORCE
            and not _is_empty(args.LOOPS_PATH)
            and not _is_empty(args.COMPARTMENT_PATH)
        ):
            logger.warning(
                "Hi-C force, loop extrusion, and compartment blocks are all enabled. "
                "This combination may over-constrain the simulation — consider whether you need all three."
            )
        # A .bed track is no longer the only way to get compartments: with Hi-C
        # data on hand, HIC_BLOCK_COPOLYMER derives A/B directly from the
        # matrix's own PC1 (see below) — so "no COMPARTMENT_PATH" alone no
        # longer means "no compartment data", and the warning should say so.
        hic_block_copolymer_requested = (
            args.HIC_BLOCK_COPOLYMER and args.HIC_USE_FORCE and not _is_empty(args.HIC_PATH)
        )
        if not _is_empty(args.CHROM) and _is_empty(args.COMPARTMENT_PATH):
            if hic_block_copolymer_requested:
                logger.info(
                    "Running chromosome-level simulation with no compartment .bed file, but "
                    "HIC_BLOCK_COPOLYMER=True will derive A/B compartments directly from the "
                    "Hi-C matrix's own PC1 instead (subject to the HIC_BLOCK_COPOLYMER_MIN_BP "
                    "size gate — see the warning below if the region turns out too small)."
                )
            else:
                logger.warning(
                    "Running chromosome-level simulation without compartment data. Consider "
                    "supplying COMPARTMENT_PATH, or — if you have a Hi-C matrix "
                    "(HIC_USE_FORCE=True + HIC_PATH) — enabling HIC_BLOCK_COPOLYMER to derive "
                    "compartments directly from it instead, with no .bed file needed."
                )
        # Any of loops, Hi-C force, or a compartment-block force (.bed-based or,
        # at large enough scale, Hi-C-derived) supplies a long-range restraint —
        # only warn when the region truly has none of them.
        has_long_range_restraint = (
            not _is_empty(args.LOOPS_PATH)
            or args.HIC_USE_FORCE
            or args.COB_USE_COMPARTMENT_BLOCKS
            or args.SCB_USE_SUBCOMPARTMENT_BLOCKS
        )
        if not _is_empty(args.LOC_START) and not has_long_range_restraint:
            logger.warning(
                "Running a TAD/region simulation without loops, Hi-C force, or compartment "
                "blocks. The polymer will lack long-range structural constraints."
            )

        # HIC_BLOCK_COPOLYMER: Hi-C-derived compartments (region-size gate, the
        # HIC_FORCE_OE warning, and the bed-vs-Hi-C conflict error all happen here).
        self._hic_block_copolymer_active = self._resolve_hic_block_copolymer(coords)

        # ── Hyperparameter summary ────────────────────────────────────────────
        self._log_hyperparameters()

        # ── Data Loading ──────────────────────────────────────────────────────
        log_section("Data Loading")

        self.hic_matrix = None
        self.hic_chrom = None
        if not _is_empty(args.HIC_PATH):
            hic_chrom = chrom if not _is_empty(chrom) else None
            self.hic_chrom = hic_chrom
            hic_start = args.LOC_START if args.LOC_START is not None else None
            hic_end   = args.LOC_END   if args.LOC_END   is not None else None
            logger.info(f"Loading Hi-C data from {args.HIC_PATH} …")
            try:
                region = (hic_start, hic_end) if hic_start is not None and hic_end is not None else None
                self.hic_matrix, _ = read_hic_matrix(
                    path=args.HIC_PATH,
                    chrom=hic_chrom,
                    N_beads=args.N_BEADS,
                    region=region,
                    normalization=args.HIC_NORMALIZATION,
                    max_gap=args.HIC_MAX_GAP,
                )
                nz  = int(np.count_nonzero(self.hic_matrix))
                tot = int(self.hic_matrix.size)
                log_table(
                    [
                        ("File",          args.HIC_PATH),
                        ("Shape",         str(self.hic_matrix.shape)),
                        ("Non-zero",      f"{nz} / {tot}  ({100*nz/tot:.1f}%)"),
                        ("Normalization", args.HIC_NORMALIZATION),
                        ("Max gap",       f"{args.HIC_MAX_GAP}%"),
                    ],
                    title="Hi-C matrix loaded",
                    log_fn=logger.info,
                )
                # HIC_BLOCK_COPOLYMER: only when no .bed compartments were
                # given (self.Cs is still None) — see _resolve_hic_block_copolymer.
                if self.Cs is None and self._hic_block_copolymer_active:
                    self._derive_hic_compartments()
            except Exception as exc:
                logger.error(f"Failed to load Hi-C data: {exc}")
                if args.HIC_USE_FORCE:
                    raise

        # Nucleosomes
        if args.NUC_DO_INTERPOLATION and args.ATACSEQ_PATH is not None:
            if args.ATACSEQ_PATH.lower().endswith(".bw") or args.ATACSEQ_PATH.lower().endswith(".bigwig"):
                self.atacseq = import_bw(
                    args.ATACSEQ_PATH,
                    self.args.N_BEADS,
                    chrom=self.args.CHROM,
                    coords=coords,
                    shuffle=args.SHUFFLE_CHROMS,
                    seed=args.SHUFFLING_SEED,
                )
            else:
                raise ValueError("ATAC-Seq file should be in .bw or .BigWig format.")

        if self.args.CHROM == "":
            write_chrom_colors(
                self.chr_ends,
                self.chrom_idxs,
                name=self.save_path + "metadata/MultiMM_chromosome_colors.cmd",
            )

        # Chromosomes
        self.chrom_spin, self.chrom_strength = np.zeros(self.args.N_BEADS), np.zeros(self.args.N_BEADS)
        if self.args.CHROM is None or self.args.CHROM == "" or str(self.args.CHROM).lower() == "none":
            for i in range(len(self.chr_ends) - 1):
                self.chrom_spin[self.chr_ends[i] : self.chr_ends[i + 1]] = self.chrom_idxs[i]
                self.chrom_strength[self.chr_ends[i] : self.chr_ends[i + 1]] = chrom_strength[i]

        log_success("Data Loading", logger)

    def _resolve_hic_block_copolymer(self, coords) -> bool:
        """Decide whether HIC_BLOCK_COPOLYMER is actually active this run.

        Combines the user's flag with the region-size gate and the two
        documented guardrails (see config.py / README's HIC_BLOCK_COPOLYMER).
        Called early (before Hi-C data loads) since the bed-vs-Hi-C conflict
        and region-size checks don't need the matrix itself.
        """
        args = self.args
        if not args.HIC_BLOCK_COPOLYMER:
            return False
        if not args.HIC_USE_FORCE or _is_empty(args.HIC_PATH):
            return False  # nothing to derive PC1 from

        if args.COB_USE_COMPARTMENT_BLOCKS or args.SCB_USE_SUBCOMPARTMENT_BLOCKS:
            raise ValueError(
                "HIC_BLOCK_COPOLYMER=True (the default) and a .bed-based compartment "
                "force (COB_USE_COMPARTMENT_BLOCKS / SCB_USE_SUBCOMPARTMENT_BLOCKS) are "
                "both enabled. Choose one source of compartments: set "
                "HIC_BLOCK_COPOLYMER=False to keep your .bed compartments, or disable "
                "COB_USE_COMPARTMENT_BLOCKS/SCB_USE_SUBCOMPARTMENT_BLOCKS to use the "
                "Hi-C-derived ones instead."
            )

        span_bp = (coords[1] - coords[0]) if coords is not None else None
        if span_bp is not None and span_bp < args.HIC_BLOCK_COPOLYMER_MIN_BP:
            logger.warning(
                "HIC_BLOCK_COPOLYMER=True but the modelled region (%.2f Mb) is below "
                "HIC_BLOCK_COPOLYMER_MIN_BP (%.2f Mb) — that's TAD-scale, not "
                "compartment-scale, so A/B compartments wouldn't be meaningful here. "
                "Disabling it for this run.",
                span_bp / 1e6, args.HIC_BLOCK_COPOLYMER_MIN_BP / 1e6,
            )
            return False

        if args.HIC_FORCE_OE:
            logger.warning(
                "HIC_FORCE_OE=True together with HIC_BLOCK_COPOLYMER=True: the "
                "Boltzmann force already strongly optimises O/E (compartment-level) "
                "enrichment on its own, so adding the block-copolymer force on top can "
                "over-weight compartments. Consider HIC_FORCE_OE=False when using "
                "HIC_BLOCK_COPOLYMER."
            )

        return True

    def _derive_hic_compartments(self) -> None:
        """Build self.Cs from the Hi-C matrix's own PC1 (HIC_BLOCK_COPOLYMER)
        instead of a .bed file: PC1 of the O/E matrix, sign-aligned to local
        contact density (dense -> B, sparse -> A), discretized to the same
        +1/-1/0 labels import_bed() uses — then add_compartment_blocks() (the
        existing block-copolymer force) consumes it exactly as it would bed
        data.
        """
        pc1 = hic_pc1(self.hic_matrix, already_oe=False, k=1)
        density = bead_contact_density(self.hic_matrix)
        pc1 = align_pc1_sign(pc1, density)
        self.Cs = discretize_compartments(pc1)
        np.save(self.save_path + "metadata/compartments_from_hic.npy", self.Cs)
        logger.info(
            "Compartments derived from Hi-C PC1 (HIC_BLOCK_COPOLYMER): "
            "%d A beads, %d B beads, %d unassigned.",
            int((self.Cs > 0).sum()), int((self.Cs < 0).sum()), int((self.Cs == 0).sum()),
        )

        # Soften the block-copolymer force for Hi-C-derived labels only — the
        # .bed path (COB_USE_COMPARTMENT_BLOCKS) never reaches this method, so
        # COB_EA/COB_EB stay untouched there. add_compartment_blocks() itself
        # is unchanged; we just scale the values it reads.
        scale = getattr(self.args, "HIC_BLOCK_COPOLYMER_STRENGTH_SCALE", 1.0)
        if scale != 1.0:
            ea0, eb0 = self.args.COB_EA, self.args.COB_EB
            self.args.COB_EA = ea0 * scale
            self.args.COB_EB = eb0 * scale
            logger.info(
                "HIC_BLOCK_COPOLYMER: scaling compartment-force strength by %.2f "
                "(Hi-C-derived labels are noisier than a curated .bed track) — "
                "Ea %.3g -> %.3g, Eb %.3g -> %.3g. Tune via "
                "HIC_BLOCK_COPOLYMER_STRENGTH_SCALE.",
                scale, ea0, self.args.COB_EA, eb0, self.args.COB_EB,
            )

    def _log_hyperparameters(self) -> None:
        """Print a compact summary table of key simulation hyperparameters."""
        a = self.args

        def _fmt(v) -> str:
            """Format a config value for display (strip Quantity units verbosely)."""
            try:
                from openmm.unit import Quantity as Q
                if isinstance(v, Q):
                    return str(v)
            except Exception:
                pass
            if v is None or (isinstance(v, str) and v.strip() == ""):
                return "—"
            return str(v)

        rows = [
            # ── Region ──────────────────────────────────────────────────────────
            "Region",
            ("Platform",          _fmt(a.PLATFORM)),
            ("N beads",           _fmt(a.N_BEADS)),
            ("Chromosome",        _fmt(a.CHROM) if not _is_empty(a.CHROM) else "genome-wide"),
            ("LOC start / end",   f"{a.LOC_START} – {a.LOC_END}" if a.LOC_START is not None else "—"),
            ("Modelling level",   _fmt(a.MODELLING_LEVEL) if not _is_empty(a.MODELLING_LEVEL) else "custom"),
            # ── Input data ──────────────────────────────────────────────────────
            "Input data",
            ("Loops",             _fmt(a.LOOPS_PATH)),
            ("Compartments",      _fmt(a.COMPARTMENT_PATH)),
            ("ATAC-seq",          _fmt(a.ATACSEQ_PATH)),
            ("Hi-C matrix",       _fmt(a.HIC_PATH)),
            # ── Polymer backbone ─────────────────────────────────────────────────
            "Polymer backbone",
            ("Bond r₀",           _fmt(a.POL_HARMONIC_BOND_R0)),
            ("Bond k",            _fmt(a.POL_HARMONIC_BOND_K)),
            ("Angle k",           _fmt(a.POL_HARMONIC_ANGLE_CONSTANT_K)),
            # ── Active forces ────────────────────────────────────────────────────
            "Active forces",
            ("Loop extrusion",    f"✓  k={_fmt(a.LE_HARMONIC_BOND_K)}  fixed={a.LE_FIXED_DISTANCES}"
                                   if a.LE_USE_HARMONIC_BOND else "—"),
            ("Hi-C force",        (f"✓  kernel={a.HIC_BOLTZMANN_KERNEL}  "
                                    f"boltzmann_alpha={a.HIC_BOLTZMANN_ALPHA}  "
                                    f"k_scale={a.HIC_K_SCALE}"
                                    + (f"  tol_frac={a.HIC_BOLTZMANN_TOL_FRAC}"
                                       if a.HIC_BOLTZMANN_TOL_FRAC else ""))
                                   if a.HIC_USE_FORCE else "—"),
            ("Compartment A/B",   (f"✓  Ea={a.COB_EA}  Eb={a.COB_EB}"
                                    + ("  (.bed)" if a.COB_USE_COMPARTMENT_BLOCKS else "  (Hi-C PC1)"))
                                   if (a.COB_USE_COMPARTMENT_BLOCKS or self._hic_block_copolymer_active) else "—"),
            ("Subcompartments",   "✓" if a.SCB_USE_SUBCOMPARTMENT_BLOCKS else "—"),
            ("Chr territories",   "✓" if a.CHB_USE_CHROMOSOMAL_BLOCKS  else "—"),
            ("Container",         f"✓  scale={a.SC_SCALE}" if a.SC_USE_SPHERICAL_CONTAINER else "—"),
            ("B-lamina",          f"✓  scale={a.IBL_SCALE}" if a.IBL_USE_B_LAMINA_INTERACTION else "—"),
            ("Central force",     f"✓  strength={a.CF_STRENGTH}" if a.CF_USE_CENTRAL_FORCE else "—"),
            # ── MD ───────────────────────────────────────────────────────────────
            "Molecular Dynamics",
            ("Run MD",            "yes" if a.SIM_RUN_MD else "no (energy minimisation only)"),
            ("Steps",             _fmt(a.SIM_N_STEPS)   if a.SIM_RUN_MD else "—"),
            ("Temperature",       _fmt(a.SIM_TEMPERATURE) if a.SIM_RUN_MD else "—"),
            ("Initial velocities", ("random (Boltzmann)" if a.SIM_SET_INITIAL_VELOCITIES else "zero")
                                   if a.SIM_RUN_MD else "—"),
            ("Integrator",        f"{a.SIM_INTEGRATOR_TYPE}  dt={_fmt(a.SIM_INTEGRATOR_STEP)}"
                                   if a.SIM_RUN_MD else "—"),
            # ── Ensemble ─────────────────────────────────────────────────────────
            "Ensemble",
            ("Ensemble",          f"yes  n={a.N_ENSEMBLE}" if a.GENERATE_ENSEMBLE else "no"),
            ("Output path",       _fmt(a.OUT_PATH)),
        ]

        log_table(rows, title="Simulation — hyperparameters", log_fn=logger.info)

    def _next_force_group_id(self) -> int:
        """Allocate the next unique OpenMM force-group id (0-31)."""
        gid = self._next_group
        if gid > 31:
            raise RuntimeError("Too many distinct force terms — OpenMM force groups are limited to 0-31.")
        self._next_group += 1
        return gid

    def _register_force(self, force, name: str) -> int:
        """Give *force* its own OpenMM force group and record it under *name*
        in self.force_groups, so run_md can later query this term's energy on
        its own (getState(groups={id})) and plot each enabled force's energy
        contribution separately over time — see plots.plot_energy_components.
        """
        gid = self._next_force_group_id()
        force.setForceGroup(gid)
        self.force_groups[name] = gid
        return gid

    def add_evforce(self):
        """Excluded volume force with optional soft-core formulations.

        Default: power-law repulsion (original model)

        Alternatives:
            - "soft_lj"
            - "gaussian_core"
        """
        mode = getattr(self.args, "EV_FORCE_TYPE", "powerlaw")

        sigma = self.args.LE_HARMONIC_BOND_R0
        if isinstance(sigma, Quantity):
            sigma_val = sigma.value_in_unit(nanometers)
        else:
            sigma_val = float(sigma)

        self.ev_force = mm.CustomNonbondedForce("0")
        self._register_force(self.ev_force, "Excluded volume")

        self.ev_force.addGlobalParameter("epsilon", self.args.EV_EPSILON)
        self.ev_force.addGlobalParameter("r_small", self.args.EV_R_SMALL)
        self.ev_force.addGlobalParameter("sigma", sigma_val)

        # add particles
        for _ in range(self.system.getNumParticles()):
            self.ev_force.addParticle()

        # 1. DEFAULT: power-law excluded volume (current model)
        if mode == "powerlaw":

            self.ev_force.setEnergyFunction("epsilon*(sigma/(r + r_small))^EV_POWER")
            self.ev_force.addGlobalParameter("EV_POWER", self.args.EV_POWER)
            logger.info(f"Excluded volume: power-law (EV_POWER={self.args.EV_POWER}, sigma={sigma_val:.4f} nm)")

        # 2. GAUSSIAN CORE (very soft polymer melt limit)
        elif mode == "gaussian_core":

            self.ev_force.setEnergyFunction("epsilon * exp(-r^2/(2*sigma^2))")
            logger.info(f"Excluded volume: Gaussian-core (sigma={sigma_val:.4f} nm)")

        else:
            logger.error(f"Unknown EV_FORCE_TYPE: {mode}")
            raise ValueError(f"Unknown EV_FORCE_TYPE: {mode}")

        self.system.addForce(self.ev_force)

    def add_compartment_blocks(self):
        """Compartment interaction model with multiple functional forms.

        Default: Gaussian A/B segregation (original model)

        Alternatives:
            - "yukawa"
            - "powerlaw"
            - "multi_gaussian"
        """
        mode = getattr(self.args, "COB_FORCE_TYPE", "gaussian")

        self.comp_force = mm.CustomNonbondedForce("0")
        self._register_force(self.comp_force, "Compartment blocks")

        # Shared parameters
        self.comp_force.addGlobalParameter("rc", self.r_comp)
        self.comp_force.addPerParticleParameter("s")

        for i in range(self.system.getNumParticles()):
            self.comp_force.addParticle([self.Cs[i]])

        # 1. DEFAULT: Gaussian compartment segregation
        if mode == "gaussian":

            logger.info("Using Gaussian compartment interaction model")

            self.comp_force.setEnergyFunction(
                "-E * exp(-r^2/(2*rc^2)); "
                "E = (Ea*(delta(s1-1)+delta(s1-2))*(delta(s2-1)+delta(s2-2)) + "
                "Eb*(delta(s1+1)+delta(s1+2))*(delta(s2+1)+delta(s2+2)))"
            )

            self.comp_force.addGlobalParameter("Ea", self.args.COB_EA)
            self.comp_force.addGlobalParameter("Eb", self.args.COB_EB)

            logger.info(f"Gaussian compartment parameters loaded: Ea={self.args.COB_EA}, Eb={self.args.COB_EB}")

        # 2. YUKAWA: screened compartment attraction
        elif mode == "yukawa":

            logger.info("Using Yukawa compartment interaction model")

            self.comp_force.setEnergyFunction(
                "-E * exp(-r/lambda) / r; "
                "E = (Ea*(delta(s1-1)+delta(s1-2))*(delta(s1-1)+delta(s1-2)) + "
                "Eb*(delta(s1+1)+delta(s1+2))*(delta(s1+1)+delta(s1+2)))"
            )

            self.comp_force.addGlobalParameter("lambda", self.r_comp)
            self.comp_force.addGlobalParameter("Ea", self.args.COB_EA)
            self.comp_force.addGlobalParameter("Eb", self.args.COB_EB)

            logger.info(f"Yukawa model: λ={self.r_comp}, Ea={self.args.COB_EA}, Eb={self.args.COB_EB}")

        # 3. THETA: hard contact compartment model
        elif mode == "theta":

            logger.info("Using theta (hard contact) compartment model")

            self.comp_force.setEnergyFunction(
                "-E * step(rc - r); "
                "E = (Ea*(delta(s1-1)+delta(s1-2))*(delta(s2-1)+delta(s2-2)) + "
                "Eb*(delta(s1+1)+delta(s1+2))*(delta(s2+1)+delta(s2+2)))"
            )

            self.comp_force.addGlobalParameter("Ea", self.args.COB_EA)
            self.comp_force.addGlobalParameter("Eb", self.args.COB_EB)

            logger.info(f"Theta model parameters: Ea={self.args.COB_EA}, Eb={self.args.COB_EB}")

        else:
            logger.error(f"Unknown COB_FORCE_TYPE: {mode}")
            raise ValueError(f"Unknown COB_FORCE_TYPE: {mode}")

        self.system.addForce(self.comp_force)

    def add_subcompartment_blocks(self):
        """Subcompartment interaction model with selectable functional forms.

        Default: Gaussian state-dependent attraction (original model)
        Alternatives:
            - "yukawa"
            - "powerlaw"
            - "gaussian_mixture"
        """
        mode = getattr(self.args, "SCB_FORCE_TYPE", "gaussian")

        self.scomp_force = mm.CustomNonbondedForce("0")
        self._register_force(self.scomp_force, "Subcompartment blocks")

        # Shared parameters
        self.scomp_force.addGlobalParameter("rsc", self.r_comp)
        self.scomp_force.addPerParticleParameter("s")

        for i in range(self.system.getNumParticles()):
            self.scomp_force.addParticle([self.Cs[i]])

        # 1. DEFAULT: Gaussian state-dependent interaction (original model)
        if mode == "gaussian":

            logger.info("Using Gaussian subcompartment interaction model")

            self.scomp_force.setEnergyFunction(
                "-E * exp(-r^2/(2*rsc^2)); "
                "E = Ea1*delta(s1-2)*delta(s2-2) + "
                "Ea2*delta(s1-1)*delta(s2-1) + "
                "Eb1*delta(s1+1)*delta(s2+1) + "
                "Eb2*delta(s1+2)*delta(s2+2)"
            )

            self.scomp_force.addGlobalParameter("Ea1", self.args.SCB_EA1)
            self.scomp_force.addGlobalParameter("Ea2", self.args.SCB_EA2)
            self.scomp_force.addGlobalParameter("Eb1", self.args.SCB_EB1)
            self.scomp_force.addGlobalParameter("Eb2", self.args.SCB_EB2)

            logger.info("Gaussian parameters loaded (Ea/Eb set)")

        # 2. YUKAWA: screened long-range attraction
        elif mode == "yukawa":

            logger.info("Using Yukawa (screened) interaction model")

            self.scomp_force.setEnergyFunction(
                "-E * exp(-r/lambda) / r; "
                "E = Ea1*delta(s1-2)*delta(s2-2) + "
                "Ea2*delta(s1-1)*delta(s2-1) + "
                "Eb1*delta(s1+1)*delta(s2+1) + "
                "Eb2*delta(s1+2)*delta(s2+2)"
            )

            self.scomp_force.addGlobalParameter("lambda", self.r_comp)
            self.scomp_force.addGlobalParameter("Ea1", self.args.SCB_EA1)
            self.scomp_force.addGlobalParameter("Ea2", self.args.SCB_EA2)
            self.scomp_force.addGlobalParameter("Eb1", self.args.SCB_EB1)
            self.scomp_force.addGlobalParameter("Eb2", self.args.SCB_EB2)

            logger.info(f"Yukawa screening length λ = {self.r_comp}")

        # 3. THETA / SQUARE-WELL (block copolymer contact model)
        elif mode == "theta":

            logger.info("Using theta (square-well) interaction model")

            self.scomp_force.setEnergyFunction(
                "-E * step(rc - r); "
                "E = Ea1*delta(s1-2)*delta(s2-2) + "
                "Ea2*delta(s1-1)*delta(s2-1) + "
                "Eb1*delta(s1+1)*delta(s2+1) + "
                "Eb2*delta(s1+2)*delta(s2+2)"
            )

            self.scomp_force.addGlobalParameter("rc", self.r_comp)

            self.scomp_force.addGlobalParameter("Ea1", self.args.SCB_EA1)
            self.scomp_force.addGlobalParameter("Ea2", self.args.SCB_EA2)
            self.scomp_force.addGlobalParameter("Eb1", self.args.SCB_EB1)
            self.scomp_force.addGlobalParameter("Eb2", self.args.SCB_EB2)

            logger.info(f"Theta cutoff rc = {self.r_comp}")

        else:
            logger.error(f"Unknown SCB_FORCE_TYPE: {mode}")
            raise ValueError(f"Unknown SCB_FORCE_TYPE: {mode}")

        self.system.addForce(self.scomp_force)

    def add_chromosomal_blocks(self):
        """Chromosome-level soft self-attraction for globule formation.

        Default: polynomial (original model)

        Alternatives:
            - "gaussian"
            - "saturating"
        """
        mode = getattr(self.args, "CHB_FORCE_TYPE", "polynomial")

        self.chrom_block_force = mm.CustomNonbondedForce("0")
        self._register_force(self.chrom_block_force, "Chromosomal blocks")

        # ----------------------------------------
        # shared parameters
        # ----------------------------------------
        self.chrom_block_force.addGlobalParameter("k_C", self.args.CHB_KC)
        self.chrom_block_force.addGlobalParameter("dE", self.args.CHB_DE)

        self.chrom_block_force.addPerParticleParameter("chrom")

        for i in range(self.system.getNumParticles()):
            self.chrom_block_force.addParticle([self.chrom_spin[i]])

        # 1. DEFAULT: polynomial
        if mode == "polynomial":

            logger.info("Using polynomial chromosomal self-attraction model")

            self.chrom_block_force.setEnergyFunction(
                "E*(k_C*r^4 - r^3 + r^2); "
                "E = dE*delta(chrom1-chrom2)"
            )

            logger.info(f"k_C={self.args.CHB_KC}, dE={self.args.CHB_DE}")

        # 2. GAUSSIAN SELF-ATTRACTION (globular collapse kernel)
        elif mode == "gaussian":

            logger.info("Using Gaussian chromosomal collapse kernel")

            self.chrom_block_force.setEnergyFunction(
                "-E * exp(-k_C*r^2); "
                "E = dE*delta(chrom1-chrom2)"
            )

            logger.info(f"k_C={self.args.CHB_KC}, dE={self.args.CHB_DE}")

        # 3. SATURATING SOFT-CORE ATTRACTION (stable clustering)
        elif mode == "saturating":

            logger.info("Using saturating chromosomal interaction model")

            self.chrom_block_force.setEnergyFunction(
                "-E / (1 + k_C*r^2); "
                "E = dE*delta(chrom1-chrom2)"
            )

            logger.info(f"k_C={self.args.CHB_KC}, dE={self.args.CHB_DE}")

        else:
            logger.error(f"Unknown CHB_FORCE_TYPE: {mode}")
            raise ValueError(f"Unknown CHB_FORCE_TYPE: {mode}")

        self.system.addForce(self.chrom_block_force)

    def add_spherical_container(self):
        self.container_force = mm.CustomExternalForce(
            "C*(max(0, r-R2)^2+max(0, R1-r)^2); r=sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2)"
        )
        self._register_force(self.container_force, "Spherical container")
        self.container_force.addGlobalParameter("C", defaultValue=self.args.SC_SCALE)
        self.container_force.addGlobalParameter("R1", defaultValue=self.radius1)
        self.container_force.addGlobalParameter("R2", defaultValue=self.radius2)
        self.container_force.addGlobalParameter("x0", defaultValue=self.mass_center[0])
        self.container_force.addGlobalParameter("y0", defaultValue=self.mass_center[1])
        self.container_force.addGlobalParameter("z0", defaultValue=self.mass_center[2])
        for i in range(self.system.getNumParticles()):
            self.container_force.addParticle(i, [])
        self.system.addForce(self.container_force)

    def add_Blamina_interaction(self):
        """B-compartment attraction to nuclear lamina with multiple functional
        forms.

        Default: sinusoidal shell (original model)

        Alternatives:
            - "gaussian_shell"
            - "harmonic_shell"
            - "logistic_shell"
        """
        mode = getattr(self.args, "BLAMINA_FORCE_TYPE", "sin")

        self.Blamina_force = mm.CustomExternalForce("0")
        self._register_force(self.Blamina_force, "B-lamina interaction")

        # Common parameters
        self.Blamina_force.addGlobalParameter("B", self.args.IBL_SCALE)
        self.Blamina_force.addGlobalParameter("R1", self.radius1)
        self.Blamina_force.addGlobalParameter("R2", self.radius2)

        self.Blamina_force.addGlobalParameter("x0", self.mass_center[0])
        self.Blamina_force.addGlobalParameter("y0", self.mass_center[1])
        self.Blamina_force.addGlobalParameter("z0", self.mass_center[2])

        self.Blamina_force.addPerParticleParameter("s")

        # radial distance
        r_expr = "r=sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2);"

        # 1. DEFAULT: sinusoidal shell (original)
        if mode == "sin":
                
            logger.info("Using sinusoidal lamina shell model")

            self.Blamina_force.setEnergyFunction(
                "B*(sin(pi*(r-R1)/(R2-R1))^8 - 1)*(delta(s+1)+delta(s+2)); " + r_expr
            )
            self.Blamina_force.addGlobalParameter("pi", np.pi)
            logger.info(f"Shell radii: R1={self.radius1}, R2={self.radius2}")

        # 2. GAUSSIAN SHELL (two lamina layers)
        elif mode == "gaussian_shell":

            logger.info("Using Gaussian lamina shell model (two-layer attraction)")
            self.Blamina_force.setEnergyFunction(
                "-B*(exp(-(r-R1)^2/(2*sigma^2)) + exp(-(r-R2)^2/(2*sigma^2)))"
                "*(delta(s+1)+delta(s+2)); " + r_expr
            )
            sigma = 0.1 * (self.radius2 - self.radius1)
            self.Blamina_force.addGlobalParameter("sigma", sigma)
            logger.info(f"sigma = {sigma}")

        # 3. HARMONIC SHELL (pull to mid-shell)
        elif mode == "harmonic_shell":
            logger.info("Using harmonic lamina shell model (mid-shell attraction)")
            self.Blamina_force.setEnergyFunction(
                "B*(r - r0)^2*(delta(s+1)+delta(s+2)); " + r_expr
            )
            r0 = 0.5 * (self.radius1 + self.radius2)
            self.Blamina_force.addGlobalParameter("r0", r0)
            logger.info(f"r0 (mid-shell) = {r0}")

        # 4. LOGISTIC WALLS (smooth boundary attraction)
        elif mode == "logistic_shell":
            logger.info("Using logistic lamina shell model (soft boundaries)")
            self.Blamina_force.setEnergyFunction(
                "-B*(1/(1+exp((r-R2)/lambda)) + 1/(1+exp(-(r-R1)/lambda)))"
                "*(delta(s+1)+delta(s+2)); " + r_expr
            )
            lam = 0.05 * (self.radius2 - self.radius1)
            self.Blamina_force.addGlobalParameter("lambda", lam)
            logger.info(f"lambda (boundary softness) = {lam}")

        else:
            logger.error(f"Unknown BLAMINA_FORCE_TYPE: {mode}")
            raise ValueError(f"Unknown BLAMINA_FORCE_TYPE: {mode}")

        # add particles
        for i in range(self.system.getNumParticles()):
            self.Blamina_force.addParticle(i, [self.Cs[i]])

        self.system.addForce(self.Blamina_force)

    def add_central_force(self):
        """
        Central nucleolar attraction with chromosome-size bias.
        """

        mode = getattr(self.args, "CENTRAL_FORCE_TYPE", "harmonic")

        self.central_force = mm.CustomExternalForce("0")
        self._register_force(self.central_force, "Central force")

        # --------------------------------------------
        # global parameters
        # --------------------------------------------
        self.central_force.addGlobalParameter("G", self.args.CF_STRENGTH)
        self.central_force.addGlobalParameter("R1", self.radius1)
        self.central_force.addGlobalParameter("x0", self.mass_center[0])
        self.central_force.addGlobalParameter("y0", self.mass_center[1])
        self.central_force.addGlobalParameter("z0", self.mass_center[2])

        self.central_force.addPerParticleParameter("chrom_s")

        # --------------------------------------------
        # radial distance (INLINE ONLY for OpenMM)
        # --------------------------------------------
        r = "(sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0)))"

        # ============================================================
        # 1. HARMONIC CENTERING
        # ============================================================
        if mode == "harmonic":
            logger.info("Using harmonic central attraction")

            self.central_force.setEnergyFunction(
                f"G*chrom_s*({r}-R1)*({r}-R1)"
            )

        # ============================================================
        # 2. GAUSSIAN CENTER
        # ============================================================
        elif mode == "gaussian":
            logger.info("Using Gaussian central enrichment")

            sigma = 0.5 * self.radius1
            self.central_force.addGlobalParameter("sigma", sigma)

            self.central_force.setEnergyFunction(
                f"-G*chrom_s*exp(-({r}*{r})/(2*sigma*sigma))"
            )

        # ============================================================
        # 3. LOGISTIC CORE
        # ============================================================
        elif mode == "logistic":
            logger.info("Using logistic central attraction")

            lam = 0.2 * self.radius1
            self.central_force.addGlobalParameter("lambda", lam)

            self.central_force.setEnergyFunction(
                f"-G*chrom_s*(1/(1+exp(({r}-R1)/lambda)))"
            )

        else:
            raise ValueError(f"Unknown CENTRAL_FORCE_TYPE: {mode}")

        # --------------------------------------------
        # add particles
        # --------------------------------------------
        for i in range(self.system.getNumParticles()):
            self.central_force.addParticle(i, [self.chrom_strength[i]])

        self.system.addForce(self.central_force)

    def add_harmonic_bonds(self):
        self.bond_force = mm.HarmonicBondForce()
        self._register_force(self.bond_force, "Harmonic bonds")
        for i in range(self.system.getNumParticles() - 1):
            if i not in self.chr_ends:
                self.bond_force.addBond(
                    i,
                    i + 1,
                    self.args.POL_HARMONIC_BOND_R0,
                    self.args.POL_HARMONIC_BOND_K,
                )
        self.system.addForce(self.bond_force)

    def add_loops(self):
        """
        Loop constraints using stable polymer bond models.

        Supported modes:
        - harmonic (default)
        - fene_safe (bounded FENE-like)
        - gaussian_tether (fully smooth bounded well)
        """

        mode = getattr(self.args, "LE_LOOP_FORCE_TYPE", "harmonic")

        # 1. HARMONIC (unchanged baseline)
        if mode == "harmonic":

            self.loop_force = mm.HarmonicBondForce()

            for i, (m, n) in enumerate(zip(self.ms, self.ns)):
                r0 = self.args.LE_HARMONIC_BOND_R0 if self.args.LE_FIXED_DISTANCES else self.ds[i]
                k = self.args.LE_HARMONIC_BOND_K
                self.loop_force.addBond(m, n, r0, k)

        # 2. SAFE FENE-LIKE (bounded, no singularity)
        elif mode == "fene_soft":

            self.loop_force = mm.CustomBondForce(
                "k * (r - r0)^2 / (1 + alpha * (r - r0)^2)"
            )

            self.loop_force.addPerBondParameter("r0")
            self.loop_force.addPerBondParameter("k")
            self.loop_force.addPerBondParameter("alpha")

            for i, (m, n) in enumerate(zip(self.ms, self.ns)):

                r0 = self.args.LE_HARMONIC_BOND_R0 if self.args.LE_FIXED_DISTANCES else self.ds[i]

                k = self.args.LE_HARMONIC_BOND_K
                alpha = 1.0 / (r0 ** 2)

                self.loop_force.addBond(m, n, [r0, k, alpha])

        # 3. GAUSSIAN TETHER (fully smooth bounded interaction)
        elif mode == "gaussian_tether":

            self.loop_force = mm.CustomBondForce(
                "k * (1 - exp(-(r - r0)^2 / sigma^2))"
            )

            self.loop_force.addPerBondParameter("r0")
            self.loop_force.addPerBondParameter("k")
            self.loop_force.addPerBondParameter("sigma")

            for i, (m, n) in enumerate(zip(self.ms, self.ns)):

                r0 = self.args.LE_HARMONIC_BOND_R0 if self.args.LE_FIXED_DISTANCES else self.ds[i]

                k = self.args.LE_HARMONIC_BOND_K
                sigma = r0 * 0.5

                self.loop_force.addBond(m, n, [r0, k, sigma])

        else:
            raise ValueError(f"Unknown loop force type: {mode}")

        # one shared group for all 3 branches above
        self._register_force(self.loop_force, "Loop extrusion")
        self.system.addForce(self.loop_force)

    def add_stiffness(self):
        self.angle_force = mm.HarmonicAngleForce()
        self._register_force(self.angle_force, "Harmonic angles")
        for i in range(self.system.getNumParticles() - 2):
            if (i not in self.chr_ends) and (i not in self.chr_ends - 1):
                self.angle_force.addAngle(
                    i,
                    i + 1,
                    i + 2,
                    self.args.POL_HARMONIC_ANGLE_R0,
                    self.args.POL_HARMONIC_ANGLE_CONSTANT_K,
                )
        self.system.addForce(self.angle_force)

    def initialize_simulation(self):
        if self.args.BUILD_INITIAL_STRUCTURE:
            logger.info("Creating initial structure...")

            structure_type = (
                "compartments"
                if (self.Cs is not None and len(np.unique(self.Cs)) <= 3)
                else "subcompartments"
            )
            logger.info(f"Detected structure type: {structure_type}")

            if self.Cs is not None:
                logger.info("Writing compartment color map (CMM file)")
                write_cmm(
                    self.Cs,
                    name=self.save_path + "metadata/MultiMM_compartment_colors.cmd",
                )

            logger.info("Building initial CIF structure")
            build_init_mmcif(
                n_dna=self.args.N_BEADS,
                chrom_ends=self.chr_ends,
                path=self.save_path + "metadata/",
                curve=self.args.INITIAL_STRUCTURE_TYPE,
                scale=(self.radius1 + self.radius2) / 2,
            )

            logger.info("Initial structure generated successfully")

        logger.info("Loading CIF structure into OpenMM system")

        self.pdb = (
            PDBxFile(self.save_path + "metadata/MultiMM_init.cif")
            if self.args.INITIAL_STRUCTURE_PATH is None or self.args.BUILD_INITIAL_STRUCTURE
            else PDBxFile(self.args.INITIAL_STRUCTURE_PATH)
        )

        self.mass_center = np.average(get_coordinates_mm(self.pdb.positions), axis=0)
        logger.info(f"Mass center computed: {self.mass_center}")

        logger.info("Creating OpenMM system from forcefield")
        forcefield = ForceField(self.args.FORCEFIELD_PATH)
        self.system = forcefield.createSystem(self.pdb.topology)

        match self.args.SIM_INTEGRATOR_TYPE:

            case "verlet":
                self.integrator = mm.VerletIntegrator(self.args.SIM_INTEGRATOR_STEP)

            case "variable_verlet":
                self.integrator = mm.VariableVerletIntegrator(self.SIM_ERROR_TOLERANCE)

            case "langevin":
                self.integrator = mm.LangevinIntegrator(
                    self.args.SIM_TEMPERATURE,
                    self.args.SIM_FRICTION_COEFF,
                    self.args.SIM_INTEGRATOR_STEP,
                )

            case "variable_langevin":
                self.integrator = mm.VariableLangevinIntegrator(
                    self.args.SIM_TEMPERATURE,
                    self.args.SIM_FRICTION_COEFF,
                    self.SIM_ERROR_TOLERANCE,
                )

            case "amd":
                self.integrator = mm.amd.AMDIntegrator(
                    self.args.SIM_INTEGRATOR_STEP,
                    self.args.SIM_AMD_ALPHA,
                    self.args.SIM_AMD_E,
                )

            case "brownian":
                self.integrator = mm.BrownianIntegrator(
                    self.args.SIM_TEMPERATURE,
                    self.args.SIM_FRICTION_COEFF,
                    self.args.SIM_INTEGRATOR_STEP,
                )

        logger.info(f"Integrator: {self.args.SIM_INTEGRATOR_TYPE}")

        logger.info("Simulation initialization complete")

    def add_hic_force(self):
        """Add Hi-C contact-guided force using the pre-loaded hic_matrix.

        The Hi-C force's contact-radius scale is fully independent of
        r_comp (the compartment/subcompartment block-copolymer force's
        interaction range, set in set_radiuses() from nucleus geometry) —
        the two forces model unrelated physics and previously shared one
        attribute, which was a bug. HIC_RC (config field) controls it:

          * HIC_RC is None (default): auto-calibrate from a percentile of
            the real initial pairwise-distance distribution (HIC_AUTO_SCALE,
            default True; see hic_force.auto_contact_scale), using a fixed
            50th-percentile/median (no longer a separate config knob — this
            is simply the right choice for "most pairs, not just the
            closest few, should start within the force's reach", and isn't
            meant to be tuned per-run). This differs from validation's own
            auto-scale, which uses a low percentile for visual contrast
            rather than force reach. If HIC_AUTO_SCALE is False, falls back
            to a nucleus-scale default computed independently of the
            compartment force's r_comp.
          * HIC_RC is set explicitly: used exactly as given — auto-scaling
            is skipped entirely, regardless of HIC_AUTO_SCALE.

        The resolved value is cached on self.hic_rc (never self.r_comp,
        which stays the compartment force's own value) for the downstream
        consumers that need it — get_heatmap's diagnostic plots and the
        validate_hic_model/validate_hic_ensemble calls in run().
        """
        if self.hic_matrix is None:
            logger.warning("add_hic_force() called but hic_matrix is None — skipping.")
            return

        # Excluded-volume floor distance: the closest two beads can
        # realistically get (same length scale add_evforce uses as its EV
        # "sigma"). Cached on self.hic_r_min so the later
        # validate_hic_model/validate_hic_ensemble calls can score against
        # the exact same target-distance calibration the force used.
        _r0 = self.args.LE_HARMONIC_BOND_R0
        self.hic_r_min = _r0.value_in_unit(nanometers) if isinstance(_r0, Quantity) else float(_r0)

        explicit_rc = getattr(self.args, "HIC_RC", None)
        calibrated = None

        if explicit_rc is not None:
            rc = float(explicit_rc)
            logger.info(
                "Hi-C force contact radius pinned explicitly: HIC_RC=%.4f nm "
                "(auto-calibration skipped).", rc,
            )
        else:
            # Nucleus-scale fallback, computed independently of the
            # compartment force's r_comp (same formula, own variable).
            rc = self.radius2 / 3.0
            if getattr(self.args, "HIC_AUTO_SCALE", True):
                # Fixed median (50th-percentile) calibration — no longer a
                # tunable config field (formerly HIC_AUTO_SCALE_PERCENTILE):
                # the force needs reach, so most contacted pairs, not just
                # the closest, should start within range.
                _AUTO_SCALE_PERCENTILE = 50.0
                try:
                    init_coords = get_coordinates_mm(self.pdb.positions)
                    calibrated = auto_contact_scale(
                        init_coords, percentile=_AUTO_SCALE_PERCENTILE,
                    )
                    logger.info(
                        "Hi-C force auto-calibrated from initial structure: rc=%.4f nm "
                        "(%.0fth percentile of initial pairwise distances; was %.4f nm)",
                        calibrated, _AUTO_SCALE_PERCENTILE, rc,
                    )
                    rc = calibrated
                except Exception as _e:
                    logger.warning(
                        "Hi-C force auto-calibration failed (%s) — falling back to rc=%.4f nm.",
                        _e, rc,
                    )

        self.hic_rc = rc

        boltzmann_alpha  = getattr(self.args, "HIC_BOLTZMANN_ALPHA", 4.0)
        boltzmann_kernel = getattr(self.args, "HIC_BOLTZMANN_KERNEL", "exponential")

        log_table(
            [
                ("Normalization", self.args.HIC_NORMALIZATION),
                ("k_scale",       f"{self.args.HIC_K_SCALE} kJ/mol"),
                ("kernel",        boltzmann_kernel),
                ("boltzmann_alpha", boltzmann_alpha),
                ("tol_frac (flat-bottom)", getattr(self.args, "HIC_BOLTZMANN_TOL_FRAC", 0.0)),
                ("r_min (EV floor)", f"{self.hic_r_min:.4f} nm"),
                ("rc (hic_rc)",   f"{rc:.4f} nm" + (
                    "  (explicit HIC_RC)" if explicit_rc is not None else
                    "  (auto-calibrated)" if calibrated is not None else
                    "  (nucleus-scale fallback)"
                )),
                ("OE normalize",  self.args.HIC_FORCE_OE),
                ("Matrix shape",  str(self.hic_matrix.shape)),
            ],
            title="Hi-C force — parameters",
            log_fn=logger.info,
        )
        use_noise = getattr(self.args, "HIC_NOISE_INTENSITY", 0.0) > 0
        hic_group_id = self._next_force_group_id()
        result = build_hic_force(
            H_raw=self.hic_matrix,
            N_beads=self.args.N_BEADS,
            rc=rc,
            alpha=boltzmann_alpha,
            k_scale=self.args.HIC_K_SCALE,
            oe_normalize=self.args.HIC_FORCE_OE,
            already_balanced=True,      # read_hic_matrix already normalises
            return_controller=use_noise,
            noise_seed=self.args.SHUFFLING_SEED,
            save_path=self.save_path,
            chrom=self.hic_chrom,
            r_min=self.hic_r_min,
            kernel=boltzmann_kernel,
            tol_frac=getattr(self.args, "HIC_BOLTZMANN_TOL_FRAC", 0.0),
            force_group=hic_group_id,
        )
        force, self.hic_noise = result if use_noise else (result, None)
        self.force_groups["Hi-C guided force"] = hic_group_id
        self.system.addForce(force)
        logger.info("Hi-C force added (Boltzmann-PMF, kernel=%s, alpha=%.2f).", boltzmann_kernel, boltzmann_alpha)
        if use_noise:
            logger.info(
                "Hi-C contact-strength noise enabled: intensity=%.3f, redrawn once per "
                "saved MD frame around each bond's original value.",
                self.args.HIC_NOISE_INTENSITY,
            )

    def add_forcefield(self):
        """Here we define the forcefield of MultiMM."""

        logger.info("Importing forcefield...")

        if self.args.EV_USE_EXCLUDED_VOLUME:
            self.add_evforce()

        if self.args.COB_USE_COMPARTMENT_BLOCKS or self._hic_block_copolymer_active:
            self.add_compartment_blocks()

        if self.args.SCB_USE_SUBCOMPARTMENT_BLOCKS:
            self.add_subcompartment_blocks()

        if self.args.CHB_USE_CHROMOSOMAL_BLOCKS:
            self.add_chromosomal_blocks()

        if self.args.SC_USE_SPHERICAL_CONTAINER:
            self.add_spherical_container()

        if self.args.IBL_USE_B_LAMINA_INTERACTION:
            self.add_Blamina_interaction()

        if self.args.CF_USE_CENTRAL_FORCE:
            self.add_central_force()

        if self.args.POL_USE_HARMONIC_BOND:
            self.add_harmonic_bonds()

        if self.args.LE_USE_HARMONIC_BOND and self.ms is not None:
            self.add_loops()

        if self.args.HIC_USE_FORCE:
            self.add_hic_force()

        if self.args.POL_USE_HARMONIC_ANGLE:
            self.add_stiffness()

        active = [
            ("Excluded volume",       "✓" if self.args.EV_USE_EXCLUDED_VOLUME else "–"),
            ("Harmonic bonds",        "✓" if self.args.POL_USE_HARMONIC_BOND else "–"),
            ("Harmonic angles",       "✓" if self.args.POL_USE_HARMONIC_ANGLE else "–"),
            ("Loop extrusion",        "✓" if (self.args.LE_USE_HARMONIC_BOND and self.ms is not None) else "–"),
            ("Compartment blocks",    "✓" if (self.args.COB_USE_COMPARTMENT_BLOCKS or self._hic_block_copolymer_active) else "–"),
            ("Subcompartment blocks", "✓" if self.args.SCB_USE_SUBCOMPARTMENT_BLOCKS else "–"),
            ("Chromosomal blocks",    "✓" if self.args.CHB_USE_CHROMOSOMAL_BLOCKS else "–"),
            ("Spherical container",   "✓" if self.args.SC_USE_SPHERICAL_CONTAINER else "–"),
            ("B-lamina interaction",  "✓" if self.args.IBL_USE_B_LAMINA_INTERACTION else "–"),
            ("Central force",         "✓" if self.args.CF_USE_CENTRAL_FORCE else "–"),
            ("Hi-C guided force",     "✓" if self.args.HIC_USE_FORCE else "–"),
        ]
        log_table(active, title="Forcefield — active terms", log_fn=logger.info)
        logger.info("Forcefield construction complete.")

    def min_energy(self):
        logger.info("Energy minimization...")
        # Try to use CUDA or OpenCL, fall back to CPU if not available
        try:
            platform = mm.Platform.getPlatformByName(self.args.PLATFORM)

            # Only check if user *wanted* GPU
            if self.args.PLATFORM in ["CUDA", "OpenCL"]:
                if platform.getName() not in ["CUDA", "OpenCL"]:
                    raise Exception(f"{self.args.PLATFORM} is not CUDA or OpenCL")
        except Exception as e:
            logger.info(f"Failed to find {self.args.PLATFORM}: {e}. Falling back to CPU.")
            platform = mm.Platform.getPlatformByName("CPU")
        if self.args.PLATFORM == "CPU" and self.args.CPU_THREADS is not None:
            platform.setPropertyDefaultValue("Threads", f"{self.args.CPU_THREADS}")

        # Run the simulation
        self.simulation = Simulation(self.pdb.topology, self.system, self.integrator, platform)
        self.simulation.context.setPositions(self.pdb.positions)

        if self.args.SIM_SET_INITIAL_VELOCITIES:
            # Random initial velocity field drawn from the Maxwell-Boltzmann
            # distribution at SIM_TEMPERATURE, so the MD run doesn't start
            # from an unphysical all-zero velocity state.
            self.simulation.context.setVelocitiesToTemperature(
                self.args.SIM_TEMPERATURE, self.args.SHUFFLING_SEED
            )
            logger.info(
                f"Initial velocities randomized at {self.args.SIM_TEMPERATURE} "
                f"(seed={self.args.SHUFFLING_SEED})."
            )
        else:
            logger.info("Initial velocities left at zero (SIM_SET_INITIAL_VELOCITIES=False).")

        # Report which platform is being used
        current_platform = self.simulation.context.getPlatform()
        logger.info(f"Simulation will run on platform: {current_platform.getName()}.")

        # Perform energy minimization
        start_time = time.time()
        self.simulation.minimizeEnergy()

        # Save the minimized structure
        self.state = self.simulation.context.getState(getPositions=True)
        PDBxFile.writeFile(
            self.pdb.topology,
            self.state.getPositions(),
            open(self.save_path + "model/MultiMM_minimized.cif", "w"),
        )

        # Cache minimised positions for the MD-mobility quality check.
        # getPositions(asNumpy=True) returns an OpenMM Quantity; extract nm values.
        try:
            pos_q = self.state.getPositions(asNumpy=True)
            self.minimized_positions = np.array(pos_q.value_in_unit(nanometers))  # (N, 3)
        except Exception as _e:
            logger.debug("Could not cache minimised positions: %s", _e)
            self.minimized_positions = None

        elapsed = time.time() - start_time
        logger.info(f"Energy minimization complete in {elapsed:.1f}s")

    def save_chromosomes(self):
        V = get_coordinates_mm(self.state.getPositions())
        for i in range(len(self.chr_ends) - 1):
            write_mmcif_chrom(
                coords=10 * V[self.chr_ends[i] : self.chr_ends[i + 1]],
                path=self.save_path + f"model/chromosomes/MultiMM_minimized_{chrs[self.chrom_idxs[i]]}.cif",
            )

    def run_md(self):
        # ── Determine effective sampling step ──────────────────────────────────
        # When TRJ_FRAMES is set, override SIM_SAMPLING_STEP so that exactly
        # TRJ_FRAMES CIF files are produced from SIM_N_STEPS total MD steps.
        # When TRJ_FRAMES is None, use SIM_SAMPLING_STEP directly and derive
        # the number of frames from SIM_N_STEPS // SIM_SAMPLING_STEP.
        if self.args.TRJ_FRAMES is not None:
            _sampling_step = max(1, self.args.SIM_N_STEPS // self.args.TRJ_FRAMES)
            _n_frames      = self.args.TRJ_FRAMES
            logger.info(
                "TRJ_FRAMES=%d → overriding SIM_SAMPLING_STEP to %d "
                "(SIM_N_STEPS=%d / TRJ_FRAMES=%d)",
                _n_frames, _sampling_step,
                self.args.SIM_N_STEPS, self.args.TRJ_FRAMES,
            )
        else:
            _sampling_step = self.args.SIM_SAMPLING_STEP
            _n_frames      = self.args.SIM_N_STEPS // _sampling_step

        logger.info(
            "Running relaxation — %d frames × %d steps = %d total MD steps …",
            _n_frames, _sampling_step, _n_frames * _sampling_step,
        )

        self.simulation.reporters.append(
            StateDataReporter(
                sys.stdout,
                _sampling_step,
                step=True,
                totalEnergy=True,
                kineticEnergy=True,
                potentialEnergy=True,
                temperature=True,
                separator="\t",
            )
        )
        self.simulation.reporters.append(
            DCDReporter(
                self.save_path + "metadata/MultiMM_annealing.dcd",
                _sampling_step,
            )
        )

        start = time.time()

        # ── Dynamics diagnostics accumulators ───────────────────────────────────
        # RMSF: per-bead positional fluctuation across the whole trajectory
        # (COM-removed each frame, two-pass-free online accumulation of the
        # mean vector and mean squared norm per bead — RMSF_i =
        # sqrt(mean(|x_i|^2) - |mean(x_i)|^2) over frames).
        N_beads_md = self.system.getNumParticles()
        _rmsf_sum = np.zeros((N_beads_md, 3))
        _rmsf_sumsq = np.zeros(N_beads_md)
        _rmsf_n = 0
        _last_velocities_nm_ps = None  # captured each frame; last one kept

        hic_noise = getattr(self, "hic_noise", None)

        # Degrees of freedom for T = 2*KE / (dof * kB), computed once —
        # matches OpenMM's own StateDataReporter convention: 3 per massive
        # particle, minus constraints, minus 3 more if COM motion is removed.
        _kB = 0.008314462618  # kJ/(mol·K)
        _num_massive = sum(
            1 for p in range(self.system.getNumParticles())
            if self.system.getParticleMass(p).value_in_unit(dalton) > 0
        )
        _dof = 3 * _num_massive - self.system.getNumConstraints()
        _has_cmm_remover = any(
            isinstance(self.system.getForce(f), mm.CMMotionRemover)
            for f in range(self.system.getNumForces())
        )
        if _has_cmm_remover:
            _dof -= 3
        _dof = max(1, _dof)

        # per-term energy history, one list per registered force (see
        # _register_force / plots.plot_energy_components)
        self.md_history["energy_components"] = {name: [] for name in self.force_groups}

        for i in range(_n_frames):

            self.simulation.step(_sampling_step)

            # Re-noise Hi-C contact strengths once per saved frame, anchored
            # back to the original data each time (see HiCNoiseController) —
            # a simple stochastic nudge to help the structure explore nearby
            # configurations instead of settling into one exact attractor.
            if hic_noise is not None:
                hic_noise.resample(self.simulation.context, self.args.HIC_NOISE_INTENSITY)

            state = self.simulation.context.getState(
                getPositions=True,
                getEnergy=True,
                getVelocities=True
            )

            # STEP (always safe)
            step = state.getStepCount()
            self.md_history["step"].append(step)

            # ENERGY (handle Quantity safely)
            pot = state.getPotentialEnergy()
            kin = state.getKineticEnergy()

            # convert to raw floats (kJ/mol in OpenMM usually)
            try:
                pot_val = pot.value_in_unit(pot.unit)
                kin_val = kin.value_in_unit(kin.unit)
            except Exception:
                # fallback if already float-like
                pot_val = float(pot)
                kin_val = float(kin)

            self.md_history["potential"].append(pot_val)
            self.md_history["kinetic"].append(kin_val)
            self.md_history["total"].append(pot_val + kin_val)

            # TEMPERATURE — always inferred from kinetic energy via the
            # equipartition theorem (T = 2*KE / (dof*kB)), never read from
            # the integrator's thermostat set-point: that's the target the
            # sim is held near, not a measurement of it.
            temp = (2.0 * kin_val) / (_dof * _kB)
            self.md_history["temperature"].append(temp)

            # ENERGY COMPONENTS — one potential-energy readout per active
            # force term (isolated via its own force group).
            for _name, _gid in self.force_groups.items():
                _comp_state = self.simulation.context.getState(getEnergy=True, groups={_gid})
                _comp_val = _comp_state.getPotentialEnergy().value_in_unit(pot.unit)
                self.md_history["energy_components"][_name].append(_comp_val)

            # RMSD vs minimised structure (COM-removed, in nm)
            if self.minimized_positions is not None:
                try:
                    pos_q = state.getPositions(asNumpy=True)
                    pos   = np.array(pos_q.value_in_unit(nanometers))   # (N, 3)
                    ref   = self.minimized_positions                     # (N, 3)
                    pos_c = pos - pos.mean(axis=0)
                    ref_c = ref - ref.mean(axis=0)
                    rmsd  = float(np.sqrt(np.mean(np.sum((pos_c - ref_c) ** 2, axis=1))))
                    self.md_history["rmsd"].append(rmsd)
                except Exception:
                    pass   # silently skip on rare Quantity conversion issues

            # RMSF accumulation (COM-removed positions, independent of the
            # minimised-structure reference above) — a global picture of how
            # much each bead moves over the whole trajectory, not just a
            # single-frame snapshot.
            try:
                pos_q = state.getPositions(asNumpy=True)
                pos_f = np.array(pos_q.value_in_unit(nanometers))   # (N, 3)
                pos_fc = pos_f - pos_f.mean(axis=0)
                _rmsf_sum += pos_fc
                _rmsf_sumsq += np.sum(pos_fc ** 2, axis=1)
                _rmsf_n += 1
            except Exception:
                pass

            # Per-bead velocities (nm/ps, OpenMM convention) — keep the last
            # frame as a representative equilibrium snapshot for the
            # Maxwell-Boltzmann velocity-distribution diagnostic.
            try:
                vel_q = state.getVelocities(asNumpy=True)
                _last_velocities_nm_ps = np.array(vel_q.value_in_unit(nanometer / picosecond))
            except Exception:
                pass

            # SAVE FRAME — every iteration saves one CIF; loop runs TRJ_FRAMES
            # times so exactly TRJ_FRAMES files are written.
            self.state = state
            PDBxFile.writeFile(
                self.pdb.topology,
                self.state.getPositions(),
                open(self.save_path + f"md_frames/frame_{i+1}.cif", "w"),
            )
        end = time.time()
        elapsed = end - start
        self.state = self.simulation.context.getState(getPositions=True)
        PDBxFile.writeFile(
            self.pdb.topology,
            self.state.getPositions(),
            open(self.save_path + "model/MultiMM_afterMD.cif", "w"),
        )
        try:
            target_temp = self.args.SIM_TEMPERATURE
            if hasattr(target_temp, "value_in_unit"):
                target_temp = target_temp.value_in_unit(kelvin)
            else:
                target_temp = float(target_temp)
        except Exception as _e:
            logger.debug("Could not resolve target temperature for plotting: %s", _e)
            target_temp = None

        # ── Bead-dynamics diagnostics (velocity/kinetic-energy Boltzmann fit
        # + per-bead fluctuation across the trajectory) ────────────────────────
        try:
            masses_amu = np.array([
                self.system.getParticleMass(p).value_in_unit(dalton)
                for p in range(N_beads_md)
            ])
            rmsf = None
            if _rmsf_n > 0:
                mean_vec = _rmsf_sum / _rmsf_n
                mean_sq = _rmsf_sumsq / _rmsf_n
                rmsf = np.sqrt(np.maximum(mean_sq - np.sum(mean_vec ** 2, axis=1), 0.0))
            if _last_velocities_nm_ps is not None:
                analyze_dynamics(
                    velocities=_last_velocities_nm_ps,
                    masses=masses_amu,
                    save_path=self.save_path,
                    name="dynamics",
                    target_temperature=target_temp,
                    rmsf=rmsf,
                )
        except Exception as _e:
            logger.warning(f"Dynamics diagnostics failed: {_e}")

        plot_md_thermo(
            self.md_history,
            self.save_path,
            target_temperature=target_temp,
        )
        plot_energy_components(
            self.md_history,
            self.save_path,
        )
        logger.info(f"MD finished in {elapsed:.1f}s — structure saved to {self.save_path}model/MultiMM_afterMD.cif")

    def nuc_interpolation(self):
        logger.info("Running nucleosome interpolation...")
        start = time.time()
        nuc_interpol = NucleosomeInterpolation(
            get_coordinates_cif(self.save_path + "model/MultiMM_minimized.cif"),
            self.atacseq,
            self.args.MAX_NUCS_PER_BEAD,
            self.args.NUC_RADIUS,
            self.args.POINTS_PER_NUC,
            self.args.PHI_NORM,
        )
        Vnuc = nuc_interpol.interpolate_structure_with_nucleosomes()
        write_mmcif_chrom(Vnuc, path=self.save_path + "model/MultiMM_minimized_with_nucs.cif")
        end = time.time()
        elapsed = end - start
        logger.info(f"Nucleosome interpolation complete in {elapsed:.1f}s")

    def set_radiuses(self):
        # Bead spacing, nucleus radius (constant-density globule: R ~ b0*N^(1/3)),
        # and nucleolus radius (20% inner volume fraction).
        b0 = self.args.POL_HARMONIC_BOND_R0
        if hasattr(b0, "value_in_unit"):
            b0 = b0.value_in_unit(nanometers)
        else:
            b0 = float(b0)

        N = float(self.args.N_BEADS)
        R2 = b0 * N ** (1.0 / 3.0)
        inner_volume_fraction = 0.20
        R1 = R2 * inner_volume_fraction ** (1.0 / 3.0)

        # r_comp: interaction range for compartment/subcompartment attraction
        # and the Hi-C force's fallback kernel scale. Must be nucleus-scale
        # (R2/3) rather than a fixed ~1-2 bead diameters (1.5*b0, kept below
        # as bead_contact_r for reference only) — a fixed microscale doesn't
        # grow with N like R2 does, so it left these long-range forces with
        # ~no gradient at realistic bead separations and no ability to fold
        # the structure toward the experimental Hi-C map. The Hi-C force
        # further refines this via its own auto-calibrated scale — see
        # add_hic_force().
        bead_contact_r = 1.5 * b0        # kept for reference/diagnostics only
        r_comp = R2 / 3.0

        self.radius2 = R2
        self.radius1 = R1
        self.bead_contact_r = bead_contact_r
        self.r_comp = r_comp

        log_table(
            [
                ("Bead spacing b0",      f"{b0:.4f} nm"),
                ("N beads",              f"{N:.0f}"),
                ("R nucleus",            f"{R2:.4f} nm"),
                ("R nucleolus",          f"{R1:.4f} nm"),
                ("r_comp (long-range)",  f"{r_comp:.4f} nm"),
                ("bead-contact scale",   f"{bead_contact_r:.4f} nm  (reference only)"),
            ],
            title="System geometry",
            log_fn=logger.info,
        )

    def make_plots(self):
        is_gw = (
            _is_empty(self.args.GENE_ID)
            and _is_empty(self.args.GENE_NAME)
            and self.args.LOC_START is None
            and self.args.LOC_END is None
            and self.chrom_idxs is not None
            and len(self.chrom_idxs) > 1
        )

        is_comp = self.Cs is not None and len(self.Cs) > 0

        def _viz_and_heat(cif_path, out_name, colors=None):
            """Unified structure + heatmap pipeline (single source of truth)."""
            V = get_coordinates_cif(cif_path)

            # 3D structure
            viz_structure(
                V,
                colors,
                r=0.2,
                cmap="coolwarm",
                save_path=self.save_path + f"plots/{out_name}.png",
            )

            # heatmap (always) — same Boltzmann-PMF parameters as the Hi-C
            # force itself (see add_hic_force / validate_hic_model call
            # sites below), so the structure-derived contact map uses
            # identical methodology to the force it is diagnosing.
            if self.args.N_BEADS<50000:
                get_heatmap(
                    cif_path,
                    viz=True,
                    save=True,
                    save_path=self.save_path + f"plots",
                    name=out_name,
                    rc=getattr(self, "hic_rc", self.radius2 / 3.0),
                    alpha=self.args.HIC_BOLTZMANN_ALPHA,
                    kernel=self.args.HIC_BOLTZMANN_KERNEL,
                )
            else:
                logger.warning("Heatmap skipped — system is too large for visualization (N_BEADS ≥ 50 000).")

            # structure analysis (NEW)
            analyze_structure(
                V,
                save_path=self.save_path,
                name=out_name,
            )

            plot_projection(
                    V,
                    self.Cs,
                    save_path=self.save_path,
                    name=out_name,
                )

            return V

        # GW MODE
        if is_gw:

            if is_comp:
                plot_projection(
                    get_coordinates_mm(self.state.getPositions()),
                    self.Cs,
                    save_path=self.save_path,
                    name="genomewide",
                )

            viz_chroms(self.save_path, r=0.2, comps=is_comp)

            for i in range(len(self.chr_ends) - 1):
                V = get_coordinates_cif(
                    self.save_path + f"model/chromosomes/MultiMM_minimized_{chrs[self.chrom_idxs[i]]}.cif"
                )
                viz_structure(
                    V,
                    r=0.2,
                    cmap="coolwarm",
                    save_path=self.save_path + f"plots/chromosomes/{chrs[self.chrom_idxs[i]]}_minimized_structure.png",
                )

            return

        # GENE / REGION MODE
        if hasattr(self, "gene_start"):

            save_chimera_cmd(
                self.gene_start,
                self.gene_end,
                self.args.N_BEADS,
                cmd_filename=self.save_path + "metadata/chimera_gene_coloring.cmd",
            )

            for tag, path in [
                ("initial_gene", "metadata/MultiMM_init.cif"),
                ("minimized_gene", "model/MultiMM_minimized.cif"),
            ]:
                V = get_coordinates_cif(self.save_path + path)
                viz_gene_structure(
                    V,
                    self.gene_start,
                    self.gene_end,
                    r=0.2,
                    cmap="coolwarm",
                    save_path=self.save_path + f"plots/{tag}.png",
                )

            if self.args.SIM_RUN_MD:
                V = get_coordinates_cif(self.save_path + "model/MultiMM_afterMD.cif")
                viz_gene_structure(
                    V,
                    self.gene_start,
                    self.gene_end,
                    r=0.2,
                    cmap="coolwarm",
                    save_path=self.save_path + "plots/structure_afterMD_gene_coloring.png",
                )

        # COMMON STRUCTURES (always executed)
        snapshots = [
            ("initial_structure", "metadata/MultiMM_init.cif"),
            ("minimized_structure", "model/MultiMM_minimized.cif"),
        ]

        for name, path in snapshots:
            _viz_and_heat(self.save_path + path, name)

        if is_comp:
            for name, path in snapshots:
                V = get_coordinates_cif(self.save_path + path)
                viz_structure(
                    V,
                    self.Cs[: len(V)],
                    r=0.2,
                    cmap="coolwarm",
                    save_path=self.save_path + f"plots/{name}_compartment_coloring.png",
                    legend_labels=("B (dense)", "A (sparse)"),
                )

        if self.args.SIM_RUN_MD:
            md_path = "model/MultiMM_afterMD.cif"

            _viz_and_heat(self.save_path + md_path, "structure_afterMD")

            if is_comp:
                V = get_coordinates_cif(self.save_path + md_path)
                viz_structure(
                    V,
                    self.Cs[: len(V)],
                    r=0.2,
                    cmap="coolwarm",
                    save_path=self.save_path + "plots/structure_afterMD_compartment_coloring.png",
                    legend_labels=("B (dense)", "A (sparse)"),
                )

    def run(self):
        """Energy minimization for GW model."""
        # ── Data Preprocessing ────────────────────────────────────────────────
        log_section("Data Preprocessing")
        self.set_radiuses()
        log_success("Data Preprocessing", logger)

        # ── Model Preparation ─────────────────────────────────────────────────
        log_section("Model Preparation")
        self.initialize_simulation()
        self.add_forcefield()
        log_success("Model Preparation", logger)

        # ── Energy Minimization ───────────────────────────────────────────────
        log_section("Energy Minimization")
        self.min_energy()
        if _is_empty(self.args.GENE_ID) and _is_empty(self.args.GENE_NAME) and self.args.LOC_START is None:
            self.save_chromosomes()
        log_success("Energy Minimization", logger)

        # ── Molecular Dynamics Relaxation ─────────────────────────────────────
        if self.args.SIM_RUN_MD:
            log_section("Molecular Dynamics Relaxation")
            self.run_md()
            log_success("Molecular Dynamics Relaxation", logger)

        # ── Validation ────────────────────────────────────────────────────────
        # Deliberately runs BEFORE Visualization/nucleosome interpolation
        # below: validation only reads the .cif files already written to
        # disk plus the Hi-C/loop/compartment data already in memory, so it
        # has no dependency on plotting succeeding. Ordering it first
        # guarantees the validation tables and metadata/*.npy files are
        # always produced — even if plotting later fails outright, including
        # a hard native crash from the offscreen PyVista/VTK renderer (an
        # environment-specific rendering-stack issue), which a Python
        # try/except cannot catch and would otherwise take the rest of the
        # run down with it.
        # Every gate below now explains itself when it skips — a validation
        # silently not running (e.g. HIC_USE_FORCE=True but the Hi-C matrix
        # never loaded) used to be indistinguishable from "nothing to do".
        hic_validation_ready = self.args.HIC_USE_FORCE and self.hic_matrix is not None
        if self.args.HIC_USE_FORCE and self.hic_matrix is None:
            logger.warning(
                "Hi-C validation skipped: HIC_USE_FORCE=True but no Hi-C matrix is loaded. "
                "Check HIC_PATH and the 'Data Loading' section above for a load error."
            )
        elif not self.args.HIC_USE_FORCE:
            logger.info("Hi-C validation skipped: HIC_USE_FORCE=False.")

        if hic_validation_ready:
            log_section("Validation")
            logger.info("Running Hi-C validation …")
            try:
                # Match n_rw to the actual frame count: if TRJ_FRAMES was set
                # use that; otherwise derive from SIM_N_STEPS / SIM_SAMPLING_STEP.
                _n_rw = (
                    self.args.TRJ_FRAMES
                    if self.args.TRJ_FRAMES is not None
                    else self.args.SIM_N_STEPS // self.args.SIM_SAMPLING_STEP
                )
                if self.args.SIM_RUN_MD:
                    # Ensemble validation: collect all saved MD frame CIFs.
                    import glob as _glob
                    frame_dir = self.save_path + "md_frames/"
                    frame_paths = sorted(
                        _glob.glob(frame_dir + "frame_*.cif"),
                        key=lambda p: int(p.rsplit("_", 1)[-1].split(".")[0]),
                    )
                    if frame_paths:
                        metrics = validate_hic_ensemble(
                            frame_paths, self.hic_matrix,
                            save_path=self.save_path, log=logger,
                            n_rw=_n_rw,
                            confine_radius_nm=self.radius2,
                            rc=self.hic_rc,
                            alpha=self.args.HIC_BOLTZMANN_ALPHA,
                            r_min=self.hic_r_min,
                            kernel=self.args.HIC_BOLTZMANN_KERNEL,
                            insulation_window=self.args.HIC_INSULATION_WINDOW,
                        )
                    else:
                        logger.warning("No MD frame files found — falling back to minimized structure.")
                        metrics = validate_hic_model(
                            self.save_path + "model/MultiMM_minimized.cif",
                            self.hic_matrix, save_path=self.save_path, log=logger,
                            n_rw=_n_rw,
                            confine_radius_nm=self.radius2,
                            rc=self.hic_rc,
                            alpha=self.args.HIC_BOLTZMANN_ALPHA,
                            r_min=self.hic_r_min,
                            kernel=self.args.HIC_BOLTZMANN_KERNEL,
                            insulation_window=self.args.HIC_INSULATION_WINDOW,
                        )
                else:
                    metrics = validate_hic_model(
                        self.save_path + "model/MultiMM_minimized.cif",
                        self.hic_matrix, save_path=self.save_path, log=logger,
                        n_rw=_n_rw,
                        confine_radius_nm=self.radius2,
                        rc=self.hic_rc,
                        alpha=self.args.HIC_BOLTZMANN_ALPHA,
                        r_min=self.hic_r_min,
                        kernel=self.args.HIC_BOLTZMANN_KERNEL,
                        insulation_window=self.args.HIC_INSULATION_WINDOW,
                    )
                np.save(self.save_path + "metadata/hic_validation.npy", metrics)
                logger.info("Hi-C validation metrics saved → %smetadata/hic_validation.npy", self.save_path)
                log_success("Validation", logger)
            except Exception as exc:
                logger.error(f"Hi-C validation failed: {exc}", exc_info=True)

        # ── Input-vs-output diagnostics (loops / compartments) ────────────────
        final_cif = self.save_path + (
            "model/MultiMM_afterMD.cif" if self.args.SIM_RUN_MD else "model/MultiMM_minimized.cif"
        )

        # ── Distance vs. experimental strength: does the force actually pull
        # high-strength (enriched) pairs closer and leave low-strength
        # (background/depleted) pairs alone? ──────────────────────────────────
        if hic_validation_ready:
            log_section("Distance vs. Strength Validation")
            try:
                dvs_metrics = validate_distance_vs_strength(
                    final_cif, self.hic_matrix, save_path=self.save_path, log=logger,
                    oe_normalize=self.args.HIC_FORCE_OE,
                )
                if dvs_metrics:
                    np.save(self.save_path + "metadata/distance_vs_strength.npy", dvs_metrics)
                log_success("Distance vs. Strength Validation", logger)
            except Exception as exc:
                logger.warning(f"Distance-vs-strength validation failed: {exc}")

        if self.ms is not None and self.ns is not None and len(self.ms) > 0:
            log_section("Loop Validation")
            logger.info("Checking input loop anchors (.bedpe) against the output structure …")
            try:
                loop_metrics = validate_loops(
                    final_cif, self.ms, self.ns, save_path=self.save_path, log=logger,
                )
                if loop_metrics:
                    np.save(self.save_path + "metadata/loop_validation.npy", loop_metrics)
                log_success("Loop Validation", logger)
            except Exception as exc:
                logger.warning(f"Loop validation failed: {exc}")
        else:
            logger.info("Loop validation skipped: no loop data (LOOPS_PATH not provided, or no valid anchors).")

        if self.Cs is not None and len(self.Cs) > 0:
            log_section("Compartment Validation")
            logger.info("Checking input compartment track (.bed) against the output structure …")
            try:
                comp_metrics = validate_compartments(
                    final_cif, self.Cs, save_path=self.save_path, log=logger,
                )
                if comp_metrics:
                    np.save(self.save_path + "metadata/compartment_validation.npy", comp_metrics)
                log_success("Compartment Validation", logger)
            except Exception as exc:
                logger.warning(f"Compartment validation failed: {exc}")

            # 3D spatial check (distinct from the 1D PC1-vs-track check above):
            # do same-compartment beads actually cluster together in space?
            # Runs for ANY source of self.Cs — a .bed track or Hi-C-derived
            # PC1 (HIC_BLOCK_COPOLYMER) alike.
            log_section("Compartment Aggregation Validation")
            try:
                agg_metrics = validate_compartment_aggregation(
                    final_cif, self.Cs, save_path=self.save_path, log=logger,
                )
                if agg_metrics:
                    np.save(self.save_path + "metadata/compartment_aggregation.npy", agg_metrics)
                log_success("Compartment Aggregation Validation", logger)
            except Exception as exc:
                logger.warning(f"Compartment aggregation validation failed: {exc}")
        else:
            logger.info(
                "Compartment validation skipped: no compartment data (COMPARTMENT_PATH not "
                "provided, and HIC_BLOCK_COPOLYMER did not derive any — see the warnings "
                "near the top of the log if you expected it to)."
            )

        save_args_to_txt(self.args, self.args.OUT_PATH + "/metadata/parameters.txt")

        # ── Quality Control ───────────────────────────────────────────────────
        log_section("Quality Control")
        logger.info("Running post-simulation quality tests …")
        try:
            final_coords = get_coordinates_mm(self.state.getPositions())
            run_quality_tests(
                coords=final_coords,
                args=self.args,
                md_history=self.md_history if self.args.SIM_RUN_MD else None,
                compartments=self.Cs,
                chr_ends=self.chr_ends,
                ms=self.ms,
                ns=self.ns,
                nucleus_radius_nm=getattr(self, "radius2", None),
                save_path=self.save_path,
            )
        except Exception as _qt_exc:
            logger.warning("Quality tests raised an exception and were skipped: %s", _qt_exc)
        log_success("Quality Control", logger)

        # ── Nucleosome interpolation ──────────────────────────────────────────
        # Writes a separate MultiMM_minimized_with_nucs.cif; does not touch
        # the files Validation above already read, so its ordering here (and
        # guarding) is independent of that concern — guarded anyway so a
        # failure here can't take out the final log_success below.
        if self.args.NUC_DO_INTERPOLATION and self.args.ATACSEQ_PATH is not None:
            try:
                self.nuc_interpolation()
            except Exception as exc:
                logger.warning(f"Nucleosome interpolation failed: {exc}")

        # ── Visualization ─────────────────────────────────────────────────────
        # Runs LAST and guarded: see the note above Validation — a plotting
        # failure must never prevent Validation/Quality Control from running,
        # and by this point they already have.
        if self.args.SAVE_PLOTS:
            log_section("Visualization")
            logger.info("Creating and saving diagnostic plots …")
            try:
                self.make_plots()
                log_success("Visualization", logger)
            except Exception as exc:
                logger.error(f"Visualization failed, continuing without it: {exc}", exc_info=True)

        log_success("MultiMM", logger)