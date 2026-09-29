import logging
import os
import sys
import time

import numpy as np
import openmm as mm
from openmm.app import DCDReporter, ForceField, PDBxFile, Simulation, StateDataReporter
from openmm.unit import Quantity, nanometers

from .initial_structure_tools import build_init_mmcif, write_cmm, write_mmcif_chrom
from .nucleosome_interpolation import NucleosomeInterpolation
from .utils import *
from .plots import *
from .read_hic import read_hic_matrix
from .hic_force import build_hic_force
from .validation import validate_hic_model, validate_hic_ensemble
from .logger import log_table

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
        }

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
        if not _is_empty(args.CHROM) and _is_empty(args.COMPARTMENT_PATH):
            logger.warning(
                "Running chromosome-level simulation without compartment data. "
                "Consider supplying COMPARTMENT_PATH for better structural accuracy."
            )
        if (
            not _is_empty(args.LOC_START)
            and _is_empty(args.LOOPS_PATH)
            and not args.HIC_USE_FORCE
        ):
            logger.warning(
                "Running a TAD/region simulation without loops or Hi-C data. "
                "The polymer will lack long-range structural constraints."
            )

        # Hi-C data loading
        self.hic_matrix = None
        if not _is_empty(args.HIC_PATH):
            hic_chrom = chrom if not _is_empty(chrom) else None
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
        self.ev_force.setForceGroup(1)

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
        self.comp_force.setForceGroup(1)

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
        self.scomp_force.setForceGroup(1)

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
        self.chrom_block_force.setForceGroup(2)

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
        self.container_force.setForceGroup(2)
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
        self.Blamina_force.setForceGroup(2)

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
        self.central_force.setForceGroup(2)

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
        self.bond_force.setForceGroup(1)
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
            self.loop_force.setForceGroup(1)

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
            self.loop_force.setForceGroup(1)

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
            self.loop_force.setForceGroup(1)

            for i, (m, n) in enumerate(zip(self.ms, self.ns)):

                r0 = self.args.LE_HARMONIC_BOND_R0 if self.args.LE_FIXED_DISTANCES else self.ds[i]

                k = self.args.LE_HARMONIC_BOND_K
                sigma = r0 * 0.5

                self.loop_force.addBond(m, n, [r0, k, sigma])

        else:
            raise ValueError(f"Unknown loop force type: {mode}")

        self.system.addForce(self.loop_force)

    def add_stiffness(self):
        self.angle_force = mm.HarmonicAngleForce()
        self.angle_force.setForceGroup(1)
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
        """Add Hi-C contact-guided force using the pre-loaded hic_matrix."""
        if self.hic_matrix is None:
            logger.warning("add_hic_force() called but hic_matrix is None — skipping.")
            return
        log_table(
            [
                ("Mode",         self.args.HIC_FORCE_MODE),
                ("Normalization",self.args.HIC_NORMALIZATION),
                ("K (SVD rank)", self.args.HIC_N_COMPONENTS),
                ("k_scale",      f"{self.args.HIC_K_SCALE} kJ/mol"),
                ("Matrix shape", str(self.hic_matrix.shape)),
            ],
            title="Hi-C force — parameters",
            log_fn=logger.info,
        )
        force = build_hic_force(
            H_raw=self.hic_matrix,
            N_beads=self.args.N_BEADS,
            r_comp=self.r_comp,
            mode=self.args.HIC_FORCE_MODE,
            K=self.args.HIC_N_COMPONENTS,
            k_scale=self.args.HIC_K_SCALE,
            already_balanced=True,   # read_hic_matrix already normalises
        )
        self.system.addForce(force)
        logger.info("Hi-C force added.")

    def add_forcefield(self):
        """Here we define the forcefield of MultiMM."""

        logger.info("Importing forcefield...")

        if self.args.EV_USE_EXCLUDED_VOLUME:
            self.add_evforce()

        if self.args.COB_USE_COMPARTMENT_BLOCKS:
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
            ("Compartment blocks",    "✓" if self.args.COB_USE_COMPARTMENT_BLOCKS else "–"),
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
        self.simulation.context.setVelocitiesToTemperature(self.args.SIM_TEMPERATURE, self.args.SHUFFLING_SEED)

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
        self.simulation.reporters.append(
            StateDataReporter(
                sys.stdout,
                self.args.SIM_SAMPLING_STEP,
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
                self.args.SIM_N_STEPS // self.args.TRJ_FRAMES,
            )
        )
        logger.info("Running relaxation...")
        start = time.time()
        for i in range(self.args.SIM_N_STEPS // self.args.SIM_SAMPLING_STEP):

            self.simulation.step(self.args.SIM_SAMPLING_STEP)

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

            # TEMPERATURE (correct OpenMM way)
            try:
                # best case: integrator exposes temperature
                temp = self.integrator.getTemperature()
                if hasattr(temp, "value_in_unit"):
                    temp = temp.value_in_unit(kelvin)
            except Exception:
                # fallback: compute from kinetic energy
                # T = 2K / (3 N k_B)
                kB = 0.008314462618  # kJ/(mol·K)
                dof = max(1, self.system.getNumParticles() * 3)
                temp = (2.0 * kin_val) / (dof * kB)

            self.md_history["temperature"].append(temp)

            # SAVE FRAME
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
        plot_md_thermo(
            self.md_history,
            self.save_path
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
        # --------------------------------------------
        # fundamental polymer scale (bead spacing)
        # --------------------------------------------
        b0 = self.args.POL_HARMONIC_BOND_R0
        if hasattr(b0, "value_in_unit"):
            b0 = b0.value_in_unit(nanometers)
        else:
            b0 = float(b0)

        N = float(self.args.N_BEADS)

        # --------------------------------------------
        # nucleus as a dense polymer globule
        # analogy: "packed ball of spaghetti"
        #
        # constant-density assumption:
        # volume ~ N * b0^3  =>  R ~ b0 * N^(1/3)
        # --------------------------------------------
        R2 = b0 * N ** (1.0 / 3.0)

        # --------------------------------------------
        # inner compartment (nucleolus-like core)
        # analogy: "denser droplet inside the globule"
        #
        # defined by volume fraction, not geometry
        # --------------------------------------------
        inner_volume_fraction = 0.20
        R1 = R2 * inner_volume_fraction ** (1.0 / 3.0)

        # --------------------------------------------
        # interaction range (NOT geometry)
        #
        # r_comp controls how far chromatin "feels"
        # compartments / lamina attraction
        #
        # analogy: interaction fuzziness around contact
        # --------------------------------------------
        r_comp = 1.5 * b0

        self.radius2 = R2
        self.radius1 = R1
        self.r_comp = r_comp

        log_table(
            [
                ("Bead spacing b0",  f"{b0:.4f} nm"),
                ("N beads",         f"{N:.0f}"),
                ("R nucleus",       f"{R2:.4f} nm"),
                ("R nucleolus",     f"{R1:.4f} nm"),
                ("r_comp",          f"{r_comp:.4f} nm"),
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

            # heatmap (always)
            if self.args.N_BEADS<50000:
                get_heatmap(
                    cif_path,
                    viz=True,
                    save=True,
                    save_path=self.save_path + f"plots",
                    name=out_name
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
                    get_coordinates_mm(self.state.getPositions()),
                    self.Cs,
                    save_path=self.save_path,
                )

            return V

        # GW MODE
        if is_gw:

            if is_comp:
                plot_projection(
                    get_coordinates_mm(self.state.getPositions()),
                    self.Cs,
                    save_path=self.save_path,
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
                )

    def run(self):
        """Energy minimization for GW model."""
        # Estimation of parameters
        self.set_radiuses()

        # Initialize simulation
        self.initialize_simulation()

        # Import forcefield
        self.add_forcefield()

        # Run simulation / Energy minimization
        self.min_energy()
        if _is_empty(self.args.GENE_ID) and _is_empty(self.args.GENE_NAME) and self.args.LOC_START is None:
            self.save_chromosomes()

        # Run molecular dynamics
        if self.args.SIM_RUN_MD:
            self.run_md()

        # Make diagnostic plots
        if self.args.SAVE_PLOTS:
            logger.info("Creating and saving plots...")
            self.make_plots()
            logger.info("Done! :)\n")
        
        # Run nucleosome interpolation
        if self.args.NUC_DO_INTERPOLATION and self.args.ATACSEQ_PATH is not None:
            self.nuc_interpolation()

        # Hi-C validation — diagonal decay, insulation score, PC1 correlation
        if self.args.HIC_USE_FORCE and self.hic_matrix is not None:
            logger.info("Running Hi-C validation …")
            if self.args.SIM_RUN_MD:
                # Ensemble validation: use all saved MD frames
                n_frames = self.args.SIM_N_STEPS // self.args.SIM_SAMPLING_STEP
                frame_paths = [
                    self.save_path + f"md_frames/frame_{i+1}.cif"
                    for i in range(n_frames)
                    if os.path.isfile(self.save_path + f"md_frames/frame_{i+1}.cif")
                ]
                if frame_paths:
                    metrics = validate_hic_ensemble(frame_paths, self.hic_matrix, save_path=self.save_path, log=logger)
                else:
                    logger.warning("No MD frame files found — falling back to minimized structure.")
                    metrics = validate_hic_model(
                        self.save_path + "model/MultiMM_minimized.cif",
                        self.hic_matrix, save_path=self.save_path, log=logger,
                    )
            else:
                metrics = validate_hic_model(
                    self.save_path + "model/MultiMM_minimized.cif",
                    self.hic_matrix, save_path=self.save_path, log=logger,
                )
            np.save(self.save_path + "metadata/hic_validation.npy", metrics)
            log_table(
                [
                    ("Diagonal decay r",   f"{metrics['diag_decay_r']:.4f}  (p={metrics['diag_decay_p']:.2e})"),
                    ("Insulation score r", f"{metrics['insulation_r']:.4f}  (p={metrics['insulation_p']:.2e})"),
                    ("|PC1| r",            f"{metrics['pc1_r']:.4f}  (p={metrics['pc1_p']:.2e})"),
                    ("Saved to",           self.save_path + "metadata/hic_validation.npy"),
                ],
                title="Hi-C validation summary",
                log_fn=logger.info,
            )

        save_args_to_txt(self.args, self.args.OUT_PATH + "/metadata/parameters.txt")

        print("\033[1;32m✅ MultiMM ran successfully!!!\033[0m\n")