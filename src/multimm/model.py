import logging
import os
import sys
import time

import numpy as np
import openmm as mm
from openmm.app import DCDReporter, ForceField, PDBxFile, Simulation, StateDataReporter
from openmm.unit import Quantity, nanometers
from tqdm import tqdm

from .initial_structure_tools import build_init_mmcif, write_cmm, write_mmcif_chrom
from .nucleosome_interpolation import NucleosomeInterpolation
from .plots import (
    plot_projection,
    save_chimera_cmd,
    viz_chroms,
    viz_gene_structure,
    viz_structure,
)
from .utils import (
    chrom_sizes,
    chrom_strength,
    chrs,
    get_coordinates_cif,
    get_coordinates_mm,
    get_gene_region,
    import_bed,
    import_bw,
    import_mns_from_bedpe,
    save_args_to_txt,
    write_chrom_colors,
)

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
                )
            else:
                raise ValueError("Compartments file should be in .bed format.")

        # Loops
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
            raise ValueError("You did not provide appropriate loop file. Loop .bedpe file is obligatory.")

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
        """
        Excluded volume force with optional soft-core formulations.

        Default: power-law repulsion (your original model)

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

        # 1. DEFAULT: power-law excluded volume (your current model)
        if mode == "powerlaw":

            self.ev_force.setEnergyFunction(
                "epsilon*(sigma/(r + r_small))^EV_POWER"
            )

            self.ev_force.addGlobalParameter("EV_POWER", self.args.EV_POWER)

        # 2. SOFT LENNARD-JONES TYPE (no divergence at r→0)
        elif mode == "soft_lj":

            self.ev_force.setEnergyFunction(
                "epsilon * (sigma^n / (r^2 + r_small^2)^(n/2))"
            )

            self.ev_force.addGlobalParameter("n", self.args.EV_POWER)

        # 3. GAUSSIAN CORE (very soft polymer melt limit)
        elif mode == "gaussian_core":

            self.ev_force.setEnergyFunction(
                "epsilon * exp(-r^2/(2*sigma^2))"
            )

        else:
            raise ValueError(f"Unknown EV_FORCE_TYPE: {mode}")

        self.system.addForce(self.ev_force)

    def add_compartment_blocks(self):
        """
        Compartment interaction model with multiple functional forms.

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

            self.comp_force.setEnergyFunction(
                "-E * exp(-r^2/(2*rc^2)); "
                "E = (Ea*(delta(s1-1)+delta(s1-2))*(delta(s2-1)+delta(s2-2)) + "
                "Eb*(delta(s1+1)+delta(s1+2))*(delta(s2+1)+delta(s2+2)))"
            )

            self.comp_force.addGlobalParameter("Ea", self.args.COB_EA)
            self.comp_force.addGlobalParameter("Eb", self.args.COB_EB)

        # 2. YUKAWA: screened compartment attraction
        elif mode == "yukawa":

            self.comp_force.setEnergyFunction(
                "-E * exp(-r/lambda) / r; "
                "E = (Ea*(delta(s1-1)+delta(s1-2))*(delta(s2-1)+delta(s2-2)) + "
                "Eb*(delta(s1+1)+delta(s1+2))*(delta(s2+1)+delta(s2+2)))"
            )

            self.comp_force.addGlobalParameter("lambda", self.r_comp)
            self.comp_force.addGlobalParameter("Ea", self.args.COB_EA)
            self.comp_force.addGlobalParameter("Eb", self.args.COB_EB)

        # 3. POWER-LAW: scale-free compartment organization
        elif mode == "powerlaw":

            self.comp_force.setEnergyFunction(
                "-E / (r^alpha + eps); "
                "E = (Ea*(delta(s1-1)+delta(s1-2))*(delta(s2-1)+delta(s2-2)) + "
                "Eb*(delta(s1+1)+delta(s1+2))*(delta(s2+1)+delta(s2+2)))"
            )

            self.comp_force.addGlobalParameter("alpha", 6.0)
            self.comp_force.addGlobalParameter("eps", 1e-3)

            self.comp_force.addGlobalParameter("Ea", self.args.COB_EA)
            self.comp_force.addGlobalParameter("Eb", self.args.COB_EB)

        # 4. MULTI-GAUSSIAN: hierarchical compartment structure
        elif mode == "multi_gaussian":

            self.comp_force.setEnergyFunction(
                "-E1*exp(-r^2/(2*s1^2)) - E2*exp(-r^2/(2*s2^2)); "
                "E1 = (Ea*(delta(s1-1)+delta(s1-2))*(delta(s2-1)+delta(s2-2)) + "
                "Eb*(delta(s1+1)+delta(s1+2))*(delta(s2+1)+delta(s2+2)))"
            )

            self.comp_force.addGlobalParameter("s1", 0.5 * self.r_comp)
            self.comp_force.addGlobalParameter("s2", 1.5 * self.r_comp)

            self.comp_force.addGlobalParameter("Ea", self.args.COB_EA)
            self.comp_force.addGlobalParameter("Eb", self.args.COB_EB)

        else:
            raise ValueError(f"Unknown COB_FORCE_TYPE: {mode}")

        self.system.addForce(self.comp_force)

    def add_subcompartment_blocks(self):
        """
        Subcompartment interaction model with selectable functional forms.

        Default: Gaussian state-dependent attraction (your original model)
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

        # 1. DEFAULT: Gaussian state-dependent interaction (your model)
        if mode == "gaussian":

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

        # 2. YUKAWA: screened long-range attraction
        elif mode == "yukawa":

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

        # 3. POWER-LAW: scale-free polymer interaction
        elif mode == "powerlaw":

            self.scomp_force.setEnergyFunction(
                "-E / (r^alpha + eps); "
                "E = Ea1*delta(s1-2)*delta(s2-2) + "
                "Ea2*delta(s1-1)*delta(s2-1) + "
                "Eb1*delta(s1+1)*delta(s2+1) + "
                "Eb2*delta(s1+2)*delta(s2+2)"
            )

            self.scomp_force.addGlobalParameter("alpha", 6.0)
            self.scomp_force.addGlobalParameter("eps", 1e-3)

            self.scomp_force.addGlobalParameter("Ea1", self.args.SCB_EA1)
            self.scomp_force.addGlobalParameter("Ea2", self.args.SCB_EA2)
            self.scomp_force.addGlobalParameter("Eb1", self.args.SCB_EB1)
            self.scomp_force.addGlobalParameter("Eb2", self.args.SCB_EB2)

        # 4. GAUSSIAN MULTI-WELL: structured chromatin contacts
        elif mode == "gaussian_mixture":

            self.scomp_force.setEnergyFunction(
                "-E1*exp(-(r-r1)^2/(2*s1^2)) - E2*exp(-(r-r2)^2/(2*s2^2)); "
                "E1 = Ea1*delta(s1-2)*delta(s2-2) + Ea2*delta(s1-1)*delta(s2-1); "
                "E2 = Eb1*delta(s1+1)*delta(s2+1) + Eb2*delta(s1+2)*delta(s2+2)"
            )

            self.scomp_force.addGlobalParameter("r1", self.r_comp * 0.8)
            self.scomp_force.addGlobalParameter("r2", self.r_comp * 1.5)
            self.scomp_force.addGlobalParameter("s1", 0.3 * self.r_comp)
            self.scomp_force.addGlobalParameter("s2", 0.5 * self.r_comp)

            self.scomp_force.addGlobalParameter("Ea1", self.args.SCB_EA1)
            self.scomp_force.addGlobalParameter("Ea2", self.args.SCB_EA2)
            self.scomp_force.addGlobalParameter("Eb1", self.args.SCB_EB1)
            self.scomp_force.addGlobalParameter("Eb2", self.args.SCB_EB2)

        else:
            raise ValueError(f"Unknown SCB_FORCE_TYPE: {mode}")

        self.system.addForce(self.scomp_force)

    def add_chromosomal_blocks(self):
        self.chrom_block_force = mm.CustomNonbondedForce("E*(k_C*r^4-r^3+r^2); E=dE*delta(chrom1-chrom2)")
        self.chrom_block_force.setForceGroup(2)
        self.chrom_block_force.addGlobalParameter("k_C", defaultValue=self.args.CHB_KC)
        self.chrom_block_force.addGlobalParameter("dE", defaultValue=self.args.CHB_DE)
        self.chrom_block_force.addPerParticleParameter("chrom")
        for i in range(self.system.getNumParticles()):
            self.chrom_block_force.addParticle([self.chrom_spin[i]])
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
        """
        B-compartment attraction to nuclear lamina with multiple functional forms.

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

        # 1. DEFAULT: sinusoidal shell (your original)
        if mode == "sin":

            self.Blamina_force.setEnergyFunction(
                "B*(sin(pi*(r-R1)/(R2-R1))^8 - 1)*(delta(s+1)+delta(s+2)); "
                + r_expr
            )
            self.Blamina_force.addGlobalParameter("pi", np.pi)

        # 2. GAUSSIAN SHELL (two lamina layers)
        elif mode == "gaussian_shell":

            self.Blamina_force.setEnergyFunction(
                "-B*(exp(-(r-R1)^2/(2*sigma^2)) + exp(-(r-R2)^2/(2*sigma^2)))"
                "*(delta(s+1)+delta(s+2)); "
                + r_expr
            )

            self.Blamina_force.addGlobalParameter("sigma", 0.1 * (self.radius2 - self.radius1))

        # 3. HARMONIC SHELL (pull to mid-shell)
        elif mode == "harmonic_shell":

            self.Blamina_force.setEnergyFunction(
                "B*(r - r0)^2*(delta(s+1)+delta(s+2)); "
                + r_expr
            )

            self.Blamina_force.addGlobalParameter("r0", 0.5 * (self.radius1 + self.radius2))

        # 4. LOGISTIC WALLS (smooth boundary attraction)
        elif mode == "logistic_shell":

            self.Blamina_force.setEnergyFunction(
                "-B*(1/(1+exp((r-R2)/lambda)) + 1/(1+exp(-(r-R1)/lambda)))"
                "*(delta(s+1)+delta(s+2)); "
                + r_expr
            )

            self.Blamina_force.addGlobalParameter("lambda", 0.05 * (self.radius2 - self.radius1))

        else:
            raise ValueError(f"Unknown BLAMINA_FORCE_TYPE: {mode}")

        # add particles
        for i in range(self.system.getNumParticles()):
            self.Blamina_force.addParticle(i, [self.Cs[i]])

        self.system.addForce(self.Blamina_force)

    def add_central_force(self):
        self.central_force = mm.CustomExternalForce(
            "G*chrom_s*(sin(r-3*R1/2)+(r-3*R1/2)^2); r=sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2)"
        )
        self.central_force.setForceGroup(2)
        self.central_force.addGlobalParameter("G", defaultValue=self.args.CF_STRENGTH)
        self.central_force.addGlobalParameter("R1", defaultValue=self.radius1)
        self.central_force.addGlobalParameter("x0", defaultValue=self.mass_center[0])
        self.central_force.addGlobalParameter("y0", defaultValue=self.mass_center[1])
        self.central_force.addGlobalParameter("z0", defaultValue=self.mass_center[2])
        self.central_force.addPerParticleParameter("chrom_s")
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
        Add loop constraints using different possible biophysical bond models.

        Supported modes:
        - "harmonic" (default original)
        - "fene" (polymer physics, soft but finite extensibility)
        - "lj_soft" (Lennard-Jones-like soft tether)
        """

        mode = getattr(self.args, "LE_LOOP_FORCE_TYPE", "harmonic")

        # -----------------------------
        # 1. Harmonic bond (original)
        # -----------------------------
        if mode == "harmonic":
            self.loop_force = mm.HarmonicBondForce()
            self.loop_force.setForceGroup(1)

            for i, (m, n) in enumerate(zip(self.ms, self.ns)):
                r0 = self.args.LE_HARMONIC_BOND_R0 if self.args.LE_FIXED_DISTANCES else self.ds[i]
                k = self.args.LE_HARMONIC_BOND_K
                self.loop_force.addBond(m, n, r0, k)

        # -----------------------------
        # 2. FENE bond (recommended for polymers)
        # -----------------------------
        elif mode == "fene":
            # U = -0.5 k R0^2 log(1 - (r/R0)^2)
            self.loop_force = mm.CustomBondForce(
                "-0.5 * k * R0^2 * log(1 - (r/R0)^2)"
            )
            self.loop_force.addPerBondParameter("R0")
            self.loop_force.addPerBondParameter("k")
            self.loop_force.setForceGroup(1)

            for i, (m, n) in enumerate(zip(self.ms, self.ns)):
                r0 = self.args.LE_HARMONIC_BOND_R0 if self.args.LE_FIXED_DISTANCES else self.ds[i]

                # interpret harmonic K as effective stiffness
                R0 = r0 * 1.5  # soft extensibility scale (can tune if needed)
                k = self.args.LE_HARMONIC_BOND_K

                self.loop_force.addBond(m, n, [R0, k])

        # -----------------------------
        # 3. Soft Lennard-Jones tether
        # -----------------------------
        elif mode == "lj_soft":
            # Soft minimum around r0 without hard constraint
            self.loop_force = mm.CustomBondForce(
                "epsilon * ((sigma/r)^12 - 2*(sigma/r)^6)"
            )
            self.loop_force.addPerBondParameter("sigma")
            self.loop_force.addPerBondParameter("epsilon")
            self.loop_force.setForceGroup(1)

            for i, (m, n) in enumerate(zip(self.ms, self.ns)):
                r0 = self.args.LE_HARMONIC_BOND_R0 if self.args.LE_FIXED_DISTANCES else self.ds[i]

                sigma = r0 / (2 ** (1/6))  # minimum at r0
                epsilon = self.args.LE_HARMONIC_BOND_K

                self.loop_force.addBond(m, n, [sigma, epsilon])

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
            logger.info("\nCreating initial structure...")
            ("compartments" if np.all(self.Cs is not None) and len(np.unique(self.Cs)) <= 3 else "subcompartments")
            if np.all(self.Cs is not None):
                write_cmm(
                    self.Cs,
                    name=self.save_path + "metadata/MultiMM_compartment_colors.cmd",
                )
            build_init_mmcif(
                n_dna=self.args.N_BEADS,
                chrom_ends=self.chr_ends,
                path=self.save_path + "metadata/",
                curve=self.args.INITIAL_STRUCTURE_TYPE,
                scale=(self.radius1 + self.radius2) / 2,
            )
            logger.info("---Done!---")
        self.pdb = (
            PDBxFile(self.save_path + "metadata/MultiMM_init.cif")
            if self.args.INITIAL_STRUCTURE_PATH is None or self.args.BUILD_INITIAL_STRUCTURE
            else PDBxFile(self.args.INITIAL_STRUCTURE_PATH)
        )
        self.mass_center = np.average(get_coordinates_mm(self.pdb.positions), axis=0)
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

    def add_forcefield(self):
        """Here we define the forcefield of MultiMM."""
        # Add forces
        logger.info("\nImporting forcefield...")
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
        if self.args.LE_USE_HARMONIC_BOND:
            self.add_loops()
        if self.args.POL_USE_HARMONIC_ANGLE:
            self.add_stiffness()

    def min_energy(self):
        logger.info("\nEnergy minimization...")
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
        logger.info(
            f"--- Energy minimization done!! Executed in {(time.time() - start_time)//3600:.0f} hours, {(time.time() - start_time)%3600//60:.0f} minutes and  {(time.time() - start_time)%60:.0f} seconds. :D ---"
        )

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
            self.state = self.simulation.context.getState(getPositions=True)
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
        logger.info(
            f"Everything is done! Simulation finished succesfully!\nMD finished in {elapsed//3600:.0f} hours, {elapsed%3600//60:.0f} minutes and  {elapsed%60:.0f} seconds. ---\n"
        )

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
        logger.info(
            f"Nucleosome interpolation finished succesfully in {elapsed//3600:.0f} hours, {elapsed%3600//60:.0f} minutes and  {elapsed%60:.0f} seconds."
        )

    def set_radiuses(self):
        raw_r1 = self.args.SC_RADIUS1
        if isinstance(raw_r1, Quantity):
            self.radius1 = raw_r1.value_in_unit(nanometers)
        else:
            self.radius1 = (self.args.N_BEADS / 50000) ** (1 / 3) if raw_r1 is None else float(raw_r1)

        raw_r2 = self.args.SC_RADIUS2
        if isinstance(raw_r2, Quantity):
            self.radius2 = raw_r2.value_in_unit(nanometers)
        else:
            self.radius2 = 3.5 * (self.args.N_BEADS / 50000) ** (1 / 3) if raw_r2 is None else float(raw_r2)

        if self.args.COB_DISTANCE is not None:
            raw_cob = self.args.COB_DISTANCE
            if isinstance(raw_cob, Quantity):
                self.r_comp = raw_cob.value_in_unit(nanometers)
            else:
                self.r_comp = float(raw_cob)
        elif self.args.SCB_DISTANCE is not None:
            raw_scb = self.args.SCB_DISTANCE
            if isinstance(raw_scb, Quantity):
                self.r_comp = raw_scb.value_in_unit(nanometers)
            else:
                self.r_comp = float(raw_scb)
        else:
            self.r_comp = (self.radius2 - self.radius1) / 20

    def make_plots(self):
        is_gw = (
            _is_empty(self.args.GENE_ID)
            and _is_empty(self.args.GENE_NAME)
            and self.args.LOC_START is None
            and self.args.LOC_END is None
        )
        is_comp = np.any(self.Cs is not None).item()
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
        else:
            if hasattr(self, "gene_start"):
                save_chimera_cmd(
                    self.gene_start,
                    self.gene_end,
                    self.args.N_BEADS,
                    cmd_filename=self.save_path + "metadata/chimera_gene_coloring.cmd",
                )
                V = get_coordinates_cif(self.save_path + "metadata/MultiMM_init.cif")
                viz_gene_structure(
                    V,
                    self.gene_start,
                    self.gene_end,
                    r=0.2,
                    cmap="coolwarm",
                    save_path=self.save_path + "plots/initial_structure_gene_coloring.png",
                )
                V = get_coordinates_cif(self.save_path + "model/MultiMM_minimized.cif")
                viz_gene_structure(
                    V,
                    self.gene_start,
                    self.gene_end,
                    r=0.2,
                    cmap="coolwarm",
                    save_path=self.save_path + "plots/minimized_structure_gene_coloring.png",
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
            V = get_coordinates_cif(self.save_path + "metadata/MultiMM_init.cif")
            viz_structure(
                V,
                r=0.2,
                cmap="coolwarm",
                save_path=self.save_path + "plots/initial_structure.png",
            )
            V = get_coordinates_cif(self.save_path + "model/MultiMM_minimized.cif")
            viz_structure(
                V,
                r=0.2,
                cmap="coolwarm",
                save_path=self.save_path + "plots/minimized_structure.png",
            )
            if is_comp:
                viz_structure(
                    V,
                    self.Cs[: len(V)],
                    cmap="coolwarm",
                    r=0.2,
                    save_path=self.save_path + "plots/minimized_structure_compartment_coloring.png",
                )
            if self.args.SIM_RUN_MD:
                V = get_coordinates_cif(self.save_path + "model/MultiMM_afterMD.cif")
                viz_structure(
                    V,
                    r=0.2,
                    cmap="coolwarm",
                    save_path=self.save_path + "plots/structure_afterMD.png",
                )
                if is_comp:
                    viz_structure(
                        V,
                        self.Cs[: len(V)],
                        cmap="coolwarm",
                        r=0.2,
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

        save_args_to_txt(self.args, self.args.OUT_PATH + "/metadata/parameters.txt")
