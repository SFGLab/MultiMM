import importlib.resources as pkg_resources
import logging
from typing import Any, Optional

from openmm.unit import Quantity
from pydantic import BaseModel, BeforeValidator, Field, model_validator
from typing_extensions import Annotated

from .enums import InitialStructureType

try:
    default_xml_path = str(pkg_resources.files("multimm.forcefields").joinpath("ff.xml"))
except Exception:
    default_xml_path = "src/multimm/forcefields/ff.xml"

try:
    default_gene_path = str(pkg_resources.files("multimm.data").joinpath("hg38_gtf_annotations.tsv"))
except Exception:
    default_gene_path = "src/multimm/data/hg38_gtf_annotations.tsv"

logger = logging.getLogger(__name__)


def parse_quantity(val: Any) -> Quantity:
    if isinstance(val, Quantity):
        return val
    if not isinstance(val, str) or val.strip() == "":
        raise ValueError("Invalid Quantity format")
    parts = val.strip().split(maxsplit=1)
    if len(parts) != 2:
        raise ValueError(f"Can't recognise Quantity format: {val}")
    value_str, unit_str = parts
    try:
        value = float(value_str)
    except ValueError:
        raise ValueError(f"Invalid float value: {value_str}")
    import openmm.unit as u

    safe_dict = {name: getattr(u, name) for name in dir(u) if not name.startswith("__")}
    try:
        unit_obj = eval(unit_str, {"__builtins__": None}, safe_dict)
    except Exception as e:
        raise ValueError(f"Can't recognise unit expression {unit_str} in {val}: {e}")
    if not isinstance(unit_obj, (u.Unit, u.BaseUnit)):
        if hasattr(unit_obj, "unit"):
            unit_obj = unit_obj.unit
        else:
            raise ValueError(f"Expression {unit_str} did not evaluate to a Unit")
    return Quantity(value=value, unit=unit_obj)


def validate_quantity(v: Any) -> Quantity:
    if isinstance(v, Quantity):
        return v
    if isinstance(v, str):
        return parse_quantity(v)
    raise ValueError(f"Cannot cast {type(v)} to Quantity")


OpenMMQuantity = Annotated[Quantity, BeforeValidator(validate_quantity)]


def validate_boolean(v: Any) -> bool:
    if isinstance(v, bool):
        return v
    if isinstance(v, (int, float)):
        return bool(v)
    if isinstance(v, str):
        val_lower = v.strip().lower()
        if val_lower in ("true", "1", "y", "yes"):
            return True
        if val_lower in ("false", "0", "n", "no", "", "none"):
            return False
    raise ValueError(f"Cannot cast {v} to boolean")


Boolean = Annotated[bool, BeforeValidator(validate_boolean)]


def validate_chrom(v: Any) -> Optional[str]:
    if v is None:
        return None
    v_str = str(v).strip()
    if not v_str or v_str.lower() == "none":
        return None
    if not v_str.startswith("chr"):
        return f"chr{v_str}"
    return v_str


ChromStr = Annotated[Optional[str], BeforeValidator(validate_chrom)]


class SimulationConfig(BaseModel):
    model_config = {
        "arbitrary_types_allowed": True,
        "populate_by_name": True,
        "validate_assignment": True,
        "validate_default": True,
        # Reject unrecognized keyword arguments outright (defense-in-depth —
        # run.py's get_config() already checks config.ini keys specifically,
        # with a clearer error naming the offending file).
        "extra": "forbid",
    }

    @model_validator(mode="after")
    def warn_hic_k_scale(self) -> "SimulationConfig":
        """Emit a runtime warning when HIC_K_SCALE is outside the recommended range."""
        k = self.HIC_K_SCALE
        if k > 200.0:
            logger.warning(
                "HIC_K_SCALE=%.1f is above the recommended maximum of 200 kJ/mol. "
                "Very high values can over-constrain the polymer, reduce conformational "
                "diversity, and cause MD instability. Consider reducing to 20–160.",
                k,
            )
        return self

    @model_validator(mode="after")
    def warn_hic_force_misuse(self) -> "SimulationConfig":
        """Warn when a Hi-C-specific field is set away from its default while
        HIC_USE_FORCE=False — the whole Hi-C force is disabled, so nothing
        under HIC_* takes effect at all. Never a hard error, since it doesn't
        make the simulation incorrect — only some configured value ends up
        unused.
        """
        hic_fields_nondefault = (
            self.HIC_BOLTZMANN_ALPHA != 4.0 or self.HIC_BOLTZMANN_KERNEL != "exponential"
        )

        if not self.HIC_USE_FORCE and hic_fields_nondefault:
            logger.warning(
                "HIC_USE_FORCE=False, but HIC_BOLTZMANN_ALPHA=%s / "
                "HIC_BOLTZMANN_KERNEL=%s (set away from "
                "default) will have no effect — the Hi-C force isn't being built at "
                "all. Set HIC_USE_FORCE=True to actually apply it, or leave these at "
                "their defaults if you don't intend to use the Hi-C force.",
                self.HIC_BOLTZMANN_ALPHA, self.HIC_BOLTZMANN_KERNEL,
            )

        return self

    @model_validator(mode="before")
    @classmethod
    def clean_fields(cls, data: Any) -> Any:
        if isinstance(data, dict):
            cleaned = {}
            for k, v in data.items():
                if isinstance(v, str):
                    v_stripped = v.strip()
                    if v_stripped == "" or v_stripped.lower() == "none":
                        if k == "LOOPS_PATH":
                            cleaned[k] = None
                            continue
                        field = cls.model_fields.get(k)
                        if field:
                            annotation = field.annotation
                            args_types = getattr(annotation, "__args__", [])
                            if type(None) in args_types or annotation is Any:
                                cleaned[k] = None
                                continue
                        cleaned[k] = ""
                        continue
                cleaned[k] = v
            return cleaned
        return data

    PLATFORM: str = Field(default="CPU", description="name of the platform. Available choices: Reference CPU")
    CPU_THREADS: Optional[int] = Field(
        default=None, description="The number of CPU threads (in case you would like to specify them)."
    )
    DEVICE: str = Field(default="", description="device index for CUDA or OpenCL (count from 0)")
    MODELLING_LEVEL: str = Field(
        default="",
        description="Choose 'GENE' or 'REGION' for gene or TAD level, 'CHROM' for chromosome leve, and 'GW' for genome level. It will setup some parameters for you and print you helpful comments.",
    )
    INITIAL_STRUCTURE_PATH: str = Field(default="", description="Path to CIF file.")
    BUILD_INITIAL_STRUCTURE: Boolean = Field(default=True, description="To build a new initial structure.")
    INITIAL_STRUCTURE_TYPE: InitialStructureType = Field(
        default=InitialStructureType.HILBERT,
        description="you can choose between: hilbert, circle, rw, confined_rw, knot, self_avoiding_rw, spiral, sphere.",
    )
    GENERATE_ENSEMBLE: Boolean = Field(
        default=False,
        description="Default value: false. True in case that you would like to have an ensemble of structures instead of one. Better to disable it for large simulations that require long computational time. Moreover it is better to start random walk initial structure in case of true value.",
    )
    COMPARTMENT_FLIP_PROB: float = Field(
        default=0.0,
        description="Probability of flipping compartment identity per bead (A↔B). Applied after BED parsing."
    )
    COMPARTMENT_NOISE_STD: float = Field(
        default=0.0,
        description="Standard deviation of Gaussian noise applied to compartment field before discretization."
    )
    N_ENSEMBLE: Optional[int] = Field(
        default=None, description="Number of samples of structures that you would like to calculate."
    )
    DOWNSAMPLING_PROB: float = Field(default=1.0, description="Probability of downsampling contacts (from 0 to 1).")
    FORCEFIELD_PATH: str = Field(
        default=default_xml_path,
        description="Path to XML file with forcefield.",
    )
    N_BEADS: int = Field(default=50000, description="Number of Simulation Beads.")
    COMPARTMENT_PATH: Optional[str] = Field(
        default=None,
        description="A .bed file with (sub)compartments from Calder (or the same format), or a .bw/.bigwig "
        "signal track (e.g. an eigenvector/PC1 track) from which A/B compartment calls are derived.",
    )
    LOOPS_PATH: Optional[str] = Field(default=None, description="A .bedpe file path with loops. Optional — simulation runs without loops if not provided.")
    HIC_PATH: Optional[str] = Field(
        default=None,
        description="Path to a Hi-C contact file (.hic, .cool, or .mcool). Required when HIC_USE_FORCE=True.",
    )
    HIC_USE_FORCE: Boolean = Field(
        default=False,
        description="Apply Hi-C contact-guided force to the simulation.",
    )
    HIC_NORMALIZATION: str = Field(
        default="KR",
        description="Hi-C matrix normalisation method. Options: KR (default), VC, VC_SQRT, NONE.",
    )
    HIC_K_SCALE: float = Field(
        default=40.0,
        description=(
            "Global energy scale for the Hi-C force [kJ/mol]. Recommended range: 20–80 "
            "(default 40); higher values (up to ~160) consistently improved insulation-"
            "score validation in testing with no diagonal-decay cost, for N_BEADS >= 300 "
            "-- below that, gains are less consistent, so lower values may suit small "
            "systems better. Above 200 triggers a runtime warning. Each pair's force is "
            "already weighted by its own observed contact strength c_ij, so weak-evidence "
            "pairs stay proportionally soft without a separate knob."
        ),
    )
    HIC_RC: Optional[float] = Field(
        default=None,
        description=(
            "Explicit contact-radius scale [nm] for the Hi-C force's own target-distance "
            "mapping — independent of r_comp (the compartment/subcompartment block-copolymer "
            "force's interaction range, set in set_radiuses() from nucleus geometry). "
            "The two forces are physically unrelated and are no longer tied to the "
            "same value. Default None: the Hi-C force picks its own scale instead — "
            "auto-calibrated from the initial structure's pairwise distances when "
            "HIC_AUTO_SCALE is True (recommended, default), or a nucleus-scale "
            "fallback otherwise. Set this explicitly to pin the Hi-C force's contact "
            "radius to a fixed value and skip auto-calibration entirely, regardless "
            "of HIC_AUTO_SCALE."
        ),
    )
    HIC_BOLTZMANN_ALPHA: float = Field(
        default=4.0,
        description=(
            "Hi-C scaling-law exponent converting contact strength to a target 3-D "
            "distance via classic Boltzmann-inversion: r_target = r_min * c_ij^(-1/alpha), "
            "clipped to [r_min, rc] (see hic_force.build_boltzmann_force), then restrained "
            "there with a harmonic well weighted by c_ij — a genuine two-sided restraint "
            "(pulls if farther than its target, pushes if closer). Typical literature "
            "range: 3-4; higher values make the strength→distance mapping steeper (weak "
            "contacts hit the rc cap sooner), for every HIC_BOLTZMANN_KERNEL."
        ),
    )
    HIC_BOLTZMANN_KERNEL: str = Field(
        default="exponential",
        description=(
            "P(r) shape assumed for the Hi-C Boltzmann-PMF's equilibrium pair-distance "
            "distribution, i.e. the c_ij -> r_target inversion (see "
            "hic_force.get_boltzmann_p_func / _boltzmann_r_target). Options: "
            "exponential (default — classic Boltzmann distribution, c ~ exp(-r/lambda), "
            "lambda = rc/alpha), power_law (Hi-C scaling law, c ~ r^-alpha), "
            "sigmoid (bounded logistic contact probability, steepness set by alpha/rc). "
            "All three are exact functional inverses of their own P(r), reuse "
            "HIC_BOLTZMANN_ALPHA as their steepness knob, and are restrained with the "
            "same harmonic well — only the strength<->distance mapping's shape changes."
        ),
    )
    HIC_FORCE_OE: Boolean = Field(
        default=False,
        description=(
            "If False (default), target raw (KR-balanced) contact frequency directly. "
            "If True, target enrichment-above-background (OE ratio minus 1, floored at "
            "0) instead: pairs at/below the expected distance-decay baseline get zero "
            "force weight ('loose'); only enriched pairs attract. OE scores higher on "
            "PC1/OE-Pearson metrics in testing, but raw frequency gave better overall "
            "results in practice (diagonal decay + insulation), hence the default."
        ),
    )
    HIC_NOISE_INTENSITY: float = Field(
        default=0.0,
        description=(
            "Std-dev of Gaussian noise (same [0, 1] scale as c_ij) added to every Hi-C "
            "bond's contact strength, redrawn once per saved MD frame around its ORIGINAL "
            "data-derived value (never drifting cumulatively). This lets different contacts "
            "take turns pulling strongest from frame to frame, nudging the structure to "
            "explore nearby configurations instead of settling into one exact attractor, "
            "while every draw stays anchored to the real data. 0 (default) disables it; "
            "try 0.05-0.2 for mild exploration. Only applied during MD (SIM_RUN_MD=True)."
        ),
    )
    HIC_AUTO_SCALE: Boolean = Field(
        default=True,
        description=(
            "If True (default, recommended), recalibrate the Hi-C force's contact-radius "
            "scale (rc) from the median (50th percentile) of the actual initial "
            "structure's pairwise-distance distribution (hic_force.auto_contact_scale) "
            "instead of a fixed default, which can leave the force with no gradient at "
            "realistic bead separations and barely change the structure from its initial "
            "state. The 50th-percentile choice is fixed, not a separate tunable knob: it's "
            "simply the right scale for most contacted pairs — not just the closest — to "
            "start within the force's reach."
        ),
    )
    HIC_MAX_GAP: int = Field(
        default=10,
        description="Maximum gap fraction (in %) tolerated when interpolating missing Hi-C bins.",
    )
    HIC_INSULATION_WINDOW: int = Field(
        default=10,
        description=(
            "Half-width (beads) of the sliding square used by Hi-C validation's "
            "insulation-score metric (validation.insulation_score) — the full window is "
            "2x this value. Should roughly match the bead-scale size of a real TAD/domain "
            "for your resolution and N_BEADS; too small or too large relative to the "
            "actual domain size weakens the insulation_r correlation reported in "
            "hic_validation.npy even when the force itself is working well, since the "
            "metric is then measuring boundaries at the wrong genomic scale."
        ),
    )
    HIC_BLOCK_COPOLYMER: Boolean = Field(
        default=False,
        description=(
            "Opt-in: derive A/B compartments straight from the Hi-C matrix's own PC1 "
            "(sign-aligned to local contact density: dense -> B, sparse -> A) and feed "
            "them into the same block-copolymer force normally built from a .bed file — "
            "the Boltzmann force alone reproduces diagonal decay/insulation well but "
            "tends to miss PC1 correlation, so this adds it back with no .bed file "
            "needed. Suggested use: enable it when you have HIC_USE_FORCE=True, no "
            "compartment .bed file, and a region large enough to actually contain "
            "compartments (see HIC_BLOCK_COPOLYMER_MIN_BP) — e.g. a whole chromosome or "
            "genome-wide run; leave it off for TAD/region-scale runs, where it has no "
            "effect anyway. Only takes effect when HIC_USE_FORCE=True and no "
            "COMPARTMENT_PATH is given; auto-disables (with a warning) when the "
            "modelled region is below HIC_BLOCK_COPOLYMER_MIN_BP. Raises an error if a "
            ".bed-based compartment force (COB_USE_COMPARTMENT_BLOCKS / "
            "SCB_USE_SUBCOMPARTMENT_BLOCKS) is also enabled — pick one source of "
            "compartments, not both."
        ),
    )
    HIC_BLOCK_COPOLYMER_MIN_BP: float = Field(
        default=5_000_000,
        description=(
            "Minimum modelled region size (bp) for HIC_BLOCK_COPOLYMER to stay enabled. "
            "A/B compartments need several alternating domains to even be visible; below "
            "this a region is TAD-scale, not compartment-scale, so deriving compartments "
            "from it doesn't make sense."
        ),
    )
    GENE_TSV: str = Field(
        default=default_gene_path,
        description="A .tsv with genes and their locations in the genome.",
    )
    GENE_NAME: str = Field(default="", description="The name of the gene of interest.")
    GENE_ID: str = Field(default="", description="The id of the gene of interest.")
    GENE_WINDOW: int = Field(default=100000, description="The window around of the area around the gene of interest.")
    ATACSEQ_PATH: Optional[str] = Field(
        default=None, description="A .bw or .BigWig file path with atacseq data. It is not required."
    )
    OUT_PATH: str = Field(default="results", description="Output folder name.")
    LOC_START: Optional[int] = Field(default=None, description="Starting region coordinate.")
    LOC_END: Optional[int] = Field(default=None, description="Ending region coordinate.")
    CHROM: ChromStr = Field(
        default=None,
        description="Chromosome that corresponds the the modelling region of interest (in case that you do not want to model the whole genome).",
    )
    SHUFFLE_CHROMS: Boolean = Field(default=False, description="Shuffle the chromosomes.")
    SHUFFLING_SEED: int = Field(default=0, description="Shuffling random seed.")
    SAVE_PLOTS: Boolean = Field(default=True, description="Save plots.")
    POL_USE_HARMONIC_BOND: Boolean = Field(default=True, description="Use harmonic bond interaction.")
    POL_HARMONIC_BOND_R0: OpenMMQuantity = Field(
        default="0.1 nanometer", description="harmonic bond distance equilibrium constant"
    )
    POL_HARMONIC_BOND_K: OpenMMQuantity = Field(
        default="300000.0 kilojoules_per_mole/nanometer**2",
        description="harmonic bond force constant (fixed unit: kJ/mol/nm^2)",
    )
    POL_USE_HARMONIC_ANGLE: Boolean = Field(default=True, description="Use harmonic angle interaction.")
    POL_HARMONIC_ANGLE_R0: OpenMMQuantity = Field(
        default="3.141592653589793 radian", description="harmonic angle distance equilibrium constant"
    )
    POL_HARMONIC_ANGLE_CONSTANT_K: OpenMMQuantity = Field(
        default="100.0 kilojoules_per_mole/radian**2",
        description="harmonic angle force constant (fixed unit: kJ/mol/radian^2)",
    )
    LE_USE_HARMONIC_BOND: Boolean = Field(
        default=True, description="Use harmonic bond interaction for long range loops."
    )
    LE_FIXED_DISTANCES: Boolean = Field(
        default=False,
        description="For fixed distances between loops. False if you want to correlate with the hatmap strength.",
    )
    LE_HARMONIC_BOND_R0: OpenMMQuantity = Field(
        default="0.1 nanometer", description="harmonic bond distance equilibrium constant"
    )
    LE_HARMONIC_BOND_K: OpenMMQuantity = Field(
        default="30000.0 kilojoules_per_mole/nanometer**2",
        description="harmonic bond force constant (fixed unit: kJ/mol/nm^2)",
    )
    EV_USE_EXCLUDED_VOLUME: Boolean = Field(default=True, description="Use excluded volume.")
    EV_EPSILON: float = Field(default=100.0, description="Epsilon parameter.")
    EV_R_SMALL: float = Field(
        default=0.05, description="Add something small in denominator to make it not exploding all the time."
    )
    EV_POWER: float = Field(default=6.0, description="Power in the exponent of EV potential.")
    SC_USE_SPHERICAL_CONTAINER: Boolean = Field(default=False, description="Use Spherical container")
    SC_RADIUS1: Optional[OpenMMQuantity] = Field(default=None, description="Spherical container radius,")
    SC_RADIUS2: Optional[OpenMMQuantity] = Field(default=None, description="Spherical container radius,")
    SC_SCALE: float = Field(default=1000.0, description="Spherical container scaling factor")
    CHB_USE_CHROMOSOMAL_BLOCKS: Boolean = Field(default=False, description="Use Chromosomal Blocks.")
    CHB_KC: float = Field(default=0.3, description="Block copolymer width parameter.")
    CHB_DE: float = Field(default=1e-04, description="Energy factor for block copolymer chromosomal model.")
    COB_USE_COMPARTMENT_BLOCKS: Boolean = Field(default=False, description="Use Compartment Blocks.")
    COB_DISTANCE: Optional[OpenMMQuantity] = Field(
        default=None, description="Block copolymer equilibrium distance for chromosomal blocks."
    )
    COB_EA: float = Field(default=1.0, description="Energy strength for A compartment.")
    COB_EB: float = Field(default=2.0, description="Energy strength for B compartment.")
    SCB_USE_SUBCOMPARTMENT_BLOCKS: Boolean = Field(default=False, description="Use Subcompartment Blocks.")
    SCB_DISTANCE: Optional[OpenMMQuantity] = Field(default=None, description="Block copolymer equilibrium distance for chromosomal blocks.")
    SCB_EA1: float = Field(default=1.0, description="Energy strength for A1 compartment.")
    SCB_EA2: float = Field(default=1.33, description="Energy strength for A2 compartment.")
    SCB_EB1: float = Field(default=1.66, description="Energy strength for B1 compartment.")
    SCB_EB2: float = Field(default=2.0, description="Energy strength for B2 compartment.")
    IBL_USE_B_LAMINA_INTERACTION: Boolean = Field(default=False, description="Interactions of B compartment with lamina.")
    IBL_SCALE: float = Field(default=400.0, description="Scaling factor for B comoartment interaction with lamina.")
    CF_USE_CENTRAL_FORCE: Boolean = Field(default=False, description="Attraction of smaller chromosomes.")
    CF_STRENGTH: float = Field(default=20.0, description="Strength of Interaction")
    NUC_DO_INTERPOLATION: Boolean = Field(default=False, description="Attraction of smaller chromosomes.")
    MAX_NUCS_PER_BEAD: int = Field(default=4, description="Maximum amount of nucleosomes per single bead.")
    NUC_RADIUS: float = Field(default=0.1, description="The radius of the single nucleosome helix.")
    POINTS_PER_NUC: int = Field(default=20, description="The number of points that consist a nucleosome helix.")
    PHI_NORM: float = Field(default=0.6283185307179586, description="Zig zag angle. ")
    SIM_RUN_MD: Boolean = Field(default=False, description="Do you want to run MD simulation?")
    SIM_N_STEPS: int = Field(default=10000, description="Number of steps in MD simulation")
    SIM_ERROR_TOLERANCE: float = Field(default=0.01, description="Error tolerance for variable MD simulation")
    SIM_AMD_ALPHA: float = Field(default=100.0, description="Alpha of AMD simulation.")
    SIM_AMD_E: float = Field(default=1000.0, description="E (energy) of AMD simulation.")
    SIM_SAMPLING_STEP: int = Field(default=100, description="It determines in t how many steps we save a structure.")
    SIM_INTEGRATOR_TYPE: str = Field(default="langevin", description="Alternative: langevin, verlet")
    SIM_INTEGRATOR_STEP: OpenMMQuantity = Field(default="1 femtosecond", description="The step of integrator.")
    SIM_FRICTION_COEFF: float = Field(
        default=0.5, description="Friction coefficient (Used only with langevin integrator)"
    )
    SIM_TEMPERATURE: OpenMMQuantity = Field(default="310 kelvin", description="Simulation temperature")
    SIM_SET_INITIAL_VELOCITIES: Boolean = Field(
        default=True,
        description=(
            "Initialize the MD simulation with a random initial velocity field, drawn from the "
            "Maxwell-Boltzmann distribution at SIM_TEMPERATURE (OpenMM's setVelocitiesToTemperature), "
            "seeded by SHUFFLING_SEED. True by default so every run starts from a physically "
            "realistic, randomized velocity field instead of the all-zero velocities OpenMM uses "
            "otherwise."
        ),
    )
    TRJ_FRAMES: int | None = Field(
        default=None,
        description=(
            "Number of CIF frames to save during MD.  When set, SIM_SAMPLING_STEP is "
            "overridden to SIM_N_STEPS // TRJ_FRAMES so that exactly TRJ_FRAMES "
            "structures are written.  When None, SIM_SAMPLING_STEP is used as-is."
        ),
    )

    EV_FORCE_TYPE: str = Field(
        default="powerlaw",
        description="Excluded volume functional form. Options: powerlaw (default), gaussian_core.",
    )

    COB_FORCE_TYPE: str = Field(
        default="gaussian",
        description="Compartment block interaction functional form. Options: gaussian (default), yukawa, theta.",
    )

    SCB_FORCE_TYPE: str = Field(
        default="gaussian",
        description="Subcompartment block interaction functional form. Options: gaussian (default), yukawa, theta.",
    )

    BLAMINA_FORCE_TYPE: str = Field(
        default="sin",
        description="B-lamina interaction functional form. Options: sin (default), gaussian_shell, harmonic_shell, logistic_shell.",
    )

    LE_LOOP_FORCE_TYPE: str = Field(
        default="harmonic",
        description="Loop extrusion bond functional form. Options: harmonic (default), fene_soft, gaussian_tether.",
    )

    CHB_FORCE_TYPE: str = Field(
        default="polynomial",
        description=(
            "Chromosome self-attraction kernel controlling global compaction into globules. "
            "Options: polynomial (default, handcrafted potential), "
            "gaussian (soft collapse kernel), "
            "saturating (soft-core bounded attraction)."
    ))

    CENTRAL_FORCE_TYPE: str = Field(
        default="harmonic",
        description=(
            "Central nucleolar attraction functional form controlling radial bias toward nucleus center. "
            "Encodes chromosome-size dependent positioning. "
            "Options: "
            "harmonic (default, quadratic confinement around R1), "
            "gaussian (soft nucleolar enrichment field), "
            "logistic (soft-core radial partitioning with smooth boundary)."
    ))