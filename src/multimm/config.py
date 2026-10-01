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
    }

    @model_validator(mode="after")
    def warn_hic_k_scale(self) -> "SimulationConfig":
        """Emit a runtime warning when HIC_K_SCALE is outside the recommended range."""
        k = self.HIC_K_SCALE
        if k > 200.0:
            logger.warning(
                "HIC_K_SCALE=%.1f is above the recommended maximum of 200 kJ/mol. "
                "Very high values can over-constrain the polymer, reduce conformational "
                "diversity, and cause MD instability. Consider reducing to 10–200.",
                k,
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
        description="It should be a .bed file with subcompartments from Calder (or something in the same format).",
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
        default=20.0,
        description=(
            "Global energy scale for the Hi-C force [kJ/mol]. Recommended range: 5–20. "
            "Values above 30 can freeze MD; above 200 triggers a runtime warning. Weak "
            "contacts are already softened independently via HIC_WEIGHT_POWER."
        ),
    )
    HIC_KERNEL: str = Field(
        default="gaussian",
        description=(
            "Distance → contact-probability kernel P(r) used by the Hi-C force.  Each kernel "
            "has its own extra parameter(s), named HIC_<KERNEL>_* below, that only take effect "
            "when that kernel is selected.  Options: "
            "'gaussian' (default, P=exp(-r²/2σ²), width HIC_GAUSSIAN_SIGMA), "
            "'power_law' / 'sigmoid' (P=1/(1+(r/r_c)^alpha), steepness HIC_POWERLAW_ALPHA), "
            "'exponential' (P=exp(-r/r_c), persistent long-range pull), "
            "'erfc' (soft step at r_c with width HIC_ERFC_SIGMA — closest to a binary Hi-C "
            "contact definition), "
            "'rouse' (separation-aware Gaussian-chain model, P=erfc(r/sqrt(2*s*b²)) with "
            "s=|i-j| in beads and Kuhn length HIC_ROUSE_KUHN_LENGTH — automatically reproduces "
            "the expected diagonal decay per genomic separation)."
        ),
    )
    HIC_POWERLAW_ALPHA: float = Field(
        default=3.0,
        description=(
            "['power_law'/'sigmoid' kernel only] Steepness of the contact-probability sigmoid: "
            "P(r) = 1 / (1 + (r/r_c)^alpha).  Has no effect unless HIC_KERNEL is 'power_law' or "
            "'sigmoid'.  "
            "Controls how sharply the force transitions between attraction and repulsion "
            "at the contact radius r_c.  "
            "Lower values (2) give a broad, gradual transition; higher values (4–6) "
            "give a sharper, TAD-like step.  Recommended range: 2–4."
        ),
    )
    HIC_GAUSSIAN_SIGMA: Optional[float] = Field(
        default=None,
        description=(
            "['gaussian' kernel only] Width σ [nm], P(r)=exp(-r²/2σ²).  Has no effect unless "
            "HIC_KERNEL='gaussian'.  Defaults to r_comp (a nucleus-scale reach, further "
            "auto-calibrated from the initial structure when HIC_AUTO_SCALE is True) when "
            "not set."
        ),
    )
    HIC_ERFC_SIGMA: Optional[float] = Field(
        default=None,
        description=(
            "['erfc' kernel only] Softening width σ_s [nm] of the step at r_c.  Has no effect "
            "unless HIC_KERNEL='erfc'.  Smaller values make the step sharper (σ_s→0 recovers a "
            "binary Hi-C contact definition).  Defaults to 0.3 * r_comp when not set."
        ),
    )
    HIC_ROUSE_KUHN_LENGTH: Optional[float] = Field(
        default=None,
        description=(
            "['rouse' kernel only] Kuhn (statistical segment) length b [nm], where "
            "<r²(s)> = s*b² for genomic separation s (in beads).  Has no effect unless "
            "HIC_KERNEL='rouse'.  Defaults to r_comp when not set."
        ),
    )
    HIC_WEIGHT_POWER: float = Field(
        default=1.0,
        description=(
            "Exponent β in the per-pair force weight w_ij = c_ij^β. With β=1 (default) "
            "force is directly proportional to contact strength c_ij, so weak-evidence "
            "pairs stay nearly inert. β<1 softens weak contacts less; β>1 suppresses them more."
        ),
    )
    HIC_THRESHOLD: float = Field(
        default=0.01,
        description=(
            "Minimum normalised contact value c_ij to include a bond at all — a sparsity "
            "cutoff only (bonds still scale continuously with c_ij via HIC_WEIGHT_POWER). "
            "Recommended range: 0.005–0.05."
        ),
    )
    HIC_FORCE_OE: Boolean = Field(
        default=True,
        description=(
            "If True (default), target enrichment-above-background (OE ratio minus 1, "
            "floored at 0) instead of raw contact frequency: pairs at or below the expected "
            "distance-decay baseline (not enriched — 'towards -1' on the displayed "
            "log2(O/E) scale) get exactly zero force weight, so they are neither pulled "
            "together nor pushed apart ('loose'); only genuinely enriched pairs attract, "
            "most strongly for the most-enriched ('towards +1'). Set False to instead "
            "target raw (KR-balanced) contact frequency directly, which also pulls "
            "short-range/background pairs together just from their high absolute count."
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
            "If True (default, recommended), recalibrate the Hi-C force's distance-kernel "
            "scale (r_comp and any unset HIC_GAUSSIAN_SIGMA/HIC_ERFC_SIGMA/"
            "HIC_ROUSE_KUHN_LENGTH) from a percentile of the actual initial structure's "
            "pairwise-distance distribution (hic_force.auto_contact_scale) instead of a "
            "fixed default, which can leave the force with no gradient at realistic bead "
            "separations and barely change the structure from its initial state."
        ),
    )
    HIC_AUTO_SCALE_PERCENTILE: float = Field(
        default=50.0,
        description=(
            "Percentile of the initial structure's pairwise distances used to calibrate "
            "the Hi-C force's scale when HIC_AUTO_SCALE is True. Higher than validation's "
            "equivalent percentile (which favors visual contrast) because the force needs "
            "reach: most contacted pairs, not just the closest, should feel a gradient. "
            "Default 50 (median); raise for a more diffuse initial structure, lower if "
            "local contact resolution suffers."
        ),
    )
    HIC_MAX_GAP: int = Field(
        default=10,
        description="Maximum gap fraction (in %) tolerated when interpolating missing Hi-C bins.",
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