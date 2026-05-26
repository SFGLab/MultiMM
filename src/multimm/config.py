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
    LOOPS_PATH: str = Field(default="", description="A .bedpe file path with loops. It is required.")
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
    CHB_DE: float = Field(default=1e-05, description="Energy factor for block copolymer chromosomal model.")
    COB_USE_COMPARTMENT_BLOCKS: Boolean = Field(default=False, description="Use Compartment Blocks.")
    COB_DISTANCE: Optional[OpenMMQuantity] = Field(
        default=None, description="Block copolymer equilibrium distance for chromosomal blocks."
    )
    COB_EA: float = Field(default=1.0, description="Energy strength for A compartment.")
    COB_EB: float = Field(default=2.0, description="Energy strength for B compartment.")
    SCB_USE_SUBCOMPARTMENT_BLOCKS: Boolean = Field(default=False, description="Use Subcompartment Blocks.")
    SCB_DISTANCE: Optional[OpenMMQuantity] = Field(
        default=None, description="Block copolymer equilibrium distance for chromosomal blocks."
    )
    SCB_EA1: float = Field(default=1.0, description="Energy strength for A1 compartment.")
    SCB_EA2: float = Field(default=1.33, description="Energy strength for A2 compartment.")
    SCB_EB1: float = Field(default=1.66, description="Energy strength for B1 compartment.")
    SCB_EB2: float = Field(default=2.0, description="Energy strength for B2 compartment.")
    IBL_USE_B_LAMINA_INTERACTION: Boolean = Field(
        default=False, description="Interactions of B compartment with lamina."
    )
    IBL_SCALE: float = Field(default=400.0, description="Scaling factor for B comoartment interaction with lamina.")
    CF_USE_CENTRAL_FORCE: Boolean = Field(default=False, description="Attraction of smaller chromosomes.")
    CF_STRENGTH: float = Field(default=10.0, description="Strength of Interaction")
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
    SIM_SET_INITIAL_VELOCITIES: Boolean = Field(
        default=False, description="Sets initial velocities based on Boltzmann distribution"
    )
    SIM_TEMPERATURE: OpenMMQuantity = Field(default="310 kelvin", description="Simulation temperature")
    TRJ_FRAMES: int = Field(default=2000, description="Number of trajectory frames to save.")
