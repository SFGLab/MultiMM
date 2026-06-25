#########################################################################
########### CREATOR: SEBASTIAN KORSAK, WARSAW 2024 ######################
#########################################################################
import argparse
import configparser
import logging
import os
import sys
import tarfile
from enum import Enum

from openmm.unit import Quantity

from .config import SimulationConfig
from .model import MultiMM
from .utils import chrom_sizes
from .logger import setup_logger

setup_logger()
logger = logging.getLogger(__name__)

STARTUP_BANNER_LINES = [
    "#########################################################################",
    "# 🧬 MultiMM Chromatin Simulation Platform 🧬",
    "#########################################################################",
    "# Creator: Sebastian Korsak (Warsaw,)",
    "# Web-server & infrastructure: Patryk Prusak",
    "# email: s.korsak@datascience.edu.pl",
    "#",
    "# 🚀 Starting simulation pipeline...",
    "# ✨ Wishing you smooth, stable and beautiful chromatin dynamics!",
    "# 🧪 May your contacts be meaningful and your loops well-formed",
    "# 🔬 Happy modeling!",
    "#########################################################################",
]

# ANSI colors (soft scientific palette)
COLORS = [
    "\033[96m",  # cyan
    "\033[95m",  # magenta
    "\033[94m",  # blue
    "\033[92m",  # green
    "\033[93m",  # yellow
    "\033[91m",  # red
]

RESET = "\033[0m"


def print_startup_banner(logger):
    """
    Print colorful MultiMM startup banner.
    """

    for i, line in enumerate(STARTUP_BANNER_LINES):
        color = COLORS[i % len(COLORS)]
        logger.info(f"{color}{line}{RESET}")

class Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for s in self.streams:
            s.write(data)
            s.flush()

    def flush(self):
        for s in self.streams:
            s.flush()

class ArgumentChanger:

    # ------------------------------------------------------------
    # ANSI colors (works in most terminals, including Linux/SSH)
    # ------------------------------------------------------------
    ORANGE = "\033[38;5;208m"
    RESET = "\033[0m"
    BOLD = "\033[1m"

    def __init__(self, args, chrom_sizes):
        self.args = args
        self.chrom_sizes = chrom_sizes

        # store original values for diff reporting
        self._original_values = {}

    def set_arg(self, name, value):
        """
        Set argument value in attribute and store change history.
        """
        if hasattr(self.args, name):
            old_value = getattr(self.args, name, None)

            # store original only once
            if name not in self._original_values:
                self._original_values[name] = old_value

            setattr(self.args, name, value)

        else:
            logger.warning(f"Argument '{name}' not found in args object.")

    def _report_changes(self):
        """
        Print all modified parameters in a readable way.
        """
        if not self._original_values:
            return

        logger.warning(
            f"{self.ORANGE}{self.BOLD}"
            "MODELLING LEVEL OVERRIDE ACTIVE: parameters have been overwritten."
            f"{self.RESET}"
        )

        print("\nChanged parameters:")
        print("-" * 60)

        for k, old_v in self._original_values.items():
            new_v = getattr(self.args, k, None)
            if old_v != new_v:
                print(f"{k:35s} : {old_v}  ->  {new_v}")

        print("-" * 60 + "\n")

    def convenient_argument_changer(self):

        self.set_arg("NUC_DO_INTERPOLATION", False)
        self.set_arg("ATACSEQ_PATH", None)

        modelling_level = self.args.MODELLING_LEVEL

        level = str(modelling_level).lower()

        if level == "gene":

            logger.warning(
                f"{self.ORANGE}{self.BOLD}"
                "Gene-level modelling activated. This will overwrite parameters."
                f"{self.RESET}"
            )

            self.set_arg("N_BEADS", 1000)
            self.set_arg("SC_USE_SPHERICAL_CONTAINER", False)
            self.set_arg("CHB_USE_CHROMOSOMAL_BLOCKS", False)
            self.set_arg("SCB_USE_SUBCOMPARTMENT_BLOCKS", False)
            self.set_arg("COB_USE_COMPARTMENT_BLOCKS", False)
            self.set_arg("IBL_USE_B_LAMINA_INTERACTION", False)
            self.set_arg("CF_USE_CENTRAL_FORCE", False)
            self.set_arg("SHUFFLE_CHROMS", False)
            self.set_arg("SIM_RUN_MD", True)
            self.set_arg("SIM_N_STEPS", 10000)

        elif level in ("region", "loc"):

            logger.warning(
                f"{self.ORANGE}{self.BOLD}"
                "Region-level modelling activated. Overwriting parameters."
                f"{self.RESET}"
            )

            self.set_arg("N_BEADS", 5000)
            self.set_arg("SC_USE_SPHERICAL_CONTAINER", False)
            self.set_arg("CHB_USE_CHROMOSOMAL_BLOCKS", False)
            self.set_arg("SCB_USE_SUBCOMPARTMENT_BLOCKS", False)
            self.set_arg("COB_USE_COMPARTMENT_BLOCKS", bool(self.args.COMPARTMENT_PATH))
            self.set_arg("IBL_USE_B_LAMINA_INTERACTION", False)
            self.set_arg("CF_USE_CENTRAL_FORCE", False)
            self.set_arg("SIM_RUN_MD", True)
            self.set_arg("SIM_N_STEPS", 10000)

        elif level in ("chromosome", "chrom"):

            logger.warning(
                f"{self.ORANGE}{self.BOLD}"
                "Chromosome-level modelling activated. Overwriting parameters."
                f"{self.RESET}"
            )

            self.set_arg("N_BEADS", 20000)
            self.set_arg("SC_USE_SPHERICAL_CONTAINER", False)
            self.set_arg("CHB_USE_CHROMOSOMAL_BLOCKS", False)
            self.set_arg("SCB_USE_SUBCOMPARTMENT_BLOCKS", False)
            self.set_arg("COB_USE_COMPARTMENT_BLOCKS", bool(self.args.COMPARTMENT_PATH))
            self.set_arg("IBL_USE_B_LAMINA_INTERACTION", False)
            self.set_arg("CF_USE_CENTRAL_FORCE", False)
            self.set_arg("SIM_RUN_MD", True)
            self.set_arg("SIM_N_STEPS", 10000)
            self.set_arg("LOC_START", 1)
            self.set_arg("LOC_END", self.chrom_sizes[self.args.CHROM])

        elif level in ("gw", "genome"):

            logger.warning(
                f"{self.ORANGE}{self.BOLD}"
                "Genome-wide modelling activated. Overwriting parameters."
                f"{self.RESET}"
            )

            self.set_arg("N_BEADS", 200000)
            self.set_arg("SC_USE_SPHERICAL_CONTAINER", True)
            self.set_arg("CHB_USE_CHROMOSOMAL_BLOCKS", False)
            self.set_arg("SCB_USE_SUBCOMPARTMENT_BLOCKS", False)
            self.set_arg("COB_USE_COMPARTMENT_BLOCKS", bool(self.args.COMPARTMENT_PATH))
            self.set_arg(
                "IBL_USE_B_LAMINA_INTERACTION",
                bool(self.args.COMPARTMENT_PATH),
            )
            self.set_arg("CF_USE_CENTRAL_FORCE", False)
            self.set_arg("SIM_RUN_MD", False)
            self.set_arg("SIM_N_STEPS", 10000)

        # final summary
        if self.args.MODELLING_LEVEL:
            self._report_changes()

def args_tests(args):

    def check_file(path, name, ext_hint=None):
        """
        Validate optional input files.
        If provided, ensure they exist.
        """
        if path is None or path == "":
            return  # optional → OK

        if not os.path.exists(path):
            ext_msg = f" (expected {ext_hint})" if ext_hint else ""
            raise ValueError(
                f"\033[91m{name} file was provided but not found: {path}{ext_msg}\033[0m"
            )

    # -----------------------------------------
    # REQUIRED INPUT
    # -----------------------------------------
    if args.LOOPS_PATH is None or args.LOOPS_PATH == "":
        raise ValueError(
            "\033[91mLoops interaction data is required to run MultiMM."
            "Please provide a valid .bedpe file via LOOPS_PATH.\033[0m"
        )

    check_file(args.LOOPS_PATH, "Loops (.bedpe)", ".bedpe")

    # -----------------------------------------
    # OPTIONAL INPUT VALIDATION (if provided)
    # -----------------------------------------
    check_file(args.COMPARTMENT_PATH, "Compartment data", ".bed")
    check_file(args.ATACSEQ_PATH, "Nucleosome/ATAC data", ".bigwig")

    if args.LOOPS_PATH is None or args.LOOPS_PATH == "":
        raise ValueError(
            "\033[91mInteraction data is required to run MultiMM."
            "Please provide a .bedpe file via LOOPS_PATH.\033[0m"
        )

    elif (args.COMPARTMENT_PATH is None or args.COMPARTMENT_PATH == "") and args.COB_USE_COMPARTMENT_BLOCKS:
        raise ValueError(
            "\033[91mCompartment modeling is enabled, but no compartment data was provided."
            "Please supply a .bed file or disable COB_USE_COMPARTMENT_BLOCKS.\033[0m"
        )

    elif args.NUC_DO_INTERPOLATION and args.ATACSEQ_PATH is None:
        raise ValueError(
            "\033[91mNucleosome interpolation is enabled, but no occupancy data was found."
            "Provide a .bigwig file via ATACSEQ_PATH or disable NUC_DO_INTERPOLATION.\033[0m"
        )

    elif (args.COMPARTMENT_PATH is None or args.COMPARTMENT_PATH == "") and args.SCB_USE_SUBCOMPARTMENT_BLOCKS:
        raise ValueError(
            "\033[91mSubcompartment modeling requires input data."
            "Please provide a .bed file or disable SCB_USE_SUBCOMPARTMENT_BLOCKS.\033[0m"
        )

    elif args.COMPARTMENT_PATH is None and args.IBL_USE_B_LAMINA_INTERACTION:
        raise ValueError(
            "\033[91mLamina interactions depend on compartment annotations."
            "Please provide a compartment .bed file or disable IBL_USE_B_LAMINA_INTERACTION.\033[0m"
        )

    elif args.IBL_USE_B_LAMINA_INTERACTION and not (
        args.SCB_USE_SUBCOMPARTMENT_BLOCKS or args.COB_USE_COMPARTMENT_BLOCKS
    ):
        raise ValueError(
            "\033[91mLamina interactions are enabled but no compartment-based forces are active."
            "Enable COB_USE_COMPARTMENT_BLOCKS or SCB_USE_SUBCOMPARTMENT_BLOCKS, or disable lamina interactions.\033[0m"
        )

    elif args.CF_USE_CENTRAL_FORCE and args.CHROM is not None:
        raise ValueError(
            "\033[91mCentral force attraction to the nucleolus is typically used for whole-genome simulations."
            "Since you are modeling a single chromosome or region, consider disabling CF_USE_CENTRAL_FORCE.\033[0m"
        )

    elif args.CHB_USE_CHROMOSOMAL_BLOCKS and args.CHROM is not None:
        logger.warning(
            "\033[93mChromosomal block interactions are more meaningful in multi-chromosome systems."
            "You may want to disable CHB_USE_CHROMOSOMAL_BLOCKS for single-chromosome simulations.\033[0m"
        )

    if args.SHUFFLE_CHROMS and (args.CHROM is not None and args.CHROM != ""):
        logger.warning(
            "\033[93mChromosome shuffling is enabled, but you are simulating a specific chromosomal region."
            "This option usually makes more sense when working with multiple chromosomes.\033[0m"
        )

    if args.CHROM is not None and args.IBL_USE_B_LAMINA_INTERACTION:
        logger.warning(
            "\033[93mLamina interactions are enabled."
            "This is not incorrect, but they are typically more relevant in whole-genome simulations.\033[0m"
        )

    if args.CHROM is not None and args.SC_USE_SPHERICAL_CONTAINER:
        logger.warning(
            "\033[93mA spherical container is being used."
            "This is fine, but it is generally more meaningful when modeling the full genome.\033[0m"
        )

    if (not args.POL_USE_HARMONIC_BOND) or (not args.POL_USE_HARMONIC_ANGLE) or (not args.EV_USE_EXCLUDED_VOLUME):
        logger.warning(
            "\033[93mSome fundamental backbone forces are disabled."
            "Make sure this is intentional, as it may strongly affect the physical behavior of the polymer.\033[0m"
        )

    if args.CHB_USE_CHROMOSOMAL_BLOCKS:
        logger.warning(
            "\033[93mChromosomal block forces are enabled."
            "These are approximate and may not always reflect biological reality."
            "Consider checking the documentation to ensure they fit your use case.\033[0m"
        )


def my_config_parser(config_parser: configparser.ConfigParser):
    """Helper function that makes flat list arg name, and it's value from
    ConfigParser object."""
    sections = config_parser.sections()
    all_nested_fields = [dict(config_parser[s]) for s in sections]
    defaults_dict = dict(config_parser.defaults())
    if defaults_dict:
        all_nested_fields.append(defaults_dict)
    args_cp = []
    for section_fields in all_nested_fields:
        for name, value in section_fields.items():
            args_cp.append((name, value))
    return args_cp


def get_config():
    """Prepare list of arguments.

    First, defaults are set. Then, optionally config file values.
    Finally, CLI arguments overwrite everything. Then internal changes
    are applied.
    """
    logger.info("Reading config...")

    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("-c", "--config_file", help="Specify config file (ini format)", metavar="FILE")

    for field_name, field in SimulationConfig.model_fields.items():
        arg_parser.add_argument(f"--{field_name.lower()}", help=field.description)

    args_ap = arg_parser.parse_args()
    args_dict = vars(args_ap)

    raw_config = {}

    if args_ap.config_file:
        config_parser = configparser.ConfigParser()
        config_parser.read(args_ap.config_file)
        args_cp = my_config_parser(config_parser)

        for cp_arg in args_cp:
            name, value = cp_arg
            raw_config[name.upper()] = value

    for name, value in args_dict.items():
        if name == "config_file":
            continue
        if value is not None:
            raw_config[name.upper()] = value

    try:
        config_obj = SimulationConfig(**raw_config)
    except Exception as e:
        logger.error(f"Configuration validation failed: {e}")
        raise e

    changer = ArgumentChanger(config_obj, chrom_sizes)
    changer.convenient_argument_changer()

    write_config(config_obj)

    return config_obj


def write_config(args):
    """Write the automatically generated config to the metadata directory."""
    metadata_dir = os.path.join(args.OUT_PATH, "metadata")
    os.makedirs(metadata_dir, exist_ok=True)
    config_path = os.path.join(metadata_dir, "config_auto.ini")

    config = configparser.ConfigParser()
    config["DEFAULT"] = {}

    for name, value in args.model_dump().items():
        if isinstance(value, Quantity):
            config["DEFAULT"][name] = f"{value._value} {value.unit.get_name()}"
        elif isinstance(value, Enum):
            config["DEFAULT"][name] = value.value
        elif value is None:
            config["DEFAULT"][name] = ""
        else:
            config["DEFAULT"][name] = str(value)

    with open(config_path, "w") as config_file:
        config.write(config_file)

    logger.info(f"Configuration saved to {config_path}")


def archive_run(run_path):
    """
    Compress a run directory and remove the original folder.
    """

    tar_path = run_path + ".tar.gz"

    logger.info(f"Creating archive: {tar_path}")

    with tarfile.open(tar_path, "w:gz") as tar:
        tar.add(run_path, arcname=os.path.basename(run_path))

    # Safety check before deleting
    if os.path.exists(tar_path) and os.path.getsize(tar_path) > 0:
        logger.info(f"Archive created successfully. Removing {run_path}")
        shutil.rmtree(run_path)
    else:
        raise RuntimeError(
            f"Archive creation failed ({tar_path}). "
            f"Original directory was NOT deleted."
        )

    logger.info(f"Archived run stored at: {tar_path}")

def main():
    try:
        print_startup_banner(logger)

        args = get_config()
        args_tests(args)

        log_dir = os.path.join(args.OUT_PATH, "metadata")
        os.makedirs(log_dir, exist_ok=True)

        log_path = os.path.join(log_dir, "output.log")

        with open(log_path, "w") as log_file:

            original_stdout = sys.stdout
            original_stderr = sys.stderr

            sys.stdout = Tee(original_stdout, log_file)
            sys.stderr = Tee(original_stderr, log_file)

            try:

                name = args.OUT_PATH

                if args.GENERATE_ENSEMBLE:

                    for i in range(args.N_ENSEMBLE):

                        args.SHUFFLING_SEED = i

                        run_path = os.path.join(name, f"run_{i:04d}")
                        args.OUT_PATH = run_path

                        os.makedirs(run_path, exist_ok=True)

                        md = MultiMM(args)
                        md.run()

                        archive_run(run_path)

                else:

                    md = MultiMM(args)
                    md.run()

            finally:
                sys.stdout = original_stdout
                sys.stderr = original_stderr

        sys.exit(0)

    except Exception as e:
        logger.error(f"ERROR: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
