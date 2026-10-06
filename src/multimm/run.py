#########################################################################
########### CREATOR: SEBASTIAN KORSAK, WARSAW 2024 ######################
#########################################################################
import argparse
import configparser
import difflib
import logging
import os
import sys
import tarfile
import shutil
from enum import Enum

from openmm.unit import Quantity

from .config import SimulationConfig
from .model import MultiMM
from .utils import chrom_sizes
from .logger import setup_logger

setup_logger()
logger = logging.getLogger(__name__)

RESET  = "\033[0m"
BOLD   = "\033[1m"
CYAN   = "\033[38;5;75m"   # steel blue – professional, readable


def print_startup_banner(logger):
    """Print a clean, single-color startup banner for MultiMM."""

    W = 72  # total banner width (inner)
    border = CYAN + BOLD + "═" * W + RESET

    def _row(text: str = "") -> str:
        pad = W - 2 - len(text)
        return CYAN + BOLD + "║ " + RESET + text + " " * max(pad, 0) + CYAN + BOLD + "║" + RESET

    lines = [
        border,
        _row(),
        _row("  MultiMM  ·  Chromatin 3D Structure Simulation Platform"),
        _row(),
        _row("  Creator : Sebastian Korsak  (Warsaw, Poland)"),
        _row("  Nucleosome interpolation : Krzysztof Banecki"),
        _row("  Web-server & infrastructure : Patryk Prusak"),
        _row(),
        _row("  Contact : s.korsak@datascience.edu.pl"),
        _row("            d.plewczynski@datascience.edu.pl"),
        _row(),
        _row("  Starting simulation pipeline … good luck!"),
        _row(),
        border,
    ]
    for line in lines:
        logger.info(line)

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

    def isatty(self) -> bool:
        # Delegate to the first stream (the real terminal).
        # This preserves ANSI colour detection after Tee wraps sys.stdout.
        return hasattr(self.streams[0], "isatty") and self.streams[0].isatty()

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
    # INPUT FILE EXISTENCE CHECKS (if provided)
    # -----------------------------------------
    check_file(args.LOOPS_PATH, "Loops (.bedpe)", ".bedpe")
    check_file(args.COMPARTMENT_PATH, "Compartment data", ".bed")
    check_file(args.ATACSEQ_PATH, "Nucleosome/ATAC data", ".bigwig")
    check_file(args.HIC_PATH, "Hi-C contact matrix", ".hic/.cool/.mcool")

    # -----------------------------------------
    # REQUIRED COMBINATIONS
    # -----------------------------------------
    if args.LE_USE_HARMONIC_BOND and (args.LOOPS_PATH is None or args.LOOPS_PATH == ""):
        raise ValueError(
            "\033[91mLE_USE_HARMONIC_BOND=True but no LOOPS_PATH provided. "
            "Please supply a .bedpe file or disable LE_USE_HARMONIC_BOND.\033[0m"
        )

    if args.HIC_USE_FORCE and (args.HIC_PATH is None or args.HIC_PATH == ""):
        raise ValueError(
            "\033[91mHIC_USE_FORCE=True but no HIC_PATH provided. "
            "Please supply a .hic / .cool / .mcool file or disable HIC_USE_FORCE.\033[0m"
        )

    if (args.COMPARTMENT_PATH is None or args.COMPARTMENT_PATH == "") and args.COB_USE_COMPARTMENT_BLOCKS:
        raise ValueError(
            "\033[91mCompartment modeling is enabled, but no compartment data was provided. "
            "Please supply a .bed file or disable COB_USE_COMPARTMENT_BLOCKS.\033[0m"
        )

    if args.NUC_DO_INTERPOLATION and (args.ATACSEQ_PATH is None or args.ATACSEQ_PATH == ""):
        raise ValueError(
            "\033[91mNucleosome interpolation is enabled, but no occupancy data was found. "
            "Provide a .bigwig file via ATACSEQ_PATH or disable NUC_DO_INTERPOLATION.\033[0m"
        )

    if (args.COMPARTMENT_PATH is None or args.COMPARTMENT_PATH == "") and args.SCB_USE_SUBCOMPARTMENT_BLOCKS:
        raise ValueError(
            "\033[91mSubcompartment modeling requires input data. "
            "Please provide a .bed file or disable SCB_USE_SUBCOMPARTMENT_BLOCKS.\033[0m"
        )

    if (args.COMPARTMENT_PATH is None or args.COMPARTMENT_PATH == "") and args.IBL_USE_B_LAMINA_INTERACTION:
        raise ValueError(
            "\033[91mLamina interactions depend on compartment annotations. "
            "Please provide a compartment .bed file or disable IBL_USE_B_LAMINA_INTERACTION.\033[0m"
        )

    if args.IBL_USE_B_LAMINA_INTERACTION and not (
        args.SCB_USE_SUBCOMPARTMENT_BLOCKS or args.COB_USE_COMPARTMENT_BLOCKS
    ):
        raise ValueError(
            "\033[91mLamina interactions are enabled but no compartment-based forces are active. "
            "Enable COB_USE_COMPARTMENT_BLOCKS or SCB_USE_SUBCOMPARTMENT_BLOCKS, or disable lamina interactions.\033[0m"
        )

    if args.CF_USE_CENTRAL_FORCE and args.CHROM is not None and args.CHROM != "":
        logger.warning(
            "\033[93mCentral force (nucleolar attraction) is enabled for a single-chromosome/region run. "
            "It is typically used in whole-genome simulations; consider disabling CF_USE_CENTRAL_FORCE.\033[0m"
        )

    # -----------------------------------------
    # ADVISORY WARNINGS
    # -----------------------------------------

    # All three data sources active simultaneously — warn about potential redundancy
    if (
        args.HIC_USE_FORCE
        and args.LE_USE_HARMONIC_BOND
        and not (args.COMPARTMENT_PATH is None or args.COMPARTMENT_PATH == "")
    ):
        logger.warning(
            "\033[93mHi-C contact force, loop extrusion, AND compartment forces are all active. "
            "This combination is valid but may over-constrain the structure. "
            "Consider using Hi-C force alone, or loops + compartments without Hi-C force.\033[0m"
        )

    # TAD / region run with no structural restraint
    if (args.CHROM is not None and args.CHROM != "") and not args.LE_USE_HARMONIC_BOND and not args.HIC_USE_FORCE:
        logger.warning(
            "\033[93mRegion/TAD simulation with neither loop extrusion nor Hi-C contact force active. "
            "The polymer will fold only under generic excluded-volume and backbone forces. "
            "Consider enabling LE_USE_HARMONIC_BOND (LOOPS_PATH) or HIC_USE_FORCE (HIC_PATH).\033[0m"
        )

    # Genome-wide or chromosome-wide run without compartment annotations
    if (args.CHROM is None or args.CHROM == "") and (args.COMPARTMENT_PATH is None or args.COMPARTMENT_PATH == ""):
        if args.SCB_USE_SUBCOMPARTMENT_BLOCKS or args.COB_USE_COMPARTMENT_BLOCKS:
            pass  # already caught above as an error
        else:
            logger.warning(
                "\033[93mChromosome-wide or genome-wide simulation without compartment data. "
                "A/B compartment organisation will not be reproduced. "
                "Supply COMPARTMENT_PATH (.bed, Calder format) to enable compartment forces.\033[0m"
            )

    if args.CHB_USE_CHROMOSOMAL_BLOCKS and args.CHROM is not None and args.CHROM != "":
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


def _check_known_config_fields(args_cp, config_path: str):
    """Raise a clear, actionable ValueError if `args_cp` (the flattened
    name/value pairs read from a config.ini file) contains any key that
    SimulationConfig doesn't define — instead of silently dropping it, which
    is exactly how a typo'd or outdated/removed field name used to go
    unnoticed. For each unrecognised key, suggests the closest actual field
    name (difflib) when one is a plausible typo match, so "HIC_BOLTZMAN_ALPHA"
    points straight at "HIC_BOLTZMANN_ALPHA" instead of a generic example.
    """
    valid_fields = sorted(SimulationConfig.model_fields)
    valid_fields_upper = set(valid_fields)
    unknown = sorted({name.upper() for name, _ in args_cp if name.upper() not in valid_fields_upper})
    if not unknown:
        return

    lines = []
    for key in unknown:
        matches = difflib.get_close_matches(key, valid_fields, n=2, cutoff=0.6)
        if matches:
            lines.append(f"  - {key}  (did you mean: {' / '.join(matches)}?)")
        else:
            lines.append(f"  - {key}  (no close match found — this field does not exist)")

    raise ValueError(
        f"Unrecognized argument(s) in config file {config_path!r} — these keys do "
        f"not exist in SimulationConfig (src/multimm/config.py):\n"
        + "\n".join(lines)
        + "\nFix the typo, remove the key, or check the README for the current "
        "field name (e.g. the old HIC_ALPHA/HIC_THRESHOLD were renamed to "
        "HIC_BOLTZMANN_ALPHA / removed)."
    )


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
        # A mistyped/missing path must not pass silently: ConfigParser.read()
        # does NOT raise for a file that doesn't exist — it just returns an
        # empty list — so without this check the run would silently fall
        # back to nothing but class defaults, with no indication the
        # requested config file was never actually loaded.
        if not os.path.isfile(args_ap.config_file):
            resolved = os.path.abspath(args_ap.config_file)
            if os.path.isdir(args_ap.config_file):
                reason = "it's a directory, not a file"
            elif os.path.exists(args_ap.config_file):
                reason = "it exists but isn't a regular file"
            else:
                reason = "no such file"
            raise FileNotFoundError(
                f"Config file not found: {args_ap.config_file!r} "
                f"(resolved to {resolved!r}) — {reason}. "
                f"Check the path passed to -c/--config_file."
            )

        # inline_comment_prefixes: without it, ConfigParser treats a trailing
        # "; comment" on the SAME line as a value as part of that value
        # (only full-line comments are stripped by default) — silently
        # corrupting any "KEY = value  ; note" style line.
        config_parser = configparser.ConfigParser(inline_comment_prefixes=(";", "#"))
        read_ok = config_parser.read(args_ap.config_file)
        if not read_ok:
            # Exists but couldn't be read as INI (permissions, encoding, …) —
            # genuine parse errors (bad syntax) already raise their own
            # configparser.Error with file/line info; this covers the rest.
            raise ValueError(
                f"Config file {args_ap.config_file!r} exists but could not be "
                f"parsed as an INI file (check permissions/encoding)."
            )
        args_cp = my_config_parser(config_parser)

        # Fail fast on any config.ini key that SimulationConfig doesn't define
        # (typo, removed/renamed field, etc.) instead of silently dropping it —
        # see SimulationConfig's `extra="forbid"` for the matching check on
        # kwargs passed directly in Python. Raises with a did-you-mean
        # suggestion per unrecognised key (see _check_known_config_fields).
        _check_known_config_fields(args_cp, args_ap.config_file)

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

        # Wire the Python logging FileHandler to the output directory so that
        # all logger.info/warning/debug calls are saved to multimm.log in
        # addition to the stdout Tee below.
        setup_logger(log_file=os.path.join(log_dir, "multimm.log"))

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
                        width = len(str(args.N_ENSEMBLE - 1))
                        run_path = os.path.join(name, f"run_{i:0{width}d}")
                        args.OUT_PATH = run_path
                        
                        os.makedirs(run_path, exist_ok=True)

                        md = MultiMM(args)
                        md.run()
                        
                        # archive_run(run_path)

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