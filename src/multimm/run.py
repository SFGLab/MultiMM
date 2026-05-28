#########################################################################
########### CREATOR: SEBASTIAN KORSAK, WARSAW 2024 ######################
#########################################################################
import argparse
import configparser
import logging
import os
import sys
from enum import Enum

from openmm.unit import Quantity

from .config import SimulationConfig
from .model import MultiMM
from .utils import chrom_sizes

logger = logging.getLogger(__name__)


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
    def __init__(self, args, chrom_sizes):
        self.args = args
        self.chrom_sizes = chrom_sizes

    def set_arg(self, name, value):
        """Set argument value in attribute."""
        if hasattr(self.args, name):
            setattr(self.args, name, value)
        else:
            logger.warning(f"Warning: Argument '{name}' not found in args object.")

    def convenient_argument_changer(self):
        self.set_arg("NUC_DO_INTERPOLATION", False)
        self.set_arg("ATACSEQ_PATH", None)

        modelling_level = self.args.MODELLING_LEVEL
        print("The modelling level is: ", modelling_level)

        # EARLY EXIT: if nothing is defined, do not modify defaults
        if modelling_level is None or str(modelling_level).strip() == "":
            logger.warning(
                "No modelling level specified. Using user-specified parameter configuration."
            )
            return
        elif str(modelling_level).lower() in ("gene"):
            logger.warning(
                "MAGIC COMMENT: For gene level it is needed to provide a loops_path, and a gene_name or gene_id to specify the target gene of interest."
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

        elif str(modelling_level).lower() in ("region", "loc"):
            logger.warning(
                "Region-level modelling selected. MultiMM will construct a reduced system "
                "of 5000 beads focused on TAD-scale dynamics, including loop interactions only. "
                "Compartment and higher-order nuclear organization terms are disabled by default. "
                "A .bed file defining the region of interest must be provided."
            )
            self.set_arg("N_BEADS", 5000)
            self.set_arg("SC_USE_SPHERICAL_CONTAINER", False)
            self.set_arg("CHB_USE_CHROMOSOMAL_BLOCKS", False)
            self.set_arg("SCB_USE_SUBCOMPARTMENT_BLOCKS", False)
            self.set_arg(
                "COB_USE_COMPARTMENT_BLOCKS",
                bool(self.args.COMPARTMENT_PATH),
            )
            self.set_arg("IBL_USE_B_LAMINA_INTERACTION", False)
            self.set_arg("CF_USE_CENTRAL_FORCE", False)
            self.set_arg("SIM_RUN_MD", True)
            self.set_arg("SIM_N_STEPS", 10000)

        elif str(modelling_level).lower() in ("chromosome", "chrom"):
            logger.warning(
                "MAGIC COMMENT: For chromosome level it is needed to provide a loops_path."
                "Do not forget to specify the beginning and end of your chromosome." 
                "You can remove the centromers or telomers that are in the boundaries."
                "You can optionally add an compartment_path to include block-copolymer forces."
            )
            self.set_arg("N_BEADS", 20000)
            self.set_arg("SC_USE_SPHERICAL_CONTAINER", False)
            self.set_arg("CHB_USE_CHROMOSOMAL_BLOCKS", False)
            self.set_arg("SCB_USE_SUBCOMPARTMENT_BLOCKS", False)
            self.set_arg(
                "COB_USE_COMPARTMENT_BLOCKS",
                bool(self.args.COMPARTMENT_PATH),
            )
            self.set_arg("IBL_USE_B_LAMINA_INTERACTION", False)
            self.set_arg("CF_USE_CENTRAL_FORCE", False)
            self.set_arg("SIM_RUN_MD", True)
            self.set_arg("SIM_N_STEPS", 10000)
            self.set_arg("LOC_START", 1)
            self.set_arg("LOC_END", self.chrom_sizes[self.args.CHROM])

        elif str(modelling_level).lower() in ("gw", "genome"):
            logger.warning(
                "MAGIC COMMENT: For gw level it is needed to provide a loops_path. You can optionally add an compartment_path to include block-copolymer forces."
            )
            self.set_arg("N_BEADS", 200000)
            self.set_arg("SC_USE_SPHERICAL_CONTAINER", True)
            self.set_arg("CHB_USE_CHROMOSOMAL_BLOCKS", False)
            self.set_arg("SCB_USE_SUBCOMPARTMENT_BLOCKS", False)
            self.set_arg(
                "COB_USE_COMPARTMENT_BLOCKS",
                bool(self.args.COMPARTMENT_PATH),
            )
            self.set_arg(
                "IBL_USE_B_LAMINA_INTERACTION",
                bool(self.args.COMPARTMENT_PATH),
            )
            self.set_arg("CF_USE_CENTRAL_FORCE", False)
            self.set_arg("SIM_RUN_MD", False)
            self.set_arg("SIM_N_STEPS", 10000)


def args_tests(args):
    if args.LOOPS_PATH is None or args.LOOPS_PATH == "":
        raise ValueError("\033[91mMultiMM cannot run without providing interactions in .bedpe format!!!\033[0m")
    elif (args.COMPARTMENT_PATH is None or args.COMPARTMENT_PATH == "") and args.COB_USE_COMPARTMENT_BLOCKS:
        raise ValueError(
            "\033[91mYou cannot model compartments without providing a file in .bed format. Either disable COB_USE_COMPARTMENT_BLOCKS or import data from some compartment caller according to the documentation.\033[0m"
        )
    elif args.NUC_DO_INTERPOLATION and args.ATACSEQ_PATH is None:
        raise ValueError(
            "\033[91mYou enabled nucleosome simulation without providing nucleosome data. Either import a .bigwig file that shows nucleosome occupancy or disable NUC_DO_INTERPOLATION.\033[0m"
        )
    elif (args.COMPARTMENT_PATH is None or args.COMPARTMENT_PATH == "") and args.SCB_USE_SUBCOMPARTMENT_BLOCKS:
        raise ValueError(
            "\033[91mYou cannot model subcompartments without providing a file in .bed format. Either disable SCB_USE_SUBCOMPARTMENT_BLOCKS or import data from some compartment caller according to the documentation.\033[0m"
        )
    elif args.COMPARTMENT_PATH is None and args.IBL_USE_B_LAMINA_INTERACTION:
        raise ValueError(
            "\033[91mLamina interactions are compartment specific, but you did not provide a .bed file for compartments. Maybe you should disable the IBL_USE_B_LAMINA_INTERACTION?\033[0m"
        )
    elif args.IBL_USE_B_LAMINA_INTERACTION and not (
        args.SCB_USE_SUBCOMPARTMENT_BLOCKS or args.COB_USE_COMPARTMENT_BLOCKS
    ):
        raise ValueError(
            "\033[91mYou have enabled lamina interactions which are compartment specific, but you did not enable compartment or subcompartment forces. Please, read the documentation and the paper to understand better the forcefield!\033[0m"
        )
    elif args.CF_USE_CENTRAL_FORCE and args.CHROM is not None:
        raise ValueError(
            "\033[91mOoo00ops! You enabled chromosome-specific attraction to the nucleolus, but you want to model only one chromosome. Maybe disable CF_USE_CENTRAL_FORCE?"
        )
    elif args.CHB_USE_CHROMOSOMAL_BLOCKS and args.CHROM is not None:
        raise ValueError("\033[91mBetter disable CHB_USE_CHROMOSOMAL_BLOCKS when you model only one chromosome.")

    if args.SHUFFLE_CHROMS and (args.CHROM is not None and args.CHROM != ""):
        logger.warning(
            "\nWarning!! You enabled chromosome shuffling, but you model only a specific region of a specific chromosome.\n"
        )
    if args.CHROM is not None and args.IBL_USE_B_LAMINA_INTERACTION:
        logger.warning(
            "\nWarning!! You enabled lamina interactions, but you want to model a specific chromosomal region. It is not imprtantly wrong, but keep in mind that it makes more sense when you model the whole genome.\n"
        )
    if args.CHROM is not None and args.SC_USE_SPHERICAL_CONTAINER:
        logger.warning(
            "\nWarning!! You enabled spherical container but you want to model a single chromosomal region. It is not importantly wrong, but it makes more sense when you model the whole genome.\n"
        )
    if (not args.POL_USE_HARMONIC_BOND) or (not args.POL_USE_HARMONIC_ANGLE) or (not args.EV_USE_EXCLUDED_VOLUME):
        logger.warning("\nWarning!! Take care when you disable fundamental forces from the backbone!.\n")
    if args.CHB_USE_CHROMOSOMAL_BLOCKS:
        logger.warning(
            "\nWarning!! You are using chromosomal block forces. Take care because they are not always very biological. Refer to the documentation to be sure that you are doing everything correctly.\n"
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


def main():
    try:
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

            name = args.OUT_PATH
            if args.GENERATE_ENSEMBLE:
                for i in range(args.N_ENSEMBLE):
                    args.SHUFFLING_SEED = i
                    args.OUT_PATH = name + f"_{i+1}"
                    md = MultiMM(args)
                    md.run()
            else:
                md = MultiMM(args)
                md.run()

            sys.stdout = original_stdout
            sys.stderr = original_stderr

        sys.exit(0)

    except Exception as e:
        logger.error(f"ERROR: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
