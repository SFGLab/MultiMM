#########################################################################
########### CREATOR: SEBASTIAN KORSAK, WARSAW 2024 ######################
#########################################################################
import argparse
import configparser
import os
import sys
from .args_definition import *
from .model import *

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
        """Set argument value in both attribute and internal argument list."""
        if hasattr(self.args, name):
            setattr(self.args, name, value)
        try:
            self.args.get_arg(name).val = value
        except AttributeError:
            print(f"Warning: Argument '{name}' not found in args object.")

    def convenient_argument_changer(self):
        self.set_arg('NUC_DO_INTERPOLATION', False)
        self.set_arg('ATACSEQ_PATH', None)

        modelling_level = self.args.MODELLING_LEVEL

        if str(modelling_level).lower() in ('gene'):
            print('\033[91m' + 'MAGIC COMMENT: For gene level it is needed to provide a loops_path, and a gene_name or gene_id to specify the target gene of interest.' + '\033[0m')
            self.set_arg('N_BEADS', 1000)
            self.set_arg('SC_USE_SPHERICAL_CONTAINER', False)
            self.set_arg('CHB_USE_CHROMOSOMAL_BLOCKS', False)
            self.set_arg('SCB_USE_SUBCOMPARTMENT_BLOCKS', False)
            self.set_arg('COB_USE_COMPARTMENT_BLOCKS', False)
            self.set_arg('IBL_USE_B_LAMINA_INTERACTION', False)
            self.set_arg('CF_USE_CENTRAL_FORCE', False)
            self.set_arg('SHUFFLE_CHROMS', False)
            self.set_arg('SIM_RUN_MD', True)
            self.set_arg('SIM_N_STEPS', 10000)

        elif str(modelling_level).lower() in ('region','loc'):
            print('\033[91m' + 'MAGIC COMMENT: For chromosome level it is needed to provide a loops_path. Do not forget to specify the beginning and end of your chromosome. You can remove the centromers or telomers that are in the boundaries. You can optionally add an compartment_path to include block-copolymer forces.' + '\033[0m')
            self.set_arg('N_BEADS', 5000)
            self.set_arg('SC_USE_SPHERICAL_CONTAINER', False)
            self.set_arg('CHB_USE_CHROMOSOMAL_BLOCKS', False)
            self.set_arg('SCB_USE_SUBCOMPARTMENT_BLOCKS', False)
            self.set_arg('COB_USE_COMPARTMENT_BLOCKS', self.args.COMPARTMENT_PATH != '' and str(self.args.COMPARTMENT_PATH).lower() != 'none')
            self.set_arg('IBL_USE_B_LAMINA_INTERACTION', False)
            self.set_arg('CF_USE_CENTRAL_FORCE', False)
            self.set_arg('SIM_RUN_MD', True)
            self.set_arg('SIM_N_STEPS', 10000)
        
        elif str(modelling_level).lower() in ('chromosome', 'chrom'):
            print('\033[91m' + 'MAGIC COMMENT: For chromosome level it is needed to provide a loops_path. Do not forget to specify the beginning and end of your chromosome. You can remove the centromers or telomers that are in the boundaries. You can optionally add an compartment_path to include block-copolymer forces.' + '\033[0m')
            self.set_arg('N_BEADS', 20000)
            self.set_arg('SC_USE_SPHERICAL_CONTAINER', False)
            self.set_arg('CHB_USE_CHROMOSOMAL_BLOCKS', False)
            self.set_arg('SCB_USE_SUBCOMPARTMENT_BLOCKS', False)
            self.set_arg('COB_USE_COMPARTMENT_BLOCKS', self.args.COMPARTMENT_PATH != '' and str(self.args.COMPARTMENT_PATH).lower() != 'none')
            self.set_arg('IBL_USE_B_LAMINA_INTERACTION', False)
            self.set_arg('CF_USE_CENTRAL_FORCE', False)
            self.set_arg('SIM_RUN_MD', True)
            self.set_arg('SIM_N_STEPS', 10000)
            self.set_arg('LOC_START', 1)
            self.set_arg('LOC_END', self.chrom_sizes[self.args.CHROM])
        
        elif str(modelling_level).lower() in ('gw', 'genome'):
            print('\033[91m' + 'MAGIC COMMENT: For gw level it is needed to provide a loops_path. You can optionally add an compartment_path to include block-copolymer forces.' + '\033[0m')
            self.set_arg('N_BEADS', 200000)
            self.set_arg('SC_USE_SPHERICAL_CONTAINER', True)
            self.set_arg('CHB_USE_CHROMOSOMAL_BLOCKS', False)
            self.set_arg('SCB_USE_SUBCOMPARTMENT_BLOCKS', False)
            self.set_arg('COB_USE_COMPARTMENT_BLOCKS', self.args.COMPARTMENT_PATH != '' and str(self.args.COMPARTMENT_PATH).lower() != 'none')
            self.set_arg('IBL_USE_B_LAMINA_INTERACTION', self.args.COMPARTMENT_PATH != '' and str(self.args.COMPARTMENT_PATH).lower() != 'none')
            self.set_arg('CF_USE_CENTRAL_FORCE', False)
            self.set_arg('SIM_RUN_MD', False)
            self.set_arg('SIM_N_STEPS', 10000)

def args_tests(args):
    if args.LOOPS_PATH==None or args.LOOPS_PATH=='':
        raise ValueError('\033[91mMultiMM cannot run without providing interactions in .bedpe format!!!\033[0m')
    elif (args.COMPARTMENT_PATH==None or args.COMPARTMENT_PATH=='') and args.COB_USE_COMPARTMENT_BLOCKS:
        raise ValueError('\033[91mYou cannot model compartments without providing a file in .bed format. Either disable COB_USE_COMPARTMENT_BLOCKS or import data from some compartment caller according to the documentation.\033[0m')
    elif args.NUC_DO_INTERPOLATION and args.ATACSEQ_PATH==None:
        raise ValueError('\033[91mYou enabled nucleosome simulation without providing nucleosome data. Either import a .bigwig file that shows nucleosome occupancy or disable NUC_DO_INTERPOLATION.\033[0m')
    elif (args.COMPARTMENT_PATH==None or args.COMPARTMENT_PATH=='') and args.SCB_USE_SUBCOMPARTMENT_BLOCKS:
        raise ValueError('\033[91mYou cannot model subcompartments without providing a file in .bed format. Either disable SCB_USE_SUBCOMPARTMENT_BLOCKS or import data from some compartment caller according to the documentation.\033[0m')
    elif args.COMPARTMENT_PATH==None and args.IBL_USE_B_LAMINA_INTERACTION:
        raise ValueError('\033[91mLamina interactions are compartment specific, but you did not provide a .bed file for compartments. Maybe you should disable the IBL_USE_B_LAMINA_INTERACTION?\033[0m')
    elif args.IBL_USE_B_LAMINA_INTERACTION and not (args.SCB_USE_SUBCOMPARTMENT_BLOCKS or args.COB_USE_COMPARTMENT_BLOCKS):
        raise ValueError('\033[91mYou have enabled lamina interactions which are compartment specific, but you did not enable compartment or subcompartment forces. Please, read the documentation and the paper to understand better the forcefield!\033[0m')
    elif args.CF_USE_CENTRAL_FORCE and args.CHROM!=None:
        raise ValueError('\033[91mOoo00ops! You enabled chromosome-specific attraction to the nucleolus, but you want to model only one chromosome. Maybe disable CF_USE_CENTRAL_FORCE?')
    elif args.CHB_USE_CHROMOSOMAL_BLOCKS and args.CHROM!=None:
        raise ValueError('\033[91mBetter disable CHB_USE_CHROMOSOMAL_BLOCKS when you model only one chromosome.')

    if args.SHUFFLE_CHROMS and (args.CHROM!=None and args.CHROM!=''):
        print('\n\033[38;5;214mWarning!! You enabled chromosome shuffling, but you model only a specific region of a specific chromosome.\033[0m\n')
    if args.CHROM!=None and args.IBL_USE_B_LAMINA_INTERACTION:
        print('\n\033[38;5;214mWarning!! You enabled lamina interactions, but you want to model a specific chromosomal region. It is not imprtantly wrong, but keep in mind that it makes more sense when you model the whole genome.\033[0m\n')
    if args.CHROM!=None and args.SC_USE_SPHERICAL_CONTAINER:
        print('\n\033[38;5;214mWarning!! You enabled spherical container but you want to model a single chromosomal region. It is not importantly wrong, but it makes more sense when you model the whole genome.\033[0m\n')
    if (not args.POL_USE_HARMONIC_BOND) or (not args.POL_USE_HARMONIC_ANGLE) or (not args.EV_USE_EXCLUDED_VOLUME):
        print('\n\033[38;5;214mWarning!! Take care when you disable fundamental forces from the backbone!.\033[0m\n')
    if args.CHB_USE_CHROMOSOMAL_BLOCKS:
        print('\n\033[38;5;214mWarning!! You are using chromosomal block forces. Take care because they are not always very biological. Refer to the documentation to be sure that you are doing everything correctly.\033[0m\n')

def my_config_parser(config_parser: configparser.ConfigParser):
    """Helper function that makes flat list arg name, and it's value from ConfigParser object."""
    sections = config_parser.sections()
    all_nested_fields = [dict(config_parser[s]) for s in sections]
    args_cp = []
    for section_fields in all_nested_fields:
        for name, value in section_fields.items():
            args_cp.append((name, value))
    return args_cp

def get_config():
    """Prepare list of arguments.
    First, defaults are set.
    Then, optionally config file values.
    Finally, CLI arguments overwrite everything. Then internal changes are applied."""

    print("Reading config...")

    # Step 1: Setup argparse
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('-c', '--config_file', help="Specify config file (ini format)", metavar="FILE")

    for arg in args:
        arg_parser.add_argument(f"--{arg.name.lower()}", help=arg.help)

    args_ap = arg_parser.parse_args()  # parse command-line arguments
    args_dict = vars(args_ap)

    # Step 2: If config file provided, parse it
    if args_ap.config_file:
        config_parser = configparser.ConfigParser()
        config_parser.read(args_ap.config_file)
        args_cp = my_config_parser(config_parser)

        # Override default args with values from config file
        for cp_arg in args_cp:
            name, value = cp_arg
            arg = args.get_arg(name)
            arg.val = value

    # Step 3: Override again with CLI arguments (if present)
    for name, value in args_dict.items():
        if name == "config_file":
            continue
        if value is not None:
            arg = args.get_arg(name.upper())
            arg.val = value
    
    # Step 4: Finalize parsing
    args.to_python()

    # Assuming chrom_sizes is already available globally
    changer = ArgumentChanger(args, chrom_sizes)
    changer.convenient_argument_changer()

    # Step 5: Save final arguments to config file
    write_config(args)

    return args

def write_config(args):
    """Write the automatically generated config to the metadata directory."""
    metadata_dir = os.path.join(args.OUT_PATH, 'metadata')
    os.makedirs(metadata_dir, exist_ok=True)
    config_path = os.path.join(metadata_dir, 'config_auto.ini')

    config = configparser.ConfigParser()
    for arg in args:
        config['DEFAULT'][arg.name] = str(arg.val)

    with open(config_path, 'w') as config_file:
        config.write(config_file)

    print(f"Configuration saved to {config_path}")
            
def main():
    try:
        # Input data
        args = get_config()
        args_tests(args)

        # Create output directory if it doesn't exist
        log_dir = os.path.join(args.OUT_PATH, 'metadata')
        os.makedirs(log_dir, exist_ok=True)

        # Redirect stdout and stderr to both terminal and file safely
        log_path = os.path.join(log_dir, 'output.log')
        with open(log_path, 'w') as log_file:
            original_stdout = sys.stdout
            original_stderr = sys.stderr
            sys.stdout = Tee(original_stdout, log_file)
            sys.stderr = Tee(original_stderr, log_file)

            # Run simulation
            name = args.OUT_PATH
            if args.GENERATE_ENSEMBLE:
                for i in range(args.N_ENSEMBLE):
                    args.SHUFFLING_SEED = i
                    args.OUT_PATH = name + f'_{i+1}'
                    md = MultiMM(args)
                    md.run()
            else:
                md = MultiMM(args)
                md.run()

            # Restore original streams before exiting 'with'
            sys.stdout = original_stdout
            sys.stderr = original_stderr

        # Normal exit
        sys.exit(0)

    except Exception as e:
        # Print error to original stderr and exit with code 1
        print(f"\033[91mERROR: {e}\033[0m", file=sys.stderr)
        sys.exit(1)

if __name__=='__main__':
    main()