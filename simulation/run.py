#########################################################################
########### CREATOR: SEBASTIAN KORSAK, WARSAW 2024 ######################
#########################################################################
import argparse
import configparser
from typing import List
from sys import stdout
from .args_definition import *
from .model import *

def args_tests(args):
    if args.COMPARTMENT_PATH==None and args.COB_USE_COMPARTMENT_BLOCKS:
        raise InterruptedError('\033[91mYou cannot model compartments without providing a file in .bed format. Either disable COB_USE_COMPARTMENT_BLOCKS or import data from some compartment caller according to the documentation.\033[0m')
    elif args.NUC_DO_INTERPOLATION and args.ATACSEQ_PATH==None:
        raise InterruptedError('\033[91mYou enabled nucleosome simulation without providing nucleosome data. Either import a .bigwig file that shows nucleosome occupancy or disable NUC_DO_INTERPOLATION.\033[0m')
    elif args.COMPARTMENT_PATH==None and args.SCB_USE_SUBCOMPARTMENT_BLOCKS:
        raise InterruptedError('\033[91mYou cannot model subcompartments without providing a file in .bed format. Either disable SCB_USE_SUBCOMPARTMENT_BLOCKS or import data from some compartment caller according to the documentation.\033[0m')
    elif args.COMPARTMENT_PATH==None and args.IBL_USE_B_LAMINA_INTERACTION:
        raise InterruptedError('\033[91mLamina interactions are compartment specific, but you did not provide a .bed file for compartments. Maybe you should disable the IBL_USE_B_LAMINA_INTERACTION?\033[0m')
    elif args.IBL_USE_B_LAMINA_INTERACTION and not (args.SCB_USE_SUBCOMPARTMENT_BLOCKS or COB_USE_COMPARTMENT_BLOCKS):
        raise InterruptedError('\033[91mYou have enabled lamina interactions which are compartment specific, but you did not enable compartment or subcompartment forces. Please, read the documentation and the paper to understand better the forcefield!\033[0m')
    elif args.CF_USE_CENTRAL_FORCE and args.CHROM!=None:
        raise InterruptedError('\033[91mOoo00ops! You enabled chromosome attraction, but you want to model only one chromosome. Maybe disable CF_USE_CENTRAL_FORCE?')
    elif args.CHB_USE_CHROMOSOMAL_BLOCKS and args.CHROM!=None:
        raise InterruptedError('\033[91mBetter disable CHB_USE_CHROMOSOMAL_BLOCKS when you model only one chromosome.')

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

def my_config_parser(config_parser: configparser.ConfigParser) -> List[tuple[str, str]]:
    """Helper function that makes flat list arg name, and it's value from ConfigParser object."""
    sections = config_parser.sections()
    all_nested_fields = [dict(config_parser[s]) for s in sections]
    args_cp = []
    for section_fields in all_nested_fields:
        for name, value in section_fields.items():
            args_cp.append((name, value))
    return args_cp

def get_config() -> ListOfArgs:
    """This function prepares the list of arguments.
    At first List of args with defaults is read.
    Then it's overwritten by args from config file (ini file).
    In the end config is overwritten by argparse options."""

    print(f"Reading config...")
    arg_parser = argparse.ArgumentParser()

    arg_parser.add_argument('-c', '--config_file', help="Specify config file (ini format)", metavar="FILE")
    for arg in args:
        arg_parser.add_argument(f"--{arg.name.lower()}", help=arg.help)
    args_ap = arg_parser.parse_args()  # args from argparse
    config_parser = configparser.ConfigParser()
    config_parser.read(args_ap.config_file)
    args_cp = my_config_parser(config_parser)
    # Override defaults args with values from config file
    for cp_arg in args_cp:
        name, value = cp_arg
        arg = args.get_arg(name)
        arg.val = value
    # Now again override args with values from command line.
    for ap_arg in args_ap.__dict__:
        if ap_arg not in ['config_file']:
            name, value = ap_arg, getattr(args_ap, ap_arg)
            if value is not None:
                arg = args.get_arg(name)
                arg.val = value
    args.to_python()
    args.write_config_file()
    return args
            
def main():
    # Input data
    args = get_config()
    args_tests(args)
    
    # Run simulation
    name = args.OUT_PATH
    if args.GENERATE_ENSEMBLE:
        for i in range(args.N_ENSEMBLE):
            args.SHUFFLING_SEED = i
            args.OUT_PATH = name+f'_{i+1}'
            md = MultiMM(args)
            md.run()
    else:
        md = MultiMM(args)
        md.run()

if __name__=='__main__':
    main()