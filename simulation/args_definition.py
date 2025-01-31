import datetime
import re
from dataclasses import dataclass
from math import pi
from typing import Union
import argparse
import importlib.resources
import openmm as mm
from openmm.unit import Quantity

# Dynamically set the default path to the XML file in the package
try:
    with importlib.resources.path('loopsage.forcefields', 'classic_sm_ff.xml') as default_xml_path:
        default_xml_path = str(default_xml_path)
except FileNotFoundError:
    # If running in a development setup without the resource installed, fallback to a relative path
    default_xml_path = 'simulation/forcefields/ff.xml'

@dataclass
class Arg(object):
    name: str
    help: str
    type: type
    default: Union[str, float, int, bool, Quantity, None]
    val: Union[str, float, int, bool, Quantity, None]

# Define custom type to parse list from string
def parse_list(s):
    try:
        return [int(x.strip()) for x in s.strip('[]').split(',')]
    except ValueError:
        raise argparse.ArgumentTypeError("Invalid list format. Must be a comma-separated list of integers.")

class ListOfArgs(list):
    quantity_regexp = re.compile(r'(?P<value>[-+]?\d+(?:\.\d+)?) ?(?P<unit>\w+)')

    def get_arg(self, name: str) -> Arg:
        """Stupid arg search in list of args"""
        name = name.upper()
        for i in self:
            if i.name == name:
                return i
        raise ValueError(f"No such arg: {name}")

    def __getattr__(self, item):
        return self.get_arg(item).val

    def parse_args(self):
        parser = argparse.ArgumentParser()
        for arg in self.arg_list:
            parser.add_argument(arg['name'], help=arg['help'], type=arg.get('type', str), default=arg.get('default', ''), val=arg.get('val', ''))

        args = parser.parse_args()
        parsed_args = {arg['name']: getattr(args, arg['name']) for arg in self.arg_list}
        return parsed_args

    def parse_quantity(self, val: str) -> Union[Quantity, None]:
        if val == '':
            return None
        match_obj = self.quantity_regexp.match(val)
        value, unit = match_obj.groups()
        try:
            unit = getattr(mm.unit, unit)
        except AttributeError:
            raise ValueError(f"I Can't recognise unit {unit} in expression {val}. Example of valid quantity: 12.3 femtosecond.")
        return Quantity(value=float(value), unit=unit)

    def to_python(self):
        """Casts string args to ints, floats, bool..."""
        for i in self:
            if i.val == '':
                i.val = None
            elif i.name == "HR_K_PARAM":  # Workaround for complex unit
                i.val = Quantity(float(i.val), mm.unit.kilojoule_per_mole / mm.unit.nanometer ** 2)
            elif i.type == str:
                continue
            elif i.type == int:
                i.val = int(i.val)
            elif i.type == float:
                i.val = float(i.val)
            elif i.type == bool:
                if i.val.lower() in ['true', '1', 'y', 'yes']:
                    i.val = True
                elif i.val.lower() in ['false', '0', 'n', 'no']:
                    i.val = False
                else:
                    raise ValueError(f"Can't convert {i.val} into bool type.")
            elif i.type == Quantity:
                try:
                    i.val = self.parse_quantity(i.val)
                except AttributeError:
                    raise ValueError(f"Can't parse: {i.name} = {i.val}")
            else:
                raise ValueError(f"Can't parse: {i.name} = {i.val}")

    def get_complete_config(self) -> str:
        w = "####################\n"
        w += "#   MultiMM Model   #\n"
        w += "####################\n\n"
        w += "# This is automatically generated config file.\n"
        w += f"# Generated at: {datetime.datetime.now().isoformat()}\n\n"
        w += "# Notes:\n"
        w += "# Some fields require units. Units are represented as objects from mm.units module.\n"
        w += "# Simple units are parsed directly. For example: \n"
        w += "# HR_R0_PARAM = 0.2 nanometer\n"
        w += "# But more complex units does not have any more sophisticated parser written, and will fail.'\n"
        w += "# In such cases the unit is fixed (and noted in comment), so please convert complex units manually if needed.\n"
        w += "# <float> and <int> types does not require any unit. Quantity require unit.\n\n"
        w += "# Default values does not mean valid value. In many places it's only a empty field that need to be filled.\n\n"

        w += '[Main]'
        for i in self:
            w += f'; {i.help}, type: {i.type.__name__}, default: {i.default}\n'
            if i.val is None:
                w += f'{i.name} = \n\n'
            else:
                if i.type == Quantity:
                    # noinspection PyProtectedMember
                    w += f'{i.name} = {i.val._value} {i.val.unit.get_name()}\n\n'
                else:
                    w += f'{i.name} = {i.val}\n\n'
        w = w[:-2]
        return w

    def write_config_file(self):
        auto_config_filename = 'config_auto.ini'
        with open(auto_config_filename, 'w') as f:
            f.write(self.get_complete_config())
        print(f"Automatically generated config file saved in {auto_config_filename}")


available_platforms = [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]

args = ListOfArgs([
    # Platform settings
    Arg('PLATFORM', help=f"name of the platform. Available choices: {' '.join(available_platforms)}", type=str, default='', val=''),
    Arg('DEVICE', help="device index for CUDA or OpenCL (count from 0)", type=str, default='', val=''),
    
    # Input data
    Arg('INITIAL_STRUCTURE_PATH', help="Path to CIF file.", type=str, default='', val=''),
    Arg('BUILD_INITIAL_STRUCTURE', help="To build a new initial structure.", type=bool, default='True', val='True'),
    Arg('INITIAL_STRUCTURE_TYPE', help="you can choose between: hilbert, circle, rw, confined_rw, knot, self_avoiding_rw, spiral, sphere.", type=str, default='hilbert', val='hilbert'),
    Arg('GENERATE_ENSEMBLE', help="Default value: false. True in case that you would like to have an ensemble of structures instead of one. Better to disable it for large simulations that require long computational time. Moreover it is better to start random walk initial structure in case of true value.", type=bool, default='False', val='False'),
    Arg('N_ENSEMBLE', help="Number of samples of structures that you would like to calculate.", type=int, default='', val=''),
    Arg('DOWNSAMPLING_PROB', help="Probability of downsampling contacts (from 0 to 1).", type=float, default='1.0', val='1.0'),
    Arg('FORCEFIELD_PATH', help="Path to XML file with forcefield.", type=str, default=default_xml_path, val=default_xml_path),
    Arg('N_BEADS', help="Number of Simulation Beads.", type=int, default='50000', val='50000'),
    Arg('COMPARTMENT_PATH', help="It should be a .bed file with subcompartments from Calder (or something in the same format).", type=str, default='', val=''),
    Arg('LOOPS_PATH', help="A .bedpe file path with loops. It is required.", type=str, default='', val=''),
    Arg('ATACSEQ_PATH', help="A .bw or .BigWig file path with atacseq data. It is not required.", type=str, default='', val=''),
    Arg('OUT_PATH', help="Output folder name.", type=str, default='results', val='results'),
    Arg('LOC_START', help="Starting region coordinate.", type=int, default='', val=''),
    Arg('LOC_END', help="Ending region coordinate.", type=int, default='', val=''),
    Arg('CHROM', help="Chromosome that corresponds the the modelling region of interest (in case that you do not want to model the whole genome).", type=str, default='', val=''),
    Arg('SHUFFLE_CHROMS', help="Shuffle the chromosomes.", type=bool, default='False', val='False'),
    Arg('SHUFFLING_SEED', help="Shuffling random seed.", type=int, default='0', val='0'),
    Arg('SAVE_PLOTS', help='Save plots.', type=bool, default='True', val='True'),
    
    # Basic Polymer Bond
    Arg('POL_USE_HARMONIC_BOND', help="Use harmonic bond interaction.", type=bool, default='True', val='True'),
    Arg('POL_HARMONIC_BOND_R0', help="harmonic bond distance equilibrium constant", type=float, default='0.1', val='0.1'),
    Arg('POL_HARMONIC_BOND_K', help="harmonic bond force constant (fixed unit: kJ/mol/nm^2)", type=float, default='300000.0', val='300000.0'),

    # Basic polymer stiffness
    Arg('POL_USE_HARMONIC_ANGLE', help="Use harmonic angle interaction.", type=bool, default='True', val='True'),
    Arg('POL_HARMONIC_ANGLE_R0', help="harmonic angle distance equilibrium constant", type=float, default=str(pi), val=str(pi)),
    Arg('POL_HARMONIC_ANGLE_CONSTANT_K', help="harmonic angle force constant (fixed unit: kJ/mol/radian^2)", type=float, default='100.0', val='100.0'),

    # Long-Range loop bonds
    Arg('LE_USE_HARMONIC_BOND', help="Use harmonic bond interaction for long range loops.", type=bool, default='True', val='True'),
    Arg('LE_FIXED_DISTANCES', help="For fixed distances between loops. False if you want to correlate with the hatmap strength.", type=bool, default='False', val='False'),
    Arg('LE_HARMONIC_BOND_R0', help="harmonic bond distance equilibrium constant", type=float, default='0.1', val='0.1'),
    Arg('LE_HARMONIC_BOND_K', help="harmonic bond force constant (fixed unit: kJ/mol/nm^2)", type=float, default='30000.0', val='30000.0'),

    # Excluded Volume
    Arg('EV_USE_EXCLUDED_VOLUME', help="Use excluded volume.", type=bool, default='True', val='True'),
    Arg('EV_EPSILON', help="Epsilon parameter.", type=float, default='100.0', val='100.0'),
    Arg('EV_R_SMALL', help="Add something small in denominator to make it not exploding all the time.", type=float, default='0.05', val='0.05'),
    Arg('EV_POWER', help="Power in the exponent of EV potential.", type=float, default='3.0', val='3.0'),
    
    # Spherical container
    Arg('SC_USE_SPHERICAL_CONTAINER', help='Use Spherical container', type=bool, default='False', val='False'),
    Arg('SC_RADIUS1', help='Spherical container radius,', type=float, default='', val=''),
    Arg('SC_RADIUS2', help='Spherical container radius,', type=float, default='', val=''),
    Arg('SC_SCALE', help='Spherical container scaling factor', type=float, default='1000.0', val='1000.0'),

    # Chromosomal Blocks
    Arg('CHB_USE_CHROMOSOMAL_BLOCKS', help='Use Chromosomal Blocks.', type=bool, default='False', val='False'),
    Arg('CHB_KC', help='Block copolymer width parameter.', type=float, default='0.3', val='0.3'),
    Arg('CHB_DE', help='Energy factor for block copolymer chromosomal model.', type=float, default='1e-5', val='1e-5'),
    
    # Compartment Blocks
    Arg('COB_USE_COMPARTMENT_BLOCKS', help='Use Compartment Blocks.', type=bool, default='False', val='False'),
    Arg('COB_DISTANCE', help='Block copolymer equilibrium distance for chromosomal blocks.', type=float, default='', val=''),
    Arg('COB_EA', help='Energy strength for A compartment.', type=float, default='1.0', val='1.0'),
    Arg('COB_EB', help='Energy strength for B compartment.', type=float, default='2.0', val='2.0'),
    
    # Subcompartment Blocks
    Arg('SCB_USE_SUBCOMPARTMENT_BLOCKS', help='Use Subcompartment Blocks.', type=bool, default='False', val='False'),
    Arg('SCB_DISTANCE', help='Block copolymer equilibrium distance for chromosomal blocks.', type=float, default='', val=''),
    Arg('SCB_EA1', help='Energy strength for A1 compartment.', type=float, default='1.0', val='1.0'),
    Arg('SCB_EA2', help='Energy strength for A2 compartment.', type=float, default='1.33', val='1.33'),
    Arg('SCB_EB1', help='Energy strength for B1 compartment.', type=float, default='1.66', val='1.66'),
    Arg('SCB_EB2', help='Energy strength for B2 compartment.', type=float, default='2.0', val='2.0'),

    # Interactions of B compartment with lamina
    Arg('IBL_USE_B_LAMINA_INTERACTION', help='Interactions of B compartment with lamina.', type=bool, default='False', val='False'),
    Arg('IBL_SCALE', help='Scaling factor for B comoartment interaction with lamina.', type=float, default='400.0', val='400.0'),
    
    # Central Force for Smaller Chromosomes Attraction to the Nucleolus
    Arg('CF_USE_CENTRAL_FORCE', help='Attraction of smaller chromosomes.', type=bool, default='False', val='False'),
    Arg('CF_STRENGTH', help='Strength of Interaction', type=float, default='100.0', val='100.0'),

    # Nucleosome interpolation
    Arg('NUC_DO_INTERPOLATION', help='Attraction of smaller chromosomes.', type=bool, default='False', val='False'),
    Arg('MAX_NUCS_PER_BEAD', help='Maximum amount of nucleosomes per single bead.', type=int, default='4', val='4'),
    Arg('NUC_RADIUS', help='The radius of the single nucleosome helix.', type=float, default='0.1', val='0.1'),
    Arg('POINTS_PER_NUC', help='The number of points that consist a nucleosome helix.', type=int, default='20', val='20'),
    Arg('PHI_NORM', help='Zig zag angle. ', type=float, default=str(pi/5), val=str(pi/5)),
    
    # Simulation parameters
    Arg('SIM_RUN_MD', help='Do you want to run MD simulation?', type=bool, default='False', val='False'),
    Arg('SIM_N_STEPS', help='Number of steps in MD simulation', type=int, default='', val=''),
    Arg('SIM_ERROR_TOLERANCE', help='Error tolerance for variable MD simulation', type=float, default='0.01', val='0.01'),
    Arg('SIM_AMD_ALPHA', help='Alpha of AMD simulation.', type=float, default='100.0', val='100.0'),
    Arg('SIM_AMD_E', help='E (energy) of AMD simulation.', type=float, default='1000.0', val='1000.0'),
    Arg('SIM_SAMPLING_STEP', help='It determines in t how many steps we save a structure.', type=int, default='100', val='100'),
    Arg('SIM_INTEGRATOR_TYPE', help='Alternative: langevin, verlet', type=str, default='langevin', val='langevin'),
    Arg('SIM_INTEGRATOR_STEP', help='The step of integrator.', type=Quantity, default='1 femtosecond', val='1 femtosecond'),
    Arg('SIM_FRICTION_COEFF', help='Friction coefficient (Used only with langevin integrator)', type=float, default='0.5', val='0.5'),
    Arg('SIM_SET_INITIAL_VELOCITIES', help='Sets initial velocities based on Boltzmann distribution', type=bool, default='False', val='False'),
    Arg('SIM_TEMPERATURE', help='Simulation temperature', type=Quantity, default='310 kelvin', val='310 kelvin'),

    # Trajectory settings
    Arg('TRJ_FRAMES', help='Number of trajectory frames to save.', type=int, default='2000', val='2000'),
])
