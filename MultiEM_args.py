import datetime
import re
from dataclasses import dataclass
from math import pi
from typing import Union

import simtk
import simtk.openmm as mm
from simtk.unit import Quantity


@dataclass
class Arg(object):
    name: str
    help: str
    type: type
    default: Union[str, float, int, bool, Quantity, None]
    val: Union[str, float, int, bool, Quantity, None]


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

    def parse_quantity(self, val: str) -> Union[Quantity, None]:
        if val == '':
            return None
        match_obj = self.quantity_regexp.match(val)
        value, unit = match_obj.groups()
        try:
            unit = getattr(simtk.unit, unit)
        except AttributeError:
            raise ValueError(f"I Can't recognise unit {unit} in expression {val}. Example of valid quantity: 12.3 femtosecond.")
        return Quantity(value=float(value), unit=unit)

    def to_python(self):
        """Casts string args to ints, floats, bool..."""
        for i in self:
            if i.val == '':
                i.val = None
            elif i.name == "HR_K_PARAM":  # Workaround for complex unit
                i.val = Quantity(float(i.val), simtk.unit.kilojoule_per_mole / simtk.unit.nanometer ** 2)
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
        w += "#   MultiEM Model   #\n"
        w += "####################\n\n"
        w += "# This is automatically generated config file.\n"
        w += f"# Generated at: {datetime.datetime.now().isoformat()}\n\n"
        w += "# Notes:\n"
        w += "# Some fields require units. Units are represented as objects from simtk.units module.\n"
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

# Every single arguments must be listed here.
# Not all of them must be provided by user.
# Invalid arguments should rise ValueError.
# Default args are overwritten by config.ini, and then they are overwritten by command line.
# Defaults value must be strings. They will be converted to python object later when ListOfArgs.to_python() will be called
args = ListOfArgs([
    # Platform settings
    Arg('PLATFORM', help=f"name of the platform. Available choices: {' '.join(available_platforms)}", type=str, default='', val=''),
    Arg('DEVICE', help="device index for CUDA or OpenCL (count from 0)", type=str, default='', val=''),

    # Input data
    Arg('INITIAL_STRUCTURE_PATH', help="Path to CIF file.", type=str, default='', val=''),
    Arg('FORCEFIELD_PATH', help="Path to XML file with forcefield.", type=str, default='', val=''),
    Arg('N_BEADS', help="Number of Simulation Beads.", type=int, default=50000, val=50000),
    Arg('COMPARTMENT_PATH', help="It can be either .bed file with subcompartments from Calder or .BigWig signal.", type=str, default='', val=''),
    Arg('LOOPS_PATH', help="A .bedpe file path with loops. It is required.", type=str, default='', val=''),

    # Basic Polymer Bond
    Arg('POL_USE_HARMONIC_BOND', help="Use harmonic bond interaction.", type=bool, default='True', val='True'),
    Arg('POL_HARMONIC_BOND_R0', help="harmonic bond distance equilibrium constant", type=Quantity, default='0.1 nanometer', val='0.1 nanometer'),
    Arg('POL_HARMONIC_BOND_K', help="harmonic bond force constant (fixed unit: kJ/mol/nm^2)", type=float, default='300000.0', val='300000.0'),

    # Basic polymer stiffness
    Arg('POL_USE_HARMONIC_ANGLE', help="Use harmonic angle interaction.", type=bool, default='True', val='True'),
    Arg('POL_HARMONIC_ANGLE_R0', help="harmonic angle distance equilibrium constant", type=float, default=str(pi), val=str(pi)),
    Arg('POL_HARMONIC_ANGLE_K', help="harmonic angle force constant (fixed unit: kJ/mol/radian^2)", type=float, default='20.0', val='20.0'),

    # Excluded Volume
    Arg('EV_USE_EXCLUDED_VOLUME', help="Use excluded volume.", type=bool, default='False', val='False'),
    Arg('EV_EPSILON', help="Epsilon parameter.", type=float, default='10.0', val='10.0'),
    Arg('EV_SIGMA', help="Sigma parameter.", type=Quantity, default='0.05 nanometer', val='0.05 nanometer'),    

    # Spherical container
    Arg('SC_USE_SPHERICAL_CONTAINER', help='Use Spherical container', type=bool, default='False', val='False'),
    Arg('SC_RADIUS1', help='Spherical container radius, fixed unit: nanometers', type=float, default=None, val=None),
    Arg('SC_RADIUS2', help='Spherical container radius, fixed unit: nanometers', type=float, default=None, val=None),
    Arg('SC_SCALE', help='Spherical container scaling factor', type=float, default='1000', val='1000'),

    # Chromosomal Blocks
    Arg('CHB_CHROMOSOMAL_BLOCKS', help='Use Chromosomal Blocks.', type=bool, default='False', val='False'),
    Arg('CHB_CHROM_DISTANCE', help='Block copolymer equilibrium distance for chromosomal blocks.', type=float, default=None, val=None),
    Arg('CHB_CHROM_dE', help='Energy factor for block copolymer chromosomal model.', type=float, default='0.1', val='0.1'),

    # Compartment Blocks
    Arg('COB_COMPARTMENT_BLOCKS', help='Use Chromosomal Blocks.', type=bool, default='False', val='False'),
    Arg('COB_DISTANCE', help='Block copolymer equilibrium distance for chromosomal blocks.', type=float, default=None, val=None),
    Arg('COB_EA', help='Energy strength for A compartment.', type=float, default='0.5', val='0.5'),
    Arg('COB_EB', help='Energy strength for B compartment.', type=float, default='2.0', val='2.0'),

    # Subcompartment Blocks
    Arg('SCB_CHROMOSOMAL_BLOCKS', help='Use Chromosomal Blocks.', type=bool, default='False', val='False'),
    Arg('SCB_CHROM_DISTANCE', help='Block copolymer equilibrium distance for chromosomal blocks.', type=float, default=None, val=None),
    Arg('COB_EA1', help='Energy strength for A1 compartment.', type=float, default='0.5', val='0.5'),
    Arg('COB_EA2', help='Energy strength for A2 compartment.', type=float, default='1.0', val='1.0'),
    Arg('COB_EB1', help='Energy strength for B1 compartment.', type=float, default='1.5', val='1.5'),
    Arg('COB_EB2', help='Energy strength for B2 compartment.', type=float, default='2.0', val='2.0'),

    # Small Chromosomes Attraction to Center
    Arg('CF_CENTRAL_FORCE', help='Attraction of smaller chromosomes.', type=bool, default='False', val='False'),
    Arg('CF_STRENGTH', help='Strength of Interaction', type=float, default='100.0', val='100.0'),
    Arg('CF_CHROM_dE', help='Spherical container scaling factor', type=float, default='1000', val='1000'),

    # Energy minimization
    Arg('MINIMIZE', help='should initial structure be minimized? - This is spring model main functionality.', type=bool, default='True', val='True'),
    Arg('MINIMIZED_FILE', help='If left empty result file will have name based on initial structure file name with _min.pdb ending.', type=str, default='', val=''),

    # Simulation parameters
    Arg('SIM_RUN_SIMULATION', help='Do you want to run MD simulation?', type=bool, default='False', val='False'),
    Arg('SIM_INTEGRATOR_TYPE', help='Alternative: langevin, verlet', type=str, default='verlet', val='verlet'),
    Arg('SIM_FRICTION_COEFF', help='Friction coefficient (Used only with langevin integrator)', type=float, default='', val=''),
    Arg('SIM_N_STEPS', help='Number of steps in MD simulation', type=int, default='', val=''),
    Arg('SIM_TIME_STEP', help='Time step (use time unit from simtk.unit module)', type=Quantity, default='', val=''),
    Arg('SIM_TEMP', help='Temperature (use temperature unit from simtk.unit module)', type=Quantity, default='', val=''),
    Arg('SIM_RANDOM_SEED', help='Random seed. Set to 0 for random seed.', type=int, default='0', val='0'),
    Arg('SIM_SET_INITIAL_VELOCITIES', help='Sets initial velocities based on Boltzmann distribution', type=bool, default='False', val='False'),
    Arg('SIM_USE_ANDERSEN_THERMOSTAT', help='Do you want to run simulation with Andersen Thermostat?', type=bool, default='False', val='False'),
    Arg('SIM_ANDERSEN_THERMOSTAT_COLLISION_FREQ_NOMINATOR', help='Number of collisions in time unit', type=int, default='1', val='1'),
    Arg('SIM_ANDERSEN_THERMOSTAT_COLLISION_FREQ_DENOMINATOR', help='Time unit in collision rate', type=Quantity, default='1 picosecond', val='1 picosecond'),

    # Trajectory settings
    Arg('TRJ_FRAMES', help='Number of trajectory frames to save.', type=int, default='2000', val='2000'),
    Arg('TRJ_FILENAME_DCD', help='Write trajectory in DCD file format, leave empty if you do not want to save.', type=str, default='', val=''),
    Arg('TRJ_FILENAME_PDB', help='Write trajectory in PDB file format, leave empty if you do not want to save.', type=str, default='', val=''),
    Arg('TRJ_LAST_FRAME_PDB', help='Write last frame of trajectory in PDB file format, leave empty if you do not want to save.', type=str, default='', val=''),

    # State reporting
    Arg('REP_STATE_N_SCREEN', help='Number of states reported on screen', type=int, default='20', val='20'),
    Arg('REP_STATE_N_FILE', help='Number of states reported to file screen', type=int, default='1000', val='1000'),
    Arg('REP_STATE_FILE_PATH', help='Filepath to save state. Leave empty if not needed.', type=str, default='state.csv', val='state.csv'),
    Arg('REP_STATE_FILE_H5_PATH', help="H5 file to save velocities. Leave empty if not needed.", type=str, default="state.h5", val='state.h5'),
    Arg('REP_PLOT_FILE_NAME', help='Filepath to save energy plot. Leave empty if not needed.', type=str, default='energy.pdf', val='energy.pdf'),
])