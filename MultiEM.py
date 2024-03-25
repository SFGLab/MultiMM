#########################################################################
########### CREATOR: SEBASTIAN KORSAK, WARSAW 2024 ######################
#########################################################################

import random as rd
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse
import configparser
from mdtraj.reporters import HDF5Reporter
from scipy import ndimage
from typing import List
from sys import stdout
import openmm as mm
from openmm.app import PDBFile, PDBxFile, ForceField, Simulation, PDBReporter, PDBxReporter, DCDReporter, StateDataReporter, CharmmPsfFile
from MultiEM_init_tools import *
from MultiEM_utils import *
from MultiEM_plots import *
from MultiEM_args import *
import time
import os

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

class MultiEM:
    def __init__(self,args):
        '''
        Input data:
        ------------
        args: list of arguments imported from config.ini file.
        '''
        ################################################################
        ################ CREATE OUTPUT FOLDER ##########################
        ################################################################
        
        self.ms, self.ns, self.ds, self.chr_ends, self.Cs = None, None, None, None, None
        
        # Make save directory
        self.save_path = args.OUT_PATH+'/'
        try:
            os.mkdir(self.save_path)
            os.mkdir(self.save_path+'ensembles')
            os.mkdir(self.save_path+'chromosomes')
            os.mkdir(self.save_path+'plots')
        except OSError as error:
            print("Folder 'chromosomes' already exists!")

        self.args  = args
        loop_path = args.LOOPS_PATH
        comp_path = args.COMPARTMENT_PATH
        nuc_path = args.NUC_PATH
        
        ################################################################
        ################ LOAD COMPARTMENTS #############################
        ################################################################
        if comp_path.endswith('.bw') or comp_path.endswith('.bigwig') or comp_path.endswith('.BigWig'):
            self.Cs, self.chr_ends = import_bw(bw_PATH=comp_path,N_beads=self.args.N_BEADS,\
                                               viz=False,binary=True,\
                                               path=self.save_PATH)
        elif comp_path.endswith('.bed'):
            self.Cs, self.chr_ends = import_compartments_from_Calder(bed_file=comp_path,N_beads=self.args.N_BEADS,\
                                                                     chrom=self.args.CHROM,coords=self.args.COORDS,\
                                                                     save_path=self.save_path)
        else:
            self.Cs = None

        ################################################################
        ################ LOAD LOOPS ####################################
        ################################################################
        # Print the value and length of LOOPS_path for debugging
        if loop_path.endswith('.bedpe'):
            self.ms, self.ns, self.ds, self.chr_ends = import_mns_from_bedpe(bedpe_file = loop_path,N_beads=self.args.N_BEADS,\
                                                                             coords = args.COORDS, chrom=args.CHROM,\
                                                                             viz=False, path=self.save_path)
        else:
            raise InterruptedError('You did not provide appropriate loop file. Loop .bedpe file is obligatory.')
        
        write_chrom_colors(self.chr_ends,name=self.save_path+'MultiEM_chromosome_colors.cmd')
        
        ################################################################
        ################ LOAD NUCLEOSOMES ##############################
        ################################################################
        if nuc_path!=None:
            self.entry_points, self.exit_points = puffin_to_array(nucs_path,args.COORDS)
            self.nuc_sim_len, self.num_nucs = int(np.max(self.exit_points+1)), len(self.entry_points)
            self.is_nuc = np.full(self.nuc_sim_len,False)
            for i,j in zip(self.entry_points,self.exit_points):
                self.is_nuc[i-1:j+1]=True
            if self.nuc_sim_len>1e7: raise InterruptedError('Too large region or structure to model nucleosomes. Please choose smaller region or simulation beads. ;)')   
        else:
            self.entry_points, self.exit_points, self.nuc_sim_len, self.num_nucs = None, None, 0, 0
        
        # if np.all(comp_path!=None): self.Cs = align_comps(self.Cs,self.ms,self.chr_ends)
        
        ################################################################
        ################ LOAD CHROMOSOMES ##############################
        ################################################################
        self.chrom_spin = np.zeros(self.args.N_BEADS)
        if self.args.CHROM==None or self.args.CHROM=='':
            print('___________________________________________________________')
            print('A region was not specified. Whole genome simulation starts.')
            print('___________________________________________________________')
            for i in range(len(self.chr_ends)-1):
                self.chrom_spin[self.chr_ends[i]:self.chr_ends[i+1]] = i

    def add_forcefield(self):
        '''
        Here we define the forcefield of MultiEM.
        '''
        # Leonard-Jones potential for excluded volume
        if self.args.EV_USE_EXCLUDED_VOLUME:
            self.ev_force = mm.CustomNonbondedForce(f'epsilon*(sigma/r)^{self.args.EV_POWER}')
            self.ev_force.addGlobalParameter('epsilon', defaultValue=self.args.EV_EPSILON)
            self.ev_force.addGlobalParameter('sigma', defaultValue=np.min(self.ds))
            for i in range(self.system.getNumParticles()):
                self.ev_force.addParticle()
            self.system.addForce(self.ev_force)

        ## Compartment force
        if self.args.COB_USE_COMPARTMENT_BLOCKS:
            self.comp_force = mm.CustomNonbondedForce('-E*exp(-r^2/(2*rc^2)); E=(Ea*delta(s1-1)*delta(s2-1)+Eb*delta(s1+1)*delta(s2+1))')
            self.comp_force.addGlobalParameter('rc',defaultValue=self.r_comp)
            self.comp_force.addGlobalParameter('Ea',defaultValue=self.args.COB_EA)
            self.comp_force.addGlobalParameter('Eb',defaultValue=self.args.COB_EB)
            self.comp_force.addPerParticleParameter('s')
            for i in range(self.system.getNumParticles()):
                self.comp_force.addParticle([self.Cs[i]])
            self.system.addForce(self.comp_force)
        
        ## Subcompartment force
        if self.args.SCB_USE_SUBCOMPARTMENT_BLOCKS:
            self.scomp_force = mm.CustomNonbondedForce('-E*exp(-r^2/(2*rsc^2)); E=Ea1*delta(s1-2)*delta(s2-2)+Ea2*delta(s1-1)*delta(s2-1)+Eb1*delta(s1+1)*delta(s2+1)+Eb2*delta(s1+2)*delta(s2+2)')
            self.scomp_force.addGlobalParameter('rsc',defaultValue=self.r_comp)
            self.scomp_force.addGlobalParameter('Ea1',defaultValue=self.args.SCB_EA1)
            self.scomp_force.addGlobalParameter('Ea2',defaultValue=self.args.SCB_EA2)
            self.scomp_force.addGlobalParameter('Eb1',defaultValue=self.args.SCB_EB1)
            self.scomp_force.addGlobalParameter('Eb2',defaultValue=self.args.SCB_EB2)
            self.scomp_force.addPerParticleParameter('s')
            for i in range(self.system.getNumParticles()):
                self.scomp_force.addParticle([self.Cs[i]])
            self.system.addForce(self.scomp_force)

        # Chromosomal blocks
        if self.args.CHB_USE_CHROMOSOMAL_BLOCKS:
            self.chrom_block_force = mm.CustomNonbondedForce('-E*exp(-r^2/(2*r_C^2)); E=dE*delta(chrom1-chrom2)')
            self.chrom_block_force.addGlobalParameter('r_C',defaultValue=self.r_chrom)
            self.chrom_block_force.addGlobalParameter('dE',defaultValue=self.args.CHB_DE)
            self.chrom_block_force.addPerParticleParameter('chrom')
            for i in range(self.system.getNumParticles()):
                self.chrom_block_force.addParticle([self.chrom_spin[i]])
            self.system.addForce(self.chrom_block_force)

        # Spherical container
        if self.args.SC_USE_SPHERICAL_CONTAINER:
            self.container_force = mm.CustomExternalForce('C*(max(0, r-R2)^2+max(0, R1-r)^2); r=sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2)')
            self.container_force.addGlobalParameter('C',defaultValue=self.args.SC_SCALE)
            self.container_force.addGlobalParameter('R1',defaultValue=self.radius1)
            self.container_force.addGlobalParameter('R2',defaultValue=self.radius2)
            self.container_force.addGlobalParameter('x0',defaultValue=self.mass_center[0])
            self.container_force.addGlobalParameter('y0',defaultValue=self.mass_center[1])
            self.container_force.addGlobalParameter('z0',defaultValue=self.mass_center[2])
            for i in range(self.system.getNumParticles()):
                self.container_force.addParticle(i, [])
            self.system.addForce(self.container_force)
        
        # Interaction of A compartment with lamina
        if self.args.IAL_USE_A_LAMINA_INTERACTION:
            self.Alamina_force = mm.CustomExternalForce('-A*sin(pi*(r-R1)/(R2-R1))^2*(delta(s-1)+delta(s-2))*step(R1)*(1-step(R2)); r=sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2)')
            self.Alamina_force.addGlobalParameter('A',defaultValue=self.args.IAL_SCALE)
            self.Alamina_force.addGlobalParameter('pi',defaultValue=np.pi)
            self.Alamina_force.addGlobalParameter('R1',defaultValue=self.radius1)
            self.Alamina_force.addGlobalParameter('R2',defaultValue=self.radius2)
            self.Alamina_force.addGlobalParameter('x0',defaultValue=self.mass_center[0])
            self.Alamina_force.addGlobalParameter('y0',defaultValue=self.mass_center[1])
            self.Alamina_force.addGlobalParameter('z0',defaultValue=self.mass_center[2])
            self.Alamina_force.addPerParticleParameter('s')
            for i in range(self.system.getNumParticles()):
                self.Alamina_force.addParticle(i, [self.Cs[i]])
            self.system.addForce(self.Alamina_force)
            
        # Interaction of B compartment with lamina
        if self.args.IBL_USE_B_LAMINA_INTERACTION:
            if self.radius1!=0.0:
                self.Blamina_force = mm.CustomExternalForce('B*(sin(pi*(r-R1)/(R2-R1))^2-1)*(delta(s+1)+delta(s+2))*step(R1)*(1-step(R2)); r=sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2)')
            else:
                self.Blamina_force = mm.CustomExternalForce('B*(sin(pi*(r-R1)/(R2-R1))^2-1)*(delta(s+1)+delta(s+2))*step(R2/2)*(1-step(R2)); r=sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2)')
            self.Blamina_force.addGlobalParameter('B',defaultValue=self.args.IBL_SCALE)
            self.Blamina_force.addGlobalParameter('pi',defaultValue=np.pi)
            self.Blamina_force.addGlobalParameter('R1',defaultValue=self.radius1)
            self.Blamina_force.addGlobalParameter('R2',defaultValue=self.radius2)
            self.Blamina_force.addGlobalParameter('x0',defaultValue=self.mass_center[0])
            self.Blamina_force.addGlobalParameter('y0',defaultValue=self.mass_center[1])
            self.Blamina_force.addGlobalParameter('z0',defaultValue=self.mass_center[2])
            self.Blamina_force.addPerParticleParameter('s')
            for i in range(self.system.getNumParticles()):
                self.Blamina_force.addParticle(i, [self.Cs[i]])
            self.system.addForce(self.Blamina_force)

        if self.args.CF_USE_CENTRAL_FORCE:
            # Force that sets smaller chromosomes closer to the center
            self.central_force = mm.CustomExternalForce('G*chrom/23*(-1/(r-R1+1)+1/(r-R1+1)^2); r=sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2)')
            self.central_force.addGlobalParameter('G',defaultValue=self.args.CF_STRENGTH)
            self.central_force.addGlobalParameter('R1',defaultValue=self.radius1)
            self.central_force.addGlobalParameter('x0',defaultValue=self.mass_center[0])
            self.central_force.addGlobalParameter('y0',defaultValue=self.mass_center[1])
            self.central_force.addGlobalParameter('z0',defaultValue=self.mass_center[2])
            self.central_force.addPerParticleParameter('chrom')
            for i in range(self.system.getNumParticles()):
                self.central_force.addParticle(i, [self.chrom_spin[i]])
            self.system.addForce(self.central_force)

        # Bond force
        if self.args.POL_USE_HARMONIC_BOND:
            self.bond_force = mm.HarmonicBondForce()
            for i in range(self.system.getNumParticles()-1):
                if (i not in self.chr_ends): self.bond_force.addBond(i,i+1,self.args.POL_HARMONIC_BOND_R0,self.args.POL_HARMONIC_BOND_K)
            self.system.addForce(self.bond_force)
        
        # Add fixed loops force
        if self.args.LE_USE_HARMONIC_BOND:
            self.loop_force = mm.HarmonicBondForce()
            counter=0
            for m,n in tqdm(zip(self.ms,self.ns),total=len(self.ms)):
                if np.any(self.ds==None):
                    self.loop_force.addBond(m,n,self.args.LE_HARMONIC_BOND_R0,self.args.LE_HARMONIC_BOND_K)
                else:
                    self.loop_force.addBond(m,n,self.ds[counter],self.args.LE_HARMONIC_BOND_K)
                counter+=1
            self.system.addForce(self.loop_force)

        # Bending potential for stiffness
        if self.args.POL_USE_HARMONIC_ANGLE:
            self.angle_force = mm.HarmonicAngleForce()
            for i in range(self.system.getNumParticles()-2):
                if (i not in self.chr_ends) and (i not in self.chr_ends-1): self.angle_force.addAngle(i, i+1, i+2, self.args.POL_HARMONIC_ANGLE_R0, self.args.POL_HARMONIC_CONSTANT_K)
            self.system.addForce(self.angle_force)

    def run_pipeline(self):
        '''
        Energy minimization for GW model.
        '''
        # Estimation of parameters
        self.radius1 = 1*(self.args.N_BEADS/50000)**(1/3) if self.args.SC_RADIUS1==None else self.args.SC_RADIUS1
        self.radius2 = 6*(self.args.N_BEADS/50000)**(1/3) if self.args.SC_RADIUS2==None else self.args.SC_RADIUS2
        if self.args.COB_DISTANCE!=None:
            self.r_comp = self.args.COB_DISTANCE
        elif self.args.SCB_DISTANCE!=None:
            self.r_comp = self.args.SCB_DISTANCE
        else:
            self.r_comp = (self.radius2-self.radius1)/20
        self.r_chrom = self.r_comp/2 if self.args.CHB_DISTANCE==None else self.args.CHB_DISTANCE

        # Initialize simulation
        if self.args.BUILD_INITIAL_STRUCTURE:
            print('\nCreating initial structure...')
            comp_mode = 'compartments' if np.all(self.Cs!=None) and len(np.unique(self.Cs))<=3 else 'subcompartments'
            if np.all(self.Cs!=None): write_cmm(self.Cs,name=self.save_path+'MultiEM_compartment_colors.cmd')
            pdb_content = build_init_mmcif(n_dna=self.args.N_BEADS,chrom_ends=self.chr_ends,path=self.save_path,hilbert=self.args.CHROM==None,scale=(self.radius1+self.radius2)/2)
            print('---Done!---')
        pdb = PDBxFile(self.save_path+'MultiEM_init.cif') if self.args.INITIAL_STRUCTURE_path==None or build_init_mmcif else PDBxFile(self.args.INITIAL_STRUCTURE_PATH)
        self.mass_center = np.average(get_coordinates_mm(pdb.positions),axis=0)
        forcefield = ForceField(self.args.FORCEFIELD_PATH)
        self.system = forcefield.createSystem(pdb.topology)
        if args.SIM_INTEGRATOR_TYPE=='verlet':
            integrator  = mm.VerletIntegrator(self.args.SIM_TEMPERATURE,self.args.SIM_INTEGRATOR_STEP)
        else:
            integrator = mm.LangevinIntegrator(self.args.SIM_TEMPERATURE, self.args.SIM_FRICTION_COEFF, self.args.SIM_INTEGRATOR_STEP)
        
        # Import forcefield
        print('\nImporting forcefield...')
        self.add_forcefield()
        print('---Done!---')

        # Run simulation / Energy minimization
        print('\nEnergy minimization...')
        platform = mm.Platform.getPlatformByName(self.args.PLATFORM)
        simulation = Simulation(pdb.topology, self.system, integrator, platform)
        
        simulation.context.setPositions(pdb.positions)
        simulation.context.setVelocitiesToTemperature(self.args.SIM_TEMPERATURE, 0)
        current_platform = simulation.context.getPlatform()
        print(f"Simulation will run on platform: {current_platform.getName()}.")
        start_time = time.time()
        simulation.minimizeEnergy()
        state = simulation.context.getState(getPositions=True)
        PDBxFile.writeFile(pdb.topology, state.getPositions(), open(self.save_path+'MultiEM_minimized.cif', 'w'))
        print(f"--- Energy minimization done!! Executed in {(time.time() - start_time)/60:.2f} minutes. :D ---\n")        
        
        start = time.time()
        V = get_coordinates_mm(state.getPositions())
        for i in range(len(self.chr_ends)-1):
            write_mmcif_chrom(coords=10*V[self.chr_ends[i]:self.chr_ends[i+1]],\
                        path=self.save_path+f'chromosomes/MultiEM_minimized_{chrs[i]}.cif')

        if self.args.SIM_RUN_MD:
            simulation.reporters.append(StateDataReporter(stdout, self.args.SIM_SAMPLING_STEP, step=True, totalEnergy=True, kineticEnergy=True ,potentialEnergy=True, temperature=True, separator='\t'))
            simulation.reporters.append(DCDReporter(self.save_path+'MultiEM_annealing.dcd', self.args.SIM_N_STEPS//self.args.TRJ_FRAMES))
            print('Running relaxation...')
            start = time.time()
            for i in range(self.args.SIM_N_STEPS//self.args.SIM_SAMPLING_STEP):
                simulation.step(self.args.SIM_SAMPLING_STEP)
                self.state = simulation.context.getState(getPositions=True)
                PDBxFile.writeFile(pdb.topology, self.state.getPositions(), open(self.save_path+f'ensembles/ens_{i+1}.cif', 'w'))
            end = time.time()
            elapsed = end - start
            state = simulation.context.getState(getPositions=True)
            PDBxFile.writeFile(pdb.topology, state.getPositions(), open(self.save_path+'MultiEM_afterMD.cif', 'w'))
            print(f'Everything is done! Simulation finished succesfully!\nMD finished in {elapsed/60:.2f} minutes.\n')

        if self.args.SAVE_PLOTS:
            print('Creating and saving plots...')
            plot_projection(get_coordinates_mm(state.getPositions()),self.Cs,save_path=self.save_path)
            print('Done! :)')

def main():
    # Input data
    args = get_config()
    
    # Run simulation
    md = MultiEM(args)
    md.run_pipeline()
    if md.args.NUC_SIM: md.run_nuc_pipeline()

if __name__=='__main__':
    main()
