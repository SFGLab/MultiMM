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
from nucs_init_struct_tools import *
from nucs_preprocessing import *
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
        if self.args.CHROM=='':
            for c, i in enumerate(range(len(self.chr_ends)-1)):
                self.chrom_spin[self.chr_ends[i]:self.chr_ends[i+1]] = c

    def add_forcefield(self):
        '''
        Here we define the forcefield of MultiEM.
        '''
        # Leonard-Jones potential for excluded volume
        if self.args.EV_USE_EXCLUDED_VOLUME:
            self.ev_force = mm.CustomNonbondedForce('epsilon*(sigma/r)^3')
            self.ev_force.addGlobalParameter('epsilon', defaultValue=self.args.EV_EPSILON)
            self.ev_force.addGlobalParameter('sigma', defaultValue=np.min(self.ds))
            for i in range(self.system.getNumParticles()):
                self.ev_force.addParticle()
            self.system.addForce(self.ev_force)
        
        # Gaussian compartmentalization potential - compartment blocks
        radius1 = (self.args.N_BEADS/50000)**(1/3) if self.args.SC_RADIUS1==None else self.args.SC_RADIUS1
        radius2 = (self.args.N_BEADS/50000)**(1/3)*6 if self.args.SC_RADIUS2==None else self.args.SC_RADIUS2
        r_comp = (radius2-radius1)/20 if self.args.COB_DISTANCE==None else self.args.COB_DISTANCE
        r_chrom = (radius2-radius1)/20 if self.args.CHB_DISTANCE==None else self.args.CHB_DISTANCE

        ## Compartment force
        if self.args.COB_USE_COMPARTMENT_BLOCKS:
            self.comp_force = mm.CustomNonbondedForce('-E*exp(-r^2/(2*r0^2)); E=(Ea*delta(s1-1)*delta(s2-1)+Eb*delta(s1+1)*delta(s2+1))')
            self.comp_force.addGlobalParameter('r0',defaultValue=r_comp)
            self.comp_force.addGlobalParameter('Ea',defaultValue=self.args.COB_EA)
            self.comp_force.addGlobalParameter('Eb',defaultValue=self.args.COB_EB)
            self.comp_force.addPerParticleParameter('s')
            for i in range(self.system.getNumParticles()):
                self.comp_force.addParticle([self.Cs[i]])
            self.system.addForce(self.comp_force)
        
        ## Subcompartment force
        if self.args.SCB_USE_SUBCOMPARTMENT_BLOCKS:
            self.comp_force = mm.CustomNonbondedForce('-E*exp(-r^2/(2*r0^2)); E=((Ea1*delta(s1-2)*delta(s2-2)+Ea2*delta(s1-1)*delta(s2-1)+Eb1*delta(s1+1)*delta(s2+1)+Eb2*delta(s1+2)*delta(s2+2))')
            self.comp_force.addGlobalParameter('r0',defaultValue=r_comp)
            self.comp_force.addGlobalParameter('Ea1',defaultValue=self.args.SCB_EA1)
            self.comp_force.addGlobalParameter('Ea2',defaultValue=self.args.SCB_EA2)
            self.comp_force.addGlobalParameter('Eb1',defaultValue=self.args.SCB_EB1)
            self.comp_force.addGlobalParameter('Eb2',defaultValue=self.args.SCB_EB2)
            self.comp_force.addPerParticleParameter('s')
            for i in range(self.system.getNumParticles()):
                self.comp_force.addParticle([self.Cs[i]])
            self.system.addForce(self.comp_force)

        # Chromosomal blocks
        if self.args.CHB_USE_CHROMOSOMAL_BLOCKS:
            self.chrom_block_force = mm.CustomNonbondedForce('-E*exp(-r^2/(2*r1^2)); E=dE*delta(chrom1-chrom2)')
            self.chrom_block_force.addGlobalParameter('r1',defaultValue=r_chrom)
            self.chrom_block_force.addGlobalParameter('dE',defaultValue=self.args.CHB_DE)
            self.chrom_block_force.addPerParticleParameter('chrom')
            for i in range(self.system.getNumParticles()):
                self.chrom_block_force.addParticle([self.chrom_spin[i]])
            self.system.addForce(self.chrom_block_force)

        # Spherical container
        if self.args.SC_USE_SPHERICAL_CONTAINER:
            self.container_force = mm.CustomExternalForce('C*(max(0, r-R2)^2+max(0, R1-r)^2); r=sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2)')
            self.container_force.addGlobalParameter('C',defaultValue=self.args.SC_SCALE)
            self.container_force.addGlobalParameter('R1',defaultValue=radius1)
            self.container_force.addGlobalParameter('R2',defaultValue=radius2)
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
            self.Alamina_force.addGlobalParameter('R1',defaultValue=radius1)
            self.Alamina_force.addGlobalParameter('R2',defaultValue=radius2)
            self.Alamina_force.addGlobalParameter('x0',defaultValue=self.mass_center[0])
            self.Alamina_force.addGlobalParameter('y0',defaultValue=self.mass_center[1])
            self.Alamina_force.addGlobalParameter('z0',defaultValue=self.mass_center[2])
            self.Alamina_force.addPerParticleParameter('s')
            for i in range(self.system.getNumParticles()):
                self.Alamina_force.addParticle(i, [self.Cs[i]])
            self.system.addForce(self.Alamina_force)
            
        # Interaction of B compartment with lamina
        if self.args.IBL_USE_B_LAMINA_INTERACTION:
            if radius1!=0.0:
                self.Blamina_force = mm.CustomExternalForce('B*(sin(pi*(r-R1)/(R2-R1))^2-1)*(delta(s+1)+delta(s+2))*step(R1)*(1-step(R2)); r=sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2)')
            else:
                self.Blamina_force = mm.CustomExternalForce('B*(sin(pi*(r-R1)/(R2-R1))^2-1)*(delta(s+1)+delta(s+2))*step(R2/2)*(1-step(R2)); r=sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2)')
            self.Blamina_force.addGlobalParameter('B',defaultValue=self.args.IBL_SCALE)
            self.Blamina_force.addGlobalParameter('pi',defaultValue=np.pi)
            self.Blamina_force.addGlobalParameter('R1',defaultValue=radius1)
            self.Blamina_force.addGlobalParameter('R2',defaultValue=radius2)
            self.Blamina_force.addGlobalParameter('x0',defaultValue=self.mass_center[0])
            self.Blamina_force.addGlobalParameter('y0',defaultValue=self.mass_center[1])
            self.Blamina_force.addGlobalParameter('z0',defaultValue=self.mass_center[2])
            self.Blamina_force.addPerParticleParameter('s')
            for i in range(self.system.getNumParticles()):
                self.Blamina_force.addParticle(i, [self.Cs[i]])
            self.system.addForce(self.Blamina_force)

        if self.args.CF_USE_CENTRAL_FORCE:
            # Force that sets smaller chromosomes closer to the center
            self.central_force = mm.CustomExternalForce('G*(chrom-1)/23*(-1/(r-R1+1)+1/(r-R1+1)^2); r=sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2)')
            self.central_force.addGlobalParameter('G',defaultValue=self.args.CF_STRENGTH)
            self.central_force.addGlobalParameter('R1',defaultValue=R1)
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
                    self.loop_force.addBond(m,n,self.ds,self.args.LE_HARMONIC_BOND_K)
            self.system.addForce(self.loop_force)

        # Bending potential for stiffness
        if self.POL_USE_HARMONIC_ANGLE:
            self.angle_force = mm.HarmonicAngleForce()
            for i in range(self.system.getNumParticles()-2):
                if (i not in self.chr_ends) and (i not in self.chr_ends-1): self.angle_force.addAngle(i, i+1, i+2, np.pi, 20)
            self.system.addForce(self.angle_force)

    def add_nuc_forcefield(self,n_wraps=2):
        '''
        Here we define the forcefield of the nucleosome model.
        '''
        # Bending potential for stiffness
        angle_force = mm.HarmonicAngleForce()
        for i in range(4*self.num_nucs,self.system_n.getNumParticles()-2):
            if self.is_nuc[i-4*self.num_nucs]:
                angle_force.addAngle(i, i+1, i+2, self.args.NAF_HARMONIC_ANGLE_R0, self.args.NAF_HARD_HARMONIC_CONSTANT_K)
            else:
                angle_force.addAngle(i, i+1, i+2, self.args.NAF_HARMONIC_ANGLE_R0, self.args.NAF_SOFT_HARMONIC_CONSTANT_K)
        angle_force.setForceGroup(2)
        self.system_n.addForce(angle_force)
        print('Angle force imported.')

        # Add DNA-DNA interactions
        print('Importing DNA-DNA interactions...')
        dnadna_force = mm.HarmonicBondForce()
        r = int(np.average(self.exit_points-self.entry_points))//n_wraps
        res = []
        for k in tqdm(range(self.num_nucs)):
            for j in range(self.entry_points[k], self.exit_points[k]-r):
                dnadna_force.addBond(4*self.num_nucs+j, 4*self.num_nucs+j+r, self.args.NDD_HARMONIC_LENGTH_R0, self.args.NDD_ARMONIC_CONSTANT_K)
                res.append(f':{self.num_nucs+j+1}\t:{self.num_nucs+j+r+1}\tred\n')
        dnadna_force.setForceGroup(3)
        self.system_n.addForce(dnadna_force)
        print('DNA-DNA force imported.')

        # Interaction between histone and DNA
        print('Import DNA-histone forcefield...')
        histone_dna_force = mm.HarmonicBondForce()
        atoms = ['HIA','HIB','HIC','HID']
        for k in tqdm(range(self.num_nucs)):
            for i in range(4*k,4*k+4):
                for w in range(n_wraps):
                    helix_factor = (self.exit_points[k]-self.entry_points[k])//n_wraps
                    wrap_factor = (self.exit_points[k]-self.entry_points[k])//n_wraps//4
                    for j in range(self.entry_points[k]+w*helix_factor+(i%4)*wrap_factor, self.entry_points[k]+w*helix_factor+(i%4+1)*wrap_factor):
                        histone_dna_force.addBond(i, 4*self.num_nucs+j, self.args.NDH_HARMONIC_LENGTH_R0, self.args.NDH_HARMONIC_CONSTANT_K)
                        res.append(f':{k+1}@{atoms[i%4]}\t:{self.num_nucs+j+1}\tgreen\n')
        histone_dna_force.setForceGroup(4)
        self.system_n.addForce(histone_dna_force)
        print('DNA-histone force imported.')

    def run_pipeline(self):
        '''
        Energy minimization for GW model.
        '''
        # Initialize simulation
        if self.args.BUILD_INITIAL_STRUCTURE:
            print('\nCreating initial structure...')
            comp_mode = 'compartments' if np.all(self.Cs!=None) and len(np.unique(self.Cs))<=3 else 'subcompartments'
            if np.all(self.Cs!=None): write_cmm(self.Cs,name=self.save_path+'MultiEM_compartment_colors.cmd')
            pdb_content = build_init_mmcif(n_dna=self.args.N_BEADS,chrom_ends=self.chr_ends,path=self.save_path,hilbert=self.args.CHROM==None)
            print('---Done!---')
        pdb = PDBxFile(self.save_path+'MultiEM_init.cif') if self.args.INITIAL_STRUCTURE_path==None or build_init_mmcif else PDBxFile(self.args.INITIAL_STRUCTURE_PATH)
        self.mass_center = np.average(get_coordinates_mm(pdb.positions),axis=0)
        forcefield = ForceField('forcefields/ff.xml')
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
        simulation.context.setVelocitiesToself.args.SIM_TEMPERATURE(self.args.SIM_TEMPERATURE, 0)
        current_platform = simulation.context.getPlatform()
        print(f"Simulation will run on platform: {current_platform.getName()}.")
        start_time = time.time()
        simulation.minimizeEnergy()
        state = simulation.context.getState(getPositions=True)
        PDBxFile.writeFile(pdb.topology, state.getPositions(), open(self.save_path+'MultiEM_minimized.cif', 'w'))
        print(f"--- Energy minimization done!! Executed in {(time.time() - start_time)/60:.2f} minutes. :D ---\n")

        if self.args.SAVE_PLOTS: plot_projection(get_coordinates_mm(state.getPositions()),self.Cs,save_path=self.save_path)
        
        # if self.args.CHROM==None:
        start = time.time()
        V = get_coordinates_mm(state.getPositions())
        for i in range(len(self.chr_ends)-1):
            write_mmcif_chrom(coords=10*V[self.chr_ends[i]:self.chr_ends[i+1]],\
                        path=self.save_path+f'chromosomes/MultiEM_minimized_{chrs[i]}.cif')

        if self.args.SIM_RUN_MD:
            simulation.reporters.append(StateDataReporter(stdout, (self.args.SIM_N_STEPS*self.args.SAMPLING_STEP)//20, step=True, totalEnergy=True, potentialEnergy=True, Temperature=True))
            simulation.reporters.append(DCDReporter(self.save_path+'MultiEM_annealing.dcd', (self.args.SIM_N_STEPS*self.args.SAMPLING_STEP)//(self.args.TRJ_FRAMES)))
            print('Running Simulated Annealing...')
            start = time.time()
            for i in range(self.args.SIM_N_STEPS//self.args.SAMPLING_STEP):
                simulation.integrator.setself.args.SIM_TEMPERATURE(self.args.SIM_TEMPERATURE)
                simulation.step(self.args.SAMPLING_STEP)
                if (i*self.args.SAMPLING_STEP)%10==0:
                    self.state = simulation.context.getState(getPositions=True)
                    PDBxFile.writeFile(pdb.topology, self.state.getPositions(), open(self.save_path+f'ensembles/ens_{i//100+1}.cif', 'w'))
            end = time.time()
            elapsed = end - start
            state = simulation.context.getState(getPositions=True)
            PDBxFile.writeFile(pdb.topology, state.getPositions(), open(self.save_path+'MultiEM_afterMD.cif', 'w'))
            print(f'Everything is done! Simulation finished succesfully!\nMD finished in {elapsed/60:.2f} minutes.\n')

    def run_nuc_pipeline(self):
        '''
        Energy minimization for nucleosome simulation.
        '''
        # Define system
        print('\nNucleosome simulation length:',self.nuc_sim_len)
        print('Number of nucleosomes',self.num_nucs)
        print('\n\nBuilding initial structure of nucleosome simulation...')
        pdb_content = build_init_mmcif_nucs(entry=self.entry_points,exit=self.exit_points,n_nucs=self.num_nucs,n_dna=self.nuc_sim_len,
                               mode='path',psf=True,path=self.save_path+'MultiEM_minimized.cif')

        print('Creating system...')
        pdb = PDBxFile('dna_histones.cif')
        forcefield = ForceField('forcefields/dna_histones_ff.xml')
        self.system_n = forcefield.createSystem(pdb.topology, nonbondedCutoff=1*mm.unit.nanometer)
        integrator = mm.LangevinIntegrator(310, 0.5, 100 * mm.unit.femtosecond)
        print('Done\n')

        # Add nucleosome forcefield
        print('Adding forcefield for nucleosomes...')
        self.add_nuc_forcefield()
        print('Done\n')
        
        # Energy minimization
        print('Energy minimizing of nucleosome potential (this may take several minutes or hours)...')
        platform = mm.Platform.getPlatformByName('CUDA')
        simulation = Simulation(pdb.topology, self.system_n, integrator, platform)
        current_platform = simulation.context.getPlatform()
        print(f"Simulation will run on platform: {current_platform.getName()}")
        simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, totalEnergy=True, potentialEnergy=True, Temperature=True))
        simulation.context.setPositions(pdb.positions)
        simulation.minimizeEnergy()
        state = simulation.context.getState(getPositions=True)
        PDBxFile.writeFile(pdb.topology, state.getPositions(), open(self.save_path+f'minimized_nucres_model.cif', 'w'))
        print('Energy minimization of nucleosome model done :D\n')

def main():
    # Input data
    args = get_config()
    
    # Run simulation
    md = MultiEM(args)
    md.run_pipeline()
    if md.args.NUC_SIM: md.run_nuc_pipeline()

if __name__=='__main__':
    main()
