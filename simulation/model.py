import numpy as np
import os
from typing import List
from sys import stdout
import openmm as mm
from openmm.app import PDBxFile, ForceField, Simulation, DCDReporter, StateDataReporter
from .args_definition import *
from .initial_structure_tools import *
from .utils import *
from .plots import *
from .nucleosome_interpolation import *
import time
import os

class MultiMM:
    def __init__(self,args):
        '''
        Input data:
        ------------
        args: list of arguments imported from config.ini file.
        '''
        # Output folder
        self.ms, self.ns, self.ds, self.chr_ends, self.Cs = None, None, None, None, None
        
        # Make save directory
        self.save_path = args.OUT_PATH+'/'
        # Create main save directory and subdirectories if they don't exist
        os.makedirs(os.path.join(self.save_path, 'ensembles'), exist_ok=True)
        os.makedirs(os.path.join(self.save_path, 'chromosomes'), exist_ok=True)
        os.makedirs(os.path.join(self.save_path, 'plots'), exist_ok=True)

        self.args  = args
        coords = [args.LOC_START,args.LOC_END] if args.LOC_START!=None else None
        
        # Compartments
        if args.COMPARTMENT_PATH!=None:
            if args.COMPARTMENT_PATH.lower().endswith('.bed'):
                self.Cs, self.chr_ends, self.chrom_idxs = import_bed(bed_file=args.COMPARTMENT_PATH,N_beads=self.args.N_BEADS,\
                                                                        chrom=self.args.CHROM,coords=coords,\
                                                                        save_path=self.save_path,\
                                                                        shuffle=args.SHUFFLE_CHROMS,seed=args.SHUFFLING_SEED)
            else:
                raise InterruptedError('Compartments file should be in .bed format.')
        
        # Loops
        if args.LOOPS_PATH.lower().endswith('.bedpe'):
            self.ms, self.ns, self.ds, self.chr_ends, self.chrom_idxs = import_mns_from_bedpe(bedpe_file = args.LOOPS_PATH,N_beads=self.args.N_BEADS,\
                                                                            coords = coords, chrom=args.CHROM,\
                                                                            path=self.save_path,\
                                                                            shuffle=args.SHUFFLE_CHROMS,seed=args.SHUFFLING_SEED,down_prob=args.DOWNSAMPLING_PROB)
        else:
            raise InterruptedError('You did not provide appropriate loop file. Loop .bedpe file is obligatory.')

        # Nucleosomes
        if args.NUC_DO_INTERPOLATION and args.ATACSEQ_PATH!=None:
            if args.ATACSEQ_PATH.lower().endswith('.bw') or args.ATACSEQ_PATH.lower().endswith('.bigwig'):
                self.atacseq = import_bw(args.ATACSEQ_PATH,self.args.N_BEADS,chrom=self.args.CHROM,coords=coords,\
                                         shuffle=args.SHUFFLE_CHROMS,seed=args.SHUFFLING_SEED)
            else:
                raise InterruptedError('ATAC-Seq file should be in .bw or .BigWig format.')
        
        if self.args.CHROM=='': write_chrom_colors(self.chr_ends,self.chrom_idxs,name=self.save_path+'MultiMM_chromosome_colors.cmd')

        # Chromosomes
        self.chrom_spin, self.chrom_strength = np.zeros(self.args.N_BEADS), np.zeros(self.args.N_BEADS)
        if self.args.CHROM==None or self.args.CHROM=='':
            for i in range(len(self.chr_ends)-1):
                self.chrom_spin[self.chr_ends[i]:self.chr_ends[i+1]] = self.chrom_idxs[i]
                self.chrom_strength[self.chr_ends[i]:self.chr_ends[i+1]] = chrom_strength[i]

    def add_evforce(self):
        sigma = self.args.LE_HARMONIC_BOND_R0
        self.ev_force = mm.CustomNonbondedForce(f'epsilon*(sigma/(r+r_small))^{self.args.EV_POWER}')
        self.ev_force.setForceGroup(1)
        self.ev_force.addGlobalParameter('epsilon', defaultValue=self.args.EV_EPSILON)
        self.ev_force.addGlobalParameter('r_small', defaultValue=self.args.EV_R_SMALL)
        self.ev_force.addGlobalParameter('sigma', defaultValue=sigma)
        for i in range(self.system.getNumParticles()):
            self.ev_force.addParticle()
        self.system.addForce(self.ev_force)
    
    def add_compartment_blocks(self):
        self.comp_force = mm.CustomNonbondedForce('-E*exp(-r^2/(2*rc^2)); E=(Ea*(delta(s1-1)+delta(s1-2))*(delta(s2-1)+delta(s2-2))+Eb*(delta(s1+1)+delta(s1+2))*(delta(s2+1)+delta(s2+2))')
        self.comp_force.setForceGroup(1)
        self.comp_force.addGlobalParameter('rc',defaultValue=self.r_comp)
        self.comp_force.addGlobalParameter('Ea',defaultValue=self.args.COB_EA)
        self.comp_force.addGlobalParameter('Eb',defaultValue=self.args.COB_EB)
        self.comp_force.addPerParticleParameter('s')
        for i in range(self.system.getNumParticles()):
            self.comp_force.addParticle([self.Cs[i]])
        self.system.addForce(self.comp_force)

    def add_subcompartment_blocks(self):
        self.scomp_force = mm.CustomNonbondedForce('-E*exp(-r^2/(2*rsc^2)); E=Ea1*delta(s1-2)*delta(s2-2)+Ea2*delta(s1-1)*delta(s2-1)+Eb1*delta(s1+1)*delta(s2+1)+Eb2*delta(s1+2)*delta(s2+2)')
        self.scomp_force.setForceGroup(1)
        self.scomp_force.addGlobalParameter('rsc',defaultValue=self.r_comp)
        self.scomp_force.addGlobalParameter('Ea1',defaultValue=self.args.SCB_EA1)
        self.scomp_force.addGlobalParameter('Ea2',defaultValue=self.args.SCB_EA2)
        self.scomp_force.addGlobalParameter('Eb1',defaultValue=self.args.SCB_EB1)
        self.scomp_force.addGlobalParameter('Eb2',defaultValue=self.args.SCB_EB2)
        self.scomp_force.addPerParticleParameter('s')
        for i in range(self.system.getNumParticles()):
            self.scomp_force.addParticle([self.Cs[i]])
        self.system.addForce(self.scomp_force)
    
    def add_chromosomal_blocks(self):
        self.chrom_block_force = mm.CustomNonbondedForce('E*(k_C*r^4-r^3+r^2); E=dE*delta(chrom1-chrom2)')
        self.chrom_block_force.setForceGroup(2)
        self.chrom_block_force.addGlobalParameter('k_C',defaultValue=self.args.CHB_KC)
        self.chrom_block_force.addGlobalParameter('dE',defaultValue=self.args.CHB_DE)
        self.chrom_block_force.addPerParticleParameter('chrom')
        for i in range(self.system.getNumParticles()):
            self.chrom_block_force.addParticle([self.chrom_spin[i]])
        self.system.addForce(self.chrom_block_force)

    def add_spherical_container(self):
        self.container_force = mm.CustomExternalForce('C*(max(0, r-R2)^2+max(0, R1-r)^2); r=sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2)')
        self.container_force.setForceGroup(2)
        self.container_force.addGlobalParameter('C',defaultValue=self.args.SC_SCALE)
        self.container_force.addGlobalParameter('R1',defaultValue=self.radius1)
        self.container_force.addGlobalParameter('R2',defaultValue=self.radius2)
        self.container_force.addGlobalParameter('x0',defaultValue=self.mass_center[0])
        self.container_force.addGlobalParameter('y0',defaultValue=self.mass_center[1])
        self.container_force.addGlobalParameter('z0',defaultValue=self.mass_center[2])
        for i in range(self.system.getNumParticles()):
            self.container_force.addParticle(i, [])
        self.system.addForce(self.container_force)

    def add_Blamina_interaction(self):
        if self.radius1!=0.0:
            self.Blamina_force = mm.CustomExternalForce('B*(sin(pi*(r-R1)/(R2-R1))^8-1)*(delta(s+1)+delta(s+2)); r=sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2)')
        else:
            self.Blamina_force = mm.CustomExternalForce('B*(sin(pi*(r-R1)/(R2-R1))^8-1)*(delta(s+1)+delta(s+2)); r=sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2)')
        self.Blamina_force.setForceGroup(2)
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
    
    def add_central_force(self):
        self.central_force = mm.CustomExternalForce('G*chrom_s*(sin(r-3*R1/2)+(r-3*R1/2)^2); r=sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2)')
        self.central_force.setForceGroup(2)
        self.central_force.addGlobalParameter('G',defaultValue=self.args.CF_STRENGTH)
        self.central_force.addGlobalParameter('R1',defaultValue=self.radius1)
        self.central_force.addGlobalParameter('x0',defaultValue=self.mass_center[0])
        self.central_force.addGlobalParameter('y0',defaultValue=self.mass_center[1])
        self.central_force.addGlobalParameter('z0',defaultValue=self.mass_center[2])
        self.central_force.addPerParticleParameter('chrom_s')
        for i in range(self.system.getNumParticles()):
            self.central_force.addParticle(i, [self.chrom_strength[i]])
        self.system.addForce(self.central_force)

    def add_harmonic_bonds(self):
        self.bond_force = mm.HarmonicBondForce()
        self.bond_force.setForceGroup(1)
        for i in range(self.system.getNumParticles()-1):
            if (i not in self.chr_ends): self.bond_force.addBond(i,i+1,self.args.POL_HARMONIC_BOND_R0,self.args.POL_HARMONIC_BOND_K)
        self.system.addForce(self.bond_force)
    
    def add_loops(self):
        self.loop_force = mm.HarmonicBondForce()
        self.loop_force.setForceGroup(1)
        counter=0
        for m,n in tqdm(zip(self.ms,self.ns),total=len(self.ms)):
            if self.args.LE_FIXED_DISTANCES:
                self.loop_force.addBond(m,n,self.args.LE_HARMONIC_BOND_R0,self.args.LE_HARMONIC_BOND_K)
            else:
                self.loop_force.addBond(m,n,self.ds[counter],self.args.LE_HARMONIC_BOND_K)
            counter+=1
        self.system.addForce(self.loop_force)

    def add_stiffness(self):
        self.angle_force = mm.HarmonicAngleForce()
        self.angle_force.setForceGroup(1)
        for i in range(self.system.getNumParticles()-2):
            if (i not in self.chr_ends) and (i not in self.chr_ends-1): 
                self.angle_force.addAngle(i, i+1, i+2, self.args.POL_HARMONIC_ANGLE_R0, self.args.POL_HARMONIC_ANGLE_CONSTANT_K)
        self.system.addForce(self.angle_force)

    def initialize_simulation(self):
        if self.args.BUILD_INITIAL_STRUCTURE:
            print('\nCreating initial structure...')
            comp_mode = 'compartments' if np.all(self.Cs!=None) and len(np.unique(self.Cs))<=3 else 'subcompartments'
            if np.all(self.Cs!=None): write_cmm(self.Cs,name=self.save_path+'MultiMM_compartment_colors.cmd')
            pdb_content = build_init_mmcif(n_dna=self.args.N_BEADS,chrom_ends=self.chr_ends,\
                                           path=self.save_path,curve=self.args.INITIAL_STRUCTURE_TYPE,scale=(self.radius1+self.radius2)/2)
            print('---Done!---')
        self.pdb = PDBxFile(self.save_path+'MultiMM_init.cif') if self.args.INITIAL_STRUCTURE_path==None or build_init_mmcif else PDBxFile(self.args.INITIAL_STRUCTURE_PATH)
        self.mass_center = np.average(get_coordinates_mm(self.pdb.positions),axis=0)
        forcefield = ForceField(self.args.FORCEFIELD_PATH)
        self.system = forcefield.createSystem(self.pdb.topology)

        match args.SIM_INTEGRATOR_TYPE:
            case 'verlet':
                self.integrator = mm.VerletIntegrator(self.args.SIM_INTEGRATOR_STEP)
            case 'variable_verlet':
                self.integrator = mm.VariableVerletIntegrator(self.SIM_ERROR_TOLERANCE)
            case 'langevin':
                self.integrator = mm.LangevinIntegrator(self.args.SIM_TEMPERATURE, self.args.SIM_FRICTION_COEFF, self.args.SIM_INTEGRATOR_STEP)
            case 'variable_langevin':
                self.integrator = mm.VariableLangevinIntegrator(self.args.SIM_TEMPERATURE, self.args.SIM_FRICTION_COEFF, self.SIM_ERROR_TOLERANCE)
            case 'amd':
                self.integrator = mm.amd.AMDIntegrator(self.args.SIM_INTEGRATOR_STEP, self.args.SIM_AMD_ALPHA, self.args.SIM_AMD_E)
            case 'brownian':
                self.integrator = mm.BrownianIntegrator(self.args.SIM_TEMPERATURE, self.args.SIM_FRICTION_COEFF, self.args.SIM_INTEGRATOR_STEP)
            
    def add_forcefield(self):
        '''
        Here we define the forcefield of MultiMM.
        '''
        # Add forces
        print('\nImporting forcefield...')
        if self.args.EV_USE_EXCLUDED_VOLUME: self.add_evforce()
        if self.args.COB_USE_COMPARTMENT_BLOCKS: self.add_compartment_blocks()
        if self.args.SCB_USE_SUBCOMPARTMENT_BLOCKS: self.add_subcompartment_blocks()
        if self.args.CHB_USE_CHROMOSOMAL_BLOCKS: self.add_chromosomal_blocks()
        if self.args.SC_USE_SPHERICAL_CONTAINER: self.add_spherical_container()
        if self.args.IBL_USE_B_LAMINA_INTERACTION: self.add_Blamina_interaction()
        if self.args.CF_USE_CENTRAL_FORCE: self.add_central_force()
        if self.args.POL_USE_HARMONIC_BOND: self.add_harmonic_bonds()
        if self.args.LE_USE_HARMONIC_BOND: self.add_loops()
        if self.args.POL_USE_HARMONIC_ANGLE: self.add_stiffness()
    
    def min_energy(self):
        print('\nEnergy minimization...')
        # Try to use CUDA or OpenCL, fall back to CPU if not available
        try:
            platform = mm.Platform.getPlatformByName(self.args.PLATFORM)
            if platform.getName() not in ["CUDA", "OpenCL"]:
                raise Exception(f"{self.args.PLATFORM} is not CUDA or OpenCL")
        except Exception as e:
            print(f"Failed to find CUDA or OpenCL: {e}. Falling back to CPU.")
            platform = mm.Platform.getPlatformByName('CPU')
        
        # Run the simulation
        self.simulation = Simulation(self.pdb.topology, self.system, self.integrator, platform)     
        self.simulation.context.setPositions(self.pdb.positions)
        self.simulation.context.setVelocitiesToTemperature(self.args.SIM_TEMPERATURE, 0)

        # Report which platform is being used
        current_platform = self.simulation.context.getPlatform()
        print(f"Simulation will run on platform: {current_platform.getName()}.")

        # Perform energy minimization
        start_time = time.time()
        self.simulation.minimizeEnergy()

        # Save the minimized structure
        self.state = self.simulation.context.getState(getPositions=True)
        PDBxFile.writeFile(self.pdb.topology, self.state.getPositions(), open(self.save_path+'MultiMM_minimized.cif', 'w'))
        print(f"--- Energy minimization done!! Executed in {(time.time() - start_time)//3600:.0f} hours, {(time.time() - start_time)%3600//60:.0f} minutes and  {(time.time() - start_time)%60:.0f} seconds. :D ---\n")
        print('---Done!---')

    def save_chromosomes(self):
        V = get_coordinates_mm(self.state.getPositions())
        for i in range(len(self.chr_ends)-1):
            write_mmcif_chrom(coords=10*V[self.chr_ends[i]:self.chr_ends[i+1]],\
                        path=self.save_path+f'chromosomes/MultiMM_minimized_{chrs[self.chrom_idxs[i]]}.cif')

    def run_md(self):
        self.simulation.reporters.append(StateDataReporter(stdout, self.args.SIM_SAMPLING_STEP, step=True, totalEnergy=True, kineticEnergy=True ,potentialEnergy=True, temperature=True, separator='\t'))
        self.simulation.reporters.append(DCDReporter(self.save_path+'MultiMM_annealing.dcd', self.args.SIM_N_STEPS//self.args.TRJ_FRAMES))
        print('Running relaxation...')
        start = time.time()
        for i in range(self.args.SIM_N_STEPS//self.args.SIM_SAMPLING_STEP):
            self.simulation.step(self.args.SIM_SAMPLING_STEP)
            self.state = self.simulation.context.getState(getPositions=True)
            PDBxFile.writeFile(self.pdb.topology, self.state.getPositions(), open(self.save_path+f'ensembles/ens_{i+1}.cif', 'w'))
        end = time.time()
        elapsed = end - start
        self.state = self.simulation.context.getState(getPositions=True)
        PDBxFile.writeFile(self.pdb.topology, self.state.getPositions(), open(self.save_path+'MultiMM_afterMD.cif', 'w'))
        print(f'Everything is done! Simulation finished succesfully!\nMD finished in {elapsed//3600:.0f} hours, {elapsed%3600//60:.0f} minutes and  {elapsed%60:.0f} seconds. ---\n')

    def nuc_interpolation(self):
        print('Running nucleosome interpolation...')
        start = time.time()
        nuc_interpol = NucleosomeInterpolation(get_coordinates_cif(self.save_path+'MultiMM_minimized.cif'),self.atacseq,\
                    self.args.MAX_NUCS_PER_BEAD, self.args.NUC_RADIUS, self.args.POINTS_PER_NUC, self.args.PHI_NORM)
        Vnuc = nuc_interpol.interpolate_structure_with_nucleosomes()
        write_mmcif_chrom(Vnuc,path=self.save_path+f'MultiMM_minimized_with_nucs.cif')
        end = time.time()
        elapsed = end - start
        print(f'Nucleosome interpolation finished succesfully in {elapsed//3600:.0f} hours, {elapsed%3600//60:.0f} minutes and  {elapsed%60:.0f} seconds.')
    
    def set_radiuses(self):
        self.radius1 = (self.args.N_BEADS/50000)**(1/3) if self.args.SC_RADIUS1==None else self.args.SC_RADIUS1
        self.radius2 = 3.5*(self.args.N_BEADS/50000)**(1/3) if self.args.SC_RADIUS2==None else self.args.SC_RADIUS2
        if self.args.COB_DISTANCE!=None:
            self.r_comp = self.args.COB_DISTANCE
        elif self.args.SCB_DISTANCE!=None:
            self.r_comp = self.args.SCB_DISTANCE
        else:
            self.r_comp = (self.radius2-self.radius1)/20
    
    def run(self):
        '''
        Energy minimization for GW model.
        '''
        # Estimation of parameters
        self.set_radiuses()
        
        # Initialize simulation
        self.initialize_simulation()
        
        # Import forcefield
        self.add_forcefield()
        
        # Run simulation / Energy minimization
        self.min_energy()
        self.save_chromosomes() 
        
        # Run molecular dynamics
        if self.args.SIM_RUN_MD: self.run_md()
        
        # Make diagnostic plots
        if self.args.SAVE_PLOTS and np.any(self.Cs!=None):
            print('Creating and saving plots...')
            plot_projection(get_coordinates_mm(self.state.getPositions()),self.Cs,save_path=self.save_path)
            print('Done! :)\n')
        
        # Run nucleosome interpolation
        if self.args.NUC_DO_INTERPOLATION and args.ATACSEQ_PATH!=None:
            self.nuc_interpolation()