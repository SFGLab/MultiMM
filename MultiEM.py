import random as rd
import matplotlib.pyplot as plt
import numpy as np
import os
from mdtraj.reporters import HDF5Reporter
from scipy import ndimage
from typing import List
from sys import stdout
import openmm as mm
from openmm.app import PDBFile, PDBxFile, ForceField, Simulation, PDBReporter, PDBxReporter, DCDReporter, StateDataReporter, CharmmPsfFile
from MultiEM_init_tools import *
from MultiEM_utils import *
from MultiEM_metrics import save_metrics
import time
import os

class MultiEM:
    def __init__(self,N_beads=None,loop_path=None,comp_path=None,n_chrom=24,comp_gw=True,out_path='results',loops_mode='k'):
        '''
        Cs (np arrays): Cs array with colors over time.
        N_beads (int): The number of beads of initial structure.
        step (int): sampling rate.
        '''
        # Create output folder only if needed
        self.save_path = out_path+'/'
        self.ms, self.ns, self.ds, self.ks, self.cs, self.chr_ends, self.Cs = None,None,None,None,None,None,None
        try:
            os.mkdir(self.save_path)
            os.mkdir(self.save_path+'chromosomes')
            os.mkdir(self.save_path+'chromosomes_info')
        except OSError as error:
            print("Folder 'chromosomes' already exists!")
        
        # Load from files
        if np.all(comp_path!=None):
            if comp_path.endswith('.bw') or comp_path.endswith('.bigwig') or comp_path.endswith('.BigWig'):
                self.Cs, self.chr_ends = import_bw(bw_path=comp_path,N_beads=N_beads,\
                                                   viz=False,binary=True,\
                                                   n_chroms=n_chrom,path=self.save_path)
            elif comp_path.endswith('.bed'):
                self.Cs, self.chr_ends = import_compartments_from_bed(bed_file=comp_path,N_beads=N_beads,\
                                                                      n_chroms=n_chrom,path=self.save_path)
            elif loop_path=='random':
                self.Cs = shuffle_blocks(self.Cs)

        elif os.path.isfile(self.save_path+'genomewide_signal.npy'):
            self.Cs = np.load(self.save_path+'genomewide_signal.npy')
            if np.all(self.Cs!=None): self.N_beads = len(self.Cs)
        
        if np.all(loop_path!=None):
            if loop_path.endswith('.bedpe'):
                self.ms, self.ns, self.ds, self.ks, self.cs, self.chr_ends = import_mns_from_bedpe(bedpe_file=loop_path,N_beads=N_beads,\
                                                                                n_chroms=n_chrom,viz=False,\
                                                                                path=self.save_path,mode=loops_mode)
            elif(loop_path.endswith('.txt')):
                self.ms, self.ns, self.ds, self.ks, self.cs, self.chr_ends = import_mns_from_txt(txt_file=loop_path,N_beads=N_beads,\
                                                                                                 n_chroms=n_chrom,path=self.save_path,mode=loops_mode)
            elif loop_path=='random':
                self.ms, self.ns, self.ks = generate_arrays(N_loops=N_beads//8, N=N_beads)
            else:
                raise InterruptedError('You did not provide appropriate loop file.')

            self.N_beads = N_beads
        else:
            raise InterruptedError('You did not provide data for loops. Check if the provided file is correct, or if your outpout path is already containing some data.')
        write_chrom_colors(self.chr_ends,name=self.save_path+'MultiEM_chromosome_colors.cmd')
        if np.all(comp_path!=None): self.Cs = align_comps(self.Cs,self.ms,self.chr_ends)

        # Define a chromsome metric
        self.chroms = np.zeros(self.N_beads)
        chr_count = 1
        for i in range(len(self.chr_ends)-1):
            if not comp_gw: self.chroms[self.chr_ends[i]:self.chr_ends[i+1]] = chr_count
            chr_count += 1

    def add_forcefield(self):
        # Leonard-Jones potential for excluded volume
        self.ev_force = mm.CustomNonbondedForce('epsilon*((sigma1+sigma2)/r)^6')
        if not self.run_MD:
            self.ev_force.addGlobalParameter('epsilon', defaultValue=100)
        else:
            self.ev_force.addGlobalParameter('epsilon', defaultValue=1)
        self.ev_force.addPerParticleParameter('sigma')
        for i in range(self.system.getNumParticles()):
            self.ev_force.addParticle([0.1])
        self.system.addForce(self.ev_force)
        
        # Gaussian compartmentalization potential
        if np.all(self.Cs!=None) and len(np.unique(self.Cs)==2):
            self.comp_force = mm.CustomNonbondedForce('E0+E*exp(-(r-r0)^2/(2*sigma^2)); E=(Ea*delta(s1-1)*delta(s2-1)+Eb*delta(s1+1)*delta(s2+1))*delta(chrom1-chrom2)')
            self.comp_force.addGlobalParameter('sigma',defaultValue=0.5)
            self.comp_force.addGlobalParameter('r0',defaultValue=0.2)
            self.comp_force.addGlobalParameter('E0',defaultValue=0.0)
            self.comp_force.addGlobalParameter('Ea',defaultValue=-1.0)
            self.comp_force.addGlobalParameter('Eb',defaultValue=-2.0)
            self.comp_force.addPerParticleParameter('s')
            self.comp_force.addPerParticleParameter('chrom')
            for i in range(self.system.getNumParticles()):
                self.comp_force.addParticle([self.Cs[i],self.chroms[i]])
            self.system.addForce(self.comp_force)
        elif np.all(self.Cs!=None) and len(np.unique(self.Cs)>=4):
            self.comp_force = mm.CustomNonbondedForce('E0+E*exp(-(r-r0)^2/(2*sigma^2)); E=(Ea1*delta(s1-2)*delta(s2-2)+Ea2*delta(s1-1)*delta(s2-1)+Eb1*delta(s1+1)*delta(s2+1)+Eb2*delta(s1+2)*delta(s2+2)')
            self.comp_force.addGlobalParameter('sigma',defaultValue=0.5)
            self.comp_force.addGlobalParameter('r0',defaultValue=0.2)
            self.comp_force.addGlobalParameter('E0',defaultValue=0.0)
            self.comp_force.addGlobalParameter('Ea1',defaultValue=-0.5)
            self.comp_force.addGlobalParameter('Ea2',defaultValue=-1.0)
            self.comp_force.addGlobalParameter('Eb1',defaultValue=-1.5)
            self.comp_force.addGlobalParameter('Eb2',defaultValue=-2.0)
            self.comp_force.addPerParticleParameter('s')
            self.comp_force.addPerParticleParameter('chrom')
            for i in range(self.system.getNumParticles()):
                self.comp_force.addParticle([self.Cs[i],self.chroms[i]])
            self.system.addForce(self.comp_force)
        
        # Spherical container
        radius = (self.N_beads/50000)**(1/3)*6
        self.container_force = mm.CustomExternalForce('C*max(0, r-R)^2; r=sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2)')
        self.container_force.addGlobalParameter('C',defaultValue=1000)
        self.container_force.addGlobalParameter('R',defaultValue=radius)
        self.container_force.addGlobalParameter('x0',defaultValue=self.mass_center[0])
        self.container_force.addGlobalParameter('y0',defaultValue=self.mass_center[1])
        self.container_force.addGlobalParameter('z0',defaultValue=self.mass_center[2])
        for i in range(self.system.getNumParticles()):
            self.container_force.addParticle(i, [])
        self.system.addForce(self.container_force)

        if np.all(self.Cs!=None):
            # Interaction of A compartment with lamina
            self.Alamina_force = mm.CustomExternalForce('-A*abs(sin(pi*r/R))^(1/4)*(delta(s-1)+delta(s-2)); r=sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2)')
            self.Alamina_force.addGlobalParameter('A',defaultValue=500)
            self.Alamina_force.addGlobalParameter('pi',defaultValue=3.14159265358979323846)
            self.Alamina_force.addGlobalParameter('R',defaultValue=radius)
            self.Alamina_force.addGlobalParameter('x0',defaultValue=self.mass_center[0])
            self.Alamina_force.addGlobalParameter('y0',defaultValue=self.mass_center[1])
            self.Alamina_force.addGlobalParameter('z0',defaultValue=self.mass_center[2])
            self.Alamina_force.addPerParticleParameter('s')
            for i in range(self.system.getNumParticles()):
                self.Alamina_force.addParticle(i, [self.Cs[i]])
            self.system.addForce(self.Alamina_force)

            # Interaction of B compartment with lamina
            self.Blamina_force = mm.CustomExternalForce('B*(sin(pi*r/R)^8-1)*(delta(s+1)+delta(s+2)); r=sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2)')
            self.Blamina_force.addGlobalParameter('B',defaultValue=500)
            self.Blamina_force.addGlobalParameter('pi',defaultValue=3.14159265358979323846)
            self.Blamina_force.addGlobalParameter('R',defaultValue=radius)
            self.Blamina_force.addGlobalParameter('x0',defaultValue=self.mass_center[0])
            self.Blamina_force.addGlobalParameter('y0',defaultValue=self.mass_center[1])
            self.Blamina_force.addGlobalParameter('z0',defaultValue=self.mass_center[2])
            self.Blamina_force.addPerParticleParameter('s')
            for i in range(self.system.getNumParticles()):
                self.Blamina_force.addParticle(i, [self.Cs[i]])
            self.system.addForce(self.Blamina_force)

        # Bond force
        self.bond_force = mm.HarmonicBondForce()
        for i in range(self.system.getNumParticles()-1):
            self.bond_force.addBond(i,i+1,0.1,300000)
        self.system.addForce(self.bond_force)
        
        # Add fixed loops force
        if np.all(self.ms!=None) and np.all(self.ns!=None):
            self.loop_force = mm.HarmonicBondForce()
            counter=0
            for m,n in tqdm(zip(self.ms,self.ns),total=len(self.ms)):
                if np.any(self.ks==None) and np.any(self.ds==None):
                    self.loop_force.addBond(m,n,0.1,3000)
                elif np.any(self.ks==None) and np.any(self.ds!=None):
                    self.loop_force.addBond(m,n,self.ds[counter],3000)
                elif np.any(self.ks!=None) and np.any(self.ds==None):
                    self.loop_force.addBond(m,n,0.1,self.ks[counter])
                else:
                    self.loop_force.addBond(m,n,self.ds[counter],self.ks[counter])
            self.system.addForce(self.loop_force)

        # Bending potential for stiffness
        self.angle_force = mm.HarmonicAngleForce()
        for i in range(self.system.getNumParticles()-2):
            self.angle_force.addAngle(i, i+1, i+2, np.pi, 40)
        self.system.addForce(self.angle_force)

    def run_pipeline(self,MD_steps=10000,run_MD=True,write_files=False,plots=False,build_init_struct=True,Temperature=300*mm.unit.kelvin,init_struct_path=None,pltf='CUDA'):
        # Initialize simulation
        self.run_MD=run_MD
        if build_init_struct:
            print('\nCreating initial structure...')
            comp_mode = 'compartments' if np.all(self.Cs!=None) and len(np.unique(self.Cs))<=3 else 'subcompartments'
            if np.all(self.Cs!=None): write_cmm(self.Cs,name=self.save_path+'MultiEM_compartment_colors.cmd',mode=comp_mode)
            pdb_content = build_init_mmcif(n_dna=self.N_beads,path=self.save_path)
            print('---Done!---')
        pdb = PDBxFile(self.save_path+'MultiEM_init.cif') if init_struct_path==None or build_init_mmcif else PDBxFile(init_struct_path)
        self.mass_center = np.average(get_coordinates_mm(pdb.positions),axis=0)
        forcefield = ForceField('ff.xml')
        self.system = forcefield.createSystem(pdb.topology)
        integrator = mm.LangevinIntegrator(Temperature, 2, 100 * mm.unit.femtosecond)
        
        # Import forcefield
        print('\nImporting forcefield...')
        self.add_forcefield()
        print('---Done!---')

        # Run simulation / Energy minimization
        print('\nEnergy minimization...')
        platform = mm.Platform.getPlatformByName(pltf)
        simulation = Simulation(pdb.topology, self.system, integrator, platform)
        simulation.context.setPositions(pdb.positions)
        current_platform = simulation.context.getPlatform()
        print(f"Simulation will run on platform: {current_platform.getName()}")
        start_time = time.time()
        simulation.minimizeEnergy()
        state = simulation.context.getState(getPositions=True)
        PDBxFile.writeFile(pdb.topology, state.getPositions(), open(self.save_path+'MultiEM_minimized.cif', 'w'))
        print(f"--- Energy minimization done!! Executed in {(time.time() - start_time)/60:.2f} minutes. :D ---")

        # Run molecular Dynamics
        if self.run_MD:
            print('\nRunning MD simulation...')
            start = time.time()
            simulation.reporters.append(DCDReporter(self.save_path+'/MultiEM_traj.dcd', 5))
            simulation.reporters.append(StateDataReporter(stdout, MD_steps//10, step=True, totalEnergy=True, potentialEnergy=True, temperature=True))
            simulation.context.setVelocitiesToTemperature(Temperature, 0)
            for i in range(0,MD_steps,10):
                self.ev_force.setGlobalParameterDefaultValue(0,1+99*i/MD_steps)
                simulation.step(10)
            end = time.time()
            elapsed = end - start
            speed = MD_steps / elapsed
            print(f"\n---MD finished in {elapsed/60:.2f} minutes ({speed:0.1f} steps/s)---")
            
            state = simulation.context.getState(getPositions=True)
            PDBxFile.writeFile(pdb.topology, state.getPositions(), open(self.save_path+'MultiEM_afterMD.cif', 'w'))

        start = time.time()
        V = get_coordinates_mm(state.getPositions())
        save_metrics(V,path_name=self.save_path+'GW_')
        for i in range(len(self.chr_ends)-1):
            write_mmcif(coords=10*V[self.chr_ends[i]:self.chr_ends[i+1]],
                        path=self.save_path+f'chromosomes/MultiEM_minimized_{chrs[i]}.cif')
            save_metrics(V[self.chr_ends[i]:self.chr_ends[i+1]],path_name=self.save_path+f'chromosomes_info/{chrs[i]}_')
        
        if plots:
            print('\nComputing heatmap...')
            heat = get_heatmap(mm_vec=state.getPositions(),viz=plots,path=self.save_path)
            end = time.time()
            elapsed = end - start
            print(f'---Heatmap computed in {elapsed/60:0.2f} minutes.---')

def main():
    # Input data
    bw_path = '/home/skorsak/Documents/data/Trios/calder_HiChIP_subcomp/CHS_m.bed'
    # loop_path = '/home/skorsak/Documents/data/Trios/ChiA-PiPE_Loops/loops_pet3+/HG00512_CHS_F_CTCF_1mb_pet3.bedpe'
    loop_path = None
    out_path_name = 'random'
    
    # Run simulation
    md = MultiEM(N_beads=50000,out_path=out_path_name,n_chrom=23,loop_path='random',comp_path=bw_path)
    md.run_pipeline(run_MD=False,build_init_struct=True,
                    init_struct_path=None,plots=False,pltf='OpenCL')

if __name__=='__main__':
    main()