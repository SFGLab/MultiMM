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
    def __init__(self,N_beads=None,loop_path=None,comp_path=None,n_chrom=24,comp_gw=False,out_path='results'):
        '''
        Cs (np arrays): Cs array with colors over time.
        N_beads (int): The number of beads of initial structure.
        step (int): sampling rate.
        '''
        # Create output folder only if needed
        self.save_path = out_path+'/'
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
                self.Cs, self.chr_ends = import_compartments_from_bed(bed_file=comp_path,N_beads=N_beads,
                                                                      n_chroms=n_chrom,path=self.save_path)
        if np.all(loop_path!=None):
            if loop_path.endswith('.bedpe'):
                self.ms, self.ns, self.ds, self.ks, self.cs, self.chr_ends = import_mns_from_bedpe(bedpe_file=loop_path,N_beads=N_beads,\
                                                                                n_chroms=n_chrom,threshold=10,viz=False,\
                                                                                path=self.save_path)
            elif(loop_path.endswith('.txt')):
                self.ms, self.ns = import_mns_from_txt(txt_file=loop_path,N_beads=N_beads,n_chroms=n_chrom,path=self.save_path)
                self.ds,self.ds,self.cs=None,None,None
            self.N_beads = N_beads
        else:
            # Import from file
            self.Cs = np.load(self.save_path+'genomewide_signal.npy')
            self.ms, self.ns = np.load(self.save_path+'ms.npy'), np.load(self.save_path+'ns.npy')
            self.ks, self.ds = np.load(self.save_path+'ks.npy'), np.load(self.save_path+'ds.npy')
            self.chr_ends = np.load(self.save_path+'chrom_lengths.npy')

            # Estimate number of beads only if needed
            if np.all(self.Cs!=None):
                self.N_beads = len(self.Cs)
            elif np.all(self.ms!=None) and np.all(self.ns!=None):
                self.N_beads = np.max(self.ns)+1
        write_chrom_colors(self.chr_ends,name=self.save_path+'MultiEM_chromosome_colors.cmd')

        # Define a chromsome metric
        self.chroms = np.zeros(self.N_beads)
        chr_count = 1
        for i in range(len(self.chr_ends)-1):
            if not comp_gw: self.chroms[self.chr_ends[i]:self.chr_ends[i+1]] = chr_count
            chr_count += 1

    def add_forcefield(self):
        # Leonard-Jones potential for excluded volume
        self.ev_force = mm.CustomNonbondedForce('epsilon*((sigma1+sigma2)/r)^12')
        self.ev_force.addGlobalParameter('epsilon', defaultValue=100)
        self.ev_force.addPerParticleParameter('sigma')
        for i in range(self.system.getNumParticles()):
            self.ev_force.addParticle([0.1])
        self.system.addForce(self.ev_force)
        
        # Gaussian compartmentalization potential
        if np.all(self.Cs!=None) and len(np.unique(self.Cs)==2):
            self.comp_force = mm.CustomNonbondedForce('E0+E*exp(-(r-r0)^2/(2*sigma^2)); E=(Ea*delta(s1+1)*delta(s2+1)+Eb*delta(s1-1)*delta(s2-1))*delta(chrom1-chrom2)')
            self.comp_force.addGlobalParameter('sigma',defaultValue=0.6)
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
            self.comp_force = mm.CustomNonbondedForce('E0+E*exp(-(r-r0)^2/(2*sigma^2)); E=(Ea1*delta(s1-2)*delta(s2-2)+Ea2*delta(s1-1)*delta(s2-1)+1/2*(Ea1+Ea2)*(delta(s1-2)*delta(s2-1)+delta(s1-1)*delta(s2-2))+Eb1*delta(s1+1)*delta(s2+1)+Eb2*delta(s1+2)*delta(s2+2)+1/2*(Eb1+Eb2)*(delta(s1+2)*delta(s2+1)+delta(s1+1)*delta(s2+2)))*delta(chrom1-chrom2)')
            self.comp_force.addGlobalParameter('sigma',defaultValue=0.6)
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
        radius = self.N_beads/50000**(1/3)*6
        self.container_force = mm.CustomExternalForce(
                '{}*max(0, r-{})^2; r=sqrt((x-{})^2+(y-{})^2+(z-{})^2)'.format(1000,radius,self.mass_center[0],self.mass_center[1],self.mass_center[2]))
        for i in range(self.system.getNumParticles()):
            self.container_force.addParticle(i, [])
        self.system.addForce(self.container_force)

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

    def run_pipeline(self,MD_steps=10000,run_MD=True,write_files=False,plots=False,build_init_struct=True,Temperature=300*mm.unit.kelvin,init_struct_path=None):
        # Initialize simulation
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
        platform = mm.Platform.getPlatformByName('CUDA')
        simulation = Simulation(pdb.topology, self.system, integrator, platform)
        simulation.context.setPositions(pdb.positions)
        start_time = time.time()
        simulation.minimizeEnergy()
        state = simulation.context.getState(getPositions=True)
        PDBxFile.writeFile(pdb.topology, state.getPositions(), open(self.save_path+'MultiEM_minimized.cif', 'w'))
        print(f"--- Energy minimization done!! Executed in {(time.time() - start_time)/60:.2f} minutes. :D ---")

        # Run molecular Dynamics
        if run_MD:
            print('\nRunning MD simulation...')
            start = time.time()
            simulation.reporters.append(DCDReporter(self.save_path+'/MultiEM_traj.dcd', 5))
            simulation.reporters.append(StateDataReporter(stdout, MD_steps//10, step=True, totalEnergy=True, potentialEnergy=True, temperature=True))
            simulation.context.setVelocitiesToTemperature(Temperature, 0)
            simulation.step(MD_steps)
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
    bw_path = '/mnt/raid/data/single_cell/coverage/k562.ATAC.merge.G1.RPGC.bw'
    loop_path = '/mnt/raid/data/single_cell/PET_cluster_with_interchr/k562.G1.all.bedpe'
    out_path_name = 'G1_k_small_structure-weak_interactions'
    
    # Run simulation
    md = MultiEM(N_beads=50000,loop_path=loop_path,comp_path=bw_path,out_path=out_path_name,n_chrom=24)
    md.run_pipeline(run_MD=False,build_init_struct=True,
                    init_struct_path=None,plots=False)

if __name__=='__main__':
    main()