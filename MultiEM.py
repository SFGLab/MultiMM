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
    def __init__(self,chrom_ends,Cs=None,ms=None,ns=None,ks=None,ds=None,comp_gw=False,path='results/'):
        '''
        Cs (np arrays): Cs array with colors over time.
        N_beads (int): The number of beads of initial structure.
        step (int): sampling rate.
        '''
        if np.all(Cs!=None):
            self.N_beads = len(Cs)
        elif np.all(ms!=None) and np.all(ns!=None):
            self.N_beads = np.max(ns)+1
        else:
            raise InterruptedError('Wrong data provided in simulation.')
        self.Cs = Cs
        self.ms, self.ns, self.ks, self.ds = ms, ns, ks, ds
        self.chr_ends = chrom_ends
        self.chroms = np.zeros(self.N_beads)
        chr_count = 1
        for i in range(len(self.chr_ends)-1):
            if comp_gw: self.chroms[self.chr_ends[i]:self.chr_ends[i+1]] = chr_count
            chr_count += 1
        self.path = path
        try: 
            os.mkdir(self.path)
            os.mkdir(self.path+'chromosomes')
        except OSError as error:
            print("Folder 'chromosomes' already exists!")

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
            self.comp_force = mm.CustomNonbondedForce('E0+E*exp(-(r-r0)^2/(2*sigma^2)); E=(Ea*delta(s1+s2-2)+Eb*delta(s1+s2+2))*delta(chrom1-chrom2)')
            self.comp_force.addGlobalParameter('sigma',defaultValue=0.5)
            self.comp_force.addGlobalParameter('r0',defaultValue=0.2)
            self.comp_force.addGlobalParameter('E0',defaultValue=0.0)
            self.comp_force.addGlobalParameter('Ea',defaultValue=-2.0)
            self.comp_force.addGlobalParameter('Eb',defaultValue=-0.5)
            self.comp_force.addPerParticleParameter('s')
            self.comp_force.addPerParticleParameter('chrom')
            for i in range(self.system.getNumParticles()):
                self.comp_force.addParticle([self.Cs[i],self.chroms[i]])
            self.system.addForce(self.comp_force)
        elif np.all(self.Cs!=None) and len(np.unique(self.Cs)==4):
            self.comp_force = mm.CustomNonbondedForce('E0+E*exp(-(r-r0)^2/(2*sigma^2)); E=(Ea1*delta(s1+s2-2)+Ea2*delta(s1+s2-4)+Eb1*delta(s1+s2+2)+Eb2*delta(s1+s2+4))*delta(chrom1-chrom2)')
            self.comp_force.addGlobalParameter('sigma',defaultValue=0.5)
            self.comp_force.addGlobalParameter('r0',defaultValue=0.2)
            self.comp_force.addGlobalParameter('E0',defaultValue=0.0)
            self.comp_force.addGlobalParameter('Ea1',defaultValue=-6.0)
            self.comp_force.addGlobalParameter('Ea2',defaultValue=-5.0)
            self.comp_force.addGlobalParameter('Eb1',defaultValue=-4.0)
            self.comp_force.addGlobalParameter('Eb2',defaultValue=-3.0)
            self.comp_force.addPerParticleParameter('s')
            self.comp_force.addPerParticleParameter('chrom')
            for i in range(self.system.getNumParticles()):
                self.comp_force.addParticle([self.Cs[i],self.chroms[i]])
            self.system.addForce(self.comp_force)
        
        # Spherical container
        radius = np.sqrt(self.N_beads/50000)*6
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
            self.angle_force.addAngle(i, i+1, i+2, np.pi, 100)
        self.system.addForce(self.angle_force)

    def run_pipeline(self,MD_steps=10000,run_MD=True,write_files=False,plots=False,build_init_struct=True,Temperature=300*mm.unit.kelvin,init_struct_path=None):
        # Initialize simulation
        if build_init_struct:
            print('\nCreating initial structure...')
            if np.all(self.Cs!=None): write_cmm(self.Cs,name=self.path+'MultiEM_compartment_colors.cmd')
            pdb_content = build_init_mmcif(n_dna=self.N_beads,path=self.path)
            print('---Done!---')
        pdb = PDBxFile(self.path+'MultiEM_init.cif') if init_struct_path==None or build_init_mmcif else PDBxFile(init_struct_path)
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
        PDBxFile.writeFile(pdb.topology, state.getPositions(), open(self.path+'MultiEM_minimized.cif', 'w'))
        print(f"--- Energy minimization done!! Executed in {(time.time() - start_time)/60:.2f} minutes. :D ---")

        # Run molecular Dynamics
        if run_MD:
            print('\nRunning MD simulation...')
            start = time.time()
            simulation.reporters.append(DCDReporter(self.path+'/MultiEM_traj.dcd', 5))
            simulation.reporters.append(StateDataReporter(stdout, MD_steps//10, step=True, totalEnergy=True, potentialEnergy=True, temperature=True))
            simulation.context.setVelocitiesToTemperature(Temperature, 0)
            simulation.step(MD_steps)
            end = time.time()
            elapsed = end - start
            speed = MD_steps / elapsed
            print(f"\n---MD finished in {elapsed/60:.2f} minutes ({speed:0.1f} steps/s)---")
            
            state = simulation.context.getState(getPositions=True)
            PDBxFile.writeFile(pdb.topology, state.getPositions(), open(self.path+'MultiEM_afterMD.cif', 'w'))

        start = time.time()
        V = get_coordinates_mm(state.getPositions())
        save_metrics(V,path_name=self.path+'GW_')
        for i in range(len(self.chr_ends)-1):
            write_mmcif(coords=10*V[self.chr_ends[i]:self.chr_ends[i+1]],
                        path=self.path+f'chromosomes/MultiEM_minimized_{chrs[i]}.cif')
            save_metrics(10*V[self.chr_ends[i]:self.chr_ends[i+1]],path_name=self.path+f'chromosomes/{chrs[i]}_')
        
        if plots:
            print('\nComputing heatmap...')
            heat = get_heatmap(mm_vec=state.getPositions(),viz=plots,path=self.path)
            end = time.time()
            elapsed = end - start
            print(f'---Heatmap computed in {elapsed/60:0.2f} minutes.---')

def main():
    N_beads, n_chrom = 200000, 24
    path = 'k562_S_combined_all_kd_big_structure/'
    try:
        os.mkdir(path)
        os.mkdir(path+'chromosomes')
    except OSError as error:
        print("Folder 'chromosomes' already exists!")

    # # Load from files
    # eigenvec, chr_ends = import_bw(bw_path="/mnt/raid/data/single_cell/coverage/k562.ATAC.merge.G2M.RPKM.bw",
    #                                N_beads=N_beads,viz=False,binary=True,n_chroms=n_chrom,path=path)
    # ms, ns, ds, ks, cs = import_mns_from_bedpe(bedpe_file='/mnt/raid/data/single_cell/PET_cluster_with_interchr/k562.G2M.all.bedpe',
    #                                    N_beads=N_beads,n_chroms=n_chrom,threshold=10,viz=False,path=path)

    # Load from numpy arrays
    ms = np.load(path+'ms.npy')#[::10]
    ns = np.load(path+'ns.npy')#[::10]
    ks = np.load(path+'ks.npy')
    ds = np.load(path+'ds.npy')
    print('Number of loops:',len(ms))
    chr_ends = np.load(path+'chrom_lengths.npy')
    eigenvec = np.load(path+'genomewide_signal.npy')

    # Run simulation
    md = MultiEM(Cs=eigenvec,chrom_ends=chr_ends,ms=ms,ns=ns,ks=ks,ds=None,path=path)
    md.run_pipeline(run_MD=False,build_init_struct=True,
                    init_struct_path=None,plots=False)#'/mnt/raid/codes/mine/MultiEM-main/init_structures/N_beads_50000/MultiEM_init_n50000.cif')

if __name__=='__main__':
    main()
