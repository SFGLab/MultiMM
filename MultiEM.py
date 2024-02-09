#########################################################################
########### CREATOR: SEBASTIAN KORSAK, WARSAW 2024 ######################
#########################################################################

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
from nucs_init_struct_tools import *
from nucs_preprocessing import *
import time
import os

class MultiEM:
    def __init__(self,N_beads,loop_path,comp_path=None,nucs_path=None,chrom=None,coords=None,comp_gw=True,out_path='results'):
        '''
        Input data:
        ------------
        N_beads (int): Number of simulation beads.
        loop_path (str): path of the .bedpe file for the looping data (required).
        comp_path (str): path of the .bw or .bed file for the compartmentalization data (optional).
        nucs_path (str): path of the PuFFIN file for the nucleosome positions (optional - it works only for small fragments of chromosome).
        out_path (str): path of saving data.
        '''
        ################################################################
        ################ CREATE OUTPUT FOLDER ##########################
        ################################################################
        self.save_path = out_path+'/'
        self.ms, self.ns, self.ks, self.chr_ends, self.Cs = None, None, None, None, None
        try:
            os.mkdir(self.save_path)
            os.mkdir(self.save_path+'chromosomes')
            os.mkdir(self.save_path+'chromosomes_info')
        except OSError as error:
            print("Folder 'chromosomes' already exists!")
        
        self.chrom, self.coords = chrom, coords
        
        ################################################################
        ################ LOAD COMPARTMENTS #############################
        ################################################################
        if comp_path.endswith('.bw') or comp_path.endswith('.bigwig') or comp_path.endswith('.BigWig'):
            self.Cs, self.chr_ends = import_bw(bw_path=comp_path,N_beads=N_beads,\
                                               viz=False,binary=True,\
                                               path=self.save_path)
        elif comp_path.endswith('.bed'):
            self.Cs, self.chr_ends = import_compartments_from_Calder(bed_file=comp_path,N_beads=N_beads,\
                                                                     chrom=self.chrom,coords=self.coords,\
                                                                     save_path=self.save_path)
        elif comp_path=='random':
            self.Cs = random_blocks()
        else:
            self.Cs = None
        
        ################################################################
        ################ LOAD LOOPS ####################################
        ################################################################
        if loop_path.endswith('.bedpe'):
            self.ms, self.ns, self.ks, self.chr_ends = import_mns_from_bedpe(bedpe_file=loop_path,N_beads=N_beads,\
                                                                             coords = self.coords, chrom=self.chrom,\
                                                                             viz=False, path=self.save_path)
        elif loop_path=='random':
            self.ms, self.ns, self.ks = generate_arrays(N_loops=N_beads//8, N=N_beads)
        else:
            raise InterruptedError('You did not provide appropriate loop file. Loop .bedpe file is obligatory.')

        self.N_beads = N_beads
        
        write_chrom_colors(self.chr_ends,name=self.save_path+'MultiEM_chromosome_colors.cmd')

        ################################################################
        ################ LOAD NUCLEOSOMES ##############################
        ################################################################
        if nucs_path!=None:
            self.entry_points, self.exit_points = puffin_to_array(nucs_path)
            self.nuc_sim_len, self.num_nucs = int(self.exit_points[-1]+1), len(self.entry_points)
            if self.nuc_sim_len>1e7: raise InterruptedError('Too large region or structure to model nucleosomes. Please choose smaller region or simulation beads. ;)')
            print('Nucleosome simulation length:',self.nuc_sim_len)
            print('Number of nucleosomes',self.num_nucs)
        else:
            self.entry_points, self.exit_points, self.nuc_sim_len, self.num_nucs = None, None, 0, 0
        
        # if np.all(comp_path!=None): self.Cs = align_comps(self.Cs,self.ms,self.chr_ends)
        
        ################################################################
        ################ LOAD CHROMOSOMES ##############################
        ################################################################
        self.chroms, self.chroms_AB = np.zeros(self.N_beads), np.zeros(self.N_beads)
        if chrom==None:
            chr_count = 1
            for i in range(len(self.chr_ends)-1):
                if not comp_gw: self.chroms_AB[self.chr_ends[i]:self.chr_ends[i+1]] = chr_count
                self.chroms[self.chr_ends[i]:self.chr_ends[i+1]] = chr_count
                chr_count += 1


    def add_forcefield(self):
        '''
        Here we define the forcefield of MultiEM.
        '''
        # Leonard-Jones potential for excluded volume
        self.ev_force = mm.CustomNonbondedForce('epsilon*((sigma1+sigma2)/r)^6')
        self.ev_force.addGlobalParameter('epsilon', defaultValue=100)
        self.ev_force.addPerParticleParameter('sigma')
        for i in range(self.system.getNumParticles()):
            self.ev_force.addParticle([0.05])
        self.system.addForce(self.ev_force)
        
        # Gaussian compartmentalization potential
        if np.all(self.Cs!=None) and len(np.unique(self.Cs)==2):
            self.comp_force = mm.CustomNonbondedForce('E0+E*exp(-(r-r0)^2/(2*sigma^2)); E=(Ea*delta(s1-1)*delta(s2-1)+Eb*delta(s1+1)*delta(s2+1))*delta(chrom1-chrom2)')
            self.comp_force.addGlobalParameter('sigma',defaultValue=0.5)
            self.comp_force.addGlobalParameter('r0',defaultValue=0.1)
            self.comp_force.addGlobalParameter('E0',defaultValue=0.0)
            self.comp_force.addGlobalParameter('Ea',defaultValue=-1)
            self.comp_force.addGlobalParameter('Eb',defaultValue=-2)
            self.comp_force.addPerParticleParameter('s')
            self.comp_force.addPerParticleParameter('chrom')
            for i in range(self.system.getNumParticles()):
                self.comp_force.addParticle([self.Cs[i],self.chroms_AB[i]])
            self.system.addForce(self.comp_force)
        elif np.all(self.Cs!=None) and len(np.unique(self.Cs)>=4):
            self.comp_force = mm.CustomNonbondedForce('E0+E*exp(-(r-r0)^2/(2*sigma^2)); E=(Ea1*delta(s1-2)*delta(s2-2)+Ea2*delta(s1-1)*delta(s2-1)+Eb1*delta(s1+1)*delta(s2+1)+Eb2*delta(s1+2)*delta(s2+2)')
            self.comp_force.addGlobalParameter('sigma',defaultValue=0.5)
            self.comp_force.addGlobalParameter('r0',defaultValue=0.1)
            self.comp_force.addGlobalParameter('E0',defaultValue=0.0)
            self.comp_force.addGlobalParameter('Ea1',defaultValue=-0.5)
            self.comp_force.addGlobalParameter('Ea2',defaultValue=-1)
            self.comp_force.addGlobalParameter('Eb1',defaultValue=-1.5)
            self.comp_force.addGlobalParameter('Eb2',defaultValue=-2.0)
            self.comp_force.addPerParticleParameter('s')
            self.comp_force.addPerParticleParameter('chrom')
            for i in range(self.system.getNumParticles()):
                self.comp_force.addParticle([self.Cs[i],self.chroms_AB[i]])
            self.system.addForce(self.comp_force)
        
        # Spherical container
        if self.chrom==None:
            radius1 = (self.N_beads/50000)**(1/3)
            radius2 = (self.N_beads/50000)**(1/3)*6
            self.container_force = mm.CustomExternalForce('C*(max(0, r-R2)^2+max(0, R1-r)^2); r=sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2)')
            self.container_force.addGlobalParameter('C',defaultValue=1000)
            self.container_force.addGlobalParameter('R1',defaultValue=radius1)
            self.container_force.addGlobalParameter('R2',defaultValue=radius2)
            self.container_force.addGlobalParameter('x0',defaultValue=self.mass_center[0])
            self.container_force.addGlobalParameter('y0',defaultValue=self.mass_center[1])
            self.container_force.addGlobalParameter('z0',defaultValue=self.mass_center[2])
            for i in range(self.system.getNumParticles()):
                self.container_force.addParticle(i, [])
            self.system.addForce(self.container_force)

            if np.all(self.Cs!=None):
                # Interaction of A compartment with lamina
                self.Alamina_force = mm.CustomExternalForce('-A*sin(pi*(r-R1)/(R2-R1))^2*(delta(s-1)+delta(s-2))*step(R1)*(1-step(R2)); r=sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2)')
                self.Alamina_force.addGlobalParameter('A',defaultValue=1000)
                self.Alamina_force.addGlobalParameter('pi',defaultValue=3.14159265358979323846)
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
                self.Blamina_force = mm.CustomExternalForce('B*(sin(pi*(r-R1)/(R2-R1))^8-1)*(delta(s+1)+delta(s+2))*step(R1)*(1-step(R2)); r=sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2)')
                self.Blamina_force.addGlobalParameter('B',defaultValue=1000)
                self.Blamina_force.addGlobalParameter('pi',defaultValue=3.14159265358979323846)
                self.Blamina_force.addGlobalParameter('R1',defaultValue=radius1)
                self.Blamina_force.addGlobalParameter('R2',defaultValue=radius2)
                self.Blamina_force.addGlobalParameter('x0',defaultValue=self.mass_center[0])
                self.Blamina_force.addGlobalParameter('y0',defaultValue=self.mass_center[1])
                self.Blamina_force.addGlobalParameter('z0',defaultValue=self.mass_center[2])
                self.Blamina_force.addPerParticleParameter('s')
                for i in range(self.system.getNumParticles()):
                    self.Blamina_force.addParticle(i, [self.Cs[i]])
                self.system.addForce(self.Blamina_force)

            # Force that sets smaller chromosomes closer to the center
            self.central_force = mm.CustomExternalForce('G*(chrom-1)/23*(-1/(r-R1+1)+1/(r-R1+1)^2); r=sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2)')
            self.central_force.addGlobalParameter('G',defaultValue=2000)
            self.central_force.addGlobalParameter('R1',defaultValue=radius1)
            self.central_force.addGlobalParameter('x0',defaultValue=self.mass_center[0])
            self.central_force.addGlobalParameter('y0',defaultValue=self.mass_center[1])
            self.central_force.addGlobalParameter('z0',defaultValue=self.mass_center[2])
            self.central_force.addPerParticleParameter('chrom')
            for i in range(self.system.getNumParticles()):
                self.central_force.addParticle(i, [self.chroms[i]])
            self.system.addForce(self.central_force)

        # Bond force
        self.bond_force = mm.HarmonicBondForce()
        for i in range(self.system.getNumParticles()-1):
            self.bond_force.addBond(i,i+1,0.1,300000)
        self.system.addForce(self.bond_force)
        
        # Add fixed loops force
        if np.all(self.ms==None) and np.all(self.ns==None):
            self.loop_force = mm.HarmonicBondForce()
            counter=0
            for m,n in tqdm(zip(self.ms,self.ns),total=len(self.ms)):
                if np.any(self.ks==None):
                    self.loop_force.addBond(m,n,0.1,30000)
                else:
                    self.loop_force.addBond(m,n,0.05,self.ks[counter])
            self.system.addForce(self.loop_force)

        # Bending potential for stiffness
        self.angle_force = mm.HarmonicAngleForce()
        for i in range(self.system.getNumParticles()-2):
            self.angle_force.addAngle(i, i+1, i+2, np.pi, 40)
        self.system.addForce(self.angle_force)

    def add_nuc_forcefield(self,n_wraps=2):
        '''
        Here we define the forcefield of the nucleosome model.
        '''
        # Bending potential for stiffness
        angle_force = mm.HarmonicAngleForce()
        for i in range(4*self.num_nucs,self.system_n.getNumParticles()-2):
            if self.is_nuc[i-4*self.num_nucs]:
                angle_force.addAngle(i, i+1, i+2, np.pi, 5000)
            else:
                angle_force.addAngle(i, i+1, i+2, np.pi, 100)
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
                dnadna_force.addBond(4*self.num_nucs+j, 4*self.num_nucs+j+r, 0.11, 3000.0)
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
                        histone_dna_force.addBond(i, 4*self.num_nucs+j, 0.175, 500)
                        res.append(f':{k+1}@{atoms[i%4]}\t:{self.num_nucs+j+1}\tgreen\n')
        histone_dna_force.setForceGroup(4)
        self.system_n.addForce(histone_dna_force)
        print('DNA-histone force imported.')

    def run_pipeline(self,MD_steps=10000,write_files=False,build_init_struct=True,Temperature=300*mm.unit.kelvin,init_struct_path=None,pltf='CUDA'):
        '''
        Energy minimization for GW model.
        '''
        # Initialize simulation
        if build_init_struct:
            print('\nCreating initial structure...')
            comp_mode = 'compartments' if np.all(self.Cs!=None) and len(np.unique(self.Cs))<=3 else 'subcompartments'
            if np.all(self.Cs!=None): write_cmm(self.Cs,name=self.save_path+'MultiEM_compartment_colors.cmd')
            pdb_content = build_init_mmcif(n_dna=self.N_beads,path=self.save_path)
            print('---Done!---')
        pdb = PDBxFile(self.save_path+'MultiEM_init.cif') if init_struct_path==None or build_init_mmcif else PDBxFile(init_struct_path)
        self.mass_center = np.average(get_coordinates_mm(pdb.positions),axis=0)
        forcefield = ForceField('ff.xml')
        self.system = forcefield.createSystem(pdb.topology)
        integrator = mm.LangevinIntegrator(Temperature, 0.5, 1 * mm.unit.femtosecond)
        
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
        print(f"Simulation will run on platform: {current_platform.getName()}.")
        start_time = time.time()
        simulation.minimizeEnergy()
        state = simulation.context.getState(getPositions=True)
        PDBxFile.writeFile(pdb.topology, state.getPositions(), open(self.save_path+'MultiEM_minimized.cif', 'w'))
        print(f"--- Energy minimization done!! Executed in {(time.time() - start_time)/60:.2f} minutes. :D ---")

        if self.chrom==None:
            start = time.time()
            V = get_coordinates_mm(state.getPositions())
            for i in range(len(self.chr_ends)-1):
                write_mmcif(coords=10*V[self.chr_ends[i]:self.chr_ends[i+1]],
                            path=self.save_path+f'chromosomes/MultiEM_minimized_{chrs[i]}.cif')
        
        if self.num_nucs>0: self.run_nuc_pipeline()

    def run_nuc_pipeline(self):
        '''
        Energy minimization for nucleosome simulation.
        '''
        # Define system
        print('\n\nBuilding initial structure of nucleosome simulation...')
        pdb_content = build_init_mmcif_nucs(entry=self.entry_points,exit=self.exit_points,n_nucs=self.num_nucs,n_dna=self.nuc_sim_len,
                               mode='path',psf=True,path=self.save_path+'MultiEM_minimized.cif')

        print('Creating system...')
        pdb = PDBxFile('dna_histones.cif')
        forcefield = ForceField('forcefields/dna_histones_ff.xml')
        self.system_n = forcefield.createSystem(pdb.topology, nonbondedCutoff=1*mm.unit.nanometer)
        integrator = mm.LangevinIntegrator(310, 5, 10 * mm.unit.femtosecond)
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
        simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, totalEnergy=True, potentialEnergy=True, temperature=True))
        simulation.context.setPositions(pdb.positions)
        simulation.minimizeEnergy()
        state = simulation.context.getState(getPositions=True)
        PDBxFile.writeFile(pdb.topology, state.getPositions(), open(self.save_path+f'/other/minimized_nucres_model.cif', 'w'))
        print('Energy minimization of nucleosome model done :D\n')

def main():
    # Input data
    bw_path = '/mnt/raid/data/Trios/calder_HiChIP_subcomp/YRB_d.bed'
    loop_path = '/mnt/raid/data/Trios/ChiA-PiPE_Loops/loops_pet3+/GM19240_YRI_C_CTCF_1mb_pet3.bedpe'
    nuc_path = None#'/mnt/raid/codes/other/PuFFIN/input/ENCFF415FEC_rep1_chr1.bed.nucs'
    out_path_name = 'region_test'
    chrom = 'chr4'
    coords = [2e7,10e7]
    
    # Run simulation
    md = MultiEM(N_beads=10000,out_path=out_path_name,loop_path=loop_path,comp_path=bw_path,nucs_path=nuc_path)
    md.run_pipeline(build_init_struct=True,init_struct_path=None,pltf='CUDA')

if __name__=='__main__':
    main()