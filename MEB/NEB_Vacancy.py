# authomatic Nudget Elastic Band Method (NEB) for vacancy migration energy barrier with LAMMPS
# Stefano Segantin
# Politecnico di Torino

# ---------- DISCUSSION ------------------------------------------------------------------------------------------------
# this routine relies on the assumption that the position in the box of all the atoms rests the same in the intial
# simulation and in the NEB simulation. It is not always like that. In order to do that you need to use the same number
# of parallelized core in both the simulations

########################################################################################################################
# import section
import os
import numpy as np
import NEB_stats as neb
import NEB_readout as nebro
import pandas as pd

myDir = os.getcwd()
os.chdir(myDir)
initialfile = r'data.***.0'
myneb = neb.NEBstats(initialfile)  # this just opens the
mynebro = nebro.NEBreadout()


microstructure = "BCC"
ao = 3.31642067489331  # lattice constant
N_rnd_atoms = 1  # number of random atoms to pick
finalfile = r"final.***"
lmpInFile = r"***_neb_py.in"
nebOutFile = r"***_neb_py.xlsx"

vacancies = []
neighbors = []
MainInfo = []
atom_coords = myneb.atomCoords()
box_size = myneb.boxSize()
num_atoms = myneb.numAtoms()
MI = []
counter = 0
for sample in range(N_rnd_atoms):
    print("NUMBER OF VACANCY SAMPLED = %s" % (sample))
    rndvacancy = myneb.randAtomID(num_atoms, atom_coords, box_size, ao, 1)
    vacancyID = int(rndvacancy[0])
    neighbors = myneb.atomNN(atom_coords, microstructure, ao, vacancyID)
    vacancies.append(rndvacancy)
    NNs = []
    for NN in neighbors:
        file = myneb.writeFinalFile(finalfile, rndvacancy, NN)
        NNs.append(NN)

        txtin = (r""" # Nudged elastic band (NEB) analysis
#Politecnico di Torino

# ---------- INITIALIZATION ----------------------------------------------------
clear
timer		full
units		metal  # A ps K bar eV
dimension	3
boundary	p p p
atom_style	atomic
atom_modify	map array

variable	ao equal %s
# ---------- GEOMETRY DEFINITION -----------------------------------------------
lattice	bcc ${ao}
read_data %s

group

# ---------- create vacancy ----------------------------------------------------
group   vacancy id %s
delete_atoms group vacancy compress no

# ---------- FORCE FIELD -------------------------------------------------------
pair_style
pair_coeff	* *
neigh_modify    every 1 delay 0 check yes

# ---------- dumb run ----------------------------------------------------------
run		0
# ---------- SETTINGS ----------------------------------------------------------
# COMPUTES
compute	csym all centro/atom bcc
compute	eng all pe/atom
compute	eatoms all reduce sum c_eng
# VARIABLES
variable	N equal count(all)
variable	Eavg equal c_eatoms/$N

# ---------- DYNAMICS ----------------------------------------------------------
# ---------- energy minimization -----------------------------------------------
# PRINT-SCREEN
thermo
thermo_style
displace_atoms all random
minimize

# --------- NEB analysis -------------------------------------------------------
variable	u uloop
reset_timestep	0
fix		1 all neb 
timestep        
min_style	quickmin
neb		 final %s
""" % (str(ao), str(initialfile), str(vacancyID), str(finalfile)))

        with open(lmpInFile, 'w') as input_file:
            input_file.write(txtin)

        os.system('mpiexec -np 10 lmp_mpi -partition 10x1 -in %s' % lmpInFile)
        # get the info about the climbing replica
        Creplica_data, Creplica_num = mynebro.climbingReplica()
        Ebarrier = mynebro.energyBarrier(Creplica_data[-1])
        MI.append([rndvacancy[0], rndvacancy[1], NN[0], NN[1], Ebarrier])

    neighbors.append(NNs)
    MainInfo.append(MI)

MI = np.array(MI)
df = pd.DataFrame(MI)
df.to_excel(excel_writer=nebOutFile)