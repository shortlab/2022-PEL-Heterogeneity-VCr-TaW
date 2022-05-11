# authomatic Nudget Elastic Band Method (NEB) for vacancy migration energy barrier with LAMMPS
# Stefano Segantin
# Politecnico di Torino

# ---------- DISCUSSION ------------------------------------------------------------------------------------------------
# this routine relies on the assumption that the position in the box of all the atoms rests the same in the intial
# simulation and in the NEB simulation. It is not always like that. In order to do that you need to use the same number
# of parallelized core in both the simulations

########################################################################################################################
# import section
import numpy as np

class NEBreadout:

    def __init__(self, filename=r"log.lammps"):
        #self.maindir = maindir
        self.filename = filename

    def is_integer(self, n):
        self.n = n
        try:
            float(n)
        except ValueError:
            return False
        else:
            return float(n).is_integer()

# function that reads the log.lammps master file
    def openLog(self):
        filename = self.filename
        with open(filename, 'r') as initial_file:
            initial = initial_file.readlines()
        initial = [line.strip("\n") for line in initial]
        initial = [line.split(" ") for line in initial]
        for line in initial:
            while '' in line:
                line.remove('')

        return initial

# function that gives the number of replica used
    def numReplicas(self):
        logfile = self.openLog()

        tot_replicas = "Cannot find the number of replicas used"
        for line in logfile:
            if "Running" in line:
                tot_replicas = int(line[2])
                break

        return tot_replicas

    # function that gives info e data about the initial replica
    def replicaData(self):
        logfile = self.openLog()

        for line in enumerate(logfile):
            if line[1][0] == "0":
                initial_line = line[0]
                break

        for line in enumerate(logfile):
            if "Climbing" in line[1]:
                final_line = line[0]
                break

        replica_data = logfile[initial_line:final_line]
        replica_data = [[float(el) for el in line] for line in replica_data]

        return replica_data

# function that gives info e data about the climbing replica
    def climbingReplica(self):
        logfile = self.openLog()

        for line in enumerate(logfile):
            if "Running" in line[1]:
                initial_line = line[0] + 2
            elif "Climbing" in line[1]:
                final_line = line[0] -1
                Creplica_num = int(line[1][3])
                break

        Creplica_data = logfile[initial_line:final_line]
        Creplica_data = [[float(el) for el in line] for line in Creplica_data]

        return Creplica_data, Creplica_num

# function that gets the PE of each replica at a given timestep
# WARNING: timestep must be chosen before casting this function, it gets just the right line
    def getPEs(self, replicadata):
        self.replicadata = replicadata

        replicadata = np.array(replicadata)
        numel = len(replicadata)
        peindex = np.arange(10,numel,2)
        PEs = replicadata[peindex]

        return PEs

# function that gives the climbing energy barrier
# WARNING: timestep must be chosen before casting this function, it gets just the right line
    def energyBarrier(self, replicadata):
        self.replicadata = replicadata

        PEs = self.getPEs(replicadata)  # gets the PEs of the given replica
        En_barrier = max(PEs) - min(PEs)

        return En_barrier
