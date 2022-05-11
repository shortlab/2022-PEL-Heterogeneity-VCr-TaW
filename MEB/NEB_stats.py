# authomatic Nudget Elastic Band Method (NEB) for vacancy migration energy barrier with LAMMPS
# Stefano Segantin
# Politecnico di Torino, May 2021

# ---------- DISCUSSION ------------------------------------------------------------------------------------------------
# this routine relies on the assumption that the position in the box of all the atoms rests the same in the intial
# simulation and in the NEB simulation. It is not always like that. In order to do that you need to use the same number
# of parallelized core in both the simulations

########################################################################################################################
# import section
import numpy as np

class NEBstats:

    def __init__(self, filename):
        #self.maindir = maindir
        self.filename = filename

# function that reads the initial data file and returns all the main informations
    def numAtoms(self):
        filename = self.filename
        num_atoms = "Couldn't get the total number of atom in the system"
        with open(filename, 'r') as initial_file:
            initial = initial_file.readlines()
        # get the index at which the list of atoms begins
        initial = [line.strip("\n") for line in initial]
        initial = [line.split(" ") for line in initial]
        # get the first part of the data: number of atoms, number of types, box size and atom masses
        linecounter = -1
        for line in initial:
            linecounter += 1
            if "atoms" in line:
                num_atoms = int(line[0])
                break
        return num_atoms

# function that gets the size of the system box in the three directions
    def boxSize(self):
        filename = self.filename
        xlo = "Couldn't find the lower x bound"
        xhi = "Couldn't find the upper x bound"
        ylo = "Couldn't find the lower y bound"
        yhi = "Couldn't find the upper y bound"
        zlo = "Couldn't find the lower z bound"
        zhi = "Couldn't find the upper z bound"
        with open(filename, 'r') as initial_file:
            initial = initial_file.readlines()
        # get the index at which the list of atoms begins
        initial = [line.strip("\n") for line in initial]
        initial = [line.split(" ") for line in initial]
        # get the first part of the data: number of atoms, number of types, box size and atom masses
        linecounter = -1
        for line in initial:
            linecounter += 1
            if "xlo" in line:
                xlo = float(line[0])
                xhi = float(line[1])
            elif "ylo" in line:
                ylo = float(line[0])
                yhi = float(line[1])
            elif "zlo" in line:
                zlo = float(line[0])
                zhi = float(line[1])
                break
        box_size = [xlo, xhi, ylo, yhi, zlo, zhi]
        return box_size

# function that gets the number of atom types in the box and also the masses of that types
    def atomTypes(self):
        filename = self.filename
        atom_types = "Couldn't find the number of atom types"
        atom_mass = "Couldn't find the mass of the atom types"
        with open(filename, 'r') as initial_file:
            initial = initial_file.readlines()
        # get the index at which the list of atoms begins
        initial = [line.strip("\n") for line in initial]
        initial = [line.split(" ") for line in initial]
        # get the first part of the data: number of atoms, number of types, box size and atom masses
        linecounter = -1
        for line in initial:
            linecounter += 1
            if "atom" in line and "types" in line:
                atom_types = int(line[0])
            elif "Masses" in line:
                masses_index = linecounter
                break
        atom_mass = []
        for types in range(atom_types):
            line = initial[masses_index+2+types]
            line[0] = int(line[0])
            line[1] = float(line[1])
            atom_mass.append(line)
        return atom_types, atom_mass

# function that gives a list of coordinates of all the atoms in the system
    def atomCoords(self):
        filename = self.filename
        with open(filename, 'r') as initial_file:
            initial = initial_file.readlines()
        # get the index at which the list of atoms begins
        for line in enumerate(initial):
            if line[1] == "Atoms # atomic\n":
                zero_line = line[0]
            if line[1] == "Velocities\n":
                end_line = line[0]
                break
        initial = [line.strip("\n") for line in initial]
        initial = [line.split(" ") for line in initial]
        # group the atom coords lines
        atom_coords = initial[zero_line + 2:end_line - 1]
        atom_coords = [[float(element) for element in line] for line in atom_coords]

        return atom_coords

    # function that gives a list of velocities of all the atoms in the system
    def atomVelocities(self):
        filename = self.filename
        with open(filename, 'r') as initial_file:
            initial = initial_file.readlines()
        # get the index at which the list of atoms begins
        for line in enumerate(initial):
            if line[1] == "Atoms # atomic\n":
                zero_line = line[0]
            if line[1] == "Velocities\n":
                end_line = line[0]
                break
        initial = [line.strip("\n") for line in initial]
        initial = [line.split(" ") for line in initial]
        # group the atom velocity lines
        atom_velocities = initial[end_line+2:]
        atom_velocities = [[float(element) for element in line] for line in atom_velocities]
        return atom_velocities

    # create functions for:
        # get the microstructure (BCC, FCC, HPC)
        # get the lattice constant (ao)

# function that gives the line of a random atom chosen inside a box that is smaller than and inside the whole system box
    def randAtomID(self, num_atoms, atom_coords, box_size, ao, aofactor=0):
        self.num_atoms = num_atoms
        self.atom_coords = atom_coords
        self.box_size = box_size
        self.ao = ao  # lattice constant in A
        self.aofactor = aofactor  # multiplication factor needed for generating the boundaries

        # construction of the box
        x_lowerbound = box_size[0] + np.multiply(aofactor, ao)
        x_upperbound = box_size[1] - np.multiply(aofactor, ao)
        y_lowerbound = box_size[2] + np.multiply(aofactor, ao)
        y_upperbound = box_size[3] - np.multiply(aofactor, ao)
        z_lowerbound = box_size[4] + np.multiply(aofactor, ao)
        z_upperbound = box_size[5] - np.multiply(aofactor, ao)
        # generate a list of random ids for generating the vacancy
        flag = 0
        while flag != 1:
            # pick a random atom for the study
            rnd_id = np.random.randint(0, num_atoms)
            rnd_line = atom_coords[rnd_id - 1]
            if rnd_line[2] >= x_lowerbound and rnd_line[2] <= x_upperbound and rnd_line[3] >= y_lowerbound \
                    and rnd_line[3] <= y_upperbound and rnd_line[4] >= z_lowerbound and rnd_line[4] <= z_upperbound:
                random_atom = rnd_line
                flag = 1

        return random_atom

    # function that gives the atoms that are near
    # WARNING: function does not mind whether atom is on the borders, use randAtom with a smaller core box to be sure your
    # atom is not on the border
    def atomNN(self, atom_coords, microstructure, ao, atomID, coordTol=0.1):
        self.atom_coords = atom_coords
        self.microstructure = microstructure  # microstructure (BCC, FCC)
        self.ao = ao  # lattice constant
        self.atomID = atomID
        self.coordTol = coordTol

        # find the atom
        for line in atom_coords:
            if line[0] == atomID:
                myatom_line = line
                myatom_id = myatom_line[0]
                myatom_type = myatom_line[1]
                myatom_coords = [myatom_line[0], myatom_line[3], myatom_line[4]]
                break

        # FCC case of NN
        if microstructure == "FCC":
            # in the FCC microstructure each atom could be considered at the center of a face of the elementary cube
            # so its nearest neighbors are 12 and are at a distance of sqrt(2)/2*ao
            # but in terms of coords (x, y, z) the distance is a0/2 in two out of the three coords
            # idea coordinates of NN atom
            # plane x = x
            NN1_coords = [myatom_line[2], myatom_line[3] - ao / 2, myatom_line[4] - ao / 2]  # atom 1: x, y-ao/2, z-ao/2
            NN2_coords = [myatom_line[2], myatom_line[3] + ao / 2, myatom_line[4] - ao / 2]  # atom 2: x, y+ao/2, z-ao/2
            NN3_coords = [myatom_line[2], myatom_line[3] + ao / 2, myatom_line[4] + ao / 2]  # atom 3: x, y+ao/2, z+ao/2
            NN4_coords = [myatom_line[2], myatom_line[3] - ao / 2, myatom_line[4] + ao / 2]  # atom 4: x, y-ao/2, z+ao/2
            # plane y = y
            NN5_coords = [myatom_line[2] - ao / 2, myatom_line[3], myatom_line[4] - ao / 2]  # atom 5: x-ao/2, y, z-ao/2
            NN6_coords = [myatom_line[2] + ao / 2, myatom_line[3], myatom_line[4] - ao / 2]  # atom 6: x+ao/2, y, z-ao/2
            NN7_coords = [myatom_line[2] + ao / 2, myatom_line[3], myatom_line[4] + ao / 2]  # atom 7: x+ao/2, y, z+ao/2
            NN8_coords = [myatom_line[2] - ao / 2, myatom_line[3], myatom_line[4] + ao / 2]  # atom 8: x-ao/2, y, z+ao/2
            # plane z = z
            NN9_coords = [myatom_line[2] - ao / 2, myatom_line[3] - ao / 2, myatom_line[4]]  # atom 9: x-ao/2, y-ao/2, z
            NN10_coords = [myatom_line[2] + ao / 2, myatom_line[3] - ao / 2, myatom_line[4]]  # atom 10: x+ao/2, y-ao/2, z
            NN11_coords = [myatom_line[2] + ao / 2, myatom_line[3] + ao / 2, myatom_line[4]]  # atom 11: x+ao/2, y+ao/2, z
            NN12_coords = [myatom_line[2] - ao / 2, myatom_line[3] + ao / 2, myatom_line[4]]  # atom 12: x-ao/2, y+ao/2, z
            # list of ideal NN coords
            ideal_NN = [NN1_coords, NN2_coords, NN3_coords, NN4_coords, NN5_coords, NN6_coords, NN7_coords, NN8_coords,
                        NN9_coords, NN10_coords, NN11_coords, NN12_coords]
            # find all the NN
            coord_tol = coordTol * ao
            NN_list = []
            for line in atom_coords:
                for NN in ideal_NN:
                    x_diff = abs(line[2] - NN[0])
                    y_diff = abs(line[3] - NN[1])
                    z_diff = abs(line[4] - NN[2])
                    if x_diff <= coord_tol and y_diff <= coord_tol and z_diff <= coord_tol:
                        NN_list.append(line)
                        break

        elif microstructure == "BCC":
            # in the BCC microstructure each atom could be considered at the center of an elementary cube
            # so its nearest neighbors are 8 and are at a distance of sqrt(3)/2*ao
            # but in terms of coords (x, y, z) the distance is a0/2 in all the three coordinates
            # idea coordinates of NN atom
            # plane -z/2
            NN1_coords = [myatom_line[2] - ao / 2, myatom_line[3] - ao / 2, myatom_line[4] - ao / 2]  # atom 1: x, y-ao/2, z-ao/2
            NN2_coords = [myatom_line[2] + ao / 2, myatom_line[3] - ao / 2, myatom_line[4] - ao / 2]  # atom 2: x, y+ao/2, z-ao/2
            NN3_coords = [myatom_line[2] + ao / 2, myatom_line[3] + ao / 2, myatom_line[4] - ao / 2]  # atom 3: x, y+ao/2, z+ao/2
            NN4_coords = [myatom_line[2] - ao / 2, myatom_line[3] + ao / 2, myatom_line[4] - ao / 2]  # atom 4: x, y-ao/2, z+ao/2
            # plane +z/2
            NN5_coords = [myatom_line[2] - ao / 2, myatom_line[3] - ao / 2, myatom_line[4] + ao / 2]  # atom 1: x, y-ao/2, z-ao/2
            NN6_coords = [myatom_line[2] + ao / 2, myatom_line[3] - ao / 2, myatom_line[4] + ao / 2]  # atom 2: x, y+ao/2, z-ao/2
            NN7_coords = [myatom_line[2] + ao / 2, myatom_line[3] + ao / 2, myatom_line[4] + ao / 2]  # atom 3: x, y+ao/2, z+ao/2
            NN8_coords = [myatom_line[2] - ao / 2, myatom_line[3] + ao / 2, myatom_line[4] + ao / 2]  # atom 4: x, y-ao/2, z+ao/2
            # list of ideal NN coords
            ideal_NN = [NN1_coords, NN2_coords, NN3_coords, NN4_coords, NN5_coords, NN6_coords, NN7_coords, NN8_coords]
            # find all the NN
            coord_tol = coordTol * ao
            NN_list = []
            for line in atom_coords:
                for NN in ideal_NN:
                    x_diff = abs(line[2] - NN[0])
                    y_diff = abs(line[3] - NN[1])
                    z_diff = abs(line[4] - NN[2])
                    if x_diff <= coord_tol and y_diff <= coord_tol and z_diff <= coord_tol:
                        NN_list.append(line)
                        break

        return NN_list

    # function that writes a final.**** file for the NEB method
    def writeFinalFile(self, filename, coordsline, moverline):
        self.filename = filename
        self.coordsline = coordsline
        self.moverline = moverline

        finalfile_text = (r"""1
%s %s %s %s""" % (int(moverline[0]), coordsline[2], coordsline[3], coordsline[4]))

        with open(filename, 'w') as input_file:
            input_file.write(finalfile_text)