# 2022-PEL-Heterogeneity-VCr-TaW
Repository for Stefano Segantin's MD manuscript.

This repository provides some example input files and scripts regarding the molecular dynamics modeling of the
Some example results are also provided for each of the models.

A folder for each of the three models presented in the paper is provided:
1) Quasi-static (system characterization)
2) MEB (Migration Energy Barriers) evaluation (NEB method)
3) primary damage (subsequent collision cascade algorithm)

In order to proceed it is necessary to have installed:
- Python 3.6 or higher
- LAMMPS molecular dynamics open source code
- Have LAMMPS installed as Python module

It is necessary to download 

The use High Performance supercomputers is strongly suggested.
It is also strongly suggested to read the [LAMMPS manual](https://docs.lammps.org/Intro.html) in order to understand all the command used and how to tune the command arguments.

Additional information and modeling parameters applied are described in the paper methodology.

## 1) Quasi-static model
The input file of this model is the "quasi_static.in" file in the relative folder. It is just necessary to complete the commands left blank with the requirement arguments with values that can be arbitrary, described in the paper or specific of the material applied. Here you can find listed some main modeling choices made in the work:
- In the work the total [box size](https://docs.lammps.org/region.html) was set equal to 40x40x40 [lattice constants](https://docs.lammps.org/lattice.html) and [created](https://docs.lammps.org/create_box.html) (i.e. 128'000 particles for BCC cristals)
- The [temperature](https://docs.lammps.org/velocity.html) was set equal to 300 K
- The hybrid Monte Carlo + Molecular Dynamics [swap method](https://docs.lammps.org/fix_atom_swap.html) was set to perform 500 swap attempts every 50 timesteps for a total of 1 million swap attempts
- The [timestep](https://docs.lammps.org/timestep.html) was set equal to 0.001 ps
- The simulation [run](https://docs.lammps.org/run.html) in the [npt ensemble](https://docs.lammps.org/fix_nh.html) (300 K, recalled every 0.1 ps) for a total of 100'000 steps
- The simulation outputs a [data](https://docs.lammps.org/write_data.html) or a [read/restart](https://docs.lammps.org/read_restart.html) file with the box configuration at the latest timestep that can be used as input file for the other two models.

It is necessary to define the [type](https://docs.lammps.org/pair_style.html) of interatomic potential adopted, according to the interatomic potential library chosen. It is also required to [point the file path](https://docs.lammps.org/pair_coeff.html) of the interatomic potential library chosen.
[Here](https://docs.lammps.org/Run_head.html) you can find how to run the LAMMPS code starting from a lammps input file.

It is suggested to run the simulation 11 times for each binary mixed system in order to cover all the concentration range with a concentration step of 10%.
It is also possibe to embed the script in a [variable-loop](https://docs.lammps.org/variable.html) command in order to authomatize such parametrization.

## 2) Migration Energy Barriers model
For this model a Python wrapper has been built as well as a couple of Phython functions. The [Nudged Elastic Band (NEB)](https://docs.lammps.org/fix_neb.html) method relies on the fact that, in addition to the classic input file, interatomic potential pair coefficient library file and the initial [data](https://docs.lammps.org/write_data.html), a "final" file is given to LAMMPS in which the IDs of one or more atoms in the box is given as well as the coordinates of such atoms. LAMMPS performs a "forced migration" of such atoms and records the path that requires the minimum energy for the migration to happen. The work studies the heterogeneity of such energy values (with focus on the migration of vacancies) for mixed systems. Hence, it was necessary to perform a huge number of different migrations. The Python files allow to do that automatically. The files work as follows:

- "NEB_stats.py": opens the initial [data](https://docs.lammps.org/write_data.html) or a [read/restart](https://docs.lammps.org/read_restart.html) file and collects some main information (box size, microstructure type, particles velocities, particles ID and coordinates and, for each particle, the ID and coordinates of all the Nearest Neighbor (NN) particles) and writes the "final" file
- "NEB_readout.py": it opens the "log.lammps" file, which is one of the simulation output files and postprocesses the results in order to identify and record the MEB values evaluated in the simulation
- "NEB_Vacancy.py": it is the main file that performs the multiple NEB simulations. It recalls several functions from the "NEB_stats.py" and "NEB_readout.py" files in order to read the model configuration and setup the input file. In generates the actual LAMMPS input file for the NEB simulation alongside the "initial" file that gets generated by a function of the "NEB_stats.py" file here recalled. It records all the important outputs (vacancy ID, NNs IDs, atom type and MEB value) and write them in a xlsx file. The whole process is embedded in a "for" cycle in order to run several simulations for the number of vacancies chosen. 

The simulation algorithm works as follows:
1) The users provides the initial [data](https://docs.lammps.org/write_data.html) or a [read/restart](https://docs.lammps.org/read_restart.html) file and choses the number of vacancies that the program has to try to migrate
2) For each vacancy a random atom ID far from the box boundaries is chosen and removed
3) Its NNs are identified
4) For each of the NNs a "final" file that points out the NN ID and its final position (i.e. the coordinates of the vacancy left by the removed atom) is generated
5) The NEB simulation is then run for each NN of each vacancy
6) A xlsx file providing main informations regarding the NEB simulation and MEB values is generated

To run the simulation it is necessary to put the three Python files and the initial [data](https://docs.lammps.org/write_data.html) or a [read/restart](https://docs.lammps.org/read_restart.html) file in the same folder. Then, it is necessary to open the "NEB_vacancy.py" file and set:
- The microstructure type ("BCC" or "FCC" are available for the moment)
- The lattice constant "ao" in Armstrong
- The number of vacancies to generate (the migration of every NN will be performed each vacancy)
- The name of the "final" file that will be generated has to be chosen
- The name of the LAMMPS input file that will be generated file has to be chosen
- The name of the outplut .xlsx file that will be generated has to be chosen
- It is necessary to modify the LAMMPS input file that the "NEB_Vacancy.py" generates directly from the Python code. It is a string-like region of the code in which it is necessary to insert typical LAMMPS input file parameters ([pair style](https://docs.lammps.org/pair_style.html), [pair coefficents](https://docs.lammps.org/pair_coeff.html), [variables](https://docs.lammps.org/variable.html), [fixes](https://docs.lammps.org/fixes.html) etc.)
- It might be necessary to modify the string that runs the simulation depending on the command line that has to be used to run LAMMPS (currently, the string is set as: "os.system('mpiexec -np 10 lmp_mpi -partition 10x1 -in %s' % lmpInFile)")

## 3) Primary damage model
The input file of this model is the "primary_damage.in" file in the relative folder. It is just necessary to put the initial [data](https://docs.lammps.org/write_data.html) or [read/restart](https://docs.lammps.org/read_restart.html) file in the same folder of the simulation. It is also necessary to complete the commands left blank with the requirement arguments with values that can be arbitrary, described in the paper or specific of the material applied. Here you can find listed some main modeling choices made in the work:
- The [temperature](https://docs.lammps.org/velocity.html) was set equal to 300 K
- A [variable-loop](https://docs.lammps.org/variable.html) command is used to simulate 1500 collisional cascades
- The [data](https://docs.lammps.org/write_data.html) or [read/restart](https://docs.lammps.org/read_restart.html) coming from the quasi-static model or from the previous cascades in the loop is read
- A [variable](https://docs.lammps.org/variable.html) is associated to each random seed number necessary and then recalled multiplied by the loop variable so as to have different seed numbers at each collisional cascade. Seed numbers are required for the primary knock on atom (PKA) choice, for the two angles of the PKA direction and for the temperature distribution
- The [pair style](https://docs.lammps.org/pair_style.html) here is [hybrid](https://docs.lammps.org/pair_hybrid.html) and the [pair coefficents](https://docs.lammps.org/pair_coeff.html) come from the chosen libraries and the [ZBL](https://docs.lammps.org/pair_zbl.html) formalism
- The [timestep](https://docs.lammps.org/timestep.html) was set equal to 0.001 ps
- The simulation is let [run](https://docs.lammps.org/run.html) to relax the box in the [npt ensemble](https://docs.lammps.org/fix_nh.html) (300 K, recalled every 0.1 ps) with 50'000 steps
- three [groups](https://docs.lammps.org/group.html) are set: the PKA group (chosen randomly at each collisional cascade), the boundary group (all the atoms at the 6 surfaces of the box) and the core group (the rest of the particles)
- [nve ensemble](https://docs.lammps.org/fix_nve.html) is applied the core group
- [Nosé-Hoover temperature rescaling](https://docs.lammps.org/fix_temp_rescale.html) at 300 K is applied to the boundary group 
- An [adaptive timestep](https://docs.lammps.org/fix_dt_reset.html) is set with tmax, tmin and xmax of 1e-6 ps, 1e-3 ps and 0.05 A, respectively
- The PKA is chosen randomly, it is given a random direction and a [velocity] according to its mass and given kinetic energy (i.e. 5 keV)
- The simulation [run](https://docs.lammps.org/run.html) for 50'000 steps
- The simulation outputs a [data](https://docs.lammps.org/write_data.html) or a [read/restart](https://docs.lammps.org/read_restart.html) file with the box configuration at the latest timestep that is used as input file for the following collisional cascade of the loop

The simulation provides about 1500 [dump](https://docs.lammps.org/dump.html) files that can be visualised and postprocessed with [Ovito](https://www.ovito.org/about/).
