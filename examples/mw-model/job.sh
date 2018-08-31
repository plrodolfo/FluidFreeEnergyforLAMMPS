#!/bin/bash

lammps="../../lmp_mpi" # Path to your LAMMPS executable.

N=2  # Number of cores in the simulation.

# Create directory structure for data output.
mkdir -p data
mkdir -p data2

mpirun -n N ${lammps} -in input_nehi.lmp -var RANDOM ${RANDOM}

mpirun -n N ${lammps} -in input_ners.lmp -var RANDOM ${RANDOM}
