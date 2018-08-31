#!/bin/bash

lammps="../../lmp_mpi" # Path to your LAMMPS executable.

N=2  # Number of cores in the simulation.

# Create directory structure for data output.
mkdir -p data

mpirun -n N ${lammps} -in input_nehi.lmp -var RANDOM ${RANDOM}

