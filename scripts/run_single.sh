#!/bin/bash
#SBATCH --job-name=r16_t6
#SBATCH --partition=nodes
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=6
#SBATCH --time=00:30:00
#SBATCH --output=scaling_r16_t6_%j.log

# Load the module
module purge
module load foss

# Compile code
mpif90 -O3 -fopenmp -o coax_hybrid coax_hybrid.f90

r=16
t=6

echo "MPI Ranks: $r"
echo "OpenMP Threads: $t"

export OMP_NUM_THREADS=$t
export OMP_PLACES=cores
export OMP_PROC_BIND=close

mpirun -np $r \
--map-by socket:PE=$t \
./coax_hybrid
