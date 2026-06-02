#!/bin/bash
#SBATCH --job-name=scaling
#SBATCH --partition=nodes
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --output=scaling_%j.log

# Load the module
module purge
module load foss

# Compile code
mpif90 -O3 -fopenmp -o coax_hybrid coax_hybrid.f90

ranks=(1 2 4 8 12)
threads=(1 2 4 8)

for r in "${ranks[@]}"; do
	for t in "${threads[@]}"; do

		echo "MPI Ranks: $r"
		echo "OpenMP Threads: $t"

		export OMP_NUM_THREADS=$t
		export OMP_PLACES=cores
		export OMP_PROC_BIND=close

		time mpirun -np $r \
		--map-by socket:PE=$t \
		 ./coax_hybrid
	done
done