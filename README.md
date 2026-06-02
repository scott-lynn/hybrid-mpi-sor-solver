# Hybrid MPI+OpenMP SOR Solver for Coaxial Potentials

This repository contains a high-performance, hybrid-parallelized (MPI + OpenMP) Successive Over-Relaxation (SOR) solver, developed in modern Fortran. The codebase was designed to perform rigorous strong and weak scaling analysis on AMD EPYC Milan architectures, scaling execution to 96 processing elements.

This project was developed as part of the graduate Scientific Supercomputing module at the University of York.

## Key Architectural Features

This solver is built around three core systems engineering principles to maximise parallel efficiency and bypass hardware latency cliffs:

### 1. Asynchronous Halo Exchange (MPI)
To manage domain decomposition across distributed memory, the solver utilizes non-blocking MPI communication (`MPI_Isend` / `MPI_Irecv`). This architecture allows the solver to actively overlap network communication latency with bulk computation, ensuring the CPU remains saturated during cross-node data transfers.

### 2. Red-Black Checkerboard Ordering (OpenMP)
Standard SOR algorithms contain inherent data dependencies that prevent thread-level parallelization. To solve this, a red-black spatial ordering scheme was implemented. By decoupling the grid into independent sub-domains, the inner computational loops are safely parallelized across shared-memory threads using OpenMP without race conditions.

### 3. Hardware & NUMA-Aware Execution
To optimize for the specific chiplet architecture of AMD EPYC processors, the execution scripts explicitly manage thread affinity and memory locality. By utilizing `OMP_PROC_BIND=close`, `OMP_PLACES=cores`, and explicit socket mapping (`--map-by socket:PE=$t`), the execution environment prevents thread migration across Core Complex Die (CCD) boundaries, minimizing L3 cache misses and cross-socket latency.

## Repository Structure

* `src/coax_hybrid.f90`: The core hybrid solver source code.
* `scripts/`: SLURM batch scripts used for executing strong/weak scaling tests on the Viking supercomputer.

## Compilation & Execution

The code is designed to be compiled with the GNU Fortran compiler, leveraging `-O3` for automatic vectorisation and loop unrolling.

```bash
# Compile the hybrid solver
mpif90 -O3 -fopenmp -o coax_hybrid src/coax_hybrid.f90

# Example execution (4 MPI ranks, 8 OpenMP threads per rank)
export OMP_NUM_THREADS=8
export OMP_PLACES=cores
export OMP_PROC_BIND=close

mpirun -np 4 --map-by socket:PE=8 ./coax_hybrid
```
