#!/bin/bash

#SBATCH --mail-user=aritram@ntnu.no
#SBATCH --mail-type=END
#SBATCH --mail-type=BEGIN
##SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=TIME_LIMIT_90

#SBATCH --account=share-iv-ept
##SBATCH --account=iv-ept

#SBATCH --partition=CPUQ

#SBATCH --job-name=test
#SBATCH --output=log.out
#SBATCH --time=03-00:00:00
#SBATCH --mem=10GB
#SBATCH --ntasks-per-node=40
#SBATCH --nodes=1

##SBATCH --test-only

# https://slurm.schedmd.com/sbatch.html

module purge
module load foss/2022b
#module load OpenMPI/4.0.5-GCC-10.2.0
#module load intel/2020b
#module load FFTW/3.3.8-intel-2020b
module list

cd ${SLURM_SUBMIT_DIR}

#module load OpenMPI/4.0.5-GCC-10.2.0
mpirun -n 40 -use-hwthread-cpus flutas> out.log 2>&1
## gfortran
##srun --mpi=pmix ./flutas
## intel
## srun ./sc_compiled/flow36
