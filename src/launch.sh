#!/bin/bash   
#SBATCH --account=NN9561K 
#SBATCH --job-name=ep2.0
#################SBATCH --qos=devel
#SBATCH --time=00-71:59:58
#SBATCH --nodes=5                                                                                                                                                                                                             
#SBATCH --ntasks-per-node=10
#SBATCH -n 50
#SBATCH -e log_%J.err
#SBATCH -o log_%J.out
#SBATCH -D .
#module load intel/2019b FFTW/3.3.8-intel-2019b
#module load OpenMPI/4.0.3-GCC-9.3.0
#module load OpenMPI/4.0.5-GCC-10.2.0

ulimit -s unlimited

module load OpenMPI/4.0.5-GCC-10.2.0
mpirun -n 50 -use-hwthread-cpus flutas> out.log 2>&1
