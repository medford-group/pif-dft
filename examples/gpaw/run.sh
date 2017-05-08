#!/bin/bash
#PBS -l nodes=2:ppn=12
#PBS -l walltime=12:00:00
#PBS -q joe-6-ge
#PBS -N optimizer
#PBS -o stdout
#PBS -e stderr
cd $PBS_O_WORKDIR

module purge
module load intel/14.0.2
module load openmpi/1.8
module load libxc/2.2.2
module load mkl/11.2
module load fftw/3.3.4
module load python/2.7

export PYTHONHOME=/usr/local/pacerepov1/python/2.7/intel-14.0.2
export PATH=/nv/hp22/amedford6/medford-shared/INSTALL/bin:$PATH
export PYTHONPATH=/nv/hp22/amedford6/medford-shared/INSTALL/lib/python2.7/site-packages:$PYTHONPATH
export LD_LIBRARY_PATH=/nv/hp22/amedford6/medford-shared/INSTALL/lib/python2.7/site-packages:$LD_LIBRARY_PATH

mpirun -np 24 gpaw-python qn_opt.py #submit using mpirun (parallel computing) to 12 cores. The 12 here must match the ppn specified in the header.
