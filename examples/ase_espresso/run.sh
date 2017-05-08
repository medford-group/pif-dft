#!/bin/bash
#PBS -l procs=16
#PBS -l walltime=24:00:00
#PBS -q joe-6-intel
#PBS -N N2
#PBS -o stdout
#PBS -e stderr
cd $PBS_O_WORKDIR

source  /gpfs/pace1/project/chbe-medford/medford-share/envs/espresso-5.1.r11289-pybeef

python qn_opt.py #submit using mpirun (parallel computing) to 12 cores. The 12 here must match the ppn specified in the header.
