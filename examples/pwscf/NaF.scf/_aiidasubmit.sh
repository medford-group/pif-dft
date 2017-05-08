#!/bin/bash

#SBATCH --no-requeue
#SBATCH --job-name="aiida-45483"
#SBATCH --get-user-env
#SBATCH --output=_scheduler-stdout.txt
#SBATCH --error=_scheduler-stderr.txt
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=8
#SBATCH --time=00:30:00


'srun' '-n' '24' '/users/mounet/Quantum_Espresso_v5.2/espresso-5.2.0_daint_rhoxml/bin/pw.x' '-nk' '24' '-in' 'aiida.in'  > 'aiida.out' 
