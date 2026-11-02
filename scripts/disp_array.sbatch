#!/bin/bash
#SBATCH --job-name=YB4t_ph_disp
#SBATCH --ntasks=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --array=1-__N__%4
#SBATCH --output=disp-%a/vasp_%A.log
#SBATCH --error=disp-%a/vasp_%A.err
export OMP_NUM_THREADS=1
module purge
# load your site VASP modules here, e.g. module load intel openmpi hdf5
cd disp-$(printf "%03d" $SLURM_ARRAY_TASK_ID)
srun --mpi=pmix --cpu-bind=none "${VASP_BIN:-vasp_std}"
