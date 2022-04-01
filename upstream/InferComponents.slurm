#!/bin/sh

#SBATCH --mem=16GB
#SBATCH --time-min=6000
#SBATCH --nice=500
#SBATCH --array=0-18

runpath=/home/lad44/davidlab/users/fsm/bacteroides_20220307/code/analysis/process_fcs_samples/
code=$runpath/InferComponents.sh

txt_out=$runpath/slurm_log/InferComponents_$SLURM_ARRAY_TASK_ID.out
txt_err=$runpath/slurm_log/InferComponents_$SLURM_ARRAY_TASK_ID.err

srun -o $txt_out -e $txt_err sh $code $SLURM_ARRAY_TASK_ID 
