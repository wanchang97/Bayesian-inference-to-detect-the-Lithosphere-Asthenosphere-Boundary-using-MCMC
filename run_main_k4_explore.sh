#!/bin/bash
#SBATCH --job-name=Matlab_main_k4_explore_Job
##Select the desired partiiton:
#SBATCH --partition=standard
##Define the maximum number of the cores you need to use:
##Keep in mind that some matlab functions use internal multithreading.
##Do some tests which different values on the number of cores to check if your application run faster.
#SBATCH -c 3
##Define the maximum memory you will use:
#SBATCH --mem=4G
##Define the maximum time your job will run.
##With format Days-hours
#SBATCH --time=5-00
echo "Start:" 'date'
echo "I ran on:"
cd $SLURM_SUBMIT_DIR
echo $SLURM_NODELIST

## load the desired matlab version
module load matlab/R2020a
##The matlab script should be called without the .m extension
matlab -nojvm -nodisplay -r "main_k4_explore;quit"
echo "End:" 'date'
