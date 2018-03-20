#!/bin/bash
#
#----------------------------------
# single GPU + single CPU example
#----------------------------------
#
#SBATCH --job-name=R2_calculations
#SBATCH --output=r2_log_hist
#
#number of CPUs to be used
#SBATCH --ntasks=1
#
#Define the number of hours the job should run. 
#Maximum runtime is limited to 10 days, ie. 240 hours
#SBATCH --time=36:00:00
#
#Define the amount of system RAM used by your job in GigaBytes
#SBATCH --mem=10G
#
#Send emails when a job starts, it is finished or it exits
#SBATCH --mail-user=ekaterina.putintseva@ist.ac.at
#SBATCH --mail-type=ALL
#
#Pick whether you prefer requeue or not. If you use the --requeue
#option, the requeued job script will start from the beginning, 
#potentially overwriting your previous progress, so be careful.
#For some people the --requeue option might be desired if their
#application will continue from the last state.
#Do not requeue the job in the case it fails.
#SBATCH --no-requeue
#
#Define the "gpu" partition for GPU-accelerated jobs
#SBATCH --partition=gpu
#
#Define the number of GPUs used by your job
#SBATCH --gres=gpu:1
#
#Define the GPU architecture (GTX980 in the example, other options are GTX1080Ti, K40)
##SBATCH --constraint=GTX980
#
#Do not export the local environment to the compute nodes
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV
#
#for single-CPU jobs make sure that they use a single thread
export OMP_NUM_THREADS=1
#
#load an CUDA software module
module load tensorflow/python-2.7/1.3.0 
#
#print out the list of GPUs before the job is started
/usr/bin/nvidia-smi
#
#run your CUDA binary through SLURM's srun
sleep 5
echo $pwd
echo $HOSTNAME

ls -l ./network.py
srun --cpu_bind=verbose python ./network.py -c S7 -n 10

