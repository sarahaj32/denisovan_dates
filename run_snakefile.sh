#!/bin/bash

#SBATCH --account=co_moorjani
#SBATCH --partition=savio4_htc
#SBATCH --ntasks-per-node=1
#SBATCH --time=30:00:00

export TMPDIR=$SCRATCH/tmp
mkdir -p $TMPDIR
echo "TMPDIR is $TMPDIR"

echo "script beginning"
source /global/scratch/users/sarahj32/software/miniconda3/bin/activate /global/scratch/users/sarahj32/software/miniconda3/envs/snakemake/

echo "conda set"
echo "starting"
module load bio/bcftools

echo "module loaded"

# make a directory to save log files
mkdir -p logs/

#unlock  snakemake
snakemake --unlock 
echo "directory unlocked"

echo "script beginning"
snakemake --rerun-incomplete --use-conda --keep-going -p --jobs 100 \
    --cluster-config $HOME/.config/snakemake/slurm.co_moorjani/settings.json \
    --default-resources "cpus=1" "time='3:00:00'" \
    --cluster "sbatch -A {cluster.account} -p 'savio4_htc' -n {cluster.ntasks} -c {resources.cpus} -t {resources.time} -J {rule} -o logs/{rule}.{jobid}.{wildcards}.out -e logs/{rule}.{jobid}.{wildcards}.err"
echo "DONE"

