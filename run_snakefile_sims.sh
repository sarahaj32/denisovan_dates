#!/bin/bash

#SBATCH --account=co_moorjani
#SBATCH --partition=savio3_htc
#SBATCH --ntasks-per-node=1
#SBATCH --time=10:00:00

sim=$1
seed=$2
nMH=$3
frag=$4
pos=$5

id=$1":"$2

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
mkdir -p logs/$id/

echo $id
echo "log files will be sent to logs/$id/"

# unlock  snakemake
#snakemake --unlock --config seed=$seed sim=$sim 
#echo "directory unlocked"

echo "script beginning"
snakemake --rerun-incomplete --use-conda --keep-going -p --jobs 100 \
    --cluster-config $HOME/.config/snakemake/slurm.co_moorjani/settings.json \
    --default-resources "cpus=1" "time='3:00:00'" \
    --cluster "sbatch -A {cluster.account} -p 'savio3_htc' -n {cluster.ntasks} -c {resources.cpus} -t {resources.time} -J {rule} -o logs/$id/{rule}.{jobid}.{wildcards}.out -e logs/$id/{rule}.{jobid}.{wildcards}.err" \
    --config seed=$seed sim=$sim nMH=$nMH frag=$frag pos=$pos
echo "DONE"

