#!/bin/bash

#SBATCH --account=co_moorjani
#SBATCH --partition=savio3
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --time=36:00:00

sim=$1
seed=$2
nMH=$3
frag=$4
pos=$5

echo "sim is $sim"
echo "seed is $seed"
echo "nMH is $nMH"
echo "frag is $frag"
echo "pos is $pos"

export TMPDIR=/tmp/$USER
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

# unlock  snakemake
#snakemake --unlock --config seed=$seed sim=$sim nMH=$nMH frag=$frag pos=$pos
#echo "directory unlocked"
#mkdir -p /global/scratch/p2p3/pl1_moorjani/sarahj32/denisova_dating/.snakemake/locks

echo "script beginning"
snakemake --rerun-incomplete --keep-going --use-conda -p -c all \
    --config seed=$seed sim=$sim nMH=$nMH frag=$frag pos=$pos
echo "DONE"


