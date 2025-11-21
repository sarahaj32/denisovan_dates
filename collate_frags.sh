#!/bin/bash
#SBATCH --job-name=collate_frags
#SBATCH --account=co_moorjani
#SBATCH --partition=savio4_htc
#SBATCH --ntasks-per-node=1
#SBATCH --time=3:00:00
#SBATCH --output=%x.%j.out

#path="/global/scratch/p2p3/pl1_moorjani/lauritsskov2/ArchaicSegments/Results/03decode/LASIDAD"
path=$1
dataset=$2

echo $path
echo $dataset

mkdir -p "data/$dataset/"

# set output file names
output="data/$dataset/${dataset}_fragments.txt"
output_08="data/$dataset/${dataset}_fragments_08.txt"
output_08_snps10="data/$dataset/${dataset}_fragments_snps10_08.txt"

# write out the header
file=$(ls $path | head -n 1)
head -n 1 "$path/${file}" | awk 'BEGIN {OFS="\t"} {print "name", "haplotype", "pop", "region", $1, $2, $3, $6, $5, $7, $8, $9, $10, $11, $12, $13}'> $output

# Then, for all files, skip their headers and append the rest
for f in ${path}/*.hap*; 
do
    echo $f
    # extract sample ID
    id=$(basename "$f" | sed 's/.hap[12].txt//')
    hap=$(basename "$f" | sed -E 's/.*\.(hap[0-9]+)\.txt/\1/')
    echo $hap
    echo $id $(basename "$f" | sed -E 's/.*\.(hap[0-9]+)\.txt/\1/')
    grep 'Archaic' "$f" | tail -n +2 | awk -v id="$id" -v hap="$hap" 'BEGIN {OFS="\t"} {print id, hap, "pop", "region", $1, $2, $3, $6, $5, $7, $8, $9, $10, $11, $12, $13}' >> $output
done

awk '$1 == "name" || ($8>=0.8)' $output > $output_08
awk '$1 == "name" || ($10 >= 10)' $output_08 > $output_08_snps10