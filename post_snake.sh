#!/bin/bash

# unlock everything
snakemake --unlock --config seed=0 sim="" nMH=0 frag="" pos=0

# first we just generate the simulations: (this generates simulations, finds true ancestry fragments, runs hmmix and ibdmix, )
for sim in  "denisovanNeanderthal_closeDiv"; #"denisovanNeanderthal_Dlow"   "denisovanNeanderthal_closeDiv"
do
    for seed in 1 2 3;
    do
        sbatch --job-name="post_${seed}:${sim}" --output=slurm-post_${sim}_${seed}_ds200k_%j.out run_snakefile_sims_bigmem.sh $sim $seed 
    done
done

# then we do all of the processing:
# all rerunning!!!
# need to do after preemption: 
slurm-post_denisovanNeanderthal_simple_nMH20_1_hmmix.0.8.200Afr.fullArchaic.1NeaMostLikely_ds200k_
slurm-post_denisovanNeanderthal_simple_nMH20_1_hmmix.0.8.200Afr.fullArchaic.mostLikely_ds200k_

# need to do after unlocking:
slurm-post_denisovanNeanderthal_simple_nMH100_3_hmmix.0.8.200Afr.fullArchaic.mostLikely.50k_ds200k_29960354.out
slurm-post_denisovanNeanderthal_simple_nMH100_3_hmmix.0.7.200Afr.fullArchaic.mostLikely_ds200k_29960351.out
slurm-post_denisovanNeanderthal_simple_nMH100_3_hmmix.0.9.200Afr.fullArchaic.mostLikely_ds200k_29960349.out

# run ds200k for most of the versions 
# here we can add on to do additional filters - I CAN DO THIS!!
for sim in "denisovanNeanderthal_simple";
do
    for seed in 1 2 3;
    do
        for nMH in 20 100;
        do
            for frag in "hmmix.0.8.200Afr.fullArchaic.1NeaMostLikely" "hmmix.0.8.200Afr.fullArchaic.mostLikely" "hmmix.0.9.200Afr.fullArchaic.mostLikely" "hmmix.0.7.200Afr.fullArchaic.mostLikely" "hmmix.0.8.200Afr.fullArchaic.mostLikely.50k" ; # "simIntro" "ibdmix" 
            do
                sbatch --job-name="post_${seed}:${sim}_nMH${nMH}_frag$frag" --output=slurm-post_${sim}_nMH${nMH}_${seed}_${frag}_ds200k_%j.out run_snakefile_sims_bigmem.sh $sim $seed $nMH $frag "ds200k"
            done
        done
    done
done


# for the close divergence time, lets just do classic hmmix and ibdmix
# I still want both 20 and 100 NAMH though
for seed in 1 2 3;
do
    for nMH in 20 100;
    do
        for frag in "hmmix.0.8.200Afr.fullArchaic.mostLikely" "ibdmix"; 
        do
            sbatch --job-name="Dclose${seed}:frag${frag}_nMH${nMH}" --output=slurm-post_denisovanNeanderthal_closeDiv_nMH20_${seed}_${frag}_ds200k_%j.out run_snakefile_sims_bigmem.sh "denisovanNeanderthal_closeDiv" $seed $nMH $frag "ds200k"
        done
    done
done

# only use "all" positions for 1 repetition of denisovanNeanderthal_simple with 20 MH, and normal hmmix/ibdmix
# since this will by far take the longest
for sim in "denisovanNeanderthal_closeDiv" "denisovanNeanderthal_simple";
do
    for nMH in 20 100;
    do
        for frag in "hmmix.0.8.200Afr.fullArchaic.mostLikely" "ibdmix";
        do
            sbatch --job-name="Dall:${sim}_nMH${nMH}_frag${frag}" --output=slurm-post_${sim}_nMH20_1_${frag}_all_%j.out run_snakefile_sims.sh $sim 1 $nMH $frag "all" 
        done
    done
done

# repeat the closest demography I have to what I did for CCB seminar - make sure that still works
### NEED TO REPEAT ALL OF THIS WITH THE CORRECTED ARCHAICS FILES (BUT MAYBE NOT NECESSARY??)
sbatch --job-name="post_${seed}:${sim}" --output=slurm-post_${sim}_${seed}_ds200k_%j.out run_snakefile_sims.sh denisovanNeanderthal_simple_previous 1 

sim="denisovanNeanderthal_simple_previous"
for nMH in 20 100;
do
    for frag in "hmmix.0.8.200Afr.fullArchaic.mostLikely" "ibdmix";
    do
        sbatch --job-name="Dall:${sim}_nMH${nMH}_frag${frag}" --output=slurm-post_${sim}_nMH${nMH}_1_${frag}_all_%j.out run_snakefile_sims.sh ${sim} 1 $nMH $frag "ds200k" 
    done
done

# want all archaic classification performance metrics
# 1 nea vs 3 
