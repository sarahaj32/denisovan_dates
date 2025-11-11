wildcard_constraints:
    rep="\d+",
    nafr="\d+",
    method="(sim_intro|hmmix.100|hmmix.500|hmmix.200)",

##### CHECK_mutSamps #####
# can now add in laurits' exact simulations
rule sims_all:
    input:
        expand("simdat/{sim}_{seed}/{archaic}/fragments/{sim}_sim_intro_{rep}.bed", rep = list(range(1,21)), sim = ["denisovanNeanderthal_simple", "denisovanNeanderthal_Dlow"], seed = [1,2,3], archaic = {"DEN", "NEA"}),
        expand("simdat/{sim}_{seed}/archaic_0.8/{archaic}/fragments/{sim}_hmmix.200_mostLikely_{rep}.bed", rep = list(range(1,21)), sim = ["denisovanNeanderthal_simple", "denisovanNeanderthal_Dlow"], seed = [1,2,3], archaic = {"DEN", "NEA"}),
        expand("simdat/{sim}_{seed}/archaic_0.8/{archaic}/fragments/{sim}_hmmix.200_ND_{rep}.bed", rep = list(range(1,21)), sim = ["denisovanNeanderthal_simple", "denisovanNeanderthal_Dlow"], seed = [1,2,3], archaic = {"DEN", "NEA"}),
        expand("simdat/{sim}_{seed}/ibdmix/fragments/{sim}_ibdmix_{archaic}_{rep}.txt", rep = list(range(1,21)), sim = ["denisovanNeanderthal_simple", "denisovanNeanderthal_Dlow"], seed = [1,2,3], archaic = {"DEN", "NEA"})

    #     expand("simdat/{sim}_{seed}/DEN/curve/results/{sim}_{frag_source}_{nafr}Afr_ingroup_{rep}_D.txt", rep = list(range(1,21)), sim = ["denisovan_simple", "denisovan_early", "denisovan_late", "denisovan_low"], seed = [1,2,3], nafr = [200], frag_source = ["sim_intro", "hmmix.200"]),
    #     #expand("simdat/{sim}_{seed}/DEN/curve/results/{sim}_hmmix.200_{nafr}Afr_ingroup_{rep}_D.txt", rep = list(range(1,21)), sim = ["denisovan_simple", "denisovan_early", "denisovan_late", "denisovan_low"], seed = [1,2,3], nafr = [100, 200]),
    #     expand("simdat/{sim}_{seed}/DEN/curve/results/{sim}_{frag_source}_{rep}_D.txt", rep = list(range(1,21)), sim = ["denisovan_simple", "denisovan_early", "denisovan_late", "denisovan_low"], frag_source = ["sim_intro", "hmmix.200"], seed = [1,2,3]),
    #     expand("simdat/{sim}_{seed}/DEN/fragments/{sim}_sim_intro_{rep}.bed", rep = list(range(1,21)), sim = ["denisovan_LS_ice", "denisovan_LS_pap"], seed = [1]), # increase to 500???? , "hmmix_LS" , "denisovan_simpleLS" "denisovan_smallAfr", "denisovan_bigNonAfr",
    #     expand("simdat/{sim}_{seed}/DEN/fragments/{sim}_hmmix.{nafr}_{rep}.bed", rep = list(range(1,21)), sim = ["denisovan_LS_ice", "denisovan_LS_pap", "denisovan_simple"], nafr = [100, 200], seed = [1]), # increase to 500???? , "hmmix_LS" , "denisovan_simpleLS" "denisovan_smallAfr", "denisovan_bigNonAfr",
        # add on neanderthal introgression, and labelling of fragments
        #expand("simdat/denisovanNeanderthal_simple_{seed}/eig/denisovanNeanderthal_simple_sim_allChr.snp", seed = [1,2,3]),
        # denisovan and Neanderthal annotations and inferences:
        # expand("simdat/{sim}_{seed}/archaic_0.8/fragments/{sim}_hmmix.200_{rep}.bed", rep = list(range(1,21)), sim = ["denisovanNeanderthal_simple"], seed = [1,2,3], method = ["sim_intro"], nafr = [200]),
        # expand("simdat/{sim}_{seed}/{archaic}/curve/results/{sim}_{frag_source}_{rep}_D.txt", rep = list(range(1,21)), seed = [1,2,3], frag_source = ["sim_intro", "hmmix.200_ND", "hmmix.200_mostLikely"], archaic = ["DEN", "NEA"], sim = ["denisovanNeanderthal_simple", "denisovan_simple", "denisovanNeanderthal_Dlow"]),
        # expand("simdat/{sim}_{seed}/{archaic}/curve/results/{sim}_{frag_source}_ds200k_{rep}_D.txt", rep = list(range(1,21)), seed = [1,2,3], frag_source = ["sim_intro", "hmmix.200_ND", "hmmix.200_mostLikely"], archaic = ["DEN", "NEA"], sim = ["denisovanNeanderthal_simple", "denisovan_simple", "denisovanNeanderthal_Dlow"]),

        # expand("simdat/{sim}_{seed}/archaics/curve/results/{sim}_annotatedArchaics_{rep}_D.txt", rep = list(range(1,21)), sim = ["denisovanNeanderthal_simple"], seed = [1,2,3]),
        # expand("simdat/{sim}_{seed}/{archaic}/curve/results/fragments/{sim}_sim_intro_{rep}_D.txt", rep = list(range(1,21)), sim = ["denisovanNeanderthal_simple"], seed = [1,2,3], archaic = {"DEN", "NEA"}),
        # expand("simdat/{sim}_{seed}/archaic_0.8/curve/results/fragments/{sim}_hmmix.200_{rep}_D.txt", rep = list(range(1,21)), sim = ["denisovanNeanderthal_simple"], seed = [1,2,3], method = ["sim_intro"], nafr = [200]),
        # expand("simdat/{sim}_{seed}/{archaic}/curve/results/fragments/{sim}_hmmix.200_mostLikely_{rep}_D.txt", rep = list(range(1,21)), sim = ["denisovanNeanderthal_simple"], seed = [1,2,3], archaic = {"DEN", "NEA"}),

        # expand("simdat/{sim}_{seed}/{archaic}/curve/results/{sim}_{frag_source}_{nafr}Afr_ingroup_{rep}_D.txt", rep = list(range(1,21)), seed = [1,2,3], nafr = [200], frag_source = ["sim_intro", "hmmix.200"], archaic = ["DEN", "NEA"], sim = ["denisovanNeanderthal_simple"]),
        # expand("simdat/{sim}_{seed}/{archaic}/curve/results/{sim}_{frag_source}_{rep}_D.txt", rep = list(range(1,21)), seed = [1,2,3], nafr = [200], frag_source = ["sim_intro", "hmmix.200"], archaic = ["DEN", "NEA"], sim = ["denisovanNeanderthal_simple"])
        
    output:
        "logs/CHECK_all_done.txt"
    shell:
        "echo 'all done' > {output}"

# I can use this as the simulation map:
# /global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/hg19/resources/genetic_maps_corrections/hg19/Shared_AfAmMap/Europeanized_AA_map.txt
rule simulate: # takes ~1 hr to simulate 125Mb
# simulate specified populations using specified demography
# produce a vcf of the populations of interest, as well as saves the tree file from the simulation (so we can access and compute on it later)
    input: 
        deme = lambda wildcards: config["simulations"]["demo_yaml"][wildcards.sim],
        json = config["simulations"]["samples_json"],
        script = "/global/scratch/p2p3/pl1_moorjani/sarahj32/superarchaic_introgression/super_clean/helper_scripts/simulate_from_deme_only.py"
    output: 
        vcf = "simdat/{sim}_{seed}/raw/{sim}_sim_{rep}.vcf",
        tree = "simdat/{sim}_{seed}/raw/{sim}_sim_{rep}.tree",
    conda: 
        "sim.yaml"
    resources:
        # increasing the cpus is the way to increase memory available when submitting with savio_htc
        # these need a lot of memory
        # this job also needs a bit more time
        cpus=3,
        time='5:00:00'
    shell: 
        # I need to update the names of the individuals to match real data so that the downstream code works properly
        """
        python {input.script} -j {input.json} -c {wildcards.rep} -y {input.deme} -o {output.vcf} -g 29 -l 125000000 -m 1.29e-8 -r 1e-8 -t {output.tree} -s {wildcards.seed}
        sed -i 's/NEA_1_1/AltaiNeandertal/g' {output.vcf}
        sed -i 's/NEA_2_1/Vindija33.19/g' {output.vcf}
        sed -i 's/NEA_3_1/Chagyrskaya-Phalanx/g' {output.vcf}
        sed -i 's/DEN_1/Denisova/g' {output.vcf}
        """

# filter to biallelic positions
rule remove_recurring_mutations:
    input: 
        "simdat/{sim}_{seed}/raw/{sim}_sim_{rep}.vcf"
    output:
        "simdat/{sim}_{seed}/clean/{sim}_sim_{rep}.vcf"
    shell:
        "bcftools view -c 1 {input} | bcftools norm -m + | awk '$3 !~ /;/' > {output}" # 


rule concatenate_sim_chroms:
    input: 
        expand("simdat/{{sim}}_{{seed}}/clean/{{sim}}_sim_{chrom}.vcf", chrom = list(range(1,21)))
    output:
        "simdat/{sim}_{seed}/vcf/{sim}_sim_allChr.vcf"
    resources:
        cpus=4
    shell:
        "bcftools concat {input} -o {output}"

# convert to eigenstrat format
rule convert_sims_to_eig:
    input:
        vcf="simdat/{sim}_{seed}/vcf/{sim}_sim_allChr.vcf",
        script="/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/bin/vcf2eigenstrat.py"
    output:
        geno="simdat/{sim}_{seed}/eig/{sim}_sim_allChr.geno",
        ind="simdat/{sim}_{seed}/eig/{sim}_sim_allChr.ind",
        snp="simdat/{sim}_{seed}/eig/{sim}_sim_allChr.snp",
    params:
        path="simdat/{sim}_{seed}/eig/{sim}_sim_allChr"
    resources:
        time='6:00:00',
        cpus=4
    shell:
        """ 
        python2 {input.script} -v {input.vcf} -o {params.path}
        """ 

# we also want these to be non-variant in africans (we'll call these ingroup)
# we only care about the positions - drop the genotypes
rule remove_outgroup_pos:
    input: 
        vcf = "simdat/{sim}_{seed}/MH_variable/{sim}_sim_{rep}.vcf",
        outgroup = "simdat/{sim}_{seed}/hmmix/outgroup/outgroup_{nafr}Afr.txt"
    output:
        "simdat/{sim}_{seed}/MH_variable/{sim}_sim_{rep}_{nafr}Afr_ingroup.vcf",
    shell:
        """
        bcftools view  -T ^{input.outgroup} -G {input.vcf} -o {output}
        """

# extract archaics
rule extract_archaics:
    # add manifesto
    input: 
        "simdat/{sim}_{seed}/clean/{sim}_sim_{rep}.vcf"
    output:
        full = "simdat/{sim}_{seed}/archaics/{sim}_sim_full_{rep}.vcf",
        mask = "simdat/{sim}_{seed}/archaics/{sim}_sim_200k_{rep}.txt",
        vcf = "simdat/{sim}_{seed}/archaics/{sim}_sim_200k_{rep}.vcf",
        gz = "simdat/{sim}_{seed}/archaics/{sim}_sim_200k_{rep}.vcf.gz",
        index = "simdat/{sim}_{seed}/archaics/{sim}_sim_200k_{rep}.vcf.gz.csi",
    shell: # Need to assign nucleotides to the archaic mutations
        """
        bcftools view -c1 -s AltaiNeandertal,Vindija33.19,Chagyrskaya-Phalanx,Denisova {input} | awk 'BEGIN{{OFS="\t"}} /^#/ {{print; next}} {{ $4="A"; $5="T"; print }}' > {output.full}
        bcftools view -H {output.full} | cut -f1,2 | shuf -n 200000 | sort -nk2 > {output.mask}
        bcftools view {output.full} -T {output.mask} -o {output.vcf}
        bcftools view {output.vcf} -Ob -o {output.gz}
        bcftools index {output.gz}
        """

# extract archaics
rule get_archaic_vcfs:
    # add manifesto
    input: 
        "simdat/{sim}_{seed}/clean/{sim}_sim_{rep}.vcf"
    output:
        den = "simdat/{sim}_{seed}/archaics/{sim}_sim_DEN_{rep}.vcf",
        nea = "simdat/{sim}_{seed}/archaics/{sim}_sim_NEA_{rep}.vcf",
    shell: # Need to assign nucleotides to the archaic mutations
        """
        bcftools view -c1 -s Vindija33.19 {input} | awk 'BEGIN{{OFS="\t"}} /^#/ {{print; next}} {{ $4="A"; $5="T"; print }}' > {output.den}
        bcftools view -c1 -s Denisova {input} | awk 'BEGIN{{OFS="\t"}} /^#/ {{print; next}} {{ $4="A"; $5="T"; print }}' > {output.nea}
        """

# extract modern humans
rule get_MH_vcfs:
    # add manifesto
    input: 
        "simdat/{sim}_{seed}/clean/{sim}_sim_{rep}.vcf"
    output:
        "simdat/{sim}_{seed}/MH20/{sim}_sim_MH20_{rep}.vcf",
    shell: 
        """
        samps=$(printf 'NAMH_%d,' {{1..20}} | sed 's/,$//')
        bcftools view -s $samps {input} | awk 'BEGIN{{OFS="\t"}} /^#/ {{print; next}} {{ $4="A"; $5="T"; print }}'> {output}
        """

rule ibdmix_generate_genotypes:
    input:
        MH="simdat/{sim}_{seed}/MH20/{sim}_sim_MH20_{rep}.vcf",
        archaic="simdat/{sim}_{seed}/archaics/{sim}_sim_{archaic}_{rep}.vcf"
    output:
        "simdat/{sim}_{seed}/ibdmix/gts/{sim}_sim_{archaic}_{rep}.txt"
    shell:
        "generate_gt -a {input.archaic} -m {input.MH} -o {output}"

rule run_ibdmix:
    input:
        "simdat/{sim}_{seed}/ibdmix/gts/{sim}_sim_{archaic}_{rep}.txt"
    output:
        "simdat/{sim}_{seed}/ibdmix/fragments/{sim}_ibdmix_{archaic}_{rep}.txt"
#"simdat/denisovanNeanderthal_simple_3/ibdmix/fragments/denisovanNeanderthal_simple_ibdmix_DEN_3.txt"
    shell:
        "ibdmix -g {input} -o {output} --write-lods --write-snps"

rule annotate_archaic_variants:
    input:
        archaics="simdat/{sim}_{seed}/archaics/{sim}_sim_200k_{rep}.vcf",
        script="add_annotations_to_archaic_observations.py"
    output:
        annotated="simdat/{sim}_{seed}/archaics/{sim}_200k_annotatedArchaics_{rep}.txt"
    shell:
        "python {input.script} -i {input.archaics} -o {output.annotated}"


# we want to get a more reasonable number of positions, equivalent to the strict mask (there are variants 61,000,000 in 1000G strict mask)
# rule downsample_to_strict_mask:
#     input:
#         "simdat/{sim}_{seed}/clean/{sim}_sim_{rep}.vcf"
#     output:
#         vcf = "simdat/{sim}_{seed}/masked/{sim}_sim_{rep}.vcf"
#         pos = temp("simdat/{sim}_{seed}/masked/{sim}_sim_{rep}_maskPos.txt")
#     shell:
#         """
#         bcftools view -H {input} | cut -f1,2 | shuf -n 6750 | sort -nk2 > {output.pos}
#         bcftools view -T {output.pos} {input} -o {output.vcf}
#         """

# we only care about positions variable within our populations of interest (and individuals of interest)      
rule remove_nonvariable_sites:
    input: 
        "simdat/{sim}_{seed}/clean/{sim}_sim_{rep}.vcf"
    output:
        "simdat/{sim}_{seed}/MH_variable/{sim}_sim_{rep}.vcf"
    shell:
        """
        samps=$(printf 'NAMH_%d,' {{1..20}} | sed 's/,$//')
        bcftools view -s $samps -c 1 {input} | awk 'BEGIN{{OFS="\t"}} /^#/ {{print; next}} {{ $4="A"; $5="T"; print }}'> {output}
        """

# A rule to downsample to 220,000 SNPs (half of what is simulated)
rule downsample:
    input:
        vcf = "simdat/{sim}_{seed}/MH_variable/{sim}_sim_{rep}.vcf",
    output:
        "simdat/{sim}_{seed}/MH_variable/{sim}_sim_{rep}_downsampled.vcf"
    shell:  
        """
        set -x
        grep -v ^# {input.vcf} | shuf -n 200000 | sort -nk2 > {output}
        """

# find simulated ND10 and ND01 positions:

rule get_NDxx_pos:
    input:
        vcf = "simdat/{sim}_{seed}/clean/{sim}_sim_{rep}.vcf",
        afr_der = "simdat/{sim}_{seed}/hmmix/outgroup/outgroup_{nafr}Afr.txt",
        script = "helper_scripts/find_NDxx_pos.py"
    output:
        vcf = "simdat/{sim}_{seed}/archaicsAfrAnc/{sim}_sim_{rep}_{nafr}AfrAnc_VinDen.vcf",
        labeled = "simdat/{sim}_{seed}/archaicsAfrAnc/{sim}_sim_{rep}_{nafr}AfrAnc_VinDenlabeled.txt",
        nd01 = "simdat/{sim}_{seed}/archaicsAfrAnc/{sim}_sim_{rep}_{nafr}AfrAnc_ND01.txt",
        nd10 = "simdat/{sim}_{seed}/archaicsAfrAnc/{sim}_sim_{rep}_{nafr}AfrAnc_ND10.txt",
    shell:
        """
        bcftools view -c1 -s Vindija33.19,Denisova -T ^{input.afr_der} {input.vcf} -o {output.vcf}
        python {input.script} -i {output.vcf} -n Vindija33.19 -d Denisova -r -o {output.labeled}
        grep -i nd01 {output.labeled} > {output.nd01}
        grep -i nd10 {output.labeled} > {output.nd10}
        """

# find the introgressed fragments from Denisovans into humans, using the simulated tree
# takes at least 30-45 minutes, depending on length and populations
rule find_intro_frags:
    input:
        tree = "simdat/{sim}_{seed}/raw/{sim}_sim_{rep}.tree",
        script="/global/scratch/p2p3/pl1_moorjani/sarahj32/superarchaic_introgression/super_clean/helper_scripts/find_intro_frags.py"
    output:
        bed = "simdat/{sim}_{seed}/{archaic}/fragments/{sim}_sim_intro_{rep}.bed"
    conda: 
        "sim.yaml"
    resources:
        time='6:00:00'
    shell:
        """
        python {input.script} -f {input.tree} -b {output.bed} -s {wildcards.archaic} -t OOA -c {wildcards.rep} -n 20
        """

#### A set of HMMix rules to call introgressed fragments ####


# I create the outgroup files myself so that I can control that 0 is always ancestral and 1 is always the derived allele
rule create_outgroup:
    input:
        vcf = expand("simdat/{{sim}}_{{seed}}/clean/{{sim}}_sim_{rep}.vcf", rep = list(range(1,21))),
    output:
        unsorted = temp("simdat/{sim}_{seed}/hmmix/outgroup/outgroup_{nafr}Afr.unsorted"),
        sorted_out = "simdat/{sim}_{seed}/hmmix/outgroup/outgroup_{nafr}Afr.txt"
    shell:
        """
        rm -f {output.unsorted} {output.sorted_out}
        samps=$(printf 'AFR_%d,' {{1..{wildcards.nafr}}} | sed 's/,$//')
        for file in {input.vcf}; 
        do
            bcftools view -s $samps $file -v snps -c 1 | bcftools query -f '%CHROM\t%POS\n' >> {output.unsorted}
        done
        sort -k1,1 -k2,2n {output.unsorted} | uniq >> {output.sorted_out}
        """

rule create_ingroup:
    input:
        # need to update the VCFs so that they have an actual reference and alternative allele
        vcfs = expand("simdat/{{sim}}_{{seed}}/MH_variable/{{sim}}_sim_{rep}.vcf", rep = list(range(1,21))),
        outgroup = "simdat/{sim}_{seed}/hmmix/outgroup/outgroup_{nafr}Afr.txt"
    output:
        "simdat/{sim}_{seed}/hmmix/obs/obs.{nafr}.NAMH_{ind}.txt"
    params:
        ingroup_path = "simdat/{sim}_{seed}/hmmix/obs/obs.{nafr}",
        vcfs = "simdat/{sim}_{seed}/MH_variable/{sim}_sim_*.vcf",
        name = "NAMH_{ind}"
    shell:
        """
        rm -f {output}
        for file in {input.vcfs}; 
        do
            echo $file
            bcftools view -s {params.name} -v snps -c 1 -T ^{input.outgroup} $file | bcftools query -f '%CHROM\t%POS\t%REF\t[%GT]\n' | awk 'BEGIN {{OFS="\t"}} {{gsub("0", "A", $4); gsub("1", "T", $4); gsub("\\\\|", "", $4); print $1,$2,$3,$4}}'  >> {output}
        done
        """

rule train:
    input:
        obs = "simdat/{sim}_{seed}/hmmix/obs/obs.{nafr}.NAMH_{ind}.txt"
    output:
        "simdat/{sim}_{seed}/hmmix/trained/trained.{nafr}.NAMH_{ind}.json"      
    shell:
        "hmmix train -obs={input.obs} -out={output} -haploid"

rule decode:
    input:
        obs = "simdat/{sim}_{seed}/hmmix/obs/obs.{nafr}.NAMH_{ind}.txt",
        trained = "simdat/{sim}_{seed}/hmmix/trained/trained.{nafr}.NAMH_{ind}.json",
        admixpop = expand("simdat/{{sim}}_{{seed}}/archaics/{{sim}}_sim_200k_{rep}.vcf.gz", rep = list(range(1,21)))
    output:
        expand("simdat/{{sim}}_{{seed}}/hmmix/decoded/decoded.{{nafr}}.NAMH_{{ind}}.hap{hap}.txt", hap = [1,2])
    params:
        out_path = "simdat/{sim}_{seed}/hmmix/decoded/decoded.{nafr}.NAMH_{ind}",
        admixpop = "simdat/{sim}_{seed}/archaics/{sim}_sim_200k_*.vcf.gz"
    shell:
        "hmmix decode -obs={input.obs} -param={input.trained} -haploid -out={params.out_path} -admixpop={params.admixpop} -extrainfo"

# combine the fragmentns from each individual, then split by chromosome
# gives us a file of fragments/chromosome
rule collate_frags:
    input:
        expand("simdat/{{sim}}_{{seed}}/hmmix/decoded/decoded.{{nafr}}.NAMH_{ind}.hap{hap}.txt", ind = list(range(1,21)), hap = [1,2])
    output:
        tmp = temp("simdat/{sim}_{seed}/archaic_{threshold}/fragments/{sim}_hmmix.{nafr}_all.bed"),
        chrs = expand("simdat/{{sim}}_{{seed}}/archaic_{{threshold}}/fragments/{{sim}}_hmmix.{{nafr}}_{rep}.bed", rep = list(range(1,21)))
    shell:
        """
        rm -f {output.tmp}
        hap_count=0
        for f in {input};
        do
        echo $f
        id=$(basename "$f" | cut -d. -f2) 
        echo $id
        hap=$(basename "$f" | cut -d. -f3 | sed 's/\.txt//')  # hap2
        echo $hap
        grep 'Archaic' "$f" | awk '($6>={wildcards.threshold})' | awk -v id="$id" -v hap="$hap" -v count="$hap_count" 'BEGIN {{OFS="\t"}} {{print count, $0, id, hap}}' >> {output.tmp}
        hap_count=$((hap_count+1))
        done
        
        for out in {output.chrs}; do
            header=$(head -n 1 $f)
            echo -e "hapID\t$header\tID\thaplotype"> $out # write out the header
            chr=$(basename "$out" | sed -E 's/.*_([0-9]+)\\.bed/\\1/')
            awk -v i="$chr" '$2 == i {{print}}' {output.tmp} >> "$out"
        done
        """

rule annotate_frags:
    input:
        frags = "simdat/{sim}_{seed}/archaic_0.8/fragments/{sim}_hmmix.{nafr}_{rep}.bed",
        annotations = "simdat/{sim}_{seed}/archaics/{sim}_200k_annotatedArchaics_{rep}.txt",
        script = "helper_scripts/annotate_frags_archaics.py"
    output:
        "simdat/{sim}_{seed}/archaic_0.8/fragments_anno/{sim}_hmmix.{nafr}_annotated_{rep}.bed",
    shell:
        "python {input.script} -b {input.frags} -i {input.annotations} -o {output}"

# need to use this to collage fragments if NEA is included
rule classify_archaic_fragments:
    input:
        bed="simdat/{sim}_{seed}/archaic_0.8/fragments_anno/{sim}_hmmix.{nafr}_annotated_{rep}.bed",
    output:
        den_full="simdat/{sim}_{seed}/archaic_0.8/DEN/fragments/{sim}_hmmix.{nafr}_mostLikely_{rep}.bed",
        nea_full="simdat/{sim}_{seed}/archaic_0.8/NEA/fragments/{sim}_hmmix.{nafr}_mostLikely_{rep}.bed",
        ND01_full="simdat/{sim}_{seed}/archaic_0.8/DEN/fragments/{sim}_hmmix.{nafr}_ND_{rep}.bed",
        ND10_full="simdat/{sim}_{seed}/archaic_0.8/NEA/fragments/{sim}_hmmix.{nafr}_ND_{rep}.bed",
    shell:
        """
        nea=$(head -1 {input.bed} | tr '\t' '\n' | grep -n -x "nea_overlap" | cut -d: -f1)
        den=$(head -1 {input.bed} | tr '\t' '\n' | grep -n -x "den_overlap" | cut -d: -f1)
        nd10=$(head -1 {input.bed} | tr '\t' '\n' | grep -n -x "ND10" | cut -d: -f1)
        nd01=$(head -1 {input.bed} | tr '\t' '\n' | grep -n -x "ND01" | cut -d: -f1)
        awk -v d="$den" -v n="$nea" 'NR==1 || ($d > $n)' {input.bed} > {output.den_full} 
        awk -v d="$den" -v n="$nea" 'NR==1 || ($d < $n)' {input.bed} > {output.nea_full}
        awk -v nd01="$nd01" -v nd10="$nd10" 'NR==1 || ($nd10 == 0 && $nd01 > 0)' {input.bed} > {output.ND01_full}
        awk -v nd01="$nd01" -v nd10="$nd10" 'NR==1 || ($nd01 == 0 && $nd10 > 0)' {input.bed} > {output.ND10_full}
        """

# this is the key of this method, create a matrix with the ancestry of each individual (0, 1, or 2 archaic haplotypes) at each position in the vcf
rule create_ancestry_matrix:
    input:
        bed = "simdat/{sim}_{seed}/{archaic}/fragments/{sim}_{frag_source}_{rep}.bed", 
        vcf = "simdat/{sim}_{seed}/MH_variable/{sim}_sim_{rep}.vcf",
        script = "helper_scripts/create_ancestry_eig.py"
    output:
        "simdat/{sim}_{seed}/{archaic}/curve/fragments/{sim}_{frag_source}_{rep}.anc"
    shell:
        "python {input.script} -b {input.bed} -i {input.vcf} -o {output}"

rule create_ancestry_matrix_ingroup: # this is for hmmix or true ancestry, using positions after removing an outgroup of nafr Africans
    input:
        bed = "simdat/{sim}_{seed}/{archaic}/fragments/{sim}_{frag_source}_{rep}.bed",
        vcf = "simdat/{sim}_{seed}/MH_variable/{sim}_sim_{rep}_{nafr}Afr_ingroup.vcf",
        script = "helper_scripts/create_ancestry_eig.py"
    output:
        "simdat/{sim}_{seed}/{archaic}/curve/fragments/{sim}_{frag_source}_{nafr}Afr_ingroup_{rep}.anc"
    shell:
        "python {input.script} -b {input.bed} -i {input.vcf} -o {output}"

rule create_ancestry_matrix_downsampled: # this is for hmmix or true ancestry, using positions after removing an outgroup of nafr Africans
    input:
        bed = "simdat/{sim}_{seed}/{archaic}/fragments/{sim}_{frag_source}_{rep}.bed",
        vcf = "simdat/{sim}_{seed}/MH_variable/{sim}_sim_{rep}_downsampled.vcf",
        script = "helper_scripts/create_ancestry_eig.py"
    output:
        "simdat/{sim}_{seed}/{archaic}/curve/fragments/{sim}_{frag_source}_ds200k_{rep}.anc"
    shell:
        "python {input.script} -b {input.bed} -i {input.vcf} -o {output}"

# this takes ~24 hrs 
rule calculate_D:
    input:
        anc = "simdat/{sim}_{seed}/{archaic}/curve/fragments/{name}_{rep}.anc",
        script = "helper_scripts/calculate_covariances.py"
    output:
        "simdat/{sim}_{seed}/{archaic}/curve/results/{name}_{rep}_D.txt"
    resources:
        time='30:00:00'
    shell:
        "python {input.script} -i {input.anc} -o {output}"
