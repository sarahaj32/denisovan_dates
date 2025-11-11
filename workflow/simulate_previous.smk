wildcard_constraints:
    rep="\d+",
    frag_source="(sim_intro|hmmix.100)"

# can now add in laurits' exact simulations
rule sims_all:
    input:
        #expand("simdat/{sim}_{seed}/archaicVar/{sim}_sim_{rep}_archaicVar_annotated.txt", rep = list(range(1,21))),
        #expand("simdat/{sim}_{seed}/vcf/{sim}_sim_{type}_allChr.vcf", type = ["ND01", "ND10"]),
        #expand("simdat/{sim}_{seed}/AfrFreq/{sim}_sim_{rep}_AfrFreqs.txt", rep = list(range(1,21)), sim = "denisovan_simple"),
        #expand("simdat/{sim}_{seed}/dating/computed/{sim}_sim_{rep}_{type}_dateCorrections.out", type = ["ND01", "ND10"], rep = list(range(1,21))),
        # expand("simdat/{sim}_{seed}/curve/results/{sim}_{source}{ds}_{rep}_D.txt", rep = list(range(1,21)), sim = ["denisovan_simple", "denisovan_early", "denisovan_late"], source = ["sim_intro", "hmmix.100"], ds = [ "_ingroup"], seed = [1,2,3]), # "", "_downsampled",  , "denisovan_low"
        expand("simdat/{sim}_{seed}/curve/results/{sim}_{source}{ds}_{rep}_D.txt", rep = list(range(1,21)), sim = ["denisovan_simple", "denisovan_early", "denisovan_late"], source = ["sim_intro"], ds = [ "_ingroup"], seed = [1,2,3]), # "", "_downsampled",  , "denisovan_low"

        # SARAH TEHSE ARE WHAT I WNAT!!!!
        # expand("simdat/{sim}_{seed}/curve/fragments/{sim}_{frag_source}{ds}_{rep}.anc", rep = list(range(1,21)), sim = ["denisovan_simple", "denisovan_early", "denisovan_late", "denisovan_low"], frag_source = ["sim_intro", "hmmix.100"], ds = ["", "_downsampled", "_ingroup"], seed = [1,2,3]),
        # expand("simdat/{sim}_{seed}/curve/results/{sim}_{source}_downsampled_{rep}_D.txt", rep = list(range(1,21)), sim = ["denisovan_simple", "denisovan_early", "denisovan_late", "denisovan_low"], source = ["sim_intro",]),
        #expand("simdat/{sim}_{seed}/curve/fragments/{sim}_{frag_source}_{rep}.anc", rep = list(range(1,21)), sim = ["denisovan_simple", "denisovan_early", "denisovan_late", "denisovan_low"], frag_source = ["sim_intro", "hmmix.100"]),
        # a test comparing the # of africans in the outgroup and the impact on transitions, emissions, and decoded/inferred fragments
        expand("simdat/{sim}_{seed}/fragments/{sim}_{method}.{nafr}_{rep}.bed", rep = list(range(1,21)), sim = ["denisovan_simple", "denisovan_smallAfr", "denisovan_bigNonAfr", "denisovan_LS_ice", "denisovan_LS_pap", "denisovan_simpleLS"], nafr = [100, 200, 500], seed = [1], method = ["hmmix"]), # increase to 500???? , "hmmix_LS"
        # expand("simdat/{sim}_{seed}/fragments/{sim}_sim_intro_{rep}.bed", rep = list(range(1,21)), sim = ["denisovan_simple", "denisovan_smallAfr", "denisovan_bigNonAfr", "denisovan_LS_ice", "denisovan_LS_pap", "denisovan_simpleLS"], seed = [1])

    output:
        "logs/sims_all_done.txt"
    shell:
        "echo 'all done' > {output}"

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

# a step to we remove recurring mutations (otherwise this causes issues in the eigenstrate files and computed analysis)
rule remove_recurring_mutations:
    input: 
        "simdat/{sim}_{seed}/raw/{sim}_sim_{rep}.vcf"
    output:
        "simdat/{sim}_{seed}/clean/{sim}_sim_{rep}.vcf"
    shell:
        "bcftools view -c 1 {input} | bcftools norm -m + | awk '$3 !~ /;/' > {output}"

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

# find the introgressed fragments from Denisovans into humans, using the simulated tree
rule find_intro_frags:
    input:
        tree = "simdat/{sim}_{seed}/raw/{sim}_sim_{rep}.tree",
        script="/global/scratch/p2p3/pl1_moorjani/sarahj32/superarchaic_introgression/super_clean/helper_scripts/find_intro_frags.py"
    output:
        bed = "simdat/{sim}_{seed}/fragments/{sim}_sim_intro_{rep}.bed"
    conda: 
        "sim.yaml"
    resources:
        time='6:00:00'
    shell:
        """
        python {input.script} -f {input.tree} -b {output.bed} -s DEN -t OOA -c {wildcards.rep} -n 20
        """

# # A rule to downsample to 220,000 SNPs (half of what is simulated)
# rule downsample:
#     input:
#         vcf = "simdat/{sim}_{seed}/MH_variable/{sim}_sim_{rep}.vcf",
#     output:
#         "simdat/{sim}_{seed}/MH_variable/{sim}_sim_{rep}_downsampled.vcf"
#     shell:  
#         """
#         set -x
#         grep -v ^# {input.vcf} | shuf -n 200000 | sort -nk2 > {output}
#         """

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
        obs = "simdat/{sim}_{seed}/{method}/obs/obs.{nafr}.NAMH_{ind}.txt"
    output:
        "simdat/{sim}_{seed}/{method}/trained/trained.{nafr}.NAMH_{ind}.json"      
    shell:
        "hmmix train -obs={input.obs} -out={output} -haploid"

rule decode:
    input:
        obs = "simdat/{sim}_{seed}/{method}/obs/obs.{nafr}.NAMH_{ind}.txt",
        trained = "simdat/{sim}_{seed}/{method}/trained/trained.{nafr}.NAMH_{ind}.json"
    output:
        expand("simdat/{{sim}}_{{seed}}/{{method}}/decoded/decoded.{{nafr}}.NAMH_{{ind}}.hap{hap}.txt", hap = [1,2])
    params:
        out_path = "simdat/{sim}_{seed}/{method}/decoded/decoded.{nafr}.NAMH_{ind}"     
    shell:
        "hmmix decode -obs={input.obs} -param={input.trained} -haploid -out={params.out_path}"

rule collate_frags:
    input:
        expand("simdat/{{sim}}_{{seed}}/{{method}}/decoded/decoded.{{nafr}}.NAMH_{ind}.hap{hap}.txt", ind = list(range(1,21)), hap = [1,2])
    output:
        tmp = temp("simdat/{sim}_{seed}/fragments/{sim}_{method}.{nafr}_all.bed"),
        chrs = expand("simdat/{{sim}}_{{seed}}/fragments/{{sim}}_{{method}}.{{nafr}}_{rep}.bed", rep = list(range(1,21)))
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
        grep 'Archaic' "$f" | awk '($6>=0.8)' | awk -v id="$id" -v hap="$hap" -v count="$hap_count" 'BEGIN {{OFS="\t"}} {{print count, $0, id, hap}}' >> {output.tmp}
        hap_count=$((hap_count+1))
        done

        for i in {{1..20}};
        do
        echo $i
        out="simdat/{wildcards.sim}_{wildcards.seed}/fragments/{wildcards.sim}_{wildcards.method}.{wildcards.nafr}_$i.bed"
        echo $out
        awk -v i="$i" '$2 == i {{print}}' {output.tmp} > "$out"
        done
        """

# this is the key of this method, create a matrix with the ancestry of each individual (0, 1, or 2 archaic haplotypes) at each position in the vcf
rule create_ancestry_matrix:
    input:
        bed = "simdat/{sim}_{seed}/fragments/{sim}_{frag_source}_{rep}.bed", 
        vcf = "simdat/{sim}_{seed}/MH_variable/{sim}_sim_{rep}.vcf",
        script = "helper_scripts/create_ancestry_eig.py"
    output:
        "simdat/{sim}_{seed}/curve/fragments/{sim}_{frag_source}_{rep}.anc"
    shell:
        "python {input.script} -b {input.bed} -i {input.vcf} -o {output}"

# rule create_ancestry_matrix_downsampled:
#     input:
#         bed = "simdat/{sim}_{seed}/fragments/{sim}_{frag_source}_{rep}.bed",
#         vcf = "simdat/{sim}_{seed}/MH_variable/{sim}_sim_{rep}_downsampled.vcf",
#         script = "helper_scripts/create_ancestry_eig.py"
#     output:
#         "simdat/{sim}_{seed}/curve/fragments/{sim}_{frag_source}_downsampled_{rep}.anc"
#     shell:
#         "python {input.script} -b {input.bed} -i {input.vcf} -o {output}"


rule create_ancestry_matrix_ingroup:
    input:
        bed = "simdat/{sim}_{seed}/fragments/{sim}_{frag_source}_{rep}.bed",
        vcf = "simdat/{sim}_{seed}/MH_variable/{sim}_sim_{rep}_ingroup.vcf",
        script = "helper_scripts/create_ancestry_eig.py"
    output:
        "simdat/{sim}_{seed}/curve/fragments/{sim}_{frag_source}_ingroup_{rep}.anc"
    shell:
        "python {input.script} -b {input.bed} -i {input.vcf} -o {output}"

# this takes ~24 hrs 
rule calculate_D:
    input:
        anc = "simdat/{sim}_{seed}/curve/fragments/{name}_{rep}.anc",
        script = "helper_scripts/calculate_covariances.py"
    output:
        "simdat/{sim}_{seed}/curve/results/{name}_{rep}_D.txt"
    resources:
        time='30:00:00'
    shell:
        "python {input.script} -i {input.anc} -o {output}"

# # a rule to find the variable archaic positions
# rule get_archaic_variable:
#     input: 
#         "simdat/{sim}_{seed}/clean/{sim}_sim_{rep}.vcf"
#     output:
#         "simdat/{sim}_{seed}/archaicVar/{sim}_sim_{rep}_archaicVar.vcf"
#     shell:
#         "bcftools view -s Chagyrskaya-Phalanx,Denisova,AltaiNeandertal,Vindija33.19 {input} | bcftools view -c 1  > {output}"

# # annotate the archaic variants
# # and separate into ND01 and ND10 files
# rule annotate_archaic_variants_sim:
#     input:
#         archaics="simdat/{sim}_{seed}/archaicVar/{sim}_sim_{rep}_archaicVar.vcf",
#         script="add_annotations_to_archaic_observations.py"
#     output:
#         annotated="simdat/{sim}_{seed}/archaicVar/{sim}_sim_{rep}_archaicVar_annotated.txt",
#         ND01="simdat/{sim}_{seed}/archaicVar/{sim}_sim_{rep}_ND01.txt",
#         ND10="simdat/{sim}_{seed}/archaicVar/{sim}_sim_{rep}_ND10.txt",
#     shell:
#         """
#         python {input.script} -i {input.archaics} -o {output.annotated}
#         awk '$12 == "ND01"' {output.annotated} | cut -f1,2 > {output.ND01}
#         awk '$12 == "ND10"' {output.annotated} | cut -f1,2 > {output.ND10}
#         """

# # get AFR frequencies at every position
# rule get_AFR_freqs:
#     input:
#         vcf="simdat/{sim}_{seed}/clean/{sim}_sim_{rep}.vcf",
#         pos="simdat/{sim}_{seed}/archaicVar/{sim}_sim_{rep}_archaicVar.vcf"
#     output:
#         "simdat/{sim}_{seed}/AfrFreq/{sim}_sim_{rep}_AfrFreqs.txt"
#     shell:
#         "bcftools view -T {input.pos} -s $(echo AFR_{{1..400}} | tr ' ' ',') {input.vcf} |  bcftools +fill-tags -Ou -- -t AC,AN,AF | bcftools query -f '%CHROM\t%POS\t%AC\t%AN\t%AF\n' > {output}"

# # filter the dataset to only ND01 or ND10 positions that are segeregating in the modern non-africans
# rule filter_sims:
#     input:
#         vcf="simdat/{sim}_{seed}/clean/{sim}_sim_{rep}.vcf",
#         pos="simdat/{sim}_{seed}/archaicVar/{sim}_sim_{rep}_{type}.txt",
#     output:
#         vcf="simdat/{sim}_{seed}/vcf/{sim}_sim_{rep}_{type}.vcf"
#     shell:
#         """
#         bcftools view -s $(echo NAMH_{{1..200}} | tr ' ' ',') -T {input.pos} -c1 {input.vcf} -o {output.vcf}
#         """

# # # concatentate all chromosomes together
# # rule concatenate_sims:
# #     input:
# #         vcf=expand("simdat/{sim}_{seed}/vcf/{sim}_sim_{rep}_{{type}}.vcf.gz", rep=list(range(1,23))),
# #     output:
# #         "simdat/{sim}_{seed}/vcf/{sim}_sim_{type}_allChr.vcf"
# #     shell:
# #         "bcftools concat {input.vcf} -o {output}"

# # rule convert_sims_to_eig:
# #     input:
# #         vcf="simdat/{sim}_{seed}/vcf/{sim}_sim_{type}_allChr.vcf",
# #         script="/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/bin/vcf2eigenstrat.py"
# #     output:
# #         geno="simdat/{sim}_{seed}/dating/eig/{sim}_sim_{type}_allChr.geno",
# #         ind="simdat/{sim}_{seed}/dating/eig/{sim}_sim_{type}_allChr.ind",
# #         snp="simdat/{sim}_{seed}/dating/eig/{sim}_sim_{type}_allChr.snp",
# #         genmap_snp="simdat/{sim}_{seed}/dating/eig/{sim}_sim_{type}_allChr_genmap.snp"
# #     params:
# #         path="simdat/{sim}_{seed}/dating/eig/{sim}_sim_{type}_allChr"
# #     shell:
# #         """ 
# #         python2 {input.script} -v {input.vcf} -o {params.path}
# #         sed -i 's/chr//g' {output.snp} 
# #         awk '{{dist = $4/100000000}} {{print $1, $2, dist, $4, $5, $6}}' {output.snp} > {output.genmap_snp}
# #         """ 
# #         # add on a step here where we create the genetic map based on positions

# rule convert_sims_to_eig:
#     input:
#         vcf="simdat/{sim}_{seed}/vcf/{sim}_sim_{rep}_{type}.vcf",
#         script="/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/bin/vcf2eigenstrat.py"
#     output:
#         geno="simdat/{sim}_{seed}/dating/eig/{sim}_sim_{rep}_{type}.geno",
#         ind="simdat/{sim}_{seed}/dating/eig/{sim}_sim_{rep}_{type}.ind",
#         snp="simdat/{sim}_{seed}/dating/eig/{sim}_sim_{rep}_{type}.snp",
#         genmap_snp="simdat/{sim}_{seed}/dating/eig/{sim}_sim_{rep}_{type}_genmap.snp"
#     params:
#         path="simdat/{sim}_{seed}/dating/eig/{sim}_sim_{rep}_{type}"
#     shell:
#         """ 
#         python2 {input.script} -v {input.vcf} -o {params.path}
#         sed -i 's/chr//g' {output.snp} 
#         awk '{{dist = $4/100000000}} {{print $1, $2, dist, $4, $5, $6}}' {output.snp} > {output.genmap_snp}
#         """ 
#         # add on a step here where we create the genetic map based on positions


# I don't actually want to use the hmmix scripts for outgroup or ingroup creation

# # make the json of ingroup and outgroup samples
# rule make_hmmix_samples:
#     input:
#         script = "helper_scripts/create_hmmix_samples.py"
#     output:
#         "simdat/{sim}_{seed}/hmmix/hmmix_samples_{nafr}Afr.json"
#     shell:
#         # for now, we'll stick with 20 ingroup and 200 outgroup Africans (this is comparable to real data)
#         "python {input.script} -o {output} -ig NAMH -ic 20 -og AFR -oc {wildcards.nafr}"

# # also do a create_outgroup LS and create_ingroup LS using Laurits' code
# # maybe i'm including too many positions here?!?!?!

# rule create_outgroup_hmmix:
#     input:
#         vcf = expand("simdat/{{sim}}_{{seed}}/clean/{{sim}}_sim_{rep}.vcf", rep = list(range(1,21))),
#         samples = "simdat/{sim}_{seed}/hmmix/hmmix_samples_{nafr}Afr.json"
#     output:
#         out = "simdat/{sim}_{seed}/hmmix_LS/outgroup/outgroup_{nafr}Afr.txt"
#     params:
#         files = "simdat/{sim}_{seed}/clean/{sim}_sim_*.vcf"
#     shell:
#         """
#         hmmix create_outgroup -ind={input.samples} -vcf={params.files} -out={output.out}
#         """

# rule create_ingroup_hmmix:
#     input:
#         vcf = expand("simdat/{{sim}}_{{seed}}/MH_variable/{{sim}}_sim_{rep}.vcf", rep = list(range(1,21))),
#         samples = "simdat/{sim}_{seed}/hmmix/hmmix_samples_{nafr}Afr.json",
#         outgroup = "simdat/{sim}_{seed}/hmmix_LS/outgroup/outgroup_{nafr}Afr.txt"
#     output:
#         expand("simdat/{{sim}}_{{seed}}/hmmix_LS/obs/obs.{{nafr}}.NAMH_{ind}.txt", ind = list(range(1,21)))
#     params:
#         ingroup_path = "simdat/{sim}_{seed}/hmmix_LS/obs/obs.{nafr}",
#         files = "simdat/{sim}_{seed}/MH_variable/{sim}_sim_*.vcf"
#     shell:
#         """
#         hmmix create_ingroup -ind={input.samples} -vcf={params.files} -outgroup={input.outgroup} -out={params.ingroup_path}
#         """
