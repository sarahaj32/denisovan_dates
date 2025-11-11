wildcard_constraints:
    rep="\d+",

##### CHECK_mutSamps #####
# can now add in laurits' exact simulations
rule sims_all:
    input: # simdat/denisovanNeanderthal_Dlow_2/raw/denisovanNeanderthal_Dlow_sim_16.vcf
        #expand("simdat/{sim}_{seed}/computed/{sim}_sim_allChr_{nafr}Afr_ingroup_{asc}_results.out", rep = list(range(1,21)), sim = ["denisovanNeanderthal_equal", "denisovanNeanderthal_equal"], seed = [1,2,3], nafr = [100], asc = ["ND10", "ND01"]), # simple", "denisovanNeanderthal_Dlow"
        expand("simdat/{sim}_{seed}/computed/{sim}_MH{nMH}sim_allChr_{nafr}Afr_{kind}_{asc}_results.out", rep = list(range(1,21)), sim = ["denisovanNeanderthal_equal", "denisovanNeanderthal_Dlow"], seed = list(range(1,11)), nafr = [100], kind = ["ingroup", "200000refFill", "500000refFill", "900000refFill", "1000000refFill"], asc = ["ND10", "ND01"], nMH = [20, 100])# , "200000refFill", "500000refFill", "900000refFill" #
    output:
        "logs/sim_D.txt"
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

rule find_fixed_AfrAnc:
    input:
        vcf = "simdat/{sim}_{seed}/clean/{sim}_sim_{rep}.vcf"
    output:
        pos = "simdat/{sim}_{seed}/AfrAnc/fixed_AfrAnc_{nafr}Afr_{rep}.txt"
    shell:
        """
        samps=$(printf 'AFR_%d,' {{1..{wildcards.nafr}}} | sed 's/,$//')
        bcftools view -s $samps -C 1 {input.vcf} | bcftools query -f '%CHROM\t%POS\n' -o {output.pos}
        """ 

# we only care about positions variable within non-African modern humans
rule remove_nonvariable_sites:
    input: 
        "simdat/{sim}_{seed}/clean/{sim}_sim_{rep}.vcf"
    output:
        "simdat/{sim}_{seed}/MH{nMH}_variable/{sim}_sim_{rep}.vcf"
    shell:
        """
        samps=$(printf 'NAMH_%d,' {{1..{wildcards.nMH}}} | sed 's/,$//')
        bcftools view -s $samps -c 1 {input} | awk 'BEGIN{{OFS="\t"}} /^#/ {{print; next}} {{ $4="A"; $5="T"; print }}'> {output}
        """

# A rule to downsample to 220,000 SNPs (half of what is simulated)
rule downsample:
    input:
        vcf ="simdat/{sim}_{seed}/archaicsAfrAnc/{sim}_sim_{rep}_{nafr}AfrAnc_{asc}.txt",
    output:
        "simdat/{sim}_{seed}/archaicsAfrAnc/{sim}_sim_{rep}_{nafr}AfrAnc_{ds}Downsampled_{asc}.txt"
        # so we'll have 200k downsampled, 500k, 900
    shell:  
        """
        set -x
        n=$(({wildcards.ds}/20))
        grep -v ^# {input.vcf} | shuf -n $n | sort -nk2 > {output}
        """

# find simulated ND10 and ND01 positions:
rule get_NDxx_pos:
    input:
        vcf = "simdat/{sim}_{seed}/clean/{sim}_sim_{rep}.vcf",
        afr_anc = "simdat/{sim}_{seed}/AfrAnc/fixed_AfrAnc_{nafr}Afr_{rep}.txt",
        script = "helper_scripts/find_NDxx_pos.py"
    output:
        vcf = "simdat/{sim}_{seed}/archaicsAfrAnc/{sim}_sim_{rep}_{nafr}AfrAnc_VinDen.vcf",
        labeled = "simdat/{sim}_{seed}/archaicsAfrAnc/{sim}_sim_{rep}_{nafr}AfrAnc_VinDenlabeled.txt",
        nd01 = "simdat/{sim}_{seed}/archaicsAfrAnc/{sim}_sim_{rep}_{nafr}AfrAnc_ND01.txt",
        nd10 = "simdat/{sim}_{seed}/archaicsAfrAnc/{sim}_sim_{rep}_{nafr}AfrAnc_ND10.txt",
    shell:
        """
        bcftools view -c1 -s Vindija33.19,Denisova -T {input.afr_anc} {input.vcf} -o {output.vcf}
        python {input.script} -i {output.vcf} -n Vindija33.19 -d Denisova -r -o {output.labeled}
        grep -i nd01 {output.labeled} | cut -f1,2 > {output.nd01}
        grep -i nd10 {output.labeled} | cut -f1,2 > {output.nd10}
        """

# these are the positions that are variable in Neanderthal or Denisovan and homozygous ancestral in Africans, and segregating in the MH ingroup
rule filter_MH_to_NDxx:
    input:
        vcf = "simdat/{sim}_{seed}/MH{nMH}_variable/{sim}_sim_{rep}.vcf",
        nd = "simdat/{sim}_{seed}/archaicsAfrAnc/{sim}_sim_{rep}_{nafr}AfrAnc_{asc}.txt"
    output:
        "simdat/{sim}_{seed}/{asc}/{sim}_MH{nMH}sim_{rep}_{nafr}Afr_ingroup_{asc}.vcf",
    shell:
        """
        bcftools view -T {input.nd} {input.vcf} -o {output}
        """

# need to add on a refFill step here - add on all other positions as reference homozygous

rule fill_in_homRef:
    input:
        vcf="simdat/{sim}_{seed}/{asc}/{sim}_MH{nMH}sim_{rep}_{nafr}Afr_ingroup_{asc}.vcf",
        pos="simdat/{sim}_{seed}/archaicsAfrAnc/{sim}_sim_{rep}_{nafr}AfrAnc_{ds}Downsampled_{asc}.txt",
        script="helper_scripts/add_homref_2vcf.py"
    output:
        unzipped = "simdat/{sim}_{seed}/{asc}/{sim}_MH{nMH}sim_{rep}_{nafr}Afr_{ds}refFill_{asc}.vcf"
    shell:
        """
        python {input.script} -v {input.vcf} -p {input.pos} -o {output.unzipped}
        """

rule concatenate_sim_chroms:
    input: 
        expand("simdat/{{sim}}_{{seed}}/{{asc}}/{{sim}}_MH{{nMH}}sim_{rep}_{{nafr}}Afr_{{kind}}_{{asc}}.vcf", rep = list(range(1,21)))
    output:
        "simdat/{sim}_{seed}/{asc}/{sim}_MH{nMH}sim_allChr_{nafr}Afr_{kind}_{asc}.vcf"
    resources:
        cpus=4
    shell:
        "bcftools concat {input} -o {output}"

# convert to eigenstrat format
rule convert_sims_to_eig:
    input:
        vcf="simdat/{sim}_{seed}/{asc}/{sim}_MH{nMH}sim_allChr_{nafr}Afr_{kind}_{asc}.vcf",
        script="/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/bin/vcf2eigenstrat.py"
    output:
        geno="simdat/{sim}_{seed}/{asc}/{sim}_MH{nMH}sim_allChr_{nafr}Afr_{kind}_{asc}.geno",
        ind="simdat/{sim}_{seed}/{asc}/{sim}_MH{nMH}sim_allChr_{nafr}Afr_{kind}_{asc}.ind",
        snp="simdat/{sim}_{seed}/{asc}/{sim}_MH{nMH}sim_allChr_{nafr}Afr_{kind}_{asc}.snp",
        genmap_snp="simdat/{sim}_{seed}/{asc}/{sim}_MH{nMH}sim_allChr_{nafr}Afr_{kind}_{asc}_genmap.snp",
    params:
        path="simdat/{sim}_{seed}/{asc}/{sim}_MH{nMH}sim_allChr_{nafr}Afr_{kind}_{asc}"
    resources:
        time='6:00:00',
        cpus=4
    shell:
        """ 
        python2 {input.script} -v {input.vcf} -o {params.path}
        awk '{{genmap=$4/1e8; print $1"\t"$2"\t"genmap"\t"$4}}' {output.snp} > {output.genmap_snp}
        """ 


rule make_computed_parfile:
    input:
        geno="simdat/{sim}_{seed}/{asc}/{sim}_MH{nMH}sim_allChr_{nafr}Afr_{kind}_{asc}.geno",
        ind="simdat/{sim}_{seed}/{asc}/{sim}_MH{nMH}sim_allChr_{nafr}Afr_{kind}_{asc}.ind",
        snp="simdat/{sim}_{seed}/{asc}/{sim}_MH{nMH}sim_allChr_{nafr}Afr_{kind}_{asc}_genmap.snp",
    params:
        out="simdat/{sim}_{seed}/computed/{sim}_MH{nMH}sim_allChr_{nafr}Afr_{kind}_{asc}_results.out",
        #jk_out="data/{dataset}/{popu}/{kind}/dating/jacknife/"
    output:
        parfile="simdat/{sim}_{seed}/computed/{sim}_MH{nMH}sim_allChr_{nafr}Afr_{kind}_{asc}_computed_parfile",
    shell:
        """
        echo "genotypename: {input.geno}" > {output.parfile}
        echo "snpname: {input.snp}" >> {output.parfile}
        echo "indivname: {input.ind}" >> {output.parfile}
        echo "maxdis: 0.01">> {output.parfile}
        echo "binsize: 0.00001">> {output.parfile}
        echo "output: {params.out}" >> {output.parfile}
        """

rule computed:
    input:
        geno="simdat/{sim}_{seed}/{asc}/{sim}_MH{nMH}sim_allChr_{nafr}Afr_{kind}_{asc}.geno",
        ind="simdat/{sim}_{seed}/{asc}/{sim}_MH{nMH}sim_allChr_{nafr}Afr_{kind}_{asc}.ind",
        snp="simdat/{sim}_{seed}/{asc}/{sim}_MH{nMH}sim_allChr_{nafr}Afr_{kind}_{asc}_genmap.snp",
        parfile="simdat/{sim}_{seed}/computed/{sim}_MH{nMH}sim_allChr_{nafr}Afr_{kind}_{asc}_computed_parfile",
    output:
        out="simdat/{sim}_{seed}/computed/{sim}_MH{nMH}sim_allChr_{nafr}Afr_{kind}_{asc}_results.out"
    resources:
        time="30:00:00",
        cpus="2"
    shell:
        "/global/home/users/moorjani/bin/computed -p {input.parfile}"
