# this is a set of snakemake rules that get ND01 and ND10 positions from modern hg38 datasets, uses my ascertained archaic SNPs
# and then filters the modern datasets to just those positions, concatenates all chromosomes together, and liftover to hg19
# this output is then processed by date_flexi.smk to get the final date estimates

# get positions of just NDxx SNPs on both builds
rule get_NDxx:
    input:
        "data/Archaics/{build}/AltDenVinPanHSanc_merged_strictArchaic_{build}_full_annotated_chimpAnc.txt"
    output:
        "data/Archaics/{build}/{type}pos_{build}_chimpAnc.txt"
    shell:
        "grep {wildcards.type} {input} > {output}"

# filter the dataset to only ND01 or ND10 positions
rule filter:
    input:
        vcf=lambda wildcards: config["datasets"][wildcards.dataset]["hg38_path"].format(chromosome=wildcards.chromosome),
        inds="data/{dataset}/inds/{popu}.inds",
        pos="data/Archaics/hg38/{type}pos_hg38_chimpAnc.txt",
        chr_map="data/resources/chr_map"
    output:
        vcf="data/{dataset}/{popu}/modernMyAsc/vcf/hg38/{dataset}_modernMyAsc_{popu}_{type}_chr{chromosome}.vcf"
    shell:
        """
        bcftools view -S {input.inds} -T {input.pos} -M2 -v snps {input.vcf} -Ou | bcftools annotate --rename-chrs {input.chr_map} -o {output.vcf}
        bcftools index {output.vcf}
        """

# concatentate all chromosomes together
rule concatenate:
    input:
        vcf=expand("data/{{dataset}}/{{popu}}/modernMyAsc/vcf/hg38/{{dataset}}_modernMyAsc_{{popu}}_{{type}}_chr{chrom}.vcf.gz", chrom=list(range(1,23))),
    output:
        "data/{dataset}/{popu}/modernMyAsc/vcf/hg38/{dataset}_modernMyAsc_{popu}_{type}_allChr.vcf.gz"
    shell:
        "bcftools concat {input.vcf} -Oz -o {output}"

# convert to snp so we can interpit and figure out which positions are not in order
rule convert_archaics_to_eig:
    input:
        vcf="data/Archaics/hg19/AltDenVinPanHSanc_merged_strictArchaic_hg19_chr{chromosome}.vcf",
        script="/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/bin/vcf2eigenstrat.py"
    output:
        geno="data/Archaics/hg19/AltDenVinPanHSanc_merged_strictArchaic_hg19_chr{chromosome}.geno",
        ind="data/Archaics/hg19/AltDenVinPanHSanc_merged_strictArchaic_hg19_chr{chromosome}.ind",
        snp="data/Archaics/hg19/AltDenVinPanHSanc_merged_strictArchaic_hg19_chr{chromosome}.snp"
    params:
        path="data/Archaics/hg19/AltDenVinPanHSanc_merged_strictArchaic_hg19_chr{chromosome}"
    shell:
        """ 
        python2 {input.script} -v {input.vcf} -o {params.path}
        sed -i 's/chr//g' {output.snp} 
        """ 
    
rule interpit_archaics:
    # cannot have "chr" in the .snp file
    input:
        snp="data/Archaics/hg19/AltDenVinPanHSanc_merged_strictArchaic_hg19_chr{chromosome}.snp",
        genmap="/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/hg19/resources/genetic_maps_corrections/hg19/Shared_AfAmMap"
    output:
        snp="data/Archaics/hg19/AltDenVinPanHSanc_merged_strictArchaic_hg19_chr{chromosome}_genmap.snp",
    shell:
        "/global/home/users/moorjani/bin/interpit -i {input.snp}  -o {output.snp}  -d {input.genmap}"

rule extract_problem_positions:
    input:
        snp="data/Archaics/hg19/AltDenVinPanHSanc_merged_strictArchaic_hg19_chr{chromosome}_genmap.snp",
        script="testing/computed_testing/find_continuous_positions.py"
    output:
        "data/Archaics/hg19/AltDenVinPanHSanc_merged_strictArchaic_hg19_chr{chromosome}_genmap_problemPositions.txt"
    shell:
        "python {input.script} -s {input.snp} -p {output}"

rule concat_problem_positions:
    input:
        expand("data/Archaics/hg19/AltDenVinPanHSanc_merged_strictArchaic_hg19_chr{chromosome}_genmap_problemPositions.txt", chromosome=list(range(1,23)))
    output:
        "data/Archaics/hg19/AltDenVinPanHSanc_merged_strictArchaic_hg19_genmap_problemPositions.txt"
    shell:
        "cat {input} > {output}"

# liftover to hg19
rule liftover:
    input:
        vcf="data/{dataset}/{popu}/modernMyAsc/vcf/hg38/{dataset}_modernMyAsc_{popu}_{type}_chr{chrom}.vcf.gz",
        chainfile="/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/resources/liftover.chainfiles/hg38ToHg19.over.chain",
        ref="/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/resources/hg19/fasta/hg19.fa",
    output:
        vcf_tmp=temp("data/{dataset}/{popu}/modernMyAsc/vcf/hg19/{dataset}_modernMyAsc_{popu}_{type}_chr{chrom}.TMP.vcf"),
        vcf="data/{dataset}/{popu}/modernMyAsc/vcf/hg19/{dataset}_modernMyAsc_{popu}_{type}_chr{chrom}_allPos.vcf",
    resources:
        cpus='6'
    shell:
        """
        CrossMap vcf {input.chainfile} {input.vcf} {input.ref} {output.vcf_tmp} --no-comp-alleles 
        bcftools sort {output.vcf_tmp} -o {output.vcf}
        """

rule remove_problem_positions:
    input:
        vcf="data/{dataset}/{popu}/modernMyAsc/vcf/hg19/{dataset}_modernMyAsc_{popu}_{type}_chr{chrom}_allPos.vcf",
        problem="data/Archaics/hg19/AltDenVinPanHSanc_merged_strictArchaic_hg19_genmap_problemPositions.txt"
    output:
        vcf="data/{dataset}/{popu}/modernMyAsc/vcf/hg19/{dataset}_modernMyAsc_{popu}_{type}_chr{chrom}.vcf"
    shell:
        """
        bcftools view -T ^{input.problem} {input.vcf}  -o {output.vcf}
        """
