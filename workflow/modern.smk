# a rule to find all the variable positions within a population
rule find_segregating:
    input:
        vcf=lambda wildcards: config["datasets"][wildcards.dataset]["hg19_path"].format(chrom=wildcards.chrom),
        inds = "data/{dataset}/inds/{popu}.inds"
    output:
        "data/{dataset}/{popu}/modern/vcf/hg19/{dataset}_snps10_08_allChr_modern_{popu}_var_chr{chrom}.vcf.gz"
        "data/{dataset}/{popu}/{kind}/vcf/hg19/{dataset}_snps10_08_allChr_{kind}_{popu}.vcf.gz", 
    shell:
        "bcftools view -S {input.inds} -M2 -v snps {input.vcf} | bcftools +fill-tags -Ou -- -t AC,AN,AF | bcftools view -i 'INFO/AF>0 && INFO/AF<1' -Oz -o {output}"

# a set of rules to find all ND01 and ND10 positions in a dataset, 
# including adding in positions that are homozygous reference

# this filters all 1000G positions to only positions in Nea or Denisova (ND11 + ND01 + ND10)
rule filter_modern_to_ND11:
    input:
        vcf=lambda wildcards: config["datasets"][wildcards.dataset]["hg19_path"].format(chrom=wildcards.chrom),
        ND11="data/hg19/snp_annotations/ND11_hg19_masked.txt"
    output:
        "data/{dataset}/ND_VCF/{dataset}_ND11pos_hg19_chr{chrom}.vcf"
    shell:
        "bcftools view {input.vcf} -T {input.ND11} -M2 -v snps -o {output}"

rule fill_in_homRef:
    input:
        vcf="data/{dataset}/ND_VCF/{dataset}_ND11pos_hg19_chr{chrom}.vcf",
        ref="/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/resources/hg19/fasta/hg19.fa",
        pos="data/hg19/snp_annotations/ND11_hg19_masked.bed",
        script="helper_scripts/add_homref_2vcf.py"
    output:
        unzipped = temp("data/{dataset}/ND_VCF_refFill_Masked/{dataset}_ND11pos_hg19_modernRefFill_chr{chrom}.vcf"),
        zipped = "data/{dataset}/ND_VCF_refFill_Masked/{dataset}_ND11pos_hg19_modernRefFill_chr{chrom}.vcf.gz"
    shell:
        """
        python {input.script} -v {input.vcf} -r {input.ref} -p {input.pos} -o {output.unzipped}
        bcftools view {output.unzipped} -Oz -o {output.zipped}
        """
# check to confirm that these are the same positions: 
# bcftools view -H 1000G_ND11pos_hg19_modernRefFill_chr17.vcf.gz | grep -v -i reffill | wc -l

rule count_lines:
    input:
        reffill=expand("data/{{dataset}}/ND_VCF_refFill_Masked/{{dataset}}_ND11pos_hg19_modernRefFill_chr{chrom}.vcf.gz", chrom = list(range(1,23))),
        orig=expand("data/{{dataset}}/ND_VCF/{{dataset}}_ND11pos_hg19_chr{chrom}.vcf", chrom = list(range(1,23))),
        script="/global/scratch/p2p3/pl1_moorjani/sarahj32/superarchaic_introgression/superarchaic_analysis/helper_scripts/count_vcf_lines.sh"
    output:
        reffill="data/{dataset}/ND_VCF_refFill_Masked/summary_counts.out",
        orig="data/{dataset}/ND_VCF/summary_counts.out"
    shell:
        """
        {input.script} {input.reffill[0]} {output.reffill}
        {input.script} {input.orig[0]} {output.orig}
        """
        
rule combine_modern:
    input:
        vcf=expand("data/{{dataset}}/ND_VCF_refFill_Masked/{{dataset}}_ND11pos_hg19_modernRefFill_chr{chrom}.vcf.gz", chrom=list(range(1,23))),
        chr_map="data/resources/chr_map_addChr"
    output:
        "data/{dataset}/ND_VCF_refFill_Masked/{dataset}_ND11pos_hg19_modernRefFill_allChr.vcf.gz"
    shell:
        "bcftools concat {input.vcf} | bcftools annotate --rename-chrs {input.chr_map} -Oz -o {output}"

# find the population frequency at each of these position
rule get_popu_ND_freqs:
    input:
        vcf="data/{dataset}/ND_VCF_refFill_Masked/{dataset}_ND11pos_hg19_modernRefFill_allChr.vcf.gz",
        inds="data/{dataset}/inds/{popu}.inds"
    output:
        "data/{dataset}/ND_VCF_refFill_Masked/{popu}/{dataset}_{popu}_ND11pos_hg19_modernRefFill_allChr_freqs.txt",
    shell:
        "bcftools view -S {input.inds} {input.vcf} | bcftools +fill-tags -Ou -- -t AC,AN,AF | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AC\t%AN\t%AF\n' > {output}" 

rule get_popu_all_freqs:
    input:
        vcf=lambda wildcards: hg19_vcf_path.format(chromosome=wildcards.chrom),
        inds="data/{dataset}/inds/{popu}.inds",
        mask="/global/scratch/p2p3/pl1_moorjani/sarahj32/data/resources/hg19/1000G/callability_masks/20141020.pilot_mask.whole_genome_noChr.bed"
    output:
        "data/{dataset}/VCF_frequencies/{popu}/{dataset}_{popu}_hg19_chr{chrom}_pilotMask_freqs.txt",
    shell:
        "bcftools view -S {input.inds} -R {input.mask} -M2 -v snps {input.vcf} | bcftools +fill-tags -Ou -- -t AC,AN,AF | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AC\t%AN\t%AF\n' > {output}" 

# filter the vcf to positions of interest
rule select_popu_modern:
    input:
        vcf="data/{dataset}/ND_VCF_refFill_Masked/{dataset}_ND11pos_hg19_modernRefFill_allChr.vcf.gz",
        inds = "data/{dataset}/inds/{popu}.inds"
    output:
        vcf = "data/{dataset}/{popu}/modernRefFill/vcf/hg19/{dataset}_snps10_08_allChr_modernRefFill_{popu}.vcf.gz"
    shell:
        "bcftools view {input.vcf} -S {input.inds} | bcftools +fill-tags -Ou -- -t AC,AN,AF| bcftools view -Oz -o {output}"

# a rule to combine all chromosomes (these should be relatively small by this point)
rule combine_chrs:
    input: 
        vcf=expand("data/{{dataset}}/{{popu}}/modern/vcf/hg19/{{dataset}}_snps10_08_allChr_modern_{{popu}}_var_chr{chrom}.vcf.gz", chrom = list(range(1,23))),
        chr_map="data/resources/chr_map_addChr"
    output:
        combined_vcf="data/{dataset}/{popu}/modern/vcf/hg19/{dataset}_snps10_08_allChr_modern_{popu}.vcf.gz"
    resources:
        cpus='4'
    shell:
        """
        mkdir -p $TMPDIR/{wildcards.dataset}_{wildcards.popu}_combine
        bcftools concat {input.vcf} | bcftools sort -T $TMPDIR/{wildcards.dataset}_{wildcards.popu}_combine | bcftools annotate --rename-chrs {input.chr_map} -o {output.combined_vcf}
        """

rule find_dataset_pos:
    input:
        vcf = "data/{dataset}/ND_VCF_refFill_Masked/{dataset}_ND11pos_hg19_modernRefFill_allChr.vcf.gz",
        pos="data/hg19/snp_annotations/{asc}_hg19.txt",
    output:
        "data/{dataset}/ND_VCF_refFill_Masked/{dataset}_{asc}_pos_hg19_modernRefFill_allChr_positions.vcf"
    shell: 
        "bcftools view {input.vcf} -T {input.pos} -o {output}"

         #"data/1000G/ND_VCF_refFill_Masked/1000G_ND01AfrAnc_pos_hg19_modernRefFill_allChr_positions.vcf"
         # "data/1000G/ND_VCF_refFill_Masked/1000G_ND10AfrAnc_pos_hg19_modernRefFill_allChr_positions.vcf"