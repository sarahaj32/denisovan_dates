rule find_all_ND_pos:
    # extract just the chrom and position for all ND positions
    # then extract nd01, nd10 and nd11 positions
    # there are several of these positions where the genetic map is not increasing with position - exclude
    input:
        "/global/scratch/users/moorjani/denisova/data/asc_Vindija/big2.count.out"
    output:
        ND_bed="data/hg19/snp_annotations/ND_hg19.bed",
        ND="data/hg19/snp_annotations/ND_hg19.txt",
        ND11="data/hg19/snp_annotations/ND11_hg19.txt",
        ND11_bed="data/hg19/snp_annotations/ND11_hg19.bed"
    shell:
        """
        tail -n +2 {input} | awk '{{print $2 "\t" $4 "\t" $4}}' > {output.ND_bed}
        tail -n +2 {input} | awk '{{ print $2 "\t" $4 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13}}' > {output.ND}
        awk '$13 == "nd01" || $13 == "nd10" || $13 == "nd11"' {input} | grep -v -E "69355282|69362336|69436333|69495068|234475007|69250201|69250566|69566677|69579433|69579696" | awk '{{ print $2 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13}}'  > {output.ND11}
        awk '{{print $1 "\t" $2 "\t" $2}}' {output.ND11} > {output.ND11_bed}
        """

# filter to those in the pilot mask (these are the ones we can reliably use)
# to do this, I need to convert the list of positions into a bed file (chrom, start, end (= start))
# I finally sort the metadata file to only include those positions
rule mask_ND_pos:
    input:
        pos_bed="data/hg19/snp_annotations/ND11_hg19.bed",
        pos="data/hg19/snp_annotations/ND11_hg19.txt",
        mask="/global/scratch/p2p3/pl1_moorjani/sarahj32/data/resources/hg19/1000G/callability_masks/20141020.pilot_mask.whole_genome_noChr.bed"
    output:
        masked="data/hg19/snp_annotations/ND11_hg19_masked.txt",
        masked_bed="data/hg19/snp_annotations/ND11_hg19_masked.bed"
    shell:
        """
        bedtools intersect -a {input.pos_bed} -b {input.mask} > {output.masked_bed}
        awk 'FNR==NR {{ keep["X:"$1"_"$2]; next }}  (("X:"$1"_"$2) in keep) {{print}}' {output.masked_bed} {input.pos} > {output.masked}
        """
    

rule liftover_ND_to_hg38:
    # extract just the chrom and position for all ND positions
    # then extract nd01, nd10 and nd11 positions
    # there are several of these positions where the genetic map is not increasing with position - exclude
    input:
        anno="data/hg19/snp_annotations/ND11_hg19_masked.txt",
        chain="/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/resources/liftover.chainfiles/hg38ToHg19.over.chain"
    output:
        to_liftover=temp("data/hg19/snp_annotations/ND11_hg19_masked_wID.bed"),
        hg38_raw=temp("data/hg38/snp_annotations/ND11_hg38_masked_raw.bed"),
        hg38="data/hg38/snp_annotations/ND11_hg38_masked.bed"
    shell:
        """
        awk '{{print $1 "\t" 2 "\t" $2 "\t" $1 ":" $2 }}' {input.anno} > {output.to_liftover}
        CrossMap bed {input.chain} {output.to_liftover} {output.hg38_raw}
        awk 'FNR==NR {{ keep[$3]; next }}  (($1":"$2) in keep) {{print keep[$3] "\t" $0}}' {output.hg38_raw} {input.anno} > {output.hg38}
        """


# rule make_final_ascertainments:
#     input:
#         "data/hg19/snp_annotations/ND11_hg19_masked.txt"
#     output:
#         ND01_afrAnc="data/hg19/snp_annotations/ND01AfrAnc_freqs_hg19.txt",
#         ND10_afrAnc="data/hg19/snp_annotations/ND10AfrAnc_freqs_hg19.txt",
#         ND01="data/hg19/snp_annotations/ND01_freqs_hg19.txt",
#         ND10="data/hg19/snp_annotations/ND10_freqs_hg19.txt",
#         ND11_afrAnc="data/hg19/snp_annotations/ND11_freqs_hg19.txt",
#     shell:
#         """
#         awk '$9 == 0.0000 && $11 == "nd10"' {input} | awk '{{print "chr" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $7 "\t" $8 "\t" $9 "\t" $11}}' > {output.ND10_afrAnc}
#         awk '$9 == 0.0000 && $11 == "nd01"' {input} | awk '{{print "chr" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $7 "\t" $8 "\t" $9 "\t" $11}}' > {output.ND01_afrAnc}
#         awk '$11 == "nd10"' {input} | awk '{{print "chr" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $7 "\t" $8 "\t" $9 "\t" $11}}'  > {output.ND10}
#         awk '$11 == "nd01"' {input} | awk '{{print "chr" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $7 "\t" $8 "\t" $9 "\t" $11}}'  > {output.ND01}
#         awk '$11 == "nd10" || $11 == "nd11"' {input} | awk '{{print "chr" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $7 "\t" $8 "\t" $9 "\t" $11}}'  > {output.ND11}
#         """


# this one below only gets positions, but we can include annotations
rule make_final_ascertainments:
    input:
        "data/hg19/snp_annotations/ND11_hg19_masked.txt"
    output:
        ND01_afrAnc="data/hg19/snp_annotations/ND01AfrAnc_hg19.txt",
        ND10_afrAnc="data/hg19/snp_annotations/ND10AfrAnc_hg19.txt",
        ND11_afrAnc="data/hg19/snp_annotations/ND11AfrAnc_hg19.txt",
        ND01="data/hg19/snp_annotations/ND01_hg19.txt",
        ND10="data/hg19/snp_annotations/ND10_hg19.txt",
        ND11="data/hg19/snp_annotations/ND11Only_hg19.txt"
    shell:
        """
        awk '$9 == 0.0000 && $11 == "nd10"' {input} | awk '{{print "chr" $1 "\t" $2}}' > {output.ND10_afrAnc}
        awk '$9 == 0.0000 && $11 == "nd01"' {input} | awk '{{print "chr" $1 "\t" $2}}' > {output.ND01_afrAnc}
        awk '$9 == 0.0000 && $11 == "nd11"' {input} | awk '{{print "chr" $1 "\t" $2}}' > {output.ND11_afrAnc}
        awk '$11 == "nd10"' {input} | awk '{{print "chr" $1 "\t" $2}}'  > {output.ND10}
        awk '$11 == "nd01"' {input} | awk '{{print "chr" $1 "\t" $2}}'  > {output.ND01}
        awk '$11 == "nd11"' {input} | awk '{{print "chr" $1 "\t" $2}}'  > {output.ND11}
        """

# rule filter_1000G_to_nd:
#     # goal: snps, biallelic, ND pos
#     # neither of these have "chr" infront of the chromosome
#     input:
#         vcf=lambda wildcards: hg19_vcf_path.format(chromosome=wildcards.chrom),
#         pos="data/hg19/snp_annotations/ND11_hg19_noChr.txt"
#     output:
#         "data/{dataset}/ND_VCF/{dataset}_ND11pos_hg19_chr{chrom}.vcf.gz"
#     shell:
#         "bcftools view {input.vcf} -T {input.pos} -v snps -M2 > {output}"

# rule fill_in_homRef:
#     input:
#         vcf="data/{dataset}/ND_VCF/{dataset}_ND11pos_hg19_chr{chrom}.vcf.gz",
#         ref="/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/resources/hg19/fasta/hg19.fa",
#         pos="data/hg19/snp_annotations/ND11_hg19_noChr.txt",
#         script="helper_scripts/add_homref_2vcf.py"
#     output:
#         "data/{dataset}/ND_VCF_refFill/{dataset}_ND11pos_hg19_refFill_chr{chrom}.vcf.gz"
#     shell:
#         "python {input.script} -v {input.vcf} -r {input.ref} -p {input.pos} -o {output}"

# rule mask_NDpos:
#     input:
#         vcf = "data/{dataset}/ND_VCF_refFill/{dataset}_ND11pos_hg19_refFill_chr{chrom}.vcf.gz",
#         mask = "/global/scratch/p2p3/pl1_moorjani/sarahj32/data/resources/hg19/1000G/callability_masks/20141020.pilot_mask.whole_genome_noChr.bed"
#     output:
#         "data/{dataset}/ND_VCF_refFill_Masked/{dataset}_ND11pos_hg19_refFill_chr{chrom}.vcf.gz",
#     shell:
#         "bcftools view -T {input.mask} {input.vcf} | bcftools view -Oz -o {output}"

# add on chimp to compare to priyas
rule add_chimp_to_ND:
    input:
        vcf = "data/{dataset}/ND_VCF_refFill/{dataset}_ND11pos_hg19_refFill_chr{chrom}.vcf.gz",
        script = "/global/scratch/p2p3/pl1_moorjani/sarahj32/superarchaic_introgression/superarchaic_analysis/helper_scripts/add_fast_to_vcf.py",
        reference = "/global/scratch/p2p3/pl1_moorjani/sarahj32/data/resources/hg19/fasta/hg19_noChr.fa",
        chimp_fasta = "/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/hg19/resources/chimp/MAF_convert/chimp_noChr_hg19.fa"
    params: 
        unzipped = temp("data/{dataset}/ND_VCF_Masked_chimp/{dataset}_ND11pos_hg19_refFill_chr{chrom}.vcf.gz")
    output:
        new_vcf = "data/{dataset}/ND_VCF_Masked_chimp/{dataset}_ND11pos_hg19_refFill_chimp_chr{chrom}.vcf.gz",
    shell:
        """
        bcftools view {input.vcf} -o {params.unzipped}
        python {input.script} -f {input.chimp_fasta} -r {input.reference} -v {params.unzipped} -n panTro5 -o {output.new_vcf}
        """

rule get_chimp_ND_freqs:
    input:
        vcf="data/{dataset}/ND_VCF_Masked_chimp/{dataset}_ND11pos_hg19_refFill_chimp_chr{chrom}.vcf.gz",
    output:
        "data/{dataset}/ND_VCF_Masked_chimp/{dataset}_ND11pos_hg19_chimp_chr{chrom}_freqs.txt",
    shell:
        "bcftools view -s panTro5 {input.vcf} | bcftools +fill-tags -Ou -- -t AC,AN,AF | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AC\t%AN\t%AF\n' > {output}" 


rule remove_chr:
    input:
        pos = "data/hg19/snp_annotations/{asc}_hg19.txt"
    output:
        pos = "data/hg19/snp_annotations/{asc}_hg19_noChr.txt"
    shell:
        "sed 's/chr//g' {input.pos} > {output.pos}"

# rule get_all_geno:
#     input:
#         vcf = "/global/scratch/p2p3/pl1_moorjani/sarahj32/superarchaic_introgression/superarchaic_analysis/data/Archaics_1000G/AltDenVinPan_1000G_merged_strictArchaic_refFilledIn_hg19_chr{chromosome}.vcf.gz",
#         pos = "data/hg19/snp_annotations/{asc}_hg19_noChr.txt",
#         inds = "data/{dataset}/inds/{popu}.inds"
#     output:
#         "data/{dataset}/{popu}/{kind}/vcf/hg19/filledIn/{dataset}_snps10_08_chr{chromosome}_{kind}_{popu}_{asc}.vcf.gz"
#     shell:
#         "bcftools view {input.vcf} -T {input.pos} -S {input.inds} | bcftools +fill-tags -Ou -- -t AC,AN,AF | bcftools view -Oz -o {output}"

rule concate_chr_all_geno:
    input:
        expand("data/{{dataset}}/{{popu}}/{{kind}}/vcf/hg19/filledIn/{{dataset}}_snps10_08_chr{chromosome}_{{kind}}_{{popu}}_{{asc}}.vcf.gz", chromosome = list(range(1,23)))
    output:
        "data/{dataset}/{popu}/{kind}/vcf/hg19/filledIn/{dataset}_snps10_08_allChr_{kind}_{popu}_{asc}.vcf.gz"
    shell:
        "bcftools concat {input} | bcftools sort -o {output}"

# rule get_all_geno_freqs:
#     input:
#         "data/{dataset}/{popu}/{kind}/vcf/hg19/filledIn/{dataset}_snps10_08_allChr_{kind}_{popu}_{asc}.vcf.gz"
#     output:
#         "data/{dataset}/{popu}/{kind}/vcf/hg19/filledIn/{dataset}_snps10_08_allChr_{kind}_{popu}_{asc}.summary.txt"
#     shell:
#         "bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AC\t%AN\t%AF\n' {input} > {output}"
