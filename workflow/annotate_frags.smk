# a set of rules for adding archaic annotations to inferred introgressed fragments
# also separates fragments by population and creates syntethetic genomes

# rule decode_GAv1_hmmix:
#     input:
#         mutationrate="/global/scratch/users/zhangyulin9806/github/ArchaicMutRate/hmmix_output/GA100K_raw/mutationrate.bed",
#         ind_json="/global/scratch/users/zhangyulin9806/github/ArchaicMutRate/hmmix_output/GA100K_raw/GA100K.samples.json",
#         # YULIN CLEARLY USED THE STRICT MASK HERE, I CAN REPEAT HER RESULTS IF I ALSO USE IT!!!!!
#         #weights="/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/resources/hg38/repeat_mask_files/hg38_all.bed",
#         weights_strict="/global/scratch/users/zhangyulin9806/github/ArchaicMutRate/hmmix_output/GA100K_raw/hg38_strictmask.bed",
#         obs="/global/scratch/p2p3/pl1_moorjani/sarahj32/GA100K_hmmix/00obs/obs.{ind}.txt",
#         archaic_bcf=expand("/global/scratch/p2p3/pl1_moorjani/lauritsskov2/ArchaicSegments/helperfiles/archaic_variants/hg38/individuals_highcov.{chrom_id}.bcf", chrom_id=range(1, 23)),
#         hmmix_trained="/global/scratch/users/zhangyulin9806/github/ArchaicMutRate/hmmix_output/GA100K_raw/01train_haploid/trained.{ind}.json",
#         hmmix_script="/global/scratch/users/sarahj32/software/miniconda3/envs/r_python/bin/hmmix"
#     params:
#         hmmix_out=lambda wildcards: f"data/GA100K/hmmix_strict/02decode/{wildcards.ind}",
#         # hmmix_out=lambda wildcards: f"data/GA100K/hmmix/02decode/{wildcards.ind}",
#         archaic_bcf="/global/scratch/p2p3/pl1_moorjani/lauritsskov2/ArchaicSegments/helperfiles/archaic_variants/hg38/",
#         admixpop="/global/scratch/p2p3/pl1_moorjani/lauritsskov2/ArchaicSegments/helperfiles/archaic_variants/hg38/individuals_highcov.*.bcf"
#     output:
#         #hmmix_out1="data/GA100K/hmmix/02decode/{ind}.hap1.txt",
#         #hmmix_out2="data/GA100K/hmmix/02decode/{ind}.hap2.txt",
#         hmmix_out1="data/GA100K/hmmix_strict/02decode/{ind}.hap1.txt",
#         hmmix_out2="data/GA100K/hmmix_strict/02decode/{ind}.hap2.txt",
#     shell:
#         """
#         {input.hmmix_script} decode -obs={input.obs} -weights={input.weights_strict} -mutrates={input.mutationrate} -param={input.hmmix_trained} -admixpop={params.admixpop} -haploid -out={params.hmmix_out} -extrainfo
#         """

# input requires a bed or text file with all fragments, filtered to 
rule annotate_fragments:
    input:
        bed = "data/{dataset}/{dataset}_fragments_snps10_08.txt",
        annotations = "data/Archaics/hg38/AltDenVinPanHSanc_merged_strictArchaic_hg38_full_annotated_hsAnc.txt",
        script = "helper_scripts/annotate_frags_archaics.py"
    output:
        "data/{dataset}/{dataset}_fragments_snps10_08_annotated.txt"
    resources:
        time="4:00:00",
        cpus="20"
    shell:
        "python {input.script} -b {input.bed} -i {input.annotations} -o {output} -a ancestral"

rule annotate_all_fragments:
    input:
        bed = "data/{dataset}/{dataset}_fragments_08.txt",
        annotations = "data/Archaics/hg38/AltDenVinPanHSanc_merged_strictArchaic_hg38_full_annotated_hsAnc.txt",
        script = "helper_scripts/annotate_frags_archaics.py"
    output:
        "data/{dataset}/{dataset}_fragments_08_annotated.txt"
    resources:
        time="4:00:00",
        cpus="20"
    shell:
        "python {input.script} -b {input.bed} -i {input.annotations} -o {output} -a ancestral"

# this assumes my input files have been ran exactly through this pipeline and have the columns
# chr, start, end = $5, $6, $7

rule convert_fragments_to_bed:
    input:
        frags = "data/{dataset}/{dataset}_fragments_08_annotated.txt",
    output:
        unsorted = temp("data/{dataset}/{dataset}_fragments_08_annotated_UNSORTED.bed"),
        bed = "data/{dataset}/{dataset}_fragments_08_annotated.bed"
    shell:
        """
        awk '{{ 
            # print columns 5, 6, 7
            printf "%s\\t%s\\t%s", $5, $6, $7;

            # print columns 1 to 4
            for (i = 1; i <= 4; i++) {{
                printf "\\t%s", $i;
            }}

            # print columns 8 to end
            for (i = 8; i <= NF; i++) {{
                printf "\\t%s", $i;
            }}

            printf "\\n";
        }}' {input.frags} > {output.unsorted}

        head -n 1 {output.unsorted} > {output.bed}
        tail -n +2 {output.unsorted} | sort -k1,1 -k2,2n >> {output.bed}
        """

# get population fragments:
# get the index lists myself (create_1000G_inds.sh or create_LASIDAD_inds.sh)
# I'll do this on my own so i have the most flexibility
# output: data/1000G_${popu}/{popu}_fragments_snps10_08_nonoverlapping.bed
# get fragments and create non-overlapping from a list of individuals
rule filter_pop_fragments:
    input:
        frags = "data/{dataset}/{dataset}_{details}.bed",
        inds = "data/{dataset}/inds/{popu}.inds"
    output:
        frags = "data/{dataset}/{popu}/fragments/{popu}_{details}.bed" # fragments_snps10_08_annotated
    shell:
        #"grep -F -f {input.inds} {input.frags} > {output.frags}"
        """
        grep -E 'chrom|'"$(paste -sd'|' {input.inds})" {input.frags} > {output.frags}
        """

rule classify_NeaDen_MostLikely:
    input:
        bed = "data/{dataset}/{popu}/fragments/{popu}_{details}.bed",
    output:
        nea_full = "data/{dataset}/{popu}/fragments/NEA/{popu}_{details}_mostLikely.bed", # fragments_snps10_08_annotated
        den_full = "data/{dataset}/{popu}/fragments/DEN/{popu}_{details}_mostLikely.bed"
    shell:
        #"grep -F -f {input.inds} {input.frags} > {output.frags}"
        """
        nea=$(head -1 {input.bed} | tr '\t' '\n' | grep -n -x "nea_overlap" | cut -d: -f1)
        den=$(head -1 {input.bed} | tr '\t' '\n' | grep -n -x "den_overlap" | cut -d: -f1)

        awk -v d="$den" -v n="$nea" 'NR==1 || ($d > $n)' {input.bed} > {output.den_full}
        awk -v d="$den" -v n="$nea" 'NR==1 || ($d < $n)' {input.bed} > {output.nea_full}
        """

rule make_pop_synthetic:
    input:
        frags = "data/{dataset}/{popu}/fragments/{popu}_fragments_snps10_08_annotated.bed",
        script = "/global/scratch/p2p3/pl1_moorjani/sarahj32/superarchaic_introgression/super_clean/helper_scripts/find_nonoverlapping.py"
    output:
        synth = "data/{dataset}/{popu}/fragments/{popu}_fragments_snps10_08_nonoverlapping_full.bed"
    shell:
        "python {input.script} -b {input.frags} -o {output.synth}"

rule classify_fragments:
    input:
        bed = "data/{dataset}/{popu}/fragments/{popu}_fragments_snps10_08_nonoverlapping_full.bed",
    output:
        high="data/{dataset}/{popu}/fragments/{popu}_fragments_snps10_08_nonoverlapping_Dhigh.txt", 
        high_bed="data/{dataset}/{popu}/fragments/{popu}_fragments_snps10_08_nonoverlapping_Dhigh.bed", 
        low="data/{dataset}/{popu}/fragments/{popu}_fragments_snps10_08_nonoverlapping_Dlow.txt",
        low_bed="data/{dataset}/{popu}/fragments/{popu}_fragments_snps10_08_nonoverlapping_Dlow.bed", 
        nea="data/{dataset}/{popu}/fragments/{popu}_fragments_snps10_08_nonoverlapping_Nea.txt",
        nea_bed="data/{dataset}/{popu}/fragments/{popu}_fragments_snps10_08_nonoverlapping_Nea.bed"
    shell:
        """
        awk '$1 == "chrom" || ($27 >= 0.7 && $26 <= 0.3)' {input.bed} > {output.high}
        cut -f1,2,3 {output.high} > {output.high_bed}
        awk '$1 == "chrom" || ($27 < 0.7 && $27 >= 0.4 && $26 <= 0.3)' {input.bed} > {output.low}
        cut -f1,2,3 {output.low} > {output.low_bed}
        awk '$1 == "chrom" || ($27 <= 0.3 && $26 >= 0.5)' {input.bed} > {output.nea}
        cut -f1,2,3 {output.nea} > {output.nea_bed}
        """
#
rule strict_classify_fragments:
    input:
        bed = "data/{dataset}/{popu}/fragments/{popu}_fragments_snps10_08_nonoverlapping_full.bed",
    output:
        high="data/{dataset}/{popu}/fragments/{popu}_fragments_snps10_08_nonoverlapping_Dhigh08.txt", 
        high_bed="data/{dataset}/{popu}/fragments/{popu}_fragments_snps10_08_nonoverlapping_Dhigh08.bed", 
        low="data/{dataset}/{popu}/fragments/{popu}_fragments_snps10_08_nonoverlapping_Dlow38.txt",
        low_bed="data/{dataset}/{popu}/fragments/{popu}_fragments_snps10_08_nonoverlapping_Dlow38.bed", 
        nea="data/{dataset}/{popu}/fragments/{popu}_fragments_snps10_08_nonoverlapping_Nea06.txt",
        nea_bed="data/{dataset}/{popu}/fragments/{popu}_fragments_snps10_08_nonoverlapping_Nea06.bed"
    shell:
        """
        awk '$1 == "chrom" || ($27 >= 0.8 && $26 <= 0.3)' {input.bed} > {output.high}
        cut -f1,2,3 {output.high} > {output.high_bed}
        awk '$1 == "chrom" || ($27 < 0.8 && $27 >= 0.3 && $26 <= 0.3)' {input.bed} > {output.low}
        cut -f1,2,3 {output.low} > {output.low_bed}
        awk '$1 == "chrom" || ($27 <= 0.3 && $26 >= 0.6)' {input.bed} > {output.nea}
        cut -f1,2,3 {output.nea} > {output.nea_bed}
        """

# here we'll need to work harder since we'll need all the hg38 paths BUT this is only for dating so I think we're good :)
# lots we can do until this point!!!
rule filter_dataset_to_fragments:
    input:
        vcf=lambda wildcards: hg38_vcf_path.format(chromosome=wildcards.chrom),
        inds="data/{dataset}/inds/{popu}.inds",
        bed="data/{dataset}/{popu}/fragments/{popu}_fragments_snps10_08_nonoverlapping_{kind}.bed"
    output:
        vcf="data/{dataset}/{popu}/{kind}/vcf/hg38/{dataset}_snps10_08_{kind}_{popu}_chr{chrom}.vcf.gz"
    shell:
        """
        bcftools view -S {input.inds} -R <(tail -n +2 {input.bed}) -M2 -v snps {input.vcf} | bcftools +fill-tags -Ou -- -t AC,AN,AF | bcftools annotate -x ^FORMAT/GT -Oz -o {output.vcf}
        bcftools index {output.vcf}
        """

rule combine_chrs:
    input: 
        vcf=expand("data/{{dataset}}/{{popu}}/{{kind}}/vcf/hg38/{{dataset}}_snps10_08_{{kind}}_{{popu}}_chr{chrom}.vcf.gz", chrom = list(range(1,23)))
    output:
        combined_vcf="data/{dataset}/{popu}/{kind}/vcf/hg38/{dataset}_snps10_08_allChr_{kind}_{popu}.vcf"
    shell:
        "bcftools concat {input.vcf} | bcftools sort -o {output.combined_vcf}"


rule prepare_to_liftover:
# combine all of the fragments
# rename the vcf to be consistent when we lift over
    input:
        vcf="data/{dataset}/{popu}/{kind}/vcf/hg38/{dataset}_snps10_08_allChr_{kind}_{popu}.vcf",
        chr_map="data/resources/chr_map"
    output:
        rename="data/{dataset}/{popu}/{kind}/vcf/hg38/{dataset}_snps10_08_allChr_{kind}_{popu}_chrRename.vcf"
    shell:
        """
        bcftools annotate --rename-chrs {input.chr_map} {input.vcf} > {output.rename}
        """

rule liftover_vcf_hg38_to_hg19:
    input:
        vcf="data/{dataset}/{popu}/{kind}/vcf/hg38/{dataset}_snps10_08_allChr_{kind}_{popu}_chrRename.vcf"
    output:
        vcf="data/{dataset}/{popu}/{kind}/vcf/hg19/{dataset}_snps10_08_allChr_{kind}_{popu}.vcf"
    params:
        chainfile="/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/resources/liftover.chainfiles/hg38ToHg19.over.chain",
        ref="/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/resources/hg19/fasta/hg19.fa",
    resources:
        cpus='6'
    shell:
        """
        CrossMap vcf {params.chainfile} {input.vcf} {params.ref} {output.vcf} --no-comp-alleles
        """

rule normalize:
    input:  
        vcf="data/{dataset}/{popu}/{kind}/vcf/hg19/{dataset}_snps10_08_allChr_{kind}_{popu}.vcf",
        ref="/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/hg19/resources/fasta/hg19.fa",
        chr_map="data/resources/chr_map_addChr"
    params:
        sorted="data/{dataset}/{popu}/{kind}/vcf/hg19/{dataset}_snps10_08_allChr_{kind}_{popu}_norm.vcf"
    output:
        compress="data/{dataset}/{popu}/{kind}/vcf/hg19/{dataset}_snps10_08_allChr_{kind}_{popu}_norm.vcf.gz",
        index="data/{dataset}/{popu}/{kind}/vcf/hg19/{dataset}_snps10_08_allChr_{kind}_{popu}_norm.vcf.gz.csi",
    shell:
        """
        bcftools sort {input.vcf} | bcftools norm -d any | bcftools annotate --rename-chrs {input.chr_map} -o {params.sorted} 
        $SCRIPTS/bgzip -f {params.sorted}
        bcftools index {output.compress}
        """

rule find_segregating:
    input:
        "data/{dataset}/{popu}/{kind}/vcf/hg19/{dataset}_snps10_08_allChr_{kind}_{popu}_norm.vcf.gz"
    output:
        "data/{dataset}/{popu}/{kind}/vcf/hg19/{dataset}_snps10_08_allChr_{kind}_{popu}_var.vcf.gz"
    resources:
        time="10:00:00"
    shell:
        "bcftools view -i 'INFO/AF>0 && INFO/AF<1' {input} -Oz -o {output}"

