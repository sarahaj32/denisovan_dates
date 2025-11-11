# a set of rules for finding the intersect of manifesto mask and 1000 G strict mask
# finding all archaic and ancestral (chimp and homo sapien ancestor) genotypes into 1 vcf
# and labelling every position

rule intersect_hg38_masks:
    input:
        mask_1000G = "/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/hg38/1000G/callability_mask/strict_mask.whole_genome.bed",
        archaic_mask = "/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/hg38/ANCIENT/ARCHAIC/manifesto_hg38/chr{chromosome}.bed"
    output:
        "data/resources/hg38/mask_intersect/strict_callability_archaic_manifesto_chr{chromosome}.bed"
    shell:
        """
        grep chr{wildcards.chromosome} {input.mask_1000G} | \
        bedtools intersect -a - -b {input.archaic_mask} > {output}
        """

# rule intersect_hg19_masks:
#     input:
#         mask_1000G="/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/hg19/1000G/callability_mask/20141020.strict_mask.whole_genome.bed",
#         archaic_mask="/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/hg19/ANCIENT/ARCHAIC/full_intersected_hg19_archaics.bed"
#     output:
#         "data/resources/hg19/mask_intersect/strict_callability_archaic_manifesto.bed"
#     shell:
#         """
#         bedtools intersect -a {input.mask_1000G} -b {input.archaic_mask} > {output}
#         """

rule mask_hg38_archaics:
    input:
        archaics = "/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/hg38/ANCIENT/ARCHAIC/individuals_highcov.{chromosome}.bcf",
        mask = "data/resources/hg38/mask_intersect/strict_callability_archaic_manifesto_chr{chromosome}.bed",
    output:
        vcf = "data/Archaics/hg38/AltDenVin_merged_strictArchaic_hg38_chr{chromosome}.vcf.gz",
        index = "data/Archaics/hg38/AltDenVin_merged_strictArchaic_hg38_chr{chromosome}.vcf.gz.csi"
    resources:
        time="200:00:00"
    shell:
    # merge with the info of chimp and AltaiNeandertal
        """
        bcftools view -g ^miss -R {input.mask} {input.archaics} -Ob -o {output.vcf}
        bcftools index {output.vcf}
        """

rule mask_hg19_archaics:
    input:
        archaics = "/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/hg19/ANCIENT/ARCHAIC/individuals_highcov.{chromosome}.bcf",
        mask = "data/resources/hg19/mask_intersect/strict_callability_archaic_manifesto_chr{chromosome}.bed",
    output:
        vcf = "data/Archaics/hg19/AltDenVin_merged_strictArchaic_hg19_chr{chromosome}.vcf.gz",
        index = "data/Archaics/hg19/AltDenVin_merged_strictArchaic_hg19_chr{chromosome}.vcf.gz.csi"
    resources:
        time="200:00:00"
    shell:
    # merge with the info of chimp and AltaiNeandertal
        """
        bcftools view -g ^miss -R {input.mask} {input.archaics} -Ob -o {output.vcf}
        bcftools index {output.vcf}
        """

rule append_chimp_fasta:
    input:
        vcf = "data/Archaics/hg38/AltDenVin_merged_strictArchaic_hg38_chr{chromosome}.vcf.gz",
        chimp_fasta = "/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/hg38/resources/chimp/chimp_hg38.chr{chromosome}.fa",
        reference = "/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/hg38/resources/fasta/hg38.fa",
        script = "/global/scratch/p2p3/pl1_moorjani/sarahj32/superarchaic_introgression/superarchaic_analysis/helper_scripts/add_fast_to_vcf.py"
    output:
        vcf = temp("data/Archaics/hg38/AltDenVin_merged_strictArchaic_hg38_chr{chromosome}.vcf"),
        new_vcf = "data/Archaics/hg38/AltDenVinPan_merged_strictArchaic_hg38_chr{chromosome}.vcf",
    
    shell:
        """
        bcftools view {input.vcf} -o {output.vcf}
        python {input.script} -f {input.chimp_fasta} -r {input.reference} -v {output.vcf} -n panTro5 -o {output.new_vcf}
        """

rule correct_chr_name_ancestor:
    input:
        anc_fasta = "/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/hg38/resources/human_ancestor/homo_sapiens_ancestor_{chromosome}.fa"
    output:
        fasta="data/resources/hg38/homo_sapiens_ancestor_reheaded_chr{chromosome}.fa"
    shell:
        "sed -E 's/^>ANCESTOR_for_chromosome:GRCh38:/>chr/' {input.anc_fasta} | sed -E 's/:.*$//' > {output.fasta}"

rule append_ancester_fasta:
    input:
        vcf="data/Archaics/hg38/AltDenVinPan_merged_strictArchaic_hg38_chr{chromosome}.vcf",
        anc_fasta="data/resources/hg38/homo_sapiens_ancestor_reheaded_chr{chromosome}.fa",
        # is this the correct ancestor to be using???????
        reference="/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/hg38/resources/fasta/hg38.fa",
        script="/global/scratch/p2p3/pl1_moorjani/sarahj32/superarchaic_introgression/superarchaic_analysis/helper_scripts/add_fast_to_vcf.py"
    output:
        vcf = "data/Archaics/hg38/AltDenVinPanHSanc_merged_strictArchaic_hg38_chr{chromosome}.vcf",
    shell:
        """
        python {input.script} -f {input.anc_fasta} -r {input.reference} -v {input.vcf} -n homo_sapien_ancestor -o {output.vcf}
        """

rule liftover_archaics_hg38_to_hg19:
    input:
        vcf="data/Archaics/hg38/AltDenVinPanHSanc_merged_strictArchaic_hg38_chr{chromosome}.vcf"
    output:
        vcf="data/Archaics/hg19/AltDenVinPanHSanc_merged_strictArchaic_hg19_chr{chromosome}.vcf",
    params:
        chainfile="/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/resources/liftover.chainfiles/hg38ToHg19.over.chain",
        ref="/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/resources/hg19/fasta/hg19.fa",
    resources:
        cpus='6'
    shell:
        """
        CrossMap vcf {params.chainfile} {input.vcf} {params.ref} {output.vcf} --no-comp-alleles
        """
        
        # compare the 2 human ancestor files (should be the same but maybe not??) 
        # THEY'RE NOT - NEED TO FIGURE OUT WHY!!!!!!!!!
rule correct_chr_name_ancestor_S:
    input:
        anc_fasta = "/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/resources/hg38/human_ancestor/homo_sapiens_ancestor_{chromosome}.fa"
    output:
        fasta="data/resources/hg38/homo_sapiens_ancestor_LS_reheaded_chr{chromosome}.fa"
    shell:
        "sed -E 's/^>ANCESTOR_for_chromosome:GRCh38:/>chr/' {input.anc_fasta} | sed -E 's/:.*$//' > {output.fasta}"

rule append_ancester_fasta_LS:
    input:
        vcf="data/Archaics/hg38/AltDenVinPanHSanc_merged_strictArchaic_hg38_chr{chromosome}.vcf",
        anc_fasta="data/resources/hg38/homo_sapiens_ancestor_LS_reheaded_chr{chromosome}.fa",
        reference="/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/hg38/resources/fasta/hg38.fa",
        script="/global/scratch/p2p3/pl1_moorjani/sarahj32/superarchaic_introgression/superarchaic_analysis/helper_scripts/add_fast_to_vcf.py"
    output:
        vcf = "data/Archaics/hg38/AltDenVinPanHSanc_LS_merged_strictArchaic_hg38_chr{chromosome}.vcf",
    shell:
        """
        python {input.script} -f {input.anc_fasta} -r {input.reference} -v {input.vcf} -n homo_sapien_ancestor -o {output.vcf}
        """
#### #HMMMM These are not the same?!?! 


rule annotate_archaic_variants_chimp:
    input:
        archaics="data/Archaics/hg38/AltDenVinPanHSanc_merged_strictArchaic_hg38_chr{chromosome}.vcf",
        script="add_annotations_to_archaic_observations.py"
    output:
        annotated="data/Archaics/hg38/AltDenVinPanHSanc_merged_strictArchaic_hg38_chr{chromosome}_annotated_chimpAnc.txt",
    shell:
        "python {input.script} -i {input.archaics} -o {output.annotated} -a panTro5"

rule annotate_archaic_variants_hsAnc:
    input:
        archaics="data/Archaics/hg{assembly}/AltDenVinPanHSanc_merged_strictArchaic_hg{assembly}_chr{chromosome}.vcf",
        script="add_annotations_to_archaic_observations.py"
    output:
        annotated="data/Archaics/hg{assembly}/AltDenVinPanHSanc_merged_strictArchaic_hg{assembly}_chr{chromosome}_annotated_hsAnc.txt",
    shell:
        "python {input.script} -i {input.archaics} -o {output.annotated} -a ancestral"

rule combine_annotated_archaics:
# we need to update everything to here for hg19
    input:
        expand("data/Archaics/hg{{assembly}}/AltDenVinPanHSanc_merged_strictArchaic_hg{{assembly}}_chr{chromosome}_annotated_{{anc}}Anc.txt", chromosome = list(range(1,23)))
    output:
        "data/Archaics/hg{assembly}/AltDenVinPanHSanc_merged_strictArchaic_hg{assembly}_full_annotated_{anc}Anc.txt"
         
    shell:
        """
        head -n 1 {input[0]} > {output}
        for file in {input};
        do
            tail -n +2 "$file" >> {output}
        done
        """

rule get_yri_inds:
    input:
        meta="/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/hg38/1000G/1000G_NYGC/1000G.ind"
    output:
        inds="data/YRI/YRI_outgroup.ind"
    shell:
        """
        grep -E "YRI" {input.meta} | cut -f1 > {output.inds}
        """

rule get_yri_freqs:
    input:
        vcf = "/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/hg38/1000G/1000G_NYGC/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr{chrom}.recalibrated_variants.vcf.gz",
        inds = "data/YRI/YRI_outgroup.ind"
    output:
        freqs = "data/YRI/YRI_frequencies_hg38_chr{chrom}.txt"
    resources:
        time="5:00:00"
    shell:
        """
        bcftools view {input.vcf} -S {input.inds} -v snps -M2 | bcftools +fill-tags -Ou -- -t AC,AN,AF | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AC\t%AN\t%AF\n' > {output.freqs}
        """

