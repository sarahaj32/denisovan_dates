rule merge_archaics: # (need to keep these unzipped so I can add on chimp manually?)
    input:
        archaics = expand("/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/hg19/ANCIENT/ARCHAIC/{pop}/chr{{chromosome}}_mq25_mapab100.vcf.gz", pop = ["Altai", "Vindija33.19", "Denisova"]),
        mask = "data/ArtificialGenome/redo/hg19/mask_intersect/strict_callability_archaic_manifesto_chr{chromosome}.bed",
    output:
        vcf = "data/Archaics/AltDenVin_merged_strictArchaic_hg19_chr{chromosome}.vcf.gz",
        index = "data/Archaics/AltDenVin_merged_strictArchaic_hg19_chr{chromosome}.vcf.gz.csi"
    resources:
        time="200:00:00"
    shell:
    # merge with the info of chimp and AltaiNeandertal
        """
        bcftools merge {input.archaics} -R {input.mask} | bcftools annotate -x ^FORMAT/GT -Ob -o {output.vcf}
        bcftools index {output.vcf}
        """

rule append_chimp_fasta:
# I confirm that this removes positions unmasked. However, it also includes positions that are not in Leo's file
# I am trying to investigate what alleles we don't have for modern human:
# file path: /global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/hg19/1000G
# there are some fixed ancestral in here, but they all have names so probably not all of them?:
# bcftools view -C0 ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | less -S
    input:
        vcf = "data/Archaics/AltDenVin_merged_strictArchaic_hg19_chr{chromosome}.vcf.gz",
        chimp_fasta = "/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/hg19/resources/chimp/MAF_convert/chimp_hg19.chr{chromosome}.fa",
        reference = "/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/hg19/resources/fasta/hg19.fa",
        script = "helper_scripts/add_fast_to_vcf.py"
    output:
        vcf = temp("data/Archaics/AltDenVin_merged_strictArchaic_hg19_chr{chromosome}.vcf"),
        new_vcf = "data/Archaics/AltDenVinPan_merged_strictArchaic_hg19_chr{chromosome}.vcf",
        #new_gvcf = "data/Archaics/AltDenVinPan_strictArchaic_merged_hg19_chr{chromosome}.vcf.gz"
    shell:
        """
        bcftools view {input.vcf} -o {output.vcf}
        python {input.script} -f {input.chimp_fasta} -r {input.reference} -v {output.vcf} -n panTro5 -o {output.new_vcf}
        """

rule get_archaic_chimp_variable:
    input: 
        "data/Archaics/AltDenVinPan_merged_strictArchaic_hg19_chr{chromosome}.vcf"
    output:
        "data/Archaics/AltDenVinPan_merged_strictArchaic_hg19_variable_chr{chromosome}.vcf.gz"
    shell:
        """
        bcftools view -c 1 {input} -Oz -o {output}
        bcftools index {output}
        """

rule get_archaic_chimp_pos:
    # these are the positions within the strict mask and archaic manifesto where panTro5 is sequenced
    input: 
        "data/Archaics/AltDenVinPan_merged_strictArchaic_hg19_chr{chromosome}.vcf"
    output:
        "data/Archaics/AltDenVinPan_merged_strictArchaicChimp_hg19_chr{chromosome}.txt"
    shell:
        """
        bcftools view {input} -g ^miss | grep -v ^# | cut -f1,2 > {output}
        """

rule add_ref_1000G: # (need to keep these unzipped so I can add on chimp manually?)
    input:
        modern = "/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/hg19/1000G/ALL.chr{chromosome}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz", 
        pos = "data/Archaics/AltDenVinPan_merged_strictArchaicChimp_hg19_chr{chromosome}.txt",
        archaics = "data/Archaics/AltDenVinPan_merged_strictArchaic_hg19_variable_chr{chromosome}.vcf.gz"
    output:
        merged = "data/Archaics_1000G/AltDenVinPan_1000G_merged_strictArchaic_refFilledIn_hg19_chr{chromosome}.vcf.gz"
    resources:
        time="190:00:00"
    shell:
    # merge with the info of chimp and AltaiNeandertal
    # I confirmed that this does fill in the reference in both merged files
        """
        bcftools merge {input.modern} {input.archaics} -R {input.pos} -0 -Oz -o {output.merged}
        bcftools index {output.merged}
        """

rule add_chimp2:
    input:
        script="helper_scripts/add_fast_to_vcf.py",
        chimp_fa="/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/hg19/resources/chimp/MAF_convert/chimp_hg19.chr{chromosome}.fa",
        reference_fa="/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/hg19/resources/fasta/hg19.fa",
        vcf="data/Archaics_1000G/AltDenVin_1000G_merged_strictArchaic_refFilledIn_hg19_chr{chromosome}.vcf.gz"
    output:
        unzipped=temp("data/Archaics_1000G/AltDenVin_1000G_merged_strictArchaic_refFilledIn_hg19_chr{chromosome}.vcf"),
        vcf_wMissing="data/ArtificialGenome/redo/hg19/archaics/chimp/AltaiNeandertal_chimp_chr{chromosome}_wMissing2.vcf",
    shell:       
        """
        bcftools view {input.vcf} -o {output.unzipped}
        echo "unzipped"
        python {input.script} -v {output.unzipped} -f {input.chimp_fa} -n chimp -r {input.reference_fa} -o {output.vcf_wMissing}
        echo "python ran"
        """
