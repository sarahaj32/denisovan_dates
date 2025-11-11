# remove the not-increasing positions from the genetic map
rule remove_problem_positions:
    input:
        genmap="/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/hg19/resources/genetic_maps_corrections/hg19/Shared_AfAmMap/genetic_map_chr{chrom}.txt",
        script="helper_scripts/find_continuous_positions.py"
    output:
        no_probs = "data/genmap/hg19/Shared_AfAmMap_continuous/Shared_AfAmMap_continuous_genmap_chr{chrom}.txt",
        probs = "data/genmap/hg19/Shared_AfAmMap_PROBLEMS/Shared_AfAmMap_PROBLEMS_genmap_chr{chrom}.txt",
    shell:
        "python {input.script} -s {input.genmap} -p {output.probs} -n {output.no_probs}"

# a rule to find all the variable positions within a population
rule find_segregating_modern:
    input:
        vcf=lambda wildcards: config["datasets"][wildcards.dataset]["hg38_path"].format(chromosome=wildcards.chrom),
        inds = "data/{dataset}/inds/{popu}.inds",
        mask="/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/hg38/1000G/callability_mask/strict_mask.whole_genome.bed"
    output:
        "data/{dataset}/{popu}/modern/vcf/hg38/{dataset}_{popu}_strictMask_var_chr{chrom}.vcf"
    shell:
        "bcftools view -S {input.inds} -T {input.mask} -M2 -v snps {input.vcf} | bcftools +fill-tags -Ou -- -t AC,AN,AF | bcftools view -i 'INFO/AF>0 && INFO/AF<1' | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AC\t%AN\t%AF\n' -o {output}"

rule downsample:
    input:
        vcf="data/{dataset}/{popu}/modern/vcf/hg38/{dataset}_{popu}_strictMask_var_chr{chrom}.vcf"
    output:
        vcf="data/{dataset}/{popu}/modern/vcf/hg38/{dataset}_{popu}_strictMask_var_Downsample{ds}_chr{chrom}.vcf"
    shell:
        """
        lines=$(wc -l {input.vcf} | awk '{{print $1}}')
        ds=$(($lines/{wildcards.ds}))
        echo -e "#CHROM\tPOS"> {output.vcf}
        cut -f1,2 {input.vcf} | shuf -n $ds | sort -nk2 >> {output.vcf}
        """

rule reformat_and_split_frags:
    input:
        "data/{dataset}/{popu}/fragments/{archaic}/{popu}_fragments_08_annotated_mostLikely.bed",
    output:
        "data/{dataset}/{popu}/fragments/{archaic}/{popu}_fragments_08_annotated_mostLikely_chr{chrom}.bed",
    shell:
        """
        header=$(head -n 1 {input})
        echo -e "hap\t$header" > {output}
        tail -n +2 {input} | \
        awk -v chrom=chr{wildcards.chrom} '$1 == chrom' | \
        sort -k4,4 -k1.4,1n -k5,5 -k2,2n | \
        awk 'BEGIN {{prev=""; id=""; n=-1}} {{if ($5 != prev || $4 != id) {{n++}} prev = $5; id = $4; print n "\t" $0}}' >> {output}
        """

rule create_ancestry_matrix:
    input:
        bed = "data/{dataset}/{popu}/fragments/{archaic}/{popu}_fragments_08_annotated_mostLikely_chr{chrom}.bed",
        vcf = "data/{dataset}/{popu}/modern/vcf/hg38/{dataset}_{popu}_strictMask_var_Downsample{ds}_chr{chrom}.vcf",
        script = "helper_scripts/create_ancestry_eig.py",
    output:
        "data/{dataset}/{popu}/hg38/{archaic}/curve/anc_matrix/{dataset}_{popu}_mostLikely_Downsample{ds}_chr{chrom}.anc"
    shell:
        """
        python {input.script} -b {input.bed} -i {input.vcf} -o {output}
        """

rule liftover_anc_matrix:
    input:
        anc = "data/{dataset}/{popu}/hg38/{archaic}/curve/anc_matrix/{dataset}_{popu}_mostLikely_Downsample{ds}_chr{chrom}.anc",
        chainfile="/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/resources/liftover.chainfiles/hg38ToHg19.over.chain",
    output:
        hg38_coords = "data/{dataset}/{popu}/hg38/{archaic}/curve/anc_matrix/{dataset}_{popu}_mostLikely_Downsample{ds}_chr{chrom}.anc.coords.bed",
        hg19_coords = "data/{dataset}/{popu}/hg19/{archaic}/curve/anc_matrix/{dataset}_{popu}_mostLikely_Downsample{ds}_chr{chrom}.anc.coords.bed",
        hg19_anc = "data/{dataset}/{popu}/hg19/{archaic}/curve/anc_matrix/{dataset}_{popu}_mostLikely_Downsample{ds}_chr{chrom}.anc",
        
    shell:
        """
        awk 'BEGIN{{OFS="\t"}} {{print $1, $2, $2, NR}}' {input.anc} > {output.hg38_coords}
        CrossMap bed {input.chainfile} {output.hg38_coords} {output.hg19_coords}
        awk 'BEGIN{{FS=OFS="\t"}} NR==FNR {{lift[$4]=$1"\t"$2; next}} {{line_num = FNR; if(line_num in lift) {{ $1=$2=""; sub(/^\t+/, ""); print lift[line_num], $0}}}}' {output.hg19_coords} {input.anc}  | grep -e "^{wildcards.chrom}" | sort -nk2 > {output.hg19_anc}
        """

# implement interpit myself - merge on the genetic map positions using the chromosome position
# then need to confirm that all positions fall within the map (they should if they're in the strict mask??)
# genmap: chrom position = 1, genetic position = 3
# anc: chrom position = 2, ancestry = 3
rule interpit_anc_ME:
    input:
        genmap="data/genmap/hg19/Shared_AfAmMap_continuous/Shared_AfAmMap_continuous_genmap_chr{chrom}.txt",
        anc="data/{dataset}/{popu}/hg19/{archaic}/curve/anc_matrix/{dataset}_{popu}_mostLikely_Downsample{ds}_chr{chrom}.anc",
    output:
        anc = "data/{dataset}/{popu}/hg19/{archaic}/curve/genmap/{dataset}_{popu}_mostLikely_Downsample{ds}_chr{chrom}_genmap_myattempt.anc",
    shell:
        """
        min=$(head -n 2 {input.genmap} | tail -n 1 | awk '{{print $2}}')
        max=$(tail -n 1 {input.genmap} | awk '{{print $2}}')
        awk 'NR==FNR {{id[$1]=$3; next}} ($2 in id) {{print $1 "\t" id[$2] "\t" $0}}' {input.genmap} {input.anc} | cut --complement -f3,4 | awk -v min=$min -v max=$max '{{$2 > $min && $2 < $max}}' > {output.anc}
        """
        # also filter to positions within the map
        # data/LASIDAD/South/hg19/DEN/curve/genmap/LASIDAD_South_mostLikely_Downsample6_chr3_genmap.anc
# rule create_snp_file:
#     input:
#         anc = "data/{dataset}/{popu}/hg19/{archaic}/curve/anc_matrix/{dataset}_{popu}_mostLikely_Downsample{ds}_chr{chrom}.anc",
#     output:
#         snp = "data/{dataset}/{popu}/hg19/{archaic}/curve/anc_matrix/{dataset}_{popu}_mostLikely_Downsample{ds}_chr{chrom}.snp",
#     shell:
#         """
#         awk '{{print $1 ":" $2 "\t" $1 "\t0.0\t" $2, "A\tT"}}' {input.anc} > {output.snp}
#         """

# rule interpit_snp:
#     input:
#         snp = "data/{dataset}/{popu}/hg19/{archaic}/curve/anc_matrix/{dataset}_{popu}_mostLikely_Downsample{ds}_chr{chrom}.snp",
#         genmap="/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/hg19/resources/genetic_maps_corrections/hg19/Shared_AfAmMap"
#     output:
#         snp = "data/{dataset}/{popu}/hg19/{archaic}/curve/genmap/{dataset}_{popu}_mostLikely_Downsample{ds}_chr{chrom}_genmap.snp",
#     shell:
#         "/global/home/users/moorjani/bin/interpit -i {input.snp}  -o {output.snp}  -d {input.genmap}"

# rule interpit_anc:
#     input:
#         genmap_snp = "data/{dataset}/{popu}/hg19/{archaic}/curve/genmap/{dataset}_{popu}_mostLikely_Downsample{ds}_chr{chrom}_genmap.snp",
#         anc="data/{dataset}/{popu}/hg19/{archaic}/curve/anc_matrix/{dataset}_{popu}_mostLikely_Downsample{ds}_chr{chrom}.anc",
#     output:
#         anc = "data/{dataset}/{popu}/hg19/{archaic}/curve/genmap/{dataset}_{popu}_mostLikely_Downsample{ds}_chr{chrom}_genmap.anc",
#     shell:
#         """
#         awk 'NR==FNR {{id[$4]=$3; next}} ($2 in id) {{print $1 "\t" id[$2] "\t" $0}}' {input.genmap_snp} {input.anc} | cut --complement -f3,4 > {output.anc}
#         """
        
# this takes >2 days for the big chromosomes 
rule calculate_data_D:
    input:
        anc = "data/{dataset}/{popu}/hg19/{archaic}/curve/genmap/{dataset}_{popu}_mostLikely_Downsample{ds}_chr{chrom}_genmap.anc",
        script = "helper_scripts/calculate_covariances_genmap_rounded.py"
    output:
        "data/{dataset}/{popu}/hg19/{archaic}_results/genmap/{dataset}_{popu}_mostLikely_Downsample{ds}_chr{chrom}_D.txt"
    resources:
        time='80:00:00'
    shell:
        "python {input.script} -i {input.anc} -o {output}"
