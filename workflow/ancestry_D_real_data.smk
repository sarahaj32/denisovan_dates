# A set of rules for analyzing choin datasets

wildcard_constraints:
    chrom="\d+",

# problem is that there currently isn't a method to make the D.txt files
rule data_all:
    input:
        expand("data/choin/curve/DEN_results/genmap/choin_Oceania_fragments_08_{classif}_Oceania{pos}_hg19_chr{chrom}_D.txt",
               classif=["ND01","mostLikely"],
               pos=["Downsample2"], # "Variable","Ingroup", 
               chrom=list(range(1,23))
        ),
        expand("data/choin/curve/NEA_results/genmap/choin_Oceania_fragments_08_{classif}_Oceania{pos}_hg19_chr{chrom}_D.txt",
               classif=["ND10","mostLikely"],
               pos=["Downsample2"], # "Variable","Ingroup", 
               chrom=list(range(1,23))
        ),
        # expand("data/choin/curve/choin_DEN_anc/genmap/choin_Oceania_fragments_08_{classif}_Oceania{pos}_hg19_chr{chrom}_genmap.anc",
        #        classif=["ND01","mostLikely"],
        #        method=["Variable","Ingroup"],
        #        chrom=list(range(1,23))
        # ),
        # expand("data/choin/curve/choin_NEA_anc/genmap/choin_Oceania_fragments_08_{classif}_Oceania{pos}_hg19_chr{chrom}_genmap.anc",
        #        classif=["ND10","mostLikely"],
        #        method=["Variable","Ingroup"],
        #        chrom=list(range(1,23))
        # ),
        # expand("data/choin/curve/{archaic}_results/genmap/choin_Oceania_fragments_08_{classif}_Oceania{pos}_hg19_chr{chrom}_D.txt",
        #        classif=["mostLikely"],
        #        method=["Downsample2"],
        #        chrom=list(range(1,23)), 
        #        archaic=["NEA", "DEN"]),
    output:
        "logs/DATA_all_done.txt",
    shell:
        "echo 'done' > {output}"

rule create_snp_file:
    input:
        anc = "data/choin/curve/choin_{archaic}_anc/choin_Oceania_fragments_08_{classif}_Oceania{pos}_hg19_chr{chrom}.anc",
    output:
        snp = "data/choin/curve/choin_{archaic}_anc/choin_Oceania_fragments_08_{classif}_Oceania{pos}_hg19_chr{chrom}.snp",
    shell:
        """
        awk '{{print $1 ":" $2 "\t" $1 "\t0.0\t" $2, "A\tT"}}' {input.anc} > {output.snp}
        """

rule interpit_snp:
    input:
        snp = "data/choin/curve/choin_{archaic}_anc/choin_Oceania_fragments_08_{classif}_Oceania{pos}_hg19_chr{chrom}.snp",
        genmap="/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/hg19/resources/genetic_maps_corrections/hg19/Shared_AfAmMap"
    output:
        snp = "data/choin/curve/choin_{archaic}_anc/genmap/choin_Oceania_fragments_08_{classif}_Oceania{pos}_hg19_chr{chrom}_genmap.snp",
    shell:
        "/global/home/users/moorjani/bin/interpit -i {input.snp}  -o {output.snp}  -d {input.genmap}"

rule interpit_anc:
    input:
        genmap_snp = "data/choin/curve/choin_{archaic}_anc/genmap/choin_Oceania_fragments_08_{classif}_Oceania{pos}_hg19_chr{chrom}_genmap.snp",
        anc="data/choin/curve/choin_{archaic}_anc/choin_Oceania_fragments_08_{classif}_Oceania{pos}_hg19_chr{chrom}.anc",
    output:
        anc = "data/choin/curve/choin_{archaic}_anc/genmap/choin_Oceania_fragments_08_{classif}_Oceania{pos}_hg19_chr{chrom}_genmap.anc",
    shell:
        """
        awk 'NR==FNR {{id[$4]=$3; next}} ($2 in id) {{print $1 "\t" id[$2] "\t" $0}}' {input.genmap_snp} {input.anc} | cut --complement -f3,4 > {output.anc}
        """
        
# this takes >2 days for the big chromosomes 
rule calculate_data_D:
    input:
        anc = "data/choin/curve/choin_{archaic}_anc/genmap/choin_Oceania_fragments_08_{classif}_Oceania{pos}_hg19_chr{chrom}_genmap.anc",
        script = "helper_scripts/calculate_covariances_genmap_rounded.py"
    output:
        "data/choin/curve/{archaic}_results/genmap/choin_Oceania_fragments_08_{classif}_Oceania{pos}_hg19_chr{chrom}_D.txt"
    resources:
        time='80:00:00'
    shell:
        "python {input.script} -i {input.anc} -o {output}"
