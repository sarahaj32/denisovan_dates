# need to liftover choins fragments to hg19
# rule liftover_beds_tohg19:
#     input:
#         fragments="data/{dataset}/{dataset}_fragments_08_annotated.txt",
#     output:
#         fragments="data/choin/choin_fragments_snps10_08_nonoverlapping_Dhigh08_hg19.bed"
#     shell:
#         """
#         Rscript {input.script} {input.fragments} {output.fragments}
#         """

rule combine_fragments:
    input:
        fragments=expand("data/{dataset}/{dataset}_fragments_08_annotated.txt", dataset = ["GAv1", "HGDP", "LASIDAD", "1000G", "choin"]) 
    output:
        combined="data/combined_fragments/HGDP_LASIDAD_1000G_choin_fragments_08_annotated.txt"
    shell:
        """
        head -n 1 {input.fragments[0]} > {output.combined}
        awk 'FNR==1 && NR!=1 {{next}} 1' {input.fragments} >> {output.combined}
        """

# make populaion files by hand, then use this and my R script to make a whole bunch of cumulative analysis files
# using shell script within this folder
rule filter_fragments_to_pop:
    input:
        fragments="data/combined_fragments/HGDP_LASIDAD_1000G_choin_fragments_08_annotated.txt",
        inds="data/combined_fragments/{population}.inds"
    output:
        data="data/combined_fragments/{population}_combined_fragments_08_annotated.txt"
    shell:
        """
        head -n 1 {input.fragments} > {output.data}
        awk 'NR==FNR {{keep[$1]; next}} $1 in keep' {input.inds} {input.fragments} >> {output.data}
        """

# get cumulative curves for each population of interest (including combined across datasets)
rule get_cumulative:
    input:
        data="data/combined_fragments/{population}_combined_fragments_08_annotated.txt",
        script="helper_scripts/get_cumulative.R"
    output:
        data="data/cumulative/{population}_cumulative_{archaic}_fragments_08_annotated.txt"
    resources:
        time="8:00:00",    
	cpus="6"
    shell:
        """
        Rscript {input.script} {input.data} {wildcards.archaic} {output.data}
        """

# also get cumulative curves for full datasets
rule get_cumulative_datasets:
    input:
        data="data/{dataset}/{dataset}_fragments_08_annotated.txt",
        script="helper_scripts/get_cumulative.R"
    output:
        data="data/cumulative/allInds_{dataset}_cumulative_{archaic}_fragments_08_annotated.txt"
    resources:
        time="8:00:00",
    shell:
        """
        Rscript {input.script} {input.data} {wildcards.archaic} {output.data}
        """
