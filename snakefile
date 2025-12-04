# snakemake --forceall --rulegraph | dot -Tpdf > dag.pdf

import json
 
configfile: "config.yaml"

#include: "workflow/get_archaics.smk"
#include: "workflow/annotate_frags.smk"
# include: "workflow/date_flexi.smk"
# include: "workflow/modern.smk"
# include: "workflow/date.smk"
# include: "workflow/cumulative_analysis.smk"
# include: "workflow/date_modern_myAsc.smk"
include: "workflow/simulate.smk"
# include: "workflow/simulated_analysis.smk"
# include: "workflow/computed_sims.smk"
# include: "workflow/ancestry_D_real_data.smk"
#include: "workflow/create_real_data_anc.smk"
# goal here is to annotate all fragments:
#GAv1_inds_json = "/global/scratch/p2p3/pl1_moorjani/zhangyulin9806/ArchaicMutRate/hmmix_output/GA100K_raw/GA100K.samples.json"

# with open(GAv1_inds_json, 'r') as f:
#     GAv1_inds = json.load(f)
#     GAv1_inds = GAv1_inds["ingroup"]

# bcftools concat data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr1.vcf.gz data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr2.vcf.gz data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr3.vcf.gz data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr4.vcf.gz data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr5.vcf.gz data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr6.vcf.gz data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr7.vcf.gz data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr8.vcf.gz data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr9.vcf.gz data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr10.vcf.gz data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr11.vcf.gz data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr12.vcf.gz data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr13.vcf.gz data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr14.vcf.gz data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr15.vcf.gz data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr16.vcf.gz data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr17.vcf.gz data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr18.vcf.gz data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr19.vcf.gz data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr20.vcf.gz data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr21.vcf.gz data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr22.vcf.gz | bcftools sort -T $TMPDIR | bcftools annotate --rename-chrs data/resources/chr_map_addChr -o data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL.vcf.gz

# bcftools concat data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr1.vcf.gz data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr2.vcf.gz data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr3.vcf.gz data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr4.vcf.gz data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr5.vcf.gz data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr6.vcf.gz data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr7.vcf.gz data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr8.vcf.gz data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr9.vcf.gz data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr10.vcf.gz data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr11.vcf.gz data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr12.vcf.gz data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr13.vcf.gz data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr14.vcf.gz data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr15.vcf.gz data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr16.vcf.gz data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr17.vcf.gz data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr18.vcf.gz data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr19.vcf.gz data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr20.vcf.gz data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr21.vcf.gz data/1000G/PJL/modern/vcf/hg19/1000G_snps10_08_allChr_modern_PJL_var_chr22.vcf.gz | bcftools sort -T $TMPDIR | bcftools annotate --rename-chrs data/resources/chr_map_addChr -o combo_sort_test.vcf

rule all:
   input:
        # expand("data/{dataset}/{dataset}_fragments_{snps}08_annotated.txt", dataset = ["1000G", "HGDP", "GAv1", "LASIDAD"], snps = ["", "snps10_"]),
        # expand("data/GA100K/hmmix_strict/02decode/{ind}.hap{hap}.txt", ind = GAv1_inds, hap = [1,2]),
        # "data/combined_fragments/HGDP_LASIDAD_1000G_choin_fragments_08_annotated.txt",
        # expand("data/1000G/{popu}/modernMyAsc/dating/computed/1000G_modernMyAsc_{popu}_{type}_chr{chrom}_results.out", popu = "EAS", type = ["ND01", "ND10"], chrom = list(range(1,23))),
        # expand("data/HGDP/{popu}/modernMyAsc/dating/computed/HGDP_modernMyAsc_{popu}_{type}_chr{chrom}_results.out", popu = "OCEANIA", type = ["ND01", "ND10"], chrom = list(range(1,23))),
        # expand("simdat/denisovan_simple/raw/denisovan_simple_sim_{rep}.vcf", rep = list(range(1,21))),
        # expand("data/Archaics/hg19/AltDenVinPanHSanc_merged_strictArchaic_hg19_chr{chromosome}_genmap.snp", chromosome = list(range(1,23))),
        # expand("data/Archaics/hg19/AltDenVinPanHSanc_merged_strictArchaic_hg19_chr{chromosome}_genmap_problemPositions.txt", chromosome = list(range(1,23))),
        # EVERYTHING ABOVE WAS PREVIOUSLY UNCOMMENTED
        # GAv1 with extra info:
        # also need the hg19 archaics list to transfer over for choin
        # expand("data/{dataset}/{popu}/fragments/{popu}_fragments_snps10_08_nonoverlapping_{thresh}.bed", dataset = ["GAv1"], popu = ["AET", "ATI", "PAP", "KOR", "JPN", "Indonesia", "Philippines", "SEA", "SAS", "NEA", "OCE", "Mongolia", "Singapore", "India"], thresh = ["Nea", "Dhigh08"]),
        # expand("data/{dataset}/{popu}/fragments/{popu}_fragments_snps10_08_nonoverlapping_{thresh}.bed", dataset = ["HGDP"], popu = ["OCEANIA"], thresh = ["Nea", "Dhigh08"]),
        # expand("data/{dataset}/{popu}/fragments/{popu}_fragments_snps10_08_nonoverlapping_{thresh}.bed", dataset = ["LASIDAD"], popu = ["North", "South", "West", "East", "Central", "North-East", "singlePulse", "multiPulse", "noMultiPulse"], thresh = ["Nea", "Dhigh08"]),
        # expand("data/{dataset}/{popu}/fragments/{popu}_fragments_snps10_08_nonoverlapping_{thresh}.bed", dataset = ["1000G"], popu = ["EAS", "CHB", "PJL"], thresh = ["Nea", "Dhigh08"]),
        # I want these, just not now!!!
        # expand("data/cumulative/{population}_cumulative_{archaic}_fragments_08_annotated.txt", archaic = ["Denisova", "Neanderthal"], population = ["choin-OCE", "choin-Agta", "choin-OCESEA", "all-PAP", "HGDP-PAP", "GAv1-Aeta", "all-Aeta", "GAv1-Indonesia", "GAv1-Philippines", "GAv1-SEA", "SAS", "EAS", "combined_metadata", "GAv1-NEA", "GAv1-SAS", "LASI", "1000G-EAS", "1000G-SAS"]),
        # FOR NICK!!!!
        #expand("data/cumulative/{population}_cumulative_{archaic}_fragments_08_annotated.txt", archaic = ["Denisova"], population = ["SEA_Agta"]),
        # expand("data/cumulative/allInds_{dataset}_cumulative_{archaic}_fragments_08_annotated.txt", archaic = ["Denisova", "Neanderthal"], dataset = ["GAv1", "HGDP", "LASIDAD", "1000G", "choin"]),
         #expand("data/LASIDAD/{popu}/modernMyAsc/dating/computed/LASIDAD_modernMyAsc_{popu}_{type}_chr{chrom}_results.out", popu = "noMultiPulse", type = ["ND01", "ND10"], chrom = list(range(1,23))),
        #"logs/sims_all_done.txt",
        #"logs/CHECK_all_done.txt",
        # get all variable sites for GBR (other populations?)
        # expand( "data/1000G/{popu}/modern/vcf/hg38/1000G_{popu}_strictMask_var_chr{chrom}.vcf", chrom = list(range(1,23)), popu = ["GBR", "TSI", "CEU", "PJL", "ITU", "SAS"]),
        # expand("data/1000G/{popu}/fragments/{popu}_fragments_08_annotated.bed", popu = ["GBR", "TSI", "CEU", "PJL", "ITU", "SAS"]), # fragments_snps10_08_annotated
        # expand("data/1000G/{popu}/fragments/NEA/{popu}_fragments_08_annotated_mostLikely.bed", popu = ["GBR", "TSI", "CEU", "PJL", "ITU", "SAS"]),
        # expand("data/1000G/{popu}/modern/vcf/hg38/1000G_{popu}_strictMask_var_Downsample{ds}_chr{chrom}.vcf", ds = [2,4], chrom = list(range(1,23)), popu = ["GBR", "TSI", "CEU", "PJL", "ITU", "SAS"]),
        # expand("data/{dataset}/{popu}/fragments/{archaic}/{popu}_fragments_08_annotated_mostLikely_chr{chrom}.bed", dataset = "1000G", archaic = ["NEA", "DEN"], chrom = list(range(1,23)), popu = ["GBR", "TSI", "CEU", "PJL", "ITU", "SAS"]),
        # expand("data/1000G/{popu}/hg38/{archaic}/curve/anc_matrix/1000G_{popu}_mostLikely_Downsample{ds}_chr{chrom}.anc", archaic = ["NEA", "DEN"], ds = [2,4], chrom = list(range(1,23)), popu = ["GBR", "TSI", "CEU", "PJL", "ITU", "SAS"]),
        # expand("data/{dataset}/{popu}/hg19/{archaic}/curve/anc_matrix/{dataset}_{popu}_mostLikely_Downsample{ds}_chr{chrom}.anc", dataset = "1000G", archaic = ["NEA", "DEN"], ds = [2,4], chrom = list(range(1,23)), popu = ["GBR", "TSI", "CEU", "PJL", "ITU", "SAS"]),
        # expand("data/{dataset}/{popu}/hg19/{archaic}_results/genmap/{dataset}_{popu}_mostLikely_Downsample{ds}_chr{chrom}_D.txt", dataset = "1000G", archaic = ["NEA", "DEN"], ds = [4], chrom = list(range(1,23)), popu = ["GBR", "CEU", "PJL", "ITU", "SAS"]),
        # "data/LASIDAD/LASIDAD_fragments_08_annotated.txt",
        #expand("data/{dataset}/{popu}/hg19/{archaic}_results/genmap/{dataset}_{popu}_mostLikely_Downsample{ds}_chr{chrom}_D.txt", dataset = "LASIDAD", archaic = ["NEA", "DEN"], ds = [6], chrom = list(range(1,23)), popu = ["North", "Central", "South", "West", "East", "North-East"]),
        # SARAH RUNNING TODAY
        # these submitted to savio4_htc
        #expand("data/1000G/{popu}/{kind}/dating/computed/{asc}_dateCorrections.out", popu = ["CHB", "BEB", "CEU", "CHS", "GBR", "IBS", "JPT", "PJL", "STU", "CDX", "GIH", "ITU", "KHV",  "TSI", "FIN", "PEL", "MXL", "CLM", "PUR", "EAS", "SAS"], asc = [ "ND10AfrAncLeo", "ND01AfrAncLeo"], kind = ["modern"]), # "ND01AfrAnc", "ND10AfrAnc", "ND01",, "modernRefFill"
        # this submitted to savio3_htc
        #expand("data/1000G/{popu}/{kind}/dating/LD_covariance/{asc}_results_chr{chrom}.txt", popu = ["CHB", "BEB", "CEU", "CHS", "GBR", "IBS", "JPT", "PJL", "STU", "CDX", "EAS", "GIH", "ITU", "KHV", "SAS", "TSI", "FIN", "PEL", "MXL", "CLM", "PUR"], asc = ["ND01AfrAnc", "ND10AfrAnc", "ND01", "ND10"], kind = ["modern", "modernRefFill"], chrom = list(range(1,23)))
        #expand("data/1000G/{popu}/{kind}/dating/eig/1000G_snps10_08_allChr_{kind}_{popu}_{asc}_genmap_chr{chrom}.snpgeno", popu = ["CHB", "BEB", "CEU", "CHS", "GBR", "IBS", "JPT", "PJL", "STU", "CDX", "EAS", "GIH", "ITU", "KHV", "SAS", "TSI", "FIN", "PEL", "MXL", "CLM", "PUR"], asc = ["ND01AfrAnc", "ND10AfrAnc", "ND01", "ND10"], kind = ["modern", "modernRefFill"], chrom = list(range(1,23))),
       # "logs/DATA_all_done.txt",
        # "logs/CHECK_all_done.txt"
        # expand("logs/sims_generated_{sim}_{seed}.txt", sim = [config["sim"]], seed = [config["seed"]])
        #expand("simdat/{sim}_{seed}/{archaic}/curve/results/{sim}_{nMH}NAMH_{frag_source}_{pos}_{rep}_D.txt", seed = [config["seed"]], nMH = [config["nMH"]], frag_source = [config["frag"]], sim = [config["sim"]], archaic = ["NEA", "DEN"], rep = list(range(1,chroms)), pos = [config["pos"]]),
        # generating simulated data:
        # expand("simdat/{sim}_{seed}/{archaic}/curve/fragments/{sim}_{nMH}NAMH_{frag_source}_{pos}_{rep}.anc", seed = [config["seed"]], nMH = [100, 20], frag_source = ["hmmix.0.8.200Afr.fullArchaic.mostLikely", "hmmix.0.9.200Afr.fullArchaic.mostLikely", "hmmix.0.7.200Afr.fullArchaic.mostLikely", "ibdmix", "hmmix.0.8.200Afr.fullArchaic.mostLikely.50k", "hmmix.0.8.200Afr.fullArchaic.1NeaMostLikely"], sim = [config["sim"]], archaic = ["NEA", "DEN"], rep = list(range(1,chroms)), pos = ["ds200k"]),
        # expand("simdat/{sim}_{seed}/{archaic}/fragments/{sim}_{nMH}NAMH_simIntro_{rep}.bed", seed = [config["seed"]], nMH = [20, 100], sim = [config["sim"]], archaic = ["DEN", "NEA"], rep = list(range(1,chroms)))
        # getting dates from simulated data
        expand("simdat/results/{sim}_{nMH}NAMH_{frag_source}_{pos}_{archaic}_{seed}_min{anal_min}_inferred_parameters.txt", seed = [config["seed"]], nMH = [config["nMH"]], frag_source =[config["frag"]], sim = [config["sim"]], archaic = ["NEA", "DEN"], rep = list(range(1,chroms)), pos = [config["pos"]], anal_min = [0, 0.01, 0.02, 0.03, 0.05, 0.1]),
        
        # collating performance results

        

# data/1000G/CDX/modern/dating/computed/ND10_dateCorrections.out
