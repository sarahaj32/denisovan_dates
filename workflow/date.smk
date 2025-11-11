# filter a vcf to only ND01 or ND10 positions, with an option to limit to only fixed ancestral in YRI (AfrAnc)
rule filter_to_nd:
    input:
        vcf = "data/{dataset}/{popu}/{kind}/vcf/hg19/{dataset}_snps10_08_allChr_{kind}_{popu}.vcf.gz", 
        pos = "data/hg19/snp_annotations/{asc}_hg19.txt"
    output:
        "data/{dataset}/{popu}/{kind}/vcf/hg19/{dataset}_snps10_08_allChr_{kind}_{popu}_{asc}.vcf"
    shell:
        "bcftools view {input.vcf} -T {input.pos} > {output}"

rule convert_to_eig:
    input:
        vcf="data/{dataset}/{popu}/{kind}/vcf/hg19/{dataset}_snps10_08_allChr_{kind}_{popu}_{asc}.vcf",
        script="/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/bin/vcf2eigenstrat.py"
    output:
        geno="data/{dataset}/{popu}/{kind}/dating/eig/{dataset}_snps10_08_allChr_{kind}_{popu}_{asc}.geno",
        ind="data/{dataset}/{popu}/{kind}/dating/eig/{dataset}_snps10_08_allChr_{kind}_{popu}_{asc}.ind",
        snp="data/{dataset}/{popu}/{kind}/dating/eig/{dataset}_snps10_08_allChr_{kind}_{popu}_{asc}.snp"
    params:
        path="data/{dataset}/{popu}/{kind}/dating/eig/{dataset}_snps10_08_allChr_{kind}_{popu}_{asc}"
    shell:
        """ 
        python2 {input.script} -v {input.vcf} -o {params.path}
        sed -i 's/chr//g' {output.snp} 
        """ 
    
rule interpit:
    # cannot have "chr" in the .snp file
    input:
        snp="data/{dataset}/{popu}/{kind}/dating/eig/{dataset}_snps10_08_allChr_{kind}_{popu}_{asc}.snp",
        genmap="/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/hg19/resources/genetic_maps_corrections/hg19/Shared_AfAmMap"
    output:
        snp="data/{dataset}/{popu}/{kind}/dating/eig/{dataset}_snps10_08_allChr_{kind}_{popu}_{asc}_genmap.snp",
    params:
        tmp="data/{dataset}/{popu}/{kind}/dating/eig/{dataset}_snps10_08_allChr_{kind}_{popu}_{asc}_genmap.snp_tmp",
    shell:
        """
        /global/home/users/moorjani/bin/interpit -i {input.snp}  -o {output.snp}  -d {input.genmap}
        grep -v -E "4:69266058|4:69280847|4:69282416|4:69283419|4:69285362" {output.snp} > {params.tmp} && mv {params.tmp} {output.snp}
        """

rule make_computed_parfile:
    input:
        geno="data/{dataset}/{popu}/{kind}/dating/eig/{dataset}_snps10_08_allChr_{kind}_{popu}_{asc}.geno",
        snp="data/{dataset}/{popu}/{kind}/dating/eig/{dataset}_snps10_08_allChr_{kind}_{popu}_{asc}_genmap.snp",
        ind="data/{dataset}/{popu}/{kind}/dating/eig/{dataset}_snps10_08_allChr_{kind}_{popu}_{asc}.ind"
    params:
        out="data/{dataset}/{popu}/{kind}/dating/computed/{asc}_results.out",
        #jk_out="data/{dataset}/{popu}/{kind}/dating/jacknife/"
    output:
        parfile="data/{dataset}/{popu}/{kind}/dating/computed/{asc}_computed_parfile",
    shell:
        """
        echo "genotypename: {input.geno}" > {output.parfile}
        echo "snpname: {input.snp}" >> {output.parfile}
        echo "indivname: {input.ind}" >> {output.parfile}
        echo "maxdis: 0.01">> {output.parfile}
        echo "binsize: 0.00001">> {output.parfile}
        echo "output: {params.out}" >> {output.parfile}
        """

        # echo "jackknife: {params.jk_out}" >> {output.parfile}
        # echo "jackblock: 10" >> {output.parfile}

rule computed:
    input:
        geno="data/{dataset}/{popu}/{kind}/dating/eig/{dataset}_snps10_08_allChr_{kind}_{popu}_{asc}.geno",
        snp="data/{dataset}/{popu}/{kind}/dating/eig/{dataset}_snps10_08_allChr_{kind}_{popu}_{asc}_genmap.snp",
        ind="data/{dataset}/{popu}/{kind}/dating/eig/{dataset}_snps10_08_allChr_{kind}_{popu}_{asc}.ind",
        parfile="data/{dataset}/{popu}/{kind}/dating/computed/{asc}_computed_parfile",
    output:
        out="data/{dataset}/{popu}/{kind}/dating/computed/{asc}_results.out"
    resources:
        time="30:00:00",
        cpus="2"
    shell:
        "/global/home/users/moorjani/bin/computed -p {input.parfile}"

rule make_snpgeno:
    input:
        geno="data/{dataset}/{popu}/{kind}/dating/eig/{dataset}_snps10_08_allChr_{kind}_{popu}_{asc}.geno",
        snp="data/{dataset}/{popu}/{kind}/dating/eig/{dataset}_snps10_08_allChr_{kind}_{popu}_{asc}_genmap.snp",
    output:
        snpgeno="data/{dataset}/{popu}/{kind}/dating/eig/{dataset}_snps10_08_allChr_{kind}_{popu}_{asc}_genmap.snpgeno"
    shell:
        """
        paste <(awk '{{print $2 "\t" $3}}' {input.snp}) {input.geno} > {output.snpgeno}
        """

rule LD_covariance:
    input:
        snpgeno="data/{dataset}/{popu}/{kind}/dating/eig/{dataset}_snps10_08_allChr_{kind}_{popu}_{asc}_genmap.snpgeno",
        script = "helper_scripts/calculate_LD.py"
    output:
        out=expand("data/{{dataset}}/{{popu}}/{{kind}}/dating/LD_covariance/{{asc}}_results_chr{chrom}.txt", chrom = list(range(1,23)))
    params:
        output="data/{dataset}/{popu}/{kind}/dating/LD_covariance/{asc}_results.txt"
    resources:
        time="20:00:00",
    shell:
        "python {input.script} -i {input.snpgeno} -o {params.output}"


rule make_date_correction_parfile:
    input:
        computed_results="data/{dataset}/{popu}/{kind}/dating/computed/{asc}_results.out",
        correction="/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/hg19/resources/genetic_maps_corrections/hg19/Shared_AfAmMap/alpha_1kg.shared/alpha.txt"
    params:
        out="data/{dataset}/{popu}/{kind}/dating/computed/{asc}_dateCorrections.out"
    output:
        parfile="data/{dataset}/{popu}/{kind}/dating/computed/{asc}_dateCorrection_parfile"
    shell:
        """
        echo "col: 1" > {output.parfile}
        echo "l: 0.02" >> {output.parfile}
        echo "h: 1.0" >> {output.parfile}
        echo "inp: {input.computed_results}" >> {output.parfile}
        echo "debug: 1" >> {output.parfile}
        echo "alpha: {input.correction}" >> {output.parfile}
        echo "output: {params.out}" >> {output.parfile}
        echo "print_interval: 10" >> {output.parfile}
        echo "glb: 25" >> {output.parfile}
        echo "gub: 33" >> {output.parfile}
        echo "mcmc_iters: 1500" >> {output.parfile}
        echo "burnin: 200" >> {output.parfile}
        echo "seed: 1" >> {output.parfile}
        """

rule get_dates:
    input:
        computed_results="data/{dataset}/{popu}/{kind}/dating/computed/{asc}_results.out",
        parfile="data/{dataset}/{popu}/{kind}/dating/computed/{asc}_dateCorrection_parfile"
    output:
        out="data/{dataset}/{popu}/{kind}/dating/computed/{asc}_dateCorrections.out"
    resources:
        time="10:00:00"
    shell:
        "/global/home/users/moorjani/bin/getposterior -p {input.parfile}"

rule make_date_2comp_correction_parfile:
    input:
        computed_results="data/{dataset}/{popu}/{kind}/dating/computed/{asc}_results.out",
        correction="/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/hg19/resources/genetic_maps_corrections/hg19/Shared_AfAmMap/alpha_1kg.shared/alpha.txt",
    params:
        out="data/{dataset}/{popu}/{kind}/dating/two_components/{asc}_dateCorrections.out"
    output:
        parfile="data/{dataset}/{popu}/{kind}/dating/two_components/{asc}_dateCorrection_parfile"
    shell:
        """
        echo "col: 1" > {output.parfile}
        echo "l: 0.02" >> {output.parfile}
        echo "h: 1.0" >> {output.parfile}
        echo "inp: {input.computed_results}" >> {output.parfile}
        echo "debug: 1" >> {output.parfile}
        echo "alpha: {input.correction}" >> {output.parfile}
        echo "output: {params.out}" >> {output.parfile}
        echo "print_interval: 10" >> {output.parfile}
        echo "glb: 25" >> {output.parfile}
        echo "gub: 33" >> {output.parfile}
        echo "mcmc_iters: 1500" >> {output.parfile}
        echo "burnin: 200" >> {output.parfile}
        echo "seed: 1" >> {output.parfile}
        echo "ncomponents: 2" >> {output.parfile}
        """
