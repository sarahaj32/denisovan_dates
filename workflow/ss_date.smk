
rule make_indroll_parfile:
    input:
        geno="data/{dataset}/{popu}/{kind}/dating/eig/{dataset}_snps10_08_allChr_{kind}_{popu}_{asc}.geno",
        snp="data/{dataset}/{popu}/{kind}/dating/eig/{dataset}_snps10_08_allChr_{kind}_{popu}_{asc}_genmap.snp",
        ind_file="data/{dataset}/{popu}/{kind}/dating/eig/{dataset}_snps10_08_allChr_{kind}_{popu}_{asc}.ind"
    params:
        out="data/{dataset}/{popu}/{kind}/dating/indroll/{indiv}/{indiv}_{asc}_results.out",
        ind_name="{indiv}"
    output:
        parfile="data/{dataset}/{popu}/{kind}/dating/indroll/{indiv}/{indiv}_{asc}_indroll_parfile",
    shell:
        """
        echo "genotypename: {input.geno}" > {output.parfile}
        echo "snpname: {input.snp}" >> {output.parfile}
        echo "indivname: {input.ind_file}" >> {output.parfile}
        echo "maxdis: 0.01">> {output.parfile}
        echo "binsize: 0.00001">> {output.parfile}
        echo "output: {params.out}" >> {output.parfile}
        echo "idname: {params.ind_name}" >> {output.parfile}
        echo "jackknife: YES" >> {output.parfile}
        """

rule indroll:
    input:
        geno="data/{dataset}/{popu}/{kind}/dating/eig/{dataset}_snps10_08_allChr_{kind}_{popu}_{asc}.geno",
        snp="data/{dataset}/{popu}/{kind}/dating/eig/{dataset}_snps10_08_allChr_{kind}_{popu}_{asc}_genmap.snp",
        ind_file="data/{dataset}/{popu}/{kind}/dating/eig/{dataset}_snps10_08_allChr_{kind}_{popu}_{asc}.ind",
        parfile="data/{dataset}/{popu}/{kind}/dating/indroll/{indiv}/{indiv}_{asc}_indroll_parfile"
    output:
        out="data/{dataset}/{popu}/{kind}/dating/indroll/{indiv}/{indiv}_{asc}_results.out"
    shell:
        "/global/home/users/moorjani/bin/indroll_v2 -p {input.parfile}"

rule make_date_indroll_correction_parfile:
    input:
        indroll_results="data/{dataset}/{popu}/{kind}/dating/indroll/{indiv}/{indiv}_{asc}_results.out",
        correction="/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/hg19/resources/genetic_maps_corrections/hg19/Shared_AfAmMap/alpha_1kg.shared/alpha.txt"
    params:
        out="data/{dataset}/{popu}/{kind}/dating/indroll/{indiv}/{indiv}_{asc}_dateCorrections.out"
    output:
        parfile="data/{dataset}/{popu}/{kind}/dating/indroll/{indiv}/{indiv}_{asc}_dateCorrection_parfile"
    shell:
        """
        echo "col: 1" > {output.parfile}
        echo "l: 0.02" >> {output.parfile}
        echo "h: 1.0" >> {output.parfile}
        echo "inp: {input.indroll_results}" >> {output.parfile}
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

rule get_ss_dates:
    input:
        indroll_results="data/{dataset}/{popu}/{kind}/dating/indroll/{indiv}/{indiv}_{asc}_results.out",
        parfile="data/{dataset}/{popu}/{kind}/dating/indroll/{indiv}/{indiv}_{asc}_dateCorrection_parfile"
    output:
        out="data/{dataset}/{popu}/{kind}/dating/indroll/{indiv}/{indiv}_{asc}_dateCorrections.out"
    resources:
        time="10:00:00"
    shell:
        "/global/home/users/moorjani/bin/getposterior -p {input.parfile}"
