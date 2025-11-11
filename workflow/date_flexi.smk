# this is the exact same as date.smk but with added flexibility to change paths and names across different datasets/analysis
rule convert_to_eig:
    input:
        vcf="data/{realDatPath}/vcf/hg19/{name}.vcf",
        script="/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/bin/vcf2eigenstrat.py"
    output:
        geno="data/{realDatPath}/dating/eig/{name}.geno",
        ind="data/{realDatPath}/dating/eig/{name}.ind",
        snp="data/{realDatPath}/dating/eig/{name}.snp"
    params:
        path="data/{realDatPath}/dating/eig/{name}"
    shell:
        """ 
        python2 {input.script} -v {input.vcf} -o {params.path}
        sed -i 's/chr//g' {output.snp} 
        """ 

rule interpit:
    # cannot have "chr" in the .snp file
    input:
        snp="data/{realDatPath}/dating/eig/{name}.snp",
        genmap="/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/hg19/resources/genetic_maps_corrections/hg19/Shared_AfAmMap"
    output:
        snp="data/{realDatPath}/dating/eig/{name}_genmap.snp",
    shell:
        "/global/home/users/moorjani/bin/interpit -i {input.snp}  -o {output.snp}  -d {input.genmap}"

rule make_computed_parfile:
    input:
        geno="{path}/dating/eig/{name}.geno",
        snp="{path}/dating/eig/{name}_genmap.snp",
        ind="{path}/dating/eig/{name}.ind"
    params:
        out="{path}/dating/computed/{name}_results.out",
    output:
        parfile="{path}/dating/computed/{name}_computed_parfile",
    shell:
        """
        echo "genotypename: {input.geno}" > {output.parfile}
        echo "snpname: {input.snp}" >> {output.parfile}
        echo "indivname: {input.ind}" >> {output.parfile}
        echo "maxdis: 0.01">> {output.parfile}
        echo "binsize: 0.00001">> {output.parfile}
        echo "output: {params.out}" >> {output.parfile}
        """

rule computed:
    input:
        geno="{path}/dating/eig/{name}.geno",
        snp="{path}/dating/eig/{name}_genmap.snp",
        ind="{path}/dating/eig/{name}.ind",
        parfile="{path}/dating/computed/{name}_computed_parfile",
    output:
        out="{path}/dating/computed/{name}_results.out"
    resources:
        time="24:00:00",
        cpus="10"
    shell:
        "/global/home/users/moorjani/bin/computed -p {input.parfile}"

rule make_date_correction_parfile:
    input:
        computed_results="{path}/dating/computed/{name}_results.out",
        correction="/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/hg19/resources/genetic_maps_corrections/hg19/Shared_AfAmMap/alpha_1kg.shared/alpha.txt"
    params:
        out="{path}/dating/computed/{name}_dateCorrections.out"
    output:
        parfile="{path}/dating/computed/{name}_dateCorrection_parfile"
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
        computed_results="{path}/dating/computed/{name}_results.out",
        parfile="{path}/dating/computed/{name}_dateCorrection_parfile"
    output:
        out="{path}/dating/computed/{name}_dateCorrections.out"
    resources:
        time="10:00:00"
    shell:
        "/global/home/users/moorjani/bin/getposterior -p {input.parfile}"

# rule make_date_2comp_correction_parfile:
#     input:
#         computed_results="{path}/dating/computed/{name}_results.out",
#         correction="/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/hg19/resources/genetic_maps_corrections/hg19/Shared_AfAmMap/alpha_1kg.shared/alpha.txt",
#     params:
#         out="{path}/dating/two_components/{name}_dateCorrections.out"
#     output:
#         parfile="{path}/dating/two_components/{name}_dateCorrection_parfile"
#     shell:
#         """
#         echo "col: 1" > {output.parfile}
#         echo "l: 0.02" >> {output.parfile}
#         echo "h: 1.0" >> {output.parfile}
#         echo "inp: {input.computed_results}" >> {output.parfile}
#         echo "debug: 1" >> {output.parfile}
#         echo "alpha: {input.correction}" >> {output.parfile}
#         echo "output: {params.out}" >> {output.parfile}
#         echo "print_interval: 10" >> {output.parfile}
#         echo "glb: 25" >> {output.parfile}
#         echo "gub: 33" >> {output.parfile}
#         echo "mcmc_iters: 1500" >> {output.parfile}
#         echo "burnin: 200" >> {output.parfile}
#         echo "seed: 1" >> {output.parfile}
#         echo "ncomponents: 2" >> {output.parfile}
#         """