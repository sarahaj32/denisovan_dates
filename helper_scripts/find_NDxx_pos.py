import argparse
from collections import defaultdict
import pandas as pd
import sys
import random

####### TESTING ########
#example with African genotypes:
# file="data/Leo_ND/Leo_chr3_NDpan.vcf"
# python helper_scripts/find_NDxx_pos.py -i $file -o helper_scripts/testing/test_afrND.csv -d DEN_3 -n Altai_1 -a AFR -r
# can use the first 17 rows to test pre-calculated african frequency
# helper functions

def parse_header(line, neanderthal, denisovan,  african, chimp):
    nea_ix = [i for i,name in enumerate(line) if neanderthal == name]
    print(f"NEA {nea_ix}")
    if not nea_ix:
        print(f"No neanderthal population of name {neanderthal} found")
        nea_ix = None
    else:
        nea_ix = nea_ix[0] # single individual

    den_ix = [i for i,name in enumerate(line) if denisovan == name]
    print(f"DEN {den_ix}")
    if not den_ix:
        print(f"No denisovan population of name {denisovan} found")
        den_ix = None
    else:
        den_ix = den_ix[0] # single individual
    
    if african:
        afr_ix = [i for i,name in enumerate(line) if african in name]
        print(f"AFR {afr_ix}")
        if not afr_ix:
            print(f"No african population of name {african} found")
            afr_ix = None
    else:
        afr_ix = None

    if chimp:
        anc_ix = [i for i,name in enumerate(line) if chimp == name]
        print(f"ANC {anc_ix}")
        if not anc_ix:
            print(f"No ancestral population of name {chimp} found")
            anc_ix = None
        else:
            anc_ix = anc_ix[0] # single individual
    else:
        anc_ix = None


    return(nea_ix, den_ix, afr_ix, anc_ix)


def get_der_count(archaic, gt, der):
    archaic_gt = gt[int(archaic[0])] + gt[int(archaic[2])]

    n_der_archaic = archaic_gt.count(der)
    if n_der_archaic == 1:
        test = random.random()
        if test >= 0.5:
            n_der_archaic += 1
        else:
            n_der_archaic -= 1
    return(n_der_archaic)

def get_NDxx(input, neanderthal, denisovan, african, chimp, random_allele, outfile):
    """
    A function that calculates the number of ND10 ND01 positions in a file
    """
    counter = 0
    with open(input, "r") as file1, open(outfile, "w") as of:
        line = ["chrom", "pos", "der", "anc", "n_der_den", "n_der_nea", "chimp", "class"]
        of.write("\t".join(line) + "\n")
        for line in file1:
            line = line.strip().split("\t")
            counter += 1
            if counter % 1000000 == 0:
                print(counter)
            if line[0].startswith("#"):  
                if "#CHROM" in line[0]:
                    nea_ix, den_ix, _, chimp_ix = parse_header(line, neanderthal, denisovan, None, chimp)
                    chrom_ix = line.index("#CHROM")
                    pos_ix = line.index("POS")
                    ref_ix = line.index("REF")
                    alt_ix = line.index("ALT")
            else:
                # get the alleles for each population
                nea = line[nea_ix]
                #print(nea)
                den = line[den_ix]
                #print(den)
                # in the case of simulations where I don't have to polarize
                if chimp_ix:
                    chimp = line[chimp_ix]
                else:
                    chimp = "0/0"
                #print(chimp)
                pos = line[pos_ix]
                chrom = line[chrom_ix]
                
                # determine alternative and reference alleles
                ref = line[ref_ix]
                #print(f"REF: {ref}")
                alt = line[alt_ix]
                #print(f"ALT: {alt}")
                if alt == ".":
                    alt = ref

                gt = [ref, alt]
                #print(gt)

                # only look at biallelic positions where there is a genotype call for all 3 genomes
                alleles = "".join([str(x)[0] + str(x)[2] for x in [nea, den, chimp]])
                if (not "." in alleles) and (all([int(a)<=1 for a in alleles])):

                    # determine which allele is derived, save as "anc, derived"
                    chimp_allele = int(chimp[0])
                    #print(chimp_allele)
                    if chimp_allele == 0:
                        anc, der = ref, alt
                    elif chimp_allele == 1:
                        anc, der = alt, ref
                    else:
                        print("NO CHIMP???")
                    #print(f"anc: {anc}, der: {der}")
                    
                    # determine how many derived alleles are in nea and den (if it's 1, randomly select)
                    # using my helper function
                    n_der_den = get_der_count(den, gt, der) 
                    n_der_nea = get_der_count(nea, gt, der) 
                    # finally get the nd type:
                    if n_der_den > 0 and n_der_nea > 0:
                        nd = "nd11"
                    elif n_der_den > 0 and n_der_nea == 0:
                        nd = "nd01"
                    elif n_der_den == 0 and n_der_nea > 0:
                        nd = "nd10"
                    elif n_der_den == 0 and n_der_nea == 0:
                        nd = "nd00"
                    else:
                        nd = "problem"
                    #print(nd)
                    line = [chrom, pos, der, anc, n_der_den, n_der_nea, 0, nd]
                    of.write("\t".join([str(x) for x in line]) + "\n")


  

 

 
parser = argparse.ArgumentParser("Process input")
parser.add_argument("-i", "--input", help = "name of input file", type = str, required = True)
parser.add_argument("-o", "--output", help = "name of output file", type = str, required = True)
parser.add_argument("-p", "--pan", help = "name of the ancestral population", type = str, default = None) # if not specified, assumes that 0 is ancestral
parser.add_argument("-d", "--denisovan", help = "name of the denisovan population", type = str, required = True) 
parser.add_argument("-a", "--african", help = "name of the african population", type = str, default = None)
parser.add_argument("-n", "--neanderthal", help = "name of the neanderthal population", type = str, default = "NEA")
parser.add_argument("-r", "--random", help = "Flag that says if we should randomly select haplotypes from the archaics", action = 'store_true') # 1 or 2

args = parser.parse_args()

# file parameters
infile = args.input
outfile = args.output

# population parameters
neanderthal = args.neanderthal
denisovan = args.denisovan
african = args.african
chimp = args.pan
# haploid parameter
random_allele = args.random


if __name__ == "__main__":
    get_NDxx(infile, chimp = chimp, african = african, denisovan = denisovan, neanderthal = neanderthal, random_allele = random_allele, outfile = outfile) 

