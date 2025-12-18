# a function that creates a geno-like file of ancestry at every position
# inputs: a vcf that is used to determine each position of interest. Note, this is only used to determine positions so no genotype info is necessary
#         a bed file with the archaic regions annotated. The first 4 positions must be haplotype (or ID if diploid), chrom, start, end
# output:
#         a text file with the ancestry label for each individual at each position. 0 = no archaic ancestry, 1 = heterozygous archaic ancestry, 2 = homozygous archaic ancestry

import argparse
from collections import defaultdict

def read_bedfile(bed):
    # read in the bedfile into a dictionary of haplotype to list of [start, end] positions
    hap_dict = defaultdict(list)
    with open(input_bed, 'r') as bed_file:
        hap_idx = 0
        start_idx = 2
        end_idx = 3
        for bed_line in bed_file:
            if "chrom" in bed_line:
                bed_line = bed_line.strip().split("\t")
                hap_idx = [i for i, val in enumerate(bed_line) if "hap" in val or "hapID" in val][0]
                start_idx = bed_line.index("start")
                end_idx = bed_line.index("end")
                continue

            bed_line = bed_line.strip().split("\t")
            # remove "hap" from the haplotype column if it exists
            hap = int(bed_line[hap_idx].replace("hap", ""))
            # iterate through all haplotypes
            start = int(bed_line[start_idx])
            end = int(bed_line[end_idx])
            hap_dict[hap].append([start, end])

    return(hap_dict)


def find_ancestry(hap_dict, input_vcf, output_txt, ploidy) -> None:
    count = 0 
    print(ploidy)
    haps = len(hap_dict)
    with open(input_vcf, 'r') as vcf_file, open(output_txt, 'w') as out_file:
        for full_line in vcf_file:
            count += 1
            if count % 10000 == 0:
                print(count)
            # hard code the indexes for POS and CHROM for a text file (these will update if a vcf)
            pos_ix = 1
            chrom_ix = 0
            if full_line.startswith('#'):
                if "#CHROM" in full_line:
                    line = full_line.strip().split("\t")
                    pos_ix = line.index("POS")
                    chrom_ix = line.index("#CHROM")
            else:
                ancestry_labs = [0 for i in range(haps)]
                line = full_line.strip().split("\t")
                pos = int(line[pos_ix])
                chrom = int(line[chrom_ix].replace("chr", ""))
                snp = f"{chrom}:{pos}"
                # check if the position is in any of the introgressed fragments
                if any([pos > hap_dict[i][0][0] for i in hap_dict]):
                    # print(pos)
                    # print("COULD BE IN FRAGMENT")
                    # for each individual, check if the position is in any of their introgressed fragments
                    # loop through every other individual * 2 (we want both haplotypes from each individual):
                    if ploidy == "haploid":
                        p = 2
                    elif ploidy == "diploid":
                        p = 1
                    else:
                        print(f"ploidy must be 'haploid' or 'diploid', you provided {ploidy}")
                        break
                    
                    # this will then iterate through either each individual or each haplotype
                    for ind in range(0, haps, p):
                        if ind in hap_dict:
                            # print(hap_dict[ind][0])
                            # if the position is after the first fragment's end, go to the next fragment:
                            while pos > hap_dict[ind][0][1]:
                                # print("TESTING")
                                # print(hap_dict)
                                # print(hap_dict[ind])
                                # print(hap_dict[ind][0])
                                # print(f"Position {pos} is greater than end of fragment {hap_dict[ind][0][1]}, removing fragment")
                                hap_dict[ind].pop(0)
                                # if we have gone through all of the fragments, set to infinity so we don't keep checking
                                if len(hap_dict[ind]) == 0:
                                    hap_dict[ind] = [[float('inf'), float('inf')]]
                                # print("NEW")
                                # print(hap_dict[ind])
                            #hap_dict[ind][0][0]
                            # if the position is in an introgressed fragment, set ancestry label to 1
                            if pos >= hap_dict[ind][0][0] and pos <= hap_dict[ind][0][1]:
                                ancestry_labs[ind] = 1
                        # in the case that we are counting up the arhcaic ancestry on both haplotypes
                        if ploidy == "haploid":
                            if ind + 1 in hap_dict:
                                while pos > hap_dict[ind + 1][0][1]:
                                    # print("TESTING")
                                    # print(hap_dict)
                                    # print(hap_dict[ind + 1])
                                    # print(hap_dict[ind + 1][0])
                                    # print(f"Position {pos} is greater than end of fragment {hap_dict[ind + 1][0][1]}, removing fragment")
                                    hap_dict[ind + 1].pop(0)
                                    if len(hap_dict[ind + 1]) == 0:
                                        hap_dict[ind + 1] = [[float('inf'), float('inf')]]
                                hap_dict[ind + 1][0][0]
                                if pos >= hap_dict[ind + 1][0][0] and pos <= hap_dict[ind + 1][0][1]:
                                # print(f"Position {pos} on chromosome {chrom} is in introgressed fragment in individual {ind} on haplotype 2")
                                    ancestry_labs[ind] += 1
                if ploidy == "haploid":
                    filtered = ancestry_labs[::2]
                else:
                    filtered = ancestry_labs
                # print(filtered)
                fields = [str(chrom), str(pos)] + [str(i) for i in filtered]
                _ = out_file.write("\t".join(fields) + "\n")
                


parser = argparse.ArgumentParser("Process input")
parser.add_argument("-i", "--input", help = "input vcf file", type = str, required = True)
parser.add_argument("-b", "--bed", help = "input bed file with introgressed regions", type = str, required = True)
parser.add_argument("-p", "--ploidy", help = "ploidy of the introgressed fragments", type = str, required = True)
parser.add_argument("-o", "--output", help = "output text file", type = str, default = "output_ancestry.txt")

args = parser.parse_args()
input_vcf = args.input
input_bed = args.bed
output_txt = args.output
ploidy = args.ploidy

if __name__ == "__main__":

    # input_vcf = "test_simulation_0.vcf"
    # input_bed = "test_introgressed_fragments_sims.bed"
    # output_txt = "output_ancestry.txt"

    hap_dict = read_bedfile(input_bed)
    print(hap_dict.keys())
    print("Hap dict created")
    #print(hap_dict)
    find_ancestry(hap_dict, input_vcf, output_txt, ploidy)
