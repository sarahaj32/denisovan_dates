# a function that creates a geno-like file of ancestry at every position
# inputs: a vcf that is used to determine each position of interest. Note, this is only used to determine positions so no genotype info is necessary
#         a bed file with the archaic regions annotated. The first 4 positions must be haplotype, chrom, start, end
# output:
#         a text file with the ancestry label for each individual at each position. 0 = no archaic ancestry, 1 = heterozygous archaic ancestry, 2 = homozygous archaic ancestry

import argparse
from collections import defaultdict

def read_bedfile(bed):
    # read in the bedfile into a dictionary of haplotype to list of [start, end] positions
    i = 0
    hap_dict = defaultdict(list)
    with open(input_bed, 'r') as bed_file:
        hap_idx = 0
        start_idx = 1
        end_idx = 2
        for bed_line in bed_file:
            if "chrom" in bed_line:
                bed_line = bed_line.strip().split("\t")
                hap_idx = bed_line.index("hap")
                start_idx = bed_line.index("start")
                end_idx = bed_line.index("end")
                continue

            bed_line = bed_line.strip().split("\t")
            # remove "hap" from the haplotype column if it exists
            hap = int(bed_line[hap_idx].replace("hap", ""))
            if hap > i:
                i += 1
            if hap == i:
                start = int(bed_line[start_idx])
                end = int(bed_line[end_idx])
                hap_dict[i].append([start, end])
# as a final step, need to make sure the lists are sorted
    return(hap_dict)


def find_ancestry(hap_dict, input_vcf, output_txt) -> None:
    count = 0 
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
                    for ind in range(0, haps, 2):
                        if ind in hap_dict:
                            # print(hap_dict[ind][0])
                            while pos > hap_dict[ind][0][1]:
                                # print("TESTING")
                                # print(hap_dict)
                                # print(hap_dict[ind])
                                # print(hap_dict[ind][0])
                                # print(f"Position {pos} is greater than end of fragment {hap_dict[ind][0][1]}, removing fragment")
                                hap_dict[ind].pop(0)
                                if len(hap_dict[ind]) == 0:
                                    hap_dict[ind] = [[float('inf'), float('inf')]]
                                # print("NEW")
                                # print(hap_dict[ind])
                            hap_dict[ind][0][0]
                            if pos >= hap_dict[ind][0][0] and pos <= hap_dict[ind][0][1]:
                                #print(f"Position {pos} on chromosome {chrom} is in introgressed fragment in individual {ind} on haplotype 1")
                                ancestry_labs[ind] = 1
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
                    # print("ANCESTRY LABS:")
                    # print(ancestry_labs)
                # print(ancestry_labs)
                filtered = ancestry_labs[::2]
                # print(filtered)
                fields = [str(chrom), str(pos)] + [str(i) for i in filtered]
                _ = out_file.write("\t".join(fields) + "\n")
                


parser = argparse.ArgumentParser("Process input")
parser.add_argument("-i", "--input", help = "input vcf file", type = str, required = True)
parser.add_argument("-b", "--bed", help = "input bed file with introgressed regions", type = str, required = True)
parser.add_argument("-o", "--output", help = "output text file", type = str, default = "output_ancestry.txt")

args = parser.parse_args()
input_vcf = args.input
input_bed = args.bed
output_txt = args.output


if __name__ == "__main__":

    # input_vcf = "test_simulation_0.vcf"
    # input_bed = "test_introgressed_fragments_sims.bed"
    # output_txt = "output_ancestry.txt"

    hap_dict = read_bedfile(input_bed)
    print("Hap dict created")
    #print(hap_dict)
    find_ancestry(hap_dict, input_vcf, output_txt)
