import argparse
from pyfaidx import Fasta

print("packages installed, reading in parameters")

parser = argparse.ArgumentParser()
parser.add_argument("-v", "--vcf", help="name of the vcf file", type = str)
parser.add_argument("-o", "--output", help="name of the output vcf file", type = str)
parser.add_argument("-r", "--reference", help="path to the reference fasta file", type = str, default = None)
parser.add_argument("-p", "--positions", help="file with all desired positions", type = str)
parser.add_argument("-i", "--intersect", help="only keep vcf positions that are also in the positions file", action = "store_true")
args = parser.parse_args()

print(f"parameters: {args}")
# vcf_file="data/1000G/ND_VCF/1000G_NDpos_refFilled_hg19_chr4.vcf.gz"
vcf_file = args.vcf
# positions_file="data/hg19/snp_annotations/ND_pos_hg19.txt"
positions_file = args.positions
# reference_fasta="/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/resources/hg19/fasta/hg19.fa"
reference_fasta = args.reference
# outfile = "test.vcf"
outfile = args.output
# if the -i flag is not included, all positions in the vcf will be output
remove = args.intersect

print('reading in reference fasta')
# read in fasta file
if reference_fasta:
    ref = Fasta(reference_fasta)
else:
    ref = None

count = 0

def advance(positions):
    try:
        pos = next(positions).strip().split("\t")
        pos_pos = int(pos[1])
        pos_chr = int(pos[0])
    except StopIteration:
        print("END of positions")
        pos_pos, pos_chr = None, None
    return(pos_pos, pos_chr) 

with open(positions_file, "r") as positions, open(vcf_file, "r") as vcf, open(outfile, "w") as of:
    print('starting to parse through')
    # get the first desired position 
    pos = next(positions).strip().split("\t")
    pos_chr = int(pos[0])
    pos_pos = int(pos[1])
    print(pos)
    # start reading the vcf line by line
    for line in vcf:
        count += 1
        # print out a status update every so often
        if count % 10000 == 0:
             print(f"{count} lines processed")
        # write out the full header
        if line.startswith("##"):
            of.write(line)
            # add on the header for RefFill to denote when we added a line with the reference
            if line.startswith("##INFO"):
                of.write("##INFO=<ID=RefFill,Number=.,Type=String,Description=\"Reference Allele Filled In.\">\n")
        # get the index of chromosome and position, and count how many individuals are present 
        # (we will fill in the homozygous reference for all of these)
        elif line.startswith("#CHROM"):
            headers = line.strip().split("\t")
            chrom_idx = headers.index("#CHROM")
            pos_idx = headers.index("POS")
            exclude = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
            indivs = [i for i in headers if i not in exclude]
            nindivs = len(indivs)
            of.write("\t".join(headers) + "\n")

        # all data lines:
        else:
            line = line.strip().split("\t")
            vcf_pos = int(line[pos_idx])
            vcf_chr = int(line[chrom_idx])
            # the first time around, make sure that we found the right chromosome in the positions file
            # advance the desired position of interest until we get to the correct chromosome
            while pos_chr and (pos_chr < vcf_chr):
                pos_pos, pos_chr = advance(positions)
            # check if it's the same chromosome
            if pos_chr and (pos_chr == vcf_chr):
                # if the vcf line is past the position of interest,
                # fill in all positions of interest until we get to the vcf position
                # print("moving forward")
                # print(f"position: {pos[1]}, vcf: {vcf_pos}")
                while pos_pos and (pos_pos < vcf_pos):
                    # print(f"position: {pos[1]}, vcf: {line[pos_idx]}")
                    # insert reference homozygous genotype for all individuals 
                    indivs_out = ["0/0"] * nindivs
                    # find the reference allele from the reference genome
                    if ref:
                        ref_allele = str(ref[f"chr{pos_chr}"][pos_pos - 1]).upper()
                    else:
                        ref_allele = "0"
                    # write out the full line (including our info flag)
                    included = [str(pos_chr), str(pos_pos), ".", ref_allele, ".", ".", ".", "RefFill", "GT"]
                    outline = included + indivs_out
                    of.write("\t".join(outline) + "\n")
                    # stop when we have reached the end of positions of interest
                    pos_pos, pos_chr = advance(positions)

                # write out any rows in the vcf that are not in the positions file 
                # (these can be filtered out later with the positions file)
                if pos_pos and (pos_pos > vcf_pos):
                    if remove:
                        pass
                    else:
                        of.write("\t".join(line) + "\n")
                # if the vcf line is the current position of interest, 
                # write it out and go to the next position of interest
                elif pos_pos == vcf_pos:

                    outline = line
                    of.write("\t".join(line) + "\n")

                    pos_pos, pos_chr = advance(positions) 
                # we should have caught every situation above, but add a catch here just incase
                elif pos_chr == None:
                    if remove:
                        pass
                    else:
                        of.write("\t".join(line) + "\n")
                else:
                    print("PROBLEM HERE:")
                    print(f"position: {pos_pos}, vcf: {vcf_pos}")
            # if we loop through all the positions of interest for this chromosome, print the rest of the vcf
            elif pos_chr == None:
                if remove:
                    pass
                else:
                    of.write("\t".join(line) + "\n")
            else:
                print("NO CHROM MATCH")
    # after we have looped through the entire VCF, add any remaining positions from the positions file:
    pos_pos, pos_chr = advance(positions)   
    while pos_chr == vcf_chr and pos_chr != None:
        indivs_out = ["0/0"] * nindivs
        # find the reference allele from the reference genome
        if ref:
            ref_allele = str(ref[f"chr{pos_chr}"][pos_pos - 1]).upper()
        else:
            ref_allele = "0"
        # write out the full line (including our info flag)
        included = [str(pos_chr), str(pos_pos), ".", ref_allele, ".", ".", ".", "RefFill", "GT"]
        outline = included + indivs_out
        of.write("\t".join(outline) + "\n")
        pos_pos, pos_chr = advance(positions) 

