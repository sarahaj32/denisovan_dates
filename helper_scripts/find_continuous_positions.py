import argparse
import re


argparser = argparse.ArgumentParser()
argparser.add_argument("-s", "--snp", type=str)
argparser.add_argument("-n", "--no_prob", type=str, default = None)
argparser.add_argument("-p", "--prob", type=str)
args = argparser.parse_args()

snp = args.snp
no_prob = args.no_prob
prob = args.prob   

# snp="../../data/Archaics/hg19/AltDenVinPanHSanc_merged_strictArchaic_hg19_chr2_genmap.snp"
# no_prob = "noprob.snp"
# prob = "prob.snp"

# intialize first position, any positive position after this will be valid
prev_pos = 0
prev_line = [""]
count = 0
if no_prob:
    npr = open(no_prob, "w")

with open(snp) as f,open(prob, "w") as p:
    for line in f:
        count += 1  
        if line.startswith("#") | line.startswith("P"):
            continue
        fields = re.split(r'[ \t]+', line.strip())
        fields = [x for x in fields if x != '']
        map_pos = float(fields[2])

        # include all positions that are greater than or equal to the previous position
        if map_pos >= prev_pos:
            prev_pos = map_pos
            prev_line = line
            if npr:
                npr.write(line)
        # write out all problem positions where the genetic map position decreases
        else:
            p.write("\t".join(fields) + "\n")
if no_prob:
    npr.close()

