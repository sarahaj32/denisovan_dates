



# goal output: 
# chromosome, position, reference, alternative, ancestral, derived, introgressed, altai, vindija, chagyrskaya, denisovan, label (ND01_all, ND10_all, ND00_all, ND11_all), comparable (T/F)
# loop through the file of variants instead - label each one, and save that to a file
import argparse
import sys

print("starting to annotate introgressed fragments")

# first read in the label for each observation into a dictionary
def read_in_pos_label(obs, ancestral):
    with open(obs) as file:
        obs_dict = {}
        count = 0
        for line in file:
            line = line.strip().split("\t")
            if "chrom" in line:
                chrom_idx = line.index("chrom")
                pos_idx = line.index("pos")
                label_idx = line.index("label")
                olabel_idx = line.index("original_label")
                if ancestral:
                    anc_idx = line.index(ancestral)
                ref_idx = line.index("ref")
                alt_idx = line.index("alt")
            else:
                ID = line[chrom_idx] + ":" + str(line[pos_idx])
                label = line[label_idx]
                olabel = line[olabel_idx]
                if ancestral:
                    anc = line[anc_idx]
                else:
                    anc = "_"
                ref = line[ref_idx]
                alt = line[alt_idx]
                dat = [ref, alt, anc, label, olabel]
                obs_dict[ID] = dat
    print(f"observations read in, with length: {len(obs_dict)}")
    return(obs_dict)

def annotate_fragments(obsdict, frags, outfile):
count = 0
with open(frags, "r") as af, open(outfile, "w") as of:
    for line in af:
        count += 1
        # print out a status update every so often
        if count % 100000 == 0:
            print(count)
        # if count > 200:
        #     break
        line = line.strip().split("\t")
        if "name" in line:
            variant_idx = line.index("variants")
            chrom_idx = line.index('chrom')
            outline = "\t".join(line + ["comp", "nocomp", "ND", "ND10", "ND01", "ND11", "ND00", "nea_overlap", "den_overlap", "nea_affinity", "den_affinity"])
            print("NAME")
        elif count == 1:
            chrom_idx = 1
            variant_idx = 15
            outline = "\t".join(line + ["comp", "nocomp", "ND", "ND10", "ND01", "ND11", "ND00", "nea_overlap", "den_overlap", "nea_affinity", "den_affinity"])
            print("FIRST LINE")
        else:
            label_count = {"ND01":0, "ND10":0, "ND11": 0, "ND00": 0, "ND":0, "comp": 0, "nocomp": 0}
            vars = line[variant_idx].split(sep = ",")
            chrom = line[chrom_idx]
            for var in vars:
                ID = f"{chrom}:{var}"
                #see if this ID is comparable to archaics
                if ID in obs_dict: 
                    ref, alt, anc, label, olabel = obs_dict[ID]
                    # count at all comparable sites (include multi-allelic for now)
                    label_count[label] += 1
                    if label != "ND":
                        label_count["comp"] += 1
                    # also get a count of how many are 'comparable'
                else:
                    label_count["nocomp"] += 1
            # calculate the number of matches to each archaic, and the overall affinity metric
            label_count["nea_overlap"] = label_count["ND11"] + label_count["ND10"]
            label_count["den_overlap"] = label_count["ND11"] + label_count["ND01"]
            if label_count["nea_overlap"] > 0:
                label_count["nea_affinity"] = label_count["nea_overlap"] / label_count["comp"]
            else:
                label_count["nea_affinity"] = 0
            if label_count["den_overlap"] > 0:
                label_count["den_affinity"] = label_count["den_overlap"] / label_count["comp"]
            else:
                label_count["den_affinity"] = 0
            outline = "\t".join(line + [str(label_count[i]) for i in ["comp", "nocomp", "ND", "ND10", "ND01", "ND11", "ND00", "nea_overlap", "den_overlap", "nea_affinity", "den_affinity"]])
        _ = of.write(outline + "\n")


                
                
parser = argparse.ArgumentParser()
parser.add_argument("-b", "--bed", help="bed file with all fragments")
parser.add_argument("-i", "--input", help="file with annotations")
parser.add_argument("-o", "--output", help="output file")
parser.add_argument("-a", "--anc_reference", help="name of individual in vcf to use as ancestral", default = None)

args = parser.parse_args()

frags = args.bed
obs = args.input
outfile = args.output
anc = args.anc_reference
print(anc)

print(f"writing to outfile: {outfile}")
if __name__ == "__main__":
    print(f"Using {anc} as ancestral")
    obs_dict = read_in_pos_label(obs, anc)
    print("labels read in, annotating fragments")
    annotate_fragments(obs_dict, frags, outfile)

#frags = "data/1000G/1000G_fragments_snps10_08.txt"
#obs = "data/1000G_snps10_08_all_observations_annotated.txt"
#outfile = f"data/1000G/summary_snps10_08_all.txt"
# anc = "ancestral" or "panTro5"

