import argparse
import random

def allele_from_gt(ref, alt, gt):
    if int(gt) == 1:
        allele = alt
    elif int(gt) == 0:
        allele = ref
    return(allele)

def get_nd(nea, den, der, rand=False):
    nea = str(nea)
    den = str(den)
    der = str(der)
    nea_der = nea.count(der) / len(nea)
    den_der = den.count(der) / len(den)
    if nea_der == den_der and nea_der == "1":
        return("ND11")
    elif nea_der == den_der and nea_der == "0":
        return("ND00")
    if rand:
        if nea_der == 0.5:
            print("NEA RANDOME")
            nea_anno = random.randint(0, 1)
            #print(nea_anno)
        else:
            assert(nea_der == 0 or nea_der == 1)
            nea_anno = int(nea_der)
        if den_der == 0.5:
            print("DEN RANDOM")
            den_anno = random.randint(0, 1)
            #print(den_anno)
        else:
            assert(den_der == 0 or den_der == 1)
            den_anno = int(den_der)

        label = f"ND{nea_anno}{den_anno}"
        #print(label)
    else:
        nea_anno = 1 if nea_der > 0 else 0
        den_anno = 1 if den_der > 0 else 0
        label = f"ND{str(nea_anno)}{str(den_anno)}"
    return(label)


def annotate_archaics(infile, outfile, anc):
    count = 0

    # read in the vcf
    with open(infile, "r") as f, open(outfile, "w") as of:
        of.write("\t".join(["chrom", "pos", "ref", "alt", "derived", "altai", "vindija", "chagyrskaya", "denisova", "panTro5", "ancestral", "label", "original_label"]))
        for line in f:
            count += 1
            if count > 200:
                break
            # parse the header
            if line.startswith("##"):
                continue
            elif line.startswith("#CHROM"):
                headers = line.strip().split("\t")
                chrom_idx = headers.index("#CHROM")
                pos_idx = headers.index("POS")
                ref_idx = headers.index("REF")
                alt_idx = headers.index("ALT")
                vin_idx = headers.index("Vindija33.19")
                altai_idx = headers.index("AltaiNeandertal")
                chag_idx = headers.index("Chagyrskaya-Phalanx")
                den_idx = headers.index("Denisova")
                if anc: 
                    anc_idx = headers.index(anc)
            # parse the data lines
            else:
                #print(line)
                info = line.strip().split("\t")
                ref_allele = info[ref_idx]
                alt_allele = info[alt_idx]
                # find the derived allele:
                # assumes that if an ancestral individual is provided, they are homozygous for the ancestral allele
                der_allele = alt_allele
                if anc:
                    anc_allele = allele_from_gt(ref_allele, alt_allele, line[anc_idx][0])
                    if anc_allele == alt_allele:
                        der_allele = ref_allele
                #print(f"der allele: {der_allele}")
                vin_gt = "".join([allele_from_gt(ref_allele, alt_allele, info[vin_idx][0]), allele_from_gt(ref_allele, alt_allele, info[vin_idx][2])])
                chag_gt = "".join([allele_from_gt(ref_allele, alt_allele, info[chag_idx][0]), allele_from_gt(ref_allele, alt_allele, info[chag_idx][2])])
                altai_gt = "".join([allele_from_gt(ref_allele, alt_allele, info[altai_idx][0]), allele_from_gt(ref_allele, alt_allele, info[altai_idx][2])])
                den_gt = "".join([allele_from_gt(ref_allele, alt_allele, info[den_idx][0]), allele_from_gt(ref_allele, alt_allele, info[den_idx][2])])
                # print(f"vin_gt: {vin_gt}, chag_gt: {chag_gt}, altai_gt: {altai_gt}, den_gt: {den_gt}")
                if not anc:
                    pan_gt = ref_allele
                    anc_indiv_gt = ref_allele
                OG_label = get_nd(vin_gt, den_gt, der_allele, rand=True)
                label = get_nd("".join([altai_gt, vin_gt, chag_gt]), den_gt, der_allele, rand=False)
                # print("LABEL:")
                # print(label)
                # print("OG_LABEL:")
                # print(OG_label)
                outline = "\t".join([info[chrom_idx], info[pos_idx], info[ref_idx], info[alt_idx], der_allele, altai_gt, vin_gt, chag_gt, den_gt, pan_gt, anc_indiv_gt, label, OG_label])
                print(outline)



parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="input archaic vcf")
parser.add_argument("-o", "--output", help="output file")
parser.add_argument("-a", "--anc_reference", help="name of individual in vcf to use as ancestral", default = None)

args = parser.parse_args()

infile = args.input
outfile = args.output
anc = args.anc_reference

print(f"writing to outfile: {outfile}")
if __name__ == "__main__":
    print(f"Using {anc} as ancestral")
    _ = annotate_archaics(infile, outfile, anc)

