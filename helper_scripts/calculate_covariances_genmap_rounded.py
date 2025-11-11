# read in the vcf:
# loop through every position, find all other positions within x distance and record ancestry
# then calculate covariance of those two vectors
# rewritten so that the start and end are always correctly rounded to 5 decimal places

import numpy as np
from collections import deque
import time
import argparse


def fast_cov(x, y):
    # Sample covariance, matching np.cov(x, y, bias=False)[0,1]
    n = len(x)
    mean_x = np.mean(x)
    mean_y = np.mean(y)
    return (np.dot(x, y) - n * mean_x * mean_y) / (n - 1)


# thanks chat gpt, this is a faster way to do many covariances at once
# i confirmed it gives the same results as looping through each row
def cov_with_many(x, Y):
    """
    Compute covariance between 1D vector x (length m)
    and each row of 2D matrix Y (shape k x m).
    Returns array of length k.
    Matches np.cov with sample covariance (N-1 in denom).
    """
    n = x.shape[0]
    mean_x = x.mean()
    mean_Y = Y.mean(axis=1)

    # dot product of x with each row of Y
    dot_xy = Y @ x  

    return (dot_xy - n * mean_x * mean_Y) / (n - 1)
    
def compute_anc_D(input_file, my_min, my_max): # output_file, 
    current_window = deque()
    count = 0
    D_sum = 0
    comps = 0
    with open(input_file, "r") as f: # , open(output_file, "w") as o
        for line in f:
            # print("NEW LINE")
            count += 1
            # if count >= 100000:
            #     break
            if count % 200000 == 0:
                print(f"On line {count}, with {comps} comparisons made so far and current D: {D_sum / comps if comps > 0 else 0}")
            line = line.strip().split("\t")
            chrom = line[0]
            pos = float(line[1])
            # all of the ancestry assignments
            ass = np.array(line[2:], dtype=int)
            # append the current position to the window
            current_window.append((pos, ass))
            # print(current_window)
            # when the first element is too far from the last, process first element
            while current_window[-1][0] - current_window[0][0] > my_max:
                # print("readjusting")
                # print(current_window)
                first_pos, first_ass = current_window[0]
                # compare first element to all other elements except itself and the last
                # for other_pos, other_ass in list(current_window)[1:-1]:
                #     if other_pos - first_pos >= my_min:
                #         # print("OUTPUT:")
                #         # print(f"{first_ass}\t{other_ass}\n")
                #         p_cov = fast_cov(first_ass, other_ass)
                #         D_sum += p_cov
                #         comps += 1
                #         _ = o.write("\t".join([str(first_pos), str(other_pos), str(other_pos - first_pos), str(p_cov)]) + "\n")
                #         # if p_cov[0,1] != 0:
                #         # print(f"{first_pos}:{other_pos} - covariance: {p_cov}")
                #     # else:
                #     #     print(f"Not far enough away: {first_pos}:{other_pos}")

                # vectorized version thanks to chat gpt (I confirmed it gives the same results - do the math though!!!!)
                # stack ancestry vectors of "other" positions into array
                others = [ass for pos, ass in list(current_window)[1:-1] if pos - first_pos >= my_min]
                other_pos = [pos for pos, ass in list(current_window)[1:-1] if pos - first_pos >= my_min]

                if others:
                    Y = np.vstack(others)  # shape (k, m)
                    covs = cov_with_many(first_ass, Y)

                    for pos, c in zip(other_pos, covs):
                        D_sum += c
                        comps += 1
                        #_ = o.write("\t".join([str(first_pos), str(pos), str(pos - first_pos), str(c)]) + "\n")

                # remove first element
                current_window.popleft()
                # print(current_window)

    print(f"D_sum: {D_sum}, comps: {comps}, pos: {count}")
    return(D_sum, comps, count)


parser = argparse.ArgumentParser("Process input")
parser.add_argument("-i", "--input", help = "input anc file", type = str, required = True)
parser.add_argument("-o", "--output", help = "output text file", type = str, default = "output_ancestry.txt")
parser.add_argument("-s", "--start", help = "starting point to continue analysis of incomplete results", type = float, default = None)

args = parser.parse_args()
input_file = args.input
output_file = args.output
restart = args.start

if restart:
    begin = float(restart)
else:
    begin = 0.0001
    
if __name__ == "__main__":
    # input_file = "simdat/denisovan_simple/curve/fragments/denisovan_simple_sim_intro_16.anc"
    # output_file = "test_covariance.txt"
    # int_path = "testing/simulate_papuans/intermediate_covs_16"
    num = round((0.01 - begin) / 0.00001)
    values = np.round(np.linspace(begin, 0.01, num=num), 5)  # 990 steps ≈ 0.00001 spacing

    if restart:
        print(f"Restarting analysis, appending to {output_file}")
        fo = open(output_file, "a")
    else:
        fo = open(output_file, "w")
    for my_min in values:
        my_max = round(my_min + 0.00001, 5)
        print("NEXT!!!")
        print(f"Calculating D for min {my_min} and max {my_max}")
        start = time.time()
        D_sum, comps, pos = compute_anc_D(input_file, my_min, my_max) # int_out,
        D = D_sum / comps
        end = time.time()
        print("COMPLETE")
        print(f"Elapsed time: {end - start:.2f} seconds")
        print(f"{my_min}\t{my_max}\t{D}\t{D_sum}\t{comps}\t{pos}\n")
        fo.write(f"{my_min}\t{my_max}\t{D}\t{D_sum}\t{comps}\t{pos}\n")
        fo.flush()
    fo.close()
    # else:
    #     with open(output_file, "w") as fo:
    #         values = np.linspace(0.0001, 0.01, num=990)  # 990 steps ≈ 0.00001 spacing

    #         for my_min in values:
                
    #             #int_out = f"{int_path}_{my_min}.txt"
    #             my_max = round(my_min + 0.00001, 5)
    #             print("NEXT!!!")
    #             print(f"Calculating D for min {my_min} and max {my_max}")
    #             start = time.time()
    #             D_sum, comps, pos = compute_anc_D(input_file,  my_min, my_max) # int_out,
    #             D = D_sum / comps
    #             end = time.time()
    #             print("COMPLETE")
    #             print(f"Elapsed time: {end - start:.2f} seconds")
    #             print(f"{my_min}\t{my_max}\t{D}\t{D_sum}\t{comps}\t{pos}\n")
    #             fo.write(f"{my_min}\t{my_max}\t{D}\t{D_sum}\t{comps}\t{pos}\n")
    #             fo.flush()