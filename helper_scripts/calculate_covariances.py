# read in the vcf:
# loop through every position, find all other positions within x distance and record ancestry
# then calculate covariance of those two vectors

import numpy as np
from collections import deque
import time
import argparse
# import os
# import shutil

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
    count = 0
    D_sum = 0
    comps = 0
    pos_window = deque()
    ass_window = deque() 
    
    max_window = 1000
    Y = None
    tmp_pos = np.empty(max_window, dtype=np.int32)

    
    with open(input_file, "r") as f: # , open(output_file, "w") as o
        for line in f:
            # print("NEW LINE")
            count += 1
            # if count >= 100000:
            #     break
            # if count % 200000 == 0:
            #     print(f"On line {count}, with {comps} comparisons made so far and current D: {D_sum / comps if comps > 0 else 0}")
            fields = line.strip().split("\t")
            chrom = fields[0]
            pos = int(fields[1])
            # all of the ancestry assignments
            ass = np.array(fields[2:], dtype=np.int8)
            # append the current position to the window
            pos_window.append(pos)
            ass_window.append(ass)
            # print(current_window)
            # when the first element is too far from the last, process first element
            while pos_window[-1] - pos_window[0] > my_max:
                # compare first element to all other elements except itself and the last

                first_pos = pos_window[0]
                first_ass = ass_window[0]

                # vectorized version thanks to chat gpt (I confirmed it gives the same results - do the math though!!!!)
                # stack ancestry vectors of "other" positions into array
                others = []
                other_pos = []

                if Y is None:
                    m = first_ass.shape[0]
                    Y = np.empty((max_window, m), dtype = np.int8)
                k = 0

                for i in range(1, len(pos_window)-1):
                    if pos_window[i] - first_pos >= my_min:
                        tmp_pos[k] = pos_window[i]
                        Y[k] = ass_window[i]
                        k += 1
                if k > 0:
                    covs = cov_with_many(first_ass, Y[:k])
                    D_sum += covs.sum()
                    comps += k                 

                        #_ = o.write("\t".join([str(first_pos), str(pos), str(pos - first_pos), str(c)]) + "\n")

                # remove first element
                pos_window.popleft()
                ass_window.popleft()

    print(f"D_sum: {D_sum}, comps: {comps}, pos: {count}")
    return(D_sum, comps, count)


parser = argparse.ArgumentParser("Process input")
parser.add_argument("-i", "--input", help = "input anc file", type = str, required = True)
parser.add_argument("-o", "--output", help = "output text file", type = str, default = "output_ancestry.txt")

args = parser.parse_args()
input_file = args.input
output_file = args.output

if __name__ == "__main__":
    # input_file = "previous_simdat/simdat/denisovan_simple/curve/fragments/denisovan_simple_sim_intro_16.anc"
    # output_file = "test_covariance.txt"
    # int_path = "testing/simulate_papuans/intermediate_covs_16"

    # tmp_dir = f"/tmp/{os.environ.get('USER')}/myjob_tmp"
    # os.makedirs(tmp_dir, exist_ok=True)
    # tmp_file = os.path.join(tmp_dir, "output.tmp")
    # tmp_file = os.path.join(tmp_dir, f"{os.path.basename(output_file)}.{os.getpid()}.tmp")


    with open(output_file, "w") as fo:
        for my_min in range(10000, 1000000, 1000):
            
            #int_out = f"{int_path}_{my_min}.txt"
            my_max = my_min + 1000
            #print("NEXT!!!")
            #print(f"Calculating D for min {my_min} and max {my_max}")
            start = time.time()
            D_sum, comps, pos = compute_anc_D(input_file, my_min, my_max) # int_out,
            D = D_sum / comps
            end = time.time()
            #print("COMPLETE")
                        # print occasionally, not every iteration
            if my_min % 10000 == 0:
                print(f"Elapsed: {end - start:.2f}s for {my_min}-{my_max}", flush=True)

            #print(f"{my_min}\t{my_max}\t{D}\t{D_sum}\t{comps}\t{pos}\n")
            fo.write(f"{my_min}\t{my_max}\t{D}\t{D_sum}\t{comps}\t{pos}\n")
            fo.flush()
    
    #shutil.move(tmp_file, output_file)