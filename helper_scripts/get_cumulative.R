# I have confirmed that this gives me the same results as bedtools merge!!!!

#to run: 
# Rscript <file.bed> <Denisova> <output.bed>

# a script that caluclates the cumulative plots
print("starting, loading packages")
library(tidyverse)
args = commandArgs(trailingOnly=TRUE)

# read in the fragments file and desired archaic from command line arguments:
fragment_file <- args[1]
#fragment_file <- "qc_checks/synthetic_EAS/EAS_fragments_snps10_08_Dhigh.bed"
archaic <- args[2]
#archaic <- "Denisova"
outfile <- args[3]
#outfile <- "test.output"

print(fragment_file)
print(archaic)
print(outfile)
# helper functions

## combine fragments:
combine_frags <- function(dat){
  # first remove fragments that are identical across haplotypes
  dat <- dat %>% select(chrom, start, end) %>% arrange(chrom, start, end) %>% unique()
  
  # create an empty dataframe to hold the combined results
  new_dat <- data.frame()
  
  # merge overlapping fragments across the two haplotypes (chromosomes)
  start <- dat$start[1]
  end <- dat$end[1]
  if(nrow(dat) > 1){  
    for (i in 2:nrow(dat)){
    if (end > dat$start[i]){
      end <- max(end, dat$end[i])
      # move on to next fragment
    }else{
      # Add the non-overlapping fragment (either individual or merged) to new dataframe
      new_frag <- data.frame("start" = start, "end" = end)
      new_dat <- rbind(new_dat, new_frag)
      # move on to next fragment, reassign starts and ends to the next fragment
      start <- dat$start[i]
      end <- dat$end[i]
      }
    }
  }

  new_frag <- data.frame("start" = start, "end" = end)
  new_dat <- rbind(new_dat, new_frag)
  return(new_dat)
}

# find the merged length for selected individuals
Archaic_segments <- function(dataset,inds,archaic){
  data_tmp <- filter(dataset,individual %in% inds & ND_type==archaic)
  amount <- 0

  for(tmp_chrom in unique(data_tmp$chrom)){
    data_tmp_chrom <- filter(data_tmp, chrom == tmp_chrom)
    combo_tmp <- combine_frags(data_tmp_chrom)
    chrom_amount <- sum(combo_tmp$end - combo_tmp$start)
    amount <- amount + chrom_amount
  }
  return(amount)
}

# a helper function to add classes to the fragments
add_classes <- function(dat){
  dat$class_1 <- ifelse(dat$nea_overlap > dat$den_overlap, "Neanderthal", 
                       ifelse(dat$nea_overlap < dat$den_overlap, "Denisova", 
                              ifelse(dat$nea_overlap == 0 & dat$den_overlap ==  0, "None",
                                     ifelse(dat$nea_overlap == dat$den_overlap, "Both", "Problem"))))
  
  print("CLASS 1 Done")
  dat %>% group_by(class_1) %>% summarize(n()) %>% print()

  return(dat)
}

############
### MAIN ###
############

# read in file
frags <- read.table(fragment_file, header = T)
print(colnames(frags))
# make sure that the headers are consistent
if(!"individual" %in% colnames(frags)){
  if("name" %in% colnames(frags)){
    colnames(frags)[colnames(frags) == "name"] <- "individual"
  }else if("ID" %in% colnames(frags)){
    print("ID")
    colnames(frags)[colnames(frags) == "ID"] <- "individual"
  }else{
    print("neither individual or name in header")
  }

  if(!"ND_type" %in% colnames(frags)){
    print("classifying fragments")
    frags <- add_classes(frags)
    print("DONE")
    colnames(frags)[colnames(frags) == "class_1"] <- "ND_type"
  }
}

# create a list of all the individuals
inds <- select(frags, individual) %>% unique() %>% pull() %>% sample(replace = F)
# randomly shuffle:
inds <- sample(inds)
print("total number of individuals: ")
print(length(inds))

results <- rep(0, length(inds))
for(i in 1:length(inds)){
  print(i)
  result <- Archaic_segments(frags, inds[1:i], archaic)
  print(result)
  results[i] <- result  
}

output <- data.frame("ind" = 1:length(inds), amount = results)
write.table(output, outfile, quote = F, row.names = F, sep = "\t")

