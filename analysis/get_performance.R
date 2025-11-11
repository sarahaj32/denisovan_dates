
library(tidyverse)
library(data.table)

# set plotting settings and variables
theme_set(theme_bw() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
            theme(axis.text = element_text(size = 10)))



# a helper function that expands the bedfile to a label for every 1000 bp
expand_rows <- function(start, end){
  # find the start and end, rounded to 1000
  # go by the majority
  start_1000 <- start - start %% 1000
  end_1000 <- end + (1000 - (end %% 1000))
  # a fix for the hmm
  end_1000 <- ifelse(end_1000 == start_1000, end_1000 + 1000, end_1000)
  
  pos_seq <- seq(from = start_1000, to = end_1000 - 1000, by = 1000)
  
  return(data.table(pos = pos_seq, pos2 = pos_seq + 1000))
}

# read in the true ancestry for a chromosome
frag_dat <- data.frame()
for (sim in c("Neanderthal_simple")){ 
  print(sim)
  for (rep in c("_1", "_2", "_3")){ #
    print(rep)
    for (arch in c("DEN", "NEA")){
      sim_frag_dat <- data.frame()
      for (i in 1:20){
        tmp_dat <- read.table(paste0("simdat/denisovan", sim, rep, "/", arch, "/fragments/denisovan", sim, "_sim_intro_", i, ".bed"))
        colnames(tmp_dat) <- c("hap", "chrom", "start", "end")
        tmp_dat$chrom <- paste0("chr", i)
        sim_frag_dat <- rbind(sim_frag_dat, tmp_dat)
      }
      sim_frag_dat$length <- sim_frag_dat$end - sim_frag_dat $start
      # calculate total for each chromosome and plot      
      sim_frag_dat$sim <- sim
      sim_frag_dat$rep <- rep
      sim_frag_dat$label <- arch
      frag_dat <- rbind(frag_dat, sim_frag_dat)
    }
  }
}

# add on the individual name
frag_dat <- frag_dat %>%  group_by(chrom, sim, rep, label) %>%
  mutate(id = cumsum(hap != lag(hap, default=first(hap)))) %>%
  mutate(id = paste0("NAMH_", ((id) %/% 2) + 1))  %>%
  ungroup()

# read in the inferred archaic fragments
inf_dat <- data.frame()
for (sim in c("Neanderthal_simple")){ # , "early", "late", "low" "early", 
  print(sim)
  for (rep in c("_1", "_2", "_3")){
  print(rep)
    sim_inf_dat <- data.frame()
    for (i in 1:20){
       tmp_dat <- read.table(paste0("simdat/denisovan", sim, rep, "/archaic_0.8/fragments_anno/denisovan", sim, "_hmmix.200_annotated_", i, ".bed"), header = T)
      colnames(tmp_dat)[colnames(tmp_dat) == "hapID"] <- "hap"
      tmp_dat$chrom <- paste0("chr", i)
      sim_inf_dat <- rbind(sim_inf_dat, tmp_dat)
    }
    
    sim_inf_dat$sim <- sim
    sim_inf_dat$rep <- rep
    colnames(sim_inf_dat)[colnames(sim_inf_dat) == "haplotype"] <- "id"
    inf_dat <- rbind(inf_dat, sim_inf_dat)
  }
}

# read in most-likely classification
inf_ML_dat <- data.frame()
for (sim in c("Neanderthal_simple")){ # , "early", "late", "low" "early", 
  print(sim)
  for (rep in c("_1", "_2", "_3")){
  print(rep)
    for (arch in c("DEN", "NEA")){
      sim_inf_ML_dat <- data.frame()
      for (i in 1:20){
         tmp_dat <- read.table(paste0("simdat/denisovan", sim, rep, "/", arch, "/fragments/denisovan", sim, "_hmmix.200_mostLikely_", i, ".bed"), header = T)
        colnames(tmp_dat)[colnames(tmp_dat) == "hapID"] <- "hap"
        tmp_dat$chrom <- paste0("chr", i)
        sim_inf_ML_dat <- rbind(sim_inf_ML_dat, tmp_dat)
      }
      
      sim_inf_ML_dat$sim <- sim
      sim_inf_ML_dat$rep <- rep
      sim_inf_ML_dat$label <- arch
      colnames(sim_inf_ML_dat)[colnames(sim_inf_ML_dat) == "haplotype"] <- "id"
      inf_ML_dat <- rbind(inf_ML_dat, sim_inf_ML_dat)
    }
  }
}

# read in ND classifications
inf_ND_dat <- data.frame()
for (sim in c("Neanderthal_simple")){ # , "early", "late", "low" "early", 
  print(sim)
  for (rep in c("_1", "_2", "_3")){
  print(rep)
    for (arch in c("DEN", "NEA")){
      sim_inf_ND_dat <- data.frame()
      for (i in 1:20){
         tmp_dat <- read.table(paste0("simdat/denisovan", sim, rep, "/", arch, "/fragments/denisovan", sim, "_hmmix.200_ND_", i, ".bed"), header = T)
        colnames(tmp_dat)[colnames(tmp_dat) == "hapID"] <- "hap"
        tmp_dat$chrom <- paste0("chr", i)
        sim_inf_ND_dat <- rbind(sim_inf_ND_dat, tmp_dat)
      }

      sim_inf_ND_dat$sim <- sim
      sim_inf_ND_dat$rep <- rep
      sim_inf_ND_dat$label <- arch
      colnames(sim_inf_ND_dat)[colnames(sim_inf_ND_dat) == "haplotype"] <- "id"

      inf_ND_dat <- rbind(inf_ND_dat, sim_inf_ND_dat)
    }
  }
}

# compare the inferred to the truth
results <- data.frame()
for (tmp_sim in c("Neanderthal_simple")){
  print(tmp_sim)
  for (tmp_ind in unique(inf_dat$id)){
    print(tmp_ind)
    # filter data to individual
    tmp_inf <- filter(inf_dat, id == tmp_ind, sim == tmp_sim)
    tmp_ML <- filter(inf_ML_dat, id == tmp_ind, sim == tmp_sim)
    tmp_frags <- filter(frag_dat, id == tmp_ind, sim == tmp_sim)
    tmp_ND <- filter(inf_ND_dat, id == tmp_ind, sim == tmp_sim)
    # rename the two haplotypes to "0" and "1"
    tmp_inf <- group_by(tmp_inf, chrom) %>% mutate(hap = dense_rank(hap) - 1)
    tmp_ML <- group_by(tmp_ML, chrom) %>% mutate(hap = dense_rank(hap) - 1)
    tmp_frags <- group_by(tmp_frags, chrom) %>% mutate(hap = dense_rank(hap) - 1)
    tmp_ND <- group_by(tmp_ND, chrom) %>% mutate(hap = dense_rank(hap) - 1)
    for (tmp_hap in c(0, 1)){
      print(tmp_hap)
      for (tmp_rep in c("_1", "_2", "_3")){
        hap_inf <- filter(tmp_inf, hap == tmp_hap, rep == tmp_rep)
        hap_frags <- filter(tmp_frags, hap == tmp_hap, rep == tmp_rep)
        hap_ND <- filter(tmp_ND, hap == tmp_hap, rep == tmp_rep)
        hap_ML <- filter(tmp_ML, hap == tmp_hap, rep == tmp_rep)
        # rename the classification titles so we can merge together
        hap_ML <- rename(hap_ML, "ML_label" = label, "inferred_ML_length" = length) %>%
          select(hap, chrom, start, end, inferred_ML_length, ML_label)
        hap_ND <- rename(hap_ND, "ND_label" = label, "inferred_ND_length" = length) %>%
          select(hap, chrom, start, end, inferred_ND_length, ND_label)

        expanded <- setDT(hap_inf)[, expand_rows(start, end), by = .(start, end)]
        expanded <- merge(expanded, hap_inf, by = c("start", "end"))  %>% dplyr::select(-c(start, end))
        expanded_ML <- setDT(hap_ML)[, expand_rows(start, end), by = .(start, end)]
        expanded_ML <- merge(expanded_ML, hap_ML, by = c("start", "end"))  %>% dplyr::select(-c(start, end))
        expanded_ND <- setDT(hap_ND)[, expand_rows(start, end), by = .(start, end)]
        expanded_ND <- merge(expanded_ND, hap_ND, by = c("start", "end"))  %>% dplyr::select(-c(start, end))
        expanded <- merge(expanded, expanded_ML, by = c("hap", "pos", "pos2", "chrom"), all.x = T) %>%
          mutate(ML_label = ifelse(is.na(ML_label), "unclassified", ML_label)) %>%
          rename("inferred_length" = length)
        
        expanded <- merge(expanded, expanded_ND, by = c("hap", "pos", "pos2", "chrom"), all.x = T) %>%
          mutate(ND_label = ifelse(is.na(ND_label) & ((ND01 + ND10) > 0), "mosaic", ND_label))
        
        expanded_true <- setDT(hap_frags)[, expand_rows(start, end), by = .(start, end)]
        expanded_true <- merge(expanded_true, hap_frags, by = c("start", "end")) %>% dplyr::select(-c(start, end, sim, rep)) %>%
          rename("true_label" = label)

        # merge and calculat the performance label of each category
        combo <- merge(expanded, expanded_true, all = T, by = c("pos", "pos2", "hap", "chrom", "id")) %>%
          mutate(arch_perf = ifelse(!is.na(state) & !is.na(true_label), "TP", 
                                ifelse(is.na(state) & !is.na(true_label), "FN", 
                                       ifelse(!is.na(state) & is.na(true_label), "FP", "Other")))) %>%
          mutate(ML_perf = ifelse(is.na(ML_label), "unclassified",
                   ifelse(ML_label == true_label, "TP", 
                                ifelse(ML_label == "unclassified", "unclassified", 
                                       ifelse(!is.na(ML_label) & !is.na(true_label) & (ML_label != true_label), "missclassified", "other"))))) %>%
          mutate(ND_perf = ifelse(is.na(ND_label), "unclassified", 
                              ifelse(ND_label == true_label, "TP", 
                                  ifelse(ND_label == "mosaic", "mosaic",
                                    ifelse(!is.na(ND_label) & !is.na(true_label) & (ND_label != true_label), "missclassified", "other")))))
        
        combo[is.na(combo)] <- ""
        combo %>% write.table(paste0("analysis/hmmix_comps/", tmp_sim, tmp_rep, "_", tmp_ind, "_", tmp_hap), quote = F, row.names = F, sep = "\t")
        tmp_results <- combo %>% group_by(arch_perf, ML_perf, ND_perf) %>% summarize("kb_blocks" = n())
        tmp_results$ind <- tmp_ind
        tmp_results$hap <- tmp_hap
        tmp_results$rep <- tmp_rep
        tmp_results$sim <- tmp_sim
        results <- rbind(results, tmp_results)
      }
    }
  }
  write.table(results, paste0("analysis/hmmix_comps/", tmp_sim, "_performanceResults.txt"), sep = "\t", quote = F, row.names = F)
}
