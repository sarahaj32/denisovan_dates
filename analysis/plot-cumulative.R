###
## Archaic accumulation
####
rm(list = ls())
library(tidyverse)

###
## Choin - all samples
###

choin_cumulative <- read.table("/global/scratch/p2p3/pl1_moorjani/sarahj32/superarchaic_introgression/Choin_analysis/choin_all_sequential_cumulative_Densiovan.tsv", header = T)
# since to match laurit's file which is by haplotype (not individual like mine is)
choin_cumulative$ind_number <- choin_cumulative$ind_number * 2
choin_cumulative$ind_frac <- choin_cumulative$ind_number / max(choin_cumulative$ind_number)
# make sure that this looks right
pdf("plots/choin_combined_all.pdf", height = 5, width = 7)
ggplot(choin_cumulative, aes(x = ind_number, y = Denisova, color = super_population)) +
  geom_line() +
  facet_wrap(~dataset, scale = "free_x") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
dev.off()

# read in laurit's results from other datasets (so we can plot together)
ww_cumulative <- read.table("/global/scratch/p2p3/pl1_moorjani/sarahj32/superarchaic_introgression/Choin_analysis/Final_out_from_Laurits.txt", sep = "\t", header = T)

ww_cumulative <- group_by(ww_cumulative, dataset, super_population) %>% mutate(ind_frac = ind_number / max(ind_number)) %>% ungroup()

# combine with choin
ww_cumulative <- choin_cumulative %>% select(ind_number, ind_frac, Denisova, Neandertal, super_population, dataset) %>% rbind(select(ww_cumulative, ind_number, ind_frac, Denisova, Neandertal, super_population, dataset))
ww_cumulative <- ww_cumulative %>% pivot_longer(cols = c(Denisova, Neandertal))

# now plot decode, choin, and lasi-dad together (select colors for priyas grant)
cols <- c(
  "LASIDAD" = "#FF4500",     # dark orange
  "deCODE" = "#00008B",      # dark blue
  "choin" = "#117733",       # dark green
  "HGDP" = "#117733"         # same as choin
)

label_data_frac <- data.frame(
  name = c("Neandertal", "Neandertal", "Neandertal", "Denisova", "Denisova", "Denisova"),
  dataset = c("LASIDAD", "deCODE", "choin", "LASIDAD", "deCODE", "choin"),
  label = c("LASI–DAD\n(n=2,679)", "deCODE\n(n=27,566)", "choin\n(n=325)",
            "LASI–DAD\n(n=2,679)", "deCODE\n(n=27,566)", "choin\n(n=325)"),
  x = rep(0.85, 6),       # adjust x position as needed
  y = c(1.5, 1.2, 0.875, 0.75, 0.275, 1.15)          # adjust y position as needed
)
label_data_frac$color <- cols[label_data_frac$dataset]


pdf("plots/choin_decode_LASI_NeaDen_all_frac.pdf", height = 7, width = 5)
ww_cumulative %>% filter(dataset %in% c("choin", "deCODE", "LASIDAD")) %>%
  mutate(name = factor(name, levels = c("Neandertal", "Denisova"), ordered = T)) %>%
  ggplot(aes(x = ind_frac, y = value / 1000, color = dataset)) +
  geom_line() +
  theme_bw() +
  scale_color_manual(values = cols, guide = "none") +
  labs(y = "Total Sequence (Gb)", x = "Proportion of Data Sequenced") +
  geom_text(data = label_data_frac,
            aes(x = x, y = y, label = label, color = dataset),
            hjust = 0, vjust = 1, size = 3, show.legend = FALSE) +
  facet_grid(as.factor(name)~., scales = "free")
dev.off()

label_data_ind <- data.frame(
  name = c("Neandertal", "Neandertal", "Neandertal", "Denisova", "Denisova", "Denisova"),
  dataset = c("LASIDAD", "deCODE", "choin", "LASIDAD", "deCODE", "choin"),
  label = c("LASI–DAD\n(n=2,679)", "deCODE\n(n=27,566)", "Choin\n(n=325)",
            "LASI–DAD\n(n=2,679)", "deCODE\n(n=27,566)", "Choin\n(n=325)"),
  x = c(4000, 45000, 1500, 4000, 45000, 1500),       # adjust x position as needed
  y = c(1.45, 1.2, 1, 0.75, 0.275, 1.175)           # adjust y position as needed
)
label_data_ind$color <- cols[label_data_ind$dataset]


pdf("plots/choin_decode_LASI_NeaDen_all_num.pdf", height = 7, width = 5)
ww_cumulative %>% filter(dataset %in% c("choin", "deCODE", "LASIDAD")) %>%
  mutate(name = factor(name, levels = c("Neandertal", "Denisova"), ordered = T)) %>%
  ggplot(aes(x = ind_number, y = value / 1000, color = dataset)) +
  geom_line() +
  theme_bw() +
  scale_color_manual(values = cols, guide = "none") +
  labs(y = "Total Sequence (Gb)", x = "# Sequenced Individuals") +
  geom_text(data = label_data_ind,
          aes(x = x, y = y, label = label, color = dataset),
          hjust = 0, vjust = 1, size = 3, show.legend = FALSE) +
  facet_grid(as.factor(name)~., scales = "free")
dev.off()


# HGDP OCEANIA
ww_cumulative %>% filter(dataset == "HGDP", super_population == "OCEANIA") %>% summarize(max(ind_number))
# 56 haplotypes
# 28 individuals

label_data_ind <- data.frame(
  name = c("Neandertal", "Neandertal", "Neandertal", "Denisova", "Denisova", "Denisova"),
  dataset = c("LASIDAD", "deCODE", "HGDP", "LASIDAD", "deCODE", "HGDP"),
  label = c("LASI–DAD\n(n=2,679)", "deCODE\n(n=27,566)", "Oceania\n(n=28)",
            "LASI–DAD\n(n=2,679)", "deCODE\n(n=27,566)", "Oceania\n(n=28)"),
  x = c(4000, 45000, 1000, 4000, 45000, 800),       # adjust x position as needed
  y = c(1.45, 0.875, 0.7, 0.55, 0.2, 0.675)           # adjust y position as needed
)
label_data_ind$color <- cols[label_data_ind$dataset]


pdf("plots/HGDP_decode_LASI_NeaDen_all_num.pdf", height = 7, width = 5)
ww_cumulative %>% filter((dataset %in% c("deCODE", "LASIDAD")) | (dataset == "HGDP" & super_population == "OCEANIA") ) %>%
  mutate(name = factor(name, levels = c("Neandertal", "Denisova"), ordered = T)) %>%
  ggplot(aes(x = ind_number, y = value / 1000, color = dataset)) +
  geom_line() +
  theme_bw() +
  scale_color_manual(values = cols, guide = "none") +
  labs(y = "Total Sequence (Gb)", x = "Individuals") +
  geom_text(data = label_data_ind,
            aes(x = x, y = y, label = label, color = dataset),
            hjust = 0, vjust = 1, size = 3.5, show.legend = FALSE) +
  theme(text = element_text(size = 10)) +
  facet_grid(as.factor(name)~., scales = "free")
dev.off()



label_data_frac <- data.frame(
  name = c("Neandertal", "Neandertal", "Neandertal", "Denisova", "Denisova", "Denisova"),
  dataset = c("LASIDAD", "deCODE", "HGDP", "LASIDAD", "deCODE", "HGDP"),
  label = c("LASI–DAD\n(n=2,679)", "deCODE\n(n=27,566)", "Oceania\n(n=28)",
            "LASI–DAD\n(n=2,679)", "deCODE\n(n=27,566)", "Oceania\n(n=28)"),
  x = rep(0.85, 6),       # adjust x position as needed
  y = c(1.45, 0.9, 0.4, 0.675, 0.2, 0.535)
)
label_data_frac$color <- cols[label_data_ind$dataset]

pdf("plots/HGDP_decode_LASI_NeaDen_all_frac.pdf", height = 7, width = 5)
ww_cumulative %>% filter((dataset %in% c("deCODE", "LASIDAD")) | (dataset == "HGDP" & super_population == "OCEANIA") ) %>%
  mutate(name = factor(name, levels = c("Neandertal", "Denisova"), ordered = T)) %>%
  ggplot(aes(x = ind_frac, y = value / 1000, color = dataset)) +
 geom_line() +
  theme_bw() +
  scale_color_manual(values = cols, guide = "none") +
  labs(y = "Total Sequence (Gb)", x = "Individuals") +
  geom_text(data = label_data_frac,
            aes(x = x, y = y, label = label, color = dataset),
            hjust = 0, vjust = 1, size = 3.35, show.legend = FALSE) +
  theme(text = element_text(size = 14), strip.text = element_text(size = 14)) +
  facet_grid(as.factor(name)~., scales = "free")
dev.off()

pdf("plots/HGDP_decode_LASI_NeaDen_all_frac.pdf", height = 7, width = 5)
ww_cumulative %>% filter((dataset %in% c("deCODE", "LASIDAD")) | (dataset == "HGDP" & super_population == "OCEANIA") ) %>%
  mutate(name = factor(name, levels = c("Neandertal", "Denisova"), ordered = T)) %>%
  ggplot(aes(x = ind_frac, y = value / 1000, color = dataset)) +
  geom_line() +
  theme_bw() +
  scale_color_manual(values = cols, guide = "none") +
  labs(y = "Total Sequence (Gb)", x = "Individuals") +
  geom_text(data = label_data_frac,
            aes(x = x, y = y, label = label, color = dataset),
            hjust = 0, vjust = 1, size = 3, show.legend = FALSE) +
  facet_grid(as.factor(name)~., scales = "free")
dev.off()

# also make plots that include 1000G
# include GAv1 (from Mara)
GA_cumulative <- read.table("GenomeAsiaV1_cumulative_regions.tsv", header = T)
GA_cumulative <- pivot_longer(GA_cumulative, cols = c("Neandertal", "Denisova")) %>% 
  select(name, value, dataset, ind_number, super_population)

ww_cumulative <- rbind(select(ww_cumulative, -c(ind_frac)), GA_cumulative)

pdf("plots/LASI_choin_GAv1_1000G_NeaDen_all_frac.pdf", height = 6, width = 9)
ww_cumulative %>% filter((dataset %in% c("LASIDAD", "choin", "GenomeAsiav1")) | (dataset == "1000genomes" & super_population %in% c("CENTRAL_SOUTH_ASIA", "EAST_ASIA") ) | (dataset == "HGDP" & super_population == "OCEANIA")) %>%
  mutate(name = factor(name, levels = c("Neandertal", "Denisova"), ordered = T)) %>%
  ggplot(aes(x = ind_number, y = value / 1000, color = paste0(dataset, "-", super_population))) +
  geom_line() +
  theme_bw() +
  labs(y = "Total Sequence (Gb)", x = "Individuals") +
  theme(text = element_text(size = 12)) +
  facet_grid(as.factor(name)~., scales = "free")
dev.off()

##########################
### CHOIN OCEANIA ONLY ###
##########################

choin_cumulative <- read.table("choin_oceania_sequential_cumulative_Densiovan.tsv", header = T)

# read in laurit's results (so we can plot together)
ww_cumulative <- read.table("Final_out_from_Laurits.txt", sep = "\t", header = T)

ww_cumulative <- group_by(ww_cumulative, dataset, super_population) %>% mutate(ind_frac = ind_number / max(ind_number)) %>% ungroup()
# combine
ww_cumulative <- choin_cumulative %>% select(ind_number, ind_frac, Denisova, Neandertal, super_population, dataset) %>% rbind(select(ww_cumulative, ind_number, ind_frac, Denisova, Neandertal, super_population, dataset))
ww_cumulative <- ww_cumulative %>% pivot_longer(cols = c(Denisova, Neandertal))

cols <- c(
  "LASIDAD" = "#B8860B",     # dark goldenrod
  "deCODE" = "#CC6677",      # dark pink
  "choin" = "#117733",       # dark green
  "HGDP" = "#117733"         # same as choin
)

label_data_frac <- data.frame(
  name = c("Neandertal", "Neandertal", "Neandertal", "Denisova", "Denisova", "Denisova"),
  dataset = c("LASIDAD", "deCODE", "choin", "LASIDAD", "deCODE", "choin"),
  label = c("LASI–DAD\n(n=2,679)", "deCODE\n(n=27,566)", "Oceania\n(n=270)",
            "LASI–DAD\n(n=2,679)", "deCODE\n(n=27,566)", "Oceania\n(n=270)"),
  x = rep(0.85, 6),       # adjust x position as needed
  y = c(1.5, 1.2, 0.925, 0.75, 0.275, 1)          # adjust y position as needed
)
label_data_frac$color <- cols[label_data_ind$dataset]


pdf("plots/choinOceania_decode_LASI_NeaDen_all_frac.pdf", height = 7, width = 5)
ww_cumulative %>% filter(dataset %in% c("choin", "deCODE", "LASIDAD")) %>%
  mutate(name = factor(name, levels = c("Neandertal", "Denisova"), ordered = T)) %>%
  ggplot(aes(x = ind_frac, y = value / 1000, color = dataset)) +
  geom_line() +
  theme_bw() +
  scale_color_manual(values = cols, guide = "none") +
  labs(y = "Total Sequence (Gb)", x = "Proportion of Data Sequenced") +
  geom_text(data = label_data_frac,
            aes(x = x, y = y, label = label, color = dataset),
            hjust = 0, vjust = 1, size = 3, show.legend = FALSE) +
  facet_grid(as.factor(name)~., scales = "free")
dev.off()

label_data_ind <- data.frame(
  name = c("Neandertal", "Neandertal", "Neandertal", "Denisova", "Denisova", "Denisova"),
  dataset = c("LASIDAD", "deCODE", "choin", "LASIDAD", "deCODE", "choin"),
  label = c("LASI–DAD\n(n=2,679)", "deCODE\n(n=27,566)", "Oceania\n(n=270)",
            "LASI–DAD\n(n=2,679)", "deCODE\n(n=27,566)", "Oceania\n(n=270)"),
  x = c(4000, 45000, 1500, 4000, 45000, 1500),       # adjust x position as needed
  y = c(1.45, 1.2, 1, 0.75, 0.275, 1)           # adjust y position as needed
)
label_data_ind$color <- cols[label_data_ind$dataset]

#  scale_color_manual(values = c("orange2", "magenta4", "darkblue")) +

pdf("plots/choinOceania_decode_LASI_NeaDen_all_num.pdf", height = 7, width = 5)
ww_cumulative %>% filter(dataset %in% c("choin", "deCODE", "LASIDAD")) %>%
  mutate(name = factor(name, levels = c("Neandertal", "Denisova"), ordered = T)) %>%
  ggplot(aes(x = ind_number, y = value / 1000, color = dataset)) +
  geom_line() +
  theme_bw() +
  scale_color_manual(values = cols, guide = "none") +
  labs(y = "Total Sequence (Gb)", x = "# Sequenced Individuals") +
  geom_text(data = label_data_ind,
            aes(x = x, y = y, label = label, color = dataset),
            hjust = 0, vjust = 1, size = 3, show.legend = FALSE) +
  facet_grid(as.factor(name)~., scales = "free")
dev.off()





#############################
### GENERATING CUMULATIVE ###
#############################

segments_choin <- read.table("summary_files/overlapping_archaic_fragments_08.tsv", header = T)


# add on summary statistics
segments_choin <- segments_choin %>% mutate(ND_type = ifelse(pmax(Altai, Vindija, Chagyrskaya) == 0 & Denisova == 0, "None",
                                                             ifelse(pmax(Altai, Vindija, Chagyrskaya) == Denisova, "Both",
                                                                    ifelse(pmax(Altai, Vindija, Chagyrskaya) > Denisova, "Neandertal",                                                           ifelse(pmax(Altai, Vindija, Chagyrskaya) < Denisova, "Denisova", "problem")))))
# summarize:
segments_choin %>% group_by(ND_type) %>% summarize(count = n())

###
## Functions
###

## combine haplotypes within an individual:
combine_frags <- function(truth_dip){
  # first remove fragments that are identical across haplotypes
  truth_dip <- truth_dip %>% select(chrom, start, end) %>% arrange(chrom, start, end) %>% unique()
  
  # create an empty dataframe to hold the combined results
  new_truth_dip <- data.frame()
  
  row <- 1
  # merge overlapping fragments across the two haplotypes (chromosomes)
  start <- truth_dip$start[row]
  end <- truth_dip$end[row]
  
  # Look ahead to combine overlapping fragments
  while (row < nrow(truth_dip)) { # end_1000 
    # If there's overlap, extend the `end` to the max `end` in overlapping rows
    if (end >= truth_dip$start[row + 1]){
      end <- max(end, truth_dip$end[row + 1])
      # move on to next fragment
      row <- row + 1
    }else{
      # Add the non-overlapping fragment (either individual or merged) to new dataframe
      new_frag <- data.frame("start" = start, "end" = end)
      new_truth_dip <- rbind(new_truth_dip, new_frag)
      # move on to next fragment, reassign starts and ends to the next fragment
      row <- row + 1
      start <- truth_dip$start[row]
      end <- truth_dip$end[row]
    }
  }
  
  return(new_truth_dip)
}

# to test:
inds=c("UV1134", "B00I2CS")
archaic <- "Denisova"
dataset <- segments_choin

Archaic_segments <- function(dataset,inds,archaic){
  data_tmp <- filter(dataset,individual %in% inds & ND_type==archaic)
  amount <- 0
  # tmp_chrom <- 3
  for(tmp_chrom in 1:22){
    data_tmp_chrom <- filter(data_tmp, chrom == tmp_chrom)
    combo_tmp <- combine_frags(data_tmp_chrom)
    chrom_amount <- sum(combo_tmp$end - combo_tmp$start)
    amount <- amount + chrom_amount
  }
  
  return(amount)
}

Archaic_segments(segments_choin, inds, archaic)
####
### RUN
####

# get constants
choin_inds <- length(unique(segments_choin$individual))
results_Den <- rep(NA, choin_inds)
results_Nea <- rep(NA, choin_inds)
choin_inds <- select(segments_choin, individual) %>% unique() %>% pull() %>% sample(replace = F)

# run (simultaneously for both denisovan and neanderthal)
for (n in 1:length(unique(segments_choin$individual))){
  print(n)
  inds <- choin_inds[1:n]
  n_amount_Den <- Archaic_segments(segments_choin, inds, "Denisova")
  n_amount_Nea <- Archaic_segments(segments_choin, inds, "Neandertal")
  
  results_Den[n] <- n_amount_Den
  results_Nea[n] <- n_amount_Nea
  print(n_amount_Den)
  print(n_amount_Nea)
}

# tranform into dataframe
choin_cumulative <- data.frame(ind_number = 1:length(choin_inds), amount_DEN = results_Den, amount_NEA = results_Nea)
choin_cumulative$Denisova <- choin_cumulative$amount_DEN / 1000000
choin_cumulative$Neandertal <- choin_cumulative$amount_NEA / 1000000

choin_cumulative$dataset <- "choin"
choin_cumulative$super_population <- "OCEANIA"

write.table(choin_cumulative, "choin_all_sequential_cumulative_Densiovan.tsv", quote = F, row.names = F, sep = "\t")

# plots to check:
# plot all populations 
pdf("plots/choin_combined_NeaDen_all.pdf", height = 5, width = 7)
ggplot(ww_cumulative, aes(x = ind_number, y = value, color = super_population)) +
  geom_line() +
  facet_grid(name~dataset, scale = "free_x") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
dev.off()

########################
## CHOIN OCEANIA ONLY ##
########################

sample_info <- read.table("sample_info/included_samples.tsv", sep = "\t", header = T, skip = 1)
segments_choin <- merge(segments_choin, sample_info, by.x = "individual", by.y = "Sample.ID")
segments_choin %>% group_by(Region) %>% summarize(length(unique(individual)))
segments_choin_oceania <- filter(segments_choin, Region == "Oceania") 

choin_inds <- length(unique(segments_choin_oceania$individual))
results_Den <- rep(NA, choin_inds)
results_Nea <- rep(NA, choin_inds)
choin_inds <- select(segments_choin_oceania, individual) %>% unique() %>% pull() %>% sample(replace = F)

# run (simultaneously for both denisovan and neanderthal)
for (n in 1:length(unique(segments_choin_oceania$individual))){
  print(n)
  inds <- choin_inds[1:n]
  n_amount_Den <- Archaic_segments(segments_choin_oceania, inds, "Denisova")
  n_amount_Nea <- Archaic_segments(segments_choin_oceania, inds, "Neandertal")
  
  results_Den[n] <- n_amount_Den
  results_Nea[n] <- n_amount_Nea
  print(n_amount_Den)
  print(n_amount_Nea)
}

# tranform into dataframe
choin_cumulative <- data.frame(ind_number = 1:length(choin_inds), amount_DEN = results_Den, amount_NEA = results_Nea)
choin_cumulative$Denisova <- choin_cumulative$amount_DEN / 1000000
choin_cumulative$Neandertal <- choin_cumulative$amount_NEA / 1000000
choin_cumulative <- choin_cumulative %>% mutate(ind_frac = ind_number / max(ind_number))
choin_cumulative$super_population <- "Oceania"
choin_cumulative$dataset <- "choin"
write.table(choin_cumulative, "choin_oceania_sequential_cumulative_Densiovan.tsv", quote = F, row.names = F, sep = "\t")



# repeat on the different sub-populations of the oceanian dataset
sample_info <- read.table("sample_info/included_samples.tsv", sep = "\t", header = T, skip = 1)
segments_choin <- merge(segments_choin, sample_info, by.x = "individual", by.y = "Sample.ID")
segments_choin %>% group_by(Region) %>% summarize(length(unique(individual)))


choin_cumulative <- data.frame()
for(region in c("East Asia", "Southeast Asia", "Oceania")){
  print(region)
  dat <- filter(segments_choin, Region == region)
  n_indivs <- length(unique(dat$individual))
  
  print("N individuals:")
  print(n_indivs)
  
  results_Den <- rep(NA, n_indivs)
  
  for (n in 1:length(unique(dat$individual))){
    print(n)
    inds <- sample(unique(dat$individual), n, replace = F)
    n_amount_Den <- Archaic_segments(dat, inds, "Denisova")
    
    results_Den[n] <- n_amount_Den
    print(n_amount_Den)
  }
  
  region_cumulative <- data.frame(ind_number = 1:n_indivs, amount = results_Den)
  region_cumulative$Denisova <- region_cumulative$amount / 1000000
  region_cumulative$dataset <- "choin"
  region_cumulative$super_population <- region
  choin_cumulative <- rbind(choin_cumulative, region_cumulative)
}

choin_cumulative$ind_number <- choin_cumulative$ind_number * 2
choin_cumulative$super_population <- ifelse(choin_cumulative$super_population == "East Asia", "EAST_ASIA", 
                                            ifelse(choin_cumulative$super_population == "Southeast Asia", "SOUTHEAST_ASIA", 
                                                   ifelse(choin_cumulative$super_population == "Oceania", "OCEANIA", 'unassigned')))
write.table(choin_cumulative, "choin_all_cumulative_Densiovan_region.tsv", quote = F, row.names = F, sep = "\t")
choin_cumulative <- read.table("choin_all_cumulative_Densiovan_region.tsv", sep = "\t", header = T)
# read in laurit's results (so we can plot together)
ww_cumulative <- read.table("Final_out_from_Laurits.txt", sep = "\t", header = T)

# combine
ww_cumulative <- choin_cumulative %>% select(ind_number, Denisova, super_population, dataset) %>% rbind(select(ww_cumulative, ind_number, Denisova, super_population, dataset))

pdf("plots/cumulative_region_all.pdf", height = 6, width = 8)
ggplot(ww_cumulative, aes(x = ind_number, y = Denisova, color = super_population)) +
  geom_line() +
  facet_wrap(~dataset, scale = "free_x") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
dev.off()

pdf("plots/cumulative_oceania.pdf", height = 4, width = 6)
ww_cumulative %>% filter(super_population == "OCEANIA") %>%
  ggplot(aes(x = ind_number, y = Denisova, color = dataset)) +
  geom_line() +
  lims(y = c(0, 1000)) +
  labs(title = "OCEANIA", y = "Denisova MB reconstructed") +
  theme_bw()
dev.off()

pdf("plots/cumulative_eastAsia.pdf", height = 4, width = 6)
ww_cumulative %>% filter(super_population == "EAST_ASIA") %>%
  ggplot(aes(x = ind_number, y = Denisova, color = dataset)) +
  geom_line() +
  labs(title = "EAST ASIA") +
  theme_bw()
dev.off()

pdf("plots/cumulative_regions_choin.pdf", height = 4, width = 6)
ww_cumulative %>% filter(dataset == "choin") %>%
  ggplot(aes(x = ind_number, y = Denisova, color = super_population)) +
  geom_line() +
  lims(y = c(0, 1000)) +
  labs(title = "Choin Dataset", y = "Denisova MB reconstructed") +
  theme_bw()
dev.off()

# repeat on Genome Asia V1:
# first generate a file with all of the fragments (collate_frags_GAv1)
GAv1_meta <- read.table("/global/scratch/users/zhangyulin9806/github/ArchaicMutRate/hmmix_output/GA100K_raw/poplabel")
colnames(GAv1_meta) <- c("individual", "population", "super_population", "dataset")

# repeat on the different sub-populations of the GAv1 dataset
segments_GA <- read.table("GAv1_fragments_08.txt", header = T)
segments_GA <- merge(segments_GA, GAv1_meta, by.x = "name", by.y = "individual") %>% rename("individual" = name)
segments_GA %>% group_by(super_population) %>% summarize(length(unique(individual)))
segments_GA$chrom <- gsub("chr", "", segments_GA$chrom)

# add on summary statistics
segments_GA <- segments_GA %>% mutate(ND_type = ifelse(pmax(AltaiNeandertal, Vindija33.19, Chagyrskaya.Phalanx) == 0 & Denisova == 0, "None",
                                                             ifelse(pmax(AltaiNeandertal, Vindija33.19, Chagyrskaya.Phalanx) == Denisova, "Both",
                                                                    ifelse(pmax(AltaiNeandertal, Vindija33.19, Chagyrskaya.Phalanx) > Denisova, "Neandertal",                                                           
                                                                           ifelse(pmax(AltaiNeandertal, Vindija33.19, Chagyrskaya.Phalanx) < Denisova, "Denisova", "problem")))))


GA_cumulative <- data.frame()
for(region in c( "NortheastAsia", "SoutheastAsia")){ #, "SouthAsia"
  print(region)
  dat <- filter(segments_GA, super_population == !!region)
  n_indivs <- length(unique(dat$individual))
  
  print("N individuals:")
  print(n_indivs)
  
  results_Den <- rep(NA, n_indivs)
  results_Nea <- rep(NA, n_indivs)
  
  for (n in 1:length(unique(dat$individual))){
    print(n)
    inds <- sample(unique(dat$individual), n, replace = F)
    n_amount_Den <- Archaic_segments(segments_GA, inds, "Denisova")
    n_amount_Nea <- Archaic_segments(segments_GA, inds, "Neandertal")
  
    results_Den[n] <- n_amount_Den
    results_Nea[n] <- n_amount_Nea
    print(n_amount_Den)
    print(n_amount_Nea)
  }
  
  region_cumulative <- data.frame(ind_number = 1:n_indivs, amountDen = results_Den, amountNea = results_Nea)
  region_cumulative$Denisova <- region_cumulative$amountDen  / 1000000
  region_cumulative$Neandertal <- region_cumulative$amountNea  / 1000000
  region_cumulative$dataset <- "GenomeAsiav1"
  region_cumulative$super_population <- region
  GA_cumulative <- rbind(GA_cumulative, region_cumulative)
}

write.table(GA_cumulative, "GenomeAsiaV1_cumulative_regions.tsv", sep = "\t", quote = F)



