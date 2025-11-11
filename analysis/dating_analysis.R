rm(list = ls())
library(MASS)
library(tidyverse)


# plot how much synthetic we recover from each population:
synth <- read.table("results/1000G_synthetic_counts.txt", header = T, fill = NA)
synth <- pivot_wider(synth, names_from = "kind", values_from = length)
synth$assigned <- synth$Dhigh + synth$Dlow + synth$Nea
synth$other <- synth$full - synth$assigned
synth$assignedStrict <- synth$Dhigh08 + synth$Dlow38 + synth$Nea06
synth$otherStrict <- synth$full - synth$assignedStrict
synth <- synth %>% pivot_longer(cols = c(other, otherStrict, Dlow, Dhigh, Nea, Dlow38, Dhigh08, Nea06), names_to = "kind", values_to = "length")
synth$dataset <- "1000G"

# our classification between the two
ggplot(synth,aes(x = pop, y = length/1000000, fill = kind)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(y = "Mb recovered") +
  theme_bw()

asc <- read.table("results/1000G_ND_positions.txt", header = T)
asc <- group_by(asc, pop, kind) %>% mutate("ND" = sum(tot))
asc$dataset <- "1000G"

pos <- read.table("results/segregating_positions.txt", header = F)
colnames(pos) <- c("pop", "kind", "count")
pos$dataset <- "1000G"

ggplot(pos, aes(x = pop, y = count / 1000, fill = pop)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(y = "Number of positions\n(in thousands)") +
  guides(fill = "none") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  facet_wrap(~kind, scales = "free")

counts <- merge(asc, pos, by = c("pop", "kind"))
counts$other <- counts$count - counts$ND
counts <- counts %>% pivot_wider(names_from = asc, values_from = tot)
counts <- pivot_longer(counts, cols = c(other, ND01AfrAnc, ND10AfrAnc))
ggplot(filter(counts, name != "other"), aes(x = pop, y = value, fill = name)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(y = "ND positions") +
  scale_fill_manual(values = c("orange", "blue"), name = "ND type") +
  facet_wrap(~kind, scales = "free") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12))

# make a plot of just EAS
ggplot(filter(counts, name != "other", pop == "EAS"), aes(x = kind, y = value, fill = name)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(y = "ND positions") +
  scale_fill_manual(values = c("orange", "blue"), name = "ND type") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12))

dates <- read.table("results/dates.txt", header = T)
dates <- group_by(dates, population, kind, asc) %>% mutate(type = seq_along(kind)) %>%
  mutate(type = ifelse(type == 1, "raw", ifelse(type == 2, "corrected", "inferred_date")))
dates$include <- (grepl("D|full|modern", dates$kind) & grepl("ND01", dates$asc)) | (grepl("Nea|full|modern", dates$kind) & grepl("ND10", dates$asc))
dates$dataset <- "1000G"

ggplot(filter(dates, type == "raw", include), aes(x = paste0(kind), y = date, color = asc)) +
  geom_point() +
  geom_errorbar(aes(ymin = low, ymax = high)) +
  labs(y = "Uncorrected tAdmix\n(in generations)", color = "ascertainment") +
  facet_wrap(~population, scales = 'free_x') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

pdf(paste0("plots/dates/SAS_EAS_pops_dates.pdf"), height = 4, width = 6)
ggplot(filter(dates, type == "raw"), aes(x = paste0(kind), y = date, color = asc)) +
  geom_point() +
  geom_errorbar(aes(ymin = low, ymax = high)) +
  labs(y = "Uncorrected tAdmix\n(in generations)", color = "ascertainment") +
  facet_wrap(~population, scales = 'free_x') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()

pdf("plots/dates/SAS.pdf", height = 4, width = 6)
dates %>% filter(type == "raw", population %in% c("SAS"), kind != "Dhigh") %>%
  mutate(label = ifelse(kind == "full", "inferred archaic", kind)) %>% 
  mutate(label = factor(label, levels = c("modern", "inferred archaic", "Nea", "Dlow"), ordered = T)) %>%
  ggplot(aes(x = label, y = date, color = asc)) +
    geom_point() +
    geom_errorbar(aes(ymin = low, ymax = high)) +
    labs(y = "Uncorrected tAdmix\n(in generations)", color = "ascertainment") +
    facet_wrap(~population, scales = 'free_x') +
    theme_bw() +
    theme(axis.title.x = element_blank())
dev.off()

get_CI <- function(ND01, ND10){
  ratio <- sort(ND10/ND01)
  estimates <- sort(47000 - (ratio - 1)*47000, decreasing = F)
  return(c(mean(ratio), ratio[5], ratio[95], mean(estimates), estimates[5], estimates[95]))
}

EAS_ND01 <- sample(1138:1233, size = 100, replace = T)
EAS_ND10 <- sample(1246:1281, size = 100, replace = T)

EASnea_ND10 <- sample(1188:1223, size = 100, replace = T)


SASmodern_ND01 <- sample(1020:1144, size = 100, replace = T)
SASmodern_ND10 <- sample(1155:1192, size = 100, replace = T)

SAS_ND01 <- sample(997:1125, size = 100, replace = T)
SAS_ND10 <- sample(1112:1149, size = 100, replace = T)

SASlocal_ND01 <- sample(962:1060, size = 100, replace = T)
SASlocal_ND10 <- sample(1078:1116, size = 100, replace = T)
SAS <- get_CI(SAS_ND01, SAS_ND10)
SASlocal <- get_CI(SASlocal_ND01, SASlocal_ND10)
SASmodern <- get_CI(SASmodern_ND01, SASmodern_ND10)

SS_ND01 <- sample(992:1008, size = 100, replace = T)
SS_ND10 <- sample((1121+16):(1121-16), size = 100, replace = T)
SS <- get_CI(SS_ND01, SS_ND10)


SAS_results <- data.frame("method" = c("inferred archaic", "local", "modern", "SS Oceania"), "date" = c(SAS[4], SASlocal[4], SASmodern[4], SS[4]), low =c(SAS[5], SASlocal[5], SASmodern[5], SS[5]), high = c(SAS[6], SASlocal[6], SASmodern[6], SS[6]))
SAS_results$method <- factor(SAS_results$method, levels = c("modern", "inferred archaic", "local", "SS Oceania"))

pdf("plots/dates/SAS_dates.pdf", height = 4, width = 6)
ggplot(SAS_results, aes(x = method, y = date, color = method == "SS Oceania")) +
  geom_point() +
  geom_hline(yintercept = 47000, linetype = "dashed", color = "blue") +
  geom_errorbar(aes(ymin = low, ymax = high)) +
  labs(y = "Inferred Date of Denisovan Gene Flow") +
  scale_color_manual(values = c("black", "red"), guide = "none") +
  theme_bw() +
  lims(y= c(35000, 50000))
dev.off()



EASlow_ND01 <- sample(1293:1497, size = 100, replace = T)
SASlow_ND01 <- sample(962:1060, size = 100, replace = T)

EAShigh_ND01 <- sample(1092:1185, size = 100, replace = T)

SAS <- get_CI(SAS_ND01, SAS_ND10)
SASlocal <- get_CI(SASlocal_ND01, SASlocal_ND10)
EAS <- get_CI(EAS_ND01, EAS_ND10)
EAShigh <- get_CI(EAShigh_ND01, EASnea_ND10)
SASlow <- get_CI(SASlow_ND01, SAS_ND10)

get_CI_higher <- function(ND01, ND10){
  ratio <- sort(ND01/ND10)
  estimates <- sort(47000 + (ratio-1)*47000, decreasing = F)
  return(c(mean(ratio), ratio[5], ratio[95], mean(estimates), estimates[5], estimates[95]))
}

EASlow <- get_CI_higher(EASlow_ND01, EASnea_ND10)

SS <- get_CI(SS_ND01, SS_ND10)


###############
### LASIDAD ###
###############

# include LASIDAD
synth_LD <- read.table("results/LASIDAD_synthetic_counts.txt", header = T)
synth_LD <- pivot_wider(synth_LD, names_from = "kind", values_from = length)
synth_LD$assigned <- synth_LD$Dhigh + synth_LD$Dlow + synth_LD$Nea
synth_LD$other <- synth_LD$full - synth_LD$assigned
synth_LD$assignedStrict <- synth_LD$Dhigh08 + synth_LD$Dlow38 + synth_LD$Nea06
synth_LD$otherStrict <- synth_LD$full - synth_LD$assignedStrict
synth_LD <- synth_LD %>% pivot_longer(cols = c(other, otherStrict, Dlow, Dhigh, Nea, Dlow38, Dhigh08, Nea06), names_to = "kind", values_to = "length")
synth_LD$dataset <- "LASI-DAD"
synth <- rbind(synth, synth_LD)
synth$pop <- factor(synth$pop, levels = c("CEU", "CHB", "CHS", "KHV", "CDX", "GIH", "PJL", "ITU", "STU", "EAS", "SAS", "South", "Central", "West", "East", "North", "North-East", "singlePulse", "multiPulse", "noMultiPulse"), ordered = T)
synth$thresh <- ifelse(synth$kind %in% c("Dhigh08", "Dlow38", "Nea06", "otherStrict"), "Strict", "")
synth$NP <- ifelse(synth$pop %in% c("multiPulse", "noMultiPulse", "singlePulse"), "LASI-DAD\nCombined", 
                   ifelse(synth$pop %in% c("EAS", "SAS"), "1000G\nCombined", synth$dataset))
# our classification between the two

pdf("plots/synthetic_summary/all_mbrecovered.pdf", height = 4, width = 8)
ggplot(filter(synth, thresh != "Strict", synth$NP != F), aes(x = pop, y = length/1000000, fill = kind)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(y = "Mb recovered") +
  theme_bw() +
  labs(x = "Population") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  facet_grid(~NP, scales = "free_x", space = "free_x")
dev.off()

asc_LD <- read.table("results/LASIDAD_NDAfrAnc_positions.txt", header = T)
asc_LD <- group_by(asc_LD, pop, kind) %>% mutate("ND" = sum(tot))
asc_LD$dataset <- "LASI-DAD"
asc <- rbind(asc, asc_LD)

pos_LD <- read.table("results/LASIDAD_segregating_positions.txt", header = T)
colnames(pos_LD) <- c("pop", "kind", "count")
pos_LD$dataset <- "LASI-DAD"
pos <- rbind(pos, pos_LD)
pos$pop <- factor(pos$pop, levels = c("CEU", "CHB", "PJL", "EAS", "SAS", "South", "Central", "West", "East", "North", "North-East", "singlePulse", "multiPulse"), ordered = T)

pdf("plots/synthetic_summary/all_segregating.pdf", height = 6, width = 8)
ggplot(pos, aes(x = as.factor(pop), y = count / 1000, fill = dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(y = "Number of Segregating Sites\n(in thousands)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  facet_wrap(~kind, scales = "free")
dev.off()

counts <- merge(asc, pos, by = c("pop", "kind", "dataset"))
counts$other <- counts$count - counts$ND
counts <- counts %>% pivot_wider(names_from = asc, values_from = tot)
counts <- pivot_longer(counts, cols = c(other, ND01AfrAnc, ND10AfrAnc))
counts <- filter(counts, (kind %in% c("Dhigh", "Dlow") & name == "ND01AfrAnc") | (kind %in% c("Nea") & name == "ND10AfrAnc") | kind %in% c("full"))
counts$pop <- factor(counts$pop, levels = c("CEU", "CHB", "PJL", "EAS", "SAS", "South", "Central", "West", "East", "North", "North-East", "singlePulse", "multiPulse"), ordered = T)

pdf("plots/synthetic_summary/all_classCounts.pdf", height = 6, width = 8)
ggplot(filter(counts, name != "other"), aes(x = as.factor(pop), y = value, fill = name)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(y = "ND positions") +
  scale_fill_manual(values = c("orange", "blue"), name = "ND type") +
  facet_wrap(~kind, scales = "free") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12))
dev.off()

dates_LD <- read.table("results/LASI_dates.txt", header = T)
dates_LD <- group_by(dates_LD, population, kind, asc) %>% mutate(type = seq_along(kind)) %>%
  mutate(type = ifelse(type == 1, "raw", ifelse(type == 2, "corrected", "inferred_date")))
dates_LD$dataset <- "LASI-DAD"

dates <- rbind(dates %>% filter(include) %>% select(-include), dates_LD)
write.table(dates, "results/LASI_1000G_dates.txt", quote = F, sep = "\t")

pdf("plots/dates/SAS_EAS_MP_SP.pdf", height = 6, width = 8)
ggplot(filter(dates, type == "raw", population %in% c("SAS", "EAS", "singlePulse", "multiPulse")), aes(x = paste0(kind), y = date, color = asc)) +
  geom_point() +
  geom_errorbar(aes(ymin = low, ymax = high)) +
  labs(y = "Uncorrected tAdmix\n(in generations)", color = "ascertainment") +
  facet_wrap(~population, scales = 'free_x') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()


pdf("plots/dates/all_pops_dates.pdf", height = 8, width = 10)
ggplot(filter(dates_LD, type == "raw"), aes(x = paste0(kind), y = date, color = asc)) +
  geom_point() +
  geom_errorbar(aes(ymin = low, ymax = high)) +
  labs(y = "Uncorrected tAdmix\n(in generations)", color = "ascertainment") +
  facet_wrap(~population, scales = 'free_x') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()

# get some LASI dates
SP_ND01 <- sample(1151:1254, size = 100, replace = T)
SP_ND10 <- sample(1153:1193, size = 100, replace = T)
SP <- get_CI_higher(SP_ND01, SP_ND10)
SP
SP <- get_CI(SP_ND01, SP_ND10)
SP

SPlocal_ND01 <- sample(834:931, size = 100, replace = T)
SPlocal_ND10 <- sample(1110:1145, size = 100, replace = T)
SP <- get_CI(SPlocal_ND01, SPlocal_ND10)

EASnea_ND10 <- sample(1188:1223, size = 100, replace = T)

SAS_ND01 <- sample(997:1125, size = 100, replace = T)
SAS_ND10 <- sample(1078:1116, size = 100, replace = T)

SS_ND01 <- sample(992:1008, size = 100, replace = T)
SS_ND10 <- sample((1121+16):(1121-16), size = 100, replace = T)

EASlow_ND01 <- sample(1293:1497, size = 100, replace = T)
SASlow_ND01 <- sample(962:1060, size = 100, replace = T)

EAShigh_ND01 <- sample(1092:1185, size = 100, replace = T)

SAS <- get_CI(SAS_ND01, SAS_ND10)
EAS <- get_CI(EAS_ND01, EAS_ND10)
EAShigh <- get_CI(EAShigh_ND01, EASnea_ND10)
SASlow <- get_CI(SASlow_ND01, SAS_ND10)

