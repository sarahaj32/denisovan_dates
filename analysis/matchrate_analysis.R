library(MASS)
library(tidyverse)

theme_set(theme_bw() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.text = element_text(size = 10)))

# countour plot parameters:
zoom=F
xlab="Match to Altai Neanderthal"
ylab="Match to Altai Denisovan"
level1=seq(0.3,0.9,0.1)
level2=seq(1,30,1)

# our hmmix results for CHB, with synthetic, non-overlapping genome

dat <- read.table("data/1000G/CHB/fragments/CHB_fragments_snps10_08_nonoverlapping_full.bed", header = T)
# for now, lets divide by the number of snps
dat$D_match <- dat$Denisova / dat$snps
dat$N_match <- dat$Vindija / dat$snps
dat <- filter(dat, snps > 10)
# how many total Gb do we recover?
sum(dat$length) / 1000000000
# 0.736

# compare dividing by snps vs accessible/comparable positions
# ok this is good that it mostly increases
ggplot(dat, aes(x = D_match, y = den_affinity)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red")

ggplot(dat, aes(x = den_overlap / snps, y = den_affinity)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red")
# this is mostly due to the denominator decreasing (by design)

ggplot(dat, aes(x = snps, y = comp)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red")


# these also increase. A larg part of this is becuase I now match to any of the 3 neanderthals
ggplot(dat, aes(x = nea_overlap / snps, y = nea_affinity)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red")


mean(dat$nea_affinity[dat$nea_affinity > 0])
mean(dat$den_affinity[dat$den_affinity > 0])

popname = "CHB"

# histogram
pdf("plots/matchrate/CHB_synthetic_snps_histogram.pdf", height = 4, width = 6)
ggplot(filter(dat,den_affinity > 0.3, nea_affinity < 0.3), aes(x = den_affinity)) +
  geom_histogram() +
  labs(title = popname) +
  theme_bw()
dev.off()

filter(dat,den_affinity > 0.3, nea_affinity < 0.3) %>% nrow()
# 398
filter(dat,den_affinity > 0.3, nea_affinity < 0.3) %>% summarize(sum(length))

# countour plots
pdf("plots/matchrate/CHB_synthetic_snps_countour.pdf", height = 4, width = 4.5)
contour(kde2d(dat$N_match,dat$D_match, lims=c(0,1,0,1),n=100),xaxs="i",yaxs="i",xlab=xlab,ylab=ylab,main=popname,levels=level1,las=1,cex.lab=1.5,cex.axis=1.4,cex.main=1.5,lty=5,labcex=0.8);abline(h=seq(0.2,0.8,0.2),v=seq(0.2,0.8,0.2),lty=3,col="gray")
contour(kde2d(dat$N_match,dat$D_match, lims=c(0,1,0,1),n=100),xaxs="i",yaxs="i",xlab=xlab,ylab=ylab,main=popname,levels=level2,las=1,cex.lab=1.5,cex.axis=1.4,cex.main=1.5,lty=5,labcex=0.8,add=T);abline(h=seq(0.2,0.8,0.2),v=seq(0.2,0.8,0.2),lty=3,col="gray")
dev.off()
# ok, i see the same pulses but the scale seems to be wrong. I think it's because I'm dividing by total derived snps,
# not those that can be compared to archaics

# MY SUMMARY STATS
pops=c("CEU","TSI","CDX","CHB","CHS","JPT","KHV","BEB","GIH","ITU","PJL","STU")
popnames=c("NW European (Utah)", "Toscani (Italy)","Chinese Dai","Han Chinese (Beijing)","Southern Han Chinese","Japanese (Tokyo)","Kinh (Vietnam)","Bengali (Bangladesh)","Gujarati Indian (Houston)","Indian Telugu (UK)","Punjabi (Pakistan)","Sri Lankan Tamil (UK)")
#i <- 4
pops=c("CEU","CHB","PJL","EAS","SAS")
popnames=c("NW European (Utah)", "Han Chinese (Beijing)","Punjabi (Pakistan)","EAS","SAS")

for(i in 1:length(pops)){
  pop = pops[i]
  popname = popnames[i]
  print(pop)
  #dat <- read.table(paste0("data/1000G_", pop, "/", pop, "_summary_snps10_08_nonoverlapping.txt"), header = T)
  dat <- read.table(paste0("data/1000G/", pop, "/fragments/", pop, "_fragments_snps10_08_nonoverlapping_full.bed"), header = T)
  dat <- filter(dat, comp > 0)
  dat <- filter(dat, comp > 5)
  dat <- filter(dat, comp >= 10)
  dat$nea_SJ <- dat$ND10 + dat$ND11
  dat$den_SJ <- dat$ND01 + dat$ND11
  dat$den_affinity <- dat$den_SJ / dat$comp #(dat$ND00 + dat$ND01 + dat$ND10)
  dat$nea_affinity <- dat$nea_SJ / dat$comp #(dat$ND00 + dat$ND01 + dat$ND10)
  
  # countour plots
  pdf(paste0("plots/matchrate/", pop, "_synthetic_comparable_contour.pdf"), height = 4, width = 4.5)
  contour(kde2d(dat$nea_affinity,dat$den_affinity, lims=c(0,1,0,1),n=100),xaxs="i",yaxs="i",xlab=xlab,ylab=ylab,main=popname,levels=level1,las=1,cex.lab=1.5,cex.axis=1.4,cex.main=1.5,lty=5,labcex=0.8);abline(h=seq(0.2,0.8,0.2),v=seq(0.2,0.8,0.2),lty=3,col="gray")
  contour(kde2d(dat$nea_affinity,dat$den_affinity, lims=c(0,1,0,1),n=100),xaxs="i",yaxs="i",xlab=xlab,ylab=ylab,main=popname,levels=level2,las=1,cex.lab=1.5,cex.axis=1.4,cex.main=1.5,lty=5,labcex=0.8,add=T);abline(h=seq(0.2,0.8,0.2),v=seq(0.2,0.8,0.2),lty=3,col="gray")
  dev.off()
  
  
  pdf(paste0("plots/matchrate/", pop, "_synthetic_comparable_histogram.pdf"), height = 4, width = 6)
  p <- ggplot(filter(dat, den_affinity > 0.3, nea_affinity < 0.3), aes(x = den_affinity)) +
    geom_histogram() +
    labs(title = popname) +
    theme_bw()
  print(p)
  dev.off()
  
  pdf(paste0("plots/matchrate/", pop, "_synthetic_comparable_histogram_ND11.pdf"), height = 4, width = 6)
  p <- ggplot(filter(dat, ND01>1, ND10 == 0), aes(x = den_affinity)) +
    geom_histogram() +
    labs(title = popname) +
    theme_bw()
  print(p)
  dev.off()
}

# get stats from some of the populations
dat <- read.table("data/1000G/CHB/fragments/CHB_fragments_snps10_08_nonoverlapping_full.bed", header = T)
dat <- filter(dat, comp >= 10)
sum(dat$length) / 1000000
# great they are the same! I can just use my previous calculations then

# how many total Gb do we recover?
sum(dat$length) / 1000000000
# 0.755 (0.659 after filtering to at least 10 comparable frags)
nrow(dat)
#6712

dat %>% dplyr::select(comp, nocomp, ND01, ND10, ND11, ND00) %>%
  pivot_longer(cols = c(comp, nocomp, ND01, ND10, ND11, ND00)) %>%
  group_by(name) %>%
  summarize(sum(value))

mean(dat$nea_affinity)
mean(dat$nea_affinity[dat$nea_affinity > 0])
mean(dat$den_affinity[dat$den_affinity > 0])

popname = "Han Chinese (Beijing)"
# countour plots
pdf("plots/matchrate/CHB_synthetic_comparable_contour.pdf", height = 4, width = 6)
contour(kde2d(dat$nea_affinity,dat$den_affinity, lims=c(0,1,0,1),n=100),xaxs="i",yaxs="i",xlab=xlab,ylab=ylab,main=popname,levels=level1,las=1,cex.lab=1.5,cex.axis=1.4,cex.main=1.5,lty=5,labcex=0.8);abline(h=seq(0.2,0.8,0.2),v=seq(0.2,0.8,0.2),lty=3,col="gray")
contour(kde2d(dat$nea_affinity,dat$den_affinity, lims=c(0,1,0,1),n=100),xaxs="i",yaxs="i",xlab=xlab,ylab=ylab,main=popname,levels=level2,las=1,cex.lab=1.5,cex.axis=1.4,cex.main=1.5,lty=5,labcex=0.8,add=T);abline(h=seq(0.2,0.8,0.2),v=seq(0.2,0.8,0.2),lty=3,col="gray")
dev.off()


pdf("plots/matchrate/CHB_synthetic_all_histogram.pdf", height = 4, width = 6)
ggplot(filter(dat), aes(x = den_affinity)) +
  geom_histogram() +
  labs(title = popname) +
  theme_bw()
dev.off()

pdf("plots/matchrate/CHB_synthetic_comparable_histogram.pdf", height = 4, width = 6)
ggplot(filter(dat, den_affinity > 0.3, nea_affinity < 0.2), aes(x = den_affinity)) +
  geom_histogram() +
  labs(title = popname) +
  theme_bw()
dev.off()
filter(dat, den_affinity > 0.3, nea_affinity < 0.3) %>% nrow()

ggplot(filter(dat, den_affinity > 0.3, nea_affinity < 0.2), aes(x = den_affinity)) +
  geom_histogram() +
  labs(title = popname) +
  theme_bw()

pdf("plots/matchrate/CHB_synthetic_comparable_histogram_ND11.pdf", height = 4, width = 6)
ggplot(filter(dat, ND01>1, ND10 == 0), aes(x = den_affinity)) +
  geom_histogram() +
  labs(title = popname) +
  theme_bw()
dev.off()

pdf("plots/matchrate/CHB_synthetic_comparable_histogram_LSscale.pdf", height = 4, width = 6)
ggplot(filter(dat,ND01>1, ND10 == 0, comp >= 30), aes(x = den_affinity)) +
  geom_histogram(bins = 40) +
  lims(y = c(0, 35)) +
  labs(title = popname) +
  theme_bw()
dev.off()

# get summary stats for this
filter(dat, ND01>1, ND10 == 0) %>% nrow()
filter(dat, ND01>1, ND10 == 0) %>% summarize(sum(length) / 1000000) 

filter(dat, den_affinity >= 0.2, nea_affinity <= 0.2) %>% nrow()
filter(dat, den_affinity >= 0.2, nea_affinity <= 0.2) %>% summarize(sum(length) / 1000000) 

filter(dat, ND01>1, ND10 == 0, den_affinity > 0.3, den_affinity < 0.75) %>% nrow()
filter(dat, ND01>1, ND10 == 0, den_affinity > 0.3, den_affinity < 0.75) %>% summarize(sum(length) / 1000000) 

filter(dat, ND01>1, ND10 == 0, den_affinity >= 0.75) %>% nrow()
filter(dat, ND01>1, ND10 == 0, den_affinity >= 0.75) %>% summarize(sum(length) / 1000000) 


# also look at potential cutoffs
# number of positions
dat %>% filter(den_affinity > 0.7) %>% select(ND01) %>% sum()
# length
dat %>% filter(den_affinity > 0.7) %>% select(length) %>% sum() / 1000000

dat %>% filter(den_affinity < 0.7 & den_affinity > 0.4) %>% select(ND01) %>% sum()
dat %>% filter(den_affinity < 0.7 & den_affinity > 0.4) %>% select(length) %>% sum() / 1000000



# REPEAT GIH
dat <- read.table("data/1000G_GIH/GIH_summary_snps10_08_nonoverlapping.txt", header = T)
dat <- filter(dat, comp > 0)
dat <- filter(dat, comp >= 10)
dat$nea_SJ <- dat$ND10 + dat$ND11
dat$den_SJ <- dat$ND01 + dat$ND11
dat$den_affinity <- dat$den_SJ / dat$comp
dat$nea_affinity <- dat$nea_SJ / dat$comp

# how many total Gb do we recover?
sum(dat$length) / 1000000000
# 0.83 after filtering to at least 10 comparable frags
nrow(dat)
#6242

dat %>% dplyr::select(comp, nocomp, ND01, ND10, ND11, ND00) %>%
  pivot_longer(cols = c(comp, nocomp, ND01, ND10, ND11, ND00)) %>%
  group_by(name) %>%
  summarize(sum(value))

mean(dat$nea_affinity)
mean(dat$nea_affinity[dat$nea_affinity > 0])
mean(dat$den_affinity[dat$den_affinity > 0])

popname = "GIH"
# countour plots
pdf(paste0("plots/matchrate/", popname, "_synthetic_comparable_contour.pdf"), height = 4, width = 4.5)
contour(kde2d(dat$nea_affinity,dat$den_affinity, lims=c(0,1,0,1),n=100),xaxs="i",yaxs="i",xlab=xlab,ylab=ylab,main=popname,levels=level1,las=1,cex.lab=1.5,cex.axis=1.4,cex.main=1.5,lty=5,labcex=0.8);abline(h=seq(0.2,0.8,0.2),v=seq(0.2,0.8,0.2),lty=3,col="gray")
contour(kde2d(dat$nea_affinity,dat$den_affinity, lims=c(0,1,0,1),n=100),xaxs="i",yaxs="i",xlab=xlab,ylab=ylab,main=popname,levels=level2,las=1,cex.lab=1.5,cex.axis=1.4,cex.main=1.5,lty=5,labcex=0.8,add=T);abline(h=seq(0.2,0.8,0.2),v=seq(0.2,0.8,0.2),lty=3,col="gray")
dev.off()


pdf(paste0("plots/matchrate/", popname, "_synthetic_comparable_histogram.pdf"), height = 4, width = 6)
ggplot(filter(dat, den_affinity > 0.3, nea_affinity < 0.3), aes(x = den_affinity)) +
  geom_histogram() +
  labs(title = popname) +
  theme_bw()
dev.off()

pdf(paste0("plots/matchrate/", popname, "_synthetic_comparable_histogram_ND11.pdf"), height = 4, width = 6)
ggplot(filter(dat, ND01>1, ND10 == 0), aes(x = den_affinity)) +
  geom_histogram() +
  labs(title = popname) +
  theme_bw()
dev.off()

pdf(paste0("plots/matchrate/", popname, "_synthetic_comparable_histogram_LSscale.pdf"), height = 4, width = 6)
ggplot(filter(dat,ND01>1, ND10 == 0, comp >= 30), aes(x = den_affinity)) +
  geom_histogram(bins = 40) +
  lims(y = c(0, 35)) +
  labs(title = popname) +
  theme_bw()
dev.off()

# get summary stats for this
filter(dat, ND01>1, ND10 == 0) %>% nrow()
filter(dat, ND01>1, ND10 == 0) %>% summarize(sum(length) / 1000000000) 


dat %>% filter(den_affinity > 0.7) %>% select(ND01) %>% sum()
dat %>% filter(den_affinity < 0.7 & den_affinity > 0.2) %>% select(ND01) %>% sum()





# look at Oceania (choin)
# focus on EAS and SAS 
ocean_dat <- read.table("data/choin/affinity_plots/Oceania_fragments_snps10_08_nonoverlapping_full.bed", header = T)
ocean_dat <- filter(ocean_dat, comp > 10)
sum(ocean_dat$length) / 1000000000
# 0.83 after filtering to at least 10 comparable frags
nrow(ocean_dat)

mean(ocean_dat$nea_affinity[ocean_dat$nea_affinity > 0])
mean(ocean_dat$den_affinity[ocean_dat$den_affinity > 0])


popname <- "Oceania"
xlab="Match to Altai Neanderthal"
ylab="Match to Altai Denisovan"

# countour plots
level1=seq(0.3,0.9,0.1)
level2=seq(1,30,1)
ocean_dat <- filter(ocean_dat, comp > 40)

pdf(paste0("plots/matchrate/oceania_synthetic_comparable_contour.pdf"), height = 4, width = 4.5)
contour(kde2d(ocean_dat$nea_affinity,ocean_dat$den_affinity, lims=c(0,1,0,1),n=100),xaxs="i",yaxs="i",xlab=xlab,ylab=ylab,main=popname,levels=level1,las=1,cex.lab=1.5,cex.axis=1.4,cex.main=1.5,lty=5,labcex=0.8);abline(h=seq(0.2,0.8,0.2),v=seq(0.2,0.8,0.2),lty=3,col="gray")
contour(kde2d(ocean_dat$nea_affinity,ocean_dat$den_affinity, lims=c(0,1,0,1),n=100),xaxs="i",yaxs="i",xlab=xlab,ylab=ylab,main=popname,levels=level2,las=1,cex.lab=1.5,cex.axis=1.4,cex.main=1.5,lty=5,labcex=0.8, add=T);abline(h=seq(0.2,0.8,0.2),v=seq(0.2,0.8,0.2),lty=3,col="gray")
dev.off()

pdf(paste0("plots/matchrate/oceania_synthetic_comparable_histogram.pdf"), height = 4, width = 4.5)
ggplot(filter(ocean_dat, den_affinity > 0.3, nea_affinity < 0.3), aes(x = den_affinity)) +
  geom_histogram() +
  labs(title = popname) +
  theme_bw()
dev.off()

pdf(paste0("plots/matchrate/oceania_synthetic_ND11_histogram.pdf"), height = 4, width = 4.5)
ggplot(filter(ocean_dat, ND01>1, ND10 == 0), aes(x = den_affinity)) +
  geom_histogram() +
  labs(title = popname) +
  theme_bw()
dev.off()

pdf(paste0("plots/matchrate/oceania_synthetic_highaffinity_histogram.pdf"), height = 4, width = 4.5)
ggplot(filter(ocean_dat, den_affinity > 0.7, nea_affinity < 0.3), aes(x = den_affinity)) +
  geom_histogram() +
  labs(title = popname) +
  theme_bw()
dev.off()

ggplot(filter(ocean_dat, den_affinity > 0.3, nea_affinity < 0.3, length > 100000), aes(x = den_affinity)) +
  geom_histogram() +
  labs(title = popname) +
  theme_bw()

pdf(paste0("plots/matchrate/oceania_synthetic_longND11_histogram.pdf"), height = 4, width = 4.5)
ggplot(filter(ocean_dat, ND01>1, ND10 == 0, length > 200000), aes(x = den_affinity)) +
  geom_histogram() +
  labs(title = popname) +
  theme_bw()
dev.off()

pdf(paste0("plots/matchrate/oceania_synthetic_long_histogram.pdf"), height = 4, width = 4.5)
ggplot(filter(ocean_dat, den_affinity > 0.3, nea_affinity < 0.3, length > 200000), aes(x = den_affinity)) +
  geom_histogram() +
  labs(title = popname) +
  theme_bw()
dev.off()

pdf(paste0("plots/matchrate/oceania_synthetic_snpsND11_histogram.pdf"), height = 4, width = 4.5)
ggplot(filter(ocean_dat, ND01>1, ND10 == 0, snps > 100), aes(x = den_affinity)) +
  geom_histogram() +
  labs(title = popname) +
  theme_bw()
dev.off()

pdf(paste0("plots/matchrate/oceania_allDen_length_histogram.pdf"), height = 4, width = 4.5)
ggplot(filter(ocean_dat, ND01>1, ND10 == 0, den_affinity == 1), aes(x = length / 1000)) +
  geom_histogram() +
  labs(title = popname) +
  theme_bw()
dev.off()

pdf(paste0("plots/matchrate/oceania_synthetic_snps_histogram.pdf"), height = 4, width = 4.5)
ggplot(filter(ocean_dat, den_affinity > 0.3, nea_affinity < 0.3, snps > 100), aes(x = den_affinity)) +
  geom_histogram() +
  labs(title = popname) +
  theme_bw()
dev.off()


# look at choin SEA and EA
SEA_dat <- read.table("data/choin/affinity_plots/SoutheastAsia_fragments_snps10_08_nonoverlapping_full.bed", header = T)
popname <- "Philippines (Agta and Cebuano)"
sum(SEA_dat$length) / 1000000000
nrow(SEA_dat)
mean(SEA_dat$nea_affinity[SEA_dat$nea_affinity > 0])
mean(SEA_dat$den_affinity[SEA_dat$den_affinity > 0])
SEA_dat <- filter(SEA_dat, comp > 40)



pdf(paste0("plots/matchrate/choinSEA_synthetic_comparable_contour.pdf"), height = 4, width = 4.5)
contour(kde2d(SEA_dat$nea_affinity,SEA_dat$den_affinity, lims=c(0,1,0,1),n=100),xaxs="i",yaxs="i",xlab=xlab,ylab=ylab,main=popname,levels=level1,las=1,cex.lab=1.5,cex.axis=1.4,cex.main=1.5,lty=5,labcex=0.8);abline(h=seq(0.2,0.8,0.2),v=seq(0.2,0.8,0.2),lty=3,col="gray")
contour(kde2d(SEA_dat$nea_affinity,SEA_dat$den_affinity, lims=c(0,1,0,1),n=100),xaxs="i",yaxs="i",xlab=xlab,ylab=ylab,main=popname,levels=level2,las=1,cex.lab=1.5,cex.axis=1.4,cex.main=1.5,lty=5,labcex=0.8, add=T);abline(h=seq(0.2,0.8,0.2),v=seq(0.2,0.8,0.2),lty=3,col="gray")
dev.off()

pdf(paste0("plots/matchrate/choinSEA_synthetic_comparable_histogram.pdf"), height = 4, width = 4.5)
ggplot(filter(SEA_dat, den_affinity > 0.3, nea_affinity < 0.3), aes(x = den_affinity)) +
  geom_histogram() +
  labs(title = popname) +
  theme_bw()
dev.off()

pdf(paste0("plots/matchrate/choinSEA_synthetic_ND11_histogram.pdf"), height = 4, width = 4.5)
ggplot(filter(SEA_dat, ND01>1, ND10 == 0), aes(x = den_affinity)) +
  geom_histogram() +
  labs(title = popname) +
  theme_bw()
dev.off()

# countour plots
level1=seq(0.1,0.9,0.1)
level2=seq(1,30,1)

# look at choin SEA and EA
EAS_dat <- read.table("data/choin/affinity_plots/EastAsia_fragments_snps10_08_nonoverlapping_full.bed", header = T)
popname <- "Taiwan (Paiwan and Atayal)"
sum(EAS_dat$length) / 1000000000
nrow(EAS_dat)
mean(EAS_dat$nea_affinity[EAS_dat$nea_affinity > 0])
mean(EAS_dat$den_affinity[EAS_dat$den_affinity > 0])
EAS_dat <- filter(EAS_dat, comp > 40)

pdf(paste0("plots/matchrate/choinEAS_synthetic_comparable_contour.pdf"), height = 4, width = 4.5)
contour(kde2d(EAS_dat$nea_affinity,EAS_dat$den_affinity, lims=c(0,1,0,1),n=100),xaxs="i",yaxs="i",xlab=xlab,ylab=ylab,main=popname,levels=level1,las=1,cex.lab=1.5,cex.axis=1.4,cex.main=1.5,lty=5,labcex=0.8);abline(h=seq(0.2,0.8,0.2),v=seq(0.2,0.8,0.2),lty=3,col="gray")
contour(kde2d(EAS_dat$nea_affinity,EAS_dat$den_affinity, lims=c(0,1,0,1),n=100),xaxs="i",yaxs="i",xlab=xlab,ylab=ylab,main=popname,levels=level2,las=1,cex.lab=1.5,cex.axis=1.4,cex.main=1.5,lty=5,labcex=0.8, add=T);abline(h=seq(0.2,0.8,0.2),v=seq(0.2,0.8,0.2),lty=3,col="gray")
dev.off()


pdf(paste0("plots/matchrate/choinEAS_synthetic_comparable_histogram.pdf"), height = 4, width = 4.5)
ggplot(filter(EAS_dat, den_affinity > 0.3, nea_affinity < 0.3), aes(x = den_affinity)) +
  geom_histogram() +
  labs(title = popname) +
  theme_bw()
dev.off()

pdf(paste0("plots/matchrate/choinEAS_synthetic_ND11_histogram.pdf"), height = 4, width = 4.5)
ggplot(filter(EAS_dat, ND01>1, ND10 == 0), aes(x = den_affinity)) +
  geom_histogram() +
  labs(title = popname) +
  theme_bw()
dev.off()









# focus on EAS and SAS 
EAS_dat <- read.table("data/1000G/EAS/fragments/EAS_fragments_snps10_08_nonoverlapping_full.bed", header = T)
popname <- "EAS"

# countour plots
level1=seq(0.3,0.9,0.1)
level2=seq(1,30,1)

pdf(paste0("plots/matchrate/EAS_synthetic_comparable_contour.pdf"), height = 4, width = 4.5)
contour(kde2d(EAS_dat$nea_affinity,EAS_dat$den_affinity, lims=c(0,1,0,1),n=100),xaxs="i",yaxs="i",xlab=xlab,ylab=ylab,main=popname,levels=level1,las=1,cex.lab=1.5,cex.axis=1.4,cex.main=1.5,lty=5,labcex=0.8);abline(h=seq(0.2,0.8,0.2),v=seq(0.2,0.8,0.2),lty=3,col="gray")
contour(kde2d(EAS_dat$nea_affinity,EAS_dat$den_affinity, lims=c(0,1,0,1),n=100),xaxs="i",yaxs="i",xlab=xlab,ylab=ylab,main=popname,levels=level2,las=1,cex.lab=1.5,cex.axis=1.4,cex.main=1.5,lty=5,labcex=0.8, add=T);abline(h=seq(0.2,0.8,0.2),v=seq(0.2,0.8,0.2),lty=3,col="gray")
dev.off()

# good to know that I can make these same plots using ggplot
#ggplot(EAS_dat, aes(x = nea_affinity, y = den_affinity)) +
#  #geom_density_2d()
#  stat_density_2d(geom = "contour", breaks = level1, linetype = "dashed", color = "black") +
#  stat_density_2d(geom = "contour", breaks = level2, linetype = "dashed", color = "black")  +
#  theme_minimal()

pdf(paste0("plots/matchrate/EAS_synthetic_comparable_histogram.pdf"), height = 4, width = 6)
EAS_dat %>% filter(den_affinity > 0.3, nea_affinity < 0.3) %>% 
  ggplot(aes(x = den_affinity)) +
  geom_histogram() +
  labs(title = "EAS", x = "Denisovan Affinity", y = "Number of Fragments")
dev.off()

nrow(EAS_dat)
sum(EAS_dat$length) / 1000000000
EAS_dat %>% filter(den_affinity > 0) %>% summarize(mean(den_affinity))
EAS_dat %>% filter(nea_affinity > 0) %>% summarize(mean(nea_affinity))

pdf(paste0("plots/matchrate/EAS_synthetic_comparable_allfrags_histogram.pdf"), height = 4, width = 6)
EAS_dat %>% 
  ggplot(aes(x = den_affinity)) +
  geom_histogram() +
  labs(title = "EAS") +
  theme_bw()
dev.off()

pdf(paste0("plots/matchrate/EAS_synthetic_comparable_histogram_ND11.pdf"), height = 4, width = 6)
ggplot(filter(EAS_dat, ND01>1, ND10 == 0), aes(x = den_affinity)) +
  geom_histogram() +
  labs(title = "EAS") +
  theme_bw()
dev.off()

# these are the initial groupings I was using:
EAS_dat$group <- ifelse(EAS_dat$nea_affinity <= 0.3 & EAS_dat$den_affinity >= 0.7, "Dhigh", 
                        ifelse(EAS_dat$nea_affinity <= 0.3 & EAS_dat$den_affinity >= 0.4 & EAS_dat$den_affinity < 0.7, "Dlow",
                               ifelse(EAS_dat$nea_affinity >= 0.5 & EAS_dat$den_affinity <= 0.3, "Nea", "other")))

# great, I confirmed that this matches what I find in my pipeline
pdf(paste0("plots/matchrate/EAS_classification_fraction.pdf"), height = 3, width = 4)
EAS_dat %>% group_by(group) %>% summarize(n = n(), mean = mean(length), sum = sum(length)) %>%
  mutate("rel_n" = n / sum(n), "rel_sum" = sum / sum(sum)) %>%
  filter(group != "other") %>%
  pivot_longer(cols = c(rel_n, rel_sum)) %>%
  ggplot(aes(x = name, y = value * 100, fill = group)) +
  theme_bw() + 
  lims(y = c(0, 70)) +
  labs(y = "Percentage of Inferred Archaic") +
  scale_fill_manual(values = c("darkviolet", "violet", "darkorange")) +
  geom_bar(stat = "identity", position = "stack")
dev.off()

# lets see how this changes if I change the thresholds
EAS_dat$groupND <- ifelse(EAS_dat$nea_affinity <= 0.3 & EAS_dat$den_affinity >= 0.7 & EAS_dat$ND10 == 0, "Dhigh", 
                        ifelse(EAS_dat$nea_affinity <= 0.3 & EAS_dat$den_affinity >= 0.4 & EAS_dat$den_affinity < 0.7 & EAS_dat$ND10 == 0, "Dlow",
                               ifelse(EAS_dat$nea_affinity >= 0.5 & EAS_dat$den_affinity <= 0.01  & EAS_dat$ND01 == 0, "Nea", "other")))

# great, I confirmed that this matches what I find in my pipeline
pdf(paste0("plots/matchrate/EAS_classificationND_fraction.pdf"), height = 3, width = 4)
EAS_dat %>% group_by(groupND) %>% summarize(n = n(), mean = mean(length), sum = sum(length)) %>%
  mutate("rel_n" = n / sum(n), "rel_sum" = sum / sum(sum)) %>%
  filter(groupND != "other") %>%
  pivot_longer(cols = c(rel_n, rel_sum)) %>%
  ggplot(aes(x = name, y = value * 100, fill = groupND)) +
  theme_bw() + 
  labs(y = "Percentage of Inferred Archaic") +
  lims(y = c(0, 75)) +
  scale_fill_manual(values = c("darkviolet", "violet", "darkorange")) +
  geom_bar(stat = "identity", position = "stack")
dev.off()

# lets see how this changes if I change the thresholds
EAS_dat$groupStrict <- ifelse(EAS_dat$nea_affinity <= 0.3 & EAS_dat$den_affinity >= 0.8, "Dhigh", 
                          ifelse(EAS_dat$nea_affinity <= 0.3 & EAS_dat$den_affinity >= 0.3 & EAS_dat$den_affinity < 0.8, "Dlow",
                                 ifelse(EAS_dat$nea_affinity >= 0.3 & EAS_dat$den_affinity <= 0.3, "Nea", "other")))

# great, I confirmed that this matches what I find in my pipeline
pdf(paste0("plots/matchrate/EAS_classificationStrict_fraction.pdf"), height = 3, width = 4)
EAS_dat %>% group_by(groupStrict) %>% summarize(n = n(), mean = mean(length), sum = sum(length)) %>%
  mutate("rel_n" = n / sum(n), "rel_sum" = sum / sum(sum)) %>%
  filter(groupStrict != "other") %>%
  pivot_longer(cols = c(rel_n, rel_sum)) %>%
  ggplot(aes(x = name, y = value * 100, fill = groupStrict)) +
  theme_bw() + 
  labs(y = "Percentage of Inferred Archaic") +
  lims(y = c(0, 75)) +
  scale_fill_manual(values = c("darkviolet", "violet", "darkorange")) +
  geom_bar(stat = "identity", position = "stack")
dev.off()


sum_th <- EAS_dat %>% group_by(groupND) %>% summarize(n = n(), mean = mean(length), sum = sum(length), archaic = sum(nea_overlap + den_overlap)) %>% mutate(thresh = "ND") %>% rename("group" = groupND)
sum_1 <- EAS_dat %>% group_by(group) %>% summarize(n = n(), mean = mean(length), sum = sum(length), archaic = sum(nea_overlap + den_overlap)) %>% mutate(thresh = "first") 
sum_strict <- EAS_dat %>% group_by(groupStrict) %>% summarize(n = n(), mean = mean(length), sum = sum(length), archaic = sum(nea_overlap + den_overlap)) %>% mutate(thresh = "strict") %>% rename("group" = groupStrict)
sums <- rbind(sum_th, sum_1) %>% rbind(sum_strict)

pdf(paste0("plots/matchrate/EAS_classification_Mb.pdf"), height = 6, width = 4)
sums %>%
  ggplot(aes(x = thresh, y = sum / 1000000, fill = thresh)) +
  geom_bar(stat = "identity") +
  labs(x = "filter", y = "Amount Recovered (Mb)", title = "EAS Synthetic Genome Classification") +
  facet_grid(group~., scales = "free")
dev.off()

pdf(paste0("plots/matchrate/EAS_classification_archaicSnps.pdf"), height = 6, width = 4)
sums %>%
  ggplot(aes(x = thresh, y = archaic, fill = thresh)) +
  geom_bar(stat = "identity") +
  labs(x = "filter", y = "# SNPs matching any Archaic", title = "EAS Synthetic Genome Classification") +
  facet_grid(group~., scales = "free")
dev.off()

ggplot(filter(EAS_dat, group == "other"), aes(x = nea_affinity, y = den_affinity)) + 
  geom_point() +
  labs(title = "Affinity of Unclassified Fragments") +
  theme_bw()

EAS_dat$groupND <- ifelse(EAS_dat$nea_affinity <= 0.3 & EAS_dat$den_affinity >= 0.8 & EAS_dat$ND10 == 0, "Dhigh", 
                          ifelse(EAS_dat$nea_affinity <= 0.3 & EAS_dat$den_affinity >= 0.4 & EAS_dat$den_affinity < 0.7 & EAS_dat$ND10 == 0, "Dlow",
                                 ifelse(EAS_dat$nea_affinity >= 0.5 & EAS_dat$den_affinity <= 0.01  & EAS_dat$ND01 == 0, "Nea", "other")))


dat %>% filter(den_affinity < 0.7 & den_affinity > 0.4) %>% select(ND01) %>% sum()
dat %>% filter(den_affinity < 0.7 & den_affinity > 0.4) %>% select(length) %>% sum() / 1000000


ggplot(EAS_dat, aes(x = group, y = length)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.1, height = 0), alpha = 0.2) +
  scale_y_log10() +
  theme_bw()

EAS_dat$super_group <- ifelse(EAS_dat$nea_affinity < 0.3 & EAS_dat$den_affinity >= 0.3, "Den", 
                        ifelse(EAS_dat$den_affinity < 0.3 & EAS_dat$nea_affinity >= 0.3, "Nea", "other"))

EAS_dat$groupLS <- ifelse(EAS_dat$ND10 == 0 & EAS_dat$ND01 + EAS_dat$ND11 >= 1, "Den", 
                          ifelse(EAS_dat$ND01 == 0 & EAS_dat$ND10 + EAS_dat$ND11 >= 1, "Nea",  "other"))
EAS_dat$super_group <- ifelse(EAS_dat$groupLS == "Den" & EAS_dat$den_affinity >= 0.7, "Dhigh",  
                              ifelse(EAS_dat$groupLS == "Den" & EAS_dat$den_affinity < 0.7, "Dlow", 
                                     ifelse(EAS_dat$groupLS == "Nea", "Nea", "other")))

filter(EAS_dat, comp >= 30) %>% group_by(super_group) %>% summarize(mean(length))

ggplot(EAS_dat, aes(x = super_group, y = length)) +
  geom_boxplot() + 
  geom_point(position = position_jitter(width = 0.1, height = 0), alpha = 0.2) +
  scale_y_log10()
  
stat <- filter(EAS_dat, den_affinity >= 0.3, nea_affinity < 0.3)
summary(lmodel <- lm(stat$length ~ stat$den_affinity))

pdf(paste0("plots/matchrate/EAS_den_length_correlation.pdf"), height = 4, width = 6)
ggplot(stat, aes(x = den_affinity, y = length)) +
  geom_point() +
  geom_abline(slope = -12080, intercept = 119684, color = "blue") +
  theme_bw()
dev.off()

stat <- filter(EAS_dat, den_affinity < 0.3, nea_affinity > 0.5)
summary(lmodel <- lm(stat$length ~ stat$nea_affinity))

pdf(paste0("plots/matchrate/EAS_nea_length_correlation.pdf"), height = 4, width = 6)
ggplot(stat, aes(x = nea_affinity, y = length)) +
  geom_point() +
  geom_abline(slope = 82865, intercept = 62179, color = "blue") +
  theme_bw()
dev.off()

# SAS
SAS_dat <- read.table("data/1000G/SAS/fragments/SAS_fragments_snps10_08_nonoverlapping_full.bed", header = T)
popname <- "SAS"
SAS_dat$length <- SAS_dat$end - SAS_dat$start

# countour plots
pdf(paste0("plots/matchrate/SAS_synthetic_comparable_contour.pdf"), height = 4, width = 4.5)
contour(kde2d(SAS_dat$nea_affinity,SAS_dat$den_affinity, lims=c(0,1,0,1),n=100),xaxs="i",yaxs="i",xlab=xlab,ylab=ylab,main=popname,levels=level1,las=1,cex.lab=1.5,cex.axis=1.4,cex.main=1.5,lty=5,labcex=0.8);abline(h=seq(0.2,0.8,0.2),v=seq(0.2,0.8,0.2),lty=3,col="gray")
contour(kde2d(SAS_dat$nea_affinity,SAS_dat$den_affinity, lims=c(0,1,0,1),n=100),xaxs="i",yaxs="i",xlab=xlab,ylab=ylab,main=popname,levels=level2,las=1,cex.lab=1.5,cex.axis=1.4,cex.main=1.5,lty=5,labcex=0.8,add=T);abline(h=seq(0.2,0.8,0.2),v=seq(0.2,0.8,0.2),lty=3,col="gray")
dev.off()

pdf(paste0("plots/matchrate/SAS_synthetic_comparable_histogram.pdf"), height = 4, width = 6)
SAS_dat %>% filter(den_affinity > 0.3, nea_affinity < 0.3) %>% 
  ggplot(aes(x = den_affinity)) +
  geom_histogram() +
  labs(title = "SAS") +
  theme_bw()
dev.off()

pdf(paste0("plots/matchrate/SAS_synthetic_comparable_allfrags_histogram.pdf"), height = 4, width = 6)
SAS_dat %>% 
  ggplot(aes(x = den_affinity)) +
  geom_histogram() +
  labs(title = "SAS") +
  theme_bw()
dev.off()

pdf(paste0("plots/matchrate/SAS_synthetic_comparable_histogram_ND11.pdf"), height = 4, width = 6)
ggplot(filter(SAS_dat, ND01>1, ND10 == 0), aes(x = den_affinity)) +
  geom_histogram() +
  labs(title = popname) +
  theme_bw()
dev.off()

nrow(SAS_dat)
sum(SAS_dat$length) / 1000000000
SAS_dat %>% filter(den_affinity > 0) %>% summarize(mean(den_affinity))
SAS_dat %>% filter(nea_affinity > 0) %>% summarize(mean(nea_affinity))

SAS_dat$super_group <- ifelse(SAS_dat$nea_affinity < 0.3 & SAS_dat$den_affinity >= 0.3, "Den", 
                              ifelse(SAS_dat$den_affinity < 0.3 & SAS_dat$nea_affinity >= 0.3, "Nea", "other"))
ggplot(SAS_dat, aes(x = super_group, y = length)) +
  geom_boxplot() + 
  geom_point(position = position_jitter(width = 0.1, height = 0), alpha = 0.2) +
  scale_y_log10()

SAS_dat$groupLS <- ifelse(SAS_dat$ND10 == 0 & SAS_dat$ND01 + SAS_dat$ND11 >= 1, "Den", 
                              ifelse(SAS_dat$ND01 == 0 & SAS_dat$ND10 + SAS_dat$ND11 >= 1, "Nea",  "other"))
SAS_dat$super_group <- ifelse(SAS_dat$groupLS == "Den" & SAS_dat$den_affinity >= 0.5, "Dhigh",  
                              ifelse(SAS_dat$groupLS == "Den" & SAS_dat$den_affinity < 0.5, "Dlow", 
                                     ifelse(SAS_dat$groupLS == "Nea", "Nea", "other")))

filter(SAS_dat, comp >= 30) %>% group_by(super_group) %>% summarize(mean(length))
#SAS_dat <- filter(SAS_dat, comp >= 30)
# wahoo!! here I can replicate what Laurits found (LASIDAD supplement page 88)
# but why is Nea still longer?? I think that there is some issue here....

contour(kde2d(SAS_dat$nea_affinity,SAS_dat$den_affinity, lims=c(0,1,0,1),n=100),xaxs="i",yaxs="i",xlab=xlab,ylab=ylab,main=popname,levels=level1,las=1,cex.lab=1.5,cex.axis=1.4,cex.main=1.5,lty=5,labcex=0.8);abline(h=seq(0.2,0.8,0.2),v=seq(0.2,0.8,0.2),lty=3,col="gray")
contour(kde2d(SAS_dat$nea_affinity,SAS_dat$den_affinity, lims=c(0,1,0,1),n=100),xaxs="i",yaxs="i",xlab=xlab,ylab=ylab,main=popname,levels=level2,las=1,cex.lab=1.5,cex.axis=1.4,cex.main=1.5,lty=5,labcex=0.8,add=T);abline(h=seq(0.2,0.8,0.2),v=seq(0.2,0.8,0.2),lty=3,col="gray")


SAS_all_dat <- read.table("data/1000G/SAS/fragments/SAS_fragments_snps10_08_annotated.bed", header = T)
SAS_all_dat$length <- SAS_all_dat$end - SAS_all_dat$start
SAS_all_dat <- filter(SAS_all_dat, comp >= 30)

SAS_all_dat$super_group <- ifelse(SAS_all_dat$nea_affinity < 0.3 & SAS_all_dat$den_affinity >= 0.3, "Den", 
                              ifelse(SAS_all_dat$den_affinity < 0.3 & SAS_all_dat$nea_affinity >= 0.3, "Nea", "other"))
SAS_all_dat$super_groupStrict <- ifelse(SAS_all_dat$nea_affinity > SAS_all_dat$den_affinity, "Nea", 
                                  ifelse(SAS_all_dat$nea_affinity < SAS_all_dat$den_affinity, "Den", "Other"))
SAS_all_dat$group <- ifelse(SAS_all_dat$nea_affinity < 0.3 & SAS_all_dat$den_affinity >= 0.7, "Dhigh", 
                        ifelse(SAS_all_dat$nea_affinity < 0.3 & SAS_all_dat$den_affinity > 0.3, "Dlow",
                               ifelse(SAS_all_dat$nea_affinity >= 0.5 & SAS_all_dat$den_affinity < 0.3, "Nea", "other")))
SAS_all_dat$groupLS <- ifelse(SAS_all_dat$ND10 == 0 & SAS_all_dat$ND01 + SAS_all_dat$ND11 >= 1, "Den", 
                                  ifelse(SAS_all_dat$ND01 == 0 & SAS_all_dat$ND10 + SAS_all_dat$ND11 >= 1, "Nea",  "other"))

SAS_all_dat %>% filter(name == "NA18942") %>% group_by(groupLS) %>% summarize(mean(length))
SAS_all_dat %>% filter(name == "NA18942") %>% group_by(super_group) %>% summarize(mean(length))
SAS_all_dat %>% group_by(super_groupStrict) %>% summarize(mean(length))
SAS_all_dat %>%  filter(comp > 10, length >= 10000, name ==) %>% group_by(ND_type) %>% summarize(mean(length))

SAS_all_dat %>% filter(comp > 10, length >= 30000) %>% group_by(super_group) %>% summarize(mean(length))

# LASIDAD (specifically multipulse and singlel pulse, but plot all) 
# no multipulse
sp_dat <- read.table("data/LASIDAD/singlePulse/fragments/singlePulse_fragments_snps10_08_nonoverlapping_full.bed", header = T)

popname <- "singlePulse"

pdf(paste0("plots/matchrate/sp_synthetic_comparable_contour.pdf"), height = 4, width = 4.5)
contour(kde2d(sp_dat$nea_affinity,sp_dat$den_affinity, lims=c(0,1,0,1),n=100),xaxs="i",yaxs="i",xlab=xlab,ylab=ylab,main=popname,levels=level1,las=1,cex.lab=1.5,cex.axis=1.4,cex.main=1.5,lty=5,labcex=0.8);abline(h=seq(0.2,0.8,0.2),v=seq(0.2,0.8,0.2),lty=3,col="gray")
contour(kde2d(sp_dat$nea_affinity,sp_dat$den_affinity, lims=c(0,1,0,1),n=100),xaxs="i",yaxs="i",xlab=xlab,ylab=ylab,main=popname,levels=level2,las=1,cex.lab=1.5,cex.axis=1.4,cex.main=1.5,lty=5,labcex=0.8,add=T);abline(h=seq(0.2,0.8,0.2),v=seq(0.2,0.8,0.2),lty=3,col="gray")
dev.off()

nrow(sp_dat)
sum(sp_dat$length) / 1000000000

sp_dat %>% filter(den_affinity > 0) %>% summarize(mean(den_affinity))
sp_dat %>% filter(nea_affinity > 0) %>% summarize(mean(nea_affinity))


mp_dat <- read.table("data/LASIDAD/multiPulse/fragments/multiPulse_fragments_snps10_08_nonoverlapping_full.bed", header = T)
popname <- "multiPulse"

# countour plots
pdf(paste0("plots/matchrate/mp_synthetic_comparable_contour.pdf"), height = 4, width = 4.5)
contour(kde2d(mp_dat$nea_affinity,mp_dat$den_affinity, lims=c(0,1,0,1),n=100),xaxs="i",yaxs="i",xlab=xlab,ylab=ylab,main=popname,levels=level1,las=1,cex.lab=1.5,cex.axis=1.4,cex.main=1.5,lty=5,labcex=0.8);abline(h=seq(0.2,0.8,0.2),v=seq(0.2,0.8,0.2),lty=3,col="gray")
contour(kde2d(mp_dat$nea_affinity,mp_dat$den_affinity, lims=c(0,1,0,1),n=100),xaxs="i",yaxs="i",xlab=xlab,ylab=ylab,main=popname,levels=level2,las=1,cex.lab=1.5,cex.axis=1.4,cex.main=1.5,lty=5,labcex=0.8,add=T);abline(h=seq(0.2,0.8,0.2),v=seq(0.2,0.8,0.2),lty=3,col="gray")
dev.off()

pdf(paste0("plots/matchrate/mp_synthetic_comparable_histogram.pdf"), height = 4, width = 6)
mp_dat %>% filter(den_affinity > 0.3, nea_affinity < 0.3) %>% 
  ggplot(aes(x = den_affinity)) +
  geom_histogram() +
  labs(title = "mp") +
  theme_bw()
dev.off()

nrow(mp_dat)
sum(mp_dat$length) / 1000000000
mp_dat %>% filter(den_affinity > 0) %>% summarize(mean(den_affinity))
mp_dat %>% filter(nea_affinity > 0) %>% summarize(mean(nea_affinity))

pdf(paste0("plots/matchrate/mp_synthetic_comparable_allfrags_histogram.pdf"), height = 4, width = 6)
mp_dat %>% 
  ggplot(aes(x = den_affinity)) +
  geom_histogram() +
  labs(title = "mp") +
  theme_bw()
dev.off()

pdf(paste0("plots/matchrate/mp_synthetic_comparable_histogram_ND11.pdf"), height = 4, width = 6)
ggplot(filter(mp_dat, ND01>1, ND10 == 0), aes(x = den_affinity)) +
  geom_histogram() +
  labs(title = "mp") +
  theme_bw()
dev.off()

#i <- 4
pops=c("Central","East","multiPulse","North","North-East", "singlePulse", "South", "West")
dat_summary <- data.frame()
for(i in 1:length(pops)){
  pop = pops[i]
  popname = pop
  print(pop)
  dat <- read.table(paste0("data/LASIDAD/", pop, "/fragments/", pop, "_fragments_snps10_08_nonoverlapping_full.bed"), header = T)
  #dat <- filter(dat, comp >= 10)
  tmp_summary = data.frame("pop" = pop, "nfrags" = nrow(dat), "length" = sum(dat$length), "nea_mean" = mean((dat[dat$nea_affinity > 0, ])$nea_affinity), "den_mean" = mean((dat[dat$den_affinity > 0, ])$den_affinity))
  dat_summary <- rbind(dat_summary, tmp_summary)
  
  # countour plots
  pdf(paste0("plots/matchrate/", pop, "_synthetic_comparable_contour.pdf"), height = 4, width = 4.5)
  contour(kde2d(dat$nea_affinity,dat$den_affinity, lims=c(0,1,0,1),n=100),xaxs="i",yaxs="i",xlab=xlab,ylab=ylab,main=popname,levels=level1,las=1,cex.lab=1.5,cex.axis=1.4,cex.main=1.5,lty=5,labcex=0.8);abline(h=seq(0.2,0.8,0.2),v=seq(0.2,0.8,0.2),lty=3,col="gray")
  contour(kde2d(dat$nea_affinity,dat$den_affinity, lims=c(0,1,0,1),n=100),xaxs="i",yaxs="i",xlab=xlab,ylab=ylab,main=popname,levels=level2,las=1,cex.lab=1.5,cex.axis=1.4,cex.main=1.5,lty=5,labcex=0.8,add=T);abline(h=seq(0.2,0.8,0.2),v=seq(0.2,0.8,0.2),lty=3,col="gray")
  dev.off()

  pdf(paste0("plots/matchrate/", pop, "_synthetic_comparable_histogram.pdf"), height = 4, width = 6)
  p <- ggplot(filter(dat, den_affinity > 0.3, nea_affinity < 0.3), aes(x = den_affinity)) +
    geom_histogram() +
    labs(title = popname) +
    theme_bw()
  print(p)
  dev.off()
  
  pdf(paste0("plots/matchrate/", pop, "_synthetic_comparable_histogram_ND11.pdf"), height = 4, width = 6)
  p <- ggplot(filter(dat, ND01>1, ND10 == 0), aes(x = den_affinity)) +
    geom_histogram() +
    labs(title = popname) +
    theme_bw()
  print(p)
  dev.off()
}


# using the data that laurits generated:
path <- "/global/scratch/p2p3/pl1_moorjani/lauritsskov2/ArchaicSegments/Results/04Summary/SegmentSummary/ArtificialGenome/datafiles/"
dat_LS <- read.table(paste0(path, "Browning_ND11_30_0.8_phased.txt"), header = T) # "Browning_ND11_30_0.8_phased.txt"
dat <- filter(dat, POP == "CHB")
mean(dat$NMATCH[dat$NMATCH > 0])
mean(dat$DMATCH[dat$DMATCH > 0])
dat$nea_affinity <- dat$Altai / dat$all_snps
mean(dat$nea_affinity[dat$nea_affinity > 0])
dat$den_affinity <- dat$Denisova / dat$all_snps
mean(dat$den_affinity[dat$den_affinity > 0])




# using the data that laurits generated:
path <- "/global/scratch/p2p3/pl1_moorjani/lauritsskov2/ArchaicSegments/Results/04Summary/SegmentSummary/ArtificialGenome/datafiles/"
dat_LS <- read.table(paste0(path, "Browning_ND11_30_0.8_phased.txt"), header = T) # "Browning_ND11_30_0.8_phased.txt"
dat <- filter(dat, POP == "CHB")
mean(dat$NMATCH[dat$NMATCH > 0])
mean(dat$DMATCH[dat$DMATCH > 0])
dat$nea_affinity <- dat$Altai / dat$all_snps
mean(dat$nea_affinity[dat$nea_affinity > 0])
dat$den_affinity <- dat$Denisova / dat$all_snps
mean(dat$den_affinity[dat$den_affinity > 0])
