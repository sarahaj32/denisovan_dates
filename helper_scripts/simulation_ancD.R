library(tidyverse)
library(DEoptim)


args <- commandArgs(trailingOnly = TRUE)
file <- args[[1]]
archaic <- args[[2]]
seed <- args[[3]]
min_dist <- args[[4]]

# file <- "../simdat/denisovanNeanderthal_simple_1/DEN/curve/results/denisovanNeanderthal_simple_100NAMH_ibdmix_all_2_D.txt"
# min_dist = 0.02
print('args read in:')
print(args)

# set plotting settings and variables
theme_set(theme_bw() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
            theme(axis.text = element_text(size = 10)))

# helper function to find the best fit parameters
opt_fit_C <- function(dat){
  # Pick optimized parameters for model 1: A*exp(m1 * -x) + C
  dist <- dat$distance / 100		# updated x values
  wcorr <- dat$D		# updated y values
  
  fm1 <- function(x) x[1]*exp(-x[2] * dist) + x[3]
  fm2 <- function(x) sum((wcorr-fm1(x))^2)
  fm3 <- DEoptim(fm2, lower=c(0.00001,500, -5), upper=c(10, 2000, 5), control=list(trace=FALSE))
  par1 <- fm3$optim$bestmem
  
  # parameters for y ~ Ae-mt + C
  names(par1) <- c("A", "m", "C")
  
  model_exp <- nls(wcorr ~ (A*exp(-m * dist) + C), start=par1, control=list(maxiter=1000, warnOnly=TRUE))
  
  
  return(coef(model_exp))
}

# a helper function to create the input to the jacknife
# calculates the fit to the data after excluding each chromosome
jacknife <- function(dat){
  # first get the fit from the data averaged across all chromosomes
  jack_means <- dat %>% group_by(distance) %>% summarize(D_sum = sum(D_sum), count = sum(count), D = D_sum / count)
  
  full_results <- opt_fit_C(jack_means)
  jack_results <- data.frame("inc" = "all", "inf_C" = full_results[3], "inf_m" = full_results[2], "inf_a" = full_results[1], "snps" = "NA")
  
  for (i in unique(dat$chrom)){
    #rm(a)
    # remove each chromosome, and calculate the mean D
    jack_dat <- filter(dat, chrom != i)
    #print(nrow(jack_dat))
    jack_means <- jack_dat %>% group_by(distance) %>% summarize(D_sum = sum(D_sum), count = sum(count), D = D_sum / count)
    
    jack <- opt_fit_C(jack_means)
    inf_C <- jack[3]
    inf_m <- jack[2]
    #print(inf_m)
    inf_a <- jack[1]
    
    snps = filter(dat, chrom == i) %>% dplyr::select(snps) %>% unique()
    
    jack_out <- data.frame("inc" = i, "inf_m" = inf_m, "inf_a" = inf_a, "inf_C" = inf_C, "snps" = snps)
    jack_results <- rbind(jack_results, jack_out)
  }
  return(jack_results)
}

# separate the file name, and extrapolate across all chromosomes
pieces <- unlist(strsplit(file, split = "_"))
file_base <- paste0(pieces[1:(length(pieces) - 2)], collapse = "_")
print(file_base)
names <- unlist(strsplit(file_base, split = "/"))
name <- names[length(names)]

# read in all the data across chromosomes
curve_dat <- data.frame()
for (chr in 1:20){
  print(chr)
  chr_file <- paste0(c(file_base, chr, "D.txt"), collapse = "_")
  tmp_curve <- read.table(chr_file)
  colnames(tmp_curve) <- c("min", "max", "D", "D_sum", "count", "snps")  
  tmp_curve$chrom <- chr
  curve_dat <- rbind(curve_dat, tmp_curve)
}

# convert distance to centimorgans
curve_dat$distance <- curve_dat$min / 1000000
# apply minimum cutoff
curve_dat <- filter(curve_dat, distance > min_dist)

# plot the curves for each chromosome
p_chrom <- ggplot(curve_dat, aes(x = distance, y = D)) +
  geom_point(size = 0.1) +
  labs(x = "Distance (in cM)", y = "Ancestry D", title = paste0("Each Chromosome Plotted Separately\n", name))

pdf(paste0("simdat/results/", name, "_", archaic, "_", seed, "_min", min_dist, "_chr.pdf"), height = 4, width = 6)
plot(p_chrom)
dev.off()


print("Jacknifing curve")
# jacknife the fits:
jack_curve_fit <- jacknife(curve_dat)
# prepare input for the weighted jacknife analysis
jin <- paste0(file_base, "_min", min_dist, ".jin")
jout <- paste0(file_base, "_min", min_dist, ".jout")
      
# write out the input for Priya's jacknife function
jack_curve_fit %>% filter(inc != "all") %>% 
  dplyr::select(inc, snps, inf_m) %>%
  write.table(jin, sep = "\t", quote = F, row.names = F, col.names = F)

full_mean <- filter(jack_curve_fit, inc == "all") %>% dplyr::select(inf_m) %>% unlist()
full_a <- filter(jack_curve_fit, inc == "all") %>% dplyr::select(inf_a) %>% unlist()
full_C <- filter(jack_curve_fit, inc == "all") %>% dplyr::select(inf_C) %>% unlist()

# run weighted jacknife
print(paste0("dowtjack -i ",jin, " -m ", full_mean,  " -o ", jout))
system(paste0("dowtjack -i ",jin, " -m ", full_mean,  " -o ", jout))
print("DONE")

jack_results <- read.table(jout)
print(jack_results)
inf_m <- jack_results[1, 1]
SE <- jack_results[1, 2]
if (inf_m == "mean:")
  inf_m <- jack_results[1, 2]
  SE <- jack_results[1, 4]

print(inf_m)
print(full_mean)
print("OK")
# plot the fit as a line
print(full_a)
print(inf_m)
print(full_C)
print(head(curve_dat))
curve_dat$fit <- full_a * exp(-inf_m * curve_dat$distance / 100) + full_C
print("FIT")

fit_p <- curve_dat  %>% group_by(distance, fit) %>% summarize(D = mean(D)) %>% 
  ggplot(aes(x = distance, y = D)) +
  geom_point() +
  theme_bw() + 
  labs(x = "distance (in cM)", y = "Ancestry D", title = paste0("Inferred Slope = ", round(inf_m), " jacknifed_SD = ", SE, "\n", name)) +
  geom_line(aes(x = distance, y = fit), color = "red", size = 0.75)

pdf(paste0("simdat/results/", name, "_", archaic, "_", seed, "_min", min_dist, "_fit.pdf"), height = 4, width = 6)
plot(fit_p)
dev.off()

write.table(data.frame(name = name, "inf_m" = inf_m, "SE" = SE, "full_mean" = full_mean, "full_a" = full_a, "full_C" = full_C), 
            paste0("simdat/results/", name, "_", archaic, "_", seed, "_min", min_dist, "_inferred_parameters.txt"), quote = F, row.names = F, sep = "\t")
# write out the inferred parameters

# also write out the average truth and fit
curve_dat %>% select(distance, D, fit) %>%
  mutate(name = name) %>% 
  write.table(paste0("simdat/results/", name, "_", archaic, "_", seed, "_min", min_dist, "_fit.txt"), quote = F, row.names = F, sep = "\t")
