library(data.table)
library(dplyr)
library(ggplot2)
library(ggsci)
library(ggjoy)
library(gridExtra)
library(RCurl)
library(pscl)

# import the data
URL <- "https://raw.githubusercontent.com/rmccoy7541/aneuploidy-analysis/master/data/aaa3337-McCoy-SM.table_S2.csv" 
data <- fread(URL, sep = ",", header = T, na.strings = c("NO_CALL", "NA"))
data[, case_embryoid := paste(case, embryoid, sep = "_")]
setorder(data, case, embryoid)

# count chromosomes per sample
count_chroms <- function(ploidy_calls) {
  return(sum(as.numeric(unlist(strsplit(gsub("H", "", unname(unlist(ploidy_calls))), ""))), na.rm = T))
}

data[, chrom_count := unlist(lapply(1:nrow(data), function(x) count_chroms(data[x, 7:29, with = F])))]

# append meiotic flag to corresponding maternal chromosome gains
matrix <- as.matrix(data[, 7:29])
index <- as.matrix(data[, 30:52] == 1)
matrix[index][!is.na(matrix[index]) & (matrix[index] %in% c("H210", "H200", "H201"))] <- paste(matrix[index][!is.na(matrix[index]) & (matrix[index] %in% c("H210", "H200", "H201"))], "mbph", sep = "_")

for (i in 7:29) {
  data[, i] <- matrix[, i - 6]
}

# append meiotic flag to corresponding paternal chromosome gains
matrix <- as.matrix(data[, 7:29])
index <- as.matrix(data[, 53:75] == 1)
matrix[index][!is.na(matrix[index]) & (matrix[index] %in% c("H120", "H020", "H002", "H102"))] <- paste(matrix[index][!is.na(matrix[index]) & (matrix[index] %in% c("H120", "H020", "H002", "H102"))], "pbph", sep = "_")
                                                                                                    
for (i in 7:29) {
  data[, i] <- matrix[, i - 6]
}

########################################################################
# Functions to count number of various chromosomal patterns per sample
########################################################################

count_m_monosomy <- function(ploidy_calls) {
  return(sum(unlist(ploidy_calls) %in% c("H010", "H001"), na.rm = T))
}

count_p_monosomy <- function(ploidy_calls) {
  return(sum(unlist(ploidy_calls) == "H100", na.rm = T))
}

count_m_trisomy <- function(ploidy_calls) {
  return(sum(unlist(ploidy_calls) %in% c("H210", "H201", "H210_mbph", "H201_mbph"), na.rm = T))
}

count_p_trisomy <- function(ploidy_calls) {
  return(sum(unlist(ploidy_calls) %in% c("H120", "H102", "H012", "H021", "H111", "H120_pbph", "H102_pbph", "H021_pbph"), na.rm = T))
}

count_nullisomy <- function(ploidy_calls) {
  return(sum(unlist(ploidy_calls) == "H000", na.rm = T))
}

count_m_upd <- function(ploidy_calls) {
  return(sum(unlist(ploidy_calls) %in% c("H200", "H200_mbph"), na.rm = T))
}

count_p_upd <- function(ploidy_calls) {
  return(sum(unlist(ploidy_calls) %in% c("H020", "H002", "H011", "H020_pbph", "H002_pbph"), na.rm = T))
}

count_disomy <- function(ploidy_calls) {
  return(sum(unlist(ploidy_calls) %in% c("H110", "H101"), na.rm = T))
}

count_meiotic <- function(ploidy_calls) {
  return(sum(grepl("bph", unlist(ploidy_calls)), na.rm = T))
}

count_na <- function(ploidy_calls) {
  return(sum(is.na(unlist(ploidy_calls)), na.rm = T))
}

# apply above functions to count various forms of aneuploidy

data[, m_monosomy := unlist(lapply(1:nrow(data), function(x) count_m_monosomy(data[x, 7:29, with = F])))]

data[, p_monosomy := unlist(lapply(1:nrow(data), function(x) count_p_monosomy(data[x, 7:29, with = F])))]

data[, m_trisomy := unlist(lapply(1:nrow(data), function(x) count_m_trisomy(data[x, 7:29, with = F])))]

data[, p_trisomy := unlist(lapply(1:nrow(data), function(x) count_p_trisomy(data[x, 7:29, with = F])))]

data[, disomy := unlist(lapply(1:nrow(data), function(x) count_disomy(data[x, 7:29, with = F])))]

data[, nullisomy := unlist(lapply(1:nrow(data), function(x) count_nullisomy(data[x, 7:29, with = F])))]

data[, na := unlist(lapply(1:nrow(data), function(x) count_na(data[x, 7:29, with = F])))]

data[, null_na := na + nullisomy]

data[, meiotic := unlist(lapply(1:nrow(data), function(x) count_meiotic(data[x, 7:29, with = F])))]

data[, m_upd := unlist(lapply(1:nrow(data),  function(x) count_m_upd(data[x, 7:29, with = F])))]

data[, p_upd := unlist(lapply(1:nrow(data), function(x) count_p_upd(data[x, 7:29, with = F])))]

# remove replicate samples, identified by Alan Handyside as consecutive samples with identical aneuploidies
handyside_dups <- fread("https://raw.githubusercontent.com/rmccoy7541/tripolar_mitosis/master/data/dups.txt") %>%
  setnames(., c("case", "embryoid"))
handyside_dups[, case_embryoid := paste(case, embryoid, sep = "_")]

data <- data[!(case_embryoid %in% handyside_dups$case_embryoid)]

# remove samples with >= 10 nullisomies and low-confidence calls
data <- data[null_na < 11]

########################################################################
# Analysis of diploid tripolar samples
########################################################################

# identify diploid tripolar samples
diploid_tripolar <- data[disomy >= 1 & nullisomy >= 1 & p_monosomy >= 1 & m_monosomy >= 1 & 
                         p_trisomy == 0 & m_trisomy == 0 & p_upd == 0 & m_upd == 0]$case_embryoid

data[, is_diploid_tripolar := (case_embryoid %in% diploid_tripolar)]

# calculate frequency at different stages
table(data[is_diploid_tripolar == T]$sample_type) / table(data$sample_type)

diploid_tripolar_summary <- data.table(group_by(data[sample_type == "blastomere"], case) %>%
  summarise(., dip_tri = sum(is_diploid_tripolar), n_dip_tri = sum(!is_diploid_tripolar), maternal_age = unique(maternal_age)))

# discovery subsample (from 2015 paper)
geno_discovery <- fread("~/tripolar_mitosis/rs2305957_f.ped")
setnames(geno_discovery, "V1", "case")
geno_discovery[, gt := as.numeric(NA)]
geno_discovery[V7 == "G" & V8 == "G", gt := 0]
geno_discovery[V7 == "A" & V8 == "G", gt := 1]
geno_discovery[V7 == "A" & V8 == "A", gt := 2]
geno_discovery <- merge(geno_discovery, diploid_tripolar_summary, "case")

# associate with PLK4 genotype
summary(glm(data = geno_discovery, formula = cbind(dip_tri, n_dip_tri) ~ gt, family = quasibinomial))
summary(glm(data = geno_discovery, formula = cbind(dip_tri, n_dip_tri) ~ gt + maternal_age, family = quasibinomial))

# validation subsample (from 2015 paper)
geno_validation <- fread("~/tripolar_mitosis/rs2305957_f_validation.ped")
setnames(geno_validation, "V1", "case")
geno_validation[, gt := as.numeric(NA)]
geno_validation[V7 == "G" & V8 == "G", gt := 0]
geno_validation[V7 == "A" & V8 == "G", gt := 1]
geno_validation[V7 == "A" & V8 == "A", gt := 2]
geno_validation <- merge(geno_validation, diploid_tripolar_summary, "case")

# associate with PLK4 genotype
summary(glm(data = geno_validation, formula = cbind(dip_tri, n_dip_tri) ~ gt, family = quasibinomial))
summary(glm(data = geno_validation, formula = cbind(dip_tri, n_dip_tri) ~ gt + maternal_age, family = quasibinomial))

# combine discovery and validation subsamples; test for association
full_model <- glm(data = rbind(geno_discovery, geno_validation), formula = cbind(dip_tri, n_dip_tri) ~ gt, family = quasibinomial)
summary(full_model)
full_model_age <- glm(data = rbind(geno_discovery, geno_validation), formula = cbind(dip_tri, n_dip_tri) ~ gt + maternal_age, family = quasibinomial)
summary(full_model_age)

# fit a hurdle model
m0 <- glm(data = rbind(geno_discovery, geno_validation), formula = dip_tri ~ gt, family = "poisson", weights = (dip_tri + n_dip_tri))
m1 <- hurdle(data = rbind(geno_discovery, geno_validation), formula = dip_tri ~ gt, dist = "poisson", zero.dist = "binomial", weights = (dip_tri + n_dip_tri))
vuong(m0, m1)

# plot the results
dt <- rbind(geno_discovery, geno_validation)[!is.na(gt)]

group_by(dt, gt) %>% summarise(., mean = mean(dip_tri / (dip_tri + n_dip_tri)))

zero_summary <- data.table(group_by(dt, gt) %>%
                             summarise(., prop_zero = sum(dip_tri == 0) / n(), n = n()))

zero_summary[, se := sqrt((prop_zero * (1 - prop_zero)) / n) ]

zero_plot <- ggplot(data = zero_summary, aes(x = as.factor(gt), y = prop_zero, fill = as.factor(gt))) + 
  geom_bar(stat = "identity") +
  geom_errorbar(aes(x = as.factor(gt), ymin = prop_zero - se, ymax = prop_zero + se), width = 0.3) +
  theme_bw() +
  xlab("Genotype") +
  ylab("Prop. w/ Zero Diploid Tripolar") +
  theme(panel.grid = element_blank(), legend.position = "none") +
  scale_fill_brewer(palette = "Dark2")

nonzero_plot <- ggplot(data = dt[dip_tri > 0 & (dip_tri + n_dip_tri) >= 3], aes(x = as.factor(gt), y = dip_tri / (dip_tri + n_dip_tri))) +
  geom_boxplot(aes(fill = as.factor(gt))) +
  theme_bw() +
  xlab("Genotype") +
  ylab("Prop. Affected Among Non-Zero Cases") +
  theme(panel.grid = element_blank(), legend.position = "none") +
  scale_fill_brewer(palette = "Dark2")

grid.arrange(nonzero_plot, zero_plot, ncol = 1)

# function to round to nearest arbitrary integer
mround <- function(x, base){ 
  base * round(x / base) 
} 

by_age_gg <- data.table(group_by(dt[gt == 0], age = mround(maternal_age, 2)) %>% 
                          summarise(., mean = mean(dip_tri / (dip_tri + n_dip_tri)), n = n(), gt = 0))
by_age_ag <- data.table(group_by(dt[gt == 1], age = mround(maternal_age, 2)) %>% 
                          summarise(., mean = mean(dip_tri / (dip_tri + n_dip_tri)), n = n(), gt = 1))
by_age_aa <- data.table(group_by(dt[gt == 2], age = mround(maternal_age, 2)) %>% 
                          summarise(., mean = mean(dip_tri / (dip_tri + n_dip_tri)), n = n(), gt = 2))
by_age <- rbind(by_age_gg, by_age_ag, by_age_aa)

by_age[, se := sqrt( ((mean * (1 - mean)) / n) )]

ggplot(data = by_age, aes(x = age, y = mean)) +
  geom_point(aes(color = as.factor(gt))) +
  geom_line(aes(color = as.factor(gt))) +
  geom_ribbon(aes(fill = as.factor(gt), x = age, ymin = mean - se, ymax = mean + se), alpha = 0.2) +
  coord_cartesian(ylim = c(0, 0.25)) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none") +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2")

########################################################################
# Tabulate the counts (and percentages) of samples falling into various categories
# also calculate the ratio of percentages: (% in blastomere samples) / (% in TE samples)
########################################################################

# euploid
nrow(data[disomy == 23]) / nrow(data)

(nrow(data[sample_type == "blastomere" & disomy == 23]) / nrow(data[sample_type == "blastomere"])) / 
(nrow(data[sample_type == "TE" & disomy == 23]) / nrow(data[sample_type == "TE"]))

euploid <- data[disomy == 23]$case_embryoid

# euploid / na
nrow(data[(disomy != 23) & ((disomy + na) == 23)]) / nrow(data)

(nrow(data[sample_type == "blastomere" & (disomy != 23) & ((disomy + na) == 23)]) / nrow(data[sample_type == "blastomere"])) / 
(nrow(data[sample_type == "TE" & (disomy != 23) & ((disomy + na) == 23)]) / nrow(data[sample_type == "TE"]))

# one to three aneuploid chromosomes
nrow(data[sample_type == "blastomere" & ((disomy + na) >= 20) & ((disomy + na) != 23)])
nrow(data[sample_type == "TE" & (disomy + na) >= 20 & (disomy + na) != 23])

nrow(data[sample_type == "blastomere" & ((disomy + na) >= 20) & ((disomy + na) != 23) & (meiotic > 0)])
nrow(data[sample_type == "TE" & ((disomy + na) >= 20) & ((disomy + na) != 23) & (meiotic > 0)])

nrow(data[sample_type == "blastomere" & ((disomy + na) >= 20) & ((disomy + na) != 23) & (meiotic == 0) & (p_monosomy > 0 | p_trisomy > 0)])
nrow(data[sample_type == "TE" & ((disomy + na) >= 20) & ((disomy + na) != 23) & (meiotic == 0) & (p_monosomy > 0 | p_trisomy > 0)])

# one to nine aneuploid chromosomes
nrow(data[(disomy + na) != 23 & (disomy + na) >= 14]) / nrow(data)

(nrow(data[sample_type == "blastomere" & (disomy + na) != 23 & (disomy + na) >= 14]) / nrow(data[sample_type == "blastomere"])) /
(nrow(data[sample_type == "TE" & (disomy + na) != 23 & (disomy + na) >=14]) / nrow(data[sample_type == "TE"])) 

aneuploid <- data[(disomy + na) != 23 & (disomy + na) >= 14]$case_embryoid

nrow(data[sample_type == "blastomere" & (disomy + na) != 23 & (disomy + na) >= 14 & meiotic > 0]) # meiotic
nrow(data[sample_type == "blastomere" & (disomy + na) != 23 & (disomy + na) >= 14 & meiotic == 0 & (p_monosomy > 0 | p_trisomy > 0)]) # mitotic only
nrow(data[sample_type == "TE" & (disomy + na) != 23 & (disomy + na) >= 14 & meiotic > 0]) # meiotic
nrow(data[sample_type == "TE" & (disomy + na) != 23 & (disomy + na) >= 14 & meiotic == 0 & (p_monosomy > 0 | p_trisomy > 0)]) # mitotic only

# karyotype-wide
nrow(data[(disomy + na) != 23 & (disomy + na) < 14]) / (nrow(data))

(nrow(data[sample_type == "blastomere" & (disomy + na) != 23 & (disomy + na) < 14]) / nrow(data[sample_type == "blastomere"])) /
(nrow(data[sample_type == "TE" & (disomy + na) != 23 & (disomy + na) < 14]) / nrow(data[sample_type == "TE"]))

karyotype_wide <- data[(disomy + na) != 23 & (disomy + na) < 14]$case_embryoid

nrow(data[sample_type == "blastomere" & (disomy + na) != 23 & (disomy + na) < 14 & meiotic > 0]) # meiotic
nrow(data[sample_type == "blastomere" & (disomy + na) != 23 & (disomy + na) < 14 & meiotic == 0 & (p_monosomy > 0 | p_trisomy > 0)]) # mitotic only
nrow(data[sample_type == "TE" & (disomy + na) != 23 & (disomy + na) < 14 & meiotic > 0]) # meiotic
nrow(data[sample_type == "TE" & (disomy + na) != 23 & (disomy + na) < 14 & meiotic == 0 & (p_monosomy > 0 | p_trisomy > 0)]) # mitotic only

# mitotic karyotype-wide
nrow(data[(disomy + na) != 23 & (disomy + na) < 14 & (p_trisomy > 0 | p_monosomy > 0 | m_trisomy > 0) & meiotic == 0])

nrow(data[(disomy + na) != 23 & (disomy + na) < 14 & (p_trisomy > 0 | p_monosomy > 0 | m_trisomy > 0) & meiotic == 0]) /
nrow(data[(disomy + na) != 23 & (disomy + na) < 14])

########################################################################
# Count various forms of karyotype-wide abnormality
########################################################################

# diploid tripolar
table(data[is_diploid_tripolar == T]$sample_type)
table(data[is_diploid_tripolar == T]$sample_type) / table(data$sample_type)

# digynic triploid / near-triploid
nrow(data[case_embryoid %in% karyotype_wide & (m_trisomy + null_na) >= 20 & !(case_embryoid %in% diploid_tripolar)])
digynic_triploid <- data[case_embryoid %in% karyotype_wide & (m_trisomy + null_na) >= 20]$case_embryoid
data[, is_digynic_triploid := case_embryoid %in% digynic_triploid]
table(data[is_digynic_triploid == T]$sample_type)
table(data[is_digynic_triploid == T]$sample_type) / table(data$sample_type)

# diandric triploid / near-triploid
nrow(data[case_embryoid %in% karyotype_wide & (p_trisomy + null_na) >= 20 & !(case_embryoid %in% diploid_tripolar)])
diandric_triploid <- data[case_embryoid %in% karyotype_wide & (p_trisomy + null_na) >= 20]$case_embryoid
data[, is_diandric_triploid := case_embryoid %in% diandric_triploid]
table(data[is_diandric_triploid == T]$sample_type)
table(data[is_diandric_triploid == T]$sample_type) / table(data$sample_type)

# maternal haploid / near-haploid
nrow(data[case_embryoid %in% karyotype_wide & (p_monosomy + null_na) >= 20 & !(case_embryoid %in% diploid_tripolar)])
maternal_haploid <- data[case_embryoid %in% karyotype_wide & (p_monosomy + null_na) >= 20 & !(case_embryoid %in% diploid_tripolar)]$case_embryoid
data[, is_maternal_haploid := case_embryoid %in% maternal_haploid]
table(data[is_maternal_haploid == T]$sample_type)
table(data[is_maternal_haploid == T]$sample_type) / table(data$sample_type)

# paternal haploid / near-haploid
nrow(data[case_embryoid %in% karyotype_wide & (m_monosomy + null_na) >= 20 & !(case_embryoid %in% diploid_tripolar)])
paternal_haploid <- data[case_embryoid %in% karyotype_wide & (m_monosomy + null_na) >= 20 & !(case_embryoid %in% diploid_tripolar)]$case_embryoid
data[, is_paternal_haploid := case_embryoid %in% paternal_haploid]
table(data[is_paternal_haploid == T]$sample_type)
table(data[is_paternal_haploid == T]$sample_type) / table(data$sample_type)

# multiple maternal gain
nrow(data[case_embryoid %in% karyotype_wide & m_monosomy == 0 & p_monosomy == 0 & m_trisomy > 1 & p_trisomy == 0 & (m_trisomy + null_na) < 20 & !(case_embryoid %in% diploid_tripolar)])
multiple_m_gain <- data[case_embryoid %in% karyotype_wide & m_monosomy == 0 & p_monosomy == 0 & m_trisomy > 1 & p_trisomy == 0 & (m_trisomy + null_na) < 20 & !(case_embryoid %in% diploid_tripolar)]$case_embryoid
data[, is_mm_gain := case_embryoid %in% multiple_m_gain]
table(data[is_mm_gain == T]$sample_type)
table(data[is_mm_gain == T]$sample_type) / table(data$sample_type)

# multiple paternal gain
nrow(data[case_embryoid %in% karyotype_wide & m_monosomy == 0 & p_monosomy == 0 & m_trisomy == 0 & p_trisomy > 1 & (p_trisomy + null_na) < 20 & !(case_embryoid %in% diploid_tripolar)])
multiple_p_gain <- data[case_embryoid %in% karyotype_wide & m_monosomy == 0 & p_monosomy == 0 & m_trisomy == 0 & p_trisomy > 1 & (p_trisomy + null_na) < 20 & !(case_embryoid %in% diploid_tripolar)]$case_embryoid
data[, is_mp_gain := case_embryoid %in% multiple_p_gain]
table(data[is_mp_gain == T]$sample_type)
table(data[is_mp_gain == T]$sample_type) / table(data$sample_type)

# maternal and paternal gains
nrow(data[case_embryoid %in% karyotype_wide & m_trisomy > 0 & p_trisomy > 0 & m_monosomy == 0 & p_monosomy == 0 & (p_trisomy + null_na) < 20 & (m_trisomy + null_na) < 20 & !(case_embryoid %in% diploid_tripolar)])
multiple_both_gain <- data[case_embryoid %in% karyotype_wide & m_trisomy > 0 & p_trisomy > 0 & m_monosomy == 0 & p_monosomy == 0 & (p_trisomy + null_na) < 20 & (m_trisomy + null_na) < 20 & !(case_embryoid %in% diploid_tripolar)]$case_embryoid
data[, is_both_gain := case_embryoid %in% multiple_both_gain]
table(data[is_both_gain == T]$sample_type)
table(data[is_both_gain == T]$sample_type) / table(data$sample_type)

# multiple maternal loss
nrow(data[case_embryoid %in% karyotype_wide & m_monosomy > 1 & p_monosomy == 0 & m_trisomy == 0 & p_trisomy == 0 & (m_monosomy + null_na) < 20 & !(case_embryoid %in% diploid_tripolar)])
multiple_m_loss <- data[case_embryoid %in% karyotype_wide & m_monosomy > 1 & p_monosomy == 0 & m_trisomy == 0 & p_trisomy == 0 & (m_monosomy + null_na) < 20 & !(case_embryoid %in% diploid_tripolar)]$case_embryoid
data[, is_mm_loss := case_embryoid %in% multiple_m_loss]
table(data[is_mm_loss == T]$sample_type)
table(data[is_mm_loss == T]$sample_type) / table(data$sample_type)

# multiple paternal loss
nrow(data[case_embryoid %in% karyotype_wide & m_monosomy == 0 & p_monosomy > 1 & m_trisomy == 0 & p_trisomy == 0 & (p_monosomy + null_na) < 20  & !(case_embryoid %in% diploid_tripolar)])
multiple_p_loss <- data[case_embryoid %in% karyotype_wide & m_monosomy == 0 & p_monosomy > 1 & m_trisomy == 0 & p_trisomy == 0 & (p_monosomy + null_na) < 20 & !(case_embryoid %in% diploid_tripolar)]$case_embryoid
data[, is_mp_loss := case_embryoid %in% multiple_p_loss]
table(data[is_mp_loss == T]$sample_type)
table(data[is_mp_loss == T]$sample_type) / table(data$sample_type)

# maternal and paternal losses
nrow(data[case_embryoid %in% karyotype_wide & m_monosomy > 0 & p_monosomy > 0 & m_trisomy == 0 & p_trisomy == 0 & (p_monosomy + null_na) < 20 & (m_monosomy + null_na) < 20 & !(case_embryoid %in% diploid_tripolar)])
multiple_both_loss <- data[case_embryoid %in% karyotype_wide & m_monosomy > 0 & p_monosomy > 0 & m_trisomy == 0 & p_trisomy == 0 & (p_monosomy + null_na) < 20 &  (m_monosomy + null_na) < 20 & !(case_embryoid %in% diploid_tripolar)]$case_embryoid
data[, is_both_loss := case_embryoid %in% multiple_both_loss]
table(data[is_both_loss == T]$sample_type)
table(data[is_both_loss == T]$sample_type) / table(data$sample_type)

# digynic triploid tripolar
nrow(data[case_embryoid %in% karyotype_wide & disomy >= 3 & m_trisomy >= 3 & p_monosomy > 0 & p_monosomy <= 8 & m_upd > 0 & m_upd <= 8])
digynic_triploid_tripolar <- data[case_embryoid %in% karyotype_wide & disomy >= 3 & m_trisomy >= 3 & p_monosomy > 0 & p_monosomy <= 8 & m_upd > 0 & m_upd <= 8]$case_embryoid
data[, is_digynic_triploid_tripolar := case_embryoid %in% digynic_triploid_tripolar]
table(data[is_digynic_triploid_tripolar == T]$sample_type)
table(data[is_digynic_triploid_tripolar == T]$sample_type) / table(data$sample_type)

# diandric triploid tripolar
nrow(data[case_embryoid %in% karyotype_wide & disomy >= 3 & p_trisomy >= 3 & m_monosomy > 0 & m_monosomy <= 8 & p_upd > 0 & p_upd <= 8])
diandric_triploid_tripolar <- data[case_embryoid %in% karyotype_wide & disomy >= 3 & p_trisomy >= 3 & m_monosomy > 0 & m_monosomy <= 8 & p_upd > 0 & p_upd <= 8]$case_embryoid
data[, is_diandric_triploid_tripolar := case_embryoid %in% diandric_triploid_tripolar]
table(data[is_diandric_triploid_tripolar == T]$sample_type)

########################################################################
# Plot examples of various forms of abnormalities
########################################################################

tdm_euploid <- data.table(melt(data[case_embryoid %in% euploid, c(1:2, 7:29), with = F], id.vars = c("case", "embryoid")))
tdm_euploid[, set := "Euploid"]
tdm_euploid[, case_embryoid := paste(case, embryoid, sep = "_")]
tdm_euploid <- tdm_euploid[(case_embryoid %in% sample(tdm_euploid$case_embryoid, 4))]

tdm_aneuploid <- data.table(melt(data[case_embryoid %in% aneuploid, c(1:2, 7:29), with = F], id.vars = c("case", "embryoid")))
tdm_aneuploid[, set := "Aneuploid"]
tdm_aneuploid[, case_embryoid := paste(case, embryoid, sep = "_")]
tdm_aneuploid <- tdm_aneuploid[(case_embryoid %in% sample(tdm_aneuploid$case_embryoid, 4))]

tdm_digynic_triploid <- data.table(melt(data[case_embryoid %in% digynic_triploid, c(1:2, 7:29), with = F], id.vars = c("case", "embryoid")))
tdm_digynic_triploid[, set := "Digynic Triploid"]
tdm_digynic_triploid[, case_embryoid := paste(case, embryoid, sep = "_")]
tdm_digynic_triploid <- tdm_digynic_triploid[(case_embryoid %in% sample(tdm_digynic_triploid$case_embryoid, 4))]

tdm_diandric_triploid <- data.table(melt(data[case_embryoid %in% diandric_triploid, c(1:2, 7:29), with = F], id.vars = c("case", "embryoid")))
tdm_diandric_triploid[, set := "Diandric Triploid"]
tdm_diandric_triploid[, case_embryoid := paste(case, embryoid, sep = "_")]
tdm_diandric_triploid <- tdm_diandric_triploid[(case_embryoid %in% sample(tdm_diandric_triploid$case_embryoid, 4))]

tdm_maternal_haploid <- data.table(melt(data[case_embryoid %in% maternal_haploid, c(1:2, 7:29), with = F], id.vars = c("case", "embryoid")))
tdm_maternal_haploid[, set := "Maternal Haploid"]
tdm_maternal_haploid[, case_embryoid := paste(case, embryoid, sep = "_")]
tdm_maternal_haploid <- tdm_maternal_haploid[(case_embryoid %in% sample(tdm_maternal_haploid$case_embryoid, 4))]

tdm_paternal_haploid <- data.table(melt(data[case_embryoid %in% paternal_haploid, c(1:2, 7:29), with = F], id.vars = c("case", "embryoid")))
tdm_paternal_haploid[, set := "Paternal Haploid"]
tdm_paternal_haploid[, case_embryoid := paste(case, embryoid, sep = "_")]
tdm_paternal_haploid <- tdm_paternal_haploid[(case_embryoid %in% sample(tdm_paternal_haploid$case_embryoid, 4))]

tdm_digynic_triploid_tripolar <- data.table(melt(data[case_embryoid %in% digynic_triploid_tripolar, c(1:2, 7:29), with = F], id.vars = c("case", "embryoid")))
tdm_digynic_triploid_tripolar[, set := "Digynic Triploid Tripolar"]
tdm_digynic_triploid_tripolar[, case_embryoid := paste(case, embryoid, sep = "_")]
tdm_digynic_triploid_tripolar <- tdm_digynic_triploid_tripolar[(case_embryoid %in% sample(tdm_digynic_triploid_tripolar$case_embryoid, 4))]

tdm_diandric_triploid_tripolar <- data.table(melt(data[case_embryoid %in% diandric_triploid_tripolar, c(1:2, 7:29), with = F], id.vars = c("case", "embryoid")))
tdm_diandric_triploid_tripolar[, set := "Diandric Triploid Tripolar"]
tdm_diandric_triploid_tripolar[, case_embryoid := paste(case, embryoid, sep = "_")]
tdm_diandric_triploid_tripolar <- tdm_diandric_triploid_tripolar[(case_embryoid %in% sample(tdm_diandric_triploid_tripolar$case_embryoid, 4))]

tdm_diploid_tripolar <- data.table(melt(data[case_embryoid %in% diploid_tripolar, c(1:2, 7:29), with = F], id.vars = c("case", "embryoid")))
tdm_diploid_tripolar[, set := "Diploid Tripolar"]
tdm_diploid_tripolar[, case_embryoid := paste(case, embryoid, sep = "_")]
tdm_diploid_tripolar <- tdm_diploid_tripolar[(case_embryoid %in% sample(tdm_diploid_tripolar$case_embryoid, 4))]

tdm <- rbind(tdm_euploid, tdm_aneuploid, tdm_digynic_triploid, tdm_diandric_triploid, tdm_maternal_haploid, tdm_paternal_haploid,
             tdm_digynic_triploid_tripolar, tdm_diandric_triploid_tripolar, tdm_diploid_tripolar)

tdm[, bph := as.character(NA)]
tdm[grepl("mbph", value), bph := "Mat. Meiotic Chrom. Gain"]
tdm[grepl("pbph", value), bph := "Pat. Meiotic Chrom. Gain"]
tdm[, value := gsub("_pbph", "", gsub("_mbph", "", value))]
tdm[, variable := gsub("c", "", variable)]

tdm$variable <- factor(tdm$variable, levels = c(1:22, "Sex"))

tdm$set <- factor(tdm$set, levels = c("Euploid", "Aneuploid", "Digynic Triploid", "Diandric Triploid", "Maternal Haploid", "Paternal Haploid",
                                      "Digynic Triploid Tripolar", "Diandric Triploid Tripolar", "Diploid Tripolar"))

ggplot(data = tdm, aes(x = variable, y = paste(case, embryoid, sep = "_"), fill = factor(value), color = bph)) +
  geom_tile(width = 0.9, height = 0.9, size = 1) +
  facet_grid(set ~ ., scales = "free", switch = "both") +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_blank(), legend.title = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  xlab("") +
  ylab("")
  
########################################################################
# Simulate tripolar mitosis and compare to observed chromosome counts
########################################################################

chaotic_melt <- melt(data[case_embryoid %in% diploid_tripolar][, c("case_embryoid", "m_monosomy", "p_monosomy", "nullisomy", "disomy")], 
                     id.vars = c("case_embryoid"), measurement.vars = c("m_monosomy", "p_monosomy", "nullisomy", "disomy"))

simulate_tripolar <- function() {
  n_maternal <- replicate(23, sample(c(0, 1), 1, prob = c(1/3, 2/3)))
  n_paternal <- replicate(23, sample(c(0, 1), 1, prob = c(1/3, 2/3)))
  m_monosomy <- sum(n_maternal == 0 & n_paternal == 1)
  p_monosomy <- sum(n_maternal == 1 & n_paternal == 0)
  nullisomy <- sum(n_maternal == 0 & n_paternal == 0)
  disomy <- sum(n_maternal == 1 & n_paternal == 1)
  return(data.table(m_monosomy = m_monosomy, p_monosomy = p_monosomy, nullisomy = nullisomy, disomy = disomy))
}

simulated_tripolar <- do.call(rbind, lapply(1:10000, function(x) simulate_tripolar()))
simulated_tripolar[, case_embryoid := .I]
simulated_tripolar_melt <- melt(simulated_tripolar, id.vars = c("case_embryoid"), measurement.vars = c("m_monosomy", "p_monosomy", "nullisomy", "disomy"))

chaotic_melt[, set := "Observed"]
simulated_tripolar_melt[, set := "Simulated"]

combined_data <- rbind(chaotic_melt, simulated_tripolar_melt)
combined_data$set <- factor(combined_data$set, levels=c("Simulated", "Observed"))

ggplot(data = combined_data, aes(x = value, y = variable, fill = variable)) +
  geom_joy(stat = "binline", binwidth = 1, scale = 0.9) +
  facet_grid(set ~ .) +
  xlim(-1, 23) +
  theme_joy() +
  theme(legend.position = "none")

# replicated genome
initial_chroms <- 46 * 2

sim_counts <- rmultinom(10000, initial_chroms, rep(1/3, 3))[1,]
sim_counts <- data.table(count = sim_counts, set = "Simulated")

obs_counts <- data[case_embryoid %in% diploid_tripolar]$chrom_count
obs_counts <- data.table(count = obs_counts, set = "Observed")

combined_counts <- rbind(sim_counts, obs_counts)
combined_counts$set <- factor(combined_counts$set, levels=c("Simulated", "Observed"))
  
ggplot(data = combined_counts, aes(x = count, y = set, fill = set)) +
  geom_joy(stat = "binline", binwidth = 2, scale = 0.9) +
  xlim(-1, 50) +
  theme_joy() +
  theme(legend.position = "none") +
  facet_grid(set ~ .) +
  scale_fill_brewer(palette = "Set2")
  
########################################
# Analysis of time-lapse phenotypes
########################################

# import time-lapse data
duc <- fread("https://raw.githubusercontent.com/rmccoy7541/tripolar_mitosis/master/data/duc.csv") %>%
  setnames(., c("DUC", "case"))
duc <- duc[case != "missing"]
duc[, case := as.integer(case)]

# fix typos, "DU4" -> "DUC"
duc[DUC == "DU4-1", DUC := "DUC-1"]
duc[DUC == "DU4-1Plus", DUC := "DUC-1Plus"]
table(duc$DUC)

# classify as DUC or non-DUC
duc[, is_DUC := as.logical(FALSE)]
duc[DUC %in% c("DUC-1", "DUC-2", "DUC-3", "DUC-1Plus", "DUC-2Plus", "DUC-3Plus"), is_DUC := TRUE]

# summarize by case
duc_summary <- data.table(group_by(duc, case) %>% 
                            summarise(., duc_t = sum(is_DUC == T), duc_f = sum(is_DUC == F),
                                      duc1 = sum(DUC == "DUC-1"), duc1_plus = sum(DUC == "DUC-1Plus"),
                                      duc2 = sum(DUC == "DUC-2"), duc2_plus = sum(DUC == "DUC-2Plus"),
                                      duc3 = sum(DUC == "DUC-3"), duc3_plus = sum(DUC == "DUC-3Plus")))

# import maternal genotype data
geno_discovery <- fread("~/tripolar_mitosis/rs2305957_f_all.ped")
geno_validation <- fread("~/tripolar_mitosis/rs2305957_f_all_validation.ped")
geno <- rbind(geno_discovery, geno_validation)
setnames(geno, "V1", "case")
geno[, gt := as.numeric(NA)]
geno[V7 == "G" & V8 == "G", gt := 0]
geno[V7 == "A" & V8 == "G", gt := 1]
geno[V7 == "A" & V8 == "A", gt := 2]

duc_summary <- merge(duc_summary, geno, "case")

# import king-inferred repeat cases
king <- fread("https://raw.githubusercontent.com/rmccoy7541/tripolar_mitosis/master/data/king.txt", header = F) %>%
  setnames(., c("case", "king"))

duc_summary <- merge(duc_summary, king, "case")

# sum DUC and non-DUC for cases inferred as same individual
duc_summary <- data.table(group_by(duc_summary, king) %>% 
                            summarise(., s_duc_t = sum(duc_t), s_duc_f = sum(duc_f), gt = unique(gt)))

summary(glm(data = duc_summary, formula = cbind(s_duc_t, s_duc_f) ~ gt, family = quasibinomial))

# plot results
duc_summary[, genotype := as.factor(paste(V7, V8, sep = ""))]
duc_summary$genotype <- factor(duc_summary$genotype, levels = c("GG", "AG", "AA"))

ggplot(data = duc_summary, aes(x = gt, y = s_duc_t / (s_duc_t + s_duc_f))) +
  geom_dotplot(aes(fill = as.factor(gt)), binwidth = 0.01, binaxis = "y", stackdir = "center", dotsize = 2) +
  theme_bw() +
  xlab("Genotype") +
  ylab("Prop. DUC") +
  theme(panel.grid = element_blank(), legend.position = "none") +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2")

