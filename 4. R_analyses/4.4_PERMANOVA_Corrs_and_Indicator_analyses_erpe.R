# load packages
library(vegan)
library(tidyverse)
library(mctoolsr)
library(GGally)
library(ecodist)
library(matrixStats)
library(viridis)

###############################################################################
# read in data and format
###############################################################################

#Loading in our data
dat=read.table("4.1.1_32_sample_geochemistry_incl_FTICR_subset.txt", header = T, sep="\t")
rownames(dat)<-dat[,1]
dat=dat[-1]

metadat=dat[,1:3]
log_dat=log(dat[,4:29])



map_pw = dat %>%
  filter(sw_pw == 'pw')

map_sw = dat %>%
  filter(sw_pw == 'sw')

mags_tab = read.table('4.3.1_MAGs_3x_depth_75AF_TMM_normalized.tsv', sep="\t", header=T,check.names=TRUE, row.names = 1) %>%
#  t() %>%
  as.data.frame()
  
vmags_tab = read.table('4.2.1_1230_vMAGs_3x_depth_75AF_TMM_normalized_32_sample.tsv', sep="\t", header=T,check.names=TRUE, row.names = 1) %>%
#  t() %>%
  as.data.frame() %>%
  select(colnames(mags_tab))

mags_tab_sw = mags_tab[1:15]
mags_tab_pw = mags_tab[16:32]

vmags_tab_sw = vmags_tab[1:15]
vmags_tab_pw = vmags_tab[16:32]

# PERMANOVA - tmm MAGs
dm.mags = calc_dm(mags_tab, method = "bray")
adonis(formula = dm.mags ~ sw_pw + timepoint, data = map, permutations = 999)

# PERMANOVA - tmm vMAGs
dm.vmags = calc_dm(vmags_tab, method = "bray")
adonis(formula = dm.vmags ~ sw_pw + timepoint, data = map, permutations = 999)

# PERMANOVA - tmm PCA
dm.chem = calc_dm(t(dat[,4:29]), method = "bray")
adonis(formula = dm.chem ~ sw_pw + timepoint, data = metadat, permutations = 999)

###############################################################################
# Next up - indicator - surf vs. pw
###############################################################################

library(indicspecies)

# TMM
set.seed(1234)
indic.swpw.getmm <- multipatt(t(vmags_tab), map$sw_pw, duleg = F, func = "r.g",
                                 permutations = 9999)

indic.swpw.getmm <- indic.swpw.getmm$sign
indic.swpw.getmm$MAG_id <- rownames(indic.swpw.getmm)
indic.swpw.getmm$p.fdr <- p.adjust(indic.swpw.getmm$p.value, method = "fdr")
indic.swpw.getmm <- indic.swpw.getmm %>%
  filter(p.fdr < 0.05)
table(indic.swpw.getmm$index)


# write out results
write.csv(indic.swpw.getmm, "indicator_out_vMAG.csv")


### AND for MAGS

# TMM
set.seed(1234)
indic.swpw.getmm <- multipatt(t(mags_tab), map$sw_pw, duleg = F, func = "r.g",
                              permutations = 9999)

indic.swpw.getmm <- indic.swpw.getmm$sign
indic.swpw.getmm$MAG_id <- rownames(indic.swpw.getmm)
indic.swpw.getmm$p.fdr <- p.adjust(indic.swpw.getmm$p.value, method = "fdr")
indic.swpw.getmm <- indic.swpw.getmm %>%
  filter(p.fdr < 0.05)

table(indic.swpw.getmm$index)

# write out results
write.csv(indic.swpw.getmm, "indicator_out_MAG")
