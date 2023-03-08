library(edgeR)

setwd('/Volumes/Macintosh HD/Users/josue.rodriguez/Library/CloudStorage/GoogleDrive-jarora2213@gmail.com/My Drive/University/wrighton_lab_phd/erpe_timeseries_viral_communities/_paper_figures/NMDS_and_PCA/NMDS_plots/MAG_NMDS/')

#Genomes as rows. Samples as columns. Columns (samples) need to be in the same exact order that the library sizes file samples are in.
counts_dat = read.csv('3.1_MAGs_counts_forTMM.csv', row.names=1, header=T)
counts_dat_vir = read.csv('3.2_vMAGs_counts_forTMM.csv', row.names=1, header=T)

counts_dat=data.matrix(counts_dat)
counts_dat_vir=data.matrix(counts_dat_vir)

#For file library sizes, it needs two columns. Samples in rows and the library size as the 2nd column. Needs to be in exact same order as abunds. No headers. No nothing. Just two columns of values.
library_sizes = read.table('/Volumes/Macintosh HD/Users/josue.rodriguez/Library/CloudStorage/GoogleDrive-jarora2213@gmail.com/My Drive/University/wrighton_lab_phd/erpe_timeseries_viral_communities/_paper_figures/NMDS_and_PCA/NMDS_plots/MAG_NMDS/MAG_coverM_output/library_sizes.tsv', row.names=1, col.names=c('sample_name', 'library_size'))

lib.size= c(library_sizes['library_size'])$library_size

counts_norm <- DGEList(counts=counts_dat, lib.size= c(library_sizes['library_size'])$library_size)
counts_norm_vir <- DGEList(counts=counts_dat_vir, lib.size= c(library_sizes['library_size'])$library_size)

counts_norm <- calcNormFactors(counts_norm)
counts_norm_vir <- calcNormFactors(counts_norm_vir)

tmms <- cpm(counts_norm)
tmms_vir <- cpm(counts_norm_vir)

write.table(tmms, "125_MAGs_3x_depth_75AF_TMM_normalized.tsv", sep='\t', quote=F)
write.table(tmms_vir, "1230_vMAGs_3x_depth_75AF_TMM_normalized.tsv", sep='\t', quote=F)
