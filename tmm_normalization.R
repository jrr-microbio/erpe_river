library(edgeR)

setwd('/Volumes/GoogleDrive/My Drive/University/Wrighton Lab; PhD/ERPE_Timeseries_viral_communities/microbial_work/read_mapping/95ID_noambig_filtered/')

#MAGs as rows. Samples as columns. Columns (samples) need to be in the same exact order that the library sizes file samples are in.
counts_dat = read.csv('/Volumes/GoogleDrive/My Drive/University/Wrighton Lab; PhD/ERPE_Timeseries_viral_communities/microbial_work/read_mapping/95ID_noambig_filtered/coverM_output/coverM_counts_125_MAGs_75AF_3xcov_forTMM_noambig_remadefile.csv', row.names=1, header=T)
counts_dat=data.matrix(counts_dat)

#Two columns. Samples in rows and the library size as the 2nd column. Needs to be in exact same order as abunds. No headers. No nothing. Just two columns of values.
library_sizes = read.table('/Volumes/GoogleDrive/My Drive/University/Wrighton Lab; PhD/ERPE_Timeseries_viral_communities/microbial_work/read_mapping/library_sizes.tsv', row.names=1, col.names=c('sample_name', 'library_size'))

lib.size= c(library_sizes['library_size'])$library_size

counts_norm <- DGEList(counts=counts_dat, lib.size= c(library_sizes['library_size'])$library_size)

counts_norm <- calcNormFactors(counts_norm)

tmms <- cpm(counts_norm)
write.table(tmms, "120_MAGs_75AF_3xcov_counts_TMM_normalized_noambig_methodsnone.tsv", sep='\t', quote=F)




########Doing this separately because it might spit out errors depending on how sparse the datasets are:

#surface

counts_dat = read.table('/Volumes/GoogleDrive/My Drive/University/Wrighton Lab; PhD/ERPE_Timeseries_viral_communities/microbial_work/read_mapping/95ID_noambig_filtered/coverM_counts_125_MAGs_75AF_3xcov_forTMM_noambig_sw.txt', row.names=1, header=T, sep="\t")

#Two columns. Samples in rows and the library size as the 2nd column. Needs to be in exact same order as abunds. No headers. No nothing. Just two columns of values.

library_sizes = read.table('/Volumes/GoogleDrive/My Drive/University/Wrighton Lab; PhD/ERPE_Timeseries_viral_communities/microbial_work/read_mapping/library_sizes_sw.tsv', row.names=1, col.names=c('sample_name', 'library_size'))

lib.size= c(library_sizes['library_size'])$library_size

counts_norm <- DGEList(counts=counts_dat, lib.size= c(library_sizes['library_size'])$library_size)

counts_norm <- calcNormFactors(counts_norm)

tmms <- cpm(counts_norm)
write.table(tmms, "125_MAGs_75AF_3xcov_counts_TMM_normalized_noambig_SW.tsv", sep='\t', quote=F)


#pore

counts_dat = read.table('/Volumes/GoogleDrive/My Drive/University/Wrighton Lab; PhD/ERPE_Timeseries_viral_communities/microbial_work/read_mapping/95ID_noambig_filtered/coverM_counts_125_MAGs_75AF_3xcov_forTMM_noambig_pw.txt', row.names=1, header=T, sep="\t")

#Two columns. Samples in rows and the library size as the 2nd column. Needs to be in exact same order as abunds. No headers. No nothing. Just two columns of values.

library_sizes = read.table('/Volumes/GoogleDrive/My Drive/University/Wrighton Lab; PhD/ERPE_Timeseries_viral_communities/microbial_work/read_mapping/library_sizes_pw.tsv', row.names=1, col.names=c('sample_name', 'library_size'))

lib.size= c(library_sizes['library_size'])$library_size

counts_norm <- DGEList(counts=counts_dat, lib.size= c(library_sizes['library_size'])$library_size)

counts_norm <- calcNormFactors(counts_norm, method= "TMM")

tmms <- cpm(counts_norm)
write.table(tmms, "125_MAGs_75AF_3xcov_counts_TMM_normalized_noambig_PW.tsv", sep='\t', quote=F)

