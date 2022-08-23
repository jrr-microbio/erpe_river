library(tidyverse)
library(vegan)
library(reshape2)

setwd("/Volumes/GoogleDrive/My Drive/University/Wrighton Lab; PhD/ERPE_Timeseries_viral_communities/microbial_work/read_mapping/95ID_noambig_filtered/clr_abunds/3samp_cutoff/")


###################################################
####################################
#Since we are splitting up the compartments for rCLR, I will first start with the SW
####################################
###################################################


#Read in the data. Rows need to be OTUs. Columns need to be samples.
tmm_sw=read.csv("/Volumes/GoogleDrive/My Drive/University/Wrighton Lab; PhD/ERPE_Timeseries_viral_communities/microbial_work/read_mapping/95ID_noambig_filtered/clr_abunds/3samp_cutoff/69_MAGs_surfacewater_erpe_3samp.csv", header=T, row.names =1)
#Make the data into a matrix that is transposed since above read in was not correct.

t_tmm_sw= as.data.frame(t(tmm_sw)) #To have a df for vegdist
t_tmm_sw_mx = as.matrix(t_tmm_sw) #To have the matrix for melt.

#Calculated the euclidean distance of non-rclr data just to compare.
noclrdist=vegdist(t_tmm_sw, method="euclidean") #non clr'd euclidean distance. Just for funsies.

#Now working on the CLR. According to pat's video (https://www.youtube.com/watch?v=ulo7WatBEAo), this needs to be in a long format. So going to use reshape for this one and a simple melt. Column headers will then be: Var1, Var2, and value where var1=genome id, var2= sample id, value = rel abund. He uses a pivot long/wide but im more used to reshape2 library.

t_tmm_sw_long = melt(t_tmm_sw_mx) %>%
  as.data.frame() #The melt needs a matrix to work, so I start it off as a matrix and then transform it into a data frame because the rCLR function below needs it that way.

#First, just calculate the total abundances for each sample in the dataframe read in above. We'll need this to compare the strategies down the line.
group_count = t_tmm_sw_long %>%
  group_by(Var1) %>%
  summarize(totalabunds=sum(value))

#Now we make the geometric mean function.
gm=function(x){
  exp(mean(log(x[x>0])))
}

gm(c(2,3,4, 0)) #this is working. Its doing "robust" clr where it is only calculating based on numbers that are greater than 0.


#Now we make the full table that contains the rCLR values.
rclr_sw = t_tmm_sw_long %>% 
  group_by(Var1) %>%
  mutate(rclr=log(value/gm(value))) %>%
  ungroup() %>%
  select(-value) %>%
  pivot_wider(names_from = Var1, values_from=rclr, values_fill = 0) %>%
  as.data.frame()

#Cleaning up the dataframe a little.
rownames(rclr_sw) = rclr_sw$Var2
rclr_sw=rclr_sw[,-1]

#My function is still computing the values for 0's and is making all of them "-Inf" which is the behavior i would expect. 

log(c(2,3,4, 0)/gm(c(2,3,4, 0)))#testing original function and can confirm it is still feeding that -inf and that its not my dataset or manipulation.

#I manually checked others in my df, and it it correctly calculating the rCLR values and ignoring those zeroes, so now i will simply just remove the "-Inf" values and put them back as 0.

rclr_sw[rclr_sw=="-Inf"]=0

#Then we just remake a data frame and transpose it so it matches the original non-clrd vegdist that I will calculate below.
rclr_sw=as.data.frame(rclr_sw)
rclr_sw_t=as.data.frame(t(rclr_sw))

#clr'd euclidean distance.
clrdist=vegdist(rclr_sw_t, method="euclidean") 

#Now i can plot and compare the clr and non-clr distances.
noclrdist_tbl = noclrdist %>%
  as.matrix %>%
  as_tibble(rownames="sample") %>%
  pivot_longer(cols=-sample) %>%
  filter(name < sample)

clrdist_tbl = clrdist %>%
  as.matrix %>%
  as_tibble(rownames="sample") %>%
  pivot_longer(cols=-sample) %>%
  filter(name < sample)

inner_join(noclrdist_tbl, clrdist_tbl, by=c("sample", "name")) %>%
  select(sample, name, noclr=value.x, rclr=value.y) %>%
  inner_join(., group_count, by=c("sample" = "Var1")) %>%
  inner_join(., group_count, by=c("name" = "Var1")) %>%
  mutate(diffs=abs(totalabunds.x - totalabunds.y)) %>%
  pivot_longer(cols=c(noclr, rclr), names_to = "method", values_to="dist") %>%
  ggplot(aes(x=diffs,y=dist)) +
  geom_point()+
  facet_wrap(~method, scales="free_y")+
  geom_smooth()

write.table(rclr_sw_t, '69_MAGs_75AF_3xcov_MAG_rclr_TMM_SW.tsv', sep='\t', col.names=NA)












###################################################
####################################
#Now i do the same as above but for the PW.
####################################
###################################################







#Read in the data. Rows need to be OTUs. Columns need to be samples.
tmm_pw=read.csv("/Volumes/GoogleDrive/My Drive/University/Wrighton Lab; PhD/ERPE_Timeseries_viral_communities/microbial_work/read_mapping/95ID_noambig_filtered/clr_abunds/3samp_cutoff/38_MAGs_porewater_erpe_3samp.csv", header=T, row.names =1)

#Make the data into a matrix that is transposed since above read in was not correct.
t_tmm_pw= as.data.frame(t(tmm_pw)) #To have a df for vegdist
t_tmm_pw_mx = as.matrix(t_tmm_pw) #To have the matrix for melt.

#Calculated the euclidean distance of non-rclr data just to compare.
noclrdist_pw=vegdist(t_tmm_pw, method="euclidean") #non clr'd euclidean distance. Just for funsies.

#Now working on the CLR. According to pat's video (https://www.youtube.com/watch?v=ulo7WatBEAo), this needs to be in a long format. So going to use reshape for this one and a simple melt. Column headers will then be: Var1, Var2, and value where var1=genome id, var2= sample id, value = rel abund. He uses a pivot long/wide but im more used to reshape2 library.

t_tmm_pw_long = melt(t_tmm_pw_mx) %>%
  as.data.frame() #The melt needs a matrix to work, so I start it off as a matrix and then transform it into a data frame because the rCLR function below needs it that way.

#First, just calculate the total abundances for each sample in the dataframe read in above. We'll need this to compare the strategies down the line.
group_count = t_tmm_pw_long %>%
  group_by(Var1) %>%
  summarize(totalabunds=sum(value))

#Now we make the geometric mean function.
gm=function(x){
  exp(mean(log(x[x>0])))
}

gm(c(2,3,4, 0)) #this is working. Its doing "robust" clr where it is only calculating based on numbers that are greater than 0.


#Now we make the full table that contains the rCLR values.
rclr_pw = t_tmm_pw_long %>% 
  group_by(Var1) %>%
  mutate(rclr=log(value/gm(value))) %>%
  ungroup() %>%
  select(-value) %>%
  pivot_wider(names_from = Var1, values_from=rclr, values_fill = 0) %>%
  as.data.frame()

#Cleaning up the dataframe a little.
rownames(rclr_pw) = rclr_pw$Var2
rclr_pw=rclr_pw[,-1]

#My function is still computing the values for 0's and is making all of them "-Inf" which is the behavior i would expect. 

log(c(2,3,4, 0)/gm(c(2,3,4, 0)))#testing original function and can confirm it is still feeding that -inf and that its not my dataset or manipulation.

#I manually checked others in my df, and it it correctly calculating the rCLR values and ignoring those zeroes, so now i will simply just remove the "-Inf" values and put them back as 0.

rclr_pw[rclr_pw=="-Inf"]=0

#Then we just remake a data frame and transpose it so it matches the original non-clrd vegdist that I will calculate below.
rclr_pw=as.data.frame(rclr_pw)
rclr_pw_t=as.data.frame(t(rclr_pw))

#clr'd euclidean distance.
clrdist_pw=vegdist(rclr_pw_t, method="euclidean") 

#Now i can plot and compare the clr and non-clr distances.
noclrdist_tbl_pw = noclrdist_pw %>%
  as.matrix %>%
  as_tibble(rownames="sample") %>%
  pivot_longer(cols=-sample) %>%
  filter(name < sample)

clrdist_tbl_pw = clrdist_pw %>%
  as.matrix %>%
  as_tibble(rownames="sample") %>%
  pivot_longer(cols=-sample) %>%
  filter(name < sample)

inner_join(noclrdist_tbl_pw, clrdist_tbl_pw, by=c("sample", "name")) %>%
  select(sample, name, noclr=value.x, rclr=value.y) %>%
  inner_join(., group_count, by=c("sample" = "Var1")) %>%
  inner_join(., group_count, by=c("name" = "Var1")) %>%
  mutate(diffs=abs(totalabunds.x - totalabunds.y)) %>%
  pivot_longer(cols=c(noclr, rclr), names_to = "method", values_to="dist") %>%
  ggplot(aes(x=diffs,y=dist)) +
  geom_point()+
  facet_wrap(~method, scales="free_y")+
  geom_smooth()

write.table(rclr_pw_t, '38_MAGs_75AF_3xcov_MAG_rclr_TMM_pw.tsv', sep='\t', col.names=NA)








###################################################
####################################
#Now i do the same as above but for the "neither" category.
####################################
###################################################







#Read in the data. Rows need to be OTUs. Columns need to be samples.
tmm_neither=read.csv("/Volumes/GoogleDrive/My Drive/University/Wrighton Lab; PhD/ERPE_Timeseries_viral_communities/viral_work/read_mapping/95ID_noambig_filtered/clr_abunds/58_vMAGs_75AF_3xcov_counts_TMM_normalized_noambig_neither.indic_in3ormoresamples.csv", header=T, row.names =1) #note, this is all samples EXCEPT SW 33hr which had 0 counts within it. So it was removed. Further, if using the "must be in at least 3 samples to be considered "neither", then we remove trimmed_erpe_2018_pw_42hr_WHONDRS.PP48.000099_95ID.sorted as it becomes empty.

#Make the data into a matrix that is transposed since above read in was not correct.
t_tmm_neither= as.data.frame(t(tmm_neither)) #To have a df for vegdist
t_tmm_neither_mx = as.matrix(t_tmm_neither) #To have the matrix for melt.

#Calculated the euclidean distance of non-rclr data just to compare.
noclrdist_neither=vegdist(t_tmm_neither, method="euclidean") #non clr'd euclidean distance. Just for funsies.

#Now working on the CLR. According to pat's video (https://www.youtube.com/watch?v=ulo7WatBEAo), this needs to be in a long format. So going to use reshape for this one and a simple melt. Column headers will then be: Var1, Var2, and value where var1=genome id, var2= sample id, value = rel abund. He uses a pivot long/wide but im more used to reshape2 library.

t_tmm_neither_long = melt(t_tmm_neither_mx) %>%
  as.data.frame() #The melt needs a matrix to work, so I start it off as a matrix and then transform it into a data frame because the rCLR function below needs it that way.

#First, just calculate the total abundances for each sample in the dataframe read in above. We'll need this to compare the strategies down the line.
group_count = t_tmm_neither_long %>%
  group_by(Var1) %>%
  summarize(totalabunds=sum(value))

#Now we make the geometric mean function.
gm=function(x){
  exp(mean(log(x[x>0])))
}

gm(c(2,3,4, 0)) #this is working. Its doing "robust" clr where it is only calculating based on numbers that are greater than 0.


#Now we make the full table that contains the rCLR values.
rclr_neither = t_tmm_neither_long %>% 
  group_by(Var1) %>%
  mutate(rclr=log(value/gm(value))) %>%
  ungroup() %>%
  select(-value) %>%
  pivot_wider(names_from = Var1, values_from=rclr, values_fill = 0) %>%
  as.data.frame()

#Cleaning up the dataframe a little.
rownames(rclr_neither) = rclr_neither$Var2
rclr_neither=rclr_neither[,-1]

#My function is still computing the values for 0's and is making all of them "-Inf" which is the behavior i would expect. 

log(c(2,3,4, 0)/gm(c(2,3,4, 0)))#testing original function and can confirm it is still feeding that -inf and that its not my dataset or manipulation.

#I manually checked others in my df, and it it correctly calculating the rCLR values and ignoring those zeroes, so now i will simply just remove the "-Inf" values and put them back as 0.

rclr_neither[rclr_neither=="-Inf"]=0

#Then we just remake a data frame and transpose it so it matches the original non-clrd vegdist that I will calculate below.
rclr_neither=as.data.frame(rclr_neither)
rclr_neither_t=as.data.frame(t(rclr_neither))

#clr'd euclidean distance.
clrdist_neither=vegdist(rclr_neither_t, method="euclidean") 

#Now i can plot and compare the clr and non-clr distances.
noclrdist_tbl_neither = noclrdist_neither %>%
  as.matrix %>%
  as_tibble(rownames="sample") %>%
  pivot_longer(cols=-sample) %>%
  filter(name < sample)

clrdist_tbl_neither = clrdist_neither %>%
  as.matrix %>%
  as_tibble(rownames="sample") %>%
  pivot_longer(cols=-sample) %>%
  filter(name < sample)

inner_join(noclrdist_tbl_neither, clrdist_tbl_neither, by=c("sample", "name")) %>%
  select(sample, name, noclr=value.x, rclr=value.y) %>%
  inner_join(., group_count, by=c("sample" = "Var1")) %>%
  inner_join(., group_count, by=c("name" = "Var1")) %>%
  mutate(diffs=abs(totalabunds.x - totalabunds.y)) %>%
  pivot_longer(cols=c(noclr, rclr), names_to = "method", values_to="dist") %>%
  ggplot(aes(x=diffs,y=dist)) +
  geom_point()+
  facet_wrap(~method, scales="free_y")+
  geom_smooth()

write.table(rclr_neither_t, '171_vMAGs_75AF_3xcov_vMAG_rclr_TMM_neither.tsv', sep='\t', col.names=NA)



#clr_getmms_sw = clr(getmms_pseudo_sw)
#This was a padckage to do clr from library compositions, but i rather do robust clr and calculate manually so i can see what is actually happening.






































#####################################################################################################################################################
#######This is some extra code that I tested. TMM vs non-TMM across distance metrics to gauge lib size impact on distances. TMM-rCLR vs TMM-noCLR to compare the impact of CLR on euclidean distances.




###################################################
####################################
#Try and do the exact same thing as above but with NON TMM normalized counts data to see how these change with/without TMM.
####################################
###################################################







tmm_swcoverm=read.csv("/Volumes/GoogleDrive/My Drive/University/Wrighton Lab; PhD/ERPE_Timeseries_viral_communities/viral_work/read_mapping/95ID_noambig_filtered/clr_abunds/testing_nonTMM/coverM_rawcounts_982_vMAGs_75AF_3xcov_forTMM_noambig_sw.indic.csv", header=T, row.names =1)
#Make the data into a matrix that is transposed since above read in was not correct.

t_tmm_swcoverm= as.data.frame(t(tmm_swcoverm)) #To have a df for vegdist
t_tmm_swcoverm_mx = as.matrix(t_tmm_swcoverm) #To have the matrix for melt.

#Calculated the euclidean distance of non-rclr data just to compare.
noclrdist_coverm=vegdist(t_tmm_swcoverm, method="euclidean") #non clr'd euclidean distance. Just for funsies.

#Now working on the CLR. According to pat's video (https://www.youtube.com/watch?v=ulo7WatBEAo), this needs to be in a long format. So going to use reshape for this one and a simple melt. Column headers will then be: Var1, Var2, and value where var1=genome id, var2= sample id, value = rel abund. He uses a pivot long/wide but im more used to reshape2 library.

t_tmm_swcoverm_long = melt(t_tmm_swcoverm_mx) %>%
  as.data.frame() #The melt needs a matrix to work, so I start it off as a matrix and then transform it into a data frame because the rCLR function below needs it that way.

#First, just calculate the total abundances for each sample in the dataframe read in above. We'll need this to compare the strategies down the line.
group_count_coverM = t_tmm_swcoverm_long %>%
  group_by(Var1) %>%
  summarize(totalabunds=sum(value))

#Now we make the geometric mean function.
gm=function(x){
  exp(mean(log(x[x>0])))
}

gm(c(2,3,4, 0)) #this is working. Its doing "robust" clr where it is only calculating based on numbers that are greater than 0.


#Now we make the full table that contains the rCLR values.
rclr_swcoverm = t_tmm_swcoverm_long %>% 
  group_by(Var1) %>%
  mutate(rclr=log(value/gm(value))) %>%
  ungroup() %>%
  select(-value) %>%
  pivot_wider(names_from = Var1, values_from=rclr, values_fill = 0) %>%
  as.data.frame()

#Cleaning up the dataframe a little.
rownames(rclr_swcoverm) = rclr_swcoverm$Var2
rclr_swcoverm=rclr_swcoverm[,-1]

rclr_swcoverm[rclr_swcoverm=="-Inf"]=0

#Then we just remake a data frame and transpose it so it matches the original non-clrd vegdist that I will calculate below.
rclr_swcoverm=as.data.frame(rclr_swcoverm)
rclr_swcoverm_t=as.data.frame(t(rclr_swcoverm))

#clr'd euclidean distance.
clrdist_coverm=vegdist(rclr_swcoverm_t, method="euclidean") 

#Now i can plot and compare the clr and non-clr distances.
noclrdist_coverm_tbl = noclrdist_coverm %>%
  as.matrix %>%
  as_tibble(rownames="sample") %>%
  pivot_longer(cols=-sample) %>%
  filter(name < sample)

clrdist_coverm_tbl = clrdist_coverm %>%
  as.matrix %>%
  as_tibble(rownames="sample") %>%
  pivot_longer(cols=-sample) %>%
  filter(name < sample)

inner_join(noclrdist_coverm_tbl, clrdist_coverm_tbl, by=c("sample", "name")) %>%
  select(sample, name, noclr=value.x, rclr=value.y) %>%
  inner_join(., group_count_coverM, by=c("sample" = "Var1")) %>%
  inner_join(., group_count_coverM, by=c("name" = "Var1")) %>%
  mutate(diffs=abs(totalabunds.x - totalabunds.y)) %>%
  pivot_longer(cols=c(noclr, rclr), names_to = "method", values_to="dist") %>%
  ggplot(aes(x=diffs,y=dist)) +
  geom_point()+
  facet_wrap(~method, scales="free_y")+
  geom_smooth()


#It's pretty clear that we need to TMM normalize our data, but that the rCLR calculation is a little less important after TMM normalization. Plot legend: I am comparing the Euclidean distance between genomes of the SW to the overall difference in number of either A) counts (no TMM) or B) total relative abundance (with TMM). A straight-ish line that is going straight on the x axis indicates that the distances between the genomes are more equal across samples, and as such are more normalized across independent samples and helps reinforce the I.I.D (independent, and identically distributed) clause for statistical analyses. We see that without TMM, this identically distributed is absolute crap, and that it would violate, but that with TMM +/- CLR we get a much more even distribution that will not violate these assumptions. (edited) 





#############################################
########################non-TMM but for pw same as above




tmm_pwcoverm=read.csv("/Volumes/GoogleDrive/My Drive/University/Wrighton Lab; PhD/ERPE_Timeseries_viral_communities/viral_work/read_mapping/95ID_noambig_filtered/clr_abunds/testing_nonTMM/coverM_counts_97_vMAGs_75AF_3xcov_forTMM_noambig_pw.indic.csv", header=T, row.names =1)
#Make the data into a matrix that is transposed since above read in was not correct.

t_tmm_pwcoverm= as.data.frame(t(tmm_pwcoverm)) #To have a df for vegdist
t_tmm_pwcoverm_mx = as.matrix(t_tmm_pwcoverm) #To have the matrix for melt.

#Calculated the euclidean distance of non-rclr data just to compare.
noclrdist_covermpw=vegdist(t_tmm_pwcoverm, method="euclidean") #non clr'd euclidean distance. Just for funsies.

#Now working on the CLR. According to pat's video (https://www.youtube.com/watch?v=ulo7WatBEAo), this needs to be in a long format. So going to use reshape for this one and a simple melt. Column headers will then be: Var1, Var2, and value where var1=genome id, var2= sample id, value = rel abund. He uses a pivot long/wide but im more used to reshape2 library.

t_tmm_pwcoverm_long = melt(t_tmm_pwcoverm_mx) %>%
  as.data.frame() #The melt needs a matrix to work, so I start it off as a matrix and then transform it into a data frame because the rCLR function below needs it that way.

#First, just calculate the total abundances for each sample in the dataframe read in above. We'll need this to compare the strategies down the line.
group_count_coverM = t_tmm_pwcoverm_long %>%
  group_by(Var1) %>%
  summarize(totalabunds=sum(value))

#Now we make the geometric mean function.
gm=function(x){
  exp(mean(log(x[x>0])))
}

gm(c(2,3,4, 0)) #this is working. Its doing "robust" clr where it is only calculating based on numbers that are greater than 0.


#Now we make the full table that contains the rCLR values.
rclr_pwcoverm = t_tmm_pwcoverm_long %>% 
  group_by(Var1) %>%
  mutate(rclr=log(value/gm(value))) %>%
  ungroup() %>%
  select(-value) %>%
  pivot_wider(names_from = Var1, values_from=rclr, values_fill = 0) %>%
  as.data.frame()

#Cleaning up the dataframe a little.
rownames(rclr_pwcoverm) = rclr_pwcoverm$Var2
rclr_pwcoverm=rclr_pwcoverm[,-1]

rclr_pwcoverm[rclr_pwcoverm=="-Inf"]=0

#Then we just remake a data frame and transpose it so it matches the original non-clrd vegdist that I will calculate below.
rclr_pwcoverm=as.data.frame(rclr_pwcoverm)
rclr_pwcoverm_t=as.data.frame(t(rclr_pwcoverm))

#clr'd euclidean distance.
clrdist_coverm_pw=vegdist(rclr_pwcoverm_t, method="euclidean") 

#Now i can plot and compare the clr and non-clr distances.
noclrdist_coverm_tbl_pw = noclrdist_covermpw %>%
  as.matrix %>%
  as_tibble(rownames="sample") %>%
  pivot_longer(cols=-sample) %>%
  filter(name < sample)

clrdist_coverm_tbl_pw = clrdist_coverm_pw %>%
  as.matrix %>%
  as_tibble(rownames="sample") %>%
  pivot_longer(cols=-sample) %>%
  filter(name < sample)

inner_join(noclrdist_coverm_tbl_pw, clrdist_coverm_tbl_pw, by=c("sample", "name")) %>%
  select(sample, name, noclr=value.x, rclr=value.y) %>%
  inner_join(., group_count_coverM, by=c("sample" = "Var1")) %>%
  inner_join(., group_count_coverM, by=c("name" = "Var1")) %>%
  mutate(diffs=abs(totalabunds.x - totalabunds.y)) %>%
  pivot_longer(cols=c(noclr, rclr), names_to = "method", values_to="dist") %>%
  ggplot(aes(x=diffs,y=dist)) +
  geom_point()+
  facet_wrap(~method, scales="free_y")+
  geom_smooth()















#############################################
########################non-TMM but for neither same as above




tmm_neithercoverm=read.csv("/Volumes/GoogleDrive/My Drive/University/Wrighton Lab; PhD/ERPE_Timeseries_viral_communities/viral_work/read_mapping/95ID_noambig_filtered/clr_abunds/testing_nonTMM/coverM_counts_58_vMAGs_75AF_3xcov_forTMM_noambig_neither.indic_ge3samples.csv", header=T, row.names =1)
#Make the data into a matrix that is transposed since above read in was not correct.

t_tmm_neithercoverm= as.data.frame(t(tmm_neithercoverm)) #To have a df for vegdist
t_tmm_neithercoverm_mx = as.matrix(t_tmm_neithercoverm) #To have the matrix for melt.

#Calculated the euclidean distance of non-rclr data just to compare.
noclrdist_covermneither=vegdist(t_tmm_neithercoverm, method="euclidean") #non clr'd euclidean distance. Just for funsies.

#Now working on the CLR. According to pat's video (https://www.youtube.com/watch?v=ulo7WatBEAo), this needs to be in a long format. So going to use reshape for this one and a simple melt. Column headers will then be: Var1, Var2, and value where var1=genome id, var2= sample id, value = rel abund. He uses a pivot long/wide but im more used to reshape2 library.

t_tmm_neithercoverm_long = melt(t_tmm_neithercoverm_mx) %>%
  as.data.frame() #The melt needs a matrix to work, so I start it off as a matrix and then transform it into a data frame because the rCLR function below needs it that way.

#First, just calculate the total abundances for each sample in the dataframe read in above. We'll need this to compare the strategies down the line.
group_count_coverM = t_tmm_neithercoverm_long %>%
  group_by(Var1) %>%
  summarize(totalabunds=sum(value))

#Now we make the geometric mean function.
gm=function(x){
  exp(mean(log(x[x>0])))
}

gm(c(2,3,4, 0)) #this is working. Its doing "robust" clr where it is only calculating based on numbers that are greater than 0.


#Now we make the full table that contains the rCLR values.
rclr_neithercoverm = t_tmm_neithercoverm_long %>% 
  group_by(Var1) %>%
  mutate(rclr=log(value/gm(value))) %>%
  ungroup() %>%
  select(-value) %>%
  pivot_wider(names_from = Var1, values_from=rclr, values_fill = 0) %>%
  as.data.frame()

#Cleaning up the dataframe a little.
rownames(rclr_neithercoverm) = rclr_neithercoverm$Var2
rclr_neithercoverm=rclr_neithercoverm[,-1]

rclr_neithercoverm[rclr_neithercoverm=="-Inf"]=0

#Then we just remake a data frame and transpose it so it matches the original non-clrd vegdist that I will calculate below.
rclr_neithercoverm=as.data.frame(rclr_neithercoverm)
rclr_neithercoverm_t=as.data.frame(t(rclr_neithercoverm))

#clr'd euclidean distance.
clrdist_coverm_neither=vegdist(rclr_neithercoverm_t, method="euclidean") 

#Now i can plot and compare the clr and non-clr distances.
noclrdist_coverm_tbl_neither = noclrdist_covermneither %>%
  as.matrix %>%
  as_tibble(rownames="sample") %>%
  pivot_longer(cols=-sample) %>%
  filter(name < sample)

clrdist_coverm_tbl_neither = clrdist_coverm_neither %>%
  as.matrix %>%
  as_tibble(rownames="sample") %>%
  pivot_longer(cols=-sample) %>%
  filter(name < sample)

inner_join(noclrdist_coverm_tbl_neither, clrdist_coverm_tbl_neither, by=c("sample", "name")) %>%
  select(sample, name, noclr=value.x, rclr=value.y) %>%
  inner_join(., group_count_coverM, by=c("sample" = "Var1")) %>%
  inner_join(., group_count_coverM, by=c("name" = "Var1")) %>%
  mutate(diffs=abs(totalabunds.x - totalabunds.y)) %>%
  pivot_longer(cols=c(noclr, rclr), names_to = "method", values_to="dist") %>%
  ggplot(aes(x=diffs,y=dist)) +
  geom_point()+
  facet_wrap(~method, scales="free_y")+
  geom_smooth()






















###################################################
####################################
######Out of curiosity, testing if different distance metrics are better at dealing with seq depth:
####################################
###################################################












t_tmm_swcovermcomp= as.data.frame(t(tmm_swcoverm)) #To have a df for vegdist
t_tmm_swcovermcomp_mx = as.matrix(t_tmm_swcovermcomp) #To have the matrix for melt.

#Calculated the euclidean distance of non-rclr non TMM data using euclidean and bray.
noclrdist_covermcompeuc=vegdist(t_tmm_swcovermcomp, method="euclidean") #non clr'd euclidean distance.
noclrdist_covermcompbray=vegdist(t_tmm_swcovermcomp, method="bray") #non clr'd bray distance.

#Melt to make long.
t_tmm_swcovermcomp_long = melt(t_tmm_swcovermcomp_mx) %>%
  as.data.frame() #The melt needs a matrix to work, so I start it off as a matrix and then transform it into a data frame because the rCLR function below needs it that way.

#Make the tibble for plot for bray
noclrdist_covermcompbray_tbl = noclrdist_covermcompbray %>%
  as.matrix %>%
  as_tibble(rownames="sample") %>%
  pivot_longer(cols=-sample) %>%
  filter(name < sample)

#Make the tibble for plot for euclidean
noclrdist_covermcompeuc_tbl = noclrdist_covermcompeuc %>%
  as.matrix %>%
  as_tibble(rownames="sample") %>%
  pivot_longer(cols=-sample) %>%
  filter(name < sample)

inner_join(noclrdist_covermcompbray_tbl, noclrdist_covermcompeuc_tbl, by=c("sample", "name")) %>%
  select(sample, name, bray=value.x, euc=value.y) %>%
  inner_join(., group_count_coverM, by=c("sample" = "Var1")) %>% #using the same group_count_coverM i calculated above since same exact "total counts".
  inner_join(., group_count_coverM, by=c("name" = "Var1")) %>% #using the same group_count_coverM i calculated above since same exact "total counts".
  mutate(diffs=abs(totalabunds.x - totalabunds.y)) %>%
  pivot_longer(cols=c(bray, euc), names_to = "method", values_to="dist") %>%
  ggplot(aes(x=diffs,y=dist)) +
  geom_point()+
  facet_wrap(~method, scales="free_y")+
  geom_smooth()





#lol. Both look equally terrible. Definitely just need to TMM no matter what distance metric I will use. Maybe slightly less noise on Bray. But veeery slightly. Just going to proceed with TMM'd aitchinsons.
