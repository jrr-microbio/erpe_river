library(tidyverse)
library(vegan)
library(reshape2)

setwd("/Volumes/GoogleDrive/My Drive/University/Wrighton Lab; PhD/ERPE_Timeseries_viral_communities/viral_work/read_mapping/95ID_noambig_filtered/clr_abunds/3samp_cutoff_threshold_clr/")


###################################################
####################################
#Since we are splitting up the compartments for rCLR, I will first start with the SW
####################################
###################################################


#Read in the data. Rows need to be OTUs. Columns need to be samples.
tmm_sw=read.csv("/Volumes/GoogleDrive/My Drive/University/Wrighton Lab; PhD/ERPE_Timeseries_viral_communities/viral_work/read_mapping/95ID_noambig_filtered/clr_abunds/3samp_cutoff_threshold_clr/1007_vMAGs_75AF_3xcov_counts_TMM_normalized_noambig_SW_3samp.csv", header=T, row.names =1)
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

write.table(rclr_sw_t, '1007_vMAGs_75AF_3xcov_vMAG_rclr_TMM_SW.tsv', sep='\t', col.names=NA)






###################################################
####################################
#Now i do the same as above but for the PW.
####################################
###################################################






#Read in the data. Rows need to be OTUs. Columns need to be samples.
tmm_pw=read.csv("/Volumes/GoogleDrive/My Drive/University/Wrighton Lab; PhD/ERPE_Timeseries_viral_communities/viral_work/read_mapping/95ID_noambig_filtered/clr_abunds/3samp_cutoff_threshold_clr/191_vMAGs_75AF_3xcov_counts_TMM_normalized_noambig_PW_3samp.csv", header=T, row.names =1)

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

write.table(rclr_pw_t, '191_vMAGs_75AF_3xcov_vMAG_rclr_TMM_pw.tsv', sep='\t', col.names=NA)

