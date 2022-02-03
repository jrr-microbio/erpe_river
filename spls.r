########VIP.R file

### VIP.R: Implementation of VIP (variable importance in projection)(*) for the
### `pls' package.
### $Id: VIP.R,v 1.2 2007/07/30 09:17:36 bhm Exp $

### Copyright © 2006,2007 Bjørn-Helge Mevik
### This program is free software; you can redistribute it and/or modify
### it under the terms of the GNU General Public License version 2 as
### published by the Free Software Foundation.
###
### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
### GNU General Public License for more details.

### A copy of the GPL text is available here:
### http://www.gnu.org/licenses/gpl-2.0.txt

### Contact info:
### Bjørn-Helge Mevik
### bhx6@mevik.net
### Rødtvetvien 20
### N-0955 Oslo
### Norway

### (*) As described in Chong, Il-Gyo & Jun, Chi-Hyuck, 2005, Performance of
### some variable selection methods when multicollinearity is present,
### Chemometrics and Intelligent Laboratory Systems 78, 103--112.

## VIP returns all VIP values for all variables and all number of components,
## as a ncomp x nvars matrix.
VIP <- function(object) {
  if (object$method != "oscorespls")
    stop("Only implemented for orthogonal scores algorithm.  Refit with 'method = \"oscorespls\"'")
  if (nrow(object$Yloadings) > 1)
    stop("Only implemented for single-response models")
  
  SS <- c(object$Yloadings)^2 * colSums(object$scores^2)
  Wnorm2 <- colSums(object$loading.weights^2)
  SSW <- sweep(object$loading.weights^2, 2, SS / Wnorm2, "*")
  sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
}


## VIPjh returns the VIP of variable j with h components
VIPjh <- function(object, j, h) {
  if (object$method != "oscorespls")
    stop("Only implemented for orthogonal scores algorithm.  Refit with 'method = \"oscorespls\"'")
  if (nrow(object$Yloadings) > 1)
    stop("Only implemented for single-response models")
  
  b <- c(object$Yloadings)[1:h]
  T <- object$scores[,1:h, drop = FALSE]
  SS <- b^2 * colSums(T^2)
  W <- object$loading.weights[,1:h, drop = FALSE]
  Wnorm2 <- colSums(W^2)
  sqrt(nrow(W) * sum(SS * W[j,]^2 / Wnorm2) / sum(SS))
}


################### SPLS run

##### Load relevant libraries
library(ggplot2)
library(pls)
library(ggdark)

setwd("/Volumes/GoogleDrive/My Drive/University/Wrighton Lab; PhD/ERPE_Timeseries_viral_communities/SPLS/vMAGs_and_MAGs") #whereve you are going to work. VIP has to be here.

source("/Volumes/GoogleDrive/My Drive/University/Wrighton Lab; PhD/ERPE_Timeseries_viral_communities/SPLS/VIP.R") ## This is a custom script (from the internet) to compute VIPs, has to be in the same directory as this script (or this path has to be changed)

##### Load OTU data. In my case, this contains all viruses with non zero abundances.

all_data_vOTU<-read.table("/Volumes/GoogleDrive/My Drive/University/Wrighton Lab; PhD/ERPE_Timeseries_viral_communities/SPLS/vMAGs_and_MAGs/1050_and_115_vMAG+MAG_25sample_75AF_1xcov_counts_TMM_normalized_nozero.txt",header=T, sep="\t")

# By default, R doesn't remove the first column if it has an ID, so we do it instead
row.names(all_data_vOTU)<-all_data_vOTU[,1]
all_data_vOTU<-all_data_vOTU[,-1] 
all_data_vOTU<-as.matrix(all_data_vOTU) # force R to see it as a matrix
dim(all_data_vOTU) # this gives the matrix dimension, so 28 vOTUs and 30 samples

## now we transpose the matrix, because sPLS wants the samples as rows
all_data_vOTU<-t(all_data_vOTU) #transpose

##### Load metabolites/geochemistry data. Mine has 31 samples and these are rows.
metabolite<-read.csv("/Volumes/GoogleDrive/My Drive/University/Wrighton Lab; PhD/ERPE_Timeseries_viral_communities/SPLS/25_samples_chemistry_for_SPLS.csv",header=T)

#same as before, write in the row names correctly.
row.names(metabolite)<-metabolite[,1]
metabolite<-metabolite[,-1]
metabolite<-as.matrix(metabolite) # force R to see it as a matrix
dim(metabolite) # this gives the matrix dimension, so 19 metabolites and 33 samples


#Step 2: Predict values

th_r2<-0.1 # We will only look at the PLS if the correlation is better than 0.1
for (i in 1:ncol(metabolite)){ # We treat each metabolite independently
  parameter<-colnames(metabolite)[i]
  obs_values<-metabolite[,i] # these are the observed values we'll try to predict
  print(paste("Trying to predict ",parameter," --- ",i,sep=""))
  # We perform the sPLS, trying to predict our specific metabolite vector (metabolite[,i]) using our whole OTU table. Validation is LOO, so "Leave one out", i.e. we train model on n-1 samples and try to predict the value for the remaining one. the "method" argument is to chose the correct type of sPLS
  pls_result<-plsr(obs_values ~ all_data_vOTU, validation="LOO",method="oscorespls") 
  # Now we check the vector of r2 (sPLS tries to use different numbers of OTUs and provides a correlation between predicted and observed for each of them, so we get a vectore of r2 and not just one r2 value)
  r2_vector<-R2(pls_result)
  max<-0
  max_comp<--1
  for (j in 1:length(r2_vector$val)){
    if(!(is.na(r2_vector$val[j]))){
      if(r2_vector$val[j]>th_r2){
        if(r2_vector$val[j]>max){
          max<-r2_vector$val[j]
          max_comp<-r2_vector$comp[j]
        }
      }
    }
  }
  print(paste(" the max r2 is ",max," corresponding to comp ",max_comp,sep=""))
}

# So here we print the highest r2 across all predictions.


#Step 3: Plotting the results

## So now we can look at the metabolites, and see if any can be predicted by the OTUs abundance
## Looks like sulfate is not bad (r2 0.36) so we set i to the corresponding column (3)
i<-13 #this number corresponds to the "iteration" of the prediction. E.g. prediction 1, prediction 2, prediction 3, etc. In my case, the value with a significant SPLS value was the 13th, "na". Change it to whatever yours is.

# And we regenerate the corresponding results
parameter<-colnames(metabolite)[i]
obs_values<-metabolite[,i] # these are the observed values we'll try to predict
pls_result<-plsr(obs_values ~ all_data_vOTU, validation="LOO",method="oscorespls") 
r2_vector<-R2(pls_result)
max<-0
max_comp<--1
for (j in 1:length(r2_vector$val)){
  if(!(is.na(r2_vector$val[j]))){
    if(r2_vector$val[j]>th_r2){
      if(r2_vector$val[j]>max){
        max<-r2_vector$val[j]
        max_comp<-r2_vector$comp[j]
      }
    }
  }
}


# Plotting predicted vs observed
df<-data.frame(x=obs_values,y=pls_result$validation$pred[,,max_comp])
colnames(df)<-c("x","y")
pdf(paste("measured_vs_predicted_",module,"-vs-",parameter,".pdf"))
ggplot(data=df, label=) + geom_point(aes(x=x,y=y), size=9) + geom_smooth(aes(x=x,y=y),method=lm) + xlab("Measured") + ylab("Predicted") + ggtitle(paste("Comparison of ",parameter," measured vs predicted -- r2=",max)) + theme(axis.text=element_text(color="black",size=10),axis.ticks=element_line(color="black"))+dark_theme_gray()
dev.off()
## So r2 is not "great", but there seems to be some common features to all the "high-sulfate" samples (the two above 300), as well as some of the "medium sulfate" (on the right side of the plot, > 300 but predicted around 100-200), while all the "low-sulfate" samples (bottom left, <100) are also well predicted.
# So next we checking the VIP (variable importance in projection), and output a table of the 100 highest values, this will tell us which OTU you need to know the abundance of to correctly predict the feature of interest
output<-paste("VIP_values_with_",parameter,".csv",sep="")
cat("Rank,OTU,VIP\n",file = output,append=FALSE)
vip_result<-VIP(pls_result)
vip_components<-sort(vip_result[max_comp,],decreasing=TRUE)[1:100]
for (k in 1:100){
  cat(paste(k,names(vip_components[k]),vip_components[k],"\n",sep=","),file=output,append=TRUE)
}
## Check the correlation between predicted and modeled (should be consistent with what plsr gave us)
cor.test(df$x,df$y)
## Alternatively, we can also use the built-in function "predplot" but I find it less pretty
## Can be good to double check the ggplot2 plot though (should be the same)
predplot(pls_result,ncomp=max_comp)

## Now we can also plot individual OTUs vs n_per based on the VIP list
for (k in 1:5){
  OTU<-unlist(names(vip_components[k]))
  vec_OTU<-all_data_vOTU[,OTU]
  df<-data.frame(x=vec_OTU,y=obs_values)
  print(ggplot(data=df) + geom_point(aes(x=x,y=y)) + xlab(OTU) + ylab(paste("Measured ",parameter,sep="")) + ggtitle(paste("Comparison of ",OTU," vs measured ",parameter)) + theme(axis.text=element_text(color="black",size=10),axis.ticks=element_line(color="black")))+dark_theme_gray()
}
