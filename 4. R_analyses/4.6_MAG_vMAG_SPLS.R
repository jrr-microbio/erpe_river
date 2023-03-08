##### Load relevant libraries
library(ggplot2)
library(pls)
library(dplyr)

source("4.6.3_VIP.R") ## This is a custom script (from the internet) to compute VIPs, has to be in the same directory as this script (or this path has to be changed)

##### Load metabolites/geochemistry data. Mine has 31 samples and these are rows.
metabolite<-t(read.csv("4.5.3_32_sample_geochemistry_incl_FTICR_subset_mod.csv",header=T, row.names=1))

metabolite<-as.matrix(metabolite) # force R to see it as a matrix

metabolite_SW=metabolite[18:32, ]
metabolite_PW=metabolite[1:17, ]

metabolite_SW<-as.matrix(metabolite_SW) # force R to see it as a matrix
dim(metabolite_SW) # this gives the matrix dimension

metabolite_PW<-as.matrix(metabolite_PW) # force R to see it as a matrix
dim(metabolite_PW) # this gives the matrix dimension

##### Load OTU data. The only difference between this file and the file used for WGCNA is that i have added in the module colors that resulted from the WGCNA analyses into each gneome. This way I can subset for doing SPLS.

sw_data<-read.table("4.6.1_WGCNA_input_MAGs-vMAGs_SW-both_modules.txt",header=T)

#sw_data<-as.matrix(sw_data) # force R to see it as a matrix
dim(sw_data) # this gives the matrix dimension

#Subset the modules
sw_data_salmon = sw_data %>%
  filter(module == "salmon")
sw_data_salmon = sw_data_salmon[ , -c(1)]
sw_data_salmon = as.matrix(t(sw_data_salmon))

sw_data_tan = sw_data %>%
  filter(module == "tan")
sw_data_tan = sw_data_tan[ , -c(1)]
sw_data_tan = as.matrix(t(sw_data_tan))

sw_data_turqoise = sw_data %>%
  filter(module == "turquoise")
sw_data_turqoise = sw_data_turqoise[ , -c(1)]
sw_data_turqoise = as.matrix(t(sw_data_turqoise))

sw_data_blue = sw_data %>%
  filter(module == "blue")
sw_data_blue = sw_data_blue[ , -c(1)]
sw_data_blue = as.matrix(t(sw_data_blue))

sw_data_brown = sw_data %>%
  filter(module == "brown")
sw_data_brown = sw_data_brown[ , -c(1)]
sw_data_brown = as.matrix(t(sw_data_brown))

sw_data_magenta = sw_data %>%
  filter(module == "magenta")
sw_data_magenta = sw_data_magenta[ , -c(1)]
sw_data_magenta = as.matrix(t(sw_data_magenta))

#Step 2: Predict values
th_r2<-0.1 # We will only look at the PLS if the correlation is better than 0.1
for (i in 1:ncol(metabolite_SW)){ # We treat each metabolite independently
  parameter<-colnames(metabolite_SW)[i]
  obs_values<-metabolite_SW[,i] # these are the observed values we'll try to predict
  print(paste("Trying to predict ",parameter," --- ",i,sep=""))
  # We perform the sPLS, trying to predict our specific metabolite vector (metabolite[,i]) using our whole OTU table. Validation is LOO, so "Leave one out", i.e. we train model on n-1 samples and try to predict the value for the remaining one. the "method" argument is to chose the correct type of sPLS
  pls_result<-plsr(obs_values ~ sw_data_brown, validation="LOO",method="oscorespls") 
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
i<-15 #this number corresponds to the "iteration" of the prediction. E.g. prediction 1, prediction 2, prediction 3, etc. In my case, the value with a significant SPLS value was the 15th, "NO3". Change it to whatever yours is.

# And we regenerate the corresponding results
parameter<-colnames(metabolite_SW)[i]
obs_values<-metabolite_SW[,i] # these are the observed values we'll try to predict
pls_result<-plsr(obs_values ~ sw_data_brown, validation="LOO",method="oscorespls") 
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
pdf(paste("measured_vs_predicted_","brown","-vs-",parameter,".pdf"))
ggplot(data=df) + geom_point(aes(x=x,y=y)) + geom_smooth(aes(x=x,y=y),method=lm) + xlab("Measured") + ylab("Predicted") + ggtitle(paste("Comparison of ",parameter," measured vs predicted -- r2=",max)) + theme(axis.text=element_text(color="black",size=10),axis.ticks=element_line(color="black"))
dev.off()

fit = lm(y ~ x, data = df)
summary(fit)
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
for (k in 1:8){
  OTU<-unlist(names(vip_components[k]))
  vec_OTU<-sw_data_brown[,OTU]
  df<-data.frame(x=vec_OTU,y=obs_values)
  print(ggplot(data=df) + geom_point(aes(x=x,y=y)) + geom_smooth(aes(x=x,y=y),method=lm) + xlab(OTU) + ylab(paste("Measured ",parameter,sep="")) + ggtitle(paste("Comparison of ",OTU," vs measured ",parameter)) + theme(axis.text=element_text(color="black",size=10),axis.ticks=element_line(color="black")))
}


######################################
######################################
################Now for the porewater.
######################################
######################################
######################################

##### Load OTU data. In my case, this contains all viruses with non zero abundances.
pw_data<-read.table("4.6.2_WGCNA_input_MAGs-vMAGs_PW-both_modules.txt",header=T)

#sw_data<-as.matrix(sw_data) # force R to see it as a matrix
dim(pw_data) # this gives the matrix dimension, so 28 vOTUs and 30 samples

#Going to subset some of the modules so I can use it on a per-module basis. Also transposing those and removing the module column to make it a matrix. There are 6 modules that are correlated to specific variables.
pw_data_blue = pw_data %>%
  filter(module == "blue")
pw_data_blue = pw_data_blue[ , -c(1)]
pw_data_blue = as.matrix(t(pw_data_blue))

pw_data_turqoise = pw_data %>%
  filter(module == "turquoise")
pw_data_turqoise = pw_data_turqoise[ , -c(1)]
pw_data_turqoise = as.matrix(t(pw_data_turqoise))

pw_data_brown = pw_data %>%
  filter(module == "brown")
pw_data_brown = pw_data_brown[ , -c(1)]
pw_data_brown = as.matrix(t(pw_data_brown))

#Step 2: Predict values

th_r2<-0.1 # We will only look at the PLS if the correlation is better than 0.1
for (i in 1:ncol(metabolite_PW)){ # We treat each metabolite independently
  parameter<-colnames(metabolite_PW)[i]
  obs_values<-metabolite_PW[,i] # these are the observed values we'll try to predict
  print(paste("Trying to predict ",parameter," --- ",i,sep=""))
  # We perform the sPLS, trying to predict our specific metabolite vector (metabolite[,i]) using our whole OTU table. Validation is LOO, so "Leave one out", i.e. we train model on n-1 samples and try to predict the value for the remaining one. the "method" argument is to chose the correct type of sPLS
  pls_result<-plsr(obs_values ~ pw_data_blue, validation="LOO",method="oscorespls") 
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
i<-26 #this number corresponds to the "iteration" of the prediction. E.g. prediction 1, prediction 2, prediction 3, etc.

# And we regenerate the corresponding results
parameter<-colnames(metabolite_PW)[i]
obs_values<-metabolite_PW[,i] # these are the observed values we'll try to predict
pls_result<-plsr(obs_values ~ pw_data_blue, validation="LOO",method="oscorespls") 
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
ggplot(data=df) + geom_point(aes(x=x,y=y)) + geom_smooth(aes(x=x,y=y),method=lm) + xlab("Measured") + ylab("Predicted") + ggtitle(paste("Comparison of ",parameter," measured vs predicted -- r2=",max)) + theme(axis.text=element_text(color="black",size=10),axis.ticks=element_line(color="black"))
dev.off()

# So next we checking the VIP (variable importance in projection), and output a table of the 100 highest values, this will tell us which OTU you need to know the abundance of to correctly predict the feature of interest

output<-paste("VIP_values_with_",parameter,".csv",sep="")
cat("Rank,OTU,VIP\n",file = output,append=FALSE)
vip_result<-VIP(pls_result)
vip_components<-sort(vip_result[max_comp,],decreasing=TRUE)[1:100]
for (k in 1:100){
  cat(paste(k,names(vip_components[k]),vip_components[k],"\n",sep=","),file=output,append=TRUE)
}


# Checking the VIP
output<-paste("VIP_values_with_",parameter,sep="")
vip_result<-VIP(pls_result)
vip_components<-sort(vip_result[max_comp,],decreasing=TRUE)[1:100]
for (i in 1:100){
  cat(paste("Rank ",i," we have ",names(vip_components[i])," with a VIP of ",vip_components[i],"\n",sep=""),file=output,append=TRUE)
}
weight_2 <- as.data.frame(weight[!is.na(weight)])
df<-data.frame(x=weight_2,y=pls_result$validation$pred[,,max_comp])
colnames(df)<-c("x","y")
pdf(paste("measured_vs_predicted_",module,"-vs-",parameter,".pdf"))
ggplot(data=df) + geom_point(aes(x=x,y=y)) + geom_smooth(aes(x=x,y=y),method=lm) + xlab("Measured") + ylab("Predicted") + ggtitle(paste("Comparison of ",parameter," measured vs predicted for module ",module)) + theme(axis.text=element_text(color="black",size=10),axis.ticks=element_line(color="black"))
dev.off()

## Check the correlation between predicted and modeled (should be consistent with what plsr gave us)
cor.test(df$x,df$y)
## Alternatively, we can also use the built-in function "predplot" but I find it less pretty
## Can be good to double check the ggplot2 plot though (should be the same)
predplot(pls_result,ncomp=max_comp)

## Now we can also plot individual OTUs vs n_per based on the VIP list
for (k in 1:5){
  OTU<-unlist(names(vip_components[k]))
  vec_OTU<-pw_data_blue[,OTU]
  df<-data.frame(x=vec_OTU,y=obs_values)
  print(ggplot(data=df) + geom_point(aes(x=x,y=y)) + xlab(OTU) + ylab(paste("Measured ",parameter,sep="")) + ggtitle(paste("Comparison of ",OTU," vs measured ",parameter)) + theme(axis.text=element_text(color="black",size=10),axis.ticks=element_line(color="black")))
}