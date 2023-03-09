## WGCNA networks co-expression data

# load 
library(WGCNA)
library(tidyverse)
library(missForest)
library(mixOmics)
library(pheatmap)
library(vegan)

###############################################################################
#Load data ###
###############################################################################

# Note: For this analysis, we had to split up the different WGCNA networks because the compartments were so different. In order to make to WGCNA networks, I selected each individual genome that was present in >3 samples for the samples in either comparmtent (SW and PW). I then made a subset of the original abundance files to use. This is why there are 2 new input files here- however they are identical to the TMM normalized abundances on Tables S2 and S3.

adjs_active = read.delim('4.5.1_WGCNA_input_MAGs-vMAGs_SW_and_both.txt',sep='\t', row.names=1) #Load the abundance matrix. Rows = samples, columns = species/genes/metabolites/etc.
adjs_active.t <- as.data.frame(t(adjs_active)) # set up data frame (e.g. rows = samples, columns = spp/genes/metabolites etc.)

#This file is a modified version of "4.1.1_32_sample_geochemistry_incl_FTICR_subset.txt" that is 1) transposed and 2) has the first two columns "diel cycle", "sw_pw", "timepoint" columns removed. 
metabolite<-t(read.csv("4.5.3_32_sample_geochemistry_incl_FTICR_subset_mod.csv",header=T, row.names=1)) # Load metabolites data / same as the metadata / if you don't have this data, this can be skipped.

metabolite<-as.matrix(metabolite) # force R to see it as a matrix
metabolite_SW=metabolite[18:32, ] #Subset for only SW
metabolite_SW = as.data.frame(metabolite_SW)  # force R to see it as a matrix

###############################################################################
# (1) Overall metaG network ####
###############################################################################

adjs_active.t_log<-log(adjs_active.t+1) # Natural log (default)

# Choosing the soft-thresholding power: analysis of network topology
allowWGCNAThreads() #this will enable multi threading.

powers = c(c(1:10), seq(from = 12, to=20, by=2)) # Choose a set of soft-thresholding powers

sft = pickSoftThreshold(adjs_active.t_log, powerVector = powers, verbose = 5,
                        networkType = "signed hybrid") # Call the network topology analysis function

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.80,col="red") # this line corresponds to using an R^2 cut-off of h


# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

power_thresh = 14 # Set power from flattening of curve above

# One-step network construction and module detection
set.seed(1234)
net = blockwiseModules(adjs_active.t_log, power = power_thresh,
                       TOMType = "signed", networkType = "signed hybrid",
                       minModuleSize = 20, ##
                       reassignThreshold = 0, mergeCutHeight = 0.3, ## 0.25
                       numericLabels = F, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE, 
                       verbose = 3)

# inspect modules
table(net$colors)
net$colors

# open a graphics window
sizeGrWindow(12, 9)

# Convert labels to colors for plotting if you selected numericLabels = T
#mergedColors = labels2colors(net$colors)
#table(mergedColors)

moduleColors = net$colors
moduleLabels = net$colors

# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

MEs = net$MEs;
geneTree = net$dendrograms[[1]];
nPop = ncol(adjs_active.t_log)
nSamples = nrow(adjs_active.t_log)

#####################################
#####################################
#####################################
#####################################
########## I feel like i can cluster these modules a little more. Using this protocol: https://fuzzyatelin.github.io/bioanth-stats/module-F21-Group1/module-F21-Group1.html#:~:text=For%20genomic%20data%20like%20this,start%20at%20a%20minClusterSize%20%3D%2030%20.
# 
# 

# MElist_toclust = moduleEigengenes(adjs_active.t_log, colors=moduleColors)
# MEs_toclust = MElist_toclust$eigengenes
# ME.dissimilarity_toclust = 1-cor(MElist_toclust$eigengenes, use="complete") #Calculate eigengene dissimilarity
# METree_toclust = hclust(as.dist(ME.dissimilarity_toclust), method = "average") #Clustering eigengenes 
# par(mar = c(0,4,2,0)) #seting margin sizes
# par(cex = 0.6);#scaling the graphic
# plot(METree_toclust)
# abline(h=.25, col = "red") #a height of .35 corresponds to correlation of .65
# merged_clusts <- mergeCloseModules(adjs_active.t_log, moduleColors, cutHeight = .25)
# # The merged module colors, assigning one color to each module
# mergedColors = merged_clusts$colors
# # Eigengenes of the new merged modules
# mergedMEs = merged_clusts$newMEs
# plotDendroAndColors(geneTree, cbind(moduleColors, mergedColors), 
#                     c("Original Module", "Merged Module"),
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05,
#                     main = "Gene dendrogram and module colors for original and merged modules")
# 
# table(net$mergedColors)
#Nothing changes so I'm just removing this out.
#####################################
#####################################
#####################################
#####################################

#I like to write out the modules and the colors of each of the organisms. So im just rerunning the net command above to get them.
write.csv(file="Module_eigen_values_surface_MAG-vMAG_signedhybrid.csv",MEs)
write.csv(file="Module_composition_surface_MAG-vMAG_signedhybrid.csv",net$colors)

# Calculate MEs with color labels
MEs0 = moduleEigengenes(adjs_active.t_log, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, metabolite_SW, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
sapply(moduleTraitPvalue, class)

#Filter out and replace any p-value that is greater than 0.05 with a zero i can use to subset.
moduleTraitCor[moduleTraitCor > -0.5 & moduleTraitCor < 0.5] <- 0
moduleTraitPvalue[moduleTraitPvalue > 0.05] <- 0
textMatrix = paste(signif(moduleTraitCor, 2), " (",signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)

# Will display correlations and their p-values
sizeGrWindow(10,10)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
par(mar = c(12, 5, 2, 2))

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(metabolite_SW),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# display with pheatmap?
moduleTraitCor2 <- as.matrix((moduleTraitCor))

pheatmap(moduleTraitCor2, 
         cellwidth = 10,
         cellheight = 10,
         color = colorRampPalette(colors = c("blue",
                                                   "#78d0e7",
                                                   "white",
                                                   "#ffd866",
                                                   "red"))(7),
                                                   border_color = "gray",
         cluster_rows = T,
         show_colnames = T,
         show_rownames = T,
         treeheight_col=T)

########################################
#Now plotting genome significant and subnetwork membership
#######################################
#Tutorial can be found here: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/a
# Calculate topological overlap
TOM = TOMsimilarityFromExpr(adjs_active.t_log, power = power_thresh,
                            networkType = "signed hybrid")

# Select modules, probes
modNames=substring(names(MEs),3)
modules = modNames
probes = colnames(adjs_active.t_log)
inModule = is.finite(match(moduleColors, modules))
modProbes = probes[inModule]

# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)

# Export the network into edge and node list files Gephi can read. Use a threshold of 0.02 to elminate everything that is poorly connected
cyt = exportNetworkToCytoscape(modTOM,
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule])
#Write out the files that Gephi needs. Each genome with its module color.
nodes <- cyt$nodeData %>% 
  rename(Id = nodeName,
         wgcna_module = `nodeAttr[nodesPresent, ]`) %>%
  dplyr::select(-altName)

edges <- cyt$edgeData %>%
  rename(Source = fromNode,
         Target = toNode,
         Weight = weight,
         Type = direction) %>%
  dplyr::select(-fromAltName, -toAltName)

# write out network results
write.csv(nodes, "SW_Gephi_network_nodes.csv", row.names = F)
write.csv(edges, "SW_Gephi_network_edges.csv", row.names = F)

###########
#After running SPLS, we saw that the brown module was highly significant for 
#I am only going to focus on the brown module and the nitrate concentrations.
adj_table_scaled <- as.data.frame(scale(adjs_active.t_log, center = TRUE, scale = TRUE))

#Trying to compare module membership vs genome significance
module = "brown"
  moduleGenes_brown = nodes %>%
    filter(wgcna_module=="brown")
  abund_brown = adj_table_scaled %>%
    dplyr::select(moduleGenes_brown$Id)
  
  modNames=substring(names(MEs),3)
  geneModuleMembership = as.data.frame(cor(adjs_active.t_log, MEs, use = "p"))
  
  #Now correlate that to the environmental variable: No3
  metabolite_SW_df = as.data.frame(metabolite_SW)
  metabolite_SW_NO3 = as.data.frame(metabolite_SW_df$NO3)
  names(metabolite_SW_NO3) = "Nitrate"
  geneTraitSignificance = as.data.frame(cor(adjs_active.t_log, metabolite_SW_NO3, use = "p"))
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
  names(geneTraitSignificance) = paste("GS.", names(metabolite_SW_NO3), sep="")
  names(GSPvalue) = paste("p.GS.", names(metabolite_SW_NO3), sep="")
  
  # Now select module and plot - NO3
  column = match(module, modNames)
  module_memb = abs(geneModuleMembership[moduleGenes_brown$Id, column])
  trait_sig = abs(geneTraitSignificance[moduleGenes_brown$Id, 1])
  brown_no3 = as.data.frame(cbind(module_memb, trait_sig))
  cor.test(brown_no3$module_memb, brown_no3$trait_sig, method = "pearson")
  brown_corr = cor(metabolite_SW_df$NO3,
                   abund_brown, use= "pairwise.complete.obs")  %>% t() 
  brown_no3$no3_corr = brown_corr
  
  
  #Assign labels for specific points of panel E:
  moduleGenes_brown$label[moduleGenes_brown$Id == "erpe_2018_sw_3hr_WHONDRS-SW48-000086_B_bin.13"] = "Proteobacteria host"
  moduleGenes_brown$label[moduleGenes_brown$Id == "erpe_2018_sw_21hr_WHONDRS-SW48-000092_A_NODE_26_length_91961_cov_18.429808__full-cat_1"] = "Linked virus"
  
  plt.brown_no3 = ggplot(brown_no3, aes(x=module_memb, y=trait_sig, color = brown_no3$color)) +
    geom_point(alpha=1, shape=21, stroke=0.25,
               fill="brown", color="black", aes(size= brown_corr)) +
    scale_size(trans="reverse") +
    theme_classic() + 
    geom_smooth(method="lm", se=T, color = "black") +
    geom_text(label = moduleGenes_brown$label, size = 0.5, hjust = 0, vjust = 0) +
    labs(x = "Subnetwork membership",
         y = "GS for nitrate",
         title = "brown module: 123 genomes\n cor=0.39, p=9.75e-06")
  
  plt.brown_no3
  
  ## WGCNA networks co-expression data
  
  # load 
  library(WGCNA)
  library(tidyverse)
  library(missForest)
  library(mixOmics)
  library(pheatmap)
  library(vegan)
  
  ###############################################################################
  #Load data ###
  ###############################################################################
  adjs_active = read.delim('4.5.2_WGCNA_input_MAGs-vMAGs_PW_and_both.txt',sep='\t', row.names=1) #Load the abundance matrix. Rows = samples, columns = species/genes/metabolites/etc.
  adjs_active.t <- as.data.frame(t(adjs_active)) # set up data frame (e.g. rows = samples, columns = spp/genes/metabolites etc.)
  
  #This file is a modified version of "4.1.1_32_sample_geochemistry_incl_FTICR_subset.txt" that is 1) transposed and 2) has the first two columns "diel cycle", "sw_pw", "timepoint" columns removed. 
  metabolite<-t(read.csv("4.5.3_32_sample_geochemistry_incl_FTICR_subset_mod.csv",header=T, row.names=1)) # Load metabolites data / same as the metadata / if you don't have this data, this can be skipped.
  metabolite<-as.matrix(metabolite) # force R to see it as a matrix
  metabolite_PW=metabolite[1:17,] #Subset for only SW
  metabolite_PW = as.data.frame(metabolite_PW)  # force R to see it as a matrix
  
  ###############################################################################
  # (1) Overall metaG network ####
  ###############################################################################
  adjs_active.t_log<-log(adjs_active.t+1) # Natural log (default)
  
  # Choosing the soft-thresholding power: analysis of network topology
  allowWGCNAThreads() #this will enable multi threading.
  
  powers = c(c(1:10), seq(from = 12, to=20, by=2)) # Choose a set of soft-thresholding powers
  
  sft = pickSoftThreshold(adjs_active.t_log, powerVector = powers, verbose = 5,
                          networkType = "signed hybrid") # Call the network topology analysis function
  
  # Plot the results:
  sizeGrWindow(9, 5)
  par(mfrow = c(1,2))
  cex1 = 0.9
  
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")
  abline(h=0.6,col="red") # this line corresponds to using an R^2 cut-off of h
  
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  
  power_thresh = 8 # Set power from flattening of curve above
  
  # One-step network construction and module detection
  set.seed(1234)
  net = blockwiseModules(adjs_active.t_log, power = power_thresh,
                         TOMType = "signed", networkType = "signed hybrid",
                         minModuleSize = 20, ##
                         reassignThreshold = 0, mergeCutHeight = 0.3, ## 0.25
                         numericLabels = F, pamRespectsDendro = FALSE,
                         saveTOMs = TRUE, 
                         verbose = 3)
  
  # inspect modules
  table(net$colors)
  net$colors
  
  # open a graphics window
  sizeGrWindow(12, 9)
  
  # Convert labels to colors for plotting if you selected numericLabels = T
  #mergedColors = labels2colors(net$colors)
  #table(mergedColors)
  
  moduleColors = net$colors
  moduleLabels = net$colors
  
  # Plot the dendrogram and the module colors underneath
  plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  
  MEs = net$MEs;
  geneTree = net$dendrograms[[1]];
  nPop = ncol(adjs_active.t_log)
  nSamples = nrow(adjs_active.t_log)
  
  # ######################################
  # ########## I feel like i can cluster these modules a little more. Using this protocol: https://fuzzyatelin.github.io/bioanth-stats/module-F21-Group1/module-F21-Group1.html#:~:text=For%20genomic%20data%20like%20this,start%20at%20a%20minClusterSize%20%3D%2030%20.
  # #####################################
  # 
  # MElist_toclust = moduleEigengenes(adjs_active.t_log, colors=moduleColors)
  # MEs_toclust = MElist_toclust$eigengenes
  # ME.dissimilarity_toclust = 1-cor(MElist_toclust$eigengenes, use="complete") #Calculate eigengene dissimilarity
  # METree_toclust = hclust(as.dist(ME.dissimilarity_toclust), method = "average") #Clustering eigengenes 
  # par(mar = c(0,4,2,0)) #seting margin sizes
  # par(cex = 0.6);#scaling the graphic
  # plot(METree_toclust)
  # abline(h=.25, col = "red") #a height of .25 corresponds to correlation of .75
  
  
  #I like to write out the modules and the colors of each of the organisms. So im just rerunning the net command above to get them.
  write.csv(file="Module_eigen_values_pore_MAG-vMAG_signedhybrid.csv",MEs)
  write.csv(file="Module_composition_pore_MAG-vMAG_signedhybrid.csv",net$colors)
  # Recalculate MEs with color labels
  MEs0 = moduleEigengenes(adjs_active.t_log, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  moduleTraitCor = cor(MEs, metabolite_PW, use = "p")
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
  sapply(moduleTraitPvalue, class)
  
  #Filter out and replace any p-value that is greater than 0.05 with a zero i can use to subset.
  moduleTraitCor[moduleTraitCor > -0.5 & moduleTraitCor < 0.5] <- 0
  moduleTraitPvalue[moduleTraitPvalue > 0.05] <- 0
  textMatrix = paste(signif(moduleTraitCor, 2), " (",signif(moduleTraitPvalue, 1), ")", sep = "")
  dim(textMatrix) = dim(moduleTraitCor)
  
  # Will display correlations and their p-values
  sizeGrWindow(10,10)
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "")
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(6, 8.5, 3, 3))
  par(mar = c(12, 5, 2, 2))
  
  # Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(metabolite_PW),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  
  # display with pheatmap?
  moduleTraitCor2 <- as.matrix((moduleTraitCor))
  
  pheatmap(moduleTraitCor2, 
           cellwidth = 10,
           cellheight = 10,
           color = colorRampPalette(colors = c("blue",
                                                     "#78d0e7",
                                                     "white",
                                                     "#ffd866",
                                                     "red"))(7),
                                                     border_color = "gray",
           cluster_rows = T,
           show_colnames = T,
           show_rownames = T,
           treeheight_col=T)
  
  ########################################
  #Now plotting genome significant and subnetwork membership
  #######################################
  
  # Recalculate topological overlap if needed
  TOM = TOMsimilarityFromExpr(adjs_active.t_log, power = power_thresh,
                              networkType = "signed hybrid")
  
  # Select modules, probes
  modules = modNames
  probes = colnames(adjs_active.t_log)
  inModule = is.finite(match(moduleColors, modules))
  modProbes = probes[inModule]
  
  # Select the corresponding Topological Overlap
  modTOM = TOM[inModule, inModule]
  dimnames(modTOM) = list(modProbes, modProbes)
  
  # Export the network into edge and node list files Gephi can read
  cyt = exportNetworkToCytoscape(modTOM,
                                 weighted = TRUE,
                                 threshold = 0.02,
                                 nodeNames = modProbes,
                                 nodeAttr = moduleColors[inModule])
  
  #Write out the files that Gephi needs. Each genome with its module color.
  nodes <- cyt$nodeData %>% 
    rename(Id = nodeName,
           wgcna_module = `nodeAttr[nodesPresent, ]`) %>%
    dplyr::select(-altName)
  
  nodes <- cyt$nodeData 
  edges <- cyt$edgeData %>%
    rename(Source = fromNode,
           Target = toNode,
           Weight = weight,
           Type = direction) %>%
    dplyr::select(-fromAltName, -toAltName)
  
  # write out network results
  write.csv(nodes, "PW_Gephi_network_nodes.csv", row.names = F)
  write.csv(edges, "PW_Gephi_network_edges.csv", row.names = F)
