library(ggplot2)
library(vegan)
library(grid)
library(MASS)
library(RColorBrewer)
library(purrr)
library(dplyr)
library(tidyr)
library(plotly)
library(htmlwidgets)
library(mds)
library(ggdark)
library(tidyverse)
library(patchwork)

setwd("/Volumes/Macintosh HD/Users/josue.rodriguez/Library/CloudStorage/GoogleDrive-jarora2213@gmail.com/My Drive/University/wrighton_lab_phd/erpe_timeseries_viral_communities/_paper_figures/NMDS_and_PCA/NMDS_plots/MAG_NMDS/")

##read in feature table with species as columns and samples as rows
genome<-t(read.table('/Volumes/Macintosh HD/Users/josue.rodriguez/Library/CloudStorage/GoogleDrive-jarora2213@gmail.com/My Drive/University/wrighton_lab_phd/erpe_timeseries_viral_communities/_paper_figures/NMDS_and_PCA/NMDS_plots/MAG_NMDS/32_sample/120_MAGs_3x_depth_75AF_TMM_normalized.tsv', sep="\t", header=T,check.names=TRUE, row.names = 1))

##read in chemistry (nona)
chem = read.table('/Volumes/Macintosh HD/Users/josue.rodriguez/Library/CloudStorage/GoogleDrive-jarora2213@gmail.com/My Drive/University/wrighton_lab_phd/erpe_timeseries_viral_communities/_paper_figures/NMDS_and_PCA/32_sample_geochemistry_incl_FTICR_subset.txt', sep = '\t', header = TRUE, check.names = T)
rownames(chem)=chem[,1]
chem = chem[,-1]

log_chem=log(chem[,4:29])
chem$sw_pw[which(chem$sw_pw=="pw")] <- "1"
chem$sw_pw[which(chem$sw_pw=="sw")] <- "2"
chem$timepoint=factor(chem$timepoint)

### Calculate species accumuation curve, richness and shannon's H (a helpful guide: https://kembellab.ca/r-workshop/biodivR/SK_Biodiversity_R.html , and https://rpubs.com/an-bui/vegan-cheat-sheet)
genome_sw = as.data.frame(genome[1:15,])
genome_pw = as.data.frame(genome[16:32,])
accumcurve_sw = specaccum(genome_sw, method = "random", permutations=999)
accumcurve_pw = specaccum(genome_pw, method = "random", permutations=999)

plot(accumcurve_pw, add=FALSE, ci.type="poly", col="black", lwd=2, ci=0.2, ci.lty=1, ci.col="#F8766D", ylim=c(0,100), xlab = "Samples", ylab = 'mag OTUs')
plot(accumcurve_sw, ci.type="poly", col="black", lwd=2, ci=0.2, ci.lty=1, ci.col="#00BFC4", add=TRUE, ylim=c(0,100))

mean(specnumber(genome_sw))
mean(specnumber(genome_pw))
sd(specnumber(genome_sw))
sd(specnumber(genome_pw))
boxplot(specnumber(genome) ~ chem$sw_pw, ylab = "# of species", col = c("#F8766D", "#00BFC4"))
t.test(specnumber(genome) ~ chem$sw_pw)


mean(diversity(genome_sw, index = "shannon"))
mean(diversity(genome_pw, index = "shannon"))
sd(diversity(genome_sw, index = "shannon"))
sd(diversity(genome_pw, index = "shannon"))
boxplot(diversity(genome) ~ chem$sw_pw, ylab = "Shannon's H'", col = c("#F8766D", "#00BFC4")) 
t.test(diversity(genome) ~ chem$sw_pw)

####Simple Ord with env factors
chem.pca_log<-prcomp(na.omit(log_chem), center=TRUE, scale.=TRUE)
ordMAG<-metaMDS(genome, distance = "bray")
fit_logMAG <- envfit(ordMAG, log_chem, perm = 999, na.rm=TRUE)
scores(fit_logMAG, "vectors")

plot(ordMAG)
plot(fit_logMAG)
plot(fit_logMAG, p.max = 0.05, col = "green")

#This establishes the coordinates for the samples and all OTUS onto an non-metric dimensional scaling
Ord_dist <-metaMDSdist(genome, distance = "bray", noshare = 0.1, trace = 1, autotransform=T)

#Run mrpp and ANOSIM on sw vs pw between samples
mrpp(genome, chem$sw_pw, permutations=999, distance="bray")
anosim(Ord_dist, chem$sw_pw, permutations = 999)

#########################
##Now to work on plotting the NMDS on ggplot in 2D (k = 3 gives 3D). Next 5 lines just making sure the plot still looks good with added features on metaMDS command.
NMDS_Bray_genome <-metaMDS(genome, distance = "bray", k =2,
                          noshare = 0.1, trace = 1, trymax = 500, autotransform=T)
ord.genome = as.data.frame(scores((NMDS_Bray_genome), display="sites"))
#ord.genome$sampleid=row.names(ord.genome)

stressplot(NMDS_Bray_genome)

#Now to start actually plotting this on ggplot.
#Pull out the scores with the new R framework
NMDS_Bray_genome <-metaMDS(genome, distance = "bray", k =2,
                              noshare = 0.1, trace = 1, trymax = 500, autotransform=T)

data.scores = as.data.frame(scores(NMDS_Bray_genome)$sites)
data.scores$chemistry = chem$sw_pw
en = envfit(NMDS_Bray_genome, log_chem, permutations = 999, na.rm = TRUE)
en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)
en_coord_cont$Compound = rownames(en_coord_cont)

###########Now do the final plot with the loadings as bars.
distance <- function(x, y, home = c(0,0)) {
  sqrt((x-home[1])^2 + (y-home[2])^2)
}

#fit_log <- envfit(NMDS_Bray_genome, log_chem, perm = 999, na.rm=TRUE)

#create loadings df for plot arrows
v_scrs <- as.data.frame(scores(en, display = "vectors"))
v_scrs <- cbind(v_scrs*ordiArrowMul(en), Compound = rownames(v_scrs))

v_scrs <- v_scrs %>%
  mutate(dist_pc1pc2 = distance(NMDS1, NMDS2), #distance in space from center
         scaled_dist = dist_pc1pc2 * 10)

#get top 10 compounds based on distance in space from center
top_dist_compounds <- v_scrs %>%
  arrange(desc(scaled_dist)) %>%
  slice_head(n = 10) %>%
  mutate(num_lab = 1:10) %>%
  select(Compound, num_lab)

#make plotting columns for loadings df
v_scrs <- v_scrs %>%
  mutate(top_c = ifelse(Compound %in% top_dist_compounds$Compound,"top","other"),
         alpha_factor = ifelse(top_c == "top", 1, 0.5),
         plot_pc1 = ifelse(top_c == "top", NMDS1, 0),
         plot_pc2 = ifelse(top_c == "top", NMDS2, 0),
         top_label = ifelse(top_c == "top", Compound, "")) %>%
  left_join(top_dist_compounds, by = c("Compound"))


#Now plot on ggplot
sw_pw <- data.scores %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_segment(data = rev(v_scrs),
               aes(x= 0, xend = NMDS1*10, 
                   y = 0, yend = NMDS2*10, 
                   alpha = alpha_factor,
                   color = Compound),
               arrow = arrow(length = unit(0.5, "cm")),
               inherit.aes = FALSE) +
  scale_color_manual(values=wes_palette("FantasticFox1", type="continuous", 26)) + 
  geom_text(data = v_scrs,
            aes(label = num_lab,
                x = NMDS1*10,
                y = NMDS2*10),
            size = 3) +
  new_scale_color()+
  geom_point(aes(fill = chemistry), size = 3, shape = 21) +
  stat_ellipse(aes(fill = chemistry, color = chemistry), alpha = 0.2, geom = "polygon", type="t") +
  theme_bw() 

sw_pw


#pc1 loadings barchart
NMDS1_loadings_bar <- v_scrs %>%
  ggplot(aes(x = reorder(Compound,NMDS1), y = abs(NMDS1))) +
  geom_bar(stat = "identity",
           aes(fill = Compound)) +
  geom_text(aes(label = Compound, angle=90, hjust=1)) +
  scale_fill_manual(values=wes_palette("FantasticFox1", type="continuous", 26)) + 
  labs(y = "Scaled NMDS1 Env.Fit Loadings",
       x = NULL) +
  labs(y = "Scaled NMDS1 Loadings",
       x = NULL) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        legend.position = "none") # theme(legend.position = "none",
#       axis.text.x = element_blank(),
#       axis.text.y = element_blank(),
#       axis.title = element_blank())

NMDS1_loadings_bar

#pc2 loadings barchart
NMDS2_loadings_bar <- v_scrs %>%
  ggplot(aes(x = abs(NMDS2), y = reorder(Compound,NMDS2))) +
  geom_bar(stat = "identity",
           aes(fill = Compound)) +
  geom_text(aes(label = Compound, hjust=1)) +
  labs(y = "Scaled NMDS2 Env.Fit Loadings",
       x = NULL) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        legend.position = "none")

NMDS2_loadings_bar

empty <- ggplot()+geom_point(aes(1,1), colour="white")+
  theme(axis.ticks=element_blank(), 
        panel.background=element_blank(), 
        axis.text.x=element_blank(), axis.text.y=element_blank(),           
        axis.title.x=element_blank(), axis.title.y=element_blank())


grid.arrange(NMDS1_loadings_bar, empty, plot, nrow = 2, ncol = 2, widths=c(4, 1), heights=c(1, 4))

(NMDS1_loadings_bar / plot)



#just with compartment colors
sw_pw<-data.scores %>%
  ggplot(aes(x=NMDS1, y=NMDS2)) +
  geom_point(mapping = aes(fill = chem$sw_pw), size=3, shape = 21)

plot=sw_pw+ggtitle('MAGs pore / surface per compartment')+theme(plot.title = element_text(hjust = 0.5))+scale_color_manual(values=c("#F8766D", "#00BFC4"))+
  stat_ellipse(aes(fill = chem$sw_pw, color=chem$sw_pw), alpha = 0.2, geom = "polygon", type="norm")+
  theme_bw()
plot