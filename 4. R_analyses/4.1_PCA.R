library(ggplot2)
library(tidyr)
library(tibble)
library(grid)
library(gridExtra)
library(ggforce)
library(ggfortify)
library(cluster)
library(vegan)
library(ggdark)
library(ggnewscale)

#Loading in our data
dat=read.table("4.1.1_32_sample_geochemistry_incl_FTICR_subset.txt", header = T, sep="\t")
rownames(dat)<-dat[,1]
dat=dat[-1]

metadat=dat[,1:3]
log_dat=log(dat[,4:29])

pca_res=prcomp(na.omit(log_dat), scale.=TRUE)
pcaresults_JRR = summary(pca_res)
scores_jrr <- as.data.frame(pcaresults_JRR$x) %>%
  rownames_to_column(var = "Compound") %>%
  separate(col = "Compound",
           into = c("compartment", "sample"),
           sep = "_")

compound_metadata = as.data.frame(colnames(log_dat))
compound_metadata$Compound = compound_metadata$`colnames(log_dat)`
compound_metadata$Group = compound_metadata$Compound
compound_metadata = compound_metadata[-1]

autoplot(pca_res, color=Compound)

###Calculating an ANOSIM and MRPP as well as Adonis PERMANOVA for timepoint and compartment.
metadat$sw_pw[which(metadat$sw_pw=="pw")] <- "1"
metadat$sw_pw[which(metadat$sw_pw=="sw")] <- "2"

mrpp(dat[,4:29], metadat$sw_pw, permutations=999, distance="bray")
anosim(dat[,4:29], metadat$sw_pw, permutations = 999)

########### Now getting plots in order
scores_jrr %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = compartment), size = 3, shape = 21) +
  geom_text(aes(label = sample), nudge_x = 1, size = 3) +
  stat_ellipse(aes(fill = compartment, color = compartment), alpha = 0.2, geom = "polygon") +
  labs(x = "PC1 (42.4%)",
       y = "PC2 (18.2%)") +
  theme_bw()

#euclidian distance funciton to get magnitude of vector arrows from 0,0
distance <- function(x, y, home = c(0,0)) {
  sqrt((x-home[1])^2 + (y-home[2])^2)
}


#loadings vectors
vl <- envfit(pca_res, log_dat, perm = 1000, p.max = 0.005)

vl
#create loadings df for plot arrows
v_scrs <- as.data.frame(scores(vl, display = "vectors"))
v_scrs <- cbind(v_scrs*ordiArrowMul(vl), Compound = rownames(v_scrs))

v_scrs <- v_scrs %>%
  mutate(dist_pc1pc2 = distance(PC1, PC2), #distance in space from center
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
         plot_pc1 = ifelse(top_c == "top", PC1, 0),
         plot_pc2 = ifelse(top_c == "top", PC2, 0),
         top_label = ifelse(top_c == "top", Compound, "")) %>%
  left_join(top_dist_compounds, by = c("Compound"))

library(wesanderson)
#Main PCA plot
pca_plot <- scores_jrr %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_segment(data = v_scrs,
               aes(x= 0, xend = PC1*10, 
                   y = 0, yend = PC2*10, 
                   alpha = alpha_factor,
                   color = Compound),
               arrow = arrow(length = unit(0.5, "cm")),
               inherit.aes = FALSE) +
  geom_text(data = v_scrs,
            aes(label = num_lab,
                x = PC1*10,
                y = PC2*10),
            size = 3) +
  scale_color_manual(values=wes_palette("FantasticFox1", type="continuous", 26)) +
  new_scale_color()+
  geom_point(aes(fill = compartment), size = 3, shape = 21) +
  stat_ellipse(aes(fill = compartment, color = compartment), alpha = 0.2, geom = "polygon", type='t') +
  theme_bw() + 
  dark_theme_gray()
pca_plot


##Basic PCA plot with just compartment and no loadings in ggplot:
#Main PCA plot
pca_plot_comp <- scores_jrr %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = compartment), size = 3, shape = 21) +
  stat_ellipse(aes(fill = compartment, color = compartment), alpha = 0.2, geom = "polygon")  +
  theme_bw()

pca_plot_comp



# pdf("figures/PCA/metab_pca_legend.pdf",
#     width = 5, height = 5)
# pca_plot
# dev.off()


limits_pc1 <- v_scrs %>%
  arrange(PC1) %>%
  pull(Compound)

#pc1 loadings barchart
pc1_loadings_bar <- v_scrs %>%
  ggplot(aes(x = reorder(Compound,PC1), y = abs(PC1))) +
  geom_bar(stat = "identity",
           aes(fill = Compound))  +
  scale_fill_manual(values=wes_palette("FantasticFox1", type="continuous", 26)) +
  geom_text(aes(label = Compound, angle=90, hjust=1.2)) +
  labs(y = "Scaled PC1 Loadings of Significant Compounds",
       x = NULL) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        legend.position = "none")
# theme(legend.position = "none",
#       axis.text.x = element_blank(),
#       axis.text.y = element_blank(),
#       axis.title = element_blank())

pc1_loadings_bar

#pc2 loadings barchart
pc2_loadings_bar <- v_scrs %>%
  ggplot(aes(x = abs(PC2), y = reorder(Compound,PC2))) +
  geom_bar(stat = "identity",
           aes(fill = Compound)) +
  geom_text(aes(label = Compound, hjust=-.1)) +
  labs(y = "Scaled PC2 Loadings of Significant Compounds",
       x = NULL) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        legend.position = "none")
# theme(legend.position = "none",
#       axis.text.x = element_blank(),
#       axis.text.y = element_blank(),
#       axis.title = element_blank())

pc2_loadings_bar

