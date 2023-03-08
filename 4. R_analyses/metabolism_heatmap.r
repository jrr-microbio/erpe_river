library(reshape2)
library(tidyverse)
library(reticulate)
library(forcats)
#path_to_python="/Users/josue.rodriguez/Library/r-miniconda/envs/r-reticulate/bin/python"
#use_python(path_to_python)
#py_install('pandas')
#np=import('numpy')
#pd=import('pandas')

metabolism_matrix=read.csv("/Volumes/GoogleDrive/My Drive/University/Wrighton Lab; PhD/ERPE_Timeseries_viral_communities/microbial_work/metabolism_per_compartment/simplified_dram_product_for_dotplot.csv",header=T)

#remove the pw/sw indicator status. I'll put this back later.
metabolism_matrix_short=metabolism_matrix[,-c(2:5)]
#Now convert this matrix into long format for ggplot2
metabolism_long = melt(data = metabolism_matrix_short)
#panda_metabolism_long=pandas.melt(metabolism_matrix_short, ignore_index=FALSE)

#index the factors that we want to join back in after the melt.
factors.join=metabolism_matrix %>%
  select("genome_id", "sumpaths.index", "indicator", "tax.index", "p.taxonomy", "taxonomy")

metabolism.long.joined = metabolism_long %>%
  left_join(factors.join)

#r was treating the columns as numeric, so I convert them to factors since it just doesn't like numbers apparently. But this fixes that.
#metabolism.long.joined<- transform(metabolism.long.joined, genome_id= reorder(genome_id, sumpaths.index))

metabolism.long.joined<- transform(metabolism.long.joined, genome_id= reorder(genome_id, tax.index))

#making a color palette because that default was very hard to look at.
pokecolorpal=c("#C5B420", "#F6E673","#FFFF00","#416A10","#7BC552", "#00FF00","#B46220","#E68329", "#392910", "#5A5A5A", "#C6DDEE", "#19ABB4", "#3C30AD", "#06525A", "#BD4141", "#FF7BAC", "#AAA6F0","#FF00FF")

metabolism.long.joined %>%
  mutate_if(is.numeric, ~1 * (.>0)) %>%
  filter(value > 0, 'Presence / Absence' > 1) %>%
  ggplot(aes(x=fct_reorder(genome_id,tax.index, .desc = TRUE), y = reorder(variable, desc(variable), size = value))) + 
  geom_point(aes(colour=factor(p.taxonomy))) +
  scale_colour_manual(values=pokecolorpal) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=0.5))

