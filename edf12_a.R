library(tidyverse)
library(data.table)
library(rtracklayer)
library(R.utils)
library(RColorBrewer)
library(ggh4x)
library(pals)
library(viridisLite)
library(ggrepel)
library(ggsignif)
library(umap)
library(dbscan)

cats <- c("Active", "Bivalent", "Poised", "Silent", "All" )
cat_clrs <- setNames(tableau20(12)[c(7, 3, 5, 1, 9)], cats)

###############################################################################################################
# load open site annotation
# cluster information
# cluster 1: active, 2:bivalnet, 3: TE-TSS?,4: poised enhancers, 0: noise 
fread("./source/human_atac.clust.v2.csv") -> human_atac.clust

human_atac.clust %>% 
  mutate(name="Cluster") %>%
  ggplot(aes(x=UMAP1, y= UMAP2, col = factor(clu, levels=c("4", "3", "0", "1", "2")))) +
  geom_point(size=.1) +
  facet_wrap(. ~ name, ncol = 2) +
  scale_color_manual(values = c("#2CA02C", "#C9CACA", "#C9CACA", "#D62728", "#FF7F0E")) +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black', size = 20),
        axis.title=element_text(size=20),
        strip.background = element_rect(fill = '#C9CACA'),
        strip.text = element_text(color = 'black', size = 22),
        legend.position = 'none') +
  labs(col = 'cluster')-> p
p
ggsave(filename = "edf12_a1.png", width = 6, height = 8)

###############################################################################################################
# IP/Input heat map
human_atac.clust %>% 
  dplyr::select(UMAP1, UMAP2, H3K4me1, H3K27ac, H3K27me3, H3K4me3) %>%
  gather(-c(UMAP1, UMAP2), key=mod, value=val) %>%
  mutate(mod = factor(mod, levels = c('H3K4me1', 'H3K4me3', 'H3K27ac', 'H3K27me3'))) %>%
  ggplot(aes(x=UMAP1, y= UMAP2, col = val)) +
  geom_point(size=.1) +
  facet_wrap(. ~ mod, ncol = 4) +
  scale_color_gradient2(low='blue', high = 'red')+
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black', size = 20),
        axis.title=element_text(size=20),
        strip.background = element_rect(fill = '#C9CACA'),
        strip.text = element_text(color = 'black', size = 22),
        legend.position = 'right') +
  labs(col = 'log2(IP/input)') -> p
p
ggsave(filename = "edf12_a2.png", width = 12, height = 4)

human_atac.clust %>% 
  dplyr::select(clu, H3K4me1, H3K27ac, H3K27me3, H3K4me3) %>%
  gather(-c(clu), key=mod, value=val) %>%
  group_by(clu, mod) %>%
  summarise(m=mean(val)) %>%
  ggplot(aes(x=clu, y = mod, fill = m)) +
  geom_tile() +
  geom_text(aes(label = round(m,2)), color = "black", size = 3) +  
  theme(plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_rect(fill = "white", colour = "white"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-4, 4), 
                       name = "Enrichment", space = "Lab") +
  labs(x = 'Cluster', y = 'Mark') -> p
p
ggsave(filename = "edf12_a3.pdf", width = 5, height = 2.5)

human_atac.clust %>% 
  dplyr::select(clu) %>%
  dplyr::count(clu)

