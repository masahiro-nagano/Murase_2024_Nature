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

caluculation <- function(data, bg, n){
  set.seed(1)
  c(1:n) %>% lapply(., function(n){
    idx <- sample(size=nrow(data), c(1:nrow(bg)))
    bg[idx,] %>% .$D -> ctrl
    data %>% .$mC-> trg
    t.test(trg, ctrl) -> res
    data.frame(diff=mean(trg) - mean(ctrl), P=res$p.value)
  }) %>% do.call(rbind, .) %>%
    summarise(d=median(diff), p=median(P))
}
###############################################################################################################

load("./source/c12.enhancer.diff.rda")
load( "./source/c12.promoter.diff.rda")
c12.enhancer.diff
bind_rows(c12.enhancer.diff %>% mutate(category = str_replace( category, "Enhancer", "")), c12.promoter.diff) %>%
  filter(cell=="c12_diff") %>%
  mutate(class=factor(class, levels = c('Promoter', 'Enhancer'))) %>%
  ungroup() %>%
  ggplot(aes(x= category, y=value)) +
  geom_violin(aes(fill=category)) +
  geom_boxplot(width =.1) +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  facet_wrap(. ~ class, nrow = 2, strip.position = "right") +
  labs(x="", y= "Log2FC(KO/WT)") +
  theme_bw() +
  scale_fill_manual(values = cat_clrs) +
  geom_hline(yintercept = 0, col= 'grey80') +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        panel.grid = element_blank(),
        axis.title=element_text(size=20),
        axis.text = element_text(color = 'black', size = 16),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        strip.background = element_rect(fill = '#C9CACA'),
        strip.text = element_text(color = 'black', size = 22),
        legend.position = 'none') -> p
p
ggsave(filename = "f5_g.pdf", width = 4, height = 4.8)