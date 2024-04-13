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


cats <- c("Active", "Bivalent", "Poised", "Silent", "All", "Genome-wide", "Early ER", "Late ER", "Imprint DMR" )
cat_clrs <- setNames(tableau20(20)[c(7, 3, 5, 1, 9, 11, 13, 15, 17)], cats)

TSS_anno_900_400 <- read.table("./source/TET1_TSS_anno_900_400.txt", sep="\t", header=T,row.names=1)

length(unique(TSS_anno_900_400$symbol)) #6217 unique genes are targeted

TSS_anno_900_400$symbol %>% unique -> tet.trg

fread('./source/promoter_classification.tsv') %>%
  filter(!category %in% c('Early', 'Late')) -> all

# read DEG
# DEG analysis
fread("./source/expression_T1KO_DEG_list2284genes_20230221.txt") -> tet.d

tet.degs <- list(
  iPSC_up =  tet.d %>% filter(diff_iPSC == 1) %>% .$Name,
  iPSC_down =  tet.d %>% filter(diff_iPSC == -1) %>% .$Name,
  iMeLC_up =  tet.d %>% filter(diff_iMeLC == 1) %>% .$Name,
  iMeLC_down =  tet.d %>% filter(diff_iMeLC == -1) %>% .$Name,
  d6_up =  tet.d %>% filter(diff_d6 == 1) %>% .$Name,
  d6_down =  tet.d %>% filter(diff_d6 == -1) %>% .$Name,
  c12_up =  tet.d %>% filter(diff_c12 == 1) %>% .$Name,
  c12_down =  tet.d %>% filter(diff_c12 == -1) %>% .$Name,
  c32_up =  tet.d %>% filter(diff_c32 == 1) %>% .$Name,
  c32_down =  tet.d %>% filter(diff_c32 == -1) %>% .$Name,
  c42_up =  tet.d %>% filter(diff_c42 == 1) %>% .$Name,
  c42_down =  tet.d %>% filter(diff_c42 == -1) %>% .$Name) 

all %>%
  dplyr::count(category) %>%
  mutate(n=n/sum(n)) -> all.dist

all %>% filter(Name %in% tet.trg) %>%
  dplyr::count(category) %>%
  mutate(n=n/sum(n)) -> tet.trg.dist

all %>%
  mutate(category = str_replace(category, 'Early', 'Early ER'), category = str_replace(category, 'Late', 'Late ER')) %>%
  mutate(tet.ol = case_when(Name %in% tet.trg ~ 'TET target',
                            TRUE ~ 'Others')) %>%
  mutate(c12_cat = case_when(Name %in% tet.degs$c12_up ~ 'Upregulated',
                             Name %in% tet.degs$c12_down ~ 'Downregulated',
                             TRUE ~ 'NonDEG')) %>%
  filter(tet.ol != 'Others') %>%
  dplyr::count(category, tet.ol, c12_cat) %>%
  group_by(tet.ol, c12_cat) %>%
  mutate(r= n/sum(n)) %>%
  dplyr::transmute(category, tet.ol,number = n,  r, c12_cat) %>%
  left_join(., tet.trg.dist, by = c(category = 'category')) %>%
  mutate(odds = r/n) %>%
  filter(c12_cat == 'Upregulated') %>%
  mutate(category = factor(category, levels = cats)) %>%
  ggplot(aes(x= category, y = odds, col = category)) +
  geom_point(size = 12) +
  geom_hline(yintercept = 1, col = 'grey') +
  scale_color_manual(values = cat_clrs) +
  theme_bw() +
  labs(y="Odds Ratio", x="", title = 'TET Target & Upregulated genes') +
  geom_hline(yintercept = 1, col = "red") +
  #coord_cartesian(ylim = c(0, 4)) +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        panel.grid = element_blank(),
        axis.title=element_text(size=20),
        axis.text = element_text(color = 'black', size = 16),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, color = 'black', size = 16),
        strip.background = element_rect(fill = '#C9CACA'),
        strip.text = element_text(color = 'black', size = 22),
        legend.position = 'none')-> p2
p2

ggsave(p2, filename = "f6_h.pdf", width = 4.3, height = 6.07)

# For labeling
all %>%
  mutate(category = str_replace(category, 'Early', 'Early ER'), category = str_replace(category, 'Late', 'Late ER')) %>%
  mutate(tet.ol = case_when(Name %in% tet.trg ~ 'TET target',
                            TRUE ~ 'Others')) %>%
  mutate(c12_cat = case_when(Name %in% tet.degs$c12_up ~ 'Upregulated',
                             Name %in% tet.degs$c12_down ~ 'Downregulated',
                             TRUE ~ 'NonDEG')) %>%
  filter(tet.ol != 'Others') %>%
  dplyr::count(category, tet.ol, c12_cat) %>%
  group_by(tet.ol, c12_cat) %>%
  filter(c12_cat == 'Upregulated') 


# Also show in c12, TET target is enriched
all %>%
  mutate(category = str_replace(category, 'Early', 'Early ER'), category = str_replace(category, 'Late', 'Late ER')) %>%
  mutate(tet.ol = case_when(Name %in% tet.trg ~ 'TET target',
                            TRUE ~ 'Others')) %>%
  mutate(c12_cat = case_when(Name %in% tet.degs$c12_up ~ 'Upregulated',
                             Name %in% tet.degs$c12_down ~ 'Downregulated',
                             TRUE ~ 'NonDEG')) %>%
  dplyr::count(tet.ol, c12_cat) %>%
  group_by(tet.ol) %>%
  mutate(r= n/sum(n)) %>%
  ungroup() %>%
  dplyr::select(-n) %>% 
  spread(key = tet.ol, value = r) %>%
  mutate(odds =  `TET target`/Others)
