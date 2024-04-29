library(tidyverse)
library(data.table)
library(patchwork)

load("./source/diff.ex.cat.rda")

fread("./source/expression_T1KO_DEG_list2284genes_20230221.txt") -> tet.d

cats <- c("Active", "Bivalent", "Poised", "Silent", "All", "Genome-wide", "Early ER", "Late ER", "Imprint DMR", "ER genes")
cat_clrs <- setNames(c("#953553", "#FFA500", "#088F8F", "#87CEEB", "#FFFF00", tableau20(20)[c(7, 3, 5, 1)], "#000000"), cats)


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



diff.ex.cat %>%
  dplyr::count(category) %>%
  mutate(class = 'background') -> bg

diff.ex.cat %>%
  filter(Name %in% c(tet.degs$c12_up ))%>%
  dplyr::count(category) %>%
  mutate(class = 'DEG')-> c12.up

diff.ex.cat %>%
  filter(Name %in% c(tet.degs$c12_down))%>%
  dplyr::count(category) %>%
  mutate(class = 'DEG')-> c12.down


bind_rows(c12.up, bg) %>%
  group_by(class) %>%
  mutate(n=n/sum(n)) %>%
  spread(key = class, value = n) %>%
  mutate(r= DEG/background) %>%
  mutate(category = factor(category, levels = cats)) %>%
  ggplot(aes(x= category, y = r, col = category)) +
  geom_point(size=6) +
  theme_bw() +
  scale_color_manual(values = cat_clrs) +
  coord_cartesian(ylim = c(0, 4)) +
  labs(y="Odds Ratio", x="", col= "Cluster", title = 'c12 Up') +
  geom_hline(yintercept = 1, col = "red") +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        panel.grid = element_blank(),
        axis.title=element_text(size=20),
        axis.text = element_text(color = 'black', size = 16),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, color = 'black', size = 16),
        strip.background = element_rect(fill = '#C9CACA'),
        strip.text = element_text(color = 'black', size = 22),
        legend.position = 'none')-> p1
p1

early <- fread('./source/231202_core_ERgene_early.txt', col.names = 'Name', header = F) %>% .$Name
late <- fread('./source/231202_core_ERgene_late.txt', col.names = 'Name', header = F) %>% .$Name

all_ER <- c(early, late)

diff.ex.cat %>%
  mutate(category = case_when(Name %in% c(early, late) ~ 'ER genes',
                              TRUE ~ category)) -> diff.ex.cat.er

diff.ex.cat.er %>%
  dplyr::count(category) %>%
  mutate(class = 'background') -> bg.er

diff.ex.cat.er %>%
  filter(Name %in% c(tet.degs$c12_up))%>%
  dplyr::count(category) %>%
  mutate(class = 'DEG') -> c12.up.er

diff.ex.cat.er %>%
  filter(Name %in% c(tet.degs$c12_down))%>%
  dplyr::count(category) %>%
  mutate(class = 'DEG') -> c12.down.er

# Down regulted genes
bind_rows(c12.down.er, bg.er) %>%
  group_by(class) %>%
  mutate(n=n/sum(n)) %>%
  mutate(category = str_replace(category, 'Early', 'Early ER'), category = str_replace(category, 'Late', 'Late ER')) %>%
  # filter(category %in% c('Early', 'Late')) %>%
  mutate(category = factor(category, levels = cats)) %>%
  spread(key = class, value = n) %>%
  mutate(r= DEG/background) %>%
  na.omit() %>%
  ggplot(aes(x= category, y = r, col = category)) +
  geom_point(size=6) +
  scale_color_manual(values = cat_clrs) +
  theme_bw() +
  labs(y="Odds Ratio", x="", col= "Cluster", title = 'c12 Down') +
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

wrap_plots(list(p1, p2), nrow = 1, ncol = 2, width = c(2.5, 2.5)) -> p
p

ggsave(filename = "f5_h.pdf", width = 4, height = 4.5, p)

# check exact numbers
bind_rows(c12.down.er, bg.er) %>%
  group_by(class)

bind_rows(c12.up, bg) %>%
  group_by(class)
