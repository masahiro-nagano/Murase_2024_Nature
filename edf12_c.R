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
cat_clrs <- c("#953553", "#FFA500", "#088F8F", "#87CEEB", "#000000")

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
# load open site annotation
# cluster information
# cluster 1: active, 2:bivalnet, 3: TE-TSS?,4: poised enhancers, 0: noise 
fread("./source/human_atac.clust.v2.csv")

# active enhancers and promoters are defined by d4hPGCLCs -> 11693 regions
fread("./source/human_atac.clust.v2.csv") %>%
  filter(clu== "1" & type == "hd4PGCLC" ) %>%
  makeGRangesFromDataFrame() %>% GenomicRanges::reduce()  -> active

# bivalent regions are defined by d4hPGCLCs -> 16245 regions
fread("./source/human_atac.clust.v2.csv") %>%
  filter(clu== "2" & type == "hd4PGCLC"  ) %>%
  makeGRangesFromDataFrame() %>% GenomicRanges::reduce() -> bivalent

# poised regions are defined by d4hPGCLCs -> 115710 regions
fread("./source/human_atac.clust.v2.csv") %>%
  filter(clu== "4" & type == "hd4PGCLC" ) %>%
  makeGRangesFromDataFrame() %>% GenomicRanges::reduce()  -> poised

# promoter information
fread("./source/GRCh38p12_wo_pseudo_v2_chrxym_PARN_promoter_900us_400ds_hcp.bed", select = c(1:3, 6), col.names = c("chr", "start", "end", "Name")) %>% makeGRangesFromDataFrame(keep.extra.columns = T)  -> proms.gr


# define active prom/enh # 8598
active %>% .[overlapsAny(., proms.gr)] -> active.prom
active %>% .[!overlapsAny(., proms.gr)] -> active.enh

# define bivalent prom/enh # 4578
bivalent %>% .[overlapsAny(., proms.gr)] -> bivalent.prom
bivalent %>% .[!overlapsAny(., proms.gr)] -> bivalent.enh

# define poised prom/enh # 4164
poised %>% .[overlapsAny(., proms.gr)] -> poised.prom
poised %>% .[!overlapsAny(., proms.gr)] -> poised.enh

# define silent promoters by no overlap with open sites #17602
proms.gr %>% .[!overlapsAny(., active.prom)] %>% .[!overlapsAny(., bivalent.prom)] %>% .[!overlapsAny(., poised.prom)] -> silent.prom



# define silent enhancers by no overlap with open sites
fread("./source/human_atac.clust.v2.csv")%>%
  filter((clu== "1"  | clu== "2"  | clu== "4" ) & type != "hd4PGCLC"  ) %>%
  makeGRangesFromDataFrame() %>% 
  .[!overlapsAny(., proms.gr)] %>%
  .[!overlapsAny(., active.enh)] %>%
  .[!overlapsAny(., bivalent.enh)] %>%
  .[!overlapsAny(., poised.enh)] %>%
  GenomicRanges::reduce() -> silent.enh

##################################################################################################################
# c42
# promoter DNA methylation level
dss <- c(
  WT = "./source/c42_WT_methyl_Prom.bedGraph",
  KO1 = "./source/c42_KO1_methyl_Prom.bedGraph",
  KO2 = "./source/c42_KO2_methyl_Prom.bedGraph"
)

samps <- c( "WT", "KO1", "KO2")

dss  %>%
  mapply(function(f, n){
    fread(f, select=c(1:4), col.names = c("chr", "start", "end", n)) 
  }, ., names(.), SIMPLIFY=F) %>%
  Reduce(function(...) left_join(..., by =c(chr="chr", start="start", end = "end"), all.x = T), .) %>% 
  data.frame() %>% na.omit()  -> all

all %>% makeGRangesFromDataFrame(keep.extra.columns = T) -> all.mc


list(ActivePromoter = active.prom, BivalentPromoter=bivalent.prom, SilentPromoter = silent.prom, PoisedPromoter = poised.prom) %>%
  mapply(function(x, n){
    all.mc %>% .[overlapsAny(., x)] %>%
      data.frame() %>% mutate(cat=n)
  }, ., names(.), SIMPLIFY = F) %>%
  do.call(rbind, .) %>% remove_rownames() %>% mutate(KO = (KO1 + KO2)/2) -> cat.mc

cat.mc %>% 
  dplyr::transmute(seqnames, start, end,D = KO- WT, cat) %>%
  gather(-c(seqnames, start, end, cat), key= cell, value=mC) -> tmp

tmp %>% mutate(cat = "AllPromoter") %>%
  bind_rows(tmp, .) %>%
  mutate(cat = str_replace(cat, "Promoter", ""), reg="Promoter") %>%
  mutate(cat = factor(cat, levels = c("Active", "Bivalent", "Poised", "Silent", "All" ))) -> prom.mc.all


# Enhancer methylation level
dss <- c(
  WT = "./source/c42_WT_methyl_OpenSites.bedGraph",
  KO1 = "./source/c42_KO1_methyl_OpenSites.bedGraph",
  KO2 = "./source/c42_KO2_methyl_OpenSites.bedGraph"
)

samps <- c( "WT", "KO1", "KO2")

dss  %>%
  mapply(function(f, n){
    fread(f, select=c(1:4), col.names = c("chr", "start", "end", n)) 
  }, ., names(.), SIMPLIFY=F) %>%
  Reduce(function(...) left_join(..., by =c(chr="chr", start="start", end = "end"), all.x = T), .) %>% 
  data.frame() %>% na.omit()  -> all

all %>% makeGRangesFromDataFrame(keep.extra.columns = T) -> all.mc

list(ActiveEnhancer = active.enh, BivalentEnhancer=bivalent.enh, SilentEnhancer = silent.enh, PoisedEnhancer = poised.enh) %>%
  mapply(function(x, n){
    all.mc %>% .[overlapsAny(., x)] %>%
      data.frame() %>% mutate(cat=n)
  }, ., names(.), SIMPLIFY = F) %>%
  do.call(rbind, .) %>% remove_rownames() %>% mutate(KO = (KO1 + KO2)/2) -> cat.mc

cat.mc %>% 
  dplyr::transmute(seqnames, start, end,D = KO- WT, cat) %>%
  gather(-c(seqnames, start, end, cat), key= cell, value=mC) -> tmp

tmp %>% mutate(cat = "AllEnhancer") %>%
  bind_rows(tmp, .) %>%
  mutate(cat = str_replace(cat, "Enhancer", ""), reg="Enhancer") %>%
  mutate(cat = factor(cat, levels = c("Active", "Bivalent", "Poised", "Silent", "All" ))) -> enh.mc.all

bind_rows(prom.mc.all, enh.mc.all) %>%
  filter(!cat %in% c("All")) %>%
  mutate(reg=factor(reg, levels = c("Promoter", "Enhancer"))) %>%
  ggplot(aes(x=cat, y= mC*100)) +
  geom_violin(aes(fill=cat), width =1) +
  geom_boxplot(width =.1) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  labs(x="Category", y= "mC (%) (KO - WT)") +
  scale_fill_manual(values = cat_clrs) +
  facet_wrap(.~ reg, nrow = 2, strip.position = "right") +
  geom_hline(yintercept = 0, col= 'grey80') +
  labs(x="") +
  coord_cartesian(ylim = c(-5, 25)) +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        panel.grid = element_blank(),
        axis.title=element_text(size=20),
        axis.text = element_text(color = 'black', size = 16),
        strip.background = element_rect(fill = '#C9CACA'),
        strip.text = element_text(color = 'black', size = 22),
        legend.position = 'none') -> p
p

ggsave(filename = "edf12_c.pdf", width = 5, height = 5)
