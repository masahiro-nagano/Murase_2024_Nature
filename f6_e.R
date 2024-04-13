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


cats <- c("Genome-wide", "Early", "Late", "Imprint DMR" )
cat_clrs <- setNames(tableau20(17)[c(11, 13, 15, 17)], cats)

# Definition of early, late core ER gene
early <- fread('./source/231202_core_ERgene_early.txt', col.names = 'Name', header = F) %>% .$Name
late <- fread('./source/231202_core_ERgene_late.txt', col.names = 'Name', header = F) %>% .$Name

# promoter information
fread("./source/GRCh38p12_wo_pseudo_v2_chrxym_PARN_promoter_900us_400ds_hcp.bed", select = c(1:3, 6), col.names = c("chr", "start", "end", "Name")) %>% makeGRangesFromDataFrame(keep.extra.columns = T)  -> proms.gr

proms.gr %>% .[.$Name %in% early]-> early.prom
proms.gr %>% .[.$Name %in% late]-> late.prom


# c12
# promoter DNA methylation level
dss <- c(
  WT = "./source/c12_WT_methyl_Prom.bedGraph",
  KO1 = "./source/c12_KO1_methyl_Prom.bedGraph",
  KO2 = "./source/c12_KO2_methyl_Prom.bedGraph"
)

samps <- c( "WT", "KO1", "KO2")

dss  %>%
  mapply(function(f, n){
    fread(f, select=c(1:4), col.names = c("chr", "start", "end", n)) 
  }, ., names(.), SIMPLIFY=F) %>%
  Reduce(function(...) left_join(..., by =c(chr="chr", start="start", end = "end"), all.x = T), .) %>% 
  data.frame() %>% na.omit()  -> all

all %>% makeGRangesFromDataFrame(keep.extra.columns = T) -> all.mc


list(Early = early.prom, Late = late.prom) %>%
  mapply(function(x, n){
    all.mc %>% .[overlapsAny(., x)] %>%
      data.frame() %>% mutate(cat=n)
  }, ., names(.), SIMPLIFY = F) %>%
  do.call(rbind, .) %>% remove_rownames() %>% mutate(KO = (KO1 + KO2)/2) %>%
  mutate_at(vars(c('WT', 'KO1', 'KO2', 'KO')), function(x){x*100})-> cat.mc

#################
# genome-wide
# 2kb analysis in c12
dss <- c(
  WT = "./source/genomewide_2kb_GRCh38_CpG_868_T1KO_WT_c12_GRCh38p12chrxymPARN.bed",
  KO1 = "./source/genomewide_2kb_GRCh38_CpG_868_T1KO_KO1_c12_GRCh38p12chrxymPARN.bed",
  KO2 = "./source/genomewide_2kb_GRCh38_CpG_868_T1KO_KO2_c12_GRCh38p12chrxymPARN.bed"
)

samps <- c( "WT", "KO1", "KO2")

dss  %>%
  mapply(function(f, n){
    fread(f, col.names = c("chr", "start", "end","depth", n)) %>%
      filter(depth>=4) %>%
      dplyr::select(-depth)
  }, ., names(.), SIMPLIFY=F) %>%
  Reduce(function(...) left_join(..., by =c(chr="chr", start="start", end = "end"), all.x = T), .) %>% 
  data.frame() %>% na.omit()   %>% makeGRangesFromDataFrame(keep.extra.columns = T)  %>%
  data.frame() %>% mutate(cat='Genome-wide') %>% mutate(KO = (KO1 + KO2)/2) -> cat.mc.genomewide


#################
# Imprint
dss <- c(
  WT = "./source/hg38_imprint_DMRs_Court_liftover_CpG_868_T1KO_parent_c12_GRCh38p12chrxymPARN.bed",
  KO1 = "./source/hg38_imprint_DMRs_Court_liftover_CpG_868_T1KO_KO1_c12_GRCh38p12chrxymPARN.bed",
  KO2 = "./source/hg38_imprint_DMRs_Court_liftover_CpG_868_T1KO_KO2_c12_GRCh38p12chrxymPARN.bed"
)

list_imprint<-c("Court_imprint_DMR_known17","Court_imprint_DMR_known12","Court_imprint_DMR_known14","Court_imprint_DMR_known29","Court_imprint_DMR_known6","Court_imprint_DMR_placenta_7","Court_imprint_DMR_near_known7","Court_imprint_DMR_placenta_16","Court_imprint_DMR_placenta_1","Court_imprint_DMR_novel3","Court_imprint_DMR_known16","Court_imprint_DMR_placenta_2","Court_imprint_DMR_placenta_17","Court_imprint_DMR_placenta_6","Court_imprint_DMR_placenta_15","Court_imprint_DMR_placenta_8","Court_imprint_DMR_placenta_10","Court_imprint_DMR_placenta_4","Court_imprint_DMR_placenta_13","Court_imprint_DMR_placenta_14","Court_imprint_DMR_placenta_5","Court_imprint_DMR_placenta_11","Court_imprint_DMR_placenta_9","Court_imprint_DMR_placenta_12","Court_imprint_DMR_placenta_3","Court_imprint_DMR_known34","Court_imprint_DMR_known35","Court_imprint_DMR_known5","Court_imprint_DMR_known26","Court_imprint_DMR_known11","Court_imprint_DMR_known31","Court_imprint_DMR_near_known1","Court_imprint_DMR_known1","Court_imprint_DMR_novel2","Court_imprint_DMR_known8","Court_imprint_DMR_known4","Court_imprint_DMR_known7","Court_imprint_DMR_known9","Court_imprint_DMR_known3","Court_imprint_DMR_novel4","Court_imprint_DMR_known10","Court_imprint_DMR_known15","Court_imprint_DMR_known25","Court_imprint_DMR_novel1","Court_imprint_DMR_near_known8","Court_imprint_DMR_near_known5","Court_imprint_DMR_known30","Court_imprint_DMR_known32","Court_imprint_DMR_novel6","Court_imprint_DMR_known28")


dss %>%
  mapply(function(f, n){
    fread(f, select=c(1:4, 8), col.names = c("chr", "start", "end", "name", n)) 
  }, ., names(.), SIMPLIFY=F) %>%
  Reduce(function(...) left_join(..., by =c(chr="chr", start="start", end = "end", name = "name"), all.x = T), .) %>% 
  data.frame() %>% na.omit()  %>% makeGRangesFromDataFrame(keep.extra.columns = T) %>%
  data.frame() %>%
  filter(name %in% list_imprint) %>% dplyr::select(-name) %>%
  mutate(cat='Imprint DMR') %>%
  mutate(KO = (KO1 + KO2)/2) -> cat.mc.imprint

cats <- c("Genome-wide", "Early ER", "Late ER", "Imprint DMR" )
cat_clrs <- setNames(tableau20(17)[c(11, 13, 15, 17)], cats)

bind_rows(cat.mc, cat.mc.genomewide, cat.mc.imprint) %>%
  mutate(cat = str_replace(cat, 'Early', 'Early ER'), cat = str_replace(cat, 'Late', 'Late ER')) %>%
  mutate(cat = factor(cat, levels = cats)) %>%
  dplyr::transmute(seqnames, start, end, D = KO- WT, cat) %>%
  gather(-c(seqnames, start, end, cat), key= cell, value=mC)%>%
  ggplot(aes(x=cat, y= mC)) +
  geom_violin(aes(fill=cat), width =1) +
  geom_boxplot(width =.1) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  labs(x="Category", y= "mC (%) (KO - WT)") +
  scale_fill_manual(values = cat_clrs) +
  geom_hline(yintercept = 0, col= 'grey80') +
  labs(x="") +
  coord_cartesian(ylim = c(-5, 30)) +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        panel.grid = element_blank(),
        axis.title=element_text(size=20),
        axis.text = element_text(color = 'black', size = 16),
        strip.background = element_rect(fill = '#C9CACA'),
        strip.text = element_text(color = 'black', size = 22),
        legend.position = 'none') -> p
p

ggsave(filename = "f6_e.pdf", width = 4, height = 4.44)
ggsave(filename = "f6_e.png", width = 4, height = 4.44)


bind_rows(cat.mc, cat.mc.genomewide, cat.mc.imprint) %>%
  dplyr::transmute(seqnames, start, end, D = KO- WT, cat) %>%
  gather(-c(seqnames, start, end, cat), key= cell, value=mC) %>%
  dplyr::count(cat)