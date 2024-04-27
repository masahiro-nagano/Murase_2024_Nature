# check 2kb bin DNAme lebel
library(tidyverse)
library(data.table)
library(rtracklayer)
library(ggsignif)

# 2kb analysis
# c12
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
  data.frame() %>% na.omit()  -> all

all %>% makeGRangesFromDataFrame(keep.extra.columns = T) -> all.mc


fread("./source/ensembl/gene.bed", col.names=c("chr", "start", "end")) %>% makeGRangesFromDataFrame() -> genes

fread("./source/ensembl/intergenic.bed", col.names=c("chr", "start", "end")) %>% makeGRangesFromDataFrame() -> intergenic

list.files("./source/peaks2") %>%
  lapply(., function(x){
    fread(sprintf("./source/peaks2/%s", x), select=c(1:3), col.names=c("chr", "start", "end"))
  }) %>% do.call(rbind, .) %>% makeGRangesFromDataFrame() %>% GenomicRanges::reduce() -> regs

regs %>% .[overlapsAny(., genes)]  -> regs.genes
regs %>% .[overlapsAny(., intergenic)]  -> regs.intergenic
intergenic %>% GenomicRanges::setdiff(., regs) -> no.regs.intergenic
genes %>% GenomicRanges::setdiff(., regs) -> no.regs.genes

list(Regs.genes=regs.genes, Genes.no.reg= no.regs.genes, Regs.intergenic = regs.intergenic, Intergenic.no.reg = no.regs.intergenic) %>%
  mapply(function(x, n){
    all.mc %>% .[overlapsAny(., x)] %>%
      data.frame() %>% mutate(cat=n)
  }, ., names(.), SIMPLIFY = F) %>%
  do.call(rbind, .) %>% remove_rownames() -> cat.mc

cat.mc %>% 
  dplyr::select(seqnames, start, end, WT, KO1, KO2, cat) %>%
  gather(-c(seqnames, start, end, cat), key= cell, value=mC) %>%
  mutate(cell= factor(cell, levels = c( "WT", "KO1", "KO2"))) %>%
  group_by(cat, cell) %>%
  summarise(Median=round(median(mC), 2), number=n(), Mean = round(mean(mC), 2)) -> stats


samps <- c("WT", "KO1", "KO2" )
samp_clrs <- setNames(tableau20(12)[c(1, 7, 3)], samps)

cat.mc %>%
  mutate(diff = (KO1 +KO2)/2 - WT) %>%
  mutate(category = case_when(cat == 'Genes.no.reg' ~ 'Gene-NRE',
                              cat == "Intergenic.no.reg" ~ 'Intergene-NRE',
                              cat == 'Regs.genes' ~ 'Gene-RE',
                              cat == "Regs.intergenic" ~ 'Intergene-RE')) %>%
  mutate(category  = factor(category, levels = c('Gene-NRE', 'Intergene-NRE', 'Gene-RE', 'Intergene-RE'))) %>%
  ggplot(aes(x=category, y = diff)) +
  geom_violin(aes(fill=category)) +
  geom_boxplot(width =.1) +
  geom_hline(yintercept = 0, col = 'yellow') +
  geom_signif(comparisons = list(c('Intergene-NRE', 'Gene-NRE'), c('Intergene-NRE','Gene-RE'), c('Intergene-NRE', 'Intergene-RE')),
              map_signif_level = T, step_increase = .1)  +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        panel.grid = element_blank(),
        axis.title=element_text(size=20),
        axis.text = element_text(color = 'black', size = 16),
        axis.text.x = element_text(angle = 90, vjust=.5),
        legend.position = 'none') +
  coord_cartesian(ylim = c(-20, 20)) +
  labs(x='', y= 'DNAme difference (KO-WT)') -> p
p
ggsave(filename = "f5_f.pdf", width=5, height = 6.75)

# number
cat.mc %>%
  mutate(diff = (KO1 +KO2)/2 - WT) %>%
  mutate(category = case_when(cat == 'Genes.no.reg' ~ 'Gene-NRE',
                              cat == "Intergenic.no.reg" ~ 'Intergene-NRE',
                              cat == 'Regs.genes' ~ 'Gene-RE',
                              cat == "Regs.intergenic" ~ 'Intergene-RE')) %>%
  mutate(category  = factor(category, levels = c('Gene-NRE', 'Intergene-NRE', 'Gene-RE', 'Intergene-RE'))) %>%
  dplyr::count(category)