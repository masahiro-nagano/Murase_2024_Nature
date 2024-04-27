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

colorscale = scale_fill_gradientn(
  colors = rev(brewer.pal(9, "YlGnBu")),
  values = c(0, exp(seq(-5, 0, length.out = 100))))

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
  data.frame() %>% na.omit()  -> all

all  %>% 
  mutate(KO=(KO1 + KO2)/2) %>%
  ggplot(aes(x=WT, y=KO)) +
  geom_hex(binwidth = c(1, 1)) + 
  colorscale +
  coord_fixed() +
  geom_abline(slope = 1, col="red") +
  geom_abline(slope = 1, intercept = 30, col="red") +
  geom_abline(slope = 1, intercept = -30, col="red") +
  theme_bw() +
  labs(x= "WT", y="KO")-> p
p
ggsave(filename = "f5_b.pdf", width=5, height = 5)