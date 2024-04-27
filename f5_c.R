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

dss <- c(
  WT = "./source/genomewide_2kb_GRCh38_CpG_868_T1KO_WT_c12_GRCh38p12chrxymPARN.bed",
  KO1 = "./source/genomewide_2kb_GRCh38_CpG_868_T1KO_KO1_c12_GRCh38p12chrxymPARN.bed",
  KO2 = "./source/genomewide_2kb_GRCh38_CpG_868_T1KO_KO2_c12_GRCh38p12chrxymPARN.bed"
)

# 2kb analysis in c12
dss  %>%
  mapply(function(f, n){
    fread(f, col.names = c("chr", "start", "end","depth", n)) %>%
      filter(depth>=4) %>%
      dplyr::select(-depth)
  }, ., names(.), SIMPLIFY=F) %>%
  Reduce(function(...) left_join(..., by =c(chr="chr", start="start", end = "end"), all.x = T), .) %>% 
  data.frame() %>% na.omit()  -> all

samps <- c( "WT", "KO1", "KO2")

dss  %>%
  mapply(function(f, n){
    fread(f, col.names = c("chr", "start", "end","depth", n)) %>%
      filter(depth>=4) %>%
      dplyr::select(-depth)
  }, ., names(.), SIMPLIFY=F) %>%
  Reduce(function(...) left_join(..., by =c(chr="chr", start="start", end = "end"), all.x = T), .) %>% 
  data.frame() %>% na.omit()  -> all
all %>% filter((KO1 - WT >30) & (KO2 - WT >30)) -> diff.mc.1kb
diff.mc.1kb %>% makeGRangesFromDataFrame() -> diff.mc.1kb.gr

diff.mc.1kb.gr %>% GRangesList() -> qSets

# Set bins of random genome background
gn <- fread('./source/hg38.chrom.sizes', header = F) 
tiles <- Seqinfo(seqnames = gn$V1, seqlengths = gn$V2, genome = "hg38") %>%
  tileGenome(tilewidth = 1000, cut.last.tile.in.chrom = T)

qSets %>% unlist %>% length -> l
set.seed(1)
sample(length(tiles), l) -> ids
tiles[ids] %>%makeGRangesFromDataFrame() -> uSet

# test regions
rSets <- list.files("./source/ensembl", pattern = ".bed") %>%
  setNames(.,sub(".bed", "", .)) %>%
  lapply(., function(x){
    fread(sprintf("./source/ensembl/%s", x), col.names = c("chr", "start", "end"))%>%
      makeGRangesFromDataFrame()
  }) %>% 
  GRangesList

# length of sets
qLen <- unlist(lapply((qSets),length))
uLen <- length(uSet)
ol.min <- 1
ol.u <- countOverlaps(rSets, uSet, minoverlap = ol.min)
ol.q <- lapply(qSets, function(x) {
  countOverlaps(rSets, x, minoverlap = ol.min)
}) %>% 
  do.call(cbind, .) %>%
  t %>%
  reshape2::melt(variable.factor = F) %>%
  dplyr::rename(userSet = Var1,
                region = Var2,
                support = value) %>%
  mutate(userSet = as.character(userSet),
         region = as.character(region),
         b = ol.u[match(region, names(ol.u))],
         c = qLen - support,
         d = uLen - b) %>%
  cbind(., apply(., 1, function(x) {
    fisher.test(matrix(as.numeric(x[3:6]), 2, 2),
                alternative = "greater")[c("p.value", "estimate")] %>%
      unlist
  }) %>%
    t %>%
    `colnames<-`(c('pValue','oddsRatio'))
  ) %>%
  left_join(sapply(rSets, length) %>%
              data.frame(size = .) %>%
              rownames_to_column("region"),
            by = "region")
ol.q %>%
  filter(pValue < 0.05) %>%
  mutate(q=-log10(pValue)) %>%
  ggplot(aes(x=oddsRatio, y=q, label=region))  +
  geom_point() +
  geom_label_repel(min.segment.length = 0.01, size=4, 
                   box.padding = 1.5) +
  theme_classic() +
  labs(title="ENCODE annotation", y="q value")-> p
p
ggsave(filename = "f5_c.pdf", width=4.5, height = 4.5)

