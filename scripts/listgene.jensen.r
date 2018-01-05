library(readr)
library(dplyr)
benchmark <- read_tsv("./human_compartment_benchmark.tsv",col_names = F)

#colnames(benchmark) <- c("ENSP","cellular_location")

experiment <- read_tsv("./human_compartment_experiments_full.tsv",col_names = F)
knowledge <- read_tsv("./human_compartment_knowledge_full .tsv",col_names = F)
predictions <- read_tsv("./human_compartment_predictions_full.tsv",col_names = F)
textmining <- read_tsv("./human_compartment_textmining_full.tsv",col_names = F)

experiment %>%
  filter(X4=="Plasma membrane",X7>=mean(X7))->mem.experiment
predictions %>%
  filter(X4=="Plasma membrane",X7>=mean(X7))->mem.predicted
knowledge%>%
  filter(X4=="Plasma membrane",X7>=mean(X7))->mem.knowledge
textmining%>%
  filter(X4=="Plasma membrane",X6>=mean(X6))->mem.textmining

symbl <- c(mem.experiment$X2,mem.knowledge$X2,mem.predicted$X2,mem.textmining$X2)
symbl <- unique(symbl)

genelist <- setdiff(symbl,membrane.protein.hpa$Gene)

write_csv(as.data.frame(genelist),"./diffgene.hpa_jensen.csv")


#explore the diffgene betwen hpa and jemsen
diffgene <- read_tsv("./diffgene.tsv")
colnames(diffgene) <- colnames(pro_reliable)
diffgene%>%
  filter(grepl("Plasma membrane",`Subcellular location`)) %>%
  rbind(filter(diffgene,is.na(`Subcellular location`)))%>%
  filter(Evidence=="Evidence at protein level",!is.na(Antibody))->diffgene.jsen

genelist.membrane <- membrane.protein.hpa%>%
  rbind(diffgene.jsen)

write_csv(genelist.membrane,"./ImmuneTCGA_geneinfo.csv")


