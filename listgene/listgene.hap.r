library(readr)

#in protein classes
protein_class_predicted <- read_tsv("./HPA/protein_class_Predicted.tsv.gz")

protein_class_predicted %>%
  filter(is.na(`Subcellular location`))-> protein.1

protein_class_predicted %>%
  filter(`Subcellular location`=="Plasma membrane") ->protein.membrn


protein_class_predicted %>%
  filter(`Subcellular location`=="Vesicles<br>Plasma membrane<br>Mitochondria") %>%
  rbind(protein.plmembrn)-> protein.plmembrn

protein.1 %>%
  rbind(protein.celljunction) %>%
  rbind(protein.membrn) %>%
  rbind(protein.plmembrn)->protein.m

write_tsv(protein.m,"./genelist_hap.tsv")

 
