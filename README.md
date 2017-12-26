# ImmuneTCGA
Immune cell surface gene list in TCGA

## Papers
[Integrating Proteomics and Transcriptomics for Systematic Combinatorial Chimeric Antigen Receptor Therapy of AML](https://www.sciencedirect.com/science/article/pii/S1535610817304087)


> fuxins prepare to wrok.

readRDS("./res/data/TCGA_expr/pancan33_expr.rds",refhook=NULL)

## Resources
[CAR T Cells: Engineering Patients’ Immune Cells to Treat Their Cancers](https://www.cancer.gov/about-cancer/treatment/research/car-t-cells)

> fuxins prepare to wrok.

library(dplyr)

#biomart_export.txt is downoaded from Ensemble:Biomart
biomart_export <- read.csv("./biomart_export.txt")
colnames(biomart_export) <- c("GeneID","Gene.name","ProteinID","CCDS.ID","UniSwi.ID")

#biomart_export$UniSwi.ID %>%
# unique() -> unipro

biomart_export$Gene.name %>%
  unique() -> symbol
  
library(readr)
library(tidyr)

#Get 3 files from HPA
celljunc<- read_tsv("./HPA/subcell_location_Cell.tsv.gz")
plasma <- read_tsv("./HPA/subcell_location_Plasma.tsv.gz")
vesicles <- read_tsv("./HPA/subcell_location_Vesicles.tsv.gz") 


membrane <- rbind(celljunc,plasma)%>%
  rbind(vesicles)%>%
  unique()
  
write_csv(membrane,"./membrn.gene.txt")