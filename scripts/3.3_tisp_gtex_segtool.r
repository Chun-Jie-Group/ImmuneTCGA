rm(list=ls())
#Lib--------------

library(dplyr)
library(stringr)

#Path----------------------

seg_path <- "/home/fux/fux/github/ImmuneTCGA/expression_analys/tisp.GTEx/segtool_result/seg_high"
expr_path <- "/home/fux/fux/github/ImmuneTCGA/expression_analys"
out_path <- "/home/fux/fux/github/ImmuneTCGA/expression_analys/tisp.GTEx/select_gene"
#-----------------------

result <- list.files(seg_path)
tisp.hpa <- readr::read_tsv(file.path(expr_path,"01_casp.tsv")) 
cancer_tissue <- xlsx::read.xlsx(file.path(expr_path,"cancer_tissue.xlsx"),1)


#Replace-------------------------------

test <- stringr::str_split(result[1],"\\.")
test[[1]][1]
levels(cancer_tissue$GTEx) <- c(levels(cancer_tissue$GTEx),"Gut")
cancer_tissue[3,3] <- "Gut"

#Function-----------------------------

segtool_high <- function(fi){
  tissue_type_00 <- stringr::str_split(fi,"\\.")
  tissue_type <- tissue_type_00[[1]][1]
  df <- readr::read_tsv(file.path(seg_path,fi))
  colnames(df) <- df[1,]
  df <- df[c(-1,-nrow(df)),]
  df%>%
    tibble::add_column(tissue_types =tissue_type, .before = 1) -> df
}
#Program------------------------------------

seg_high <- segtool_high(result[1])

for(i in 2:length(result)){
  df <- segtool_high(result[i])
  seg_high%>%
    rbind(df)->seg_high
}

seg_high_01 <- unique(seg_high)

#Overlap---------------------------------

colnames(seg_high_01) <- c("tissue_types","symbol")
tissue_type <- unique(seg_high$tissue_types)


overlap <- function(i)
{
  cancer_tissue%>%
    filter(cancer_tissue$GTEx==tissue_type[i])->tisp_00
  tisp.hpa%>%
    filter(cancer_types %in% tisp_00$cancer_type) %>%
    select(symbol)->tisp.hpa_01
  seg_high_01%>%
    filter(seg_high_01$tissue_types==tissue_type[i])%>%
    select(symbol)->tisp.gtex_01
  
  merge_00 <- intersect(tisp.hpa_01,tisp.gtex_01)
  merge_00%>%
    tibble::add_column(tissue_type =tissue_type[i], .before = 1) -> merge_01
}

tisp.gtex <- overlap(1)

for (i in 2:length(tissue_type)) {
  tisp.gtex_01 <- overlap(i)
  tisp.gtex%>%
    rbind(tisp.gtex_01)->tisp.gtex
}


