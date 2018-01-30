rm(list=ls())
#Lib--------------

library(dplyr)

#Path----------------------

seg_path <- "/home/fux/fux/github/ImmuneTCGA/expression_analys/tissue.specif.tcga/segtool_result_high"
expr_path <- "/home/fux/fux/github/ImmuneTCGA/expression_analys"
#-----------------------

result <- list.files(seg_path)
tisp.hpa <- readr::read_tsv(file.path(expr_path,"01_casp.tsv")) 
#Function-----------------------------

segtool_high <- function(fi){
  cancer_type <- stringr::str_sub(fi,1,4)
  df <- readr::read_tsv(file.path(seg_path,fi))
  colnames(df) <- df[1,]
  df <- df[c(-1,-nrow(df)),]
  df%>%
    tibble::add_column(cancer_types =cancer_type, .before = 1) -> df
}

#Program------------------------------------

seg_high <- segtool_high(result[1])

for(i in 2:length(result)){
  df <- segtool_high(result[i])
  seg_high%>%
    rbind(df)->seg_high
}

#Overlap--------------------------------------

cancer_type <- unique(tisp.hpa$cancer_types)
colnames(seg_high) <- c("cancer_types","symbol")

seg_high%>%
  filter(cancer_types==cancer_type[1])->seg_high_01
tisp.hpa%>%
  filter(cancer_types==cancer_type[1])%>%
  select("cancer_types","symbol")->tisp.hpa_01
merge_00 <- intersect(seg_high_01,tisp.hpa_01)

for (i in 2:length(cancer_type)) {
  seg_high%>%
    filter(cancer_types==cancer_type[i])->seg_high_01
  tisp.hpa%>%
    filter(cancer_types==cancer_type[i])%>%
    select("cancer_types","symbol")->tisp.hpa_01
  merge_00 <- merge_00%>%
    rbind(intersect(seg_high_01,tisp.hpa_01))
    
}

#-------------------------------------------------

tisp.hpa%>%
  filter(symbol %in% merge_00$symbol)->merge_00

readr::write_tsv(merge_00,file.path(seg_path,"01_merge_hpa.taga.tsv"))

