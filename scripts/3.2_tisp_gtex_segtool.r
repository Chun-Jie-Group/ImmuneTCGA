rm(list=ls())

#Lib--------------------------------------------------------------------
library(dplyr)
#Path-------------------------------------------------------------------

GTEx_path <- "/data/GTEx/V7"
diff_path <- "/data/fux/github/ImmuneTCGA/diff_expr"
out_path <- "/home/fux/fux/github/ImmuneTCGA/expression_analys/tisp.GTEx"
exp_path <- "/home/fux/fux/github/ImmuneTCGA/expression_analys"

#Input file--------------------------------------------------------------

GTEx_expr <- readr::read_rds(file.path(GTEx_path,"gtex_gene_tmp_annotation_phenotype_v7.rds.gz"))
GTEx_expr <- GTEx_expr%>%
  select("SMTS","expr")
gene_list <- readr::read_tsv(file.path(diff_path,"tsv_gene_fc_pvalue_filter.tsv"))
cancer_tissue <- readxl::read_xlsx(file.path(exp_path,"cancer_tissue.xlsx"))

#Tissue name consensus-------------------------------------------------------------------------
  
GTEx_expr$SMTS %>%
  stringr::str_replace_all("Colon","Gut") %>%
  stringr::str_replace_all("Small Intestine","Gut") -> GTEx_expr$SMTS

cancer_tissue%>%
  mutate(fixed_GTEx=cancer_tissue$GTEx) -> cancer_tissue
  cancer_tissue[3,4] <- "Gut"

GTEx_expr%>%
  filter(GTEx_expr$SMTS %in% cancer_tissue$fixed_GTEx)->GTEx_expr

#Exclude------------------------------------------------------------------------

doubt_gene <-c("CSF2RA","IL3RA","IL9R","JAG1","LSP1","P2RY8") 
dupli_gene <- function(tb){
  tb %>%
    filter(symbol %in% doubt_gene)->tb_01
  tb %>%
    anti_join(tb_01)->tb_02
  return(tb_02)
}

#First part--------------------------------------------------------------------------------------
i <- 1
GTEx_expr[[i,2]]%>%
  filter(.$symbol %in% gene_list$symbol) -> gene_list_expr

gene_list_expr %>%
  select(-1,-2) %>%
  mutate(average=apply(., 1,mean))->gene_list_expr_00


gene_list_expr %>%
  select(symbol) %>%
  mutate(average=gene_list_expr_00$average)->gene_list_expr_01

gene_list_expr_02 <- dupli_gene(gene_list_expr_01)

colnames(gene_list_expr_02)[2] <- GTEx_expr[[i,1]]

GTEx_matrix <- gene_list_expr_02

#Loop--------------------------------------------------------------------------

for (i in 2:nrow(GTEx_expr)) {
  GTEx_expr[[i,2]]%>%
    filter(.$symbol %in% gene_list$symbol) -> gene_list_expr
  
  gene_list_expr %>%
    select(-1,-2) %>%
    mutate(average=apply(., 1,mean))->gene_list_expr_00
  
  gene_list_expr %>%
    select(symbol) %>%
    mutate(average=gene_list_expr_00$average)->gene_list_expr_01
  
  gene_list_expr_02 <- dupli_gene(gene_list_expr_01)
  
  colnames(gene_list_expr_02)[2] <- GTEx_expr[[i,1]]
  
  GTEx_matrix%>%
    left_join(gene_list_expr_02,by="symbol")->GTEx_matrix
}



#Test------------------------------------------------------------------

#GTEx_matrix %>%
 # filter(symbol=="IL3RA")
  
#Write---------------------------------------------------------------
row.names(GTEx_matrix) <- GTEx_matrix$symbol
GTEx_matrix <- GTEx_matrix %>%
  select(-1)
write.table(GTEx_matrix,
                 file.path(out_path,"GTEx_matrix_02.tsv"))
