rm(list=ls())

library(stringr)
library(tidyr)

data_path.hpa="/home/fux/fux/data/ImmuneTCGA/HPA"
data_path.hpm="/home/fux/fux/data/ImmuneTCGA/HPM"
out_path="/home/fux/fux/github/ImmuneTCGA/expression_analys"
diffgene_path="/home/fux/fux/github/ImmuneTCGA/diff_expr"


expr_hpa <- readr::read_tsv(file.path(data_path.hpa,"normal_tissue.tsv.zip"))%>%
  filter(`Gene name` %in% diffgene$symbol)

diffgene <- readr::read_tsv(file.path(diffgene_path,"tsv_gene_fc_pvalue_filter.tsv"))
expr_hpm <- readr::read_csv(file.path(data_path.hpm,"HPM_gene_level_epxression_matrix_Kim_et_al_052914.csv"))
tissue_nomen <- readxl::read_excel(file.path(out_path,"tissue_nomenclature.xlsx"))
HK_genes <- read.table("/home/fux/fux/github/ImmuneTCGA/expression_analys/HK_genes.txt")

#------------------------------------------

diffgene%>%
  filter(symbol %in% HK_genes$V1)->HK_diffgene
diffgene%>%
  anti_join(HK_diffgene)%>%
  mutate(target=0)->diffgene

#------------------------------

colnames(tissue_nomen) <- c("X","tissue","HPA","HPM","cancer_type")
tissue_nomen$HPA%>%
  stringr::str_replace_all("'","")%>%
  stringr::str_replace_all("\\\r\\\n"," ")->tissue_nomen$HPA
tissue_nomen$HPM%>%
  stringr::str_replace_all("'","")%>%
  stringr::str_replace_all("\\\r\\\n"," ")%>%
  stringr::str_replace_all("\\."," ")->tissue_nomen$HPM

#------------------------------------------
#tissue nomenlature correction

expr_hpa$Tissue %>%
  str_replace_all("adrenal gland","adrenal")%>%
  str_replace_all("urinary bladder","bladder") %>%
  str_replace_all("bone marrow","bone") %>%
  str_replace_all("cerebellum","brain") %>%
  str_replace_all("hippocampus","brain")%>%
  str_replace_all("lateral ventricle","brain")%>%
  str_replace_all("cerebral cortex","brain")%>%
  str_replace_all("hypothalamus","brain")%>%
  str_replace_all("pituitary gland","brain")%>%
  str_replace_all("retina","brain")%>%
  
  str_replace_all("cervix, uterine","cervix") %>%
  str_replace_all("colon","gut") %>%
  str_replace_all("duodenum","gut")%>%
  str_replace_all("small intestine","gut")%>%
  
  str_replace_all("heart muscle","heart") %>%
  str_replace_all("oral mucosa","oropharynx") %>%
  str_replace_all("salivary gland","oropharynx")%>%
  
  str_replace_all("parathyroid gland","parathyroid") %>%
  str_replace_all("seminal vesicle","seminal") %>%
  str_replace_all("skin 1","skin") %>%
  str_replace_all("skin 2","skin")%>%
  
  str_replace_all("soft tissue 1","soft tissue") %>%
  str_replace_all("soft tissue 2","soft tissue")%>%
  
  str_replace_all("stomach 1","stomach") %>%
  str_replace_all("stomach 2","stomach")%>%
  
  str_replace_all("thyroid gland","thyroid") %>%
  str_replace_all("endometrium","uterus") %>%
  str_replace_all("endometrium 1","uterus")%>%
  str_replace_all("endometrium 2","uterus")%>%
  str_replace_all("uterus 1","uterus")%>%
  str_replace_all("uterus 2","uterus")%>%
  
  str_replace_all("lactating breast","breast")->expr_hpa$Tissue
  
#---------------------------------------------------------

tissue_nomen%>%
  filter(!is.na(cancer_type))%>%
  select(tissue,cancer_type)->cancer_tissue

#------------------------

fetal_tissue <- c("liver","heart","brain","ovary","gut","testis")
diffgene%>%
  mutate(target=0)

#---------------------------------------------------------
for(i in 1:nrow(diffgene)){
  
  diffgene[i,]->gene
  tissue.rela <- cancer_tissue$tissue[which(grepl(gene$cancer_types,cancer_tissue$cancer_type))]
  
  expr.hpa <- expr_hpa%>%
    filter(`Gene name`==gene$symbol)
  expr.hpa.fetal <- expr.hpa%>%
    filter(Tissue %in% fetal_tissue)
  
  expr.level <- expr.hpa[which(expr.hpa$Tissue==tissue.rela),]$Level
  
  
  gene_hp <- sum(expr.hpa$Level=="High")/nrow(expr.hpa)
  gene_mp <- sum(expr.hpa$Level=="Medium")/nrow(expr.hpa)
  gene_np <- sum(expr.hpa$Level=="Not detected")/nrow(expr.hpa)
  
  fetal_hp <- sum(expr.hpa.fetal$Level=="High")/nrow(expr.hpa.fetal)
  
 if(length(expr.level)!=0){
   
   if((gene_mp+gene_hp)>0.5){
     diffgene[i,11] <- -1
     next
   }#we don't want gene expressed so much broad
   
   if("High" %in% expr.level){
     if(gene_hp<0.05)diffgene[i,11] <- 1
   }else{
     if(fetal_hp>0.1)diffgene[i,11] <- -1
   } 
   
 }else{diffgene[i,11] <- -100}
  
}

#--hpm-----------------------------------------------------

diffgene%>%
  filter(target==1)%>%
  count(cancer_types)
