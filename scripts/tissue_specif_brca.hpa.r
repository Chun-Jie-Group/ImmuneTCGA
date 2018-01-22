
#-------------------------------------------

tissue_group <- function(expres,symbol){
  
  expres%>%
    filter(`Gene name`==symbol)%>%
    select(`Gene name`,Tissue,Level) %>%
    unique()->expr.hpa_01
  
  if(nrow(expr.hpa_01)==0){
    return(0)
  }else{
    expr.hpa_01%>%
      count(Tissue) %>%
      filter(n>1)->tissue_00
    
    expr.hpa_01%>%
      count(Tissue)%>%
      filter(n==1)->tissue_01
    
    expr.hpa_01%>%
      filter(Tissue %in% tissue_01$Tissue)->tissue_01
    
    for (t in tissue_00$Tissue) {
      tissue_03 <- expr.hpa_01%>%
        filter(Tissue==t)
      if("High" %in% tissue_03$Level){
        tissue_01 <- rbind(tissue_01,tissue_03[which(tissue_03$Level=="High"),])
        next
      }
      if("Medium" %in% tissue_03$Level){
        tissue_01 <- rbind(tissue_01,tissue_03[which(tissue_03$Level=="Medium"),])
        next
      }
      if("Low" %in% tissue_03$Level){
        tissue_01 <- rbind(tissue_01,tissue_03[which(tissue_03$Level=="Low"),])
        next
      }
      
    }
  }
  return(tissue_01)
}


#--------------------------------------------

#In BRCA

diffgene%>%
  filter(cancer_types=="BRCA")->diffgene.brca
expr_hpa%>%
  filter(Tissue=="breast",`Gene name` %in% diffgene.brca$symbol,Level=="High"|Level=="Medium")->expr.hpa.brca

brca_info <- data.frame(gene=unique(expr.hpa.brca$`Gene name`),level=character(82),pro=numeric(82),stringsAsFactors = F)

for (i in 1:nrow(brca_info)) {
  brca_selec <- tissue_group(expr.hpa.brca,brca_info[i,1])
  brca_info[i,2] <- brca_selec$Level
  gene_expr <- tissue_group(expr_hpa,brca_info[i,1])
  brca_info[i,3] <- (sum(gene_expr$Level=="High")+sum(gene_expr$Level=="Medium"))/nrow(gene_expr)
  
}

brca_tar <- brca_info%>%
  filter(pro %in% sort(brca_info$pro,decreasing = F)[1:7])

#-------------------------------------------------------
 