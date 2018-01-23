
out_path <- "/data/fux/github/ImmuneTCGA/expression_analys"
cancer_tissue <- readxl::read_xlsx("/data/fux/github/ImmuneTCGA/expression_analys/cancer_tissue.xlsx")

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


#---------------------------------------------------------------------

tissue_specif <- function(can){
  
  cancer_tissue%>%
    filter(cancer_type==can)->cancer_tissue_00
  
  diffgene.cancer.spe <- diffgene%>%
    filter(cancer_types==can)
  
  expr_hpa%>%
    filter(Tissue %in% cancer_tissue_00$tissue,
           `Gene name` %in% diffgene.cancer.spe$symbol,
           Level=="High"|Level=="Medium")->expr.hpa
  
 if(nrow(expr.hpa)==0){
   return(0)
   }
  
  n.gene <- length(unique(expr.hpa$`Gene name`))
  
  can_info <- data.frame(gene=unique(expr.hpa$`Gene name`),
                         level.HPA=character(n.gene),Ratio.HPA=numeric(n.gene),
                         stringsAsFactors = F)
  
  for (i in 1:nrow(can_info)) {
    casp_selec <- tissue_group(expr.hpa,can_info[i,1])
    can_info[i,2] <- casp_selec$Level
    gene_expr <- tissue_group(expr_hpa,can_info[i,1])
    can_info[i,3] <- (sum(gene_expr$Level=="High")+sum(gene_expr$Level=="Medium"))/nrow(gene_expr)
    
  }
  
  diffgene%>%
    select(cancer_types,symbol,Tumor,Normal,fc,fdr)%>%
    filter(cancer_types==can,symbol %in% can_info$gene)%>%
    left_join(can_info,by=c("symbol"="gene"))->diffgene.cancer.spe_00
 
   diffgene.cancer.spe_00[order(diffgene.cancer.spe_00$Ratio.HPA,decreasing=F),]->diffgene.cancer.spe_01
   diffgene.cancer.spe_01%>%
     filter(Ratio.HPA<0.3)->diffgene.cancer.spe_02
  return(diffgene.cancer.spe_02)
}

#----------------------------------------------------------------------------------------
casp.blca <- tissue_specif("BLCA")
casp.hnsc<- tissue_specif("HNSC")
casp.brca <- tissue_specif("BRCA")
casp.coad<- tissue_specif("COAD")
casp.kich <- tissue_specif("KICH")
casp.kirp<- tissue_specif("KIRP")

casp.kicr<- tissue_specif("KICR")
casp.esca <- tissue_specif("ESCA")

casp.lihc <- tissue_specif("LIHC")
casp.lusc <- tissue_specif("LUSC")
casp.luad <- tissue_specif("LUAD")
casp.prad <- tissue_specif("PRAD")
casp.stad <- tissue_specif("STAD")
casp.thca <- tissue_specif("THCA")
#----------------------------------------------------------------------------
casp <- data.frame()
casp%>%
  rbind(casp.blca)%>%
  rbind(casp.brca)%>%
  rbind(casp.coad)%>%
  rbind(casp.hnsc)%>%
  rbind(casp.kich)%>%
  rbind(casp.kirp)%>%
  rbind(casp.lihc)%>%
  rbind(casp.luad)%>%
  rbind(casp.lusc)%>%
  rbind(casp.prad)%>%
  rbind(casp.stad)%>%
  rbind(casp.thca)->casp   
#-----------------------------

casp%>%
  filter(level.HPA=="High")->casp_01
casp_01%>%
  count(cancer_types)
readr::write_tsv(casp,
                 file.path(out_path,"3.0_tissue.specifi_hpa.tsv"))
