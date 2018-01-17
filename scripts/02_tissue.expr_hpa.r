#the objections are consistent with 01
out_path <- "/home/fux/fux/github/ImmuneTCGA/expression_analys"

#hk_gene had been filtered out
diffgene$target <- 0

tissue_group <- function(symbol){
 
   expr_hpa%>%
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

#-------------------------------------------------------------------------

for (i in 1:nrow(diffgene)) {
  
    expr.hpa <- tissue_group(symbol = diffgene[i,]$symbol)
  if(expr.hpa!=0){
    sum_high <- 0
    
    sum_high <- sum(expr.hpa$Level=="High")
    mp <- sum(expr.hpa$Level=="Medium")/nrow(expr.hpa)
    lp <- sum(expr.hpa$Level=="Low")/nrow(expr.hpa)
    
  if(sum_high>0){
    
    if(sum_high<=5&&(mp+lp)<=0.1)diffgene[i,11] <- 1
  }else{
    diffgene[i,11] <- -1
  }
  
  
  }else{
  diffgene[i,11] <- -100
}
}
#-------------------------------


