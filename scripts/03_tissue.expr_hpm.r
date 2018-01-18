#the objects are consistent with 01

out_path <- "/home/fux/fux/github/ImmuneTCGA/expression_analys"

diffgene$target <- 0


for (i in 1:nrow(diffgene)) {
  expr_hpm %>%
    filter(Gene==diffgene[i,]$symbol)%>%
    as.matrix()->expr_hpm_00
  
  l <- length(expr_hpm_00)
  
  if(l!=0){
    
    expr_hpm_00 <- expr_hpm_00%>%
      .[1,-1]%>%
      as.numeric()
    
    expr_hpm.sort <- sort(expr_hpm_00,decreasing = T)
    zero.p <- sum(expr_hpm_00==0)/l
    
    if(zero.p<=0.5){
      diffgene[i,11] <- -100
      next
    }
    
    if(expr_hpm.sort[1]==0){
      diffgene[i,11] <- -1
      next
    }
    if(expr_hpm.sort[2]==0){
      diffgene[i,11] <- 1
    }else if(expr_hpm.sort[1]/expr_hpm.sort[2]>=10){
      diffgene[i,11] <- 1
    }
    
    
  }
}

#-----------------
readr::write_tsv(target.hpam_01,file.path(out_path,"target.hpm_01.tsv"))
