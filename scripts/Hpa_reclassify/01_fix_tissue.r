data_path.hpa <- "/home/fux/fux/data/ImmuneTCGA/HPA"
diffgene_path <- "/home/fux/fux/github/ImmuneTCGA/diff_expr"
work_path <- "/home/fux/fux/github/ImmuneTCGA/expression_analys/HPA_reclassify_090518"

diffgene_01 <-  readr::read_tsv(file.path(diffgene_path,"diffgene_HKexclu.tsv"))
hpa_tissue_classify <- readxl::read_xlsx("/home/fux/fux/github/ImmuneTCGA/expression_analys/HPA_reclassify_090518/HPA_tissue_classify.xlsx",
                                         1,
                                         col_names = T)

expr_hpa <- readr::read_tsv(file.path(data_path.hpa,"normal_tissue.tsv.zip"))%>%
  dplyr::filter(`Gene name` %in% diffgene_01$symbol)

get_tmark <- function(ti)
  {
  hpa_tissue_classify%>%
    dplyr::filter(Organ == ti[3])->htc
  
  if(nrow(htc) != 0)
  {
    return(htc$TissueMark) }
  else
  {
    return(ti[3])}
  
}


expr_hpa%>%
  dplyr::mutate(tmark=apply(expr_hpa,1,get_tmark))->expr_hpa_fixed
expr_hpa_fixed$tmark %>%
  str_replace_all("mucle","muscle")->expr_hpa_fixed$tmark

readr::write_tsv(expr_hpa_fixed,file.path(work_path,"hpa_expr_tissuefixed.tsv"))
