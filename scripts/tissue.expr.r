
source_path <- "/home/fux/fux/github/ImmuneTCGA/expression_analys"

target_hpa <- readr::read_tsv(file.path(source_path,"02_tissue.expr_hpa.tsv"))%>%
  filter(target==1)

target_hpm <- readr::read_tsv(file.path(source_path,"target.hpm_01.tsv"))

target_temp <- intersect(target_hpa,target_hpm)                              
