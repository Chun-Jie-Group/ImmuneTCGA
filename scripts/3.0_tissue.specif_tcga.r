
library(dplyr)
#path--------------------
expr_path <- "/data/TCGA/TCGA_data"
out_path <- "/data/fux/github/ImmuneTCGA/expression_analys/tissue.specif.tcga"
diff_path <- "/data/fux/github/ImmuneTCGA/diff_expr"
expr.analy_path <- "/data/fux/github/ImmuneTCGA/expression_analys"

#infile--------------------------
expr <- readr::read_rds(file.path(expr_path, "pancan33_expr.rds.gz"))
gene_list <- readr::read_tsv(file.path(diff_path,"tsv_gene_fc_pvalue_filter.tsv"))


casp <- readr::read_tsv()

# filter out genes ---------------
filter_gene_list <- function(.x, gene_list) {
  gene_list %>%
    dplyr::select(symbol) %>%
    dplyr::left_join(.x, by = "symbol")
}

expr %>%
  dplyr::mutate(filter_expr = purrr::map(expr, filter_gene_list, gene_list = gene_list)) %>%
  dplyr::select(-expr) -> gene_list_expr


#----------------------------

  
get_samples <-function(.x,.y){
    .y %>%
      tibble::add_column(cancer_types = .x, .before = 1) -> df
 
  samples <-  tibble::tibble(barcode = colnames(df)[-c(1:3)]) %>%
    dplyr::mutate(
      sample = stringr::str_sub(
        string = barcode,
        start = 1,
        end = 12
      ),#sample is the patient id
      type = stringr::str_split(barcode, pattern = "-", simplify = T)[, 4] %>% stringr::str_sub(1, 2)
    ) %>%
    #type is 01/11
    dplyr::filter(type %in% c("01", "11")) %>%
    dplyr::mutate(type = plyr::revalue(
      x = type,
      replace = c("01" = "Tumor", "11" = "Normal"),
      warn_missing = F
    )) %>%
    dplyr::group_by(sample) %>%
    dplyr::filter(n() >= 2, length(unique(type)) == 2) %>%#for every patient,he'll be included in if there are 01&11 samples sample
    dplyr::ungroup()
  
  sample_type_summary <- table(samples$type) %>% as.numeric()
  if (gtools::invalid(sample_type_summary) ||
      any(sample_type_summary < c(10, 10))) {
    return(NULL)
  }
  
  df_f <-
    df %>%
    dplyr::select(c(1, 2, 3), samples$barcode) %>%
    tidyr::gather(key = barcode, value = expr, -c(1, 2, 3)) %>%
    dplyr::left_join(samples, by = "barcode")
  
  df_f %>%
    dplyr::group_by(cancer_types, symbol, entrez_id, type) %>%
    tidyr::drop_na(expr) %>%
    dplyr::summarise(mean = mean(expr)) %>%
    dplyr::ungroup()  -> df_mean
  
  }
  


purrr::map2(.x = gene_list_expr$cancer_types,
            .y = gene_list_expr$filter_expr,
            .f = get_samples)-> sample.pair


names(sample.pair) <- gene_list_expr$cancer_types

sample.pair %>% dplyr::bind_rows() -> sample_mean


sample_mean%>%
  filter(type=="Normal")%>%
  reshape2::melt(id=c("cancer_types","symbol"),measure.vars="mean")%>%
  reshape2::dcast(symbol~cancer_types)->sample.df

#rm(sample_melt)

#write-----------------------------------------------------
readr::write_tsv(sample.df,file.path(out_path,"tisp_tcga.tsv"))


