library(magrittr)
library(ggplot2)

# processed path
tcga_path = "/data/TCGA/TCGA_data"

#output path
expr_path <- "/home/fux/fux/github/ImmuneTCGA/diff_expr"
#expr_path_a <- file.path(expr_path, "")

# Read gene list
# Gene list was compress as rds
expr <- readr::read_rds(file.path(tcga_path, "pancan33_expr.rds.gz"))

gene_list <- readr::read_csv(file.path(expr_path, "ImmuneTCGA_genelist.csv")) 
colnames(gene_list) <- "symbol"

#######################
# filter out genes ---------------
#######################
filter_gene_list <- function(.x, gene_list) {
  gene_list %>%
    dplyr::select(symbol) %>%
    dplyr::left_join(.x, by = "symbol")
}

expr %>%
  dplyr::mutate(filter_expr = purrr::map(expr, filter_gene_list, gene_list = gene_list)) %>%
  dplyr::select(-expr) -> gene_list_expr

readr::write_rds(x = gene_list_expr, path = file.path(expr_path, "ImmuneTCGA_genelist_expr.rds.gz"), compress = "gz")

#################################
# Caculate p-value and fold-change.
##################################
calculate_fc_pvalue <- function(.x, .y) { #what is .x,.y?
  .y %>%
    tibble::add_column(cancer_types = .x, .before = 1) -> df #what is .before for?
  
  # get cancer types and get # of smaple >= 10
  samples <-
    tibble::tibble(barcode = colnames(df)[-c(1:3)]) %>%
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
  }#It means the number of patient should be more than 10?
  #In short,here we get the patient in 14 types of cancer with tumor/normal samples
  #---------------------------------------------------------------------------------------------------
  
  # filter out cancer normal pairs
  df_f <-
    df %>%
    dplyr::select(c(1, 2, 3), samples$barcode) %>%
    tidyr::gather(key = barcode, value = expr, -c(1, 2, 3)) %>%
    dplyr::left_join(samples, by = "barcode")
  
  # pvalue & fdr
  df_f %>%
    dplyr::group_by(cancer_types, symbol, entrez_id) %>%
    tidyr::drop_na(expr) %>%
    dplyr::do(broom::tidy(t.test(expr ~ type, data = .))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
    dplyr::select(cancer_types, symbol, entrez_id, p.value, fdr) -> df_pvalue
  
  # log2 fold change mean
  df_f %>%
    dplyr::group_by(cancer_types, symbol, entrez_id, type) %>%
    tidyr::drop_na(expr) %>%
    dplyr::summarise(mean = mean(expr)) %>%
    tidyr::spread(key = type, mean) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(fc = (Tumor + 0.1) / (Normal + 0.1))-> df_fc
  
  df_fc %>%
    dplyr::inner_join(df_pvalue, by = c("cancer_types", "symbol", "entrez_id")) %>%
    dplyr::mutate(n_normal = sample_type_summary[1], n_tumor = sample_type_summary[2]) -> res
  return(res)
}

purrr::map2(.x = gene_list_expr$cancer_types,
            .y = gene_list_expr$filter_expr,
            .f = calculate_fc_pvalue) -> gene_list_fc_pvalue

names(gene_list_fc_pvalue) <- gene_list_expr$cancer_types

gene_list_fc_pvalue %>% dplyr::bind_rows() -> gene_list_fc_pvalue_simplified

#we get p.value,fdr.fc of the 14 cancer with expression info. of target genes. 

#-----------------------------------------------------------------------------------write.part
readr::write_rds(
  x = gene_list_fc_pvalue_simplified,
  path = file.path(expr_path, "03_ImmuneTCGA_fc_pvalue_simplified.rds.gz"),
  compress = "gz"
)
readr::write_tsv(
  x = gene_list_fc_pvalue_simplified,
  path = file.path(expr_path, "tsv_ImmuneTCGA_fc_pvalue_simplified.tsv")
)
#------------------------------------------------------------------------------------------------------
# check autophagy gene expression distribution-----------------------------
expr_path_a <- "/home/fux/fux/github/ImmuneTCGA/diff_expr/expression_distribution"

gene_list_fc_pvalue_simplified %>% 
  dplyr::left_join(gene_list, by = "symbol") %>% 
  tidyr::drop_na(Tumor, Normal) %>% 
  # dplyr::filter(status %in% c( "p")) %>%
  tidyr::gather(key = nt, value = expr, c(Tumor, Normal)) %>% 
  ggplot(aes(x = log2(expr + 1), color = nt)) +
  geom_density(aes(y = ..density..)) +
  ggthemes::scale_color_gdocs(name = "Type") +
  ggthemes::theme_gdocs() +
  facet_wrap( ~ cancer_types) +
  labs(
    x = "Expression (log2)",
    y = "Frenquncy"
  ) -> p
ggsave(filename = "03_a_expression_distribution_promotion.pdf", plot = p, device = "pdf", path = expr_path_a, width = 6, height = 6)
#--------------------------------------------------------------------------------------------

#------------- filter -----------------------------------------------------------------------
gene_list_fc_pvalue_simplified %>%
dplyr::left_join(gene_list, by = "symbol") -> gene_fc_pvalue #what difference?


filter_fc_pval <- function(.x){
  .x %>% 
    # dplyr::filter(Normal > 10 & Tumor > 10) %>%
    dplyr::filter(fc>2)%>% #only high expressor included
    dplyr::filter(abs(log2(fc)) >= log2(3 / 2), fdr <= 0.05) %>%
    dplyr::filter(p.value<0.05) %>%
    dplyr::mutate(p.value = -log10(p.value)) %>% 
    dplyr::mutate(p.value = ifelse(p.value > 15, 15, p.value)) %>% 
    dplyr::mutate(fdr = -log10(fdr)) %>% 
    dplyr::mutate(fdr = ifelse(fdr > 15, 15, fdr)) %>% 
    dplyr::mutate(fc = ifelse(fc < 1/8, 1/8, ifelse(fc > 8, 8, fc)))
}

gene_fc_pvalue %>% filter_fc_pval() -> gene_fc_pvalue_filter

readr::write_tsv(
  x = gene_fc_pvalue_filter,
  path = file.path(expr_path, "tsv_gene_fc_pvalue_filter.tsv")
)

