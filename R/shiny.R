write_rds_name <- function(obj) {
  path <- file.path("shiny", "data")
  if (!dir.exists(path)) dir.create(path)
  obj_name <- deparse(substitute(obj))
  file_name <- file.path(path, str_glue("{obj_name}.rds"))
  write_rds(obj, file_name, compress = "xz")
}

save_data_for_shiny <- function(bm_genes, star, edger, terms, fterms, gse, clock_genes) {
  de <- edger |> 
    rename(id = gene_id, p_value = PValue, fdr = FDR, log_fc = logFC, log_exp = logCPM)
  
  data <- star$dat |> rename(id = gene_id)
  metadata <- star$metadata |> rename(group = condition)
  features <- star$genes |> rename(id = gene_id, name = gene_symbol)
  
  write_rds_name(data)
  write_rds_name(metadata)
  write_rds_name(features)
  write_rds_name(de)
  write_rds_name(fterms)
  write_rds_name(terms)
  write_rds_name(gse)
  write_rds_name(clock_genes)
}



make_de_shiny <- function(edger_pabpc1_e1, edger_total_e1, edger_pair_e1, de_norm) {
  bind_rows(edger_pabpc1_e1 |> mutate(base = "Base PABPC1"),
            edger_total_e1 |> mutate(base = "Base Total"),
            edger_pair_e1 |> add_column(base = "Pairwise"),
            de_norm |> add_column(base = "Normalised to Total")
  ) |>
    select(gene_id, gene_symbol, description, logFC, logCPM, PValue, FDR, contrast, base)
}


make_star_shiny <- function(star_e1, star_norm) {
  star_e1$dat <- star_e1$dat |> 
    left_join(star_norm$dat |> select(gene_id, sample, logratio_Total))
  star_e1
}
