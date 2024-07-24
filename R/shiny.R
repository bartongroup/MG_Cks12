write_shiny_file <- function(obj, method = c("qs", "rds")) {
  method <- match.arg(method)
  path <- file.path("shiny", "data")
  if (!dir.exists(path)) dir.create(path)
  obj_name <- deparse(substitute(obj))
  cat(paste("  Writing", obj_name))
  if(method == "qs") {
    file_name <- file.path(path, str_glue("{obj_name}.qs"))
    qs::qsave(obj, file_name)
  } else if(method == "rds") {
    file_name <- file.path(path, str_glue("{obj_name}.rds"))
    write_rds(obj, file_name, compress = "xz")
  }
  cat("\n")
}

 
save_data_for_shiny <- function(star, de1, de2, terms, fterms) {
  de1 <- de1 |>
    select(id, log_fc = logFC, expr = logCPM, p_value = PValue, contrast)
  de2 <- de2 |>
    select(id, log_fc = logFC, expr = logCPM, p_value = PValue, contrast)
  data <- star$dat |>
    filter(!bad) |>
    mutate(value = count_norm) |>
    select(id, sample, value)
  metadata <- star$metadata |>
    filter(!bad) |>
    select(-c(bad, cell_line, time_point, raw_sample, experiment))
  features <- star$genes |>
    select(id, name = gene_symbol, description)
  
  write_shiny_file(de1)
  write_shiny_file(de2)
  write_shiny_file(data)
  write_shiny_file(metadata)
  write_shiny_file(features)
  write_shiny_file(fterms)
  write_shiny_file(terms)
}


