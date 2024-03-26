make_metadata <- function() {
  read_tsv("info/metadata.txt", show_col_types = FALSE) |> 
    mutate(
      condition = as_factor(condition) |> fct_relevel("PABPC1")
    ) |> 
    arrange(condition, replicate) |> 
    mutate(sample = as_factor(sample))
}


GSE_EXAMPLE <- list(
  ontology = "go",
  term_id = "GO:0016055",
  contrast = "conditionYTHDF2"
)
