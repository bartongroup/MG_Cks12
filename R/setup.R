make_metadata <- function(file) {
  read_excel(file) |> 
    clean_names() |>
    mutate(raw_sample = if_else(experiment == "EK_042",
      str_extract(novogene_sample_name, "EK_\\d+_\\d+"),
      str_extract(novogene_sample_name, "^A\\d\\d")
    )) |> 
    mutate(raw_sample = if_else(experiment == "Pilot",
      str_c("A_", str_extract(raw_sample, "\\d+") |> as.numeric() |> as.character()),
      raw_sample
    )) |> 
    select(experiment, cell_line, treatment, time_point = timepoint, raw_sample) |> 
    mutate(treatment = treatment |> str_remove("_[AB]") |> str_replace("siCks1_2", "siCks1.2")) |> 
    distinct() |> 
    group_by(treatment, cell_line, time_point) |> 
    mutate(replicate = row_number() |> as.character()) |> 
    ungroup() |> 
    unite(sample, c(treatment, replicate), remove = FALSE) |>
    mutate(
      Cks1 = if_else(str_detect(treatment, "1"), "+", "-") |> fct_relevel("-"),
      Cks2 = if_else(str_detect(treatment, "2"), "+", "-") |> fct_relevel("-"),
      experiment = fct_relevel(experiment, "Pilot"),
      time_point = fct_relevel(time_point, "24h"),
      treatment = fct_relevel(treatment, "siNT"),
      bad = cell_line == "HeLa" | time_point == "48h"
    ) |> 
    arrange(bad, treatment, time_point, replicate, experiment) |> 
    mutate(sample = as_factor(sample))
}

