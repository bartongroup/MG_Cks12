top_expressed_genes <- function(set, gns, group_var = "treatment", n_top = 10) {
  genes <- gns |> 
    select(chr, start, end, id, gene_symbol, description) |> 
    distinct()
  
  mn <- star_mean_expr(set, what = "count_norm", group_var = group_var)
  mn |> 
    arrange(-M) |> 
    mutate(group = get(group_var)) |> 
    group_split(group) |> 
    map_dfr(\(w) head(w, n_top)) |> 
    left_join(genes, by = "id")
}

remove_mt <- function(set, gns) {
  mt <- gns |> 
    filter(chr == "MT") |> 
    pull(gene_id)
  set$sel <- setdiff(set$sel, mt)
  set
}

coding_only <- function(set, gns) {
  cod <- gns |> 
    filter(gene_biotype == "protein_coding") |> 
    pull(gene_id)
  set$sel <- intersect(set$sel, cod)
  set
}


read_fastq_lengths <- function(config, meta) {
  trimmed_path <- file.path(config$data_path, config$dir$fastq_trimmed)
  f_trim <- read_fastq_length(trimmed_path, meta) |> 
    add_column(type = "Trimmed")
  clean_path <- file.path(config$data_path, config$dir$clean_trimmed)
  f_clean <- read_fastq_length(clean_path, meta) |> 
    add_column(type = "Clean")
  bind_rows(f_trim, f_clean)
}




make_upset_contrasts <- function(de) {
  contrasts <- levels(de$contrast)
  map(contrasts, function(ctr) {
    de |> 
      filter(contrast == ctr & sig) |> 
      pull(gene_id)
  }) |> 
    set_names(contrasts)
}

make_upset <- function(edger_rest, edger_full, edger_cond) {
  list(
    rest = edger_rest |> filter(contrast == "YTHDF-Total" & sig) |> pull(gene_id),
    full = edger_full |> filter(contrast == "conditionYTHDF" & sig) |> pull(gene_id),
    cond = edger_cond |> filter(contrast == "conditionYTHDF" & sig) |> pull(gene_id)
  )
}

# Minimum count per sample in each experiment
filter_per_experiment <- function(set, min_count) {
  md <- set$metadata |> 
    select(sample, experiment)
  set$dat <- set$dat |> 
    left_join(md, by = join_by(sample))  |> 
    group_by(id, experiment) |> 
    mutate(bad = max(count) < min_count)  |> 
    ungroup() |> 
    group_by(id) |> 
    mutate(bad = any(bad)) |> 
    ungroup() |> 
    select(-experiment)
  set$sel <- set$dat |> 
    filter(!bad) |> 
    pull(id) |> 
    unique()
  
    set  
}