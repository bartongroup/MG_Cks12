one_count_file <- function(config, meta, n) {
  path <- file.path(config$data_path, "readcount")
  rs <- meta |> 
    slice(n) |> 
    pull(raw_sample) 
  file.path(path, str_glue("{rs}.txt"))
}


read_and_process_star <- function(config, meta, gene_info, suffix = ".txt", min_count = 10,
                                  fix_names_fun = NULL) {
  readcount_path <- file.path(config$data_path, config$dir$read_count)
  starmap_path <- file.path(config$data_path, config$dir$mapping)
  
  star_column <- star_strand_column(readcount_path, meta, suffix)
  cat(paste0("Detected strand - ", names(star_column), ". Reading STAR column ", star_column, ".\n"))
  
  set <- parse_star_counts(readcount_path, meta, star_column, suffix, fix_names_fun) |> 
    add_star_log(starmap_path, meta) |> 
    add_gene_names(gene_info) |> 
    normalize_star_counts(gene_info) |> 
    filter_star_min_count(min_count)

  return(set)
}

add_star_log <- function(set, path, meta) {
  set$star_log <- parse_star_logs(path, meta)
  return(set)
}

parse_one_star_log <- function(file, smpl) {
  if (!file.exists(file)) {
    warning(paste("File", file, "not found"))
    return(NULL)
  }
  tibble(x = readLines(file)) |> 
    separate(x, c("key", "value"), sep = "\\t", fill = "right") |> 
    drop_na() |> 
    mutate(key = str_remove(key, " \\|")) |> 
    mutate(key = str_remove(key, "^\\s+")) |> 
    mutate(raw_sample = smpl)
}

# parse all star logs
parse_star_logs <- function(path, meta) {
  s2n <- set_names(meta$sample, meta$raw_sample)
  meta$raw_sample |> 
    map(~parse_one_star_log(file.path(path, paste0(.x, "_Log.final.out")), .x)) |> 
    list_rbind() |> 
    mutate(sample = as.character(s2n[raw_sample]))
}

# tabulate star log for printing
tabulate_star_log <- function(slog) {
  slog |> 
    select(-raw_sample) |> 
    pivot_wider(names_from = "sample", values_from = "value") |> 
    rename(Description = key) 
}

parse_one_star_count <- function(file, smpl, column = 2, fix_names_fun = NULL) {
  if (!file.exists(file)) {
    warning(str_glue("File {file} not found."))
    return(NULL)
  }
  d <- read_tsv(file, col_names = FALSE, skip = 4, col_types = "ciii") |> 
    select(1, all_of(column)) |> 
    set_names("gene_id", "count") |> 
    add_column(raw_sample = smpl, .after = "gene_id")
  if (!is.null(fix_names_fun)) {
    d <- d |> mutate(gene_id = fix_names_fun(gene_id))
  }
  d
}

# Parse all star count files. Returns a list of tab (wide matrix) and dat (long tibble).
parse_star_counts <- function(path, meta, column = 2, suffix = ".txt", fix_names_fun) {
  s2n <- set_names(meta$sample, meta$raw_sample)
  dat <- meta$raw_sample |> 
    map(~parse_one_star_count(file.path(path, paste0(.x, suffix)), .x, column, fix_names_fun)) |> 
    list_rbind() |> 
    mutate(sample = as.character(s2n[raw_sample])) |> 
    group_by(gene_id) |>
    mutate(gene_count = sum(count)) |>
    ungroup() |>
    filter(gene_count > 0) |>
    select(gene_id, sample, count) |> 
    mutate(sample = factor(sample, levels = meta$sample))
  tab <- dat |> 
    pivot_wider(id_cols = "gene_id", names_from = "sample", values_from = "count") |> 
    column_to_rownames("gene_id") |> 
    as.matrix()
  list(dat = dat, tab = tab, metadata = meta, sel = rownames(tab))
}


add_gene_names <- function(set, gene_info) {
  set$genes <- tibble(gene_id = rownames(set$tab)) |> 
    left_join(gene_info |> select(gene_id, gene_symbol) |> distinct(), by = "gene_id") |> 
    mutate(gene_symbol = if_else(is.na(gene_symbol), gene_id, gene_symbol))
  set
}

# used in normalise_to_library
normalise_to_size <- function(set, libsize = NULL) {
  if (is.null(libsize)) {
    libsize <- set$dat |> 
      group_by(sample) |> 
      summarise(size = sum(count, na.rm = TRUE)) 
  }
  
  libsize <- libsize |> mutate(normfac = size / mean(size))
  
  # look-up tables are much faster than left_join
  list(
    norm_fac = set_names(libsize$normfac, libsize$sample),
    size_fac = set_names(libsize$size / 1e6, libsize$sample),
    normfac_tab = libsize
  )
  
}

# normalize to library size and RPKM
normalise_to_library <- function(set, gene_info, input_size) {
  
  len_fac = set_names(gene_info$length / 1e3, gene_info$gene_id)
  
  # normalised to total mapped and counted reads
  mapped <- normalise_to_size(set)
  # normalised to total input reads
  input <- normalise_to_size(set, input_size)
  
  
  set$dat <- set$dat |> 
    mutate(
      count_norm = count / mapped$norm_fac[sample] |> unname(),
      rpkm = (count + 1) / (mapped$size_fac[sample] * len_fac[gene_id] |> unname())
      #count_inputnorm = count / input$norm_fac[sample] |> unname,
      #rpkm_inputnorm = (count + 1) / (input$size_fac[sample] * len_fac[gene_id] |> unname),
    )
  
  set$mapped_normfac <- mapped$normfac_tab
  set$input_normfac <- input$normfac_tab
  
  return(set)
}


normalise_edger <- function(set) {
  ed <- edgeR::DGEList(set$tab) |> 
    edgeR::calcNormFactors() |> 
    pluck("samples") |> 
    rownames_to_column("sample") |> 
    rename(normfac = norm.factors) |> 
    select(sample, normfac) |> 
    as_tibble()
  
  # look-up tables are much faster than left_join
  norm_fac <- set_names(ed$normfac, ed$sample)
  
  set$dat <- set$dat |> 
    mutate(count_tmm = count / norm_fac[sample] |> unname())
  set$edger_normfac <- ed
  
  return(set)
}

# regularised log from DESeq2 (originally log2, we convert to log10)
regularised_log <- function(set) {
  rdat <- DESeq2::rlog(set$tab) |> 
    as_tibble(rownames = "gene_id") |>
    pivot_longer(-gene_id, names_to = "sample", values_to = "rlog") |> 
    mutate(
      rlog = rlog / log2(10),
      sample = factor(sample, levels = set$metadata$sample)
    )
  set$dat <- set$dat |> 
    left_join(rdat, by = c("gene_id", "sample"))
 
  return(set) 
}


# edit: added DESeq2's regularised logarithm
normalize_star_counts <- function(set, gene_info) {
  libsize <- set$star_log |> 
    filter(key == "Number of input reads") |> 
    mutate(size = as.numeric(value)) |> 
    select(sample, size)
  
  set |> 
    normalise_to_library(gene_info, libsize) |> 
    #normalise_edger() |> 
    regularised_log()
}

# Normalise to one condition assuming replicates are matched
normalize_to_condition <- function(set, ref_cond, min_count = 10) {
  # Select genes where each replicate in reference is > min_count
  sel_genes <- set$dat |> 
    left_join(set$metadata, by = "sample") |> 
    filter(condition == ref_cond) |> 
    group_by(gene_id) |> 
    summarise(m = min(count)) |> 
    filter(m >= min_count) |> 
    pull(gene_id)
  
  d <- set$dat |> 
    filter(gene_id %in% sel_genes)
  
  nrm <- set$metadata |> 
    filter(condition == ref_cond) |> 
    select(sample, replicate) |> 
    left_join(d, by = "sample") |> 
    select(gene_id, replicate, ref = count_norm)
  
  lr <- d |> 
    left_join(set$metadata |> select(sample, replicate, condition), by = "sample") |> 
    filter(condition != ref_cond) |> 
    droplevels() |> 
    left_join(nrm, by = c("gene_id", "replicate")) |> 
    mutate(logratio = log2((count_norm + 0.5) / ref)) |> 
    select(gene_id, sample, logratio) |> 
    set_names("gene_id", "sample", str_c("logratio_", ref_cond))
  
  set$metadata <- set$metadata |> 
    filter(condition != ref_cond) |> 
    droplevels()
  set$dat <- d |> 
    right_join(lr, by = c("gene_id", "sample")) |> 
    droplevels()
  set$sel <- set$dat$gene_id |> unique()
  
  set
}


plot_star_log <- function(slog, meta, group_var = "condition",
                        descs = c(
                          "Number of input reads",
                          "Uniquely mapped reads %",
                          "% of reads mapped to multiple loci",
                          "% of reads mapped to too many loci",
                          "% of reads unmapped: too short")
                        ) {
  s <- slog |> 
    filter(key %in% descs) |>
    mutate(value = str_remove(value, "%") |> as.numeric()) |> 
    mutate(key = factor(key, levels = descs)) |> 
    select(-raw_sample) |> 
    left_join(meta, by = "sample") |> 
    mutate(group = get(group_var)) |> 
    mutate(value = if_else(key ==  "Number of input reads", value / 1e6, value)) |> 
    mutate(key = recode(key,  "Number of input reads" =  "Number of input reads (millions)",))
  su <- s |> 
    filter(key == "Uniquely mapped reads %") |> 
    arrange(value)
  
  s |> 
    mutate(sample = factor(sample, levels = meta$sample)) |> 
  ggplot(aes(x = sample, y = value, fill = group)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.position = "none"
    ) +
    #geom_segment(aes(xend = sample, yend = 0), colour = "grey70") +
    geom_col(colour = "grey30") +
    facet_wrap(~key, scales = "free_y", ncol = 1) +
    scale_y_continuous(expand = c(0, 0), labels = scientific_10) +
    scale_fill_manual(values = okabe_ito_palette) +
    geom_text(aes(y = value * 1.05, label = "")) + # blank geom to expand axis
    labs(x = NULL, y = "Value")
}

plot_star_log_map <- function(slog, meta) {
  slog |> 
    filter(str_detect(key, "%")) |> 
    mutate(value = str_remove(value, "%") |> as.numeric()) |> 
    mutate(value = na_if(value, 0)) |> 
    mutate(sample = factor(sample, levels = meta$sample)) |> 
    mutate(Description = factor(key)) |> 
  ggplot(aes(x = sample, y = Description, fill = value)) + 
    theme_bw() +
    geom_tile() +
    scale_fill_viridis_c(option = "cividis", trans = "log10", limits = c(0.1, 100), breaks = c(0.1,1,10,100), labels = c(0.1,1,10,100)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) +
    labs(x = NULL, y = NULL, fill = "Percentage")
}


plot_star_sense <- function(file) {
  read_tsv(file, skip = 4, col_names = c("gene", "Unstranded", "First", "Second"), col_types = cols()) |> 
    pivot_longer(-gene, names_to = "column", values_to = "count") |>
    filter(count > 0) |>
  ggplot(aes(x = log10(count), fill = column, group = column)) +
    theme_bw() +
    geom_density(alpha = 0.3)
}

star_col_count <- function(file, smpl) {
  read_tsv(file, skip = 4, col_names = c("gene", "unstranded", "first", "second"), col_types = cols()) |> 
    filter(unstranded > 10) |> 
    summarise(unstranded = sum(unstranded), first = sum(first), second = sum(second)) |> 
    mutate(raw_sample = smpl)
}

star_strand <- function(path, meta, suffix = ".txt", prop.limit = 0.8) {
  s2n <- set_names(meta$sample, meta$raw_sample)
  m <- meta$raw_sample |> 
    map_dfr(~star_col_count(file.path(path, paste0(.x, suffix)), .x)) |> 
    mutate(sample = as.character(s2n[raw_sample])) |>
    mutate(r1 = first / unstranded, r2 = second / unstranded) |> 
    mutate(strand = if_else(r1 > prop.limit, "first", if_else(r2 > prop.limit, "second", "unstranded")))
  ms <- m |> distinct(strand)
  if (nrow(ms) == 1) {
    return(ms$strand)
  } else {
    print("Different stranding detected")
    print(m)
    stop()
  }
}

# wrapper around star_strand, returns STAR column
star_strand_column <- function(...) {
  str2col <- set_names(
    c(3, 4, 2),
    c("first","second", "unstranded")
  )
  
  strand <- star_strand(...)
  str2col[strand]
}

# Select genes with at least min_count count in at least one sample
filter_star_min_count <- function(set, min_count = 10, count_column = "count_norm"){
  set$dat <- set$dat |> 
    group_by(gene_id) |> 
    mutate(good = max(get(count_column)) >= min_count) |> 
    ungroup()
    
  set$sel <- set$dat |> 
    filter(good) |> 
    pull(gene_id) |> 
    unique()

  set  
}


filter_star_samples <- function(set, expr) {
  meta <- set$metadata |> 
    filter(eval(rlang::parse_expr(expr))) |> 
    droplevels()
  smpl_sel <- as.character(meta$sample)
  
  set$dat <- set$dat |> filter(sample %in% smpl_sel)
  set$tab <- set$tab[, smpl_sel]
  set$metadata <- meta  
  
  set
}

# returns genes with zero in at least one group
# gene_id and sum of counts per group in each group
find_zeroes <- function(set, group_var = "group") {
  set$dat |>
    left_join(set$metadata, by = c("sample", "raw_sample")) |>
    mutate(group = get(group_var)) |> 
    group_by(gene_id, group) |>
    summarise(S = sum(count)) |>
    ungroup() |>
    group_by(gene_id) |>
    mutate(m = min(S)) |>
    ungroup() |>
    filter(m == 0) |> 
    pivot_wider(id_cols = gene_id, names_from = group, values_from = S)
}


merge_star_dat <- function(dat, meta, genes, columns = c("group", "replicate")) {
  dat |> 
    left_join(meta, by = "sample") |> 
    left_join(genes, by = "gene_id") |> 
    mutate(gene_symbol = na_if(gene_symbol, "")) |> 
    select(gene_id, gene_symbol, description, sample, all_of(columns), count, count_norm, rpkm)
}


plot_qualities <- function(qc) {
  qc |> 
    ggplot(aes(x = Base, y = Mean, group = root, colour = pair)) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    geom_line(alpha = 0.3, linewidth = 0.6) +
    scale_colour_manual(values = okabe_ito_palette) +
    labs(x = "Base", y = "Mean quality score")
}

plot_cluster_qualities <- function(qc, text.size = 10) {
  tab <- qc |> 
    pivot_wider(id_cols = c(pair, Base), names_from = sample, values_from = Mean) |> 
    select(-c(pair, Base)) |> 
    as.matrix()
  
  hc <- t(tab) |> 
    dist() |> 
    hclust()
  
  dendr <- ggdendro::dendro_data(hc)
  seg <- ggdendro::segment(dendr)
  theme.d <- ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(size = text.size),
    axis.line.y = ggplot2::element_blank(),
    axis.line.x = ggplot2::element_line(size = 0.5),
    axis.ticks.y = ggplot2::element_blank()
  )
  ggplot() +
    theme.d +
    coord_flip() +
    geom_segment(data = seg, aes(x = x, y = y, xend = xend, yend = yend)) +
    scale_x_continuous(breaks = seq_along(dendr$labels$label), labels = dendr$labels$label) +
    scale_y_continuous(expand = c(0,0), limits = c(0, max(seg$y) * 1.03)) +
    scale_colour_manual(values = okabe_ito_palette) +
    labs(x = NULL, y = "Distance")
  
}


plot_map_qual <- function(qc, slog, base = 20) {
  qcf <- qc |> 
    filter(Base == base) |> 
    group_by(raw_sample) |> 
    summarise(qual = mean(Mean))
  slogf <- slog |> 
    filter(key == "Uniquely mapped reads %") |> 
    mutate(mapped = str_remove(value, "%") |> as.numeric())
  qcf |> 
    left_join(slogf, by = "raw_sample") |> 
  ggplot(aes(x = mapped, y = qual)) +
    theme_bw() +
    geom_point() +
    labs(y = paste("Mean read quality at", base), x = "Uniquely mapped reads")
}







despike <- function(set) {
  set$sel <- set$sel |> 
    str_subset("SPIKE", negate = TRUE)
  
  set
}


plot_mapped_count <- function(set) {
  sm <- colSums(set$tab)
  cnt <- tibble(
    sample = names(sm),
    value = sm,
    key = "Reads counted in genes"
  )
  counts <- set$star_log |> 
    filter(key %in% c("Number of input reads", "Uniquely mapped reads number")) |> 
    mutate(key = recode(key, "Number of input reads" = "Input reads", "Uniquely mapped reads number" = "Uniquely mapped reads")) |> 
    select(sample, value, key) |> 
    mutate(value = as.numeric(value)) |> 
    bind_rows(cnt) |> 
    mutate(value = value / 1e6) |> 
    mutate(key = fct_relevel(key, c("Input reads", "Uniquely mapped reads")))
  
  sample_input <- counts |> 
    filter(key == "Input reads") |>
    rename(input = value) |> 
    arrange(input) |> 
    mutate(sample = sample |> as_factor()) |> 
    select(-key)
  
  perc <- counts |> 
    left_join(sample_input, by = "sample") |> 
    mutate(value = 100 * value / input) |> 
    select(-input)
  
  rank_by <- function(w, sel) {
    r <- w |> 
      filter(key == sel) |> 
      mutate(rank = rank(value)) |> 
      select(sample, rank)
    w |> 
      left_join(r, by = "sample")
  }
  
  dat <- bind_rows(
      counts |>
        rank_by("Input reads") |> 
        add_column(what = "Count"),
      perc |>
        rank_by("Uniquely mapped reads") |> 
        add_column(what = "Percentage") |> 
        mutate(rank = rank + 1000)
    )
  
  labs <- dat |> 
    select(sample, rank) |> 
    distinct() |> 
    mutate(rank = as.character(rank))
  
  get_labels <- function(rnk) {
    tibble(rank = rnk) |> 
      left_join(labs, by = "rank") |> 
      pull(sample)
  }
  
  dat |> 
    ggplot(aes(x = as_factor(rank), y = value, colour = key)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    ) +
    geom_segment(aes(xend = as_factor(rank), yend = 0), colour = "grey90") +
    geom_point() +
    scale_colour_manual(values = okabe_ito_palette) +
    scale_x_discrete(labels = get_labels) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    facet_wrap(~ what, scales = "free") +
    labs(x = NULL, y = "Read count (millions)", colour = "Legend")
}


select_strong_genes <- function(set, limit = 100) {
  set$dat |> 
    select(gene_id, sample, count) |> 
    left_join(set$metadata, by = "sample") |> 
    group_by(gene_id, group) |> 
    summarise(min_group_count = min(count)) |> 
    ungroup() |> 
    group_by(gene_id) |> 
    summarise(best_count = max(min_group_count)) |> 
    filter(best_count > limit) |> 
    pull(gene_id) |> 
    unique()
}


save_count_data <- function(set, clock_genes, file, what = "count_norm") {
  set$dat |> 
    mutate(val = get(what)) |> 
    left_join(set$genes, by = "gene_id") |> 
    pivot_wider(id_cols = c(gene_id, gene_symbol), names_from = sample, values_from = val) |> 
    mutate(across(where(is.numeric), ~signif(.x, 4))) |> 
    left_join(clock_genes |> select(-gene_symbol), by = "gene_id") |> 
    relocate(class, .after = "gene_symbol") |> 
    mutate(class = as.character(class)) |> 
    mutate(class = replace_na(class, "-")) |> 
    rename(clock_class = class) |> 
    write_tsv(file)
}

dat2mat <- function(dat, what = "count", id_col = "gene_id", sample_var = "sample") {
  dat |> 
    pivot_wider(id_cols = !!id_col, names_from = !!sample_var, values_from = !!what) |> 
    column_to_rownames(id_col) |> 
    as.matrix()
}


star_mean_expr <- function(set, what = "count", id_var = "gene_id", group_var = "condition") {
  set$dat |> 
    left_join(set$metadata, by = "sample") |> 
    mutate(val = get(what), group = get(group_var), id = get(id_var)) |> 
    group_by(id, group) |> 
    summarise(M = mean(val)) |> 
    ungroup() |> 
    rename(!!id_var := id, !!group_var := group) 
}


read_fastq_length <- function(path, meta) {
  map2_dfr(meta$raw_sample, meta$sample, function(rsam, sam) {
    cfile <- file.path(path, str_glue("{rsam}_count.txt"))
    if(file.exists(cfile)) {
      read_tsv(cfile, col_names = "read_count", show_col_types = FALSE) |> 
        mutate(sample = sam, raw_sample = rsam, .before = 1)
    }
  })
}

mat2dat <- function(tab, to_what = "abu_norm", names = "sample", id_var = "gene_id") {
  tab |> 
    as.data.frame() |> 
    rownames_to_column(id_var) |> 
    pivot_longer(-{{ id_var }}, names_to = names, values_to = to_what)
}

remove_batch_effects <- function(set, what = "count_norm", names = "sample", id_var = "gene_id",
                                 batch_var = "day", formula = "~ condition + day",
                                 filt = "TRUE") {
  
  file_dat <- tempfile("data")
  file_desc <- tempfile("description")
  
  meta <- set$metadata |> 
    filter(!!rlang::parse_expr(filt)) |> 
    droplevels() |> 
    mutate(batch = get(batch_var)) |> 
    mutate(ibatch = as.integer(batch)) |> 
    arrange(ibatch)
  design_mat <- model.matrix(as.formula(formula), data = meta)
  
  dat <- set$dat |> 
    pivot_wider(id_cols = {{ id_var }}, names_from = !!names, values_from = !!what)
  dat <- dat[, c(id_var, as.character(meta$sample))]
  tab <- dat |> 
    column_to_rownames(id_var) |> 
    as.matrix()
  
  desc <- meta |> 
    rename(ID = sample) |> 
    mutate(sample = row_number(), batch = ibatch) |> 
    select(ID, sample, batch)
  
  write_tsv(dat, file_dat, na = "NaN")
  write_csv(desc, file_desc)
  
  # hr_combat <- HarmonizR::harmonizR(file_dat, file_desc, algorithm = "ComBat")
  hr_limma <- HarmonizR::harmonizR(file_dat, file_desc, algorithm = "limma")
  # bat <- limma::removeBatchEffect(tab, batch = meta[[batch_var]], design = design_mat)
  
  # dat_bat <- mat2dat(bat, "count_batch", names)
  # dat_combat <- mat2dat(hr_combat, "count_combat", names)
  dat_limma <- mat2dat(hr_limma, "count_limma", names)
  set$dat <- set$dat |>   
    # inner_join(dat_combat, by = c(id_var, names)) |> 
    inner_join(dat_limma, by = c(id_var, names)) 
    # inner_join(dat_bat, by = c(id_var, names))
  set$metadata <- meta
  
  set
}


