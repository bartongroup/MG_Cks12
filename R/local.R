top_expressed_genes <- function(set, gns, n_top = 10) {
  genes <- gns |> 
    select(chr, start, end, gene_id, gene_symbol, description) |> 
    distinct()
  
  mn <- star_mean_expr(set, what = "count_norm")
  mn |> 
    arrange(-M) |> 
    group_split(condition) |> 
    map_dfr(\(w) head(w, n_top)) |> 
    left_join(genes, by = "gene_id")
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


read_fastq_lengths <- function(meta) {
  f_trim <- read_fastq_length("rna_seq/fastq_trimmed", meta) |> 
    add_column(type = "Trimmed")
  f_clean <- read_fastq_length("rna_seq/fastq_clean", meta) |> 
    add_column(type = "Clean")
  bind_rows(f_trim, f_clean)
}


plot_ribo_proportion <- function(fl) {
  fl |> 
    pivot_wider(id_cols = sample, names_from = type, values_from = read_count) |> 
    mutate(prop = 1 - Clean / Trimmed) |> 
  ggplot(aes(x = sample, y = prop)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    ) +
    geom_col() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    labs(x = NULL, y = "Ribosomal proportion")
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



get_clock_genes <- function(cl_file) {
  readxl::read_excel(cl_file) |> 
    janitor::clean_names() |> 
    select(gene_symbol = gene_name, class) |> 
    mutate(gene_symbol = case_match(
      gene_symbol,
      "KIAA0922" ~ "TMEM131L",
      "FAM122A" ~ "PABIR1",
      .default = gene_symbol
    )) |> 
    mutate(class = tolower(class) |> fct_relevel("phase", "anti-phase", "phase (ambiguous)"))
}


plot_significant_clock_genes <- function(de, clock_genes) {
  d <- de |> 
    filter(str_detect(contrast, "YTHDF")) |> 
    mutate(contrast = str_remove(contrast, "condition")) |> 
    inner_join(clock_genes, by = join_by(gene_id, gene_symbol)) |> 
    select(gene_symbol, logFC, logCPM, PValue, FDR, contrast, class, sig) |> 
    group_by(gene_symbol) |> 
    mutate(
      any_sig = any(sig),
      mean_logFC = mean(logFC)
    ) |> 
    ungroup() |> 
    filter(any_sig)
  
  mx <- max(abs(d$logFC))
    
  classes <- levels(d$class)
  sizes <- d |> arrange(class) |> dplyr::count(class) |> pull(n)
  pl <- map(classes, function(cl) {
    d |>
      filter(class == cl) |> 
      ggplot(aes(x = logFC, y = fct_reorder(gene_symbol, mean_logFC), fill = sig)) +
      th +
      theme(legend.position = "none") +
      geom_col() +
      scale_x_continuous(limits = c(-mx, mx)) +
      scale_fill_manual(values = c("grey70", "black")) +
      facet_wrap(~ contrast, nrow = 1) +
      labs(x = expression(log[2]~FC), y = NULL, title = cl)
  })
  plot_grid(plotlist = pl, ncol = 1, align = "v", rel_heights = sizes + 25)

}


prepare_clock_genes_for_gsea <- function(clock_genes) {
  mapping <- clock_genes |> 
    select(gene_id, term_id = class)
  terms <- clock_genes |> 
    select(term_id = class) |> 
    distinct() |> 
    mutate(term_name = term_id)
  list(
    clock = prepare_for_gsea(terms, mapping)
  )
}

de_clock_genes <- function(de, clock_genes) {
  de |> 
    filter(contrast %in% c("conditionYTHDF1", "conditionYTHDF2")) |> 
    mutate(contrast = str_remove(contrast, "condition")) |> 
    inner_join(clock_genes, by = join_by(gene_id, gene_symbol)) |> 
    select(gene_id, logFC, logCPM, PValue, FDR, contrast, class, sig) |> 
    group_by(gene_id) |> 
    mutate(
      any_sig = any(sig),
      mean_logFC = mean(logFC)
    ) |> 
    ungroup() |> 
    filter(any_sig) |> 
    select(gene_id, class, mean_logFC) |> 
    distinct() |> 
    arrange(class, -mean_logFC)
}


plot_volcano_clock <- function(de, clock_genes) {
  d <- de |> 
    filter(str_detect(contrast, "condition")) |> 
    left_join(clock_genes, by = join_by(gene_symbol, gene_id)) |> 
    mutate(
      contrast = str_remove(contrast, "condition"),
      sel = !is.na(class)
    )
    
  classes <- levels(clock_genes$class)
  pl <- map(classes, function(cl) {
    d_sel <- d |> filter(class == cl)
    
    ggplot(mapping = aes(x = logFC, y = -log10(PValue))) +
      th +
      geom_point(data = d, colour = "grey80", size = 0.6, alpha = 0.3) +
      geom_point(data = d_sel, colour = "black", size = 1.3) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
      geom_vline(xintercept = 0, linewidth = 0.1, alpha = 0.5) +
      facet_wrap( ~ contrast) +
      labs(x = expression(log[2]~FC), y = expression(-log[10]~P), title = cl)
  })
  
  plot_grid(plotlist = pl, ncol = 1, align = "v")
}