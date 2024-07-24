


fgsea_run <- function(term_data, res, min_size = 3, n_perm = 1000) {
  res <- res |>
    filter(!is.na(value) & !is.na(feature_id)) |> 
    mutate(first_id = str_remove(feature_id, ";.*$"))
  ranks <-  set_names(res$value, res$first_id)
  fgsea::fgsea(pathways = term_data$term2feature, stats = ranks, nproc = 6, minSize = min_size, eps = 0, nPermSimple = n_perm) |>
    as_tibble() |> 
    rename(term_id = pathway) |> 
    mutate(term_name = term_data$term2name[term_id]) |> 
    arrange(NES) |>
    select(term_id, term_name, p_value = pval, fdr = padj, nes = NES, size, leading_edge = leadingEdge)
}

fgsea_groups <- function(d, term_data, feature_var, value_var, group_var, min_size = 3, n_perm = 1000) {
  d |>
    mutate(value = get(value_var), feature_id = get(feature_var)) |>
    group_split(!!sym(group_var)) |>
    map_dfr(function(w) {
      fgsea_run(term_data, w, min_size, n_perm) |>
        mutate(!!group_var := dplyr::first(w[[group_var]]))
    })
}

fgsea_all_terms <- function(d, terms, feature_var = "id", value_var = "logFC", group_var = "contrast",
                            min_size = 3, n_perm = 1000) {
  ontologies <- names(terms)
  map(ontologies, function(ont) {
    cat(str_glue("  Computing fgsea for {ont}\n\n"))
    fgsea_groups(d, terms[[ont]], feature_var, value_var, group_var, min_size, n_perm)
  }) |>
    set_names(ontologies)
}





plot_fgsea_enrichment <- function(de, gterms, info, value_var = "logFC") {
  de <- de |> filter(contrast == info$contrast)
  lst <- gterms[[info$ontology]]$term2feature[[info$term_id]]
  rnks <- set_names(de[[value_var]], de$id)
  fgsea::plotEnrichment(lst, rnks)
}

split_genes_fgsea <- function(se, fg, groupvar = "contrast") {
  fg |> 
    filter(padj < 0.05) |> 
    group_split(term, !!sym(groupvar)) |> 
    map_dfr(function(w) {
      term <- as.character(w$term)
      gr <- as.character(w[[groupvar]])
      genes <- w$leading_edge[[1]]
      se |> 
        filter(id %in% genes & !!sym(groupvar) == gr) |> 
        mutate(term_id = term, .before = "id")
    })
}


gsea_de <- function(gse, res, fdr_limit = 0.01) {
  contrasts <- unique(res$contrast) |> 
    as.character()
  ontologies <- names(gse)
  
  map_dfr(ontologies, function(ont) {
    map_dfr(contrasts, function(ctr) {
      res_sel <- res |> 
        filter(contrast == ctr)
      gse_sel <-  gse[[ont]] |> 
        filter(contrast == ctr)
      sig_genes <- res_sel |> 
        filter(sig) |> 
        pull(id)
      
      gse_sel |> 
        filter(fdr < fdr_limit) |> 
        unnest(leading_edge) |> 
        filter(leading_edge %in% sig_genes) |> 
        rename(id = leading_edge) |> 
        left_join(res_sel, by = "id") |> 
        select(term_id, term_name, nes, fdr, id, gene_symbol, logFC, logCPM, PValue, FDR) |> 
        add_column(ontology = ont, .before = 1) |> 
        add_column(contrast = ctr)
    })
  })
}


get_terms_str <- function(gso, query, fdr_limit = 0.05) {
  gso |> 
    filter(str_detect(term_name, query) & padj < fdr_limit) |>
    pull(term_id)
}



gsa_text_selection <- function(gse, texts, fdr_limit = 0.05) {
  ontologies <- names(gse)
  map(ontologies, function(ont) {
    g <- gse[[ont]]
    map(texts, function(txt) {
      g |> 
        filter(str_detect(term_name, regex(txt, ignore_case = TRUE)) & fdr < fdr_limit) |> 
        add_column(search = txt)
    }) |> 
      list_rbind()
  }) |> 
    list_rbind()
}



make_gse_example <- function(de, terms, gse, example) {
  genes_used <- de |>
    pull(id) |>
    unique()
  n_genes_used <- length(genes_used)
  n_genes_in_term <- terms[[example$ontology]]$mapping |>
    filter(term_id == example$term_id) |>
    filter(id %in% genes_used) |> 
    nrow()
  n_leading_edge <- gse[[example$ontology]] |>
    filter(term_id == example$term_id & contrast == example$contrast) |>
    pull(leading_edge) |>
    unlist() |>
    length()
  term_name <- terms[[example$ontology]]$terms |> 
    filter(term_id == example$term_id) |> 
    pull(term_name)
  list(
    term_id = example$term_id,
    term_name = term_name,
    contrast = example$contrast,
    n_genes_used = n_genes_used,
    n_genes_in_term = n_genes_in_term,
    n_leading_edge = n_leading_edge
  )
}

write_gse <- function(gse, gene_info, file, n_top = 10, fdr_limit = 0.01) {
  gene_info <- gene_info |> 
    select(id, gene_symbol) |> 
    distinct()

  list_rbind(gse, names_to = "ontology") |> 
    filter(fdr < fdr_limit) |> 
    rename(id = leading_edge) |> 
    unnest(id) |> 
    left_join(gene_info, by = join_by(id)) |> 
    select(-id) |> 
    group_by(across(-gene_symbol)) |> 
    summarise(leading_edge_top_10 = str_c(head(gene_symbol, n_top), collapse = ", ")) |> 
    arrange(nes) |> 
    mutate(across(where(is.numeric), ~signif(.x, 4))) |> 
    write_csv(file)
}
