targets_main <- function() {
  
  get_annotations <- tar_plan(
    mart = new_mart(config),
    chromosomes = biomart_fetch_chromosomes(mart),
    gns = biomart_fetch_genes(mart),
    terms = get_functional_terms(config, gns),
    fterms = prepare_terms_fenr(terms, star$genes$id)
  )
  
  setup_experiment <- tar_plan(
    tar_target(config_file, "config.json", format = "file"),
    config = jsonlite::read_json(config_file),
    metadata = make_metadata(config$metadata_file)
  )
  
  qc <- tar_plan(
    fastq_len = read_fastq_lengths(config, metadata),
    fig_ribo_prop = plot_ribo_proportion(fastq_len),
    
    fscrn = parse_fscreens(config, metadata),
    qcs = parse_qcs(config, metadata, paired = TRUE),
    idxstats = parse_idxstats(config, metadata),
    gene2name = rlang::set_names(star$genes$gene_symbol, star$genes$gene_id),
    scaffolds = idxstats |> filter(count > 1000 & !(chr %in%  chromosomes)) |> pull(chr) |> unique(),
    
    fig_fscreen = plot_fscreen_map(fscrn),
    fig_read_qual = plot_qualities(qcs),
    fig_read_qual_clust = plot_cluster_qualities(qcs),
    fig_chrom_proportion = plot_chrom_proportion(idxstats, chromosomes, per_length = TRUE),
    fig_chrom_count = plot_chrom_proportion(idxstats, chromosomes, per_length = FALSE),
    fig_scaffold_count = plot_chrom_proportion(idxstats, scaffolds, per_length = FALSE),
    
    png_sample_dist = plot_sample_quasirandom(star, colour_var = "experiment") |> gs("sample_quasirandom", 8, 3.5),
    fig_mean_var = plot_mean_var(star, group_var = "experiment"),
    fig_distance_mat = plot_distance_matrix(star),
    fig_clustering = plot_clustering(star, colour_var = "experiment"),
    fig_pca = plot_pca(star, colour_var = "experiment", shape_var = "treatment")
  )
  
  star <- tar_plan(
    star = read_and_process_star(config, metadata, gns, min_count = 0) |> 
      filter_per_experiment(min_count = 10),
    #bg = read_bedgraphs(config, metadata),
    
    tab_star_log = star$star_log,
    fig_star_log = plot_star_log(star, group_var = "experiment"),
    fig_star_log_map = plot_star_log_map(star),
    fig_map_count = plot_mapped_count(star, nrow = 2),
    
    example_count_file = one_count_file(config, metadata, 1),
    fig_star_sense = plot_star_sense(example_count_file),

    all_ids = star$dat |> filter(!bad) |> pull(id) |> unique()
  )

  de <- tar_plan(
    def1 = edger_de_formula(star, "~ treatment"),
    def2 = edger_de_formula(star, "~ siCks1 + siCks2"),
    
    figs_def1 = plot_vmp(def1, fdr_limit = config$fdr_limit),
    figs_def2 = plot_vmp(def2, fdr_limit = config$fdr_limit),

    def1_ids = select_de_ids(def1, config$fdr_limit),
    def2_ids = select_de_ids(def2, config$fdr_limit)
  )

  get_environment <- tar_plan(
    tar_target(version_file, file.path(config$data_path, "environment.txt"), format = "file"),
    software_versions = get_software_versions(
      version_file,
      selection = c("fastqc", "fastq-screen", "fastp", "ribodetector", "multiqc", "star", "samtools", "bedtools", "snakemake", "drmaa")
    )
  )

  gene_set_enrichment <- tar_plan(
    gse1 = fgsea_all_terms(def1, filt = "logCPM > 0", fterms = fterms),
    gse2 = fgsea_all_terms(def2, filt = "logCPM > 0", fterms = fterms),
    gse_example = list(ontology = "go", term_id = "GO:0007059", contrast = "treatmentsiCks1.2"),
    gse_random = list(ontology = "bp", term_id = "bioplanet_1340", contrast = "treatmentsiCks1.2"),
    gse_example_stats = make_gse_example(def1, terms, gse1, gse_example),
    gse_random_stats = make_gse_example(def1, terms, gse1, gse_random),
    
    fig_fg_example = plot_fgsea_enrichment(def1, fterms, gse_example),
    fig_fg_example_random = plot_fgsea_enrichment(def1, fterms, gse_random)
  )

  for_report <- tar_plan(
    cks_ids = gns |> filter(gene_symbol %in% c("CKS1B", "CKS2")) |> pull(id) |> unique(),
    fig_cks = plot_gene_groups(star, cks_ids, "count_norm", group_var = "treatment")
  )
  
  for_shiny <- tar_plan(
    sav_shiny = save_data_for_shiny(star, def1, def2, terms, fterms, gse1, gse2)
  )

  downloads <- tar_plan(
    tbl_cnt = write_counts(star, "tab/normalized_counts.csv"),
    tbl_def1 = write_de(def1, star$genes, "tab/de_treatment.csv"),
    tbl_def2 = write_de(def2, star$genes, "tab/de_separate.csv"),
    tbl_gse = write_gse(gse1, gns, "tab/gse.csv")
  )

  c(
    get_annotations,
    setup_experiment,
    get_environment,
    qc,
    star,
    de,
    gene_set_enrichment,
    for_report,
    for_shiny,
    downloads
  )
}