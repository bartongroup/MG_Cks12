targets_main <- function() {
  
  get_annotations <- tar_plan(
    mart = new_mart(config),
    chromosomes = biomart_fetch_chromosomes(mart),
    gns = biomart_fetch_genes(mart),
    terms = get_functional_terms(config, gns)
    #gterms = prepare_terms_gsea(terms, star$genes$gene_id),
    #fterms = prepare_terms_fenr(terms, star$genes$gene_id)
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
    
    fig_sample_dist = plot_sample_distributions(star, x_breaks = c(0, 1, 2, 3, 4, 5), x_lim = c(-1, 5), ncol = 4, colour_var = "condition"),
    fig_mean_var = plot_mean_var(star, group_var = "condition"),
    fig_distance_mat = plot_distance_matrix(star),
    fig_clustering = plot_clustering(star, colour_var = "condition"),
    fig_pca = plot_pca(star, colour_var = "condition", shape_var = "experiment")
  )
  
  star <- tar_plan(
    star = read_and_process_star(config, metadata, bm_genes, min_count = 10),
    bg = read_bedgraphs(config, metadata, "50000"),
    
    tab_star_log = star$star_log,
    fig_star_log = plot_star_log(star$star_log, metadata),
    fig_star_log_map = plot_star_log_map(star$star_log, metadata),
    fig_map_count = plot_mapped_count(star),
    
    example_count_file = one_count_file(config, metadata, 1),
    fig_star_sense = plot_star_sense(example_count_file)
  )

  get_environment <- tar_plan(
    tar_target(version_file, file.path(config$data_path, "environment.txt"), format = "file"),
    software_versions = get_software_versions(
      version_file,
      selection = c("fastqc", "fastq-screen", "fastp", "ribodetector", "multiqc", "star", "samtools", "bedtools", "snakemake", "drmaa")
    )
  )

  
  c(
    get_annotations,
    setup_experiment,
    get_environment,
    qc,
    star
  )
}