# Get "proper" chromosome names, without patches and haplotypes
biomart_fetch_chromosomes <- function(mart) {
  chrnames <- getBM(attributes = "chromosome_name", mart = mart, uniqueRows = TRUE)
  nums <- chrnames |> 
    as_tibble() |> 
    filter(str_detect(chromosome_name, "^\\d{1,2}$")) |>
    pull(chromosome_name) |> 
    as.numeric() |> 
    sort() |> 
    as.character()
  c(nums, "X", "Y", "MT")
}

biomart_fetch_genes <- function(mart, use_cache = TRUE, with_ncbi = FALSE) {
  attr <- c(
    chr = "chromosome_name",
    start = "start_position",
    tss = "transcription_start_site",
    end = "end_position",
    strand = "strand",
    gene_biotype = "gene_biotype",
    gc_content = "percentage_gene_gc_content",
    gene_id = "ensembl_gene_id",
    gene_symbol = "external_gene_name",
    description = "description"
  )
  if(with_ncbi)
    attr <- c(attr, c(ncbi_id = "entrezgene_id"))
  
  getBM(attributes = attr, mart = mart, useCache = use_cache) |> 
  dplyr::rename(!!attr) |>
    dplyr::mutate(
      description = str_remove(description, "\\s\\[.*\\]"),
      gene_symbol = na_if(gene_symbol, "")
    ) |> 
    tibble::as_tibble() |> 
    dplyr::mutate(gene_symbol = if_else(is.na(gene_symbol), gene_id, gene_symbol))
}

biomart_fetch_structures <- function(mart, gene_ids, use_cache = TRUE) {
  attr <-  c(
    chr = "chromosome_name",
    gene_id = "ensembl_gene_id",
    transcript_id = "ensembl_transcript_id",
    exon_id = "ensembl_exon_id",
    gene_symbol = "external_gene_name",
    start = "start_position",
    end = "end_position",
    tss = "transcription_start_site",
    exon_start = "exon_chrom_start",
    exon_end = "exon_chrom_end",
    start_5_utr = "5_utr_start",
    end_5_utr = "5_utr_end",
    start_3_utr = "3_utr_start",
    end_3_utr = "3_utr_end",
    transcript_start = "transcript_start",
    transcript_end = "transcript_end",
    strand = "strand"
  )
  
  getBM(attributes = attr, mart = mart, filters = "ensembl_gene_id", values = gene_ids, useCache = use_cache) |> 
    dplyr::rename(!!attr) |>
    dplyr::mutate(
      gene_symbol = na_if(gene_symbol, "")
    ) |> 
    tibble::as_tibble() |> 
    dplyr::mutate(gene_symbol = if_else(is.na(gene_symbol), gene_id, gene_symbol))
}

biomart_fetch_proteins <- function(mart, gene_ids, use_cache = TRUE) {
  getBM(attributes = c(
    "ensembl_gene_id",
    "external_gene_name",
    "ensembl_peptide_id"
  ), mart = mart, filters = "ensembl_gene_id", values = gene_ids, useCache = use_cache) |> 
    dplyr::rename(
      gene_id = ensembl_gene_id,
      gene_symbol = external_gene_name,
      protein_id = ensembl_peptide_id    
    ) |> 
    tibble::as_tibble() |> 
    filter(protein_id != "")
}

# this is inspired by getGeneLengthAndGCContent function in EDASeq
biomart_fetch_gene_lengths <- function(mart, gene_ids, use_cache = TRUE) {
  exons <- getBM(attributes = c(
    "ensembl_gene_id",
    "ensembl_exon_id",
    "chromosome_name",
    "exon_chrom_start",
    "exon_chrom_end"),
    mart = mart, filters = "ensembl_gene_id", values = gene_ids, useCache = use_cache
  )
  # find all exons, reduce overlap
  coords <- exons |>
    group_split(ensembl_gene_id) |> 
    map_dfr(function(w) {
      IRanges::IRanges(w$exon_chrom_start, w$exon_chrom_end) |>
        IRanges::reduce() |>
        as.data.frame() |> 
        as_tibble() |> 
        mutate(ensembl_gene_id = w$ensembl_gene_id[1], chr = w$chromosome_name[1])
    }) |> 
    rename(gene_id = ensembl_gene_id, length = width)
  coords |> 
    group_by(gene_id, chr) |> 
    summarise(length = sum(length))
}


biomart_fetch_gene_sequences <- function(mart, gids, batch_size = 1000) {
  attr <- c(
    "ensembl_gene_id",
    "external_gene_name",
    "gene_exon_intron"
  )
  
  n_batch <- ceiling(length(gids) / batch_size)
  pb <- progress_bar$new(format = "  Fetching gene sequences [:bar] :percent eta: :eta", total = n_batch)
  
  tibble(
    gene_id = gids,
    batch = ntile(gene_id, n_batch)
  ) |> 
    group_split(batch) |> 
    map(function(w) {
      tb <- biomaRt::getBM(attributes = attr, mart = mart, filters = "ensembl_gene_id", values = w$gene_id) |> 
        dplyr::rename(
          gene_id = ensembl_gene_id,
          gene_symbol = external_gene_name,
          sequence = gene_exon_intron
        ) |>
       tibble::as_tibble()
      pb$tick()
      tb
    }) |> 
    list_rbind()
}

