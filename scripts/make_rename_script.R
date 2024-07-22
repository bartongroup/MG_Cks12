# Create a script called "merge_rename_fastq.sh" that merges lanes and copies
# original FASTQ files to a new directory. The names are simplified.

library(tidyverse)

read_csv("rna_seq/fastq_files.txt", col_names = "fq_file", show_col_types = FALSE) |> 
  filter(str_detect(fq_file, "^X")) |> 
  separate_wider_delim(fq_file, delim = "/", names = c("top_dir", "next_dir", "sample", "file_name"), cols_remove = FALSE) |> 
  mutate(
    lane = str_extract(file_name, "L[0-9]+"),
    pair = str_extract(file_name, "_[12](?=\\.fq\\.gz$)") %>% str_replace_all("_", "")
  ) |> 
  group_by(sample, pair) |> 
  summarise(files = str_c(fq_file, collapse = " "), n = n()) |> 
  ungroup() |> 
  mutate(
    new_name = str_glue("fastq/{sample}_R{pair}.fastq.gz"),
    cmd = if_else(n == 1,
      str_glue("cp {files} {new_name}"),
      str_glue("cat {files} > {new_name}")
    )
  ) |> 
  select(cmd) |> 
  write_csv("merge_rename_fastq.sh", col_names = FALSE)
