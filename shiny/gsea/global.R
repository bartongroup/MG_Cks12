# GENE SET ENRICHMENT EXPLORER
#
# Global settings and variables
#


library(shiny)
library(shinydashboard)
library(jsonlite)

# datatable needs this option to display Inf correctly

options(htmlwidgets.TOJSON_ARGS = list(na = 'string'))

options(dplyr.summarise.inform = FALSE)
select <- dplyr::select
rename <- dplyr::rename

# Source all R files

files_R <- list.files(c("../R", "modules"), pattern = "*.R$", full.names = TRUE)
sr_ <- sapply(files_R, source)

# Read config json file

config <- jsonlite::read_json("../shiny_config.json")

# Read all data, make it available to all modules

shd <- sh_read_data(c("data", "metadata", "features", "de", "gse"), config, with_progress = FALSE)
dexset <- list(
  de = shd$de,
  data = shd$data,
  metadata = shd$metadata,
  name = "Treatments"
) 
dexset <- list(dexset) |> rlang::set_names(dexset$name)

  # Mutliple test corrections
  for(nm in names(dexset)) {
    dexset[[nm]]$de <- dexset[[nm]]$de |>
      dplyr::group_by(contrast) |>
      dplyr::mutate(fdr = p.adjust(p_value, method = "BH")) |>
      dplyr::ungroup()
  }


DATA <- list(
  dex = dexset,
  names = names(dexset),
  features = shd$features,
  gse = shd$gse
)

DATA$id2name <- set_names(DATA$features$name, DATA$features$id)
DATA$name2id <- set_names(DATA$features$id, DATA$features$name)

