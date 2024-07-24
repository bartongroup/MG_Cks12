library(dexdash)
library(readr)

# datatable needs this option to display Inf correctly

options(htmlwidgets.TOJSON_ARGS = list(na = 'string'))
options(dplyr.summarise.inform = FALSE)

# Source all R files

files_R <- list.files(c("../R"), pattern = "*.R$", full.names = TRUE)
sr_ <- sapply(files_R, source)


# Read config json file

config <- jsonlite::read_json("../shiny_config.json")

# Read all data, make it available to all modules

shd <- sh_read_data(c("de1", "de2", "data", "metadata", "features", "fterms"), config, with_progress = FALSE)
dexset1 <- dexdash::dexdash_set(shd$de1, shd$data, shd$metadata, name = "Treatments")
dexset2 <- dexdash::dexdash_set(shd$de2, shd$data, shd$metadata, name = "Separate Cks")
dexset <- dexdash::dexdash_list(dexset1, dexset2)
dexdash::run_app(dexset, shd$features, shd$fterms, title = "Impact of Cks1/2 on transcriptome")
