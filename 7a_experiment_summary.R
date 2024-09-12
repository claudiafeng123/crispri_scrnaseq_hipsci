
library(data.table)

## set variabless
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
date <- "2022-05-31"

args <- commandArgs(trailingOnly = TRUE)
HomeFolder=args[[1]]
date=args[[2]]
print(date)

section_name <- "7a_experiment_summary"
# set relevent i/o paths
source(paste0(HomeFolder, "scripts/io/Magpie_io.R"))
source(paste0(HomeFolder, "scripts/io/Magpie_Utils.R"))

## read lane metadata
LaneMetadata <- fread(LaneMetadataPath)


#execute R markdown

rmarkdown::render(input = paste0(CodeFolder, ExperimentName, "/pipeline/", section_name, ".Rmd"), 
                  output_file = paste0(HTMLFolder, "pipeline/", section_name, "/", section_name, "_", date, ".html"), 
                  params = list(home_folder = HomeFolder, date = date, section_name = section_name, experiment_name = ExperimentName))

sessionInfo()