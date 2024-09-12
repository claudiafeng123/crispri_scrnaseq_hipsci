date <- "2022-08-15"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
ExperimentName <- "Magpie"
section_name <- "10a_number_of_DEGs"

library(Seurat)
library(tidyverse)
library(data.table)
library(reshape2)
library(ggpointdensity)
library(magrittr)


# options
args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
  if (arg == "--rmd_file"){ rmd_file <- args[[ind + 1]] }
}

# i/o
source(file.path(HomeFolder, paste0("scripts/io/", ExperimentName, "_io.R")))
source(file.path(HomeFolder, "scripts/io/Magpie_Utils.R"))
outdir <- file.path(OutFolder, section_name)
plotsdir <- file.path(HTMLFolder, "pipeline/", section_name, '/')
if(!dir.exists(outdir)) dir.create(outdir)
if(!dir.exists(plotsdir)) dir.create(plotsdir)

# define parameters for Rmd
params_list <- list(utils_path = file.path(HomeFolder, "scripts/io/Magpie_Utils.R"),
                    io_path = file.path(HomeFolder, "scripts/io/Magpie_io.R"),
                    date = date,
                    section_name = section_name)
print(params_list)

## ---- LoadData

target_df <- read.csv(file.path(OutFolder, "7b_target_summary", paste0(date, "_target_meta_data.csv")))
target_df <- dplyr::filter(target_df, gene != "NonTarget")

## add in K562 and RPE1


gprofiler_annotations <- readRDS("/lustre/scratch123/hgi/projects/crispri_scrnaseq/resources/gprofiler/gprofiler_full_hsapiens.ENSG.RDS")
gprofiler_terms <- fread("/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Magpie/Preprocess_External_Datasets/gprofiler/go_terms_2023-10-18.tsv")

## ----

# run Rmd and create htmls
fnm_html <- paste0(section_name, ".html")
rmarkdown::render(paste0(CodeFolder, 'Magpie/pipeline/', section_name, '.Rmd'),
                  envir = new.env(),
                  output_file = file.path(plotsdir, fnm_html),
                  params =  params_list)
