#xvfb-run -a /software/R-4.1.3/bin/R

HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
section_name <- "11d_get_clusters_to_aggregate"
date <- "2022-08-15"
ExperimentName <- "Magpie"

library(data.table)
library(dplyr)
library(doParallel)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggExtra)
library(ggpubr)
library(tidyverse)

# options
args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
  if (arg == "--experiment_name"){ ExperimentName <- args[[ind + 1]] }
}

# i/o
source(file.path(HomeFolder, paste0("scripts/io/", ExperimentName, "_io.R")))
source(file.path(HomeFolder, "scripts/io/Magpie_Utils.R"))
outdir <- file.path(OutFolder, section_name)
plotsdir <- file.path(HTMLFolder, "pipeline/", section_name)
if(!dir.exists(outdir)) dir.create(outdir, recursive = T)
if(!dir.exists(plotsdir)) dir.create(plotsdir, recursive = T)

## ---- LoadData

target_meta <- fread(file.path(OutFolder, "7b_target_summary", paste0(date, "_target_meta_data.csv")))
target_target_cor <- fread(file.path(OutFolder, "11a_target_target_cor", paste0(date, "_target_target_cor.tsv.gz")))

## ---- Visualize

rmarkdown::render(
  file.path(CodeFolder, ExperimentName, 'pipeline', paste0(section_name, '.Rmd')),
  output_file = file.path(plotsdir, paste0(section_name, "_", date, ".html")),
  params = list(
    date = date,
    home_folder = HomeFolder,
    section_name = section_name
  )
)

## ---- scratch
