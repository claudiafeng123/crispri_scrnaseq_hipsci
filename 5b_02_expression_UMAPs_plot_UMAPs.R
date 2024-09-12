#xvfb-run -a /software/R-4.1.3/bin/R
date <- "2022-08-15"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
section_name <- "5b_expression_UMAPs"
subsection_name <- "5b_02_expression_UMAPs_plot_UMAPs"
correction <- "none"
ExperimentName <- "Pica"
## assigned - do for "Magpie", "Pica" and "both"
## control/unassigned - do for "Pica" and "both"
#correction <- "technical" # optionas are "none", "technical" and "all"
#
perturbation_status <- "control" # options are: "perturbed", "control", "unassigned", "assigned" (control + perturbed), "unperturbed" (control + unassigned), "all"
#magpie_date <- "2022-08-15"; pica_date <- '2022-10-05'


library(ggpubr)
library(data.table)
library(Seurat)
library(dplyr)
library(reshape2)
library(ggplot2)
library(doParallel)
library(ggrastr)
library(ggpointdensity)
library(tidyverse)
library(viridis)

# options
args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
  if (arg == "--correction"){ correction <- args[[ind + 1]] }
  if (arg == "--experiment_name"){ ExperimentName <- args[[ind + 1]] }
  if (arg == "--perturbation_status"){ perturbation_status <- args[[ind + 1]] }
}



# i/o
print(paste0("Experiment name: ", ExperimentName))
print(paste0("Perturbation status: ", perturbation_status))
print(paste0("Correction: ", correction))
og_experiment_name <- ExperimentName


ExperimentName <- og_experiment_name
if (ExperimentName == 'both'){
  ExperimentName <- c("Magpie", "Pica")
}

og_perturbation_status <- perturbation_status
if (perturbation_status == 'control'){
  perturbation_status = c('nontarget', 'control')
} else if (perturbation_status == 'assigned'){
  perturbation_status = c("perturbed","nontarget", "control")
} else if (perturbation_status == 'unperturbed'){
  perturbation_status = c("unassigned", "control", "nontarget")
} else if (perturbation_status == 'all'){
  perturbation_status = c("unassigned", "control", "nontarget", "perturbed")
}

config <-paste0( og_experiment_name, "_combined_", og_perturbation_status, "_cells_correction_", correction)

# i/o
print(paste0("date: ", date))
source(file.path(HomeFolder, "scripts/io/", paste0(ifelse(og_experiment_name == "Magpie", "Magpie", "Pica"), "_io.R")))
outdir <- file.path(OutFolder, section_name)
plotsdir <- file.path(HTMLFolder, "pipeline", section_name, config)
if(!dir.exists(outdir)) dir.create(outdir)
if(!dir.exists(plotsdir)) dir.create(plotsdir)

## ---- Load Data

inlets_combined <- readRDS(paste0(outdir, '/', date, '_', config, '.RDS'))


## ---- General UMAPs
# take a look at what we have before going into more detail

df4umap <- data.frame(
  id = row.names(inlets_combined@meta.data),
  UMAP_1 = Embeddings(inlets_combined, reduction = 'umap')[,1],
  UMAP_2 = Embeddings(inlets_combined, reduction = 'umap')[,2]
)
df4umap <- left_join(df4umap,
                     inlets_combined@meta.data %>% rownames_to_column('id'))
fwrite(df4umap, paste0(outdir, '/', date, '_', config, '_umap_coords.tsv.gz'), sep = '\t', compress = 'gzip')

set.seed(0)
sample_inds <- sample(1:dim(df4umap)[1], min(dim(df4umap)[1], 100000), replace = F)
df4plot <- df4umap[sample_inds,]

# run Rmd and create htmls
rmarkdown::render(file.path(CodeFolder, "Magpie", "pipeline", paste0(subsection_name, ".Rmd")),
                  envir = new.env(),
                  output_file = paste0(plotsdir, "/", date, '_', config, ".html"),
                  params =  list(correction = correction,
                                 experiment_name = og_experiment_name,
                                 perturbation_status = og_perturbation_status,
                                 section_name = section_name,
                                 subsection_name = subsection_name))

sessionInfo()