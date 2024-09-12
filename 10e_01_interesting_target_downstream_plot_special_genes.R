
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))
suppressMessages(library(glmnet))
suppressMessages(library(doParallel))
suppressMessages(library(ggpubr))
suppressMessages(library(tidyverse))
suppressMessages(library(doParallel))


section_name <- "10e_interesting_target_downstream"
subsection_name <- "10e_01_interesting_target_downstream_plot_special_genes"
analysis_name <- "pipeline"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
ExperimentName <- "Magpie"
date <- "2022-08-15" #(magpie)


args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  ##needed arguments
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--experiment_name"){ ExperimentName <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
}

# set relevent i/o paths
source(file.path(HomeFolder, paste0("scripts/io/", ExperimentName, "_io.R")))
source(file.path(HomeFolder, "scripts/io/Magpie_Utils.R"))
outdir <- file.path(OutFolder, "", section_name)
plotsdir <- file.path(HTMLFolder, analysis_name, section_name)
if(!dir.exists(outdir)) dir.create(outdir)
if(!dir.exists(plotsdir)) dir.create(plotsdir)
print(section_name)


## ---- LoadData

gprofiler_annotations <- readRDS(paste0(ResourcesFolder, "gprofiler/gprofiler_full_hsapiens.ENSG.RDS"))
gprofiler_terms <- fread(paste0(OutFolder, "Preprocess_External_Datasets/gprofiler/go_terms_2022-10-10.tsv"))

target_go_enrichment <- fread(paste0(OutFolder, "10d_go_term_enrichment/", date, "_target_go_enrichment.tsv.gz"))
target_go_enrichment$pval_adj <- p.adjust(target_go_enrichment$pval, method = 'BH')
downstream_go_enrichment <- fread(paste0(OutFolder, "10d_go_term_enrichment/", date, "_target_go_enrichment.tsv.gz"))
downstream_go_enrichment$pval_adj <- p.adjust(downstream_go_enrichment$pval, method = "BH")
lfcs <- fread(paste0(OutFolder, "6b_calc_lfcs_transcriptome_wide_by_gene/", date, "_combined_with_adj_pval.tsv.gz"))

target_meta <- fread(paste0(OutFolder, "/7b_target_summary/", date, "_target_meta_data.csv"))
downstream_meta <- fread(paste0(OutFolder, "10b_number_of_perturbing_knockdowns/", date, "_downstream_meta_data.csv"))

## ---- Rmarkdown

rmarkdown::render(file.path(CodeFolder, "Magpie", 'pipeline', paste0(subsection_name, ".Rmd")),
                  output_file = paste0(plotsdir, '/', subsection_name, '/', date, ".html"),
                  params = list(
                    date = date,
                    home_folder = HomeFolder,
                    section_name = section_name,
                    analysis_name = subsection_name
                  ))

## ---- SessionInfo()

sessionInfo()
