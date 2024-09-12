

library(tidyverse)
library(dplyr)
library(knitr) 
library(ggplot2)
library(data.table)



HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
section_name <- "7b_target_summary"
date <- "2022-10-05"
ExperimentName <- "Pica"

# options
args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
  if (arg == "--experiment_name"){ ExperimentName <- args[[ind + 1]] }
  if (arg == "--padj_th"){ padj_th <- as.numeric(args[[ind + 1]]) }
  if (arg == "--abs_lfc_th"){ abs_lfc_th <- as.numeric(args[[ind + 1]]) }
}

# i/o
source(file.path(HomeFolder, "scripts/io/", paste0(ExperimentName, "_io.R")))
outdir <- file.path(OutFolder, section_name)
plotsdir <- file.path(HTMLFolder, "pipeline/", section_name)
if(!dir.exists(outdir)) dir.create(outdir)
if(!dir.exists(plotsdir)) dir.create(plotsdir)

if (!exists("padj_th")){padj_th <- sig_pval_thresh}
if (!exists("abs_lfc_th")){abs_lfc_th <- sig_abs_lfc_thresh}

## LoadData
MagpieGuideMetadata <- fread(GuideMetadataPath)
lfcs <- fread(file.path(OutFolder, "6b_calc_lfcs_transcriptome_wide_by_gene",
                        paste0(date, "_combined_with_adj_pval.tsv.gz")))
#filter out crispr sequences
crispr_sequences <- gsub(list.files(file.path(ResourcesFolder, 'reference_sequences', 'added_sequences'), pattern = '[.]gtf'), pattern = '.gtf', replacement = '')
lfcs <- dplyr::filter(lfcs, !(downstream_gene_name %in% crispr_sequences))
lfcs <- lfcs %>% 
  mutate(target = gsub(target, pattern = "_", replacement = '-')) %>%
  mutate(is_on_target = (target == downstream_gene_name))


# define parameters for Rmd
params_list <- list(home_folder = HomeFolder,
                    date = date,
                    section_name = section_name,
                    experiment_name = ExperimentName,
                    padj_th = sig_pval_thresh,
                    abs_lfc_th = sig_abs_lfc_thresh)

# run Rmd and create htmls
fnm_html <- paste0(section_name,".html")
rmarkdown::render(paste0(CodeFolder, "/Magpie/pipeline/", section_name, ".Rmd"),
                  envir = new.env(),
                  output_file = file.path(plotsdir, fnm_html),
                  params =  params_list)
