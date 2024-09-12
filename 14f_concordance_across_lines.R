suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(doParallel))
suppressMessages(library(umap))

# set relevent i/o paths
section_name <- "14f_concordance_across_lines"


# set relevent i/o paths
date <- "2024-05-17"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
ExperimentName <- "Pica"

# options for testing
io_path <- paste0("scripts/io/", ExperimentName, "_io.R")
utils_path <- "scripts/io/Magpie_Utils.R"
REDO <- F


args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  ##needed arguments
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--experiment_name"){ ExperimentName <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--utils_path"){ utils_path <- args[[ind + 1]] }
  if (arg == "--io_path"){ io_path <- args[[ind + 1]] }
  if (arg == "--rnd_seed"){ rnd_seed <- as.numeric(args[[ind + 1]]) }
  if (arg == "--REDO"){ REDO <- as.logical(args[[ind + 1]]) }
}

print(paste0(section_name))
print(paste("date:", date))
print(paste("home folder:", HomeFolder))
if(exists("donors2include")){print('including donors: '); print(donors2include)}



setwd(HomeFolder)
source(io_path)
source(utils_path)

GuideMetadata <- fread(GuideMetadataPath)
LineMetadata <- fread(LineMetadataPath)

outdir <- file.path(file.path(OutFolder, section_name))
if(!dir.exists(outdir)) dir.create(outdir, recursive = T)
plotsdir <- file.path(file.path(HTMLFolder, 'pipeline', section_name))
if(!dir.exists(file.path(plotsdir))) dir.create(file.path(plotsdir), recursive = T)

## ---- Choose Genes

fnm <- sort(list.files(file.path(OutFolder, '7b_target_summary'), pattern = '_target_meta_data.csv'), decreasing = T )[1]
pica_target_meta <- fread(file.path(OutFolder, '7b_target_summary', fnm))
pica_lfcs_by_line <- fread("/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Pica/6f_calc_lfcs_transcriptome_wide_by_gene_per_line/2024-05-17_all_lines_with_adj_pval.tsv.gz")
pica_lfcs_by_line_split <- split.data.frame(pica_lfcs_by_line, f = pica_lfcs_by_line$cell_line)


## ---- Correlation of Trans effects

df4cor <- pica_lfcs_by_line %>%
  mutate(to_match = paste0(target, ":", downstream_gene_name)) %>%
  dplyr::select(c("to_match", 'lfc', 'cell_line')) %>%
  reshape2::dcast(to_match ~ cell_line, value.var = 'lfc') %>%
  column_to_rownames('to_match')
cor_across_lines <- cor(df4cor, use="complete.obs")
fwrite(cor_across_lines, paste0(outdir, '/', date, "_cor_by_line_all_trans.tsv"), sep = '\t')
my_heatmap(cor_across_lines, 
           treeheight_row = 0, treeheight_col = 0,
           border_col = NA,
           min_c = -0.6, max_c = 0.6)

