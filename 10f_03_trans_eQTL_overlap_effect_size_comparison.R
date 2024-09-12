
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))
suppressMessages(library(doParallel))
suppressMessages(library(ggpubr))
suppressMessages(library(tidyverse))
library(plotly)

section_name <- "10f_trans_eQTL_overlap"
subsection_name <- "10f_03_trans_eQTL_overlap_effect_size_comparison"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
ExperimentName <- "Magpie"
date <- "2022-08-15" 

args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  ##needed arguments
  if (arg == "--date"){ magpie_version <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--start_ind"){ start_ind <- as.numeric(args[[ind + 1]]) }
  if (arg == "--no_genes2test"){ no_genes2test <- as.numeric(args[[ind + 1]]) }
  
}

# set relevent i/o paths
source(paste0(HomeFolder, "scripts/io/", ExperimentName, "_io.R"))
source(paste0(HomeFolder, "scripts/io/Magpie_Utils.R"))
outdir <- file.path(OutFolder, section_name, subsection_name)
plotsdir <- file.path(HTMLFolder, 'pipeline/', section_name, subsection_name)
if(!dir.exists(outdir)) dir.create(outdir, recursive = T)
if(!dir.exists(plotsdir)) dir.create(plotsdir, recursive = T)
print(section_name)
control_tag <- c("unassigned", NonTargetGeneName) # use all unassigned cells as control


## ---- LoadData

cor_fnm <- paste0(outdir, '/', 'cis_effect_cor_df.tsv')
if (file.exists(cor_fnm)){
  control_cis_effects <- fread(cor_fnm)
} else {
  fnms <- list.files(file.path(OutFolder, section_name, '10f_02_trans_eQTL_overlap_compute_cis_effect_size'), full.names = T)
  control_cis_effects <- as.data.frame(bind_rows(lapply(fnms, readRDS)))
  fwrite(control_cis_effects, , sep = '\t')
}
target_meta <- fread("/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Magpie/7b_target_summary/2022-08-15_target_meta_data.csv")
magpie_lfcs <- fread(paste0(OutFolder, '6b_calc_lfcs_transcriptome_wide_by_gene/', date, '_combined_with_adj_pval.tsv.gz'))
magpie_lfcs_split <- split.data.frame(magpie_lfcs, f = magpie_lfcs$target)

cis_eQTL_meta <- fread("/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Magpie/10f_trans_eQTL_overlap/10f_01_trans_eQTL_overlap_overall_signal/trans_eQTL_overlap.tsv.gz")
effect_cor_df <- fread("/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Magpie/10f_trans_eQTL_overlap/10f_01_trans_eQTL_overlap_overall_signal/trans_effect_info.tsv")
trans_eQTL_list <- fread("/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Magpie/Preprocess_External_Datasets/eQTLs/eQTL_interaction_pairs_magpie-2022-08-15-version.tsv")


expression_variance <- fread("/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Magpie/Preprocess_External_Datasets/coexpression/sc_regressed_variance_per_gene_magpie-2022-08-15-version.tsv")
names(expression_variance) <- c('gene_name', 'sc_var')
expression_variance$gene_name <- unlist(lapply(strsplit(expression_variance$gene_name, ":"), "[[", 2))


## ---- Preprocess Data

trans_eQTL_df <- left_join(control_cis_effects, target_meta %>% dplyr::select('gene_name' = 'gene', 'target_lfc', 'pval_adj'))


## ---- Compare Effect Sizes

### interactive plot
df4plot <- left_join(trans_eQTL_df, expression_variance)
p <- ggplot(df4plot, aes(x = abs(control_beta), y = abs(target_lfc), label = gene_name, col = sc_var)) + 
  xlab("Natural effects") + ylab("CRISPR down-regulation") + 
  geom_point(size = 0.4) + geom_abline(slope = c(1,2), lty = 3, col = 'red') + theme_bw()
ggplotly(p)



## ---- SessionInfo

sessionInfo()


