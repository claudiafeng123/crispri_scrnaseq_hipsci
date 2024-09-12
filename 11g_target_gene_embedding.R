library(data.table)
library(dplyr)
library(reshape2)
library(ggplot2)
library(umap)
library(scales)
library(ggrepel)
library(tidyverse)
library(factoextra)
library("grid")
library("gridExtra")
library(ggpubr)
library(ggtext)


# set relevent i/o paths
section_name <- "11g_target_gene_embedding"
#correlation within the complex

# set relevent i/o paths
date <- "2022-08-15"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
ExperimentName <- "Magpie"

args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  ##needed arguments
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--experiment_name"){ ExperimentName <- args[[ind + 1]] }
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
}

setwd(HomeFolder)
source(paste0(HomeFolder, "scripts/io/", ExperimentName, "_io.R"))
source(paste0(HomeFolder, "/scripts/io/Magpie_Utils.R"))

print(paste("date:", date))
print(paste("home folder:", HomeFolder))
print(paste("Experiment Name:", ExperimentName))


## set io
outdir <- file.path(OutFolder, section_name)
if(!dir.exists(outdir)) dir.create(outdir, recursive = T)
plotsdir <- file.path(HTMLFolder, "pipeline", section_name)
if(!dir.exists(plotsdir)) dir.create(plotsdir, recursive = T)


## ---- Load Data

target_complex_cor <- fread(paste0(OutFolder, "11f_target_complex_cor/", date, "_high_target_complex_cor_ALLY_DAVID_validated.csv"))
target_meta <- fread(file.path(OutFolder, "7b_target_summary", paste0(date, "_target_meta_data.csv")))
fnm <- sort(list.files(file.path(OutFolder, 'Preprocess_External_Datasets', 'omnipath'), full.names = T, pattern = "omnipath_complex_meta"), decreasing = T)[1]
complex_meta <- fread(fnm)

target_target_cor_anno <- fread( paste0(OutFolder, "11c_deconvoluting_signal_from_glm", "/", date, "_target_target_cor_with_adj_pval.tsv.gz"))
high_cor_annotated_fnm <- paste0(outdir, '/', date, '_high_cor_annotated.tsv')
if (file.exists(high_cor_annotated_fnm)){
  high_cor_annotated <- fread(high_cor_annotated_fnm)
} else {
  ## target-target correlation
  ##determine the new pairs
  high_cor_annotated <- target_target_cor_anno %>% 
    dplyr::select(c('target_gene_1' = 'gene_1',
                    'target_gene_2' = 'gene_2',
                    'iPSC_cor' = 'cor_all',
                    'K562_gwps_cor',
                    'K562_essential_cor',
                    "RPE1_raw_cor")) %>%
    dplyr::filter(abs(iPSC_cor) > sig_cor_thresh_target) %>%
    dplyr::filter(target_gene_1 < target_gene_2) %>%
    dplyr::mutate(status = ifelse(
      ## missing in everything
      is.na(K562_gwps_cor) & is.na(K562_essential_cor) & is.na(RPE1_raw_cor), 'Missing in other cell types',
      ## not missing
      ifelse(
        # present in K562 but missing in one or more of the essentiality screens
        !is.na(K562_gwps_cor) & (is.na(K562_essential_cor) | is.na(RPE1_raw_cor)), 
        #"Missing in essentiality screens"
        ifelse(
          # correlated in gwps K562, missing in essentiality
          abs(K562_gwps_cor) > sig_cor_thresh_target, "Correlated in Genome-wide K562 Screen, Missing in Essentiality Screens",
          ifelse(
            # not correlated/borderline
            abs(K562_gwps_cor) > insig_cor_thresh_target, "Uncategorized in Genome-wide K562 Screen, Missing in Essentiality Screen",
            "Not Correlated in K562 Screen, Missing in Essentiality Screens"
          ) 
        ),
        #"Present in RPE1 and K562 screens"
        ifelse(
          abs(K562_gwps_cor) > sig_cor_thresh_target & abs(K562_essential_cor) > sig_cor_thresh_target & abs(RPE1_raw_cor) > 0.3,
          "Correlated in K562 & RPE1 Cells",
          ifelse(
            ## not correlated in anything
            abs(K562_gwps_cor) < insig_cor_thresh_target & abs(K562_essential_cor) < insig_cor_thresh_target & abs(RPE1_raw_cor) < insig_cor_thresh_target,
            "Not correlated in K562 or RPE1 cells", "Uncategorized, but data available for all three cell types"
          )
        )
      )
    ))
  fwrite(high_cor_annotated, high_cor_annotated_fnm, sep='\t')
}

new_correlations <- high_cor_annotated %>% 
  dplyr::filter(status %in% c("Not Correlated in K562 Screen, Missing in Essentiality Screens", 'Missing in other cell types', "Not correlated in K562 or RPE1 cells"))
fwrite(new_correlations, paste0(outdir, '/', date, '_high_cor_annotated_iPSC_specific.tsv'), sep='\t')

## magpie
fnm <- sort(list.files(file.path(OutFolder, '6b_calc_lfcs_transcriptome_wide_by_gene'), full.names = T, pattern = "combined_with_adj_pval"), decreasing = T)
magpie_res <- fread(fnm)
magpie_res_split <- split.data.frame(magpie_res, f = magpie_res$target)

#protein_half_life <- fread('/lustre/scratch123/hgi/projects/crispri_scrnaseq/resources/protein_complexes/hela_half_lives_2023-03-22_download.csv')
off_target_fnm <- sort(list.files(file.path(OutFolder, "/Preprocess_External_Datasets/off_target/"), pattern = 'off_target_pairs_magpie-', full.names = T), decreasing = T)[1]
off_target <- fread(off_target_fnm)
## remove anything that causes signifcant down-regulation of the off-target
off_target <- left_join(off_target %>% dplyr::select(c("target", "off_target")), 
                        magpie_res %>% dplyr::select(c("target", "off_target" = 'downstream_gene_name', 'lfc', 'pval_adj'))) %>%
  dplyr::filter(!(is.na(lfc)) & pval_adj < sig_pval_thresh & lfc < 0) %>% unique()
genes2exclude <- sort(unique(off_target$target, off_target$off_target))

tf_list <- fread('/lustre/scratch123/hgi/mdt2/teams/parts/jh47/claudia_files/DATA/DatabaseExtract_v_1.01.csv')
predicted_tfs <- tf_list %>% dplyr::filter(`Is TF?` == "Yes")

## ----

rmarkdown::render(paste0(CodeFolder, "Magpie/pipeline/", section_name, ".Rmd"),
                  output_file = paste0(plotsdir, '/', date, ".html"),
                  params = list(date= date,
                                home_folder = HomeFolder,
                                section_name = section_name))

#plotly::ggplotly(p)
