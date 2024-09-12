
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))
suppressMessages(library(igraph))
suppressMessages(library(ggpubr))
suppressMessages(library(tidyverse))


section_name <- "10e_find_interesting_target_downstream_pairs"
analysis_name <- "pipeline"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
ExperimentName <- "Magpie"
date <- "2022-08-15" #(magpie)
replogle_date <- '2022-09-12'


args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  ##needed arguments
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--experiment_name"){ ExperimentName <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
  if (arg == "--target_gene"){ target_gene <- args[[ind + 1]] }
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

target_meta <- fread(file.path(OutFolder, "7b_target_summary", paste0(date, "_target_meta_data.csv")))

cell_type_comparison_fnm <- file.path(paste0(outdir, '/', date, '_target_downstream_replogle_comparison.tsv'))
magpie_lfcs <- fread(file.path(OutFolder, "6b_calc_lfcs_transcriptome_wide_by_gene", paste0(date, "_combined_with_adj_pval.tsv.gz")))

if (!(file.exists(cell_type_comparison_fnm))){
  lfcs_by_cell_type <- magpie_lfcs %>%
    dplyr::select(c('target', 'downstream_gene_name', "iPSC_lfc" = 'lfc', "iPSC_pval_adj" = 'pval_adj'))
  
  ## load replogle as well
  for (cell_type in c("K562_essential", "K562_gwps", "RPE1_raw")){
    cell_type_lfcs <- fread(file.path(ProjectFolder, "outs", "Other", "Replogle", "6a_analysis_heatmap", 
                                      paste0(replogle_date, "_", cell_type, "_target_downstream_lfc.tsv.gz")))
    #targets_in_screen <- unique(cell_type_lfcs)
    eval(parse(text = paste0("lfcs_by_cell_type <- left_join(lfcs_by_cell_type, ", 
                             "cell_type_lfcs %>% dplyr::select(", 
                             "c('target', 'downstream_gene_name', '", paste0(cell_type, '_lfc' ), "' = 'lfc', '", paste0(cell_type, '_pval_adj'),"' = 'pval_adj')))"
    )))
  }
  fwrite(lfcs_by_cell_type, cell_type_comparison_fnm, sep ='\t')
} else {
  lfcs_by_cell_type <- fread(cell_type_comparison_fnm)
}
all_targets <- unique(lfcs_by_cell_type$target); all_downstream <- unique(lfcs_by_cell_type$downstream_gene_name); targets_expressed <- intersect(all_targets, all_downstream)
expressed_targets <- intersect(all_targets, all_downstream)

## load external data
gprofiler_annotations <- readRDS(paste0(ResourcesFolder, "gprofiler/gprofiler_full_hsapiens.ENSG.RDS"))
gprofiler_terms <- fread(paste0(OutFolder, "Preprocess_External_Datasets/gprofiler/go_terms_2022-10-10.tsv"))

fnm <- sort(grep(list.files(file.path(OutFolder, "Preprocess_External_Datasets", "dorothea"),
                            pattern = "dorothea_interaction_pairs_"),
                 pattern = paste0(tolower(ExperimentName), "-", date, "-version.tsv"), value = T), decreasing = T)[1]
tf_pairs <-  fread(file.path(OutFolder, "Preprocess_External_Datasets", "dorothea", fnm))
tf_list <- fread("/lustre/scratch123/hgi/mdt2/teams/parts/jh47/claudia_files/DATA/DatabaseExtract_v_1.01.csv")

## ---- Preprocess Data

new_de_pairs <- lfcs_by_cell_type %>%
  dplyr::filter(iPSC_pval_adj < sig_pval_thresh & abs(iPSC_lfc) > sig_abs_lfc_thresh) %>%
  dplyr::mutate(sig_in_K562 = ifelse(!is.na(K562_gwps_pval_adj), ifelse(K562_gwps_pval_adj < 0.5, 'DE', 'not_DE'), 'missing'),
                sig_in_other_cell_types = ifelse(
                  !is.na(K562_essential_pval_adj) & 
                    !is.na(K562_gwps_pval_adj) & !is.na(RPE1_raw_pval_adj),  
                  ifelse(K562_essential_pval_adj > 0.5 & 
                           K562_gwps_pval_adj > 0.5 & 
                           RPE1_raw_pval_adj > 0.5, "not_DE_in_all", "DE_in_at_least_one" ), 'missing')
  )%>%
  dplyr::filter(target != downstream_gene_name & ((sig_in_other_cell_types == "not_DE_in_all") | 
                                                    (sig_in_other_cell_types == "missing" & sig_in_K562 == "not_DE")))
fwrite(new_de_pairs, paste0(outdir, '/', date, '_iPSC_de_pairs.tsv'), sep= '\t')



## ---- SessionInfo()

sessionInfo()

## ---- scratch



