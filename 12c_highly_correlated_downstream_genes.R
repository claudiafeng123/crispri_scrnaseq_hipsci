#xvfb-run -a /software/R-4.1.3/bin/R

suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(reshape2))
suppressMessages(library(ggpointdensity))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
suppressMessages(library(lattice))
suppressMessages(library(igraph))

## get protein complexes

section_name <- "12c_highly_correlated_downstream_genes"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
date <- "2022-08-15"
# can probably make a verynice heatmap

args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  ##needed arguments
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--experiment_name"){ ExperimentName <- args[[ind + 1]] }
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
  
  ### not needed
  if (arg == "--utils_path"){ utils_path <- args[[ind + 1]] }
  if (arg == "--io_path"){ io_path <- args[[ind + 1]] }
  if (arg == "--REDO"){ REDO <- as.logical(args[[ind + 1]]) }
  
}

# set relevent i/o paths
# i/o
source(file.path(HomeFolder, "scripts/io/Magpie_io.R"))
source(file.path(HomeFolder, "scripts/io/Magpie_Utils.R"))
outdir <- file.path(OutFolder, section_name)
plotsdir <- file.path(HTMLFolder, "pipeline/", section_name)
if(!dir.exists(outdir)) dir.create(outdir)
if(!dir.exists(plotsdir)) dir.create(plotsdir)

## ---- LoadData

GuideMetadata <- fread(GuideMetadataPath)

target_meta <- fread(file.path(OutFolder, "7b_target_summary", paste0(date, "_target_meta_data.csv")))
downstream_downstream_cor_anno <- fread(paste0(OutFolder, "12a_downstream_gene_embedding", "/", date, "_downstream_downstream_cor_with_adj_pval.tsv.gz"))
#downstream_genes <- unique(downstream_downstream_cor$downstream_gene_1)
#target_target_cor <- fread("/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Magpie/11d_target_target_cor/2022-08-15_high_correlations_all_conf.tsv")

## target lfcs
magpie_lfcs <- fread(file.path(OutFolder, "6b_calc_lfcs_transcriptome_wide_by_gene", paste0(date, "_combined_with_adj_pval.tsv.gz")))
magpie_lfcs_by_target <- split.data.frame(magpie_lfcs, f = magpie_lfcs$target)
magpie_lfcs_by_downstream <- split.data.frame(magpie_lfcs, f = magpie_lfcs$downstream_gene_name)
all_downstream_genes <- names(magpie_lfcs_by_downstream)

#umap_coords <- fread(paste0(OutFolder, '12a_downstream_gene_embedding/', date, '_min_deg-1_min_targets-1_min_cells-10/', date, '_umap_coords.tsv'))
data4umap <- readRDS('/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Magpie/12a_downstream_gene_embedding/2022-08-15_data4umap.RDS')
umap_coords <- data4umap$umap_coords

## sc coexpression
sc_coexpression <- fread(file.path(OutFolder, "Preprocess_External_Datasets", "coexpression", paste0("sc_regressed_correlation_pairs_magpie-", date, "-version.tsv")))
sc_coexpression$to_match <- paste0(sc_coexpression$gene_1, "_", sc_coexpression$gene_2)

## complexes
omnipath_complexes_pairs <- fread("/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Magpie/Preprocess_External_Datasets/omnipath/omnipath_protein_complex_pairs_data-from-2022-09-07_magpie-2022-08-15-version.tsv")
omnipath_complex_meta_fnm <- sort(grep(list.files(paste0(OutFolder, "Preprocess_External_Datasets/omnipath/"), pattern = "complex_meta"), pattern = date, value = T), decreasing  = T)[1]
omnipath_complex_meta <- fread(paste0(OutFolder, "Preprocess_External_Datasets/omnipath/", omnipath_complex_meta_fnm))


pathways <- readRDS(paste0(ResourcesFolder, "/gene_sets/msigdb.v7.4.symbols.RDS"))

##tf
tfs_fnm <- sort(list.files(file.path(OutFolder, "Preprocess_External_Datasets", "dorothea"), pattern = "dorothea_interaction_pairs"), decreasing = T)[1]
tfs <- fread(file.path(OutFolder, "Preprocess_External_Datasets", "dorothea", tfs_fnm))
#tfs <- filter(tfs, target %in% all_downstream_genes)

off_targets <- fread(paste0(OutFolder, "Preprocess_External_Datasets/off_target/off_target_pairs_magpie-", date, "-version.tsv.gz"), header = T) %>%
  filter(target_gene != off_target_gene)
genes2exclude <- sort(unique(c(off_targets$target_gene, off_targets$off_target_gene)))

gprofiler_annotations <- readRDS("/lustre/scratch123/hgi/projects/crispri_scrnaseq/resources/gprofiler/gprofiler_full_hsapiens.ENSG.RDS")
gprofiler_terms <- fread("/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Magpie/Preprocess_External_Datasets/gprofiler/go_terms_2023-10-18.tsv")



## ---- Correlated Downstream Genes
## heatmaps

params_list <- list(
  date = date,
  home_folder = HomeFolder,
  section_name = section_name
)
rmarkdown::render(file.path(CodeFolder, ExperimentName, "pipeline", paste0(section_name, ".Rmd")),
                  output_file = file.path(plotsdir, paste0(date, "_highly_correlated_downstream_genes.html")),
                  params = params_list)


## ---- Double_Heatmap
# heatmap of single-cell co-expression + downstream-downstream correlation

if (T){
  ## just anything with signal
  dg_with_signal <- lapply(magpie_lfcs_by_downstream, FUN = function(df){
    dim(df %>% filter(pval_adj < sig_pval_thresh))[1]
  }) %>% unlist()
  dg_with_signal <- names(dg_with_signal)[which(dg_with_signal > 10)]
  dg_with_signal <- dg_with_signal[which(dg_with_signal %in% sc_coexpression$gene_1)]
  df4heatmap_downstream_downstream_cor <- magpie_lfcs_by_downstream[dg_with_signal] %>%
    bind_rows() %>% 
    as.data.frame() %>%
    dplyr::select(c("target", "downstream_gene_name", 'lfc')) %>%
    reshape2::dcast(target ~ downstream_gene_name, value.var = 'lfc') %>%
    column_to_rownames('target') %>%
    cor()
  df4heatmap_sc_coexpression <- sc_coexpression %>%
    dplyr::filter(gene_1 %in% dg_with_signal & gene_2 %in% dg_with_signal) %>%
    dplyr::select(c("gene_1", "gene_2", "correlation")) %>%
    reshape2::dcast(gene_1 ~ gene_2, value.var = 'correlation') %>%
    column_to_rownames("gene_1")
  
  pdf(paste0(plotsdir, "/", date, '_downstream_heatmap.pdf'), width = 40, height = 40)
  print(double_heatmap(cor_mat1 = df4heatmap_downstream_downstream_cor, 
                       cor_mat2 = df4heatmap_sc_coexpression, min_c = -0.7, max_c = 0.7))
  dev.off()
  
  
  ## slightly more refined
  targets_of_interest <- list(
    "Cytosolic ribosome" = grep(names(magpie_lfcs_by_downstream), pattern = "^RPL|^RPS", value = T)[which(!(grep(names(magpie_lfcs_by_downstream), pattern = "^RPL|^RPS", value = T) %in% c(c("RPL4", "RPL36AL", "RPL22L1", "RPS6KA1", "RPL39L", "RPS4Y1", "RPS27L", "RPS19BP1", "RPL26L1", "RPS6KB2", "RPS6KC1", "RPL7L1"))))],
    "Electron transport chain" = c("MT-CO2", "MT-CO3", "MT-CYB", "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L"),
    "Apoptosis via the p53 pathway" = c('IKBIP', "RRM2B", "TNFRSF10B", "BBC3", "MDM2", "BAX", "TM7SF3", "TIGAR", "TRAF4"),
    "Amino acid metabolic process" = c("EIF2S2", "YARS", "BCAT1", "MARS", 'SARS', "AARS", "SLC7A11", "SLC1A5", "CARS", "GARS", "PYCR1", "SESN2", 'SHMT2', 'EIF4EBP1', "SLC3A2", 'CHAC1',  'PSAT1', 'PCK2', 'STC2', 'MTHFD2'),
    "Calcium regulation" = c("HSP90B1", "PDIA3", "CALR", "HSPA5"),
    "Cholesterol biosynthesis" = c("MSMO1", "MVD", "HMGCS", "FDFT", "ACLY", "FDPS", "HMGCP", "CYP51A1", "INSIG1", "LSS", "DHCR7", "HSD17B7"),
    "Glycolysis" = c("ALDOA", "TPI1", "GPI", "PKM", "ENO1", "PGK1", "PGAM1", "LDHA"),
    "Differentiation to ectoderm" = c("UTF1", "MDK", "CCND1", "TMSB4X", "DSP", "RBP1", "IFITM1", "IFITM3", "NEFL", "NMU", "CRABP1")
  )
  gene_anno <- data.frame(
    gene = unlist(targets_of_interest),
    pathway = rep(names(targets_of_interest), unlist(lapply(targets_of_interest, length)))
  )
  row.names(gene_anno) <- gene_anno$gene
  gene_anno <- gene_anno %>%
    filter(gene %in% sc_coexpression$gene_1)
  gene_anno$gene <- NULL
  
  
  df4heatmap_downstream_downstream_cor <- magpie_lfcs_by_downstream[unlist(targets_of_interest)] %>%
    bind_rows() %>% 
    as.data.frame() %>%
    dplyr::select(c("target", "downstream_gene_name", 'lfc')) %>%
    reshape2::dcast(target ~ downstream_gene_name, value.var = 'lfc') %>%
    column_to_rownames('target') %>%
    cor()
  df4heatmap_downstream_downstream_cor <- df4heatmap_downstream_downstream_cor[row.names(gene_anno), row.names(gene_anno)]
  
  df4heatmap_sc_coexpression <- sc_coexpression %>%
    dplyr::filter(gene_1 %in% row.names(gene_anno) & gene_2 %in% row.names(gene_anno)) %>%
    dplyr::select(c("gene_1", "gene_2", "correlation")) %>%
    reshape2::dcast(gene_1 ~ gene_2, value.var = 'correlation') %>%
    column_to_rownames("gene_1")
  df4heatmap_sc_coexpression <- df4heatmap_sc_coexpression[row.names(gene_anno), row.names(gene_anno)]
  
  data4figure <- list(
    df4heatmap_downstream_downstream_cor = df4heatmap_downstream_downstream_cor,
    df4heatmap_sc_coexpression = df4heatmap_sc_coexpression,
    gene_anno = gene_anno
  )
  saveRDS(data4figure, paste0(outdir, "/", date, "_data4double_heatmap.RDS"))
  
  pdf(paste0(plotsdir, "/", date, '_downstream_sc_cor_heatmap.pdf'), width =25, height = 20)
  print(double_heatmap(cor_mat1 = df4heatmap_downstream_downstream_cor, 
                       cor_mat2 = df4heatmap_sc_coexpression,
                       min_c = -0.7, max_c = 0.7,
                       treeheight_row = 0, treeheight_col = 0,
                       annotation_row = gene_anno))
  dev.off()
  
}

## ---- ER Stress


params_list <- list(
  date = date,
  home_folder = HomeFolder,
  section_name = section_name,
  complex_of_interest = "nmd"
)

if (!(dir.exists(paste0(outdir, '/', tolower(gsub( x = params$complex_of_interest, pattern = " ", replacement = '_')))))){
  dir.create(paste0(outdir, '/', tolower(gsub( x = params$complex_of_interest, pattern = " ", replacement = '_'))))
}
rmarkdown::render(file.path(CodeFolder, ExperimentName, "pipeline", paste0(section_name, "_", params_list$complex_of_interest, ".Rmd")),
                  output_file = file.path(plotsdir, paste0(date, "_highly_correlated_downstream_genes_", params_list$complex_of_interest, ".html")),
                  params = params_list)

## ---- ER Stress


params_list <- list(
  date = date,
  home_folder = HomeFolder,
  section_name = section_name,
  complex_of_interest = "histones_and_apoptosis"
)

if (!(dir.exists(paste0(outdir, '/', tolower(gsub( x = params$complex_of_interest, pattern = " ", replacement = '_')))))){
  dir.create(paste0(outdir, '/', tolower(gsub( x = params$complex_of_interest, pattern = " ", replacement = '_'))))
}
rmarkdown::render(file.path(CodeFolder, ExperimentName, "pipeline", paste0(section_name, "_", params_list$complex_of_interest, ".Rmd")),
                  output_file = file.path(plotsdir, paste0(date, "_highly_correlated_downstream_genes_", params_list$complex_of_interest, ".html")),
                  params = params_list)




## ---- Hypoxia Genes


params_list <- list(
  date = date,
  home_folder = HomeFolder,
  section_name = section_name,
  complex_of_interest = "hypoxia_pathway"
)
rmarkdown::render(file.path(CodeFolder, ExperimentName, "pipeline", paste0(section_name, "_", params_list$complex_of_interest, ".Rmd")),
                  output_file = file.path(plotsdir, paste0(date, "_highly_correlated_downstream_genes_", params_list$complex_of_interest, ".html")),
                  params = params_list)



## ---- Apoptosis Genes


params_list <- list(
  date = date,
  home_folder = HomeFolder,
  section_name = section_name,
  complex_of_interest = "apoptosis_pathway"
)
rmarkdown::render(file.path(CodeFolder, ExperimentName, "pipeline", paste0(section_name, "_", params_list$complex_of_interest, ".Rmd")),
                  output_file = file.path(plotsdir, paste0(date, "_highly_correlated_downstream_genes_", params_list$complex_of_interest, ".html")),
                  params = params_list)



## ---- Cholesterol biosynthesis


params_list <- list(
  date = date,
  home_folder = HomeFolder,
  section_name = section_name,
  complex_of_interest = "cholesterol_biosynthesis_pathway"
)

if (!(dir.exists( paste0(outdir, "/", tolower(gsub( x = params$complex_of_interest, pattern = " ", replacement = '_')))))){
  dir.create( paste0(outdir, "/", tolower(gsub( x = params$complex_of_interest, pattern = " ", replacement = '_'))))
}
rmarkdown::render(file.path(CodeFolder, ExperimentName, "pipeline", paste0(section_name, "_", params_list$complex_of_interest, ".Rmd")),
                  output_file = file.path(plotsdir, paste0(date, "_highly_correlated_downstream_genes_", params_list$complex_of_interest, ".html")),
                  params = params_list)


