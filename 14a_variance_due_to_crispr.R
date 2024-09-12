suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(doParallel))
suppressMessages(library(umap))

# set relevent i/o paths
section_name <- "14a_variance_due_to_crispr"


# set relevent i/o paths
date <- "2022-10-05"
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

targets_of_interest <- pica_target_meta %>%
  dplyr::filter(n_downstream_excl_target > 100 & control_norm_expr > 0.5) %>%
  dplyr::select(c('gene', 'n_downstream_excl_target', 'control_norm_expr' ))

## ---- Coordinates for UMAP

## get list of deg to care about
fnm <- sort(list.files(file.path(OutFolder, '6b_calc_lfcs_transcriptome_wide_by_gene'), pattern = '_combined_with_adj_pval.tsv.gz', full.names = T), decreasing = T )[1]
mean_effect <- fread(fnm)
deg <- mean_effect %>%
  dplyr::filter(target %in% targets_of_interest$gene & pval_adj < sig_pval_thresh & 
                  !(downstream_gene_name %in% targets_of_interest$gene)) %>%
  .$downstream_gene %>% table() %>% as.data.frame() %>%
  slice_max(Freq, n = 100) %>% .$`.` %>% as.character()

control_expression <- get_control_expression()
control_deg_expression <- control_expression[, deg]

knockdown_expression_all <- mapply(tg = targets_of_interest$gene, FUN = function(tg){
  rtn <- get_knockdown_expression(target = tg)
  # get some special gene
  dCas9_counts <- as.vector(rtn %>% select(contains('dCas9-KRAB-MeCP2')))[,1]
  # on target
  on_target_expr <- as.vector(rtn %>% select(contains(paste0(":", tg, ":"))))[,1]
  cell_meta <- rtn %>% dplyr::select(c("target", 'cell_line', 'donor', 'id'))
  rtn <- rtn[,deg ]
  rtn$on_target_expr <- on_target_expr
  rtn$dCas9_counts <- dCas9_counts
  rtn <- as.data.frame(bind_cols(cell_meta, rtn))
  return(rtn)
}, SIMPLIFY = F) %>% bind_rows() %>% as.data.frame()


df4umap <- bind_rows(knockdown_expression_all[, deg],
                     control_deg_expression[, deg])
## calculate umap coordinates
set.seed(0)
umap_coords <- umap(as.matrix(df4umap))
umap_coords <- data.frame(
  UMAP_1 = umap_coords$layout[,1],
  UMAP_2 = umap_coords$layout[,2]
)

## ---- CellMetadata
## 
most_recent <- sort(gsub(list.files(file.path(OutFolder, "5a_combine"), pattern = "perturbed_cells_metadata.tsv|control_cells_metadata.tsv"), pattern = "_.*", replacement = ""), decreasing = T)[1]
pica_cell_metadata_fnms <- grep(list.files(file.path(OutFolder, "5a_combine"), pattern = "perturbed_cells_metadata.tsv|control_cells_metadata.tsv"),
                                pattern = most_recent, value = T)
pica_cell_metadata <- mapply(pica_cell_metadata_fnms, FUN = function(fnm){
  fread(file.path(OutFolder, "5a_combine", fnm))
}, SIMPLIFY = F)
pica_cell_metadata <- as.data.frame(bind_rows(pica_cell_metadata))

## guide freq df
most_recent <- sort(gsub(list.files(file.path(OutFolder, "4a_add_guides"), pattern = "_guide_freqs.txt"), pattern = "_.*", replacement = ""), decreasing = T)[1]
guide_freq_fnms <- grep(list.files(file.path(OutFolder, "4a_add_guides"), pattern = "_guide_freqs.txt"), pattern = most_recent, value = T)
guide_freqs <- lapply(guide_freq_fnms, FUN = function(fnm){
  rtn <- fread(file.path(OutFolder, "4a_add_guides", fnm)) %>%
    mutate(inlet = gsub(fnm, pattern = paste0(most_recent, "_|_guide_freqs.txt"), replacement = ""))
  return(rtn)
}) %>% bind_rows() %>% as.data.frame()
guide_freqs$V1 <- paste0(guide_freqs$inlet, "_", guide_freqs$cell_barcode)

cell_meta <- as.data.frame(bind_rows(
  knockdown_expression_all %>% dplyr::select(c('target', "cell_line", "donor", 'id', 'on_target_expr', 'dCas9_counts')),
  control_expression %>% 
    dplyr::select(c( "cell_line", "donor", 'id')) %>%
    mutate(on_target_expr = 0,
           target = NonTargetGeneName,
           dCas9_counts = (control_expression  %>% select(contains('dCas9-KRAB-MeCP2')))[,1]
)))

### add in guide umis
cell_meta <- left_join(cell_meta, guide_freqs %>%
                         dplyr::select(c('id' = "V1", 'top_guide_freq' = "top_guide", 'top_guide_nUMI' = 'n_umi_top_guide')))


## add some 'effect size' metrics
cell_meta$abs_lfc_deg_mean <- c(apply(abs(knockdown_expression_all[, deg]), 1, mean), apply(abs(control_deg_expression[, deg]), 1, mean))
cell_meta$abs_lfc_deg_var <- c(apply(abs(knockdown_expression_all[, deg]), 1, var), apply(abs(control_deg_expression[, deg]), 1, var))


## ---- Save

df2save <- as.data.frame(bind_cols(cell_meta, umap_coords))
fwrite(df2save, paste0(outdir, '/', date, '_umap_coords.tsv'), sep = '\t')

## ---- Effect Size for all

expressed_targets <- pica_target_meta %>%
  dplyr::filter(!(is.na(control_norm_expr)) & n_cells_total > 100)
knockdown_expression_all <- mapply(tg = targets_of_interest$gene, FUN = function(tg){
  rtn <- get_knockdown_expression(target = tg)
  on_target_expr <- as.vector(rtn %>% select(contains(paste0(":", tg, ":"))))[,1]
  cell_meta <- rtn %>% dplyr::select(c("target", 'cell_line', 'donor', 'id'))
  rtn <- rtn[,deg ]
  rtn$on_target_expr <- on_target_expr
  rtn$dCas9_counts <- dCas9_counts
  rtn <- as.data.frame(bind_cols(cell_meta, rtn))
  return(rtn)
}, SIMPLIFY = F) %>% bind_rows() %>% as.data.frame()



## ---- MakePlots

rmarkdown::render(input = paste0(CodeFolder, 'Magpie/pipeline/', section_name, '.Rmd'),
                  output_file = paste0(plotsdir, '/', date, '.html'))
