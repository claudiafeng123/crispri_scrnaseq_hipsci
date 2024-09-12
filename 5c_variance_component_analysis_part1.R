library(tidyverse)
library(lme4)
library(reshape2)
library(Seurat)
library(doParallel)
library(data.table)
library(variancePartition)

HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
section_name <- "5c_variance_component_analysis"
date <- "2022-08-15"
date_out <- "2024-05-23"

#settings:
n_genes_total <- 10
start_ind <- 0
redo <- F

# options
args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--date_out"){ date_out <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
  if (arg == "--n_genes"){ n_genes_total <- as.numeric(args[[ind + 1]]) }  
  if (arg == "--start_ind"){ start_ind <- as.numeric(args[[ind + 1]]) }  
  if (arg == "--redo"){ redo <- as.logical(args[[ind + 1]])}
}

# i/o
source(file.path(HomeFolder, "scripts/io/Magpie_io.R"))
outdir <- file.path(OutFolder, section_name, date_out, "by_ind")
plotsdir <- file.path(HTMLFolder, "pipeline/", section_name)
if(!dir.exists(outdir)) dir.create(outdir, recursive = T)
if(!dir.exists(plotsdir)) dir.create(plotsdir)
date_new <- date
inlet_strengths <- c("Moderate", "Strong")


genes2test <- fread(file.path(OutFolder, "5a_combine", paste0(date, "_gene_list.txt")))$claudia
genes2test <- genes2test[(start_ind + 1):(min(length(genes2test), start_ind + n_genes_total))]


## load data

inlets_combined <- mapply(inlet_strength = paste0(rep(inlet_strengths, rep(2, length(inlet_strengths))), rep(c("_perturbed", "_nontarget", "_unassigned"), length(inlet_strengths))), FUN = function(inlet_strength){
  readRDS(file.path(OutFolder, "5a_combine",
                    paste0(date, "_", inlet_strength, "_cells_combined.RDS")))
}, SIMPLIFY = F)

## cell counts
cell_metadata <- lapply(1:length(inlets_combined), FUN = function(i){
  seurat_object <- inlets_combined[[i]]
  rtn <- seurat_object@meta.data %>%
    select(c("orig.ident", "nCount_RNA", "percent_MT", "target_gene_name", "donor", "cell_line", "Phase", "Guide_Strength", "Batch")) %>%
    mutate(status = unlist(lapply(strsplit(names(inlets_combined)[i], "_"), "[[", 2))) %>%
    mutate(status_2 = ifelse(target_gene_name %in% c("unassigned", NonTargetGeneName), "unperturbed", target_gene_name))
  return(rtn)
}) %>% bind_rows() %>% as.data.frame()
cell_metadata$Batch <- as.character(cell_metadata$Batch)
## rescale the continuous variables
cell_metadata$percent_MT <- (cell_metadata$percent_MT - mean(cell_metadata$percent_MT))/sd(cell_metadata$percent_MT)
cell_metadata$nCount_RNA <- (cell_metadata$nCount_RNA - mean(cell_metadata$nCount_RNA))/sd(cell_metadata$nCount_RNA)

expr_df <- lapply(1:length(inlets_combined), FUN = function(i){
  rtn <- data.frame(
    expr = inlets_combined[[i]][["RNA"]]@data[(start_ind + 1):(start_ind + length(genes2test)), ]
  )
  #rtn <- data.frame(
  #  expr = inlets_combined[[i]][["RNA"]]@data[(start_ind + 1):min(length(genes2test), start_ind + n_genes_total)), ]
  #)## 41...10
}) %>% bind_cols() %>% as.data.frame()
expr_df <- as.data.frame(t(expr_df))

rm(inlets_combined); gc()

registerDoParallel(10)
varianceExplained <- foreach (gene_ind = 1:length(genes2test)) %dopar% {
  df4lmm <- as.data.frame(bind_cols(cell_metadata, expr_df[, gene_ind]))
  names(df4lmm)[length(names(df4lmm))] <- "expr"
  fit <- tryCatch(
    lmer(formula = expr ~ (1|donor) + (1|donor:cell_line) + (1|Batch) + (1|Batch:orig.ident) + 
           (1|Phase) + (1|status_2) + percent_MT + nCount_RNA, data = df4lmm),
    error = function(e) {
      return(NULL)
    }) 
  
  if (is.null(fit)){
    varExplained <- NULL
  } else {
    varExplained <- tryCatch(
      calcVarPart(fit, showWarnings = F), 
      error = function(e) {
        return(NULL)
      })
  }
  
  return(varExplained)
}
names(varianceExplained) <- genes2test
varianceExplained <- varianceExplained[!sapply(varianceExplained,is.null)]
gene_names <- names(varianceExplained)

varianceExplained <- as.data.frame(bind_rows(varianceExplained))
varianceExplained$downstream_gene_name <- unlist(lapply(strsplit(gene_names, ":"), "[[", 2))

## write this as a new file
fwrite(varianceExplained, paste0(outdir, '/', date, '_', start_ind, '.tsv'), sep= '\t')

## session info
sessionInfo()