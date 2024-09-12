library(Matrix)
library(Seurat)
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(R.utils)
library(ggpubr)
library(SeuratObject)
library(doParallel)

HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
section_name <- "12e_matrix_for_DCDI"
date <- "2022-08-15"
ExperimentName <- "Magpie"
#gene_set <- "hypoxia_pathway"
gene_set <- "apoptosis_genes"
scale_factor <- 1
lfc_base <- 10


args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  ##needed arguments
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--gene_set"){ gene_set <- args[[ind + 1]] }
  if (arg == "--experiment_name"){ ExperimentName <- args[[ind + 1]] }
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
}


# set relevent i/o paths
source(file.path(HomeFolder, "scripts/io/", paste0(ExperimentName, "_io.R")))
source(paste0(HomeFolder, "scripts/io/Magpie_Utils.R"))
outdir <- file.path(OutFolder, section_name, gene_set, '/')
plotsdir <- file.path(HTMLFolder, "pipeline/", section_name)
if(!dir.exists(outdir)) dir.create(outdir, recursive = T)
if(!dir.exists(file.path(outdir, 'control_lm'))) dir.create(file.path(outdir, 'control_lm'), recursive = T)
if(!dir.exists(file.path(outdir, 'lfcs'))) dir.create(file.path(outdir, 'lfcs'), recursive = T)
if(!dir.exists(plotsdir)) dir.create(plotsdir)

LaneMetadata <- fread(LaneMetadataPath)
GuideMetadata <- fread(GuideMetadataPath)
IDs <- LaneMetadata$ID
control_tag <- c("unassigned", NonTargetGeneName)

## genes to consider
co_perturbed_genes <- fread(paste0(OutFolder, "12c_highly_correlated_downstream_genes", '/', date, "_", gene_set, '.tsv'))


# read all data
files <- list.files(paste0(OutFolder, "4b_add_donors/"))
files <- files[grepl("_with_donor.RDS", files) & !grepl("combined", files)]
files <- files[gsub(files, pattern = paste0(date, "_|_with_donor.RDS"), replacement = "") %in% IDs]
#use only IDs where data already exists
IDs_available <- gsub(files, pattern = paste0(date, "_|_with_donor.RDS"), replacement = "")

print("IDs without data:")
print(IDs[which(!(IDs %in% IDs_available))])
print("IDs going into combined object:")
IDs <- IDs[which(IDs %in% IDs_available)]
print(paste(length(IDs), "inlets total"))

#load data
print("Loading data...")
inlet_list <- lapply(files, function(fnm){
    single_inlet <- readRDS(file.path(OutFolder, "4b_add_donors", fnm))
    ID <- gsub(fnm, pattern = paste0(date, "_|_with_donor.RDS"), replacement = "")
  
  gene_inds <- which(unlist(lapply(strsplit(row.names(single_inlet[["RNA"]]@data), ":"), "[[", 2)) %in% co_perturbed_genes$gene)
  
  single_inlet[["RNA"]]@data <- single_inlet[["RNA"]]@data[gene_inds, ]
  single_inlet[["RNA"]]@counts <- single_inlet[["RNA"]]@counts[gene_inds, ]
  single_inlet[["RNA"]]
  single_inlet <- subset(single_inlet, target_gene_name %in% c(co_perturbed_genes$gene, NonTargetGeneName, "unassigned"))
  return(single_inlet)
})
names(inlet_list) <- IDs_available


## combine all inlets (strong, moderate or all)
relevant_IDs <- LaneMetadata$ID
relevant_IDs <- intersect(relevant_IDs, IDs_available)
relevant_inlets <- inlet_list[relevant_IDs]

inlets_combined <- merge(x = relevant_inlets[[1]],
                         y = relevant_inlets[-1],
                         add.cell.ids = names(relevant_inlets),
                         project = paste(ExperimentName, "combined", sep = "_"),
                         merge.data = TRUE)
nms <- rownames(inlets_combined@meta.data)
inlets_combined@meta.data <- inlets_combined@meta.data %>% select(c("orig.ident", "nCount_RNA", "nFeature_RNA", "nCount_guides", "nFeature_guides", "percent_MT", "num_umis",
                                                                    "feature_call", "guide_call", "target_gene_name", "num_features",
                                                                    "donor", "cell_line", "hipsci_line_name",
                                                                    "S.Score", "G2M.Score", "Phase"
))
inlets_combined@meta.data <- left_join(inlets_combined@meta.data,
                                       select(LaneMetadata, ID, Guide_Group, Guide_Strength, Pool, Timepoint, Batch),
                                       by = c("orig.ident" = "ID"))
rownames(inlets_combined@meta.data) <- nms

print("seurat object combined successfully!")
print(inlets_combined)
saveRDS(inlets_combined, file = paste0(outdir,'/', date, "_", gene_set, "_combined.RDS"))
write.csv(inlets_combined@meta.data,
          file =paste0(outdir, '/', date, "_", gene_set, "_combined_meta.csv"))
print("seurat object written")

## ---- control lm

control_cells <- subset(inlets_combined, target_gene_name %in% c(NonTargetGeneName, "unassigned"))
normalized_expression <- control_cells[["RNA"]]@data
cell_meta4lm <- data.frame(cell_line = control_cells$cell_line,
                      inlet = c(control_cells$orig.ident),
                      nCount_RNA = c(control_cells$nCount_RNA),
                      percent_MT = c(control_cells$percent_MT),
                      s.score = c(control_cells$S.Score),
                      g2m.score = c(control_cells$G2M.Score)
)
registerDoParallel(10)
foreach (downstream_gene_ind = 1:dim(normalized_expression)[1]) %dopar% {
  gene_name <- row.names(normalized_expression)[downstream_gene_ind]
  
  df4lm <- cell_meta4lm %>%
    mutate(y4lm = normalized_expression[downstream_gene_ind,])
  fit <- lm(y4lm ~ cell_line + inlet + nCount_RNA + percent_MT + s.score + g2m.score, data = df4lm)
  
  #some betas are NA [change so that they're the mean effect across other values with the same variable]
  rtn <- list(residuals_var = var(fit$residuals), coefficients = fit$coefficients)
  rtn$coefficients[which(is.na(rtn$coefficients))] <- 0
  
  n_control = length(inlets_combined$orig.ident) + length(inlets_combined$orig.ident)
  rtn$n_control <- n_control
  #control_norm_expr
  rtn$control_num_reads <- sum(inlets_combined[["RNA"]]@counts[downstream_gene_ind,])
  rtn$control_norm_expr <- sum(inlets_combined[["RNA"]]@data[downstream_gene_ind,])/n_control
  
  rtn$av_log_expression_control <- log(sum(expm1(as.matrix(inlets_combined[["RNA"]]@data[downstream_gene_ind, ]))/n_control) + 1/scale_factor, base = lfc_base) 
  
  saveRDS(rtn, file = file.path(outdir, '/control_lm/', paste0(date, "_gene-", gsub(gene_name, pattern = "/", replacement = "-"), "_coefs.RDS")))
  
  ## save the residuals
  saveRDS(residuals(fit), file = paste0(outdir,'/control_lm/', date, "_gene-", gsub(gene_name, pattern = "/", replacement = "-"), "_residuals.RDS"))
}
                     
## ---- LM

perturbed_cells <- subset(inlets_combined, target_gene_name %in% co_perturbed_genes$gene)
targets <- unique(perturbed_cells$target_gene_name)
targets <- targets[which(targets %in% unlist(lapply(strsplit(row.names(perturbed_cells), ":"), "[[", 2)))]

registerDoParallel(10)
foreach (target_ind = 1:length(targets)) %dopar% {
  cells_by_target <- subset(perturbed_cells, target_gene_name == targets[target_ind])
  normalized_expression <- cells_by_target[["RNA"]]@data
  perturbed_cell_metadata <- data.frame(cell_line = cells_by_target$cell_line,
                                        inlet = c(cells_by_target$orig.ident),
                                        nCount_RNA = c(cells_by_target$nCount_RNA),
                                        percent_MT = c(cells_by_target$percent_MT),
                                        s.score = c(cells_by_target$S.Score),
                                        g2m.score = c(cells_by_target$G2M.Score)
  )
  
  lm_out <- lapply(1:dim(cells_by_target[["RNA"]]@data)[1], FUN = function(downstream_gene_ind){
    
    downstream_gene <- row.names(normalized_expression)[downstream_gene_ind]
    control_lm <- readRDS(file.path(outdir, 'control_lm', paste0(date, "_gene-", gsub(downstream_gene, pattern = "/", replacement = "-"), "_coefs.RDS")))
    
    reformatted_perturbed_cell_metadata <- reformat_metadata_matrix(in_matrix = perturbed_cell_metadata, coefficients = control_lm$coefficients)
    lm_res <- run_lm(perturbed_normalized_expression = normalized_expression[downstream_gene_ind,], 
                     perturbed_cell_metadata = reformatted_perturbed_cell_metadata,
                     control_fit = control_lm)
    lm_res$n_control = control_lm$n_control
    lm_res$control_num_reads = control_lm$control_num_reads
    lm_res$control_norm_expr = control_lm$control_norm_expr
    lm_res$av_log_expression_control = control_lm$av_log_expression_control
    
    return(lm_res)
  })
  names(lm_out) <- row.names(normalized_expression)
  lm_fit <- as.data.frame(t(bind_rows(lapply(lm_out, 
                                             FUN = function(x){
                                               return(c(x$lfc, x$v, x$pval, x$n_control, x$control_num_reads, x$control_norm_expr, x$av_log_expression_control))
                                             }))))
  names(lm_fit) <- c("lfc", "v", "pval", "n_control", "control_num_reads", "control_norm_expr", "av_log_expression_control")
  
  lm_fit$perturbed_norm_expr <- rowMeans(normalized_expression)
  lm_fit$av_log_expression_perturbed <- log(rowMeans(expm1(normalized_expression)) + 1/scale_factor, base = lfc_base)
  lm_fit$lfc_unadjusted <- lm_fit$av_log_expression_perturbed - lm_fit$av_log_expression_control
  lm_fit$perturbed_num_reads <- rowSums(cells_by_target@assays$RNA@counts)
  rtn <- data.frame(downstream_gene = row.names(normalized_expression),
                    n_control = lm_fit$n_control,
                    n_perturbed = dim(normalized_expression)[2],
                    lfc = lm_fit$lfc,
                    lfc_unadjusted = lm_fit$lfc_unadjusted,
                    v = lm_fit$v,
                    pval_lm = lm_fit$pval,
                    control_norm_expr = lm_fit$control_norm_expr,
                    perturbed_norm_expr = lm_fit$perturbed_norm_expr,
                    control_num_reads = lm_fit$control_num_reads,
                    perturbed_num_reads = lm_fit$perturbed_num_reads,
                    target = targets[target_ind])
  
  
  saveRDS(rtn, file = file.path(outdir, 'lfcs', paste0(date, "_", targets[target_ind], "_coefs.RDS")))
  
  ## save the residuals
  
  #print residuals
  cnms <- c("id", names(lm_out))
  lm_res <- data.frame(
    id = row.names(perturbed_cell_metadata),
    as.data.frame(bind_cols(mapply(lm_out, FUN = function(x){return(x$res)})))
  )
  colnames(lm_res) <- cnms
  
  saveRDS(lm_res, file = paste0(outdir,'/lfcs/', date, "_", targets[target_ind], "_residuals.RDS"))
  
}

## ---- Combine

## residuals
fnms <- list.files(file.path(outdir, 'lfcs'), pattern = 'residuals')
resid <- lapply(fnms, FUN = function(f){
  rtn <- readRDS(file.path(outdir, 'lfcs', f))
}) %>% bind_rows() %>% as.data.frame()
fwrite(resid, paste0(outdir, '/', date, "_resid.tsv.gz"), sep = '\t', compress = 'gzip')


## lfcs
fnms <- list.files(file.path(outdir, 'lfcs'), pattern = 'coefs')
lfcs <- lapply(fnms, FUN = function(f){
  rtn <- readRDS(file.path(outdir, 'lfcs', f))
}) %>% bind_rows() %>% as.data.frame()
lfcs$downstream_gene_name <- unlist(lapply(strsplit(lfcs$downstream_gene, ":"), "[[", 2))
lfcs$pval_adj <- p.adjust(lfcs$pval_lm, method = "BH")
fwrite(lfcs, paste0(outdir, '/', date, "_lfcs.tsv.gz"), sep = '\t', compress = 'gzip')



## ---- SessionInfo

sessionInfo()
      