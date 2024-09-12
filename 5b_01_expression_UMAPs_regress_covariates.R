#xvfb-run -a /software/R-4.1.3/bin/R
date <- "2022-10-05"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
section_name <- "5b_expression_UMAPs"
subsection_name <- "5b_01_UMAPs_regress_covariates"
REDO <- F
ExperimentName <- "Pica" # options are: "Magpie", "Pica" or "both"
## assigned - do for "Magpie", "Pica" and "both"
## control/unassigned - do for "Pica" and "both"
correction <- "none" # optionas are "none", "technical" and "all"
perturbation_status <- "control" # options are: "perturbed", "control", "unassigned", "assigned" (control + perturbed), "unperturbed" (control + unassigned), "all"
magpie_date <- "2022-08-15"; pica_date <- '2022-10-05'

library(ggpubr)
library(data.table)
library(Seurat)
library(dplyr)
library(reshape2)
#library(ggplot2)
#library(doParallel)
#library(ggrastr)
#library(ggpointdensity)
library(tidyverse)
#library(viridis)

# options
args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--correction"){ correction <- args[[ind + 1]] }
  if (arg == "--experiment_name"){ ExperimentName <- args[[ind + 1]] }
  if (arg == "--perturbation_status"){ perturbation_status <- args[[ind + 1]] }
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
}

print(paste0("Experiment name: ", ExperimentName))
print(paste0("Perturbation status: ", perturbation_status))
print(paste0("Correction: ", correction))
og_experiment_name <- ExperimentName

# i/o
print(paste0("date: ", date))
source(file.path(HomeFolder, "scripts/io/", paste0(ifelse(og_experiment_name == "Magpie", "Magpie", "Pica"), "_io.R")))
outdir <- file.path(OutFolder, section_name)
plotsdir <- file.path(HTMLFolder, "pipeline/", section_name)
if(!dir.exists(outdir)) dir.create(outdir)
if(!dir.exists(plotsdir)) dir.create(plotsdir)

## ---- LoadData
# figure out which files you want to load

ExperimentName <- og_experiment_name
if (ExperimentName == 'both'){
  ExperimentName <- c("Magpie", "Pica")
}

og_perturbation_status <- perturbation_status
if (perturbation_status == 'control'){
  perturbation_status = c('nontarget', 'control')
} else if (perturbation_status == 'assigned'){
  perturbation_status = c("perturbed","nontarget", "control")
} else if (perturbation_status == 'unperturbed'){
  perturbation_status = c("unassigned", "control", "nontarget")
} else if (perturbation_status == 'all'){
  perturbation_status = c("unassigned", "control", "nontarget", "perturbed")
}


all_possible_seurat_objects <- c(grep(list.files(paste0(ProjectFolder, 'outs/Magpie/5a_combine'), pattern = '[.]RDS', full.names = T), pattern = magpie_date, value = T),
                                 grep(list.files(paste0(ProjectFolder, 'outs/Pica/5a_combine'), pattern = '[.]RDS', full.names = T), pattern = pica_date, value = T))
all_possible_seurat_objects <- grep(all_possible_seurat_objects, pattern = "_ALL", invert = T, value = T)

seurat_objects_to_load <- c()
for (e in ExperimentName){
  for (p in perturbation_status){
    fnms2load <- grep(grep(all_possible_seurat_objects, pattern = paste0(p, "_cells"), value = T),
                      pattern = e, value = T)
    seurat_objects_to_load <- c(seurat_objects_to_load, fnms2load)
  }
}
seurat_objects_to_load <- unique(seurat_objects_to_load)

inlets_combined <- mapply(seurat_objects_to_load, FUN = readRDS, SIMPLIFY = F)
expressed_in_both <- table(unlist(lapply(inlets_combined, FUN = function(rtn){
  row.names(rtn[["RNA"]]@data)
})))
expressed_in_both <- names(expressed_in_both)[which(expressed_in_both== max(expressed_in_both))]
downstream_genes_of_interest <- sort(unique(unlist(lapply(inlets_combined, FUN = function(rtn){
  row.names(rtn[["RNA"]]@scale.data)
}))))
downstream_genes_of_interest <- intersect(downstream_genes_of_interest, expressed_in_both)
inlets_combined <- lapply(inlets_combined, FUN = function(rtn){
  counts <- rtn[["RNA"]]@counts[downstream_genes_of_interest,]
  new_object <- CreateSeuratObject(counts, meta.data = rtn@meta.data)
  return(new_object)
})
#inlets_combined = inlets_combined[-which(sapply(inlets_combined, is.null))]
#
## focus only on the highly variable genes since it is just a UMAP

inlets_combined <- merge(x = inlets_combined[[1]],
                         y = inlets_combined[-1],
                         add.cell.ids = names(inlets_combined),
                         project = paste(og_experiment_name, og_perturbation_status, sep = "_"),
                         merge.data = TRUE)
inlets_combined <- JoinLayers(inlets_combined)
inlets_combined <- NormalizeData(inlets_combined)
print('combined object:')
inlets_combined



## ---- Regression

if (correction == "technical"){
  
  inlets_combined <- ScaleData(inlets_combined,
                               vars.to.regress = c("S.Score", "G2M.Score", "orig.ident",
                                                   "nCount_RNA", "percent_MT"),
                               features = rownames(inlets_combined))
  #table(inlets_combined@meta.data$seurat_clusters)
  
}

if (correction == "all"){
  inlets_combined <- ScaleData(inlets_combined,
                               vars.to.regress = c('cell_line', "S.Score", "G2M.Score", "orig.ident",
                                                   "nCount_RNA", "percent_MT"),
                               features = rownames(inlets_combined))
  
}

if (correction == "none"){
  inlets_combined <- ScaleData(inlets_combined)
}

inlets_combined <- FindVariableFeatures(inlets_combined)
inlets_combined <- RunPCA(inlets_combined, npcs = min(50, ncol(inlets_combined)/2))
inlets_combined <- RunUMAP(inlets_combined, dims = 1:20)

inlets_combined <- FindNeighbors(inlets_combined, dims = 1:10)
## for technical correction, there should be two large clusters
inlets_combined <- FindClusters(inlets_combined, resolution = 0.5)

#DimPlot(inlets_combined, group.by = 'donor')

## ---- SaveRDS

fnm_out <- paste0(outdir, "/", date, "_", og_experiment_name, "_combined_", og_perturbation_status, "_cells_correction_", correction, ".RDS")

print('writing to')
print(fnm_out)
saveRDS(inlets_combined, fnm_out )


## ----

sessionInfo()