#xvfb-run -a /software/R-4.1.3/bin/R
date <- "2022-08-15"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
section_name <- "5b_expression_UMAPs"
subsection_name <- "5b_03_expression_UMAPs_technical_clusters"
correction <- "technical"
perturbation_status <- 'all'
ExperimentName <- "Magpie"

library(ggpubr)
library(data.table)
library(Seurat)
library(dplyr)
library(reshape2)
library(ggplot2)
library(doParallel)
library(ggrastr)
library(ggpointdensity)
library(tidyverse)
library(viridis)

# options
args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
  if (arg == "--correction"){ correction <- args[[ind + 1]] }
  if (arg == "--experiment_name"){ ExperimentName <- args[[ind + 1]] }
  if (arg == "--perturbation_status"){ perturbation_status <- args[[ind + 1]] }
}



# i/o
print(paste0("Experiment name: ", ExperimentName))
print(paste0("Perturbation status: ", perturbation_status))
print(paste0("Correction: ", correction))
og_experiment_name <- ExperimentName


ExperimentName <- og_experiment_name
if (ExperimentName == 'both'){
  ExperimentName <- "Pica"
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

config <-paste0( og_experiment_name, "_combined_", og_perturbation_status, "_cells_correction_", correction)

# i/o
print(paste0("date: ", date))
source(file.path(HomeFolder, "scripts/io/", paste0(ifelse(og_experiment_name == "Magpie", "Magpie", "Pica"), "_io.R")))
outdir <- file.path(OutFolder, section_name)
plotsdir <- file.path(HTMLFolder, "pipeline", section_name, config)
if(!dir.exists(outdir)) dir.create(outdir)
if(!dir.exists(plotsdir)) dir.create(plotsdir)

## ---- Load Data

inlets_combined <- readRDS(paste0(outdir, '/', date, '_', config, '.RDS'))
crispr_counts <- fread(file.path(OutFolder, "5a_combine", paste0(date, "_crispr_sequence_counts.tsv")))
line_metadata <- fread(paste0(FilesFolder, "LineMetadata.tsv"))
line_metadata <- line_metadata[!(duplicated(line_metadata$name)),] %>%
  dplyr::select(c('name', 'sex', 'disease_status', 'age', 'pluripotency_score', 'novelty_score'))

fnm <- sort(list.files(paste0(ResourcesFolder, "/ipsc_marker_genes_from_sunay/"),
                       pattern = "Jerber_iNeuron_Differentiation_Efficiency_by_Line", full.names = T), decreasing = T)[1]
jerber_differentiation_efficiency <- fread(fnm)
jerber_differentiation_efficiency$cell_line <- unlist(lapply(strsplit(jerber_differentiation_efficiency$donor_id, "-"), "[[", 2))
jerber_differentiation_efficiency$donor <- unlist(lapply(strsplit(unlist(lapply(strsplit(jerber_differentiation_efficiency$donor_id, "-"), "[[", 2)), "_"), "[[", 1))
jerber_differentiation_efficiency <- jerber_differentiation_efficiency %>%
  dplyr::filter(cell_line %in% line_metadata$name)

fnm <- sort(list.files(paste0(ResourcesFolder, "/ipsc_marker_genes_from_sunay/"),
                       pattern = "Jerber_iNeuron_Differentiation_Efficiency_downloaded", full.names = T), decreasing = T)[1]
jerber_de_eff_deg <- fread(fnm)

## ---- DEG
# should be done later, but to save time

if (F){
  high_novelty_clusters <- c("1", "2", "12")
  low_novelty_clusters <- c("0", "3", "4", "5", "6")
  
  
  cell_meta <- inlets_combined@meta.data %>%
    mutate(Status = ifelse(as.character(seurat_clusters) %in% low_novelty_clusters, 'low_novelty', ifelse(as.character(seurat_clusters) %in% high_novelty_clusters, 'high_novelty', 'unassigned'))) %>%
    dplyr::filter(Status != 'unassigned') %>%
    dplyr::select(c("Status", 'orig.ident', 'nCount_RNA', 'percent_MT', 'S.Score', 'G2M.Score'))
  
  downstream_genes2test <- rownames(inlets_combined[["RNA"]]@data)
  registerDoParallel(min(length(downstream_genes2test), 50))
  regressed_lfc <- foreach(downstream_gene2test = downstream_genes2test) %dopar% {
    df4lm <- cell_meta %>%
      mutate(y = as.vector(inlets_combined[["RNA"]]@data[downstream_gene2test, row.names(cell_meta)]))
    fit <- lm(y ~ Status + orig.ident + nCount_RNA + percent_MT + S.Score + G2M.Score, data = df4lm)
    rtn <- c(lfc = summary(fit)$coefficients['Statuslow_novelty', 1],
             pval = summary(fit)$coefficients['Statuslow_novelty', 4])
    return(rtn)
  }
  regressed_lfc <- as.data.frame(bind_rows(regressed_lfc)) %>% mutate(
    downstream_gene = downstream_genes2test,
    downstream_gene_name = unlist(lapply(strsplit(downstream_genes2test, ":"), "[[", 2))
  )
  fwrite(regressed_lfc, paste0(outdir, '/', date, '_', config, '_cluster_deg.tsv'), sep = '\t')
}


## ---- General UMAPs
# take a look at what we have before going into more detail

df4umap <- data.frame(
  id = row.names(inlets_combined@meta.data),
  UMAP_1 = Embeddings(inlets_combined, reduction = 'umap')[,1],
  UMAP_2 = Embeddings(inlets_combined, reduction = 'umap')[,2]
)
df4umap <- left_join(df4umap,
                     inlets_combined@meta.data %>% rownames_to_column('id')) %>%
  mutate(id = gsub(id, pattern = ".*[.]RDS_", replacement = ''))
df4umap <- left_join(df4umap, line_metadata, by = c('cell_line' = 'name'))
df4umap <- left_join(df4umap, crispr_counts %>%
                       mutate(id = gsub(id, pattern = "Moderate_|Strong_", replacement = "")))
fwrite(df4umap, paste0(outdir, '/', config, '_umap_coords_with_line_meta.tsv.gz'), sep = '\t', compress = 'gzip')

set.seed(0)
df4umap <- fread(paste0(outdir, '/', config, '_umap_coords_with_line_meta.tsv.gz'))
df4plot <- df4umap[sample(1:dim(df4umap)[1], 100000, replace = F),]


# run Rmd and create htmls
rmarkdown::render(file.path(CodeFolder, ExperimentName, "pipeline", paste0(subsection_name, ".Rmd")),
                  envir = new.env(),
                  output_file = paste0(plotsdir, "/", date, "_", config, ".html"),
                  params =  list(correction = correction,
                                 experiment_name = og_experiment_name,
                                 perturbation_status = og_perturbation_status,
                                 section_name = section_name,
                                 subsection_name = subsection_name))

sessionInfo()