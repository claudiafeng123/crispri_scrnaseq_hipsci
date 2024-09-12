library(tidyverse)
library(data.table)
library(parallel)
library(VennDiagram)

# options
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
section_name="9a_target_gene_downregulation"
min_cells=10 # minimum number of cells for a target to be considered
min_degs=2 # minimum number of DEGs for a target to be considered
ncores=2
ExperimentName <- "Magpie"
date <- "2022-08-15"
#ExperimentName <- "Pica"
#date <- "2022-10-05"


args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--experiment_name"){ ExperimentName <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
  if (arg == "--min_cells"){ min_cells <- as.numeric(args[[ind + 1]])}
  if (arg == "--min_degs"){ min_degs <- as.numeric(args[[ind + 1]])}
  if (arg == "--min_degs"){ ncores <- as.numeric(args[[ind + 1]]) }
}

source(file.path(HomeFolder, "scripts/io/", paste0(ExperimentName, "_io.R")))
source(file.path(HomeFolder, "scripts/io/Magpie_Utils.R"))
outdir <- file.path(OutFolder, section_name)
plotsdir <- file.path(HTMLFolder, "pipeline/", section_name)
if(!dir.exists(outdir)) dir.create(outdir)
if(!dir.exists(plotsdir)) dir.create(plotsdir)

## ---- TargetDownregulation

if (!(exists("min_cells"))){min_cells <- min_cells_per_gene}
if (!(exists("min_degs"))){min_degs <- min_deg_per_target}

LaneMetadata <- fread(LaneMetadataPath)
GuideMetadata <- fread(GuideMetadataPath)

# which targets to include
target_df <- read.csv(file.path(OutFolder, "7b_target_summary", paste0(date, "_target_meta_data.csv")))
target_df <- dplyr::filter(target_df, gene != NonTargetGeneName)
target_df <- mutate(target_df, group = ifelse(group == "Strong", "fitness genes", "other genes"))
names(cols4Pools) <- c("fitness genes", "other genes")
#target_df <- dplyr::filter(target_df, n_cells >= min_cells, n_downstream_excl_target >= min_degs)
#targets_to_include <- target_df$gene

lfcs <- fread(file.path(OutFolder, "6b_calc_lfcs_transcriptome_wide_by_gene", paste0(date, "_combined_with_adj_pval.tsv.gz")))


## ---- CRISPR Marker Genes

## crispr metadata
crispr_sequences <- gsub(list.files(file.path(ResourcesFolder, "reference_sequences", "Magpie_added_sequences"), pattern = ".reference.fa"),
                         pattern = "[.]reference[.]fa", replacement = "")
crispr_counts_df <- fread(paste0(OutFolder, "5a_combine/", date, "_crispr_sequence_counts.tsv"))

## ---- Guide Assignments

if (ExperimentName == "Pica"){
  cell_metadata_fnms <- grep(list.files(file.path(OutFolder, "5a_combine"), pattern = "_cells_metadata.tsv"),
                             pattern = date, value = T)
} else if (ExperimentName == "Magpie"){
  cell_metadata_fnms <- grep(list.files(file.path(OutFolder, "5a_combine"), pattern = "_cells_combined_meta.csv"),
                             pattern = date, value = T)
  
}
cell_metadata <- mapply(cell_metadata_fnms, FUN = function(fnm){
  fread(file.path(OutFolder, "5a_combine", fnm))
}, SIMPLIFY = F)
cell_metadata <- as.data.frame(bind_rows(cell_metadata))
cell_metadata <- select(cell_metadata, c("V1", "orig.ident", "nCount_guides", "nFeature_guides", "num_umis", "target_gene_name", "cell_line"))
if (ExperimentName == "Magpie"){
  cell_metadata <- cell_metadata %>%
    mutate(V1 = ifelse(substr(V1, 4, 5) == "MD", 
                                     paste0("Moderate_", V1), 
                                     paste0("Strong_", V1))
    )
}

## ---- Target Gene Expression

## also have a random assignment for benchmarking
set.seed(0)
guide_assignment_df <- data.frame(
  cell_metadata %>%
  dplyr::select(c("V1", "orig.ident", "target_gene_name", "num_umis")),
  random_assignment = sample(cell_metadata$target_gene_name, replace = F)
) %>% 
  mutate(V1 = gsub(V1, pattern = "Moderate_|Strong_", replacement = ""))

seurat_object_fnms <- grep(grep(list.files(file.path(OutFolder, "3_Load_Data_QC"), pattern = ".RDS"),
                                pattern = date, value = T),
                           pattern = "cutoffs", invert = T, value = T)
on_target_expr <- mapply(seurat_object_fnms, FUN = function(f){
  
  seurat_object <- readRDS(file.path(OutFolder, "3_Load_Data_QC", f))
  
  target_inds <- which(unlist(lapply(strsplit(row.names(seurat_object[["RNA"]]@data), ":"), "[[", 2)) %in% GuideMetadata$gene)
  target_expr <- as.data.frame(t(as.matrix(seurat_object[["RNA"]]@data[target_inds, ])))
  colnames(target_expr) <-  paste0(unlist(lapply(strsplit(colnames(target_expr), ":"), "[[", 2)),  "_expr")
  
  
  ## expression of 
  rtn <- data.frame(
    V1 = paste0(gsub(f, pattern = paste0(date, "_|[.]RDS"), replacement = ""), "_", colnames(seurat_object[["RNA"]]@data)),
    target_expr
  )
  
  return(rtn)
  
}, SIMPLIFY = F)
on_target_expr <- as.data.frame(bind_rows(on_target_expr))
on_target_expr <- left_join(guide_assignment_df, on_target_expr)
on_target_expr$true_assignment_expr <- lapply(1:dim(on_target_expr)[1], FUN = function(i){
  if (paste0(on_target_expr$target_gene_name[i], '_expr') %in% colnames(on_target_expr)){
    rtn <- on_target_expr[i, paste0(on_target_expr$target_gene_name[i], '_expr')]
  } else{
    rtn <- -Inf
  }
  return(rtn)
}) %>% unlist()
on_target_expr$random_assignment_expr <- lapply(1:dim(on_target_expr)[1], FUN = function(i){
  if (paste0(on_target_expr$random_assignment[i], '_expr') %in% colnames(on_target_expr)){
    rtn <- on_target_expr[i, paste0(on_target_expr$random_assignment[i], '_expr')]
  } else{
    rtn <- -Inf
  }
  return(rtn)
}) %>% unlist()

on_target_expr$true_assignment_control_mean_expr <- target_df$control_norm_expr[match(on_target_expr$target_gene_name, target_df$gene)]
on_target_expr$random_assignment_control_mean_expr <- target_df$control_norm_expr[match(on_target_expr$random_assignment, target_df$gene)]
on_target_expr <- dplyr::select(on_target_expr, c("V1", "orig.ident", "num_umis",
                                                  "target_gene_name", "true_assignment_expr", "true_assignment_control_mean_expr",
                                                  "random_assignment", "random_assignment_expr", "random_assignment_control_mean_expr"))
on_target_expr$dCas9_expression <- crispr_counts_df$`dCas9-KRAB-MeCP2`[match(on_target_expr$V1, gsub(crispr_counts_df$id, pattern = "Moderate_|Strong_", replacement = ""))]
on_target_expr <- on_target_expr %>%
  mutate(dCas9_expression = ifelse(dCas9_expression > 0, "Expressed", "Not Expressed"))

fwrite(on_target_expr, paste0(outdir, '/', date, "_on_target_expr_by_cells.tsv"), sep = '\t')

## ---- Bad Guides

## use pica data to check out the guides that we thought were bad??



## ---- Subsample
## For LFC vs. Effect Size

## Choose genes to subsample pica
if (ExperimentName == "Pica"){
  targets2subsample <- filter(target_df, n_cells > 80 & !is.na(target_lfc) & control_norm_expr > 0.5) %>%
    dplyr::select(c("gene", "n_cells", "n_downstream_excl_target", "target_lfc", "control_norm_expr"))
  lots_of_cells <- target_df %>% filter(n_cells > 80) %>% .$gene
  high_expression <- target_df %>% filter(control_norm_expr > 0.5) %>% .$gene
  target_expressed <-target_df %>% filter(!is.na(target_lfc)) %>% .$gene
  venn.diagram(
    x = list(lots_of_cells,
             high_expression,
             target_expressed),
    category.names = c("> 80 Cells", "High Expression of Target Gene", "Target Gene Expressed"),
    filename = paste0(plotsdir, "/", date, "_genes_subsampled.png"), width = 1080, height = 800, imagetype =  'png',
  )
  } else if (ExperimentName == "Magpie"){
  targets2subsample <- filter(target_df, n_cells > 80 & !is.na(target_lfc) & control_norm_expr > 0.1 & n_downstream_excl_target > 20) %>%
    dplyr::select(c("gene", "group", "n_cells", "n_downstream_excl_target", "target_lfc", "control_norm_expr")) %>%
    mutate(group = ifelse(group == "other genes", "Moderate", "Strong"))
}
fwrite(targets2subsample, paste0(outdir, "/", date, "_targets2subsample.tsv"), sep = '\t')


## ---- MakePlots

rmarkdown::render(
  input = file.path(CodeFolder, "Magpie", 'pipeline', paste0(section_name, ".Rmd")), 
  output_file = paste0(plotsdir, '/', date, "_", section_name, ".html")
)


## ---- SessionInfo

sessionInfo()