suppressWarnings(suppressMessages(library(data.table, warn.conflicts = F)))
suppressWarnings(suppressMessages(library(dplyr, warn.conflicts = F)))
suppressWarnings(suppressMessages(library(ggplot2, warn.conflicts = F)))
suppressWarnings(suppressMessages(library(dplyr, warn.conflicts = F)))
suppressWarnings(suppressMessages(library(ggpubr, warn.conflicts = F)))
suppressWarnings(suppressMessages(library(ggExtra, warn.conflicts = F)))
suppressWarnings(suppressMessages(library(doParallel, warn.conflicts = F)))

# set relevent i/o paths
section_name <- "11f_target_complex_cor_calc_cor"

# set relevent i/o paths
date <- "2022-08-15"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
ExperimentName <- "Magpie"
#REDO <- T
complex_ind <- 1

args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  ##needed arguments
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--experiment_name"){ ExperimentName <- args[[ind + 1]] }
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
  if (arg == "--complex_ind"){ complex_ind <- as.numeric(args[[ind + 1]]) }
  
  
  ### not needed
  if (arg == "--utils_path"){ utils_path <- args[[ind + 1]] }
  if (arg == "--io_path"){ io_path <- args[[ind + 1]] }
  
}

setwd(HomeFolder)
if (!(exists("io_path"))){io_path <- paste0(HomeFolder, "scripts/io/", ExperimentName, "_io.R")}
if (!(exists("utils_path"))){ utils_path <- paste0(HomeFolder, "/scripts/io/Magpie_Utils.R")}
source(io_path)
source(utils_path)
control_tag <- c("unassigned", "NonTarget") # use all unassigned cells as control

print(paste("date:", date))
print(paste("home folder:", HomeFolder))
print(paste("Experiment Name:", ExperimentName))
print(paste("min_cells_per_gene:", min_cells_per_gene))
print(paste("section_name:", section_name))


## set io
#config <- paste0(date, "_min_degs-", min_deg_per_target)
outdir <- file.path(OutFolder, section_name, "by_complex/")
if(!dir.exists(outdir)) dir.create(outdir, recursive = T)
plotsdir <- file.path(HTMLFolder, "pipeline", section_name)
if(!dir.exists(plotsdir)) dir.create(plotsdir, recursive = T)

#if (min_deg_per_target <= 10){min_deg_per_target <- 10}


## ---- CalcCorr


complex_meta <- fread(paste0(OutFolder, "11d_get_clusters_to_aggregate/", date, "_complex_meta.tsv"))#fread(complex_meta_path)

fnm <- paste0(OutFolder, "11e_calc_lfcs_transcriptome_wide_by_complex/by_gene/", complex_ind, "_", date, ".tsv.gz")
if(file.exists(fnm)){
  complex_lfc <- fread(fnm)
  complex_lfc <- mutate(complex_lfc, complex_name = complex_meta$complex_name[complex_ind],
                downstream_gene_name = unlist(lapply(strsplit(downstream_gene, ":"), "[[", 2)))
} else {complex_lfc <- NULL}


lfcs_by_target_split <- fread(paste0(OutFolder, "6b_calc_lfcs_transcriptome_wide_by_gene", "/", date, "_combined_with_adj_pval.tsv.gz"))
lfcs_by_target_split <- split.data.frame(lfcs_by_target_split, f = lfcs_by_target_split$target)
target_meta <- fread(file.path(OutFolder, '7b_target_summary', paste0(date,"_target_meta_data.csv")))

all_downstream_genes <- lfcs_by_target_split[[1]]$downstream_gene_name

## ---- Target_Complex_Cor

fnm_out <- paste0(outdir, "/", date, "_", complex_ind, "_target_complex_cor.tsv")

if (!(is.null(complex_lfc))){
  rtn <- lapply(1:length(lfcs_by_target_split), FUN = function(i){
    #print(paste0("target: ", lfcs_by_target_split[[i]]$target[1]))
    #print(paste0("complex: ", complex_meta$complex_name[complex_ind]))
    cor_by_target <- calc_cor(df1 = lfcs_by_target_split[[i]],
                              df2 = complex_lfc, by_pval = F,
                              is_2_complex = T,
                              xlab1 = lfcs_by_target_split[[i]]$target[1],
                              ylab1 = complex_meta$complex_name[complex_ind], anno = T)
  })
  target_complex_cor_df <- lapply(rtn, FUN = function(l){
    rtn <- l$df
  }) %>% bind_rows() %>% as.data.frame()
  fwrite(target_complex_cor_df, fnm_out, sep = "\t")
  
  ## plot top hits
  pdf(paste0(plotsdir, "/",date, "_", complex_ind, "_top_correlates.pdf"), width = 5, height = 5)
  top_correlates <- filter(target_complex_cor_df, cor_all > 0.3) %>% .$downstream_gene_name
  plots <- lapply(which(abs(target_complex_cor_df$cor_all) > 0.3), FUN = function(x){rtn[[x]]$p})
  print(plots)
  dev.off()
} 



