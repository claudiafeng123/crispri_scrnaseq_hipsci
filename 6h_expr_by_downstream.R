suppressWarnings(suppressMessages(library(Matrix, warn.conflicts = F)))
suppressWarnings(suppressMessages(library(Seurat, warn.conflicts = F)))
suppressWarnings(suppressMessages(library(data.table, warn.conflicts = F)))
suppressWarnings(suppressMessages(library(dplyr, warn.conflicts = F)))
suppressWarnings(suppressMessages(library(ggplot2, warn.conflicts = F)))
suppressWarnings(suppressMessages(library(tidyverse, warn.conflicts = F)))
suppressWarnings(suppressMessages(library(parallel, warn.conflicts = F)))
suppressWarnings(suppressMessages(library(doParallel, warn.conflicts = F)))



## for pica
date <- "2024-05-17"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
ExperimentName <- "Pica"
no_genes2test <- 10
start_ind <- 0

args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  ##needed arguments
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--experiment_name"){ ExperimentName <- args[[ind + 1]] }
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
  
  if (arg == "--start_ind"){ start_ind <- as.numeric(args[[ind + 1]]) }
  if (arg == "--no_genes2test"){ no_genes2test <- as.numeric(args[[ind + 1]]) }
  if (arg == "--max_unperturbed"){ max_unperturbed <- as.numeric(args[[ind + 1]]) }
  
}

#scale_factor <- scale_factor*10^5

if (!(exists('out_date'))){out_date <- date}
print(paste("date:", date))
print(paste("home folder:", HomeFolder))

setwd(HomeFolder)
if (!(exists("io_path"))){io_path <- paste0(HomeFolder, "scripts/io/", ExperimentName, "_io.R")}
if (!(exists("utils_path"))){ utils_path <- paste0(HomeFolder, "/scripts/io/Magpie_Utils.R")}
source(io_path)
source(utils_path)

outdir <- paste0(OutFolder, "6b_calc_lfcs_transcriptome_wide_by_gene/by_downstream/")
if (!(dir.exists(outdir))){dir.create(outdir, recursive = T)}

## ---- LoadData


perturbed_fnms <- grep(list.files(file.path(OutFolder, "6b_calc_lfcs_transcriptome_wide_by_gene", 'by_target'), full.names = T),
                        pattern = NonTargetGeneName, invert = T, value = T)
nontarget_fnms <- list.files(file.path(OutFolder, '6f_calc_lfcs_transcriptome_wide_by_gene_per_line'), 
                             pattern = paste0(date, "_", NonTargetGeneName, "_knockdowns_residuals"), recursive = T, full.names = T)
unassigned_fnms <- grep(list.files(file.path(OutFolder, '6f_calc_lfcs_transcriptome_wide_by_gene_per_line'), 
                                   pattern = paste0("unassigned"), recursive = T, full.names = T),
                        pattern = "knockdowns_residuals", value = T)
fnms <- c(perturbed_fnms, nontarget_fnms, unassigned_fnms)
lfcs <- lapply(fnms, FUN = function(fnm){
  #print(fnm)
  if (!grepl(fnm, pattern = 'unassigned') & !(grepl(fnm, pattern = NonTargetGeneName))){
    rtn <- fread(fnm)
  } else {
    cell_line <- gsub(fnm, pattern = ".*/6f_calc_lfcs_transcriptome_wide_by_gene_per_line/|/by_target_by_line/.*", replacement = '')
    rtn <- fread(fnm) %>%
      mutate(cell_line = cell_line,
             id = paste0(cell_line, "_", id))
    if ("V1" %in% colnames(rtn)){
      rtn <- select(rtn, !(contains("V1")))
    }
  }
  return(rtn)
})
cell_ids <- lapply(lfcs, FUN = function(df){df[,1]}) %>% unlist()
genes2test <- grep(colnames(lfcs[[1]]), pattern = paste0(':', "PRKD3", ':Gene-Expression'), value = T)
#genes2test <- colnames(lfcs[[1]])[-1][(start_ind + 1):(min(length(colnames(lfcs[[1]])) - 1, (start_ind+no_genes2test)))]

pica_cell_meta_fnms <- list.files(paste0(OutFolder, '5a_combine'), pattern = 'cells_metadata', full.names = T)
pica_cell_meta <- lapply(pica_cell_meta_fnms, FUN = fread) %>% bind_rows() %>% as.data.frame()
pica_cell_meta <- pica_cell_meta %>%
  dplyr::filter(paste0(cell_line, "_", V1) %in% cell_ids)


## ---- ExpressionFiles

registerDoParallel(10)

foreach (downstream_gene_name =genes2test) %dopar% {
  expr <- lapply(lfcs, FUN = function(df){
    eval(parse(text = paste0('df$`', downstream_gene_name, '`')))
  }) %>% unlist()
  control_lm <- readRDS(paste0(OutFolder, "6a_run_control_lm/", date,
                               "/gene-", gsub(gsub(downstream_gene_name, pattern = ':', replacement = '-'),
                                              pattern = '/', replacement = '-'), ".RDS"))
  rtn <- pica_cell_meta %>%
    dplyr::select(!(contains(c("CellRanger", 'guide_call', 'snn', 'cluster', 'num_features')))) %>%
    mutate(regressed_lfc = expr[match(paste0(cell_line, "_", V1), cell_ids)],
           cell_line_coef = ifelse(paste0('cell_line', cell_line) %in% names(control_lm$coefficients), 
                                   control_lm$coefficients[paste0('cell_line', cell_line)], 0 )) %>%
    mutate(regressed_expr = regressed_lfc + cell_line_coef)
  fwrite(rtn, paste0(outdir, date, "_", gsub(gsub(downstream_gene_name, pattern = ':', replacement = '-'),
                                  pattern = '/', replacement = '-'), '.tsv.gz' ),
         sep = '\t', compress = 'gzip' )
}


#sessionInfo()

