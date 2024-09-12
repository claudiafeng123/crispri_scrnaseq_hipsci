suppressMessages(library(Matrix))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))


# set relevent i/o paths
section_name <- "15b_downstream_heritability"
subsection_long_name <- "15b_03_downstream_heritability_write_dfs"

# set relevent i/o paths
date <- "2024-05-17"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
ExperimentName <- "Pica"

# options for testing
io_path <- paste0("scripts/io/", ExperimentName, "_io.R")
utils_path <- "scripts/io/Magpie_Utils.R"
control_tag <- c("unassigned", "NonTarget") # use all unassigned cells as control
keep_singleton_lines <- FALSE
make_negatives_zero <- TRUE
line_set_ind <- 10

args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  ##needed arguments
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--line_set_ind"){ line_set_ind <- as.numeric(args[[ind + 1]]) }
  
  ## not needed
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
  if (arg == "--utils_path"){ utils_path <- args[[ind + 1]] }
  if (arg == "--io_path"){ io_path <- args[[ind + 1]] }
  if (arg == "--make_negatives_zero"){ make_negatives_zero <- as.logical(args[[ind + 1]]) }
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

fnm <- sort(list.files(file.path(OutFolder, "15a_get_line_perms"), pattern = "_lines_per_target.tsv", full.names = T), decreasing = T)[1]
donor_sets <- fread(fnm )
if (exists("donors2include")){
  d <- sort(unlist(strsplit(donors2include, '_')))
} else {
  d <- sort(unlist(strsplit(donor_sets$lines[which(donor_sets$line_set_name == line_set_ind)], ";")))
}

outdir <- file.path(file.path(OutFolder, section_name))
if(!dir.exists(outdir)) dir.create(outdir, recursive = T)
plotsdir <- file.path(HTMLFolder, "pipeline", section_name, subsection_long_name, paste0('line_set_ind_', line_set_ind)) 
if(!dir.exists(plotsdir)) dir.create(plotsdir, recursive = T)


## ---- load data

fnm <- sort(list.files(file.path(OutFolder, section_name, 
                                   '15b_01_downstream_heritability_calc_downstream_meta'), 
                                   pattern = paste0(paste0(d, collapse = '_'), ".tsv.gz"), full.names = T),
            decreasing = T)[1]
downstream_meta <- fread(fnm)

downstream_gene_names <- downstream_meta$downstream_gene_name

ctrl_heritability_pvals <- mapply(1:length(downstream_gene_names), FUN = function(i){
  #print(i)
  downstream_gene_name <- downstream_gene_names[i]
  fnm <- list.files(file.path(OutFolder, 
                              section_name, 
                              "15b_02_downstream_heritability_calc_control_heritability", 
                              paste0('line_set_ind_', line_set_ind)),
                    pattern = paste0("_", downstream_gene_name, "-", '[0-9]'), full.names = T)
  
  if (length(fnm) == 1){
    h_perms <- fread(fnm)
    h <- downstream_meta$lfc_dLL_donor[i]
    
    if (!(is.na(h)) & !(is.infinite(h))){
      if (make_negatives_zero == T){
        if (h < 0){h <- 0}
        h_perms$lfc_dLL_donor[which(h_perms$lfc_dLL_donor < 0)] <- 0
      } 
      perms_beating_true_value <- length(which(h_perms$lfc_dLL_donor >= h))
      n_perms_tested <- dim(h_perms)[1]
      p_val <- (perms_beating_true_value + 1)/(n_perms_tested + 1)
    } else {
      p_val <- -Inf
    }
  } else {
    p_val <- -Inf
  }
  return(p_val)
})

rtn <- downstream_meta %>%
  mutate(pval=unlist(ctrl_heritability_pvals))
rtn$pval_adj <- p.adjust(rtn$pval, method = 'BH')

## ---- Write

fwrite(rtn, paste0(outdir, '/', date,  "_", paste0(d, collapse = '_'), ".tsv.gz"), sep = '\t', compress = 'gzip')

## ---- SessionInfo

sessionInfo()

## AARS had error at 840
## EXOSC2 had error at 2690

#fit1 <- lmer(lfc ~ on_target_expr + (1|cell_line) + (1|donor), df4lmm)
#fit2 <- lmer(lfc ~ on_target_expr + (1|cell_line), df4lmm)
#lrtest(fit1, fit2)

## ----

## add in single-cell variance
var_fnms <- list.files('/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Pica/6a_run_control_lm/2024-05-17/', full.names = T)
var_resid <- lapply(var_fnms, FUN = function(f){readRDS(f)$residuals_var})

var_resid <- data.frame(
  gene = unlist(lapply(strsplit(gsub(var_fnms, pattern = '-Gene-Expression.RDS|.*gene-', replacement = ''), '-'), FUN = function(x){paste(x[-1], collapse = '-')})),
  var_resid = unlist(var_resid)
)
fwrite(var_resid, '/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Pica/15b_downstream_heritability/2024-05-17_sc_resid_var.tsv', sep = '\t')

