
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(tidyverse))



# set relevent i/o paths
section_name <- "15c_calc_var_expl_by_target"

# set relevent i/o paths
date <- "2024-01-10"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
ExperimentName <- "Pica"
subsection_name <- "15c_04_calc_var_expl_by_target_get_permutation_status"

# options for testing
io_path <- paste0("scripts/io/", ExperimentName, "_io.R")
utils_path <- "scripts/io/Magpie_Utils.R"
control_tag <- c("unassigned", "NonTarget") # use all unassigned cells as control
keep_singleton_lines <- FALSE
include_line_quality <- T
include_sex <- T
num_perms_to_beat <- 10
target_gene <- "PAF1"; dg <- "PRKD3"
make_negatives_zero <- T
val2calc <- "lfc" #'lfc' or 'post_kd_expr'

args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  ##needed arguments
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--target_gene"){ target_gene <- args[[ind + 1]] }
  if (arg == "--downstream_gene_name"){ downstream_gene_name <- args[[ind + 1]] }
  if (arg == "--val2calc"){ val2calc <- args[[ind + 1]] }
  
  ## not needed
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
  if (arg == "--subsection_name"){ subsection_name <- args[[ind + 1]] }
  if (arg == "--utils_path"){ utils_path <- args[[ind + 1]] }
  if (arg == "--io_path"){ io_path <- args[[ind + 1]] }
  
  ## linear model parameters
  if (arg == "--keep_singleton_lines"){ keep_singleton_lines <- as.logical(args[[ind + 1]]) }
  if (arg == "--include_line_quality"){ include_line_quality <- as.logical(args[[ind + 1]]) }
  if (arg == "--include_sex"){ include_sex <- as.logical(args[[ind + 1]]) }
  
  ## parameters
  if (arg == "--num_perms_to_beat"){ num_perms_to_beat <- as.numeric(args[[ind + 1]]) }
  if (arg == "--make_negatives_zero"){ make_negatives_zero <- as.logical(args[[ind + 1]]) }
  
}

print(paste0(section_name))
print(paste("date:", date))
print(paste("home folder:", HomeFolder))


setwd(HomeFolder)
source(io_path)
source(utils_path)

GuideMetadata <- fread(GuideMetadataPath)
LineMetadata <- fread(LineMetadataPath)
sexes <- LineMetadata$sex; names(sexes) <- LineMetadata$name


outdir <- file.path(file.path(OutFolder, section_name, subsection_name, val2calc, target_gene))
if(!dir.exists(outdir)) dir.create(outdir, recursive = T)
plotsdir <- file.path(HTMLFolder, "pipeline", section_name) 
if(!dir.exists(plotsdir)) dir.create(plotsdir, recursive = T)

fnm_out <- paste0(outdir, '/', date, "_", target_gene, "_dg_permutation_status.tsv")

## ---- LoadData

folder_name <- list.files(file.path(OutFolder, section_name), pattern = "15c_02")
most_recent <- sort(unlist(lapply(strsplit(list.files(file.path(OutFolder, section_name, folder_name)), "_"), "[[", 1)), decreasing = T)[1]
dLL <- fread(file.path(OutFolder, section_name, folder_name, paste0(most_recent, "_var_expl_combined_no_pval.tsv.gz"))) %>%
  dplyr::select(c("target", "downstream_gene_name", paste0(val2calc, '_dLL_donor'))) %>%
  dplyr::filter(target == target_gene)

# all pairs permuted
folder_name <- list.files(file.path(OutFolder, section_name), pattern = "15c_00")
most_recent <- sort(unlist(lapply(strsplit(list.files(file.path(OutFolder, section_name, folder_name)), "_"), "[[", 1)), decreasing = T)[1]
pairs_permuted <- fread( paste0(OutFolder, section_name, '/', folder_name, '/', most_recent, "_pairs_to_permute.tsv"))

## ---- Consolidate Resultsx


most_recent <- sort(unlist(lapply(strsplit(list.files(file.path(OutFolder, section_name, "15c_03_calc_var_expl_by_target_permute_lines_all_pairs", val2calc, target_gene)), "_"), "[[", 1)), decreasing = T)[1]
completed_downstream_genes_fnms <- grep(list.files(file.path(OutFolder, section_name, "15c_03_calc_var_expl_by_target_permute_lines_all_pairs", val2calc, target_gene), 
                                                   pattern = '[.]tsv[.]gz'), pattern = most_recent, value = T)
completed_downstream_genes <- unlist(lapply(strsplit(
  gsub(completed_downstream_genes_fnms, pattern = paste0(most_recent, "_", target_gene, "_"), replacement = ""), "_"), "[[", 1))
if (length(completed_downstream_genes) < dim(pairs_permuted %>% filter(target == target_gene))[1]){print("some pairs missing")}

#pdf(paste0(plotsdir, '/', date, "_", target_gene, "_lfc_permutation_dLL_dist.pdf"), width = 4.5, height = 3)
#i <- grep(completed_downstream_genes_fnms, pattern = dg)
permutation_status <- mapply(1:length(completed_downstream_genes_fnms), FUN = function(i){
  #print(i)
  dg <- completed_downstream_genes[i]
  #print(dg)
  eval(parse(text = paste0("h <- dLL %>% filter(downstream_gene_name == dg) %>% .$", val2calc, "_dLL_donor")))
  
  df2write <- fread(paste0(file.path(OutFolder, section_name, "15c_03_calc_var_expl_by_target_permute_lines_all_pairs", val2calc, target_gene, completed_downstream_genes_fnms[i])))
  names(df2write) <- gsub(names(df2write), pattern = "var_expl", replacement = "dLL")
  
  num_non_converging <- length(which(is.infinite(df2write$LL_m2)))
  df2write <- df2write %>% filter(!is.infinite(lfc_dLL_donor))
  
  
  if (is.na(h) | dim(df2write)[1] == 0){
    rtn <- c(
      permutations_tested = 0, 
      max_dLL = -Inf,
      perms_with_higher_dLL = 0, 
      perms_with_equal_dLL = 0, 
      perms_with_negative_dLL = 0,
      perms_with_zero_dLL = 0,
      perms_that_didnt_converge = 0,
      num_repeated_permuted_hs = 0,
      dLL_mode = -Inf,
      num_mode_appearances = 0,
      lfc_dLL_donor = -Inf,
      adjusted_h = -Inf,
      perms_beating_true_value = -Inf,
      pval = -Inf
    )
    return(rtn)
  }
  
  ## plot
  #p_nonzero <- ggplot(df2write, aes(x = lfc_dLL_donor)) + 
  #  xlab("Change in Log-likelihood") + ylab("# of Permutations") + ggtitle(paste0(dg, " Expression in ", target_gene, " Knockdowns")) + 
  #  xlim(c(min(c(df2write$lfc_dLL_donor, h, -0.1)), max(c(0.1, df2write$lfc_dLL_donor, h)*1.1))) + 
  #  geom_histogram(bins = 50) + 
  #  geom_vline(xintercept = h, col = 'red', lty = 3) + 
  #  theme_bw()
  #print(p_nonzero)
  
  freq_table <-  table(df2write$lfc_dLL_donor)
  freq_table <- freq_table[which(!names(freq_table) == "0")]
  num_repeated_dLLs <- length(freq_table[which(freq_table > 1)])
  if (num_repeated_dLLs >= 1){
    freq_table <- freq_table[which.max(freq_table)]
    freq_table <- freq_table[which.max(abs(as.numeric(names(freq_table))))]
  } else{
    freq_table <- 0
    names(freq_table) <- 0
  }
  
  rtn <- c(
    permutations_tested = dim(df2write)[1], 
    max_dLL = max(df2write$lfc_dLL_donor)[1],
    perms_with_higher_dLL = length(which(df2write$lfc_dLL_donor > h)), 
    perms_with_equal_dLL = length(which(df2write$lfc_dLL_donor == h)), 
    perms_with_negative_dLL = length(which(df2write$lfc_dLL_donor < 0)),
    perms_with_zero_dLL = length(which(df2write$lfc_dLL_donor == 0)),
    perms_that_didnt_converge = num_non_converging,
    num_repeated_permuted_hs = num_repeated_dLLs,
    dLL_mode = names(freq_table),
    num_mode_appearances = as.numeric(freq_table),
    lfc_dLL_donor = h)
  if (val2calc == "post_kd_expr"){names(rtn)[length(rtn)] <- "post_kd_expr_dLL_donor"}
  
  if (make_negatives_zero == T){
    if (h < 0){h <- 0}
    df2write$lfc_dLL_donor[which(df2write$lfc_dLL_donor < 0)] <- 0
  } 
  rtn['adjusted_h'] <- h
  rtn['perms_beating_true_value'] <- length(which(df2write$lfc_dLL_donor >= h))
  rtn['pval'] <- (as.numeric(rtn['perms_beating_true_value']) + 1)/(as.numeric(rtn['permutations_tested']) + 1)
  
  return(rtn)
}, SIMPLIFY = F)
#dev.off()

permutation_status <- data.frame(
  to_match =paste0(target_gene, "_", completed_downstream_genes),
  target = target_gene,
  downstream_gene_name = completed_downstream_genes,
  permutations_tested = as.numeric(unlist(lapply(permutation_status, "[[", 1))),
  max_dLL = as.numeric(unlist(lapply(permutation_status, "[[", 2))),
  perms_with_higher_dLL = as.numeric(unlist(lapply(permutation_status, "[[", 3))),
  perms_with_equal_dLL = as.numeric(unlist(lapply(permutation_status, "[[", 4))),
  perms_with_negative_dLL = as.numeric(unlist(lapply(permutation_status, "[[", 5))),
  perms_with_zero_dLL = as.numeric(unlist(lapply(permutation_status, "[[", 6))),
  perms_that_didnt_converge = as.numeric(unlist(lapply(permutation_status, "[[", 7))),
  num_repeated_permuted_hs = as.numeric(unlist(lapply(permutation_status, "[[", 8))),
  dLL_mode = as.numeric(unlist(lapply(permutation_status, "[[", 9))),
  num_mode_appearances = as.numeric(unlist(lapply(permutation_status, "[[", 10))),
  lfc_dLL_donor = as.numeric(unlist(lapply(permutation_status, "[[", 11))),
  adjusted_h = as.numeric(unlist(lapply(permutation_status, "[[", 12))),
  perms_beating_true_value = as.numeric(unlist(lapply(permutation_status, "[[", 13))),
  pval = as.numeric(unlist(lapply(permutation_status, "[[", 14)))
)
if (val2calc == "post_kd_expr" ){names(permutation_status) <- gsub(names(permutation_status), pattern = "lfc", replacement = val2calc)}

## ---- How much more to do?

fwrite(permutation_status, paste0(outdir, '/', date, "_", target_gene, "_permutation_status.tsv"), sep = '\t')

## ---- SessionInfo

sessionInfo()

