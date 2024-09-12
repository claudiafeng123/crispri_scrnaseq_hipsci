suppressMessages(library(umap))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(tidyverse))
suppressMessages(library( org.Hs.eg.db))

# set relevent i/o paths
section_name <- "15c_calc_var_expl_by_target"
subsection_name <- "15c_02_calc_var_expl_by_target_variance_decomposition"
## combine the by-target results
## produce p-values

# set relevent i/o paths
date <- "2024-06-03"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
ExperimentName <- "Pica"

# options for testing
io_path <- paste0("scripts/io/", ExperimentName, "_io.R")
utils_path <- "scripts/io/Magpie_Utils.R"
REDO <- F
donors2include <- 'eipl_iudw_jejf_kolf_paab'

args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  ##needed arguments
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--target_gene"){ target_gene <- args[[ind + 1]] }
  if (arg == "--donors2include"){ donors2include <- args[[ind + 1]] }
  if (arg == "--REDO"){ REDO <- as.logical(args[[ind + 1]]) }
  
  ## not needed
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
  if (arg == "--utils_path"){ utils_path <- args[[ind + 1]] }
  if (arg == "--io_path"){ io_path <- args[[ind + 1]] }
  
  ## permuting
  if (arg == "--keep_singleton_lines"){ keep_singleton_lines <- as.logical(args[[ind + 1]]) }
  if (arg == "--include_line_quality"){ include_line_quality <- as.logical(args[[ind + 1]]) }
  if (arg == "--include_sex"){ include_sex <- as.logical(args[[ind + 1]]) }
}

print(paste0(section_name))
print(paste("date:", date))
print(paste("home folder:", HomeFolder))
if(exists("donors2include")){print('including donors: '); print(donors2include)}



setwd(HomeFolder)
source(io_path)
source(utils_path)


outdir <- file.path(file.path(OutFolder, section_name, subsection_name))
if(!dir.exists(outdir)) dir.create(outdir, recursive = T)
plotsdir <- file.path(HTMLFolder, "pipeline", section_name, subsection_name) 
if(!dir.exists(plotsdir)) dir.create(plotsdir, recursive = T)

## ---- Load Other Data

#no_fixed_effects <- fread(paste0(outdir, '/', date, '_var_expl_no_fixed_effects_combined_no_pval.tsv.gz'))
fnm <- sort(list.files(file.path(OutFolder, "6f_calc_lfcs_transcriptome_wide_by_gene_per_line"), pattern = "_all_lines_with_adj_pval.tsv.gz", full.names = T), decreasing = T)[1]
lfc_by_line <- fread(fnm )



## ---- true values with fixed effects

fnm_out <- paste0(outdir, '/', date, '_var_expl_combined_no_pval.tsv.gz')
if (REDO == T | !(file.exists(fnm_out))){
  
  
  
  var_expl_fnms <- list.files(file.path(OutFolder, section_name, '15c_01_calc_var_expl_by_target_run_lme'), 
                              pattern = paste0(date, "_permute-cells-FALSE_permute-lines-FALSE_include-line-quality-TRUE_include-sex-FALSE"), 
                              full.names = T, recursive = T)
  var_expl <- lapply(var_expl_fnms, fread) %>%  bind_rows() %>% as.data.frame()
  var_expl <- var_expl %>% 
    dplyr::select(c('target', 'downstream_gene_name', 
                    "bulk_mean_lfc", "bulk_lfc_var", "max_abs_lfc", "bulk_range_lfc",
                    'lfc_dLL_donor', 'lfc_dLL_line', 'lfc_LL_m0',
                    'lfc_varExpl_cell_line', 'lfc_varExpl_donor', 'lfc_varExpl_guide', 'lfc_varExpl_on_target_expr', 'lfc_varExpl_Residual',
                    'lfc_coef_on_target_expr', 'lfc_tval_on_target_expr', 
                    'post_kd_expr_dLL_donor', 'post_kd_expr_dLL_line', 'post_kd_expr_LL_m0'))
  
  
  ## add in control expression heritability
  fnm <- sort(list.files(file.path(OutFolder, "15a_get_line_perms"), pattern = '_target_line_mapping', full.names = T), decreasing = T)[1]
  lines_per_target <- fread(fnm)
  if (exists('donors2include')){
    ds <- unlist(strsplit(donors2include, '_'))
    fnm <- sort(list.files(file.path(OutFolder, "15a_get_line_perms"), pattern = '_lines_per_target', full.names = T), decreasing = T)[1]
    line_ind_table <- fread(fnm)
    lines_per_target$paired_lines <- unlist(lapply(1:length(lines_per_target$paired_lines), FUN = function(i){
      well_powered_lines <- unlist(strsplit(lines_per_target$paired_lines[i], ";"))
      well_powered_lines <- well_powered_lines[which(well_powered_lines %in% ds)]
      well_powered_lines <- paste(well_powered_lines, collapse = ';')
      return(well_powered_lines)
    }))
    lines_per_target$n_paired_lines <- unlist(lapply(strsplit(lines_per_target$paired_lines, ";"), length))
    lines_per_target$line_set_name <- line_ind_table$line_set_name[match(lines_per_target$paired_lines, line_ind_table$lines)]
  }
  fnm <- sort(list.files(file.path(OutFolder, "7b_target_summary"), pattern = '_target_meta_data', full.names = T), decreasing = T)[1]
  target_meta <- fread(fnm)
  ## split by line set index
  target_sets <- split.data.frame(lines_per_target, f = lines_per_target$line_set_name)
  var_expl <- mapply(1:length(target_sets), FUN = function(line_set_ind){
    targets_in_set <- target_sets[[line_set_ind]]$gene
    n_paired_lines <- target_sets[[line_set_ind]]$n_paired_lines[1]
    
    ## control expression
    fnm <- sort(list.files(file.path(OutFolder, "15b_downstream_heritability"), pattern = paste0(paste(unlist(strsplit(target_sets[[line_set_ind]]$paired_lines[1], ';')), collapse = '_'),'.tsv.gz'), full.names = T), decreasing = T)[1]
    downstream_meta <- fread(fnm)
    df2merge_downstream_meta <- downstream_meta %>%
      mutate(n_paired_lines = n_paired_lines) %>%
      dplyr::select(c('downstream_gene_name',
                      n_paired_lines,
                      'bulk_mean_control_expr', 'bulk_control_expr_var', 
                      'control_expr_dLL_line' = 'lfc_dLL_line', 'control_expr_dLL_donor' = 'lfc_dLL_donor', 'control_expr_LL_m0' = 'LL_m0',
      ))
    
    df2merge_heritability_df <- var_expl %>%
      dplyr::filter(target %in% targets_in_set)
    rtn <- left_join(df2merge_heritability_df, df2merge_downstream_meta)
    
  }, SIMPLIFY = F) %>% bind_rows() %>% as.data.frame()
  
  
  ## add in genomic coordinates
  downstream_gene_entrez_ids <- select(org.Hs.eg.db, 
                                       keys = unique(var_expl$downstream_gene_name),
                                       columns = c("ENTREZID", "SYMBOL"),
                                       keytype = "SYMBOL") %>%
    dplyr::filter(!(is.na(ENTREZID)))
  chr_meta <- as.list(org.Hs.egCHR[downstream_gene_entrez_ids$ENTREZID])
  downstream_gene_entrez_ids$CHROM <- unlist(lapply(chr_meta[match(downstream_gene_entrez_ids$ENTREZID, names(chr_meta))], FUN = function(x){paste(x, collapse = "_")}))
  var_expl$chrom <- downstream_gene_entrez_ids$CHROM[match(var_expl$downstream_gene_name, downstream_gene_entrez_ids$SYMBOL)]
  
  fwrite(var_expl, fnm_out, sep = '\t', compress = 'gzip')
} else {
  var_expl <- fread(fnm_out)
  
}

## ---- MakePlots

## load some more data
if (F){
  rmarkdown::render(file.path(CodeFolder, "Magpie", "pipeline", paste0(subsection_name, ".Rmd") ),
                    output_file = paste0(plotsdir, '/', date, '.html'),
                    params = list(
                      date =date,
                      section_name = section_name,
                      subsection_name = subsection_name
                    ))
  
}

## ---- SessionInfo

sessionInfo()
#mean_effect <- fread(file.path(OutFolder, "6b_calc_lfcs_transcriptome_wide_by_gene", paste0(date, "_combined_with_adj_pval.tsv.gz")))
if (F){
  lines_per_target_fnm <- sort(list.files(paste0(OutFolder, "/15a_get_line_perms/"), pattern = "_lines_per_target.tsv", full.names = T), decreasing = T)[1]
  lines_per_target <- fread(lines_per_target_fnm)
  target_meta_fnm <- sort(list.files(file.path(OutFolder, '7b_target_summary'), pattern = "_target_meta_data.csv", full.names = T), decreasing = T)[1]
  target_meta <- fread(target_meta_fnm)
  target_meta_by_line_fnm <- sort(list.files(file.path(OutFolder, '9b_target_gene_downregulation_by_line'), pattern = "_target_meta_by_line.csv", full.names = T), decreasing = T)[1]
  target_meta_by_line <- fread(target_meta_by_line_fnm)
  
}

