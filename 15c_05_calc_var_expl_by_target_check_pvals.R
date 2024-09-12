## mainly exploration
## adjusted p-values can be computed in the next section. this is to decide a cutoff
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(ggpointdensity))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(tidyverse))
suppressMessages(library(lme4))


# set relevent i/o paths
section_name <- "15c_calc_var_expl_by_target"
subsection_name <- "15c_05_calc_var_expl_by_target_check_pvals"

# set relevent i/o paths
date<- "2024-09-04"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
ExperimentName <- "Pica"

# options for testing
io_path <- paste0("scripts/io/", ExperimentName, "_io.R")
utils_path <- "scripts/io/Magpie_Utils.R"
control_tag <- c("unassigned", "NonTarget") # use all unassigned cells as control
keep_singleton_lines <- FALSE
include_line_quality <- T
include_sex <- T
num_perms_to_beat <- 10
max_perms2compute <- 10^4
REDO <- T
donors2include <- 'eipl_iudw_jejf_kolf_paab'

args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  ##needed arguments
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--date_new"){ date_new <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--donors2include"){ donors2include <- args[[ind + 1]] }
  if (arg == "--REDO"){ REDO <- as.logical(args[[ind + 1]]) }
  
  ## not needed
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
  if (arg == "--utils_path"){ utils_path <- args[[ind + 1]] }
  if (arg == "--io_path"){ io_path <- args[[ind + 1]] }
  
  ## linear model parameters
  if (arg == "--keep_singleton_lines"){ keep_singleton_lines <- as.logical(args[[ind + 1]]) }
  if (arg == "--include_line_quality"){ include_line_quality <- as.logical(args[[ind + 1]]) }
  if (arg == "--include_sex"){ include_sex <- as.logical(args[[ind + 1]]) }
  
  ## parameters
  if (arg == "--num_perms_to_beat"){ num_perms_to_beat <- as.numeric(args[[ind + 1]]) }
  if (arg == "--max_perms2compute"){ max_perms2compute <- as.numeric(args[[ind + 1]]) }
  
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

outdir <- file.path(file.path(OutFolder, section_name, subsection_name,'/'))
if(!dir.exists(outdir)) dir.create(outdir, recursive = T)
plotsdir <- file.path(HTMLFolder, "pipeline", section_name, subsection_name,'/') 
if(!dir.exists(plotsdir)) dir.create(plotsdir, recursive = T)

## ---- LoadData

crispr_sequences <- gsub(list.files(file.path(ResourcesFolder, 'reference_sequences', 'added_sequences'), pattern = '[.]gtf'), pattern = '.gtf', replacement = '')

fnm <- sort(grep(list.files(file.path(OutFolder, "7b_target_summary"), pattern = "_target_meta_data.csv", full.names = T), pattern = "params", invert = T, value = T), decreasing =  T)[1]
target_meta <- fread(fnm)

## heritability
folder_name <- list.files(file.path(OutFolder, section_name), pattern = "15c_02")
most_recent <- sort(unlist(lapply(strsplit(list.files(file.path(OutFolder, section_name, folder_name)), "_"), "[[", 1)), decreasing = T)[1]
heritability_df <- fread(file.path(OutFolder, "15c_calc_var_expl_by_target", folder_name, paste0(most_recent, "_var_expl_combined_no_pval.tsv.gz")))

## add in p-values
## LFC
perm_status_fnms <- list.files(file.path(OutFolder, section_name, '15c_04_calc_var_expl_by_target_get_permutation_status', 'lfc'), pattern = "permutation_status.tsv", full.names = T, recursive = T)
lfc_pvals <- lapply(perm_status_fnms, FUN = fread) %>% bind_rows() %>% as.data.frame()
lfc_pvals <- left_join(lfc_pvals , heritability_df %>% 
                            dplyr::select(c("target", "downstream_gene_name", 'n_paired_lines')))
#post_kd_expr
most_recent <- sort(gsub(unlist(lapply(strsplit(
  list.files(file.path(OutFolder, section_name, '15c_04_calc_var_expl_by_target_get_permutation_status', 'post_kd_expr'), pattern = "permutation_status.tsv", recursive = T), "_"), "[[", 1)),
  pattern = '.*/', replacement = ""), decreasing = T)[1]
perm_status_fnms <- grep(list.files(file.path(OutFolder, section_name, '15c_04_calc_var_expl_by_target_get_permutation_status', 'post_kd_expr'), pattern = "permutation_status.tsv", full.names = T, recursive = T), pattern = most_recent, value = T)
post_kd_expr_pvals <- lapply(perm_status_fnms, FUN = fread) %>% bind_rows() %>% as.data.frame()
post_kd_expr_pvals  <- left_join(post_kd_expr_pvals  , heritability_df %>% 
                         dplyr::select(c("target", "downstream_gene_name", 'n_paired_lines')))

off_target <- fread(paste0(OutFolder, "/../Magpie/Preprocess_External_Datasets/off_target/", "off_target_pairs_magpie-2024-05-17-version.tsv.gz"), header = T) %>%
  dplyr::filter(target != off_target)

## crossmapping
cross_mapped_pairs <- fread(paste0(ResourcesFolder, 'eQTLs/crossmapping_reads_iPSC_PairedEndData_with_gene_name.txt'))

## ---- DataPreprocessing

heritability_fnm_out <- paste0(outdir, '/', date, "_heritability_all_no_pval_adj.tsv.gz")


## add in target metadata
target_meta2add <- target_meta %>% dplyr::select(c('target' = 'gene',
                                                   'magpie_deg', 'n_downstream_excl_target',
                                                   contains('is_')))
heritability_df <- left_join(heritability_df,  target_meta2add)

## add a column for off-target
heritability_df$potential_off_target <- ifelse(
  paste0(heritability_df$target, ":", heritability_df$downstream_gene_name) %in% 
    paste0(off_target$target, ":", off_target$off_target), "y", "n")
heritability_df$potential_cross_map <- ifelse(
  paste0(heritability_df$target, ":", heritability_df$downstream_gene_name) %in% 
    paste0(cross_mapped_pairs$gene_1, ":", cross_mapped_pairs$gene_2), "y", "n")

## lfc and post-kd expression
heritability_df <- left_join(left_join(heritability_df,  
                                       lfc_pvals %>% dplyr::select(c("target", "downstream_gene_name",
                                                                     'lfc_dLL_donor_pval' = 'pval'))), 
                                       post_kd_expr_pvals %>% 
                                         dplyr::select(c("target", "downstream_gene_name",
                                                         'post_kd_expr_dLL_donor_pval' = 'pval')))
                             
## control expression
## add in wild-type expression
## add in control expression heritability
fnm <- sort(list.files(file.path(OutFolder, "15a_get_line_perms"), pattern = '_target_line_mapping', full.names = T), decreasing = T)[1]
line_ind_by_target <- fread(fnm)
if (exists('donors2include')){
  ds <- unlist(strsplit(donors2include, '_'))
  fnm <- sort(list.files(file.path(OutFolder, "15a_get_line_perms"), pattern = '_lines_per_target', full.names = T), decreasing = T)[1]
  line_ind_table <- fread(fnm)
  line_ind_by_target$paired_lines <- unlist(lapply(1:length(line_ind_by_target$paired_lines), FUN = function(i){
    well_powered_lines <- unlist(strsplit(line_ind_by_target$paired_lines[i], ";"))
    well_powered_lines <- well_powered_lines[which(well_powered_lines %in% ds)]
    well_powered_lines <- paste(well_powered_lines, collapse = ';')
    return(well_powered_lines)
  }))
  line_ind_by_target$n_paired_lines <- unlist(lapply(strsplit(line_ind_by_target$paired_lines, ";"), length))
  line_ind_by_target$line_set_name <- line_ind_table$line_set_name[match(line_ind_by_target$paired_lines, line_ind_table$lines)]
}


targets_by_donors_with_data <- split.data.frame(line_ind_by_target, f = line_ind_by_target$line_set_name)
wt_expr_df <- lapply(targets_by_donors_with_data, FUN = function(df){
  line_set_ind <- df$line_set_name[1]
  #print(line_set_ind)
  lines_in_set <- sort(unlist(strsplit(df$paired_lines[1], ';')))
  most_recent <- sort(unique(unlist(lapply(strsplit(list.files(file.path(OutFolder, "15b_downstream_heritability"), pattern = '.tsv'), "_"), "[[", 1))), decreasing = T)
  wt_downstream_meta <- fread(file.path(OutFolder, "15b_downstream_heritability", paste0(most_recent, '_', paste(lines_in_set, collapse = '_'), '.tsv.gz')))
  rtn <- left_join(heritability_df %>%
                     dplyr::filter(target %in% df$gene),
                   wt_downstream_meta %>% dplyr::select(c('downstream_gene_name', 
                                                          'control_expr_dLL_donor_pval' = 'pval', 'control_expr_dLL_donor_pval_adj' = 'pval_adj'))
  )
})
heritability_df <- as.data.frame(bind_rows(wt_expr_df))


## do a full list
fwrite(heritability_df, heritability_fnm_out, sep = '\t', compress = 'gzip')


## ---- Calculate Adjusted P-values
## do a shortened list

heritability_df <- fread(heritability_fnm_out)

pairs2compute <- heritability_df %>%
  dplyr::filter(!(is.na(control_expr_dLL_donor_pval)) & !(is.na(lfc_dLL_donor_pval)) & !(is.na(post_kd_expr_dLL_donor_pval)) & 
                  !(is.infinite(control_expr_dLL_donor_pval)) & !(is.infinite(lfc_dLL_donor_pval)) & !(is.infinite(post_kd_expr_dLL_donor_pval)) &
                  bulk_range_lfc > sig_abs_lfc_thresh &
                  !(downstream_gene_name %in% crispr_marker_sequences) &
                  potential_cross_map == 'n' &
                  potential_off_target == 'n') 
##  
pairs2compute$lfc_dLL_donor_pval_adj <- p.adjust(pairs2compute$lfc_dLL_donor_pval, method = 'BH')
pairs2compute$post_kd_expr_dLL_donor_pval_adj <- p.adjust(pairs2compute$post_kd_expr_dLL_donor_pval, method = 'BH')
pairs2compute <- pairs2compute %>%
  dplyr::select(c("target", 'downstream_gene_name', 
                  "bulk_mean_control_expr", "bulk_control_expr_var", ## target meta
                  "bulk_mean_lfc", "bulk_lfc_var", ## bulk lfc meta
                  'control_expr_dLL_donor', 'control_expr_dLL_donor_pval', "control_expr_dLL_donor_pval_adj",## control expression
                  'lfc_dLL_donor', 'lfc_dLL_donor_pval', "lfc_dLL_donor_pval_adj", ## lfc
                  'post_kd_expr_dLL_donor', 'post_kd_expr_dLL_donor_pval', "post_kd_expr_dLL_donor_pval_adj", ## post knockdown,
                  contains('varExpl')
  )) 
fwrite(pairs2compute, paste0(outdir, '/', date, "_heritability_with_pval_adj.tsv.gz"), sep = '\t', compress = 'gzip')

significant_pairs <- pairs2compute %>%
  dplyr::filter(lfc_dLL_donor_pval_adj < sig_pval_thresh)
fwrite(significant_pairs, paste0(outdir, '/', date, "_significant_pairs.tsv.gz"), sep = '\t', compress = 'gzip')


## ---- Check P-values

## plot the significant hits
date <- '2023-06-03'
#pairs2plot <- fread(paste0(outdir, '/', date, "_significant_pairs.tsv.gz"))
pairs2plot <- fread('/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Pica/15c_calc_var_expl_by_target/15c_05_calc_var_expl_by_target_check_pvals/2024-06-03_significant_pairs.tsv.gz')
#pairs2plot <- fread('/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Pica/15c_calc_var_expl_by_target/2024-05-17/2024-05-17_significant_pairs.tsv.gz')

i <- 1
#tg <- pairs2plot$target[i]; dg <- pairs2plot$downstream_gene_name[i]
dgs <- sort(unique(pairs2plot$downstream_gene_name))

pdf(paste0(plotsdir, '/', date, '_sig_pairs.pdf'), width= 12, height = 3)
lapply(dgs, FUN = function(dg){
  print(dg)
  expr_df_fnm <- list.files(paste0(OutFolder, "6b_calc_lfcs_transcriptome_wide_by_gene/by_downstream/"),
                            pattern = paste0('-', gsub( dg, pattern = '/|:', replacement= '-'), '-Gene-Expression'), full.names = T)
  expr_df <- fread(expr_df_fnm)
  if (exists('donors2include')){
    d <- unlist(strsplit(donors2include, "_"))
    expr_df <- expr_df %>% dplyr::filter(unlist(lapply(strsplit(cell_line, "_"), "[[", 1)) %in% d)
  }
    
  tgs <- pairs2plot %>% dplyr::filter(downstream_gene_name == dg) %>% .$target %>% sort() %>% unique()
  
  lapply(tgs, FUN = function(tg){
    
    i <- which(pairs2plot$target == tg & pairs2plot$downstream_gene_name == dg)
    expr_df_subset <- expr_df %>%
      dplyr::filter(target_gene_name %in% c(tg, NonTargetGeneName, "unassigned"))
    control_expr <- expr_df_subset %>%
      dplyr::filter(target_gene_name %in% c(NonTargetGeneName, "unassigned")) %>%
      dplyr::select(c("cell_line", 'y' = 'regressed_expr'))%>%
      group_by(cell_line) %>%
      summarize(y = mean(y),
                status = 'control_expr')
    knockdown_lfc <- expr_df_subset %>%
      dplyr::filter(target_gene_name ==tg) %>%
      dplyr::select(c("cell_line", "y" = 'regressed_lfc'))%>%
      group_by(cell_line) %>%
      summarize(y = mean(y),
                status = 'lfc')
    post_kd_expr <- expr_df_subset %>%
      dplyr::filter(target_gene_name == tg) %>%
      dplyr::select(c("cell_line", "y" = 'regressed_expr')) %>%
      group_by(cell_line) %>%
      summarize(y = mean(y),
                status = 'post_kd_expr')
    df4plot <- as.data.frame(bind_rows(control_expr, knockdown_lfc, post_kd_expr)) %>%
      mutate(donor = unlist(lapply(strsplit(cell_line, "_"), "[[", 1)))
    
    p_lfc <- ggplot(df4plot, aes(x = cell_line, y = y, fill = donor)) + 
      facet_wrap(~status, ncol = 3) + 
      ggtitle(paste0(dg, ' expression in ', tg, ' knockdowns')) + 
      geom_bar(stat = 'identity')+theme_bw()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    
    
    df4plot <- pairs2plot[i,] %>% dplyr::select(contains('pval_adj')) %>%
      reshape2::melt() %>%
      mutate(variable = gsub(variable, pattern = "_dLL.*", replacement = ''))
    p_h <- ggplot(df4plot, aes(x= variable, y= value)) + 
      ggtitle("heritability (pval)") + 
      geom_bar(stat = 'identity') + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    
    p <- ggarrange(p_lfc, p_h, ncol = 2, widths = c(2.5, 1))
    print(p)
  })
})
dev.off()

## gain in heritability 
heritability_gain <- pairs2plot %>%
  dplyr::select(c("target", 'downstream_gene_name', contains('pval_adj'))) %>%
  dplyr::filter(control_expr_dLL_donor_pval_adj > 0.5 & lfc_dLL_donor_pval_adj < sig_pval_thresh & post_kd_expr_dLL_donor_pval_adj < sig_pval_thresh)
dgs <- sort(unique(heritability_gain$downstream_gene_name))

pdf(paste0(plotsdir, '/', date, '_gained_heritability.pdf'), width= 12, height = 3)
plotlist <- lapply(dgs, FUN = function(dg){
  print(dg)
  expr_df_fnm <- list.files(paste0(OutFolder, "6b_calc_lfcs_transcriptome_wide_by_gene/by_downstream/"),
                            pattern = paste0('-', gsub( dg, pattern = '/|:', replacement= '-'), '-Gene-Expression'), full.names = T)
  expr_df <- fread(expr_df_fnm)
  if (exists('donors2include')){
    d <- unlist(strsplit(donors2include, "_"))
    expr_df <- expr_df %>% dplyr::filter(unlist(lapply(strsplit(cell_line, "_"), "[[", 1)) %in% d)
  }
  
  tgs <- heritability_gain %>% dplyr::filter(downstream_gene_name == dg) %>% .$target %>% sort() %>% unique()
  
  lapply(tgs, FUN = function(tg){
    
    i <- which(heritability_gain$target == tg & heritability_gain$downstream_gene_name == dg)
    expr_df_subset <- expr_df %>%
      dplyr::filter(target_gene_name %in% c(tg, NonTargetGeneName, "unassigned"))
    control_expr <- expr_df_subset %>%
      dplyr::filter(target_gene_name %in% c(NonTargetGeneName, "unassigned")) %>%
      dplyr::select(c("cell_line", 'y' = 'regressed_expr'))%>%
      group_by(cell_line) %>%
      summarize(y = mean(y),
                status = 'control_expr')
    knockdown_lfc <- expr_df_subset %>%
      dplyr::filter(target_gene_name ==tg) %>%
      dplyr::select(c("cell_line", "y" = 'regressed_lfc'))%>%
      group_by(cell_line) %>%
      summarize(y = mean(y),
                status = 'lfc')
    post_kd_expr <- expr_df_subset %>%
      dplyr::filter(target_gene_name == tg) %>%
      dplyr::select(c("cell_line", "y" = 'regressed_expr')) %>%
      group_by(cell_line) %>%
      summarize(y = mean(y),
                status = 'post_kd_expr')
    df4plot <- as.data.frame(bind_rows(control_expr, knockdown_lfc, post_kd_expr)) %>%
      mutate(donor = unlist(lapply(strsplit(cell_line, "_"), "[[", 1)))
    
    p_lfc <- ggplot(df4plot, aes(x = cell_line, y = y, fill = donor)) + 
      facet_wrap(~status, ncol = 3) + 
      ggtitle(paste0(dg, ' expression in ', tg, ' knockdowns')) + 
      geom_bar(stat = 'identity')+theme_bw()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    
    
    df4plot <- heritability_gain[i,] %>% dplyr::select(contains('pval_adj')) %>%
      reshape2::melt() %>%
      mutate(variable = gsub(variable, pattern = "_dLL.*", replacement = ''))
    p_h <- ggplot(df4plot, aes(x= variable, y= value)) + 
      ggtitle("heritability (pval)") + 
      geom_bar(stat = 'identity') + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    
    p <- ggarrange(p_lfc, p_h, ncol = 2, widths = c(2.5, 1))
    #print(p)
    return(p)
  })
})
plotlist
dev.off()

## ---- SessionInfo

sessionInfo()



