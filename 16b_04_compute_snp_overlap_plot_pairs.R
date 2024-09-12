suppressMessages(library(Matrix))
suppressMessages(library(Seurat))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(tidyverse))
suppressMessages(library(parallel))
suppressMessages(library(doParallel))
suppressMessages(library(lme4))
suppressMessages(library(lattice))
suppressMessages(library(variancePartition))
suppressMessages(library(gtools))
suppressMessages(library(lmtest))


# set relevent i/o paths
section_name <- "16b_compute_snp_overlap"
subsection_name <- "16b_04_compute_snp_overlap_plot_pairs"

# set relevent i/o paths
date <- "2024-06-03"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
ExperimentName <- "Pica"
donors2include <- 'eipl_iudw_jejf_kolf_paab'

# options for testing
io_path <- paste0("scripts/io/", ExperimentName, "_io.R")
utils_path <- "scripts/io/Magpie_Utils.R"

args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  ##needed arguments
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  ## not needed
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
  if (arg == "--subsection_name"){ subsection_name <- args[[ind + 1]] }
  if (arg == "--utils_path"){ utils_path <- args[[ind + 1]] }
  if (arg == "--io_path"){ io_path <- args[[ind + 1]] }
}

print(paste0(section_name))
print(paste("date:", date))
print(paste("home folder:", HomeFolder))


setwd(HomeFolder)
source(io_path)
source(utils_path)

GuideMetadata <- fread(GuideMetadataPath)
LineMetadata <- fread(LineMetadataPath)

outdir <- file.path(file.path(OutFolder, section_name, subsection_name))
if(!dir.exists(outdir)) dir.create(outdir, recursive = T)
plotsdir <- file.path(HTMLFolder, "pipeline", section_name) 
if(!dir.exists(plotsdir)) dir.create(plotsdir, recursive = T)
if(length(donors2include) == 1){
  donors2include <- unlist(strsplit(donors2include, "_"))
}

## ---- load data

pica_heritability_df <- fread(paste0(OutFolder, '15c_calc_var_expl_by_target/15c_05_calc_var_expl_by_target_check_pvals/', date, '_heritability_with_pval_adj.tsv.gz'))

pica_var_expl <- fread(paste0(OutFolder, '15c_calc_var_expl_by_target/2024-05-17/2024-05-17_heritability_all_no_pval_adj.tsv.gz'))

donor_eQTL_status <- fread(paste0(OutFolder, section_name, '/', "16b_02_compute_snp_overlap_identify_snps_of_interest", '/', date, '_downstream_cis_eQTL_gts.tsv'))
cis_eQTL_list <- fread(paste0(OutFolder, section_name, '/', "16b_02_compute_snp_overlap_identify_snps_of_interest", '/', date, '_all_downstream_gene_cis_eQTLs.tsv'))
cis_eQTL_list$beta4comparison <- ifelse(unlist(lapply(strsplit(cis_eQTL_list$snp_id, "_"), "[[", 4)) == cis_eQTL_list$assessed_allele, cis_eQTL_list$beta, 
                                        ifelse(unlist(lapply(strsplit(cis_eQTL_list$snp_id, "_"), "[[", 3)) == cis_eQTL_list$assessed_allele, -cis_eQTL_list$beta, NA))
var_expl_snp_fnms <- list.files(paste0(OutFolder, section_name, '/16b_03_compute_snp_overlap_run_snp_model/by_gene/'), pattern = 'gt_effect', full.names = T)
var_expl_snp <- lapply(var_expl_snp_fnms, fread) %>% bind_rows() %>% as.data.frame()
var_expl_snp$gt_pval_adj <- p.adjust(var_expl_snp$gt_pval, method = "BH")
var_expl_snp <- left_join(var_expl_snp,
                          pica_heritability_df %>%
                            dplyr::select(c("target", 'snp_gene' = 'downstream_gene_name',
                                            'control_expr_dLL_donor_pval_adj',
                                            'lfc_dLL_donor_pval_adj',
                                            'post_kd_expr_dLL_donor_pval_adj')))
var_expl_snp$status <- ifelse(var_expl_snp$lfc_dLL_donor_pval_adj > sig_pval_thresh, 'Insignificant',
                              ifelse(var_expl_snp$control_expr_dLL_donor_pval_adj > 0.5 & var_expl_snp$post_kd_expr_dLL_donor_pval_adj < sig_pval_thresh, 'Gain in heritability',
                                     ifelse(var_expl_snp$control_expr_dLL_donor_pval_adj < sig_pval_thresh & var_expl_snp$post_kd_expr_dLL_donor_pval_adj < sig_pval_thresh, 'Heritability maintained',
                                            ifelse(var_expl_snp$control_expr_dLL_donor_pval_adj < sig_pval_thresh & var_expl_snp$post_kd_expr_dLL_donor_pval_adj > 0.5,'Loss of heritability', 'Undefined' ))))
var_expl_snp <- left_join(var_expl_snp %>%
                            mutate(snp_id = paste(snp_chrom, snp_pos, snp_ref, snp_alt, sep= '_')),
                          cis_eQTL_list %>%
                            dplyr::select('snp_id', 'beta' = 'beta4comparison', 'p_value'))




fwrite(var_expl_snp, paste0(outdir, '/', date, '_gt_effect.tsv'), sep='\t')


## ---- Plot

var_expl_snp <- fread(paste0(outdir, '/', date, '_gt_effect.tsv'))
#pica_lfcs_by_line <- fread('/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Pica/6f_calc_lfcs_transcriptome_wide_by_gene_per_line/2024-05-17_all_lines_with_adj_pval.tsv.gz')
#pica_lfcs_by_line_split <- split.data.frame(pica_lfcs_by_line)
gain_in_h <- pica_heritability_df %>%
  dplyr::filter(varExpl)

tg <- 'EP400'; dg <- 'SPP1'
snp_id2match <- "4_88823621_C_T"
pairs2plot <- var_expl_snp %>%
  dplyr::filter(gt_pval_adj < 0.1) %>%
  arrange(snp_gene)

cis_eQTL_info <- cis_eQTL_list %>%
  dplyr::filter(gene_name == dg &  snp_id == snp_id2match)



snp_jitters <- lapply(1:dim(pairs2plot)[1], FUN = function(i){
  tg <- pairs2plot$target[i]; dg <- pairs2plot$snp_gene[i]
  snp_id <- paste0(pairs2plot$snp_chrom[i], "_", pairs2plot$snp_pos[i], "_", pairs2plot$snp_ref[i], "_", pairs2plot$snp_alt[i])
  snp_vec <- donor_eQTL_status %>%
    dplyr::filter(gene == dg & paste0(CHROM, "_", POS, "_", REF, "_", ALT) == snp_id)%>% dplyr::select(contains("_"))
  
  snp_id2 <- snp_id
  beta <- cis_eQTL_list %>%
    dplyr::filter(snp_id == snp_id2) %>% .$beta
  snp_vec <- data.frame(
    cell_line = names(snp_vec),
    donor = unlist(lapply(strsplit(names(snp_vec), "_"), "[[", 1)),
    gt = gsub(gsub(snp_vec, pattern = '1[|]0', replacement = '0|1'), pattern = '/', replacement = '|')
  ) %>%
    dplyr::filter(gt != '.|.' & donor %in% donors2include) %>%
    group_by(donor) %>%
    summarize(gt = unique(gt)) 
  regressed_lfc_fnm <- list.files(paste0(OutFolder, '/6b_calc_lfcs_transcriptome_wide_by_gene/by_downstream/'), pattern = paste0('-', dg, '-Gene-Expression'), full.names = T)
  regressed_lfc <- fread(regressed_lfc_fnm)
  
  df4plot <-regressed_lfc %>%
    dplyr::filter(target_gene_name == tg & donor %in% donors2include) %>%
    dplyr::select(c('V1', 'guide' = 'feature_call', 'lfc' = 'regressed_lfc', 'donor', 'cell_line')) %>%
    mutate(gt = snp_vec$gt[match(donor, snp_vec$donor)])
  
  
  p1 <- ggplot(df4plot, aes(x = cell_line, y = lfc, col = gt)) + 
    ggtitle(paste0(dg, ' expression in ', tg, ' knockdowns')) + 
    geom_jitter() + 
    stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.5, color="red") + 
    annotate('text', label = paste0('beta = ', beta), x = -Inf, y = Inf, vjust = 2, hjust = -0.5) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  #p1
  
  df4bar_plot_control <- regressed_lfc %>%
    dplyr::filter(target_gene_name == NonTargetGeneName & donor %in% donors2include) %>%
    group_by(cell_line) %>%
    summarize(
      donor = unique(donor),
      status = 'control',
      expr = mean(regressed_expr)
    )
  df4bar_plot_kd <- regressed_lfc %>%
    dplyr::filter(target_gene_name == tg & donor %in% donors2include) %>%
    group_by(cell_line) %>%
    summarize(
      donor = unique(donor),
      lfc = mean(regressed_lfc),
      post_kd_expr = mean(regressed_expr)
    ) %>%
    reshape2::melt(id = c("cell_line", "donor"))
  names(df4bar_plot_kd) <- names(df4bar_plot_control)
  df4bar_plot <- as.data.frame(bind_rows(df4bar_plot_control, df4bar_plot_kd))
  
  p2 <- ggplot(df4bar_plot, aes(x = cell_line, y = expr, fill = donor))+ 
    facet_wrap(~status, ncol = 1) + 
    geom_bar(stat = 'identity') + theme_bw() +
    theme(legend.position = 'none', 
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  p <- ggarrange(p1, p2, ncol = 1, heights = c(1,3))
  return(p)
})

pdf('../../../scratch/2024-07-30_jitters2.pdf', width = 5, height = 3)
snp_jitters
dev.off()



## ---- 


rmarkdown::render(input = paste0(CodeFolder, 'Magpie/pipeline/16b_04_compute_snp_overlap_plot_pairs.Rmd'),
                  output_file = paste0(plotsdir, '/', date, '.html'))







