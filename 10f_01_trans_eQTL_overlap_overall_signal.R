
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))
suppressMessages(library(doParallel))
suppressMessages(library(ggpubr))
suppressMessages(library(tidyverse))

section_name <- "10f_trans_eQTL_overlap"
subsection_name <- "10f_01_trans_eQTL_overlap_overall_signal"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
ExperimentName <- "Magpie"
magpie_version <- "2022-08-15" 

args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  ##needed arguments
  if (arg == "--date"){ magpie_version <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  ## not needed
  if (arg == "--utils_path"){ utils_path <- args[[ind + 1]] }
  if (arg == "--io_path"){ io_path <- args[[ind + 1]] }
}

# set relevent i/o paths
source(paste0(HomeFolder, "scripts/io/", ExperimentName, "_io.R"))
source(paste0(HomeFolder, "scripts/io/Magpie_Utils.R"))
outdir <- file.path(OutFolder, section_name, subsection_name)
plotsdir <- file.path(HTMLFolder, 'pipeline/', section_name)
if(!dir.exists(outdir)) dir.create(outdir, recursive = T)
if(!dir.exists(plotsdir)) dir.create(plotsdir, recursive = T)
print(section_name)

## ---- LoadData

eQTL_list <- fread(paste0(OutFolder, "Preprocess_External_Datasets/eQTLs/", date, "_cis_eQTL_list_with_cis_beta.tsv.gz"))
magpie_lfcs <- fread(paste0(OutFolder, '6b_calc_lfcs_transcriptome_wide_by_gene/', date, '_combined_with_adj_pval.tsv.gz'))

sig_trans_eQTLs <- fread(paste0(OutFolder, 'Preprocess_External_Datasets/eQTLs/eQTL_interaction_pairs_magpie-', date, '-version.tsv'))
magpie_target_meta <- fread(paste0(OutFolder, '7b_target_summary/', date, '_target_meta_data.csv'))

## ---- Preprocess Data


df <- left_join(eQTL_list %>%
                  dplyr::select(c("trans_gene_name", "cis_gene_name", 
                                  'snp_beta', 'snp_p_value', 'snp_id',
                                  'trans_eQTL_p_value', 'trans_eQTL_beta')),
                magpie_lfcs %>%
                  dplyr::select(c("trans_gene_name" = 'downstream_gene_name', 'cis_gene_name' = 'target', 'crispr_lfc' = 'lfc', 'crispr_pval' = 'pval_lm', 'crispr_pval_adj' = 'pval_adj'))
) %>%
  dplyr::filter(!is.na(crispr_lfc)) 
df <- left_join(df,
                magpie_target_meta %>%
                  dplyr::select('cis_gene_name' = 'gene',
                                "crispr_lfc" = 'target_lfc',
                                'crispr_pval_adj' = 'pval_adj',
                                'crispr_n_trans' = 'n_downstream_excl_target'))
fwrite(df, paste0(outdir, '/trans_eQTL_overlap.tsv.gz' ), sep = '\t', compress = 'gzip')

## add in on-target lfc
## also add in # of degs

## ---- Compute Correlation

df <- fread(paste0(outdir, '/trans_eQTL_overlap.tsv.gz' ))
snps2plot <- split.data.frame(df, f = df$snp_id)

## correlation per snp
trans_effect_info <- lapply(snps2plot, FUN = function(df2plot){
  trans_genes <- sig_trans_eQTLs %>%
    dplyr::filter(cis_gene == df2plot$cis_gene_name[1]) %>% .$trans_gene
  p_effect <- ggplot(df2plot, aes(x = trans_eQTL_beta, y = crispr_lfc,
                                  col = ifelse(crispr_pval_adj < 0.1, 'y', 'n'),
                                  label = ifelse(trans_gene_name %in% trans_genes, trans_gene_name, ''))) + 
    ggrepel::geom_text_repel() + 
    scale_color_manual('', values = c('y' = 'red', 'n' = 'gray')) + 
    ggtitle(paste0(df2plot$cis_gene_name, " - ", df2plot$snp_id)) + 
    geom_point()+ theme_bw() + theme(legend.position = 'none') 
  p_qq <- ggplot(df2plot, aes(x = trans_eQTL_p_value, y = crispr_pval)) + 
    ggtitle(paste0(df2plot$cis_gene_name, " - ", df2plot$snp_id)) + 
    geom_point()+ theme_bw() + theme(legend.position = 'none') 
  rtn <- list()
  rtn$plot <- ggarrange(p_effect, p_qq)
  rtn$cor_effect <- cor(df2plot$crispr_lfc, df2plot$trans_eQTL_beta)
  rtn$cor_pval <- cor(df2plot$crispr_pval, df2plot$trans_eQTL_p_value)
  return(rtn)
})

df2save <- lapply(1:length(trans_effect_info), FUN = function(i){
  #print(i)
  r <- trans_effect_info[[i]]
  snp_info <- snps2plot[[i]]
  rtn <- data.frame(
    cis_gene = snp_info$cis_gene_name[1],
    snp_beta = snp_info$snp_beta[1],
    snp_p_value = snp_info$snp_p_value[1],
    snp_id = snp_info$snp_id[1],
    beta_cor = r$cor_effect,
    pval_cor = r$cor_pval
  )
  return(rtn)
}) %>% bind_rows() %>% as.data.frame()
fwrite(df2save, paste0(outdir, '/trans_effect_info.tsv'), sep = '\t')


## ---- MakePlots

df <- fread(paste0(outdir, '/trans_eQTL_overlap.tsv.gz' ))
trans_effect_info <- fread(paste0(outdir, '/trans_effect_info.tsv')) %>%
  mutate(has_trans_eQTL = ifelse(cis_gene %in% sig_trans_eQTLs$cis_gene, 'y', 'n'))

pdf(paste0(plotsdir, '/', date, '_correlation_per_snp.pdf'), width = 8, height = 4)
lapply(trans_effect_info, FUN = function(r){r$plot})
dev.off()

## one-by-one

ggplot(df2save, aes(x = abs(snp_beta), y = abs(beta_cor))) + 
  geom_point()

ggplot(df2save, aes(x = -log10(snp_p_value), y = abs(beta_cor))) + 
  geom_point()

ggplot(df2plot, aes(x = ))
