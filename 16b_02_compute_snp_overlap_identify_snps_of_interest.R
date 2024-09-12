

suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))
suppressMessages(library(glmnet))
suppressMessages(library(doParallel))
suppressMessages(library(ggpubr))
suppressMessages(library(tidyverse))
suppressMessages(library(doParallel))
suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(vcfR))


section_name <- "16b_compute_snp_overlap"
subsection_name <- "16b_02_compute_snp_overlap_identify_snps_of_interest"
analysis_name <- "pipeline"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
ExperimentName <- "Pica"
date <- "2024-06-03" #(magpie)
donors2include='eipl_iudw_jejf_kolf_paab'


args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  ##needed arguments
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--experiment_name"){ ExperimentName <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
  if (arg == "--donors2include"){ donors2include <- args[[ind + 1]] }
}

if (length(donors2include) == 1){donors2include <- unlist(strsplit(donors2include, "_"))}

# set relevent i/o paths
source(file.path(HomeFolder, paste0("scripts/io/", ExperimentName, "_io.R")))
source(file.path(HomeFolder, "scripts/io/Magpie_Utils.R"))
outdir <- file.path(OutFolder, "", section_name, '/', subsection_name)
plotsdir <- file.path(HTMLFolder, section_name)
if(!dir.exists(outdir)) dir.create(outdir)
if(!dir.exists(plotsdir)) dir.create(plotsdir)
print(section_name)

## ---- LoadData

all_var_expl <- fread("/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Pica/15c_calc_var_expl_by_target/15c_02_calc_var_expl_by_target_variance_decomposition/2024-06-03_var_expl_combined_no_pval.tsv.gz")


pica_heritability_df <- fread('/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Pica/15c_calc_var_expl_by_target/15c_05_calc_var_expl_by_target_check_pvals/2024-06-03_heritability_with_pval_adj.tsv.gz')
pica_heritability_df <- pica_heritability_df %>%
    dplyr::filter(control_expr_dLL_donor_pval_adj > 0.5 & 
                    lfc_dLL_donor_pval_adj < 0.1 & 
                    post_kd_expr_dLL_donor_pval_adj < 0.1)

cis_eQTL_list <- fread("/lustre/scratch123/hgi/projects/crispri_scrnaseq/resources/eQTLs/mjb_eQTL_list_downloaded_2023-02-27.csv", skip = 1, header = T)
gene_names <- ensembldb::select(EnsDb.Hsapiens.v86, keys= cis_eQTL_list$feature_id, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
cis_eQTL_list$gene_name <- gene_names$SYMBOL[match(cis_eQTL_list$feature_id, gene_names$GENEID)]

## all downstream genes
cis_eQTL_list <- cis_eQTL_list %>%
  dplyr::filter(gene_name %in% all_var_expl$downstream_gene_name)
fwrite(cis_eQTL_list, paste0(outdir, '/', date, '_all_downstream_gene_cis_eQTLs.tsv'), sep = '\t')


## all gain in heritability
#cis_eQTL_list <- cis_eQTL_list %>%
#  dplyr::filter(gene_name %in% pica_heritability_df$downstream_gene_name)
#fwrite(cis_eQTL_list, paste0(outdir, '/', date, '_gain_in_heritability_downstream_gene_cis_eQTLs.tsv'), sep = '\t')


donor_gts <- fread('/lustre/scratch123/hgi/projects/crispri_scrnaseq/resources/hipsci_genotypes/pica_genotypes_hg19/all_lines_gts.tsv')
hipsci_line_ids <- unlist(lapply(strsplit(gsub(list.files(paste0(ResourcesFolder, 'hipsci_genotypes/pica_genotypes_hg19'), pattern = 'genotypes.vcf.gz.tbi'), 
                        pattern = "[.]we.*", replacement = ''), "-"), "[[", 2))
names(donor_gts) <- c("CHROM", "POS", "REF", "ALT", hipsci_line_ids , 'na')

snps_of_interest <- data.frame(
  CHROM = unlist(lapply(strsplit(cis_eQTL_list$snp_id, "_"), '[[', 1)),
  POS = as.integer(unlist(lapply(strsplit(cis_eQTL_list$snp_id, "_"), '[[', 2))),
  REF = unlist(lapply(strsplit(cis_eQTL_list$snp_id, "_"), '[[', 3)),
  ALT = unlist(lapply(strsplit(cis_eQTL_list$snp_id, "_"), '[[', 4)),
  gene = cis_eQTL_list$gene_name,
  maf = cis_eQTL_list$maf
)
donor_eQTL_status <- left_join(snps_of_interest, donor_gts)
fwrite(donor_eQTL_status, paste0(outdir, '/', date, '_downstream_cis_eQTL_gts.tsv'), sep = '\t')






