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
subsection_name <- "16b_03_compute_snp_overlap_run_snp_model"

# set relevent i/o paths
date <- "2024-06-03"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
ExperimentName <- "Pica"

# options for testing
io_path <- paste0("scripts/io/", ExperimentName, "_io.R")
utils_path <- "scripts/io/Magpie_Utils.R"
control_tag <- c("unassigned", "NonTarget") # use all unassigned cells as control
#downstream_gene <- 'PRKD3'
#downstream_gene_fnm <- '/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Pica/6b_calc_lfcs_transcriptome_wide_by_gene/by_downstream/2024-05-17_ENSG00000104722-NEFM-Gene-Expression.tsv.gz'
downstream_gene_fnm <- '/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Pica/6b_calc_lfcs_transcriptome_wide_by_gene/by_downstream/2024-05-17_ENSG00000185885-IFITM1-Gene-Expression.tsv.gz'
donors2include <- 'eipl_iudw_jejf_kolf_paab'

args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  ##needed arguments
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--downstream_gene_fnm"){ downstream_gene_fnm <- args[[ind + 1]] }
  if (arg == "--donors2include"){ donors2include <- args[[ind + 1]] }
  
  ## not needed
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
  if (arg == "--subsection_name"){ subsection_name <- args[[ind + 1]] }
  if (arg == "--utils_path"){ utils_path <- args[[ind + 1]] }
  if (arg == "--io_path"){ io_path <- args[[ind + 1]] }
}

print(paste0(section_name))
print(paste("date:", date))
print(paste("home folder:", HomeFolder))
if(exists("donors2include")){print('including donors: '); print(donors2include)}
if(length(donors2include) == 1){
  donors2include <- unlist(strsplit(donors2include, "_"))
}


setwd(HomeFolder)
source(io_path)
source(utils_path)

GuideMetadata <- fread(GuideMetadataPath)
LineMetadata <- fread(LineMetadataPath)

outdir <- file.path(file.path(OutFolder, section_name, subsection_name, 'by_gene'))
if(!dir.exists(outdir)) dir.create(outdir, recursive = T)
plotsdir <- file.path(HTMLFolder, "pipeline", section_name) 
if(!dir.exists(plotsdir)) dir.create(plotsdir, recursive = T)

downstream_gene <- gsub(gsub(gsub(downstream_gene_fnm, pattern = ".*by_downstream/", replacement = ''),
                        pattern = "-Gene-Expression.tsv.gz", replacement = ''),
                        pattern = '.*_', replacement = '')
downstream_gene <- unlist(lapply(strsplit(downstream_gene, '-'), FUN = function(x){paste(x[-1], collapse = '-')}))
print(paste0("computing cis eqtl effects for ", downstream_gene))

## ---- load data


pica_heritability_df <- fread(paste0(OutFolder, '15c_calc_var_expl_by_target/15c_05_calc_var_expl_by_target_check_pvals/', date, '_heritability_with_pval_adj.tsv.gz'))
#pica_heritability_df <- pica_heritability_df %>%
#  dplyr::filter(#control_expr_dLL_donor_pval_adj > 0.5 & 
#                  lfc_dLL_donor_pval_adj < 0.1)# & 
                  #post_kd_expr_dLL_donor_pval_adj < 0.1)

regressed_lfc_fnm <- list.files(paste0(OutFolder, '/6b_calc_lfcs_transcriptome_wide_by_gene/by_downstream/'), pattern = paste0('-', downstream_gene, '-Gene-Expression'), full.names = T)
regressed_lfc <- fread(regressed_lfc_fnm)

donor_eQTL_status <- fread(paste0(OutFolder, section_name, '/', "16b_02_compute_snp_overlap_identify_snps_of_interest", '/', date, '_downstream_cis_eQTL_gts.tsv'))
snps2test <- donor_eQTL_status %>%
  dplyr::filter(gene == downstream_gene)

## ---- Calculate Heritability

## decide which models need to be fun


targets2test <- pica_heritability_df %>%
  dplyr::filter(downstream_gene_name == downstream_gene) %>% 
  .$target
if (length(targets2test) >= 1 & dim(snps2test)[1] >= 1){
  res <- lapply(targets2test, FUN = function(tg){
    print(tg)
    on_target_expr_fnm <- list.files(paste0(OutFolder, '/6b_calc_lfcs_transcriptome_wide_by_gene/by_downstream/'), pattern = paste0('-', tg, '-Gene-Expression'), full.names = T)
    print(on_target_expr_fnm)
    if (length(on_target_expr_fnm) == 1){
      on_target_expr <- fread(on_target_expr_fnm)
      on_target_expr <- on_target_expr %>%
        dplyr::filter(target_gene_name == tg)
      
      #gt_ind <- 4
      varExpl_snps <- lapply(1:dim(snps2test)[1], FUN = function(gt_ind){
        print(gt_ind)
        snp_vec <- snps2test[gt_ind,] %>% dplyr::select(contains("_"))
        snp_vec <- data.frame(
          cell_line = names(snp_vec),
          donor = unlist(lapply(strsplit(names(snp_vec), "_"), "[[", 1)),
          gt = gsub(gsub(snp_vec, pattern = '1[|]0', replacement = '0|1'), pattern = '/', replacement = '|')
        ) %>%
          dplyr::filter(gt != '.|.' & donor %in% donors2include) %>%
          group_by(donor) %>%
          summarize(gt = unique(gt))
        
        if (length(unique(snp_vec$gt)) > 1 & all(table(snp_vec$donor)==1)){
          df4lmm <- regressed_lfc %>%
            dplyr::filter(target_gene_name == tg & donor %in% donors2include) %>%
            dplyr::select(c('V1', 'guide' = 'feature_call', 'lfc' = 'regressed_lfc', 'donor', 'cell_line')) %>%
            mutate(on_target_expr = on_target_expr$regressed_expr[match(V1, on_target_expr$V1)]) %>%
            mutate(gt = snp_vec$gt[match(donor, snp_vec$donor)])
          
          ## model 1: variance explained
          fit <- suppressMessages(suppressWarnings(
            tryCatch(expr = lme4::lmer(lfc ~ 1 + on_target_expr + (1|guide) + (1|cell_line) + (1|gt) + (1|gt:donor),
                                       data = df4lmm), 
                     warning = function(x){
                       return(NULL)
                     }, 
                     error = function(e) {
                       return(NULL)
                     })
          ))
          
          
          if (!(is.null(fit))){
            
            ## variance explained model
            varExplained <- tryCatch(
              calcVarPart(fit, showWarnings = F),
              warning = function(x){
                return(rep(NULL, length(c(fixed_effects, random_effects))))
              }, 
              error = function(e) {
                return(rep(NULL, length(c(fixed_effects, random_effects))))
              }
            )
            
            
            
          } else {varExplained <- NULL}
          
          
          ## model 2: effect of # of alt alleles
          
          df4lm <-regressed_lfc %>%
            dplyr::filter(target_gene_name == tg & donor %in% donors2include) %>%
            dplyr::select(c('V1', 'guide' = 'feature_call', 'lfc' = 'regressed_lfc', 'donor', 'cell_line')) %>%
            mutate(gt = snp_vec$gt[match(donor, snp_vec$donor)]) %>%
            mutate(n_alt = unlist(lapply(strsplit(gt, "[|]"), FUN = function(x){length(which(as.numeric(x) == 1))})))
          
          fit <- lm(df4lm, formula= lfc ~ n_alt)
          
          rtn <- list()
          if (!(is.null(fit))){
            rtn$gt_effect <- snps2test[gt_ind,] %>%
              dplyr::select(c("snp_chrom" = "CHROM",
                              'snp_pos' = 'POS',
                              'snp_ref'= 'REF',
                              'snp_alt' = 'ALT',
                              'snp_gene' = 'gene')) %>%
              mutate(target = tg,
                     gt_coef = summary(fit)$coefficients['n_alt', 1],
                     gt_pval = summary(fit)$coefficients['n_alt', 4])
          } else {rtn$gt_effect <- NULL}
          
          if (!(is.null(varExplained))){
            rtn$varExpl <- data.frame(
              snp_chrom = snps2test$CHROM[gt_ind],
              snp_pos = snps2test$POS[gt_ind],
              snp_ref = snps2test$REF[gt_ind],
              snp_alt = snps2test$ALT[gt_ind],
              snp_gene = snps2test$gene[gt_ind],
              target = tg, 
              t(varExplained)
            )
          } else{rtn$varExpl <- NULL}
          
        } else {rtn <- list(varExpl = NULL,
                            gt_effect = NULL)}
        
        return(rtn)
        
      })
      
    } else {varExpl_snps <- NULL}
    
    return(varExpl_snps)
  })
  
  
  varExpl_res <- lapply(res, FUN = function(x){lapply(x, FUN = function(y){y$varExpl}) %>% bind_rows() %>% as.data.frame()}) %>% bind_rows() %>% as.data.frame()
  if (dim(varExpl_res)[1] > 0){
    fwrite(varExpl_res, paste0(outdir, '/', date, "_", downstream_gene, "_var_expl.tsv.gz"), sep = '\t')
  }
  gt_effect_res <-  lapply(res, FUN = function(x){lapply(x, FUN = function(y){y$gt_effect}) %>% bind_rows() %>% as.data.frame()}) %>% bind_rows() %>% as.data.frame()
  if (dim(gt_effect_res)[1] > 0){
    fwrite(gt_effect_res, paste0(outdir, '/', date, "_", downstream_gene, "_gt_effect.tsv.gz"), sep = '\t')
  }
  
} else {'no targets to test'}



## ---- SessionInfo

sessionInfo()
 