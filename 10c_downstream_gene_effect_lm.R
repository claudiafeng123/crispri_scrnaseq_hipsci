
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))
suppressMessages(library(glmnet))
suppressMessages(library(doParallel))
suppressMessages(library(ggpubr))
suppressMessages(library(tidyverse))

section_name <- "10c_downstream_gene_effect_lm"
analysis_name <- "pipeline"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
ExperimentName <- "Magpie"
magpie_version <- "2022-08-15" #(magpie)
out_date <- "2024-05-17"

io_path <- paste0(HomeFolder, "scripts/io/", ExperimentName, "_io.R")
utils_path <- paste0(HomeFolder, "scripts/io/Magpie_Utils.R")
REDO <- TRUE


args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  ##needed arguments
  if (arg == "--date"){ magpie_version <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
  ## not needed
  if (arg == "--utils_path"){ utils_path <- args[[ind + 1]] }
  if (arg == "--io_path"){ io_path <- args[[ind + 1]] }
}

# set relevent i/o paths
source(io_path)
source(utils_path)
outdir <- file.path(OutFolder, "", section_name)
plotsdir <- file.path(HTMLFolder, analysis_name, section_name)
if(!dir.exists(outdir)) dir.create(outdir)
if(!dir.exists(plotsdir)) dir.create(plotsdir)
print(section_name)
print(paste0("REDO: ", REDO))


annotations <- c("genomic_distance", 'wge_off_targets', 
                 'dorothea', 'eQTL',
                 "omnipath_protein_complex",
                 "sc_coexpressed_gene",
                 "coessentiality_modules",
                 "gprofiler",
                 "msigdb_hallmark_gene")



glm_out_fnm <- paste0(outdir, "/", out_date,  "_glm_out.RDS")
lfc_res_out_fnm <- paste0(outdir, '/', out_date, '_lfc_res.tsv.gz')
if (!file.exists(glm_out_fnm) | !(file.exists(lfc_res_out_fnm )) | REDO == TRUE){
  
  ## ---- LoadData
  #external datasets
  for (annotation in annotations){
    print(annotation)
    if (annotation == "omnipath_interaction"){
      fnm <- sort(grep(list.files(file.path(OutFolder, "Preprocess_External_Datasets", "omnipath"),
                                  pattern = "omnipath_interaction_pairs_data-from"),
                       pattern = paste0(tolower(ExperimentName), "-", magpie_version, "-version.tsv"), value = T), decreasing = T)[1]
      assign("omnipath_interaction_pairs",
             fread(file.path(OutFolder, "Preprocess_External_Datasets", "omnipath", fnm)))
    } else if (annotation == "omnipath_protein_complex"){
      fnm <- sort(grep(list.files(file.path(OutFolder, "Preprocess_External_Datasets", "omnipath"),
                                  pattern = "omnipath_protein_complex_pairs_data-from"),
                       pattern = paste0(tolower(ExperimentName), "-", magpie_version, "-version.tsv"), value = T), decreasing = T)[1]
      assign("omnipath_protein_complex_pairs",
             fread(file.path(OutFolder, "Preprocess_External_Datasets", "omnipath", fnm)))
    } else if (annotation == "sc_coexpressed_gene") {
      fnm <- sort(grep(list.files(file.path(OutFolder, "Preprocess_External_Datasets", "coexpression"),
                                  pattern = "sc_regressed_correlation_pairs_"),
                       pattern = paste0(tolower(ExperimentName), "-", magpie_version, "-version.tsv"), value = T), decreasing = T)[1]
      assign("sc_coexpressed_gene_pairs",
             fread(file.path(OutFolder, "Preprocess_External_Datasets", "coexpression", fnm)))
    } else if (annotation == "bulk_coexpressed_gene"){
      fnm <- sort(grep(list.files(file.path(OutFolder, "Preprocess_External_Datasets", "coexpression"),
                                  pattern = "bulk_correlation_pairs_"),
                       pattern = paste0(tolower(ExperimentName), "-", magpie_version, "-version.tsv"), value = T), decreasing = T)[1]
      assign("bulk_coexpressed_gene_pairs",
             fread(file.path(OutFolder, "Preprocess_External_Datasets", "coexpression", fnm)))
    } else if (annotation == "gprofiler"){
      fnm <- sort(grep(list.files(file.path(OutFolder, "Preprocess_External_Datasets", "gprofiler"),
                                  pattern = "num-common-gprofiler-annotations-pairs_"),
                       pattern = paste0(tolower(ExperimentName), "-", magpie_version, "-version.tsv"), value = T), decreasing = T)[1]
      assign("gprofiler_pairs",
             fread(file.path(OutFolder, "Preprocess_External_Datasets", "gprofiler", fnm)))
    } else if (annotation == "coessentiality_modules") {
      fnm <- sort(grep(list.files(file.path(OutFolder, "Preprocess_External_Datasets", "coessentiality"),
                                  pattern = "num_common_coessentiality_modules_"),
                       pattern = paste0(tolower(ExperimentName), "-", magpie_version, "-version.tsv"), value = T), decreasing = T)[1]
      assign(paste0(annotation, "_pairs"),
             fread(file.path(OutFolder, "Preprocess_External_Datasets", "coessentiality", fnm)))
    } else if (annotation == "msigdb_hallmark_gene"){
      fnm <- sort(grep(list.files(file.path(OutFolder, "Preprocess_External_Datasets", "msigdb_gene_sets"),
                                  pattern = "hallmark_gene_set_pairs_"),
                       pattern = paste0(tolower(ExperimentName), "-", magpie_version, "-version.tsv"), value = T), decreasing = T)[1]
      assign(paste0(annotation, "_pairs"),
             fread(file.path(OutFolder, "Preprocess_External_Datasets", "msigdb_gene_sets", fnm)))
    } else if (annotation == "dorothea"){
      fnm <- sort(grep(list.files(file.path(OutFolder, "Preprocess_External_Datasets", annotation),
                                  pattern = "dorothea_interaction_pairs_"),
                       pattern = paste0(tolower(ExperimentName), "-", magpie_version, "-version.tsv"), value = T), decreasing = T)[1]
      assign(paste0(annotation, "_pairs"),
             fread(file.path(OutFolder, "Preprocess_External_Datasets", annotation, fnm)))
    } else if (annotation == "eQTL"){
      fnm <- sort(grep(list.files(file.path(OutFolder, "Preprocess_External_Datasets", "eQTLs"),
                                  pattern = "eQTL_interaction_pairs_"),
                       pattern = paste0(tolower(ExperimentName), "-", magpie_version, "-version.tsv"), value = T), decreasing = T)[1]
      assign(paste0(annotation, "_pairs"),
             fread(file.path(OutFolder, "Preprocess_External_Datasets", "eQTLs", fnm)))
    } else if (annotation == "off_target"){
      fnm <- sort(grep(list.files(file.path(OutFolder, "Preprocess_External_Datasets", annotation),
                                  pattern = "off_target_pairs_"),
                       pattern = paste0(tolower(ExperimentName), "-", magpie_version, "-version.tsv"), value = T), decreasing = T)[1]
      assign(paste0(annotation, "_pairs"),
             fread(file.path(OutFolder, "Preprocess_External_Datasets", annotation, fnm), header = T))
    } else if (annotation == "gwas_colocalization"){
      fnm <- sort(grep(list.files(file.path(OutFolder, "Preprocess_External_Datasets", annotation),
                                  pattern = paste0(annotation, "_pairs_")),
                       pattern = paste0(tolower(ExperimentName), "-", magpie_version, "-version.tsv"), value = T), decreasing = T)[1]
      assign(paste0(annotation, "_pairs"),
             fread(file.path(OutFolder, "Preprocess_External_Datasets", annotation, fnm), header = T))
    }  else if (annotation == "genomic_distance"){
      fnm <- sort(grep(list.files(file.path(OutFolder, "Preprocess_External_Datasets", 'off_target'),
                                  pattern = paste0(annotation, "_1000_nt")),
                       pattern = magpie_version, value = T), decreasing = T)[1]
      assign(paste0(annotation, "_pairs"),
             fread(file.path(OutFolder, "Preprocess_External_Datasets", 'off_target', fnm), header = T))
    } else if (annotation == "wge_off_targets"){
      fnm <- sort(grep(list.files(file.path(OutFolder, "Preprocess_External_Datasets", 'off_target'),
                                  pattern = paste0(annotation, '_1000_nt_3')),
                       pattern = magpie_version, value = T), decreasing = T)[1]
      assign(paste0(annotation, "_pairs"),
             fread(file.path(OutFolder, "Preprocess_External_Datasets", 'off_target', fnm), header = T))
    }
  }
  
  ## ---- LoadMagpieData
  
  print("loading magpie data...")
  lfc_res <- fread(file.path(OutFolder, "6b_calc_lfcs_transcriptome_wide_by_gene", paste0(magpie_version, '_combined_with_adj_pval.tsv.gz')))
  target_genes <- unique(lfc_res$target)
  lfc_res <- select(lfc_res, c("target", "downstream_gene_name", "lfc", "pval_adj")) %>%
    mutate(to_match_1 = paste0(target, ":", downstream_gene_name),
           to_match_2 = paste0(downstream_gene_name, ":", target)) %>%
    dplyr::filter(target != downstream_gene_name)
  
  ## ---- MakeTable4LM
  
  print("Adding external data...")
  ## add to interaction info
  
  vars2regress <- vars2regress_pretty <- c()
  for (annotation in annotations){
    #coessentiality_module_pairs 
    if (annotation == "coessentiality_modules"){
      coessentiality_modules_pairs$to_match <- paste(coessentiality_modules_pairs$gene_1, coessentiality_modules_pairs$gene_2, sep = ":")
      lfc_res$is_coessentiality_module_pair <- (lfc_res$to_match_1 %in% coessentiality_modules_pairs$to_match | lfc_res$to_match_2 %in% coessentiality_modules_pairs$to_match)
      vars2regress <- c(vars2regress, 'is_coessentiality_module_pair')
      vars2regress_pretty <- c(vars2regress_pretty, 'Co-essential Gene')
    }
    
    #msigdb_hallmark_gene_pairs
    if (annotation == "msigdb_hallmark_gene"){
      msigdb_hallmark_gene_pairs$to_match <- paste(msigdb_hallmark_gene_pairs$gene_1, msigdb_hallmark_gene_pairs$gene_2, sep = ":")
      lfc_res$is_msigdb_hallmark_gene_set_pair <- (lfc_res$to_match_1 %in% msigdb_hallmark_gene_pairs$to_match | lfc_res$to_match_2 %in% msigdb_hallmark_gene_pairs$to_match)
      vars2regress <- c(vars2regress, 'is_msigdb_hallmark_gene_set_pair')
      vars2regress_pretty <- c(vars2regress_pretty, 'Common Pathway')
    }
    
    #coexpression
    #single-cell co-expression
    if (annotation == "sc_coexpressed_gene"){
      sc_coexpressed_gene_pairs$to_match <- paste(sc_coexpressed_gene_pairs$gene_1, sc_coexpressed_gene_pairs$gene_2, sep = ":")
      lfc_res$sc_expression_correlation <- sc_coexpressed_gene_pairs$correlation[match(lfc_res$to_match_1, sc_coexpressed_gene_pairs$to_match)]
      lfc_res$sc_expression_correlation[which(is.na(lfc_res$sc_expression_correlation))] <- 0
      lfc_res$is_sc_coexpressed_pair <- (abs(lfc_res$sc_expression_correlation) > sig_cor_thresh_target)#bulk co-expression
      vars2regress <- c(vars2regress, 'is_sc_coexpressed_pair')
      vars2regress_pretty <- c(vars2regress_pretty, 'Co-expressed Gene')
    }
    
    if (annotation == "bulk_coexpressed_gene"){
      bulk_coexpressed_gene_pairs$to_match <- paste(bulk_coexpressed_gene_pairs$gene_1, bulk_coexpressed_gene_pairs$gene_2, sep = ":")
      lfc_res$bulk_expression_correlation <- bulk_coexpressed_gene_pairs$bulk_cor[match(lfc_res$to_match_1, bulk_coexpressed_gene_pairs$to_match)]
      lfc_res$bulk_expression_correlation[which(is.na(lfc_res$bulk_expression_correlation))] <- 0
      lfc_res$is_bulk_coexpressed_pair <- (abs(lfc_res$bulk_expression_correlation) > sig_cor_thresh_target)
      vars2regress <- c(vars2regress, 'is_bulk_coexpressed_pair')
      vars2regress_pretty <- c(vars2regress_pretty, 'Co-expressed Gene (Bulk across Donors)')
    }
    
    #omnipath protein complexes
    if (annotation == "omnipath_protein_complex"){
      omnipath_protein_complex_pairs$to_match <- paste(omnipath_protein_complex_pairs$gene_1, omnipath_protein_complex_pairs$gene_2, sep = ":")
      lfc_res$is_omnipath_protein_complex_pair <- (lfc_res$to_match_1 %in% omnipath_protein_complex_pairs$to_match | lfc_res$to_match_2 %in% omnipath_protein_complex_pairs$to_match)
      vars2regress <- c(vars2regress, 'is_omnipath_protein_complex_pair')
      vars2regress_pretty <- c(vars2regress_pretty, 'Common Protein Complex')
    }
    
    #omnipath interactions
    if (annotation == "omnipath_interaction"){
      omnipath_interaction_pairs$to_match <- paste(omnipath_interaction_pairs$gene_1, omnipath_interaction_pairs$gene_2, sep = ":")
      lfc_res$is_omnipath_interaction_pair <- (lfc_res$to_match_1 %in% omnipath_interaction_pairs$to_match | lfc_res$to_match_2 %in% omnipath_interaction_pairs$to_match)
      vars2regress <- c(vars2regress, 'is_omnipath_interaction_pair')
      vars2regress_pretty <- c(vars2regress_pretty, 'Known Protein-protein Interaction')
    }
    #GO terms (gprofiler)
    if (annotation == "gprofiler"){
      #lfc_res$num_common_gprofiler_annotations <- gprofiler_pairs$num_common_annots[match(lfc_res$to_match_1, paste(gprofiler_pairs$gene_1, gprofiler_pairs$gene_2, sep = ":"))]
      #lfc_res$num_common_gprofiler_annotations[is.na(lfc_res$num_common_gprofiler_annotations)] <- 0
      #lfc_res$is_gprofiler_pair <- (lfc_res$num_common_gprofiler_annotations > 0)
      #vars2regress <- c(vars2regress, 'is_gprofiler_pair')
      #vars2regress_pretty <- c(vars2regress_pretty, "Common GO annotation")
      for (term_type in c("bp", "cc", "mf")){
        
        eval(parse(text = paste0(
          'num_common_annotations <- gprofiler_pairs$common_', term_type, '_annotations[match(lfc_res$to_match_1, paste(gprofiler_pairs$gene_1, gprofiler_pairs$gene_2, sep = ":"))]'
        )))
        num_common_annotations[is.na(num_common_annotations)] <- 0
        eval(parse(text = paste0(
          'lfc_res$is_gprofiler_', term_type, '_pair <- (num_common_annotations > 0)'
        )))
        vars2regress <- c(vars2regress, paste0('is_gprofiler_', term_type, '_pair'))
        vars2regress_pretty <- c(vars2regress_pretty, ifelse(term_type == "bp",
                                                             "Common biological process", ifelse(term_type == "cc", "Common cellular component", "Common molecular function")))
        
      }
    }
    
    #dorothea
    if (annotation == "dorothea"){
      dorothea_pairs$to_match <- paste0(dorothea_pairs$source, ":", dorothea_pairs$target)
      lfc_res$is_tf_pair <- (lfc_res$to_match_1 %in% dorothea_pairs$to_match)
      vars2regress <- c(vars2regress, 'is_tf_pair')
      vars2regress_pretty <- c(vars2regress_pretty, 'Known Transcription Factor-Target Pair')
    }
    
    #eqtl pairs
    if (annotation == "eQTL"){
      eQTL_pairs$to_match <- paste0(eQTL_pairs$trans_gene, ":", eQTL_pairs$cis_gene)
      lfc_res$is_eQTL_pair <- (lfc_res$to_match_2 %in% eQTL_pairs$to_match)
      vars2regress <- c(vars2regress, 'is_eQTL_pair')
      vars2regress_pretty <- c(vars2regress_pretty, 'Known trans-/cis-eQTL Pair')
    }
    
    #off-target
    if (annotation == "off_target"){
      off_target_pairs$to_match <- paste(off_target_pairs$target_gene, off_target_pairs$off_target_gene, sep = ":")
      lfc_res$is_off_target_pair <- (lfc_res$to_match_1 %in% off_target_pairs$to_match | lfc_res$to_match_2 %in% off_target_pairs$to_match)
      vars2regress <- c(vars2regress, 'is_off_target_pair')
      vars2regress_pretty <- c(vars2regress_pretty, 'Off-target')
    }
    
    
    if (annotation == "genomic_distance"){
      genomic_distance_pairs$to_match <- paste(genomic_distance_pairs$gene, genomic_distance_pairs$nearby_gene, sep = ":")
      lfc_res$is_genomic_distance_pair <- (lfc_res$to_match_1 %in% genomic_distance_pairs$to_match | lfc_res$to_match_2 %in% genomic_distance_pairs$to_match)
      vars2regress <- c(vars2regress, 'is_genomic_distance_pair')
      vars2regress_pretty <- c(vars2regress_pretty, "Genomic Proximity")
    }
    
    
    if (annotation == "wge_off_targets"){
      wge_off_targets_pairs$to_match <- paste(wge_off_targets_pairs$gene, wge_off_targets_pairs$off_target_genes, sep = ":")
      lfc_res$is_wge_off_target_pair <- (lfc_res$to_match_1 %in% wge_off_targets_pairs$to_match | lfc_res$to_match_2 %in% wge_off_targets_pairs$to_match)
      vars2regress <- c(vars2regress, 'is_wge_off_target_pair')
      vars2regress_pretty <- c(vars2regress_pretty, "Sequence Similarity")
    }
    
    
  }
  df4lm <- lfc_res %>% 
    dplyr::select(c("target", "downstream_gene_name", 
                    'lfc',
                    'pval_adj',
                    all_of(vars2regress))) %>%
    mutate('var2regress' = ifelse(pval_adj < sig_pval_thresh, 1, 0))
  for (v in vars2regress){
    eval(parse(text = paste0(
      'df4lm$', v, ' <- ifelse(df4lm$', v, ', 1, 0)'
    )))
  }
  
  
  ## ---- LogisticRegression
  
  
  lm4model <- paste('var2regress', '~', paste(vars2regress, collapse = ' + '))
  fit_elastic_net <- glmnet::cv.glmnet(y =df4lm$var2regress ,
                                       x = as.matrix(df4lm %>% dplyr::select(all_of(vars2regress))), 
                                       family = "poisson")
  df4lm <- df4lm %>%
    mutate(var2regress = ifelse(pval_adj < sig_pval_thresh, 1, 0))
  fit_glm <- glm(as.formula(lm4model), family=binomial(link='logit'),
                 data = df4lm)
  
  
  lm_result <- data.frame(
    var = vars2regress,
    label = vars2regress_pretty,
    elastic_net_coef = coef(fit_elastic_net)[-1],
    glm_coef = coefficients(fit_glm)[-1],
    stderr = summary(fit_glm)$coefficients[-1,2],
    z_value = summary(fit_glm)$coefficients[-1,3],
    p_val = summary(fit_glm)$coefficients[-1,4]
  ) 
  # save
  saveRDS(lm_result, file = glm_out_fnm)
  
  df4lm <- df4lm %>%
    mutate(y4glm = ifelse(pval_adj < sig_pval_thresh, 1, 0),
           y4elastic_net = abs(lfc)) %>%
    dplyr::select('lfc',
                  'pval_adj',
                  'y4elastic_net',
                  'y4glm',
                  all_of(vars2regress))
  df4lm$glm_res <- residuals(fit_glm)
  predicted_values <- predict(fit_elastic_net, as.matrix(df4lm %>% dplyr::select(all_of(vars2regress))))
  df4lm$elastic_net_res <- df4lm$y4elastic_net-predicted_values
  ## save the residuals as well
  fwrite(lfc_res, lfc_res_out_fnm, compress = "gzip", sep = "\t")
  
} else {
  lm_result <- readRDS(glm_out_fnm)
}


## actually do things

#plot
pdf(paste0(plotsdir, "/", out_date, '_elastic_net_coef.pdf'), width = 6, height = 6)
pretty_plot <- ggplot(lm_result, aes(x = label, y = elastic_net_coef)) + 
  ylab('GLM Coefficient') + xlab('') + ggtitle("Functional Enrichment of\nTarget Gene-Downstream Gene Effects") + 
  geom_bar(stat="identity", fill = target_downstream_col) + theme_bw()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
        plot.margin = margin(1,1,0, 2, "cm"), 
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5))
print(pretty_plot)
dev.off()

pdf(paste0(plotsdir, "/", magpie_version, '_glm.pdf'), width = 6, height = 6)
pretty_plot <- ggplot(lm_result, aes(x = label, y = glm_coef)) + 
  ylab('GLM Coefficient') + xlab('') + ggtitle("Functional Enrichment of\nTarget Gene-Downstream Gene Effects") + 
  geom_bar(stat="identity", fill = target_downstream_col) + theme_bw()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
        plot.margin = margin(1,1,0, 2, "cm"), 
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5))
print(pretty_plot)

## ---- SessionInfo()

sessionInfo()

## ---- Scratch

## just rerunning the model without loading the data
if (F){
  df4lm <- fread(lfc_res_out_fnm)
  vars2regress <- grep(colnames(df4lm), pattern = 'is_', value = T)
  df4lm$var2regress <- ifelse(abs(df4lm$lfc), 1, 0)
  df4lm <- filter(df4lm, target != downstream_gene_name)
  
  lm4model <- paste('var2regress', '~', paste(vars2regress, collapse = ' + '))
  fit <- glmnet::cv.glmnet(y =df4lm$var2regress ,
                                       x = as.matrix(df4lm %>% dplyr::select(all_of(vars2regress))), 
                                       family = "binomial")
  
  
  df4lm <- df4lm %>%
    mutate(var2regress = ifelse(pval_adj < sig_pval_thresh, 1, 0))
  fit_glm <- glm(as.formula(lm4model), family=binomial(link='logit'),
                 data = df4lm)
  
}

