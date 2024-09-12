
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))
suppressMessages(library(glmnet))
suppressMessages(library(doParallel))
suppressMessages(library(ggpubr))
suppressMessages(library(tidyverse))

section_name <- "11b_target_target_correlation_lm"
analysis_name <- "pipeline"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
ExperimentName <- "Magpie"
magpie_version <- "2022-08-15" #(magpie)
date_out <- '2024-05-17'

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
  if (arg == "--REDO"){ REDO <- as.logical(args[[ind + 1]]) }
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

sc_cor_thresh <- sig_cor_thresh_target


annotations <- c("genomic_distance", 'wge_off_targets', #off targets 
                 "omnipath_protein_complex",
                 "sc_coexpressed_gene",
                 "coessentiality_modules",
                 "gprofiler",
                 "msigdb_hallmark_gene")



## ---- LogisticRegression

#run lm
glm_out_fnm <- paste0(outdir, "/", date_out,  "_glm_out.RDS")
lfc_res_out_fnm <- paste0(outdir, '/', date_out, '_lfc_res.tsv.gz')
if (!file.exists(glm_out_fnm) | !(file.exists(lfc_res_out_fnm )) | REDO == TRUE){
  
  ## ---- LoadData
  
  
  #external datasets
  for (annotation in annotations){
    print(annotation)
    if (annotation == "omnipath_protein_interaction"){
      fnm <- sort(grep(list.files(file.path(OutFolder, "Preprocess_External_Datasets", "omnipath"),
                                  pattern = "omnipath_interaction_pairs_data-from"),
                       pattern = paste0(tolower(ExperimentName), "-", magpie_version, "-version.tsv"), value = T), decreasing = T)[1]
      assign("omnipath_protein_interaction_pairs",
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
    } else if (annotation == "genomic_distance"){
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
  #config_name <- paste0(magpie_version, "_min_deg-", min_degs, "_min_targets-", min_targets, "_min_cells-", min_cells)
  outdir <- file.path(OutFolder, "", section_name)#, config_name)
  if(!dir.exists(outdir)) dir.create(outdir)
  lfc_res <- fread(file.path(OutFolder, "11a_target_target_cor", paste0(magpie_version, "_target_target_cor.tsv.gz")))
  lfc_res <- lfc_res %>%
    dplyr::select(
      'target_1' = 'gene_1',
      'target_2' = 'gene_2',
      'cor' = 'cor_all'
    ) %>%
    dplyr::filter(target_1 > target_2) %>%
    as.data.frame()
  print("number of targets/target pairs:")
  #print(dim(lfc_res)[1])
  lfc_res$to_match_1 <- paste(lfc_res$target_1, lfc_res$target_2, sep = ":")
  lfc_res$to_match_2 <- paste(lfc_res$target_2, lfc_res$target_1, sep = ":")
  
  
  
  ## ---- MakeTable4LM
  
  
  print("Adding external data...")
  ## add to interaction info
  
  vars2regress <- c()
  vars2regress_pretty <- c()
  
  #coessentiality_module_pairs 
  if ("coessentiality_modules" %in% annotations){
    coessentiality_modules_pairs$to_match <- paste(coessentiality_modules_pairs$gene_1, coessentiality_modules_pairs$gene_2, sep = ":")
    lfc_res$is_coessentiality_module_pair <- (lfc_res$to_match_1 %in% coessentiality_modules_pairs$to_match | lfc_res$to_match_2 %in% coessentiality_modules_pairs$to_match)
    vars2regress <- c(vars2regress, 'is_coessentiality_module_pair')
    vars2regress_pretty <- c(vars2regress_pretty, 'Co-essential')
  }
  
  #msigdb_hallmark_gene_pairs
  if ("msigdb_hallmark_gene" %in% annotations){
    msigdb_hallmark_gene_pairs$to_match <- paste(msigdb_hallmark_gene_pairs$gene_1, msigdb_hallmark_gene_pairs$gene_2, sep = ":")
    lfc_res$is_msigdb_hallmark_gene_set_pair <- (lfc_res$to_match_1 %in% msigdb_hallmark_gene_pairs$to_match | lfc_res$to_match_2 %in% msigdb_hallmark_gene_pairs$to_match)
    vars2regress <- c(vars2regress, 'is_msigdb_hallmark_gene_set_pair')
    vars2regress_pretty <- c(vars2regress_pretty, 'Common pathway')
  }
  
  #coexpression
  #single-cell co-expression
  if ("sc_coexpressed_gene" %in% annotations){
    sc_coexpressed_gene_pairs$to_match <- paste(sc_coexpressed_gene_pairs$gene_1, sc_coexpressed_gene_pairs$gene_2, sep = ":")
    lfc_res$sc_expression_correlation <- sc_coexpressed_gene_pairs$correlation[match(lfc_res$to_match_1, sc_coexpressed_gene_pairs$to_match)]
    lfc_res$sc_expression_correlation[which(is.na(lfc_res$sc_expression_correlation))] <- 0
    lfc_res$is_sc_coexpressed_pair <- (abs(lfc_res$sc_expression_correlation) > sc_cor_thresh)#bulk co-expression
    vars2regress <- c(vars2regress, 'is_sc_coexpressed_pair')
    vars2regress_pretty <- c(vars2regress_pretty, "Co-expressed")
  }
  
  
  #omnipath protein complexes
  if ("omnipath_protein_complex" %in% annotations){
    omnipath_protein_complex_pairs$to_match <- paste(omnipath_protein_complex_pairs$gene_1, omnipath_protein_complex_pairs$gene_2, sep = ":")
    lfc_res$is_omnipath_protein_complex_pair <- (lfc_res$to_match_1 %in% omnipath_protein_complex_pairs$to_match | lfc_res$to_match_2 %in% omnipath_protein_complex_pairs$to_match)
    vars2regress <- c(vars2regress, 'is_omnipath_protein_complex_pair')
    vars2regress_pretty <- c(vars2regress_pretty, "Common protein complex")
  }
  
  #omnipath interactions
  if ("omnipath_interaction" %in% annotations){
    omnipath_interaction_pairs$to_match <- paste(omnipath_interaction_pairs$gene_1, omnipath_interaction_pairs$gene_2, sep = ":")
    lfc_res$is_omnipath_interaction_pair <- (lfc_res$to_match_1 %in% omnipath_interaction_pairs$to_match | lfc_res$to_match_2 %in% omnipath_interaction_pairs$to_match)
    print("omnipath interactions summary:")
    #table(lfc_res$is_omnipath_interaction_pair)
    vars2regress <- c(vars2regress, 'is_omnipath_interaction_pair')
    vars2regress_pretty <- c(vars2regress_pretty, 'Known Protein Interaction')
  }
  
  #GO terms (gprofiler)
  if ("gprofiler" %in% annotations){
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
  
  
  #off-target
  if ("off_target" %in% annotations){
    off_target_pairs$to_match <- paste(off_target_pairs$target_gene, off_target_pairs$off_target_gene, sep = ":")
    lfc_res$is_off_target_pair <- (lfc_res$to_match_1 %in% off_target_pairs$to_match | lfc_res$to_match_2 %in% off_target_pairs$to_match)
    print("off-target summary:")
    #table(lfc_res$is_off_target_pair)
    vars2regress <- c(vars2regress, 'is_off_target_pair')
    vars2regress_pretty <- c(vars2regress_pretty, "Off-target")
  }
  
  if ("genomic_distance" %in% annotations){
    genomic_distance_pairs$to_match <- paste(genomic_distance_pairs$gene, genomic_distance_pairs$nearby_gene, sep = ":")
    lfc_res$is_genomic_distance_pair <- (lfc_res$to_match_1 %in% genomic_distance_pairs$to_match | lfc_res$to_match_2 %in% genomic_distance_pairs$to_match)
    vars2regress <- c(vars2regress, 'is_genomic_distance_pair')
    vars2regress_pretty <- c(vars2regress_pretty, "Genomic Proximity")
  }
  
  
  if ("wge_off_targets" %in% annotations){
    wge_off_targets_pairs$to_match <- paste(wge_off_targets_pairs$gene, wge_off_targets_pairs$off_target_genes, sep = ":")
    lfc_res$is_wge_off_target_pair <- (lfc_res$to_match_1 %in% wge_off_targets_pairs$to_match | lfc_res$to_match_2 %in% wge_off_targets_pairs$to_match)
    vars2regress <- c(vars2regress, 'is_wge_off_target_pair')
    vars2regress_pretty <- c(vars2regress_pretty, "Sequence Similarity")
  }
  
  
  # remove anything without enough power
  num_distinct <- sapply(lfc_res, n_distinct)
  cols2remove <- names(num_distinct)[which(num_distinct == 1)]
  
  
  ## load lfc_res from here
  lfc_res$var2regress <- (abs(lfc_res$cor) > sig_cor_thresh_target)
  print("Number of target-target pairs:")
  print(dim(lfc_res)[1])
  print("Number of targets:")
  print(length(unique(lfc_res$target_1)))
  
  lm4model <- paste('var2regress', '~', paste(vars2regress, collapse = ' + '))
  df4lm <- lfc_res %>% dplyr::select(c('var2regress', all_of(vars2regress)))
  ## convert T/F to 1/0
  for (v in colnames(df4lm)){
    eval(parse(text = paste0(
      'df4lm$', v, ' <- ifelse(df4lm$', v, ', 1, 0)'
    )))
  }
  
  print(paste0("running lm..."))
  fit_elastic_net <- glmnet::cv.glmnet(y =df4lm$var2regress ,
                           x = as.matrix(df4lm %>% dplyr::select(all_of(vars2regress))), 
                           family = "binomial", alpha = 0.1)
  ## sorry attachment to the without the penalty
  fit_glm <- glm(as.formula(lm4model), family=binomial(link='logit'),
                 data = df4lm)
  lm_result <- data.frame(
    var = vars2regress,
    label = vars2regress_pretty,
    elastic_net_coef = coef(fit_elastic_net)[-1],
    glm_coef = coefficients(fit_glm)[-1],
    glm_stderr = summary(fit_glm)$coefficients[-1,2],
    glm_z_value = summary(fit_glm)$coefficients[-1,3],
    glm_p_val = summary(fit_glm)$coefficients[-1,4]
  ) 
  # save
  saveRDS(lm_result, file = glm_out_fnm)
  
  df2save <- lfc_res %>%
    mutate(
      glm_res = residuals(fit_glm)
    )
  predicted_values <- predict(fit_elastic_net, as.matrix(df4lm %>% dplyr::select(all_of(vars2regress))))
  df2save$elastic_net_res <- df4lm$var2regress-predicted_values
  ## save the residuals as well
  fwrite( df2save, lfc_res_out_fnm, compress = "gzip", sep = "\t")
  
  
} else {
  lm_result <- readRDS(glm_out_fnm)
}

#plot
pdf(paste0(plotsdir, "/", magpie_version, '_glm_coef.pdf'), width = 6, height = 6)
df4plot <- lm_result
yrange <- max(df4plot$elastic_net_) - min(c(df4plot$elastic_net_coef, 0))
pretty_plot <- ggplot(df4plot, aes(x = label, y = elastic_net_coef)) + 
  ylab('GLM Coefficient') + xlab('') + 
  geom_bar(stat="identity", fill = target_col) + theme_bw()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
        plot.margin = margin(1,1,0, 2, "cm"), 
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5))
print(pretty_plot)
dev.off()

## ---- Converse



df4roc <- target_target_cor_anno %>%
  dplyr::filter(gene_1 > gene_2) %>%
  dplyr::select(c("gene_1", "gene_2", "cor_all", 'is_protein_complex_pair'))

claudia <- calc_roc(tbl = df4roc, 
                    true_col = is_protein_complex_pair,
                    var_col = cor_all)
ggplot(claudia, aes(x = fpr, y = tpr)) + 
  geom_point() + geom_line() + 
  geom_abline(slope = 1, col = 'gray', lty = 3) + theme_bw()
fwrite(claudia, '/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Magpie/11b_target_target_correlation_lm/2022-08-15_roc_complexes.tsv', sep = '\t')

p <- ggplot(claudia, aes(x = precision, y = tpr, label = thresh)) + 
  geom_path()
+ geom_line() + 
  geom_abline(slope = 1, col = 'gray', lty = 3) + theme_bw()
plotly::ggplotly(p)

ggplot(claudia, aes(x = precision, y = tpr, label = thresh)) + 
  geom_path()
  geom_point() + geom_line() + 
  geom_abline(slope = 1, col = 'gray', lty = 3) + theme_bw()

## ---- SessionInfo()

sessionInfo()
