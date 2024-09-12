
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))
suppressMessages(library(glmnet))
suppressMessages(library(doParallel))
suppressMessages(library(ggpubr))
suppressMessages(library(tidyverse))

section_name <- "12b_downstream_downstream_lm"
analysis_name <- "pipeline"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
ExperimentName <- "Magpie"
magpie_version <- "2022-08-15" #(magpie)

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

#sc_cor_thresh <- sig_cor_thresh_downstream


annotations <- c("omnipath_protein_complex",
                 "dorothea",
                 "sc_coexpressed_gene",
                 "gprofiler",
                 "msigdb_hallmark_gene")



## ---- LogisticRegression

#run lm
glm_out_fnm <- paste0(outdir, "/", magpie_version,  "_glm_out.RDS")
lfc_res_out_fnm <- paste0(outdir, '/', magpie_version, '_lfc_res.tsv.gz')
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
      fnm <- sort(grep(list.files(file.path(OutFolder, "Preprocess_External_Datasets", annotation),
                                  pattern = paste0(annotation, "_pairs_")),
                       pattern = paste0(tolower(ExperimentName), "-", magpie_version, "-version.tsv"), value = T), decreasing = T)[1]
      assign(paste0(annotation, "_pairs"),
             fread(file.path(OutFolder, "Preprocess_External_Datasets", annotation, fnm), header = T))
    }
  }
  
  ## ---- LoadMagpieData
  
  print("loading magpie data...")
  #config_name <- paste0(magpie_version, "_min_deg-", min_degs, "_min_targets-", min_targets, "_min_cells-", min_cells)
  lfc_res <- fread(file.path(OutFolder, "12a_downstream_gene_embedding", paste0(magpie_version, "_downstream_downstream_cor_with_adj_pval.tsv.gz")))
  lfc_res <- lfc_res %>%
    dplyr::select(
      'downstream_gene_1', 'downstream_gene_2',
      'cor' = 'iPSC_cor'
    ) %>%
    dplyr::filter(downstream_gene_1 > downstream_gene_2) %>%
    as.data.frame()
  #print("number of downstream-downstream pairs:")
  #print(dim(lfc_res)[1])
  lfc_res$to_match_1 <- paste(lfc_res$downstream_gene_1, lfc_res$downstream_gene_2, sep = ":")
  lfc_res$to_match_2 <- paste(lfc_res$downstream_gene_2, lfc_res$downstream_gene_1, sep = ":")
  
  
  
  ## ---- MakeTable4LM
  
  
  print("Adding external data...")
  ## add to interaction info
  
  vars2regress <- c()
  vars2regress_pretty <- c()
  
  #coessentiality_module_pairs 
  if ("dorothea" %in% annotations){
    shared_tf <- split.data.frame(dorothea_pairs, f = dorothea_pairs$source)
    shared_tf <- lapply(shared_tf, FUN = function(df){
      if (dim(df)[1] >= 2){
        rtn <- as.data.frame(t(combn(df$target, 2)))
      } else{
        rtn <- NULL
      }
    })
    shared_tf <- as.data.frame(bind_rows(shared_tf))
    names(shared_tf) <- c('gene_1', 'gene_2')
    shared_tf$to_match <- paste(shared_tf$gene_1, shared_tf$gene_2, sep = ":")
    lfc_res$is_tf_downstream_pair <- (lfc_res$to_match_1 %in% shared_tf$to_match | lfc_res$to_match_2 %in% shared_tf$to_match)
    vars2regress <- c(vars2regress, 'is_tf_downstream_pair')
    vars2regress_pretty <- c(vars2regress_pretty, 'Shared transcription factor')
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
    lfc_res$is_sc_coexpressed_pair <- (abs(lfc_res$sc_expression_correlation) > sig_cor_thresh_downstream)#bulk co-expression
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
  
  # remove anything without enough power
  num_distinct <- sapply(lfc_res, n_distinct)
  cols2remove <- names(num_distinct)[which(num_distinct == 1)]
  
  
  ## load lfc_res from here
  lfc_res$var2regress <- (abs(lfc_res$cor) > sig_cor_thresh_downstream)
  print("Number of target-target pairs:")
  print(dim(lfc_res)[1])
  print("Number of targets:")
  print(length(unique(c(lfc_res$downstream_gene_1, lfc_res$downstream_gene_2))))
  
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
yrange <- max(df4plot$elastic_net_coef) - min(c(df4plot$elastic_net_coef, 0))
pretty_plot <- ggplot(df4plot, aes(x = label, y = elastic_net_coef)) + 
  ylab('GLM Coefficient') + xlab('') + 
  geom_bar(stat="identity", fill = downstream_col) + theme_bw()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
        plot.margin = margin(1,1,0, 2, "cm"), 
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5))
print(pretty_plot)
dev.off()


## ---- SessionInfo()

sessionInfo()
