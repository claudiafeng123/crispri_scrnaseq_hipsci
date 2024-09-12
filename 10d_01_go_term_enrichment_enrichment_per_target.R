
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))
suppressMessages(library(glmnet))
suppressMessages(library(doParallel))
suppressMessages(library(ggpubr))
suppressMessages(library(tidyverse))
suppressMessages(library(doParallel))


section_name <- "10d_go_term_enrichment"
analysis_name <- "pipeline"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
ExperimentName <- "Magpie"
date <- "2022-08-15" #(magpie)
#by <- 'target'; target_gene <- "CD72"
#by <- 'downstream_gene_name'; downstream_gene <- "ZNF483"

args <- commandArgs(trailingOnly = TRUE)
#print(args)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  ##needed arguments
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--experiment_name"){ ExperimentName <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
  if (arg == "--by"){ by <- args[[ind + 1]] }
  if (arg == "--downstream_gene"){ downstream_gene <- args[[ind + 1]] }
  if (arg == "--target_gene"){ target_gene <- args[[ind + 1]] }
}

# set relevent i/o paths
source(file.path(HomeFolder, paste0("scripts/io/", ExperimentName, "_io.R")))
source(file.path(HomeFolder, "scripts/io/Magpie_Utils.R"))
outdir <- file.path(OutFolder, "", section_name, paste0('by_', by))
plotsdir <- file.path(HTMLFolder, analysis_name, section_name)
if(!dir.exists(outdir)) dir.create(outdir)
if(!dir.exists(plotsdir)) dir.create(plotsdir)
print(section_name)


## ---- LoadData

gprofiler_annotations <- readRDS(paste0(ResourcesFolder, "gprofiler/gprofiler_full_hsapiens.ENSG.RDS"))
gprofiler_terms <- fread(paste0(OutFolder, "Preprocess_External_Datasets/gprofiler/go_terms_2022-10-10.tsv"))

magpie_lfcs <- fread(file.path(OutFolder, "6b_calc_lfcs_transcriptome_wide_by_gene", paste0(date, "_combined_with_adj_pval.tsv.gz")))
if (by == "target"){
  print(paste0("target gene: ", target_gene))
  magpie_lfcs <- dplyr::filter(magpie_lfcs, target == target_gene)
} else if (by == "downstream_gene_name"){
  print(paste0("downstream gene: ", downstream_gene))
  dg <- downstream_gene
  magpie_lfcs <- dplyr::filter(magpie_lfcs, downstream_gene_name == dg)
}

## ---- t-test

if (dim(magpie_lfcs)[1] > 0){
  print(paste0('testing ', dim(magpie_lfcs)[1], ' genes'))
  # only test the terms with at least 10 genes in our data
  gprofiler_terms2test <- lapply(1:length(gprofiler_annotations), FUN = function(i){
    if (by == "target"){
      gprofiler_annotations[[i]][gprofiler_annotations[[i]] %in% magpie_lfcs$downstream_gene_name]
    } else if (by == "downstream_gene_name"){
      gprofiler_annotations[[i]][gprofiler_annotations[[i]] %in% magpie_lfcs$target]
    }
    
  })
  names(gprofiler_terms2test) <- names(gprofiler_annotations)
  gprofiler_terms2test <- gprofiler_terms2test[sapply(gprofiler_terms2test, FUN = function(l){length(l) >= 10})]
  
  registerDoParallel(50)
  res <- foreach (i = 1:length(gprofiler_terms2test)) %dopar% {
    if (by == "target"){
      t_test_res <- t.test(
        magpie_lfcs$lfc[which(magpie_lfcs$downstream_gene_name %in% gprofiler_terms2test[[i]])],
        magpie_lfcs$lfc[which(!(magpie_lfcs$downstream_gene_name %in% gprofiler_terms2test[[i]]))]
      )
      } else if (by == "downstream_gene_name"){
        t_test_res <- t.test(
          magpie_lfcs$lfc[which(magpie_lfcs$target %in% gprofiler_terms2test[[i]])],
          magpie_lfcs$lfc[which(!(magpie_lfcs$target %in% gprofiler_terms2test[[i]]))]
        )}
    
    rtn <- c(t_stat = as.numeric(t_test_res$statistic), pval = t_test_res$p.value, 
             mean_in = as.numeric(t_test_res$estimate[1]), mean_out = as.numeric(t_test_res$estimate[2]))
    return(rtn)
  }
  
  if (by == "target"){
    res <- as.data.frame(bind_rows(res)) %>%
      mutate(target = target_gene,
             go_id = names(gprofiler_terms2test)) %>%
      mutate(go_term = gprofiler_terms$Term[match(go_id, gprofiler_terms$GO_ID)],
             go_source = gprofiler_terms$Ontology[match(go_id, gprofiler_terms$GO_ID)])
    
    ## ---- write
    
    fwrite(res, paste0(outdir, '/', date, "_", target_gene, "_go_enrichment.tsv"), sep = '\t')
  } else {
    res <- as.data.frame(bind_rows(res)) %>%
      mutate(downstream_gene_name = downstream_gene,
             go_id = names(gprofiler_terms2test)) %>%
      mutate(go_term = gprofiler_terms$Term[match(go_id, gprofiler_terms$GO_ID)],
             go_source = gprofiler_terms$Ontology[match(go_id, gprofiler_terms$GO_ID)])
    
    ## ---- write
    
    print(paste0("writing to..."))
    print(paste0(outdir, '/', date, "_", downstream_gene, "_go_enrichment.tsv"))
    fwrite(res, paste0(outdir, '/', date, "_", downstream_gene, "_go_enrichment.tsv"), sep = '\t')
    
  }
  
}


## ---- SessionInfo()

sessionInfo()
