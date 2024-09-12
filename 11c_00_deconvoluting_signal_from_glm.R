#xvfb-run -a /software/R-4.1.3/bin/R

HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
section_name <- "11c_deconvoluting_signal_from_glm"
date <- "2022-08-15"
ExperimentName <- "Magpie"

library(data.table)
library(dplyr)
library(doParallel)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggExtra)
library(ggpubr)

# options
args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
  if (arg == "--experiment_name"){ ExperimentName <- args[[ind + 1]] }
}

# i/o
source(file.path(HomeFolder, paste0("scripts/io/", ExperimentName, "_io.R")))
source(file.path(HomeFolder, "scripts/io/Magpie_Utils.R"))
outdir <- file.path(OutFolder, section_name)
plotsdir <- file.path(HTMLFolder, "pipeline/", section_name)
if(!dir.exists(outdir)) dir.create(outdir, recursive = T)
if(!dir.exists(plotsdir)) dir.create(plotsdir, recursive = T)
#if(!dir.exists(paste0(plotsdir, "/downstream_scatters/"))) dir.create(paste0(plotsdir, "/downstream_scatters/"), recursive = T)
if (!exists("min_degs")){min_degs <- min_deg_per_target}
if (!exists("min_cells")){min_cells <- min_cells_per_gene}

## ---- LoadData

target_meta <- fread(file.path(OutFolder, "7b_target_summary", paste0(date, "_target_meta_data.csv")))

anno_target_target_fnm <- paste0(OutFolder, section_name, "/", date, "_target_target_cor_with_adj_pval.tsv.gz")
if (file.exists(paste0(OutFolder, section_name, "/", date, "_target_target_cor_with_adj_pval.tsv.gz"))){
  target_target_cor_anno <- fread(anno_target_target_fnm)
} else {
  target_target_cor <- fread(paste0(OutFolder, '11a_target_target_cor', "/", date, "_target_target_cor.tsv.gz"))
  ## add in adjusted p-values
  target_target_cor$cor_all_pval_adj <- p.adjust(target_target_cor$cor_all_pval, method = "BH")
  target_target_cor$cor_all_deg_pval_adj <- p.adjust(target_target_cor$cor_all_deg_pval, method = "BH")
  target_target_cor$cor_all_deg_excl_target_pval_adj <- p.adjust(target_target_cor$cor_all_deg_excl_target_pval, method = "BH")
  target_target_cor$cor_common_deg_pval_adj <- p.adjust(target_target_cor$cor_common_deg_target_pval, method = "BH")
  target_target_cor$cor_common_deg_excl_target_pval_adj <- p.adjust(target_target_cor$cor_common_deg_excl_target_pval, method = "BH")
  
  pair_anno <- fread(paste0(OutFolder, "11b_target_target_correlation_lm/", date, "_lfc_res.tsv.gz"))
  target_target_cor_anno <- left_join(target_target_cor, pair_anno,
                                      by = c("gene_1" = "target_1", "gene_2" = "target_2"))
  
  ## add in target-target correlation from rpe1 and k562 cells
  for (cell_type in c("K562_essential", "K562_gwps", "RPE1_raw")){
    cell_type_target_cor <- fread(paste0(ProjectFolder, "outs/Other/Replogle/6a_analysis_heatmap/2022-09-12_", cell_type, "_target_target_cor.tsv.gz"))
    eval(parse(text = paste0("target_target_cor_anno <- left_join(target_target_cor_anno, ", 
                                       " cell_type_target_cor %>% ", 
                                          "dplyr::select(c('gene_1' = 'target_1', 'gene_2' = 'target_2', '", paste0(cell_type, '_cor'), "' = 'cor_all')))"
    )))
  }
  
  
  fwrite(target_target_cor_anno, anno_target_target_fnm)
}


## target lfcs
magpie_lfcs_split <- fread(file.path(OutFolder, "6b_calc_lfcs_transcriptome_wide_by_gene", paste0(date, "_combined_with_adj_pval.tsv.gz")))
magpie_lfcs_split <- split.data.frame(magpie_lfcs_split, f = magpie_lfcs_split$target)


#off-target
off_target <- fread(paste0(OutFolder, "Preprocess_External_Datasets/off_target/off_target_pairs_magpie-", date, "-version.tsv.gz"),
                    header = T) %>% filter(target_gene != off_target_gene)

#complexes
complex_pair_fnm <- sort(grep(list.files(file.path(OutFolder, 'Preprocess_External_Datasets', 'omnipath'), 
                                    pattern = "complex_pairs"),
                         pattern = date, value = T))[1]
complex_pairs <- fread(file.path(OutFolder, 'Preprocess_External_Datasets', 'omnipath', complex_pair_fnm))
complex_meta_fnm <- sort(grep(list.files(file.path(OutFolder, "Preprocess_External_Datasets", "omnipath"), pattern = date), pattern = "complex_meta", value = T), decreasing = T)[1]
complex_meta <- fread(file.path(OutFolder, "Preprocess_External_Datasets/omnipath/", complex_meta_fnm))

# coessentiality
coessentiality_modules <- fread(paste0(ResourcesFolder, 'coessentiality/Wainberg2021_Supplementary_Data_2.csv'))
module_meta <- select(coessentiality_modules, c("Module #", "Most-enriched GO term", "Most-enriched p"))
coessentiality_modules <- mapply(1:dim(coessentiality_modules)[1], FUN = function(i){
  genes_in_module <- as.character(coessentiality_modules[i,-(1:13)])
  genes_in_module <- genes_in_module[which(genes_in_module != "NA")]
  genes_in_module <- genes_in_module[which(genes_in_module %in% target_meta$gene)]
  if (length(genes_in_module) > 1){
    return(genes_in_module)
  } else {return(NULL)}
}, SIMPLIFY = F)
names(coessentiality_modules) <- paste("Module:", module_meta$`Module #`)
modules2consider <- which(unlist(lapply(coessentiality_modules, length)) > min_genes_per_complex)
module_meta <- module_meta[modules2consider,]
coessentiality_modules <- coessentiality_modules[modules2consider]

gprofiler_annotations <- readRDS(paste0(ResourcesFolder, "gprofiler/gprofiler_full_hsapiens.ENSG.RDS"))
gprofiler_terms <- fread(paste0(OutFolder, "Preprocess_External_Datasets/gprofiler/go_terms_2023-10-18.tsv"))

## ---- Protein Complexes

subsection_name <- '11c_01_deconvoluting_signal_from_glm_protein_complexes'
params_list <- list(home_folder = HomeFolder,
                    date = date,
                    section_name = section_name,
                    analysis_name = subsection_name)
print("settings:")
print(params_list)
if(!dir.exists(paste0(plotsdir, "/", subsection_name, "/"))) dir.create(paste0(plotsdir, "/", subsection_name, "/"), recursive = T)

complexes_cor <- filter(target_target_cor_anno, 
                             !(gene_1 %in% c(off_target$target_gene, off_target$off_target_gene)) & !(gene_2 %in% c(off_target$target_gene, off_target$off_target_gene)) & 
                               gene_1 < gene_2)

rmarkdown::render(paste0(CodeFolder, "/", ExperimentName, "/", "pipeline", "/", subsection_name, ".Rmd"),
                  output_file = file.path(plotsdir, paste0( date, "_", subsection_name, ".html")),
                  params =  params_list)

## ---- Co-essentiality Modules

## this isn't working and i'm not sure that i care

if (F){
  analysis_name <- '11c_02_deconvoluting_signal_from_glm_coessentiality_modules'
  params_list <- list(home_folder = HomeFolder,
                      date = date,
                      section_name = section_name,
                      analysis_name = analysis_name)
  print("settings:")
  print(params_list)
  if(!dir.exists(paste0(plotsdir, "/", analysis_name, "/"))) dir.create(paste0(plotsdir, "/", analysis_name, "/"), recursive = T)
  
  coessentiality_cor <- filter(target_target_cor_anno, 
                               !(gene_1 %in% c(off_target$target_gene, off_target$off_target_gene)) &
                                 !(gene_2 %in% c(off_target$target_gene, off_target$off_target_gene)) & 
                                 gene_1 < gene_2 & is_protein_complex_pair == F)
  
  rmarkdown::render(paste0(CodeFolder, "/", ExperimentName, "/", "pipeline", "/", analysis_name, ".Rmd"),
                    output_file = file.path(plotsdir, paste0(date, "_", analysis_name, ".html")),
                    params =  params_list)
}


## ---- Highly Correlated Pairs
#find all highly correlated pairs

analysis_name <- "11c_03_deconvoluting_signal_from_glm_get_highly_correlated_pairs"

## all high conf correlations



params_list <- list(home_folder = HomeFolder,
                    date = date,
                    section_name = section_name,
                    analysis_name = analysis_name)
print("settings:")
print(params_list)
if(!dir.exists(paste0(plotsdir, "/", analysis_name, "/"))) dir.create(paste0(plotsdir, "/", analysis_name, "/"), recursive = T)

rmarkdown::render(paste0(CodeFolder, "/", ExperimentName, "/", "pipeline", "/", analysis_name, ".Rmd"),
                  output_file = file.path(plotsdir, paste0( date, "_", analysis_name, ".html")),
                  params =  params_list)

## plot everything
high_cor_conf <- fread(file.path(OutFolder, section_name, paste0(date, "_all_high_correlations.tsv")))
pdf(paste0(plotsdir, "/", date, "_", "all_high_cor.pdf"), width = 4, height = 4)
for (i in 1:dim(high_cor_conf)[1]){
  t1 <- high_cor_conf$gene_1[i]; t2 <- high_cor_conf$gene_2[i]
  df4plot <- data.frame(
    tg1 = magpie_lfcs_split[[t1]]$lfc,
    tg2 = magpie_lfcs_split[[t2]]$lfc,
    gene = magpie_lfcs_split[[t1]]$downstream_gene_name
  ) %>%
    mutate(is_de = ifelse(
      magpie_lfcs_split[[t2]]$pval_adj < sig_pval_thresh & magpie_lfcs_split[[t1]]$pval_adj < sig_pval_thresh, "both", 
      ifelse(
        magpie_lfcs_split[[t2]]$pval_adj < sig_pval_thresh | magpie_lfcs_split[[t1]]$pval_adj < sig_pval_thresh, "one", "not"
      )
    ))
  p <- ggplot(df4plot, aes(x = tg1, y = tg2, 
                           col = is_de, size = is_de, alpha = is_de)) + 
    xlab(t1) + ylab(t2) + 
    scale_color_manual("DEG", values = c("both" = 'darkred', 'one' = 'red', 'not' = 'gray')) + 
    scale_alpha_manual("DEG", values = c("both" = 1, 'one' = 1, 'not' = 0.4)) + 
    scale_size_manual("DEG", values = c("both" = 1, 'one' = 1, 'not' = 0.5)) + 
    geom_point() + annotate("text", x = min(df4plot$tg1)*0.9, y = max(df4plot$tg2)*0.8, 
                            hjust = 0,
                            label = paste0("r_all = ", round(high_cor_conf$cor_all[i], 4), ", p-val = ", round(high_cor_conf$cor_all_pval_adj[i], 4), "\n",
                                           "r_deg = ", round(high_cor_conf$cor_all_deg_excl_target[i], 4), ", p-val = ", round(high_cor_conf$cor_all_deg_pval_adj[i], 4), "\n",
                                           high_cor_conf$common_deg[i], " common DEGs")) + 
    theme_bw() + theme(legend.position = 'none')
  print(p)
}
dev.off()



### ----- scratch

