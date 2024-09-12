date <- "2022-08-15"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
ExperimentName <- "Magpie"
section_name <- "10b_number_of_perturbing_knockdowns"



library(Seurat)
library(tidyverse)
library(data.table)
library(reshape2)
library(ggpointdensity)
library(magrittr)


# options
args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
  if (arg == "--rmd_file"){ rmd_file <- args[[ind + 1]] }
}

# i/o
source(file.path(HomeFolder, paste0("scripts/io/", ExperimentName, "_io.R")))
source(file.path(HomeFolder, "scripts/io/Magpie_Utils.R"))
outdir <- file.path(OutFolder, section_name)
plotsdir <- file.path(HTMLFolder, "pipeline/", section_name)
if(!dir.exists(outdir)) dir.create(outdir)
if(!dir.exists(plotsdir)) dir.create(plotsdir)

# define parameters for Rmd
params_list <- list(utils_path = file.path(HomeFolder, "scripts/io/Magpie_Utils.R"),
                    io_path = file.path(HomeFolder, "scripts/io/Magpie_io.R"),
                    date = date,
                    section_name = section_name)
print(params_list)

## ---- LoadData

crispr_sequences <- gsub(list.files(file.path(ResourcesFolder, 'reference_sequences', 'added_sequences'), pattern = '[.]gtf'), pattern = '.gtf', replacement = '')
magpie_lfcs <- fread(file.path(OutFolder, "6b_calc_lfcs_transcriptome_wide_by_gene", paste0(date, "_combined_with_adj_pval.tsv.gz")))
magpie_lfcs <- magpie_lfcs %>%
  dplyr::filter(!(downstream_gene_name %in% crispr_sequences))
#off_targets <- fread("/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Magpie/Preprocess_External_Datasets/off_target/off_target_pairs_magpie-2022-08-15-version.tsv.gz",
#                     header = T) %>% dplyr::filter(target_gene != off_target_gene)
#off_target_genes <- sort(unique(c(off_targets$target_gene, off_targets$off_target_gene)))


## ---- DataPreprocessing

downstream_meta <- magpie_lfcs %>%
  dplyr::filter(target != downstream_gene_name) %>%
  group_by(downstream_gene_name) %>%
  summarize(n_targets_that_de = sum(pval_adj < sig_pval_thresh & abs(lfc) > sig_abs_lfc_thresh),
            n_targets_that_upregulate = sum(pval_adj < sig_pval_thresh & lfc > sig_abs_lfc_thresh),
            n_targets_that_downregulate = sum(pval_adj < sig_pval_thresh & lfc < -sig_abs_lfc_thresh),
            mean_downstream_expr = control_norm_expr[1])

## add in variance of expression
sc_variance_fnm <- sort(grep(list.files(file.path(OutFolder, "Preprocess_External_Datasets", "coexpression"), pattern = "sc_regressed_variance_per_gene"), pattern = date, value = T), decreasing = T)[1]
sc_variance <- as.data.frame(fread(file.path(OutFolder, "Preprocess_External_Datasets", "coexpression", sc_variance_fnm)))
names(sc_variance) <- c("gene", "sc_variance")
sc_variance$gene <- unlist(lapply(strsplit(sc_variance$gene, ":"), "[[", 2))
downstream_meta$sc_variance <- sc_variance[,2][match(downstream_meta$downstream_gene_name, sc_variance[,1])]

## # of coexpressed genes
coexpression_fnm <- sort(grep(list.files(file.path(OutFolder, "Preprocess_External_Datasets", "coexpression"), pattern = 'sc_regressed_correlation_pairs', full.names = T), pattern = date, value = T), decreasing = T)[1]
coexpression_df <- fread(coexpression_fnm)
n_coexpressed_genes <- coexpression_df %>%
  dplyr::filter(abs(correlation) > sig_cor_thresh_downstream)
n_coexpressed_genes <- as.data.frame(table(c(n_coexpressed_genes$gene_1, n_coexpressed_genes$gene_2)))
downstream_meta <- downstream_meta %>%
  mutate(n_coexpressed = n_coexpressed_genes$Freq[match(downstream_gene_name, n_coexpressed_genes$Var1)]) %>%
  mutate(n_coexpressed = ifelse(is.na(n_coexpressed), 0, n_coexpressed))


## add heritability of expression
h <- read.csv(file.path(ResourcesFolder, "heritability/expr_variance_expl_by_donor.csv"))
downstream_meta <- left_join(downstream_meta, h, by = c("downstream_gene_name" = 'gene'))

## conservation score
cons_scores <- read.csv(file.path(ResourcesFolder,
                                  "conservation_scores",
                                  "conservation_scores_phastCons100way.UCSC.hg38.csv"))
downstream_meta <- left_join(downstream_meta, cons_scores, by = c("downstream_gene_name" = "GeneSymbol"))


## dorothea, msigdb, eQTLs
other_data_sources <- c("dorothea", "msigdb_gene_sets", "eQTLs")

for (other_data_source in other_data_sources){
  fnm <- sort(grep(list.files(file.path(OutFolder, "Preprocess_External_Datasets", other_data_source), pattern = paste0("annotated-guides")), pattern = paste0(date, "-version.tsv.gz"), value = T), decreasing = T)[1]
  annotations <- fread(file.path(OutFolder, "Preprocess_External_Datasets", other_data_source, fnm))
  if (other_data_source == "dorothea"){
    downstream_meta$num_dorothea_source <- annotations$num_dorothea_source[match(downstream_meta$downstream_gene_name, annotations$gene)]
    downstream_meta$num_dorothea_target <- annotations$num_dorothea_target[match(downstream_meta$downstream_gene_name, annotations$gene)]
    downstream_meta$num_dorothea_genes <- annotations$num_dorothea_genes[match(downstream_meta$downstream_gene_name, annotations$gene)]
  } else if (other_data_source == "msigdb_gene_sets"){
    downstream_meta$num_msigdb_gene_sets <- annotations$num_hallmark_gene_sets[match(downstream_meta$downstream_gene_name, annotations$gene)]
  } else if (other_data_source == "eQTLs"){
    downstream_meta$num_eQTLs_trans_gene <- annotations$num_trans_gene_pairs[match(downstream_meta$downstream_gene_name, annotations$gene)]
    downstream_meta$num_eQTLs_cis_gene <- annotations$num_cis_gene_pairs[match(downstream_meta$downstream_gene_name, annotations$gene)]
    downstream_meta$num_eQTLs_genes <- annotations$num_eQTL_gene_pairs[match(downstream_meta$downstream_gene_name, annotations$gene)]
  }
}

## protein complexes (omnipath)

#interactions
omnipath_complexes_fnm <- sort(grep(list.files(file.path(OutFolder, "Preprocess_External_Datasets", "omnipath"), pattern = paste0("omnipath-interactions_annotated-guides-")), pattern = paste0(date, "-version.tsv.gz"), value = T), decreasing = T)[1]
omnipath_complexes <- fread(file.path(OutFolder, "Preprocess_External_Datasets", "omnipath", omnipath_complexes_fnm))
downstream_meta$num_omnipath_interactions <- omnipath_complexes$num_omnipath_interactions[match(downstream_meta$downstream_gene_name, omnipath_complexes$gene)]

#complexes
omnipath_complexes_fnm <- sort(grep(list.files(file.path(OutFolder, "Preprocess_External_Datasets", "omnipath"), pattern = paste0("protein-complexes_annotated-guides-")), pattern = paste0(date, "-version.tsv.gz"), value = T), decreasing = T)[1]
omnipath_complexes <- fread(file.path(OutFolder, "Preprocess_External_Datasets", "omnipath", omnipath_complexes_fnm))
downstream_meta$num_omnipath_protein_complexes <- omnipath_complexes$num_omnipath_interactions[match(downstream_meta$downstream_gene_name, omnipath_complexes$gene)]


## add number of gprofiler terms
gprofiler_annotations <- readRDS("/lustre/scratch123/hgi/projects/crispri_scrnaseq/resources/gprofiler/gprofiler_full_hsapiens.ENSG.RDS")
gprofiler_terms <- fread("/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Magpie/Preprocess_External_Datasets/gprofiler/go_terms_2023-10-18.tsv")
num_go_anno <- table(unlist(gprofiler_annotations[grep(names(gprofiler_annotations), pattern = "GO:")]))
downstream_meta <- downstream_meta %>% 
  mutate(num_go_terms = as.numeric(num_go_anno)[match(downstream_gene_name, names(num_go_anno))])%>%
  mutate(num_go_terms = ifelse(is.na(num_go_terms), 0, num_go_terms))

fwrite(downstream_meta, paste0(outdir, '/', date, "_downstream_meta_data.csv"))

## ----

crispr_sequences <- gsub(list.files(file.path(ResourcesFolder, 'reference_sequences', 'added_sequences'), pattern = '[.]gtf'), pattern = '.gtf', replacement = '')
downstream_meta <- fread(paste0(outdir, '/', date, "_downstream_meta_data.csv"))
downstream_meta <- downstream_meta %>%
  dplyr::filter(!(downstream_gene_name %in% crispr_marker_sequences))

# run Rmd and create htmls
fnm_html <- paste0(section_name, ".html")
rmarkdown::render(paste0(CodeFolder, "Magpie/pipeline/", section_name, ".Rmd"),
                  output_file = file.path(plotsdir, fnm_html),
                  params =  params_list)
