library(data.table)
library(dplyr)
library(reshape2)
library(ggplot2)
library(tidyverse)


# set relevent i/o paths
section_name <- "13a_get_pluripotency_genes"
#correlation within the complex

# set relevent i/o paths
date <- "2022-08-15"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
ExperimentName <- "Magpie"

args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  ##needed arguments
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--experiment_name"){ ExperimentName <- args[[ind + 1]] }
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
}

setwd(HomeFolder)
source(paste0(HomeFolder, "scripts/io/", ExperimentName, "_io.R"))
source(paste0(HomeFolder, "/scripts/io/Magpie_Utils.R"))

print(paste("date:", date))
print(paste("home folder:", HomeFolder))
print(paste("Experiment Name:", ExperimentName))


## set io
outdir <- file.path(OutFolder, section_name)
if(!dir.exists(outdir)) dir.create(outdir, recursive = T)
plotsdir <- file.path(HTMLFolder, "pipeline", section_name)
if(!dir.exists(plotsdir)) dir.create(plotsdir, recursive = T)

## ---- LoadData

guide_meta <- fread(paste0(FilesFolder, ExperimentName, "Library.tsv"))
guide_meta <- guide_meta[!(duplicated(guide_meta$gene)),] %>%
  dplyr::select(c('gene', 'group', contains("pre_day")))

downstream_meta <- fread(paste0(OutFolder, "10b_number_of_perturbing_knockdowns/", date, "_downstream_meta_data.csv"))
target_meta <- fread(paste0(OutFolder, "7b_target_summary/", date, "_target_meta_data.csv"))
magpie_lfcs <- fread(paste0(OutFolder, "6b_calc_lfcs_transcriptome_wide_by_gene/", date, "_combined_with_adj_pval.tsv.gz"))
downstream_downstream_cor_anno <- fread(paste0(OutFolder, "12a_downstream_gene_embedding/", date, "_downstream_downstream_cor_with_adj_pval.tsv.gz"))
target_target_cor_anno <-fread(paste0(OutFolder, "11a_target_target_cor/", date, "_target_target_cor.tsv.gz"))

pluripotency_clusters_df <- fread('/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Magpie/17a_choose_pluripotency_genes/2022-08-15_tf_clusters.tsv')
pluripotency_lfc_fnms <- lapply(1:2,
                                FUN = function(i){
                                  fread(paste0('/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Magpie/17c_calc_lfcs_transcriptome_wide_by_pluripotency_cluster/by_gene/', i, "_", "2022-08-15.tsv.gz"))
                                })
pluripotency_lfcs <- pluripotency_lfc_fnms %>%
  bind_rows() %>%
  mutate(downstream_gene_name = unlist(lapply(strsplit(downstream_gene, ':'), "[[", 2))) %>%
  dplyr::select(c('genes_in_complex', 'downstream_gene_name', 'lfc', 'pval_lm'))
pluripotency_lfcs$pval_adj <- p.adjust(pluripotency_lfcs$pval_lm, method = 'BH')


gprofiler_annotations <- readRDS(paste0(ProjectFolder, "resources//gprofiler/gprofiler_full_hsapiens.ENSG.RDS"))
fnm <- sort(list.files(paste0(OutFolder, "Preprocess_External_Datasets/gprofiler/"), pattern = 'go_terms', full.names = T), decreasing = T)[1]
gprofiler_terms <- fread(fnm)
fnm <- sort(grep(list.files(paste0(OutFolder, "Preprocess_External_Datasets/omnipath/"), pattern = 'omnipath_complex_meta', full.names = T), pattern = date, value = T), decreasing = T)[1]
omnipath_protein_complexes <- fread(fnm)
off_targets <- fread(paste0(OutFolder, "Preprocess_External_Datasets/off_target/off_target_pairs_magpie-", date, "-version.tsv.gz"), header = T) %>%
  dplyr::filter(target_gene != off_target_gene)
targets2exclude <- unique(c(off_targets$target_gene, off_targets$off_target_gene))
fnm <- sort(grep(list.files(paste0(OutFolder, "Preprocess_External_Datasets/dorothea/"), pattern = 'dorothea_interaction_pairs', full.names = T), pattern = date, value = T), decreasing = T)[1]
dorothea_tfs <- fread(fnm)
msigdb_hallmark_pathways <- readRDS(paste0(ProjectFolder, "resources/gene_sets/h.all.v7.4.symbols.RDS"))



## ---- Run R Markdown

fnm_html <- paste0(section_name, ".html")
rmarkdown::render(rmd_file,
                  envir = new.env(),
                  output_file = file.path(plotsdir, fnm_html),
                  params =  params_list)


## subset to just the genes with lots of signal
targets2consider <- magpie_lfcs %>%
  dplyr::filter(pval_adj < 0.1 & abs(lfc) > 0.25 & target != downstream_gene_name) %>%
  group_by(target) %>%
  summarize(n_deg = n())%>%
  dplyr::filter(n_deg >= 10) %>%
  .$target
downstreams2consider <- magpie_lfcs %>%
  dplyr::filter(pval_adj < 0.1 & abs(lfc) > 0.25 & target != downstream_gene_name) %>%
  group_by(downstream_gene_name) %>%
  summarize(n_reg = n())%>%
  dplyr::filter(n_reg >= 5)
df4heatmap <- magpie_lfcs %>%
  dplyr::filter(target %in% targets2consider$target  & 
                  downstream_gene_name %in% downstreams2consider$downstream_gene_name ) %>%
  dplyr::select(c('target', 'downstream_gene_name', 'lfc')) %>%
  reshape2::dcast(target ~ downstream_gene_name, value.var = 'lfc') %>%
  column_to_rownames('target') %>%
  as.matrix()
my_heatmap(df4heatmap,
           treeheight_row=0, treeheight_col = 0,
           min_c = -0.5, max_c = 0.5,
           border_col = NA)

## pluripotency genes
col_anno <- data.frame(
  gene = downstreams2consider$downstream_gene_name,
  ipsc_gene = factor(downstreams2consider$downstream_gene_name %in% ips_genes)
) %>%
  column_to_rownames('gene')
my_heatmap(df4heatmap,
           annotation_col = col_anno,
           treeheight_row=0, treeheight_col = 0,
           min_c = -0.5, max_c = 0.5,
           border_col = NA)

ips_genes <- unique(c(gprofiler_annotations[['GO:0048863']], gprofiler_annotations[['GO:0019827']], gprofiler_annotations[['GO:0048864']]))

## cluster things
downstream_clust <- hclust(dist(df4heatmap))
## cluster the genes
fviz_nbclust(df4heatmap, kmeans, method = "silhouette", k.max = 20)+
  labs(subtitle = "Silhouette method")
cluster_res <- hclust(dist(t(df4heatmap)))







## ---- Pick Genes

targets2consider <- target_target_cor_anno %>%
  dplyr::filter(gene_1 != gene_2 & cor_all > 0.4) %>%
  group_by(gene_1) %>%
  summarize(n_high_correlates = n()) %>%
  dplyr::filter(n_high_correlates >= 3) %>% .$gene_1
downstreams2consider <- magpie_lfcs %>%
  dplyr::filter(target %in% targets2consider & pval_adj < 0.1 & abs(lfc) > 0.25 & target != downstream_gene_name) %>%
  .$downstream_gene_name %>% 
  table() %>% as.data.frame() %>%
  dplyr::filter(Freq >= 3) %>%
  .$`.`

df4heatmap <- magpie_lfcs %>%
  dplyr::filter(target %in% targets2consider  & 
                  downstream_gene_name %in% downstreams2consider) %>%
  dplyr::select(c('target', 'downstream_gene_name', 'lfc', 'pval_adj')) %>%
  #mutate(lfc = ifelse(pval_adj < 0.1, lfc, 0)) %>%
  reshape2::dcast(target ~ downstream_gene_name, value.var = 'lfc') %>%
  column_to_rownames('target') %>%
  as.matrix()

row_anno <- target_meta %>%
  dplyr::filter(gene %in% row.names(df4heatmap)) %>%
  dplyr::select(c("gene", "pre_day10_lfc_mean")) %>%
  column_to_rownames('gene')
my_heatmap(df4heatmap,
           treeheight_row=0, treeheight_col = 0,
           min_c = -0.5, max_c = 0.5,
           border_col = NA, 
           annotation_row=row_anno)


## stick to the umap genes
targets2consider <- target_umap_coords$gene
downstreams2consider <- downstream_umap_coords$gene

df4heatmap <- magpie_lfcs %>%
  dplyr::filter(target %in% targets2consider  & 
                  downstream_gene_name %in% downstreams2consider) %>%
  dplyr::select(c('target', 'downstream_gene_name', 'lfc', 'pval_adj')) %>%
  #mutate(lfc = ifelse(pval_adj < 0.1, lfc, 0)) %>%
  reshape2::dcast(target ~ downstream_gene_name, value.var = 'lfc') %>%
  column_to_rownames('target') %>%
  as.matrix()

## cluster rows and columns
downstream_clust <- hclust(dist(df4heatmap))
## cluster the genes
fviz_nbclust(df4heatmap, kmeans, method = "silhouette", k.max = 20)+
  labs(subtitle = "Silhouette method")
cluster_res <- hclust(dist(t(df4heatmap)))
col_anno <- data.frame(
  gene = colnames(df4heatmap),
  cluster = factor(as.numeric(cutree(cluster_res, k = 16)))
) %>%column_to_rownames('gene')


gene_anno <- data.frame(
  gene = row.names(df4heatmap),
  annotation = ifelse(row.names(df4heatmap) %in% unlist(strsplit(omnipath_protein_complexes %>% dplyr::filter( complex_name ==  "Mediator complex") %>% .$genes_in_complex, "_")), "Mediator complex", 
                      ifelse(row.names(df4heatmap) %in% unlist(strsplit(omnipath_protein_complexes %>% dplyr::filter( complex_name %in% c("TFIID complex","TFIIF","TFIIA", "TFIIE complex")) %>% .$genes_in_complex, "_")), "Pre-initiation complexes", 
                             ifelse(row.names(df4heatmap) %in% unlist(strsplit(omnipath_protein_complexes %>% dplyr::filter( complex_name ==  "RNA polymerase II core complex") %>% .$genes_in_complex, "_")), "RNA polymerase II core complex", 
                                    ifelse(row.names(df4heatmap) %in% unlist(strsplit(omnipath_protein_complexes %>% dplyr::filter( complex_name ==  "Spliceosome") %>% .$genes_in_complex, "_")), "Spliceosome", 
                                           ifelse(row.names(df4heatmap) %in% unlist(strsplit(omnipath_protein_complexes %>% dplyr::filter( complex_name ==  "CCR4-NOT complex") %>% .$genes_in_complex, "_")), "CCR4-NOT complex", 
                                                  ifelse(row.names(df4heatmap) %in% unlist(strsplit(omnipath_protein_complexes %>% dplyr::filter( complex_name ==  "INO80 chromatin remodeling complex") %>% .$genes_in_complex, "_")), "INO80 chromatin remodeling complex", 
                                                         ifelse(row.names(df4heatmap) %in% unlist(strsplit(omnipath_protein_complexes %>% dplyr::filter( complex_name ==  "NuA4/Tip60-HAT complex") %>% .$genes_in_complex, "_")), "NuA4/Tip60-HAT complex", 
                                                                ifelse(row.names(df4heatmap) %in% unlist(strsplit(omnipath_protein_complexes %>% dplyr::filter( complex_name ==  "Paf complex") %>% .$genes_in_complex, "_")), "Paf complex", 
                                                                       ifelse(row.names(df4heatmap) %in% gprofiler_annotations[['GO:0036396']], "RNA N6-methyladenosine methyltransferase complex",
                                                                              ifelse(row.names(df4heatmap) %in% unlist(strsplit(omnipath_protein_complexes %>% dplyr::filter( complex_name ==  "EIF2B2-EIF2B3-EIF2B4-EIF2B5 complex") %>% .$genes_in_complex, "_")), "EIF2 complex", 
                                                                                     ifelse(row.names(df4heatmap) %in% gprofiler_annotations[['GO:0006418']], "tRNA aminoacylation for protein translation", 
                                                                                            ifelse(row.names(df4heatmap) %in% unlist(strsplit(omnipath_protein_complexes %>% dplyr::filter( complex_name ==  "EIF3_complex") %>% .$genes_in_complex, "_")), "EIF3 complex", 
                                                                                                   ' ' )
                                                                                     )
                                                                              )
                                                                       )
                                                                )
                                                         )
                                                  )
                                           )
                                    ))))
)
row_anno <- gene_anno %>% 
  mutate(essentiality = target_meta$pre_day10_lfc_mean[match(gene, target_meta$gene)]) %>%
  column_to_rownames('gene')
row_anno$annotation <- factor(row_anno$annotation)




## complex colors
histone_complexes <- c("CCR4-NOT complex", "INO80 chromatin remodeling complex", "NuA4/Tip60-HAT complex", "RNA N6-methyladenosine methyltransferase complex")
histone_colors <- unlist(lapply(1:length(histone_complexes), FUN = function(i){hue_pal(h = c(-15, 15), l = i*80/length(histone_complexes))(1)})); names(histone_colors) <- histone_complexes
transcription_complexes <- c("Pre-initiation complexes","Mediator complex" , "RNA polymerase II core complex", "Paf complex")
transcription_colors <- unlist(lapply(1:length(transcription_complexes), FUN = function(i){hue_pal(h =90 + c(-15, 15), l = i*80/length(transcription_complexes))(1)})); names(transcription_colors) <- transcription_complexes
splicing_complexes <- c("Spliceosome")
splicing_colors <- hue_pal(h = 180 + c(-15, 15), l = 40)(length(splicing_complexes)); names(splicing_colors) <- splicing_complexes
translation_complexes <- c("EIF2 complex", 'tRNA aminoacylation for protein translation', "EIF3 complex")
translation_colors <- unlist(lapply(1:length(translation_complexes), FUN = function(i){hue_pal(h =270 + c(-15, 15), l = i*80/length(translation_complexes))(1)})); names(transcription_colors) <- transcription_complexes; names(translation_colors) <- translation_complexes
complex_colors <- c(histone_colors, transcription_colors, translation_colors, splicing_colors, c(' ' = NA))


my_heatmap(df4heatmap,
           treeheight_row=0, treeheight_col = 0,
           min_c = -0.25, max_c = 0.25,
           border_col = NA, 
           annotation_row = row_anno,
           annotation_colors = list(' ' = complex_colors ),
           annotation_col = col_anno)
## ---- try again (based on the aggregated things)

downstreams2consider <- downstream_umap_coords$gene
targets2consider <- magpie_lfcs %>% 
  dplyr::filter(downstream_gene_name %in% downstreams2consider &
                  pval_adj < 0.1 & abs(lfc) > 0.4 &
                  target != downstream_gene_name) %>%
  .$target %>%
  table() %>% as.data.frame() %>%
  dplyr::filter(Freq >= 3) %>%
  .$`.`

df4heatmap <- magpie_lfcs %>%
  dplyr::filter(target %in% targets2consider  & 
                  downstream_gene_name %in% downstreams2consider) %>%
  dplyr::select(c('target', 'downstream_gene_name', 'lfc', 'pval_adj')) %>%
  #mutate(lfc = ifelse(pval_adj < 0.1, lfc, 0)) %>%
  reshape2::dcast(target ~ downstream_gene_name, value.var = 'lfc') %>%
  column_to_rownames('target') %>%
  as.matrix()

row_anno <- target_meta %>%
  dplyr::filter(gene %in% row.names(df4heatmap)) %>%
  dplyr::select(c("gene", "pre_day10_lfc_mean")) %>%
  column_to_rownames('gene')

pdf('scratch/2024-02-06_heatmap.pdf', width = 30, height = 20)
my_heatmap(df4heatmap,
           treeheight_row=0, treeheight_col = 0,
           min_c = -0.5, max_c = 0.5,
           border_col = NA, 
           annotation_row = row_anno)#,
           #annotation_colors = list(' ' = complex_colors ))
dev.off()




## possible ipsc genes
sort(gprofiler_annotations[['GO:0072089']])
sort(gprofiler_annotations[['GO:0048863']])
sort(gprofiler_annotations[['GO:0048864']])
sort(gprofiler_annotations[['GO:0001824']])
sort(gprofiler_annotations[['GO:0090263']])
