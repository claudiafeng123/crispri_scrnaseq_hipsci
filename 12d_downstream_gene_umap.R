library(data.table)
library(dplyr)
library(reshape2)
library(ggplot2)
library(umap)
library(scales)
library(ggrepel)
library(Seurat)

# set relevent i/o paths
section_name <- "12d_downstream_gene_umap"
#correlation within the complex

# set relevent i/o paths
date <- "2022-08-15"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
ExperimentName <- "Magpie"
#REDO <- T
min_deg <- 1 #min_deg is min_deg in trans and strict
min_target <- 1
min_cells <- 10

args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  ##needed arguments
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--experiment_name"){ ExperimentName <- args[[ind + 1]] }
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
  
  ## parameters
  if (arg == "--min_deg"){ min_deg <- as.numeric(args[[ind + 1]]) }
  if (arg == "--min_target"){ min_target <- as.numeric(args[[ind + 1]]) }
  if (arg == "--min_cells"){min_cells <- as.numeric(args[[ind + 1]]) }
  
  
  ### not needed
  if (arg == "--utils_path"){ utils_path <- args[[ind + 1]] }
  if (arg == "--io_path"){ io_path <- args[[ind + 1]] }
  if (arg == "--outdir"){ outdir1 <- args[[ind + 1]] }
  if (arg == "--REDO"){ REDO <- as.logical(args[[ind + 1]]) }
  
}

setwd(HomeFolder)
if (!(exists("io_path"))){io_path <- paste0(HomeFolder, "scripts/io/", ExperimentName, "_io.R")}
if (!(exists("utils_path"))){ utils_path <- paste0(HomeFolder, "/scripts/io/Magpie_Utils.R")}
source(io_path)
source(utils_path)
control_tag <- c("unassigned", "NonTarget") # use all unassigned cells as control

print(paste("date:", date))
print(paste("home folder:", HomeFolder))
print(paste("Experiment Name:", ExperimentName))
print(paste("min_cells_per_gene:", min_cells_per_gene))
print(paste("section_name:", section_name))


## set io
#
config <- paste0(date, "_min_deg-", min_deg, "_min_targets-", min_target, "_min_cells-", min_cells)
outdir <- file.path(OutFolder, section_name)
if(!dir.exists(outdir)) dir.create(outdir, recursive = T)
plotsdir <- file.path(HTMLFolder, "pipeline", section_name)
if(!dir.exists(plotsdir)) dir.create(plotsdir, recursive = T)

#if (min_deg_per_target <= 10){min_deg_per_target <- 10}


## ---- Load Data

#complex_meta <- fread(file.path(OutFolder, '13a_calc_lfcs_transcriptome_wide_by_complex', paste0(date, '_complex_meta.tsv')))
target_meta <- fread(file.path(OutFolder, "7b_target_summary", paste0(date, "_target_meta_data.csv")))
umap_coords <- fread(file.path(OutFolder, "12a_downstream_gene_embedding", config, paste0(date, "_umap_coords.tsv")))

#magpie_lfcs_by_downstream <- split.data.frame(magpie_lfcs_by_downstream, f = magpie_lfcs_by_downstream$downstream_gene_name)
#all_downstream_genes <- names(magpie_lfcs_by_downstream)


## go terms
gprofiler_terms <- fread(file.path(OutFolder, "Preprocess_External_Datasets", "gprofiler", "go_terms_2022-09-17.tsv"))
gprofiler_annotations <- readRDS(file.path(ResourcesFolder, "gprofiler", "gprofiler_full_hsapiens.ENSG.RDS"))

msigdb_pathways <- readRDS(file.path(ResourcesFolder, "gene_sets", "h.all.v7.4.symbols.RDS"))

off_target <- fread(paste0(OutFolder, "/Preprocess_External_Datasets/off_target/off_target_pairs_magpie-", date, "-version.tsv.gz"), header = T) %>%
  filter(target_gene != off_target_gene)

tf_fnm <- sort(grep(list.files(paste0(OutFolder, "/Preprocess_External_Datasets/dorothea/"), pattern = date), pattern =  "interaction_pairs", value = T), decreasing = T)[1]
tfs <- fread(paste0(OutFolder, "/Preprocess_External_Datasets/dorothea/", tf_fnm))
#tfs <- unique(tfs$source)

protein_ixn_fnm <- sort(grep(list.files(paste0(OutFolder, "/Preprocess_External_Datasets/omnipath/"), pattern = date), pattern =  "interaction_pairs", value = T), decreasing = T)[1]
protein_ixns <- fread(paste0(OutFolder, "/Preprocess_External_Datasets/omnipath/", protein_ixn_fnm))

gene_anno <- fread('/lustre/scratch123/hgi/projects/crispri_scrnaseq/resources/gene_sets/gene_complex_links_2023-03-13.tsv')
highly_correlated_modules <- fread("/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Magpie/11d_target_target_cor/2022-08-15_high_cor_anno.tsv")

#remove anything that might be an off-target
umap_coords <- umap_coords[which(!(umap_coords$gene %in% off_target$off_target_gene)), ]

## ---- Data Exploration - UMAP by complex/other manual annotations

data_exploration <- F
if (data_exploration == T){
  
  ## go terms
  pdf(paste0(plotsdir, "/", date, "_umap_by_go_term.pdf"), width = 6, height = 6)
  complexes2plot <- gprofiler_annotations[which(unlist(lapply(gprofiler_annotations, length)) > 25 & 
                                                  unlist(lapply(gprofiler_annotations, length)) < 100 )]
  #complexes2plot <- complexes2plot[1:10]
  for (i in 1:length(complexes2plot)){
    #complex_name <- complexes2plot[i]
    genes_of_interest <- complexes2plot[[i]]
    p <- ggplot(umap_coords, aes(x = umap_1, y = umap_2, 
                                 col = as.character(gene %in% genes_of_interest),
                                 label = ifelse(gene %in% genes_of_interest, gene , ""))) + 
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = 'gray')) + 
      geom_point() + ggrepel::geom_text_repel(max.overlaps = Inf) + 
      ggtitle(gprofiler_terms$Term[match(names(complexes2plot)[i], gprofiler_terms$GO_ID)]) + 
      theme_bw() + theme(legend.position = "none")
    print(p)
  }
  dev.off()
  
  ##by knockdown
  ## go terms
  targets2plot <- filter(target_meta, n_downstream_excl_target > 1)
  magpie_lfcs <- fread(file.path(OutFolder, "6b_calc_lfcs_transcriptome_wide_by_gene", paste0(date, "_combined_with_adj_pval.tsv.gz")))
  magpie_lfcs <- dplyr::filter(magpie_lfcs, 
                               downstream_gene_name %in% umap_coords$gene & target %in% targets2plot$gene)
  magpie_lfcs <- split.data.frame(magpie_lfcs, f = magpie_lfcs$target)
  pdf(paste0(plotsdir, "/", date, "_umap_by_target_knockdown.pdf"), width = 6, height = 6)
  #complexes2plot <- complexes2plot[1:10]
  for (i in 1:length(magpie_lfcs)){
    knockdown <- names(magpie_lfcs)[i]
    df4plot <- umap_coords %>%
      mutate(lfc = magpie_lfcs[[knockdown]]$lfc[match(gene, magpie_lfcs[[knockdown]]$downstream_gene_name)])
    p <- ggplot(df4plot, aes(x = umap_1, y = umap_2, 
                                 col = ifelse(abs(lfc) > 1, sign(lfc), lfc))) + 
      scale_colour_gradient2(low = "red",
                             mid = "white",
                             high = "darkgreen", ) + 
      xlab("UMAP_1") + ylab("UMAP_2") + 
      geom_point(size = 0.4) + 
      ggtitle(knockdown) + 
      theme_bw() + theme(legend.position = "none")
    print(p)
  }
  dev.off()
  
  
  ## thomas's list
  pdf(paste0(plotsdir, "/", date, "_umap_by_function_manual_annotated.pdf"), width = 6, height = 6)
  complexes2plot <- unique(gene_anno$complex)
  for (i in 1:length(complexes2plot)){
    complex_name <- complexes2plot[i]
    genes_of_interest <- gene_anno %>% filter(complex == complex_name) %>% .$gene
    p <- ggplot(umap_coords, aes(x = umap_1, y = umap_2, 
                                 col = as.character(gene %in% genes_of_interest),
                                 label = ifelse(umap_coords$gene %in% genes_of_interest, gene , ""))) + 
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = 'gray')) + 
      ggtitle(complex_name) + 
      geom_point() + geom_text_repel(max.overlaps = Inf) + 
      theme_bw() + theme(legend.position = "none")
    print(p)
  }
  dev.off()
  
  p <- ggplot(umap_coords, aes(x = umap_1, y = umap_2, label = gene)) + 
    geom_point(size = 0.4) + theme_bw()
  plotly::ggplotly(p)
}

## ----


## copy and paste
if (F){
  complex_name1 <- "DNA replication"
  genes_of_interest <- c("PRIM2", "POLO", "POLE4", "POLE2", "POLA1", "POLG")
  genes2label <- unique(c(genes2label, genes_of_interest))
  complex_of_interest <- filter(umap_coords, gene %in% genes_of_interest)
  label_info <- data.frame(
    gene = "",
    umap_1 = mean(complex_of_interest$umap_1),
    umap_2 = mean(complex_of_interest$umap_2),
    color = "gray",
    size = 1,
    transparency = '1',
    genes2label = paste(complex_of_interest$gene, collapse = "_"),
    label =  complex_name1,
    label_x = mean(complex_of_interest$umap_1) + 0.25,
    label_y = mean(complex_of_interest$umap_2)
  )
  label_df <- as.data.frame(bind_rows(label_info, label_df))
  segments2draw <- data.frame(
    label = complex_name1,
    x_start = mean(complex_of_interest$umap_1),
    x_end = complex_of_interest$umap_1,
    y_start = mean(complex_of_interest$umap_2),
    y_end = complex_of_interest$umap_2
  )
  segment_df <- as.data.frame(bind_rows(segments2draw, segment_df))
  
}


umap_by_complex <- T
if (umap_by_complex == T){
  ## label df
  label_df <- data.frame(
    gene = character(0),
    umap_1 = numeric(0),
    umap_2 = numeric(0),
    color = character(0),
    size = numeric(0),
    transparency = character(0),
    genes2label = character(0),
    label = character(0),
    label_x = numeric(0),
    label_y = numeric(0)
  )
  segment_df <- data.frame(
    label = character(0),
    x_start = numeric(0),
    x_end = numeric(0),
    y_start = numeric(0),
    y_end = numeric(0)
  )
  genes2label <- c()
  #genes_considered <- c()
  
  
  
  complex_name1 <- "Mitochondrial ribosome"
  genes_of_interest <- c("MRPL20", "MRPS34", "MRPS24", "MRPL48", "MRPL21", "MRPS2", "MRPL52", "MRPL28", "MRPL37")
  genes2label <- unique(c(genes2label, genes_of_interest))
  complex_of_interest <- filter(umap_coords, gene %in% genes_of_interest)
  label_info <- data.frame(
    gene = "",
    umap_1 = mean(complex_of_interest$umap_1),
    umap_2 = mean(complex_of_interest$umap_2),
    color = "gray",
    size = 1,
    transparency = '1',
    genes2label = paste(complex_of_interest$gene, collapse = "_"),
    label =  complex_name1,
    label_x = mean(complex_of_interest$umap_1),
    label_y = mean(complex_of_interest$umap_2)+ 0.15
  )
  label_df <- as.data.frame(bind_rows(label_info, label_df))
  
  
  complex_name1 <- "ADP/ATP translocator"
  genes_of_interest <- c("SLC25A6", "SLC25A5")
  enes2label <- unique(c(genes2label, genes_of_interest))
  complex_of_interest <- filter(umap_coords, gene %in% genes_of_interest)
  label_info <- data.frame(
    gene = "",
    umap_1 = mean(complex_of_interest$umap_1),
    umap_2 = mean(complex_of_interest$umap_2),
    color = "gray",
    size = 1,
    transparency = '1',
    genes2label = paste(complex_of_interest$gene, collapse = "_"),
    label =  complex_name1,
    label_x = mean(complex_of_interest$umap_1),
    label_y = mean(complex_of_interest$umap_2)+ 0.15
  )
  label_df <- as.data.frame(bind_rows(label_info, label_df))
  segments2draw <- data.frame(
    label = complex_name1,
    x_start = mean(complex_of_interest$umap_1),
    x_end = complex_of_interest$umap_1,
    y_start = mean(complex_of_interest$umap_2)+ 0.1,
    y_end = complex_of_interest$umap_2 
  )
  segment_df <- as.data.frame(bind_rows(segments2draw, segment_df))
  
  
  complex_name1 <- "TOMM complex"
  genes_of_interest <- grep(umap_coords$gene, pattern = "TOMM", value = T)
  genes_of_interest <- genes_of_interest[which(!(genes_of_interest %in% 
                                                   c("TOMM34")))]
  genes2label <- unique(c(genes2label, genes_of_interest))
  complex_of_interest <- filter(umap_coords, gene %in% genes_of_interest)
  label_info <- data.frame(
    gene = "",
    umap_1 = mean(complex_of_interest$umap_1),
    umap_2 = mean(complex_of_interest$umap_2),
    color = "gray",
    size = 1,
    transparency = '1',
    genes2label = paste(complex_of_interest$gene, collapse = "_"),
    label =  complex_name1,
    label_x = mean(complex_of_interest$umap_1),
    label_y = mean(complex_of_interest$umap_2)
  )
  label_df <- as.data.frame(bind_rows(label_info, label_df))
  
  
  complex_name1 <- "ETC complex I"
  genes_of_interest <- grep(umap_coords$gene, pattern = "MT-ND|NDUF", value = T)
  genes_of_interest <- genes_of_interest[which(!(genes_of_interest %in% 
                                                   c("NDUFA6", "NDUFS2", "NDUFA8", "NDUFS7", "NDUFAF3", "MT-ND6", "MT-ND5",
                                                     "NDUFV2", "NDUFB10", "NDUFAF8")))]
  genes2label <- unique(c(genes2label, genes_of_interest))
  complex_of_interest <- filter(umap_coords, gene %in% genes_of_interest)
  label_info <- data.frame(
    gene = "",
    umap_1 = mean(complex_of_interest$umap_1),
    umap_2 = mean(complex_of_interest$umap_2),
    color = "gray",
    size = 1,
    transparency = '1',
    genes2label = paste(complex_of_interest$gene, collapse = "_"),
    label =  complex_name1,
    label_x = mean(complex_of_interest$umap_1),
    label_y = mean(complex_of_interest$umap_2) - 0.05
  )
  label_df <- as.data.frame(bind_rows(label_info, label_df))
  segments2draw <- data.frame(
    label = complex_name1,
    x_start = mean(complex_of_interest$umap_1),
    x_end = complex_of_interest$umap_1,
    y_start = mean(complex_of_interest$umap_2),
    y_end = complex_of_interest$umap_2
  )
  segment_df <- as.data.frame(bind_rows(segments2draw, segment_df))
  
  
  complex_name1 <- "ETC complex IV"
  genes_of_interest <- c("COX5B", "COX5A", "COX4I1", "COX7B", "COX8A", "COX6C", "COX7A2", "COX7C")
  genes2label <- unique(c(genes2label, genes_of_interest))
  complex_of_interest <- filter(umap_coords, gene %in% genes_of_interest)
  label_info <- data.frame(
    gene = "",
    umap_1 = mean(complex_of_interest$umap_1),
    umap_2 = mean(complex_of_interest$umap_2),
    color = "gray",
    size = 1,
    transparency = '1',
    genes2label = paste(complex_of_interest$gene, collapse = "_"),
    label =  complex_name1,
    label_x = mean(complex_of_interest$umap_1) + 0.25,
    label_y = mean(complex_of_interest$umap_2) + 0.05
  )
  label_df <- as.data.frame(bind_rows(label_info, label_df))
  
  
  complex_name1 <- "ATP synthase"
  genes_of_interest <- grep(umap_coords$gene, pattern = "ATP5", value = T)
  genes_of_interest <- genes_of_interest[which(!(genes_of_interest %in% c("ATP5IF1", "ATP5PB", "ATP5F1B")))]
  genes2label <- unique(c(genes2label, genes_of_interest))
  complex_of_interest <- filter(umap_coords, gene %in% genes_of_interest)
  label_info <- data.frame(
    gene = "",
    umap_1 = mean(complex_of_interest$umap_1),
    umap_2 = mean(complex_of_interest$umap_2),
    color = "gray",
    size = 1,
    transparency = '1',
    genes2label = paste(complex_of_interest$gene, collapse = "_"),
    label =  complex_name1,
    label_x = mean(complex_of_interest$umap_1) + 0.25,
    label_y = mean(complex_of_interest$umap_2) - 0.1
  )
  label_df <- as.data.frame(bind_rows(label_info, label_df))
  
  
  complex_name1 <- "Cytosolic ribosome"
  genes_of_interest <- grep(umap_coords$gene, pattern = "^RPL|^RPS", value = T)
  genes_of_interest <- genes_of_interest[which(!(genes_of_interest %in% c("RPS6KB2", "RPL39L", "RPS6KC1", "RPL7L1", "RPS27L", "RPL26L1", "RPL19")))]
  genes2label <- unique(c(genes2label, genes_of_interest))
  complex_of_interest <- filter(umap_coords, gene %in% genes_of_interest)
  label_info <- data.frame(
    gene = "",
    umap_1 = mean(complex_of_interest$umap_1),
    umap_2 = mean(complex_of_interest$umap_2),
    color = "gray",
    size = 1,
    transparency = '1',
    genes2label = paste(complex_of_interest$gene, collapse = "_"),
    label =  complex_name1,
    label_x = mean(complex_of_interest$umap_1) + 0.25,
    label_y = mean(complex_of_interest$umap_2) + 0.1
  )
  label_df <- as.data.frame(bind_rows(label_info, label_df))
  segments2draw <- data.frame(
    label = complex_name1,
    x_start = mean(complex_of_interest$umap_1),
    x_end = complex_of_interest$umap_1,
    y_start = mean(complex_of_interest$umap_2),
    y_end = complex_of_interest$umap_2
  )
  segment_df <- as.data.frame(bind_rows(segments2draw, segment_df))
  
  
  
  
  complex_name1 <- "Amino acid metabolic process"
  genes_of_interest <- c("EIF2S2", "YARS", "BCAT1", "MARS", 'SARS', "AARS", "SLC7A11", "SLC1A5", "CARS", "GARS", "PYCR1", "SESN2", 'SHMT2', 'EIF4EBP1', "SLC3A2", 'CHAC1',  'PSAT1')
  genes2label <- unique(c(genes2label, genes_of_interest))
  complex_of_interest <- filter(umap_coords, gene %in% genes_of_interest)
  label_info <- data.frame(
    gene = "",
    umap_1 = mean(complex_of_interest$umap_1),
    umap_2 = mean(complex_of_interest$umap_2),
    color = "gray",
    size = 1,
    transparency = '1',
    genes2label = paste(complex_of_interest$gene, collapse = "_"),
    label =  complex_name1,
    label_x = mean(complex_of_interest$umap_1) + 0.25,
    label_y = mean(complex_of_interest$umap_2)
  )
  label_df <- as.data.frame(bind_rows(label_info, label_df))
  segments2draw <- data.frame(
    label = complex_name1,
    x_start = mean(complex_of_interest$umap_1),
    x_end = complex_of_interest$umap_1,
    y_start = mean(complex_of_interest$umap_2),
    y_end = complex_of_interest$umap_2
  )
  segment_df <- as.data.frame(bind_rows(segments2draw, segment_df))
  
  
  complex_name1 <- "Apoptosis via p53"
  genes_of_interest <- c("RPS27L", "IKIP", "RRM2B", "TNFRSF10B", "BBC3", "MDM2", "BAX", "TM7SF3", "TIGAR", "TRAF4")
  genes2label <- unique(c(genes2label, genes_of_interest))
  complex_of_interest <- filter(umap_coords, gene %in% genes_of_interest)
  label_info <- data.frame(
    gene = "",
    umap_1 = mean(complex_of_interest$umap_1),
    umap_2 = mean(complex_of_interest$umap_2),
    color = "gray",
    size = 1,
    transparency = '1',
    genes2label = paste(complex_of_interest$gene, collapse = "_"),
    label =  complex_name1,
    label_x = mean(complex_of_interest$umap_1) + 0.25,
    label_y = mean(complex_of_interest$umap_2)
  )
  label_df <- as.data.frame(bind_rows(label_info, label_df))
  segments2draw <- data.frame(
    label = complex_name1,
    x_start = mean(complex_of_interest$umap_1),
    x_end = complex_of_interest$umap_1,
    y_start = mean(complex_of_interest$umap_2),
    y_end = complex_of_interest$umap_2
  )
  segment_df <- as.data.frame(bind_rows(segments2draw, segment_df))
  

  complex_name1 <- "Glycolysis"
  genes_of_interest <- c("ALDOA", "TPI1", "GPI", "PKM", "ENO1", "PGK1", "PGAM1", "LDHA")
  genes2label <- unique(c(genes2label, genes_of_interest))
  complex_of_interest <- filter(umap_coords, gene %in% genes_of_interest)
  label_info <- data.frame(
    gene = "",
    umap_1 = mean(complex_of_interest$umap_1),
    umap_2 = mean(complex_of_interest$umap_2),
    color = "gray",
    size = 1,
    transparency = '1',
    genes2label = paste(complex_of_interest$gene, collapse = "_"),
    label =  complex_name1,
    label_x = mean(complex_of_interest$umap_1) + 0.25,
    label_y = mean(complex_of_interest$umap_2)
  )
  label_df <- as.data.frame(bind_rows(label_info, label_df))
  segments2draw <- data.frame(
    label = complex_name1,
    x_start = mean(complex_of_interest$umap_1),
    x_end = complex_of_interest$umap_1,
    y_start = mean(complex_of_interest$umap_2),
    y_end = complex_of_interest$umap_2
  )
  segment_df <- as.data.frame(bind_rows(segments2draw, segment_df))
  
  
  complex_name1 <- "Cholesterol biosynthesis"
  genes_of_interest <- c("MSMO1", "MVD", "HMGCS", "FDFT", "ACLY", "FDPS", "HMGCP", "CYP51A1", "INSIG1", "LSS", "DHCR7", "HSD17B7")
  genes2label <- unique(c(genes2label, genes_of_interest))
  complex_of_interest <- filter(umap_coords, gene %in% genes_of_interest)
  label_info <- data.frame(
    gene = "",
    umap_1 = mean(complex_of_interest$umap_1),
    umap_2 = mean(complex_of_interest$umap_2),
    color = "gray",
    size = 1,
    transparency = '1',
    genes2label = paste(complex_of_interest$gene, collapse = "_"),
    label =  complex_name1,
    label_x = mean(complex_of_interest$umap_1) + 0.25,
    label_y = mean(complex_of_interest$umap_2)
  )
  label_df <- as.data.frame(bind_rows(label_info, label_df))
  segments2draw <- data.frame(
    label = complex_name1,
    x_start = mean(complex_of_interest$umap_1),
    x_end = complex_of_interest$umap_1,
    y_start = mean(complex_of_interest$umap_2),
    y_end = complex_of_interest$umap_2
  )
  segment_df <- as.data.frame(bind_rows(segments2draw, segment_df))
  
  
  
  ## plot
  df4plot <- dplyr::select(umap_coords, c("gene", "umap_1", "umap_2")) %>%
    mutate(color = ifelse(gene %in% genes2label, "steelblue", "gray"), 
           size = 1, 
           transparency = ifelse(gene %in% unique(genes2label), 'highlighted', 'not highlighted'), 
           label = "", 
           label_x = NA, 
           label_y = NA)
  
  # color some of the closer complexes together
  #complexes2highlight <- c("RNA polymerase II core complex", "Cytoplasmic ribosome", "TOM40 complex")
  complexes2highlight <- c("")
  df4plot <- mutate(df4plot,
                    color = ifelse(gene %in% unlist(strsplit( label_df$genes2label[label_df$label %in% complexes2highlight], "_" )), "navy", color))
  
  df4plot <- bind_rows(df4plot, 
                       label_df)
  pdf(paste0(plotsdir, "/", date, "_umap_annotated.pdf"), width = 12, height = 12)
  p <- ggplot(df4plot[is.na(df4plot$genes2label),], 
              aes(x = umap_1, y = umap_2, alpha = transparency, col = color)) + 
    xlab("") + ylab("") + 
    geom_point() + 
    scale_color_manual(values = c('gray' = 'gray', 'red' = 'red', 'black' = 'black',
                                  'navy' = 'navy', 'steelblue' = 'steelblue',
                                  'darkgoldenrod' = 'darkgoldenrod')) + 
    scale_alpha_manual(values = c('highlighted' = 0.6, 'not highlighted' = 0.2, '1' = 1)) + 
    geom_segment(data = segment_df, aes(x = x_start, y = y_start, xend = x_end, yend = y_end), alpha = 0.4, col = 'gray') + 
    geom_text(data = df4plot[label != "",], 
              aes(x = label_x, y = label_y, label = label), nudge_x = -0.1, nudge_y = -0.1, col = 'black') + 
    theme_bw() + theme(legend.position = 'none')
  print(p)
  dev.off()
  
  complex_umap_anno <- list(
    complex_segment_df = segment_df,
    complex_label_df = label_df
  )
  
  saveRDS(complex_umap_anno , paste0(outdir, "/", date, "_complex_umap_anno.RDS"))
  
}








## ---- scratch

if (FALSE){
  genes_of_interest <- grep(umap_coords$gene, pattern = "SP", value = T)
  genes_of_interest <-c("HSP90B1", "PDIA3", "CALR", "HSPA5")
  genes_of_interest <- gprofiler_annotations[["GO:1904293"]]
  ##targets of PPRC1
  genes_of_interest <- magpie_lfcs_split[["PPRC1"]] %>%
    filter(pval_adj < sig_pval_thresh) %>%
    .$downstream_gene_name
  ggplot(umap_coords, aes(x = umap_1, y = umap_2,
                          label = ifelse(gene %in% genes_of_interest, gene, ""),
                          col = ifelse(gene %in% genes_of_interest, "y", 'no'))) + 
    scale_color_manual("", values = c("y" = 'red', 'no' = 'gray')) +
    geom_point(size = 0.4, alpha = 0.6) + ggrepel::geom_text_repel(max.overlaps = Inf) + theme_bw()
  
  
  
  df4heatmap <- magpie_lfcs_by_downstream[genes_of_interest] %>%
    bind_rows() %>% 
    as.data.frame() %>%
    dplyr::select(c("target", "downstream_gene_name", 'lfc')) %>%
    reshape2::dcast(target ~ downstream_gene_name, value.var = 'lfc') %>%
    column_to_rownames('target') %>%
    cor()
  my_heatmap(df4heatmap, min_c = -0.5, max_c = 0.5,
             treeheight_row = 0, treeheight_col = 0)
  
  dg_with_signal <- lapply(magpie_lfcs_by_downstream, FUN = function(df){
    dim(df %>% filter(pval_adj < sig_pval_thresh))[1]
  }) %>% unlist()
  dg_with_signal <- names(dg_with_signal)[which(dg_with_signal > 10)]
  df4heatmap <- magpie_lfcs_by_downstream[dg_with_signal] %>%
    bind_rows() %>% 
    as.data.frame() %>%
    dplyr::select(c("target", "downstream_gene_name", 'lfc')) %>%
    reshape2::dcast(target ~ downstream_gene_name, value.var = 'lfc') %>%
    column_to_rownames('target') %>%
    cor()
  df4plot <- data.frame(
    dg1 = magpie_lfcs_by_downstream[["IFNGR2"]]$lfc,
    dg2 = magpie_lfcs_by_downstream[["IFNGR1"]]$lfc
  )
  ggplot(df4plot, aes(x = dg1, y = dg2)) + 
    geom_point()
  pdf('scratch/2023-03-30_downstream_heatmap.pdf', width = 30, height = 30)
  print(my_heatmap(df4heatmap, min_c = -0.5, max_c = 0.5,
                   treeheight_row = 0, treeheight_col = 0))
  dev.off()
  
  claudia1 <- gprofiler2::gost(c("TNFRSF12A", "TNFRSF1A",  "STAT3", "STAT1", "SOCS1", "CD9"))$result
  claudia2 <- gprofiler2::gost(c("IFNGR2", "IFNGR1", "HMOSX1", "JUN"))$result
  gprofiler_annotations[["GO:0035976"]]
  
  
  ## some clusters i don't want
  
  complex_name1 <- "Calcium regulation"
  genes_of_interest <- c("HSP90B1", "PDIA3", "HSPA5") #skip out calr
  genes2label <- unique(c(genes2label, genes_of_interest))
  complex_of_interest <- filter(umap_coords, gene %in% genes_of_interest)
  label_info <- data.frame(
    gene = "",
    umap_1 = mean(complex_of_interest$umap_1),
    umap_2 = mean(complex_of_interest$umap_2),
    color = "gray",
    size = 1,
    transparency = '1',
    genes2label = paste(complex_of_interest$gene, collapse = "_"),
    label =  complex_name1,
    label_x = mean(complex_of_interest$umap_1),
    label_y = mean(complex_of_interest$umap_2)+ 0.15
  )
  label_df <- as.data.frame(bind_rows(label_info, label_df))
  
  
  
  complex_name1 <- "Protein export"
  genes_of_interest <- c("SPCS1", "SRP19", "SPCS2", "SRP68", "HSPA5", "SPCS3")
  genes2label <- unique(c(genes2label, genes_of_interest))
  complex_of_interest <- filter(umap_coords, gene %in% genes_of_interest)
  label_info <- data.frame(
    gene = "",
    umap_1 = mean(complex_of_interest$umap_1),
    umap_2 = mean(complex_of_interest$umap_2),
    color = "gray",
    size = 1,
    transparency = '1',
    genes2label = paste(complex_of_interest$gene, collapse = "_"),
    label =  complex_name1,
    label_x = mean(complex_of_interest$umap_1),
    label_y = mean(complex_of_interest$umap_2)+ 0.15
  )
  label_df <- as.data.frame(bind_rows(label_info, label_df))
  
  
  complex_name1 <- "Mediator complex"
  genes_of_interest <- grep(umap_coords$gene, pattern = "^MED", value = T)
  genes2label <- unique(c(genes2label, genes_of_interest))
  complex_of_interest <- filter(umap_coords, gene %in% genes_of_interest)
  label_info <- data.frame(
    gene = "",
    umap_1 = mean(complex_of_interest$umap_1),
    umap_2 = mean(complex_of_interest$umap_2),
    color = "gray",
    size = 1,
    transparency = '1',
    genes2label = paste(complex_of_interest$gene, collapse = "_"),
    label =  complex_name1,
    label_x = mean(complex_of_interest$umap_1),
    label_y = mean(complex_of_interest$umap_2)+ 0.15
  )
  label_df <- as.data.frame(bind_rows(label_info, label_df))
  
  
  
  complex_name1 <- "Eukaryotic translation elongation"
  genes_of_interest <- c("EEF1A1", "EEF1G", "EEF1B2", "EEF1D", "EEF2")
  genes2label <- unique(c(genes2label, genes_of_interest))
  complex_of_interest <- filter(umap_coords, gene %in% genes_of_interest)
  label_info <- data.frame(
    gene = "",
    umap_1 = mean(complex_of_interest$umap_1),
    umap_2 = mean(complex_of_interest$umap_2),
    color = "gray",
    size = 1,
    transparency = '1',
    genes2label = paste(complex_of_interest$gene, collapse = "_"),
    label =  complex_name1,
    label_x = mean(complex_of_interest$umap_1) - 0.1,
    label_y = mean(complex_of_interest$umap_2) + 0.1
  )
  label_df <- as.data.frame(bind_rows(label_info, label_df))
  
  
  complex_name1 <- "Translation initiation (EIF4/5)"
  genes_of_interest <- c("EIF4G3", "EIF4A3", "EIF3A", "EIF4A2", "EIF4G1", "EIF4EBP2", "EIF2AK4", "EIF5")
  genes2label <- unique(c(genes2label, genes_of_interest))
  complex_of_interest <- filter(umap_coords, gene %in% genes_of_interest)
  label_info <- data.frame(
    gene = "",
    umap_1 = mean(complex_of_interest$umap_1),
    umap_2 = mean(complex_of_interest$umap_2),
    color = "gray",
    size = 1,
    transparency = '1',
    genes2label = paste(complex_of_interest$gene, collapse = "_"),
    label =  complex_name1,
    label_x = mean(complex_of_interest$umap_1),
    label_y = mean(complex_of_interest$umap_2)
  )
  label_df <- as.data.frame(bind_rows(label_info, label_df))
  
  
  
  
  complex_name1 <- "Translation initiation (EIF1/3/6)"
  genes_of_interest <- c("EIF4B", "EIF2S3", "EIF3L", "EIF1", "EIF3K", "EIF3F", "EIF6", "EIF3M", "EIF3E", "EIF4A1")
  genes2label <- unique(c(genes2label, genes_of_interest))
  complex_of_interest <- filter(umap_coords, gene %in% genes_of_interest)
  label_info <- data.frame(
    gene = "",
    umap_1 = mean(complex_of_interest$umap_1),
    umap_2 = mean(complex_of_interest$umap_2),
    color = "gray",
    size = 1,
    transparency = '1',
    genes2label = paste(complex_of_interest$gene, collapse = "_"),
    label =  complex_name1,
    label_x = mean(complex_of_interest$umap_1),
    label_y = mean(complex_of_interest$umap_2) - 0.05
  )
  label_df <- as.data.frame(bind_rows(label_info, label_df))
  
  
  
  complex_name1 <- "Mitochondrial ribosome"
  genes_of_interest <- c("MRPS17", "MRPL38", "MRPL42", "MRPS16", "MRPL2", "MRPS33", "MRPS26", "MRPL40", "MRPL44", "MRPL17", "MRPS5", "MRPS25", "MRPL38", "MRPL2", "MRPL42", "MRPS17", "MRPL58", "MRPL16")
  genes2label <- unique(c(genes2label, genes_of_interest))
  complex_of_interest <- filter(umap_coords, gene %in% genes_of_interest)
  label_info <- data.frame(
    gene = "",
    umap_1 = mean(complex_of_interest$umap_1),
    umap_2 = mean(complex_of_interest$umap_2),
    color = "gray",
    size = 1,
    transparency = '1',
    genes2label = paste(complex_of_interest$gene, collapse = "_"),
    label =  complex_name1,
    label_x = mean(complex_of_interest$umap_1),
    label_y = mean(complex_of_interest$umap_2)+ 0.15
  )
  label_df <- as.data.frame(bind_rows(label_info, label_df))
  
  
  complex_name1 <- "PA700 complex"
  genes_of_interest <- c("PSMD7", "PSMD2", "PSMD11", "PSMC1", "PSMD1", "PSMD4", "PSMD13", "PSMC2", "PSMC4")
  genes2label <- unique(c(genes2label, genes_of_interest))
  complex_of_interest <- filter(umap_coords, gene %in% genes_of_interest)
  label_info <- data.frame(
    gene = "",
    umap_1 = mean(complex_of_interest$umap_1),
    umap_2 = mean(complex_of_interest$umap_2),
    color = "gray",
    size = 1,
    transparency = '1',
    genes2label = paste(complex_of_interest$gene, collapse = "_"),
    label =  complex_name1,
    label_x = mean(complex_of_interest$umap_1) + 0.25,
    label_y = mean(complex_of_interest$umap_2)
  )
  label_df <- as.data.frame(bind_rows(label_info, label_df))
  segments2draw <- data.frame(
    label = complex_name1,
    x_start = mean(complex_of_interest$umap_1),
    x_end = complex_of_interest$umap_1,
    y_start = mean(complex_of_interest$umap_2),
    y_end = complex_of_interest$umap_2
  )
  segment_df <- as.data.frame(bind_rows(segments2draw, segment_df))
  
  ##targets of PPRC1
  PPRC1_targets <- magpie_lfcs_split[["PPRC1"]] %>%
    filter(pval_adj < sig_pval_thresh) %>%
    .$downstream_gene_name
  
}



























































