#xvfb-run -a /software/R-4.1.3/bin/R

library(Seurat)
library(tidyverse)
library(data.table)
library(reshape2)
library(ggpointdensity)
#library(MOFAdata)
library(umap)
library(scales)
library(ggrepel)
library(ggpubr)
#library(factoextra)
library(ggExtra)
library(plotly)

HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
section_name <- "12a_downstream_gene_embedding"
date <- "2022-08-15"
ExperimentName <- "Magpie"
#min_cells <- 10
#min_target <- 25
#min_deg <- 1

# options
args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
  if (arg == "--experiment_name"){ ExperimentName <- args[[ind + 1]] }
  if (arg == "--min_cells"){ min_cells <- as.numeric(args[[ind + 1]])}
  if (arg == "--min_target"){ min_target <- as.numeric(args[[ind + 1]])}
  if (arg == "--min_deg"){ min_deg <- as.numeric(args[[ind + 1]])}
}

# i/o
print(paste0("date: ", date))
source(file.path(HomeFolder, "scripts/io/", paste0(ExperimentName, "_io.R")))
source(file.path(HomeFolder, "scripts/io/Magpie_Utils.R"))
outdir <- file.path(OutFolder, section_name)
plotsdir <- file.path(HTMLFolder, "pipeline/", section_name, '/')
if(!dir.exists(outdir)) dir.create(outdir)
if(!dir.exists(plotsdir)) dir.create(plotsdir)


## ---- Load Data

target_df <- fread(file.path(OutFolder, "7b_target_summary", paste0(date, "_target_meta_data.csv")))
downstream_meta <- fread(paste0(OutFolder, "10b_number_of_perturbing_knockdowns/", date, "_downstream_meta_data.csv"))

magpie_res <- fread(file.path(OutFolder, "6b_calc_lfcs_transcriptome_wide_by_gene",
                              paste0(date, "_combined_with_adj_pval.tsv.gz")))
magpie_res_split <- split.data.frame(magpie_res, f = magpie_res$target)


fnm <- sort(list.files(file.path(OutFolder, "Preprocess_External_Datasets", "gprofiler"), pattern= 'go_terms_', full.names = T), decreasing = T)[1]
gprofiler_terms <- fread(fnm)
gprofiler_annotations <- readRDS(file.path(ResourcesFolder, "gprofiler", "gprofiler_full_hsapiens.ENSG.RDS"))
omnipath_protein_complexes <- fread("/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Magpie/Preprocess_External_Datasets/omnipath/omnipath_complex_meta_from-2022-09-12_magpie-2022-08-15-version.txt")

tf_list <- fread('/lustre/scratch123/hgi/mdt2/teams/parts/jh47/claudia_files/DATA/DatabaseExtract_v_1.01.csv')
predicted_tfs <- tf_list %>% dplyr::filter(`Is TF?` == "Yes")
fnm <- sort(list.files(file.path(OutFolder, "Preprocess_External_Datasets", "dorothea"), pattern= 'dorothea_interaction_pairs', full.names = T), decreasing = T)[1]
dorothea_tfs <- fread(fnm)


fnm <- sort(list.files(file.path(OutFolder, "5b_expression_UMAPs"), pattern= 'Magpie_combined_all_cells_correction_technical_cluster_deg.tsv', full.names = T), decreasing = T)[1]
iPSC_subtype_deg <- fread(fnm)
low_novelty_genes <- iPSC_subtype_deg %>% dplyr::filter(avg_log2FC > 0) %>% .$downstream_gene_name
high_novelty_genes <- iPSC_subtype_deg %>% dplyr::filter(avg_log2FC < 0) %>% .$downstream_gene_name

fnm <- sort(list.files(file.path(OutFolder, 'Preprocess_External_Datasets', 'omnipath'), pattern = 'omnipath_complex_meta', full.names = T), decreasing = T)[1]
complex_meta <- fread(fnm)
complex_meta <- complex_meta %>% dplyr::filter(complex_name != "")
protein_complexes <- strsplit(complex_meta$genes_in_complex, "_")
names(protein_complexes) <- complex_meta$complex_name

fnm <- sort(list.files(paste0(ResourcesFolder, "/ipsc_marker_genes_from_sunay/"),
                       pattern = "Jerber_iNeuron_Differentiation_Efficiency_downloaded", full.names = T), decreasing = T)[1]
jerber_de_eff_deg <- fread(fnm)
msig_db <- readRDS("/lustre/scratch123/hgi/projects/crispri_scrnaseq/resources/gene_sets/h.all.v7.4.symbols.RDS")
canonical_gene_sets <- readRDS("/lustre/scratch123/hgi/projects/crispri_scrnaseq/resources/gene_sets/c2.cp.v2023.2.Hs.symbols_downloaded-on-2022-08-15.RDS")
perturbation_gene_sets <- readRDS("/lustre/scratch123/hgi/projects/crispri_scrnaseq/resources/gene_sets/h.all.v7.4.symbols.RDS")

## sc coexpression
sc_coexpression <- fread(file.path(OutFolder, "Preprocess_External_Datasets", "coexpression", paste0("sc_regressed_correlation_pairs_magpie-", date, "-version.tsv")))
#sc_coexpression$to_match <- paste0(sc_coexpression$gene_1, "_", sc_coexpression$gene_2)

off_target <- fread(paste0(OutFolder, "/Preprocess_External_Datasets/off_target/off_target_pairs_magpie-", date, "-version.tsv.gz"), header = T) %>%
  dplyr::filter(target_gene != off_target_gene)
genes2exclude <- sort(unique(off_target$target_gene, off_target$off_target_gene))

## ---- Compute Downstream-downstream correlation

anno_downstream_downstream_fnm <- paste0(OutFolder, section_name, "/", date, "_downstream_downstream_cor_with_adj_pval.tsv.gz")
if (!file.exists(anno_downstream_downstream_fnm)){
  
  cor_mat <- magpie_res %>%
    dplyr::select(c('target', 'downstream_gene_name', 'lfc')) %>%
    reshape2::dcast(target ~ downstream_gene_name, value.var = 'lfc') %>%
    column_to_rownames('target') %>%
    as.matrix() %>%
    cor() %>%
    as.data.frame() %>%
    rownames_to_column('downstream_gene_1') %>%
    reshape2::melt(id = 'downstream_gene_1')
  cor_mat$variable <- as.character(cor_mat$variable)
  names(cor_mat) <- c("downstream_gene_1", "downstream_gene_2", "iPSC_cor")
  
  ## add in target-target correlation from rpe1 and k562 cells
  for (cell_type in c("K562_essential", "K562_gwps", "RPE1_raw")){
    cell_type_downstream_cor <- fread(paste0(ProjectFolder, "outs/Other/Replogle/6a_analysis_heatmap/2022-09-12_", cell_type, "_downstream_downstream_cor.tsv.gz"))
    eval(parse(text = paste0("cor_mat <- left_join(cor_mat, ", 
                             " cell_type_downstream_cor %>% ", 
                             "dplyr::select(c('downstream_gene_1', 'downstream_gene_2', '", paste0(cell_type, '_cor'), "' = 'cor_all')))"
    )))
  }
  
  
  fwrite(cor_mat, anno_downstream_downstream_fnm, sep = '\t', compress = 'gzip')
}
downstream_downstream_cor_anno <- fread(anno_downstream_downstream_fnm)

## ---- Plots4Intuition

if (T){
  
  ## umap coords
  umap_coords <- fread(paste0(outdir, '/', date, "_downstream_umap_coords.tsv"))
  
  
  ## plot signal per target
  tgs <- target_df %>%
    dplyr::filter(n_cells >= min_cells_per_gene & n_downstream_excl_target > min_deg_per_target) %>% .$gene %>% sort()
  pdf(paste0(plotsdir, '/', date, "_downstream_umap_by_target.pdf"),
      width = 8, height = 6)
  for (tg in tgs){
    lfc_df <- magpie_res_split[[tg]]
    
    df4plot <- umap_coords
    df4plot$lfc <- lfc_df$lfc[match(df4plot$gene, lfc_df$downstream_gene_name)]
    df4plot$pval_adj <- lfc_df$pval_adj[match(df4plot$gene, lfc_df$downstream_gene_name)]
    df4plot$is_de <- ifelse(df4plot$pval_adj < sig_pval_thresh, 'y', 'n')
    if (length(which(df4plot$is_de == 'y')) > 3){
      p <- ggplot(df4plot %>% arrange(is_de), aes(x = umap_1, y = umap_2, 
                                                  col = ifelse(lfc < 0 & pval_adj < sig_pval_thresh, 'downregulated', 
                                                               ifelse(lfc > 0 & pval_adj < sig_pval_thresh, 'upregulated', 'neither')))) + 
        ggtitle(tg) + 
        scale_color_manual('', values = c('downregulated' = down_regulated_blue, 'upregulated' = up_regulated_red, 'neither' = 'gray')) + 
        xlab("UMAP_1") + ylab("UMAP_2") + 
        geom_point(size = 1) + theme_bw() + theme(legend.position = 'none')
      print(p)}
    
  }
  dev.off()
  
  ## plot signal per msigdb pathway
  msigdb_plots <- mapply(pathway_ind = 1:length(msig_db), FUN = function(pathway_ind){
    genes_of_interset <- intersect(msig_db[[pathway_ind]], umap_coords$gene)
    df4plot <- umap_coords
    df4plot$in_pathway <- ifelse(df4plot$gene %in% genes_of_interset, 'y', 'n')
    if (length(which(df4plot$in_pathway == 'y')) > 3){
      
      df4heatmap_iPSCs <- downstream_downstream_cor_anno %>%
        dplyr::select(c('downstream_gene_1', 'downstream_gene_2', 'iPSC_cor')) %>%
        dplyr::filter(downstream_gene_1 %in% genes_of_interset & downstream_gene_2 %in% genes_of_interset) %>%
        reshape2::dcast(downstream_gene_1 ~ downstream_gene_2, value.var = 'iPSC_cor') %>%
        column_to_rownames('downstream_gene_1') %>%
        as.matrix()
      df4heatmap_K562 <- downstream_downstream_cor_anno %>%
        dplyr::select(c('downstream_gene_1', 'downstream_gene_2', 'K562_gwps_cor')) %>%
        dplyr::filter(downstream_gene_1 %in% genes_of_interset & downstream_gene_2 %in% genes_of_interset) %>%
        reshape2::dcast(downstream_gene_1 ~ downstream_gene_2, value.var = 'K562_gwps_cor') %>%
        column_to_rownames('downstream_gene_1') %>%
        as.matrix()
      
      p_heatmap <- double_heatmap(cor_mat1 = df4heatmap_iPSCs, cor_mat2 = df4heatmap_K562,
                          min_c = -0.6, max_c = 0.6,
                          treeheight_row = 0, treeheight_col = 0, border_col = NA)
      
      p_umap <- ggplot(df4plot %>% arrange(in_pathway), aes(x = umap_1, y = umap_2, 
                                                  col = in_pathway,
                                                  label = ifelse(in_pathway == 'y', gene, ''))) + 
        ggtitle(names(msig_db)[pathway_ind]) + 
        ggrepel::geom_text_repel(max.overlaps = Inf) + 
        scale_color_manual('', values = c('y' = up_regulated_red, 'n' = 'gray')) + 
        xlab("UMAP_1") + ylab("UMAP_2") + 
        geom_point(size = 1) + theme_bw() + theme(legend.position = 'none')
      rtn <- ggarrange(p_heatmap[[4]], p_umap)
    }
    
  }, SIMPLIFY = F)
  pdf(paste0(plotsdir, '/', date, "_downstream_umap_by_msigdb_hallmark_pathway.pdf"),
      width = 20, height = 8)
  msigdb_plots
  dev.off()
  
  pdf(paste0(plotsdir, '/', date, "_downstream_umap_by_msigdb_canonical_pathway.pdf"),
      width = 8, height = 6)
  for (pathway_ind in 1:length(canonical_gene_sets)){
    genes_of_interset <- canonical_gene_sets[[pathway_ind]]
    df4plot <- umap_coords
    df4plot$in_pathway <- ifelse(df4plot$gene %in% genes_of_interset, 'y', 'n')
    if (length(which(df4plot$in_pathway == 'y')) > 3){
      p <- ggplot(df4plot %>% arrange(in_pathway), aes(x = umap_1, y = umap_2, 
                                                       col = in_pathway,
                                                       label = ifelse(in_pathway == 'y', gene, ''))) + 
        ggtitle(names(canonical_gene_sets)[pathway_ind]) + 
        ggrepel::geom_text_repel(max.overlaps = Inf) + 
        scale_color_manual('', values = c('y' = up_regulated_red, 'n' = 'gray')) + 
        xlab("UMAP_1") + ylab("UMAP_2") + 
        geom_point(size = 1) + theme_bw() + theme(legend.position = 'none')
      print(p)
    }
    
  }
  dev.off()
  
  pdf(paste0(plotsdir, '/', date, "_downstream_umap_by_msigdb_perturbation_gene_set.pdf"),
      width = 8, height = 6)
  for (pathway_ind in 1:length(perturbation_gene_sets)){
    genes_of_interset <- perturbation_gene_sets[[pathway_ind]]
    df4plot <- umap_coords
    df4plot$in_pathway <- ifelse(df4plot$gene %in% genes_of_interset, 'y', 'n')
    if (length(which(df4plot$in_pathway == 'y')) > 3){
      p <- ggplot(df4plot %>% arrange(in_pathway), aes(x = umap_1, y = umap_2, 
                                                       col = in_pathway,
                                                       label = ifelse(in_pathway == 'y', gene, ''))) + 
        ggtitle(names(perturbation_gene_sets)[pathway_ind]) + 
        ggrepel::geom_text_repel(max.overlaps = Inf) + 
        scale_color_manual('', values = c('y' = up_regulated_red, 'n' = 'gray')) + 
        xlab("UMAP_1") + ylab("UMAP_2") + 
        geom_point(size = 1) + theme_bw() + theme(legend.position = 'none')
      print(p)
    }
    
  }
  dev.off()
  
  pdf(paste0(plotsdir, '/', date, "_downstream_umap_by_msigdb_cell_type_gene_set.pdf"),
      width = 8, height = 6)
  for (pathway_ind in 1:length(cell_type_gene_sets)){
    genes_of_interset <- cell_type_gene_sets[[pathway_ind]]
    df4plot <- umap_coords
    df4plot$in_pathway <- ifelse(df4plot$gene %in% genes_of_interset, 'y', 'n')
    if (length(which(df4plot$in_pathway == 'y')) > 3){
      p <- ggplot(df4plot %>% arrange(in_pathway), aes(x = umap_1, y = umap_2, 
                                                       col = in_pathway,
                                                       label = ifelse(in_pathway == 'y', gene, ''))) + 
        ggtitle(names(cell_type_gene_sets)[pathway_ind]) + 
        ggrepel::geom_text_repel(max.overlaps = Inf) + 
        scale_color_manual('', values = c('y' = up_regulated_red, 'n' = 'gray')) + 
        xlab("UMAP_1") + ylab("UMAP_2") + 
        geom_point(size = 1) + theme_bw() + theme(legend.position = 'none')
      print(p)
    }
    
  }
  dev.off()
  
  ## subset to transcription factors
  ## plot signal per target
  tgs <- target_df %>%
    dplyr::filter(n_cells >= min_cells_per_gene & n_downstream_excl_target > min_deg_per_target) %>% .$gene
  tgs <- sort(tgs[which(tgs %in% predicted_tfs$`HGNC symbol`)])
  pdf(paste0(plotsdir, '/', date, "_downstream_umap_tfs_only.pdf"),
      width = 8, height = 6)
  for (tg in tgs){
    lfc_df <- magpie_res_split[[tg]]
    
    df4plot <- umap_coords
    df4plot$lfc <- lfc_df$lfc[match(df4plot$gene, lfc_df$downstream_gene_name)]
    df4plot$pval_adj <- lfc_df$pval_adj[match(df4plot$gene, lfc_df$downstream_gene_name)]
    df4plot$is_de <- ifelse(df4plot$pval_adj < sig_pval_thresh, 'y', 'n')
    if (length(which(df4plot$is_de == 'y')) > 2){
      p <- ggplot(df4plot %>% arrange(is_de), aes(x = umap_1, y = umap_2, 
                                                  col = ifelse(lfc < 0 & pval_adj < sig_pval_thresh, 'downregulated', 
                                                               ifelse(lfc > 0 & pval_adj < sig_pval_thresh, 'upregulated', 'neither')))) + 
        ggtitle(tg) + 
        scale_color_manual('', values = c('downregulated' = down_regulated_blue, 'upregulated' = up_regulated_red, 'neither' = 'gray')) + 
        xlab("UMAP_1") + ylab("UMAP_2") + 
        geom_point(size = 1) + theme_bw() + theme(legend.position = 'none')
      print(p)}
    
  }
  dev.off()
  
  ## transcription factor targets
  pdf(paste0(plotsdir, '/', date, "_downstream_umap_by_tf.pdf"),
      width = 8, height = 6)
  tfs <- sort(unique(dorothea_tfs$source))
  for (tf in tfs){
    lfc_df <- magpie_res_split[[tf]]
    
    df4plot <- umap_coords
    df4plot$is_target <- ifelse(df4plot$gene %in% 
                                  (dorothea_tfs %>%
                                     dplyr::filter(source == tf) %>%
                                     .$target), "y", "n")
    if (length(which(df4plot$is_target == 'y') > 2)){
      p <- ggplot(df4plot %>%
                    arrange(is_target), aes(x = umap_1, y = umap_2, 
                                            col = is_target,
                                            label = ifelse(is_target == 'y', gene, ''))) + 
        ggtitle(tf) + 
        scale_color_manual('', values = c('y' = up_regulated_red, 'n' = 'gray')) + 
        ggrepel::geom_text_repel(max.overlaps = Inf) + 
        xlab("UMAP_1") + ylab("UMAP_2") + 
        geom_point(size = 1) + theme_bw() + theme(legend.position = 'none')
      print(p)
    }
    
    
  }
  dev.off()
  
  ## ipscs
  if (T){
    pdf(paste0(plotsdir, '/', date, "_downstream_umap_ipsc_status.pdf"),
        width = 8, height = 6)
    
    df4plot <- umap_coords
    df4plot <- left_join(umap_coords,
                         jerber_de_eff_deg %>%
                           dplyr::select(c('hgnc_symbol', 'coef')),
                         by = c('gene' = 'hgnc_symbol'))
    df4plot <- df4plot %>%
      mutate(coef = ifelse(is.na(coef), 'aaa_neither', ifelse(coef < 0, 'team_utf1', 'pro-differentiation' )))
    p <- ggplot(df4plot %>% arrange(coef), aes(x = umap_1, y = umap_2, 
                                                    col = coef,
                                                    label = ifelse(coef!='aaa_neither', gene, ''))) + 
      ggtitle("DEG in Jerber") + 
      ggrepel::geom_text_repel(max.overlaps = Inf) +
      scale_color_manual('', values = c('team_utf1' = down_regulated_blue, 'pro-differentiation' = up_regulated_red, 'aaa_neither' = 'gray')) + 
      xlab("UMAP_1") + ylab("UMAP_2") + 
      geom_point(size = 1) + theme_bw() + theme(legend.position = 'none')
    print(p)
    
    ## pou5f1 knockdowns
    df4plot <- umap_coords
    df4plot$ipsc_de = ifelse(iPSC_subtype_deg$avg_log2FC[match(umap_coords$gene, iPSC_subtype_deg$downstream_gene_name)] > 0, "upregulated", "downregulated")
    df4plot$ipsc_de[is.na(df4plot$ipsc_de)] <- 'aaa_neither'
    p <- ggplot(df4plot %>% arrange(ipsc_de), aes(x = umap_1, y = umap_2, 
                                                  col = ipsc_de,
                                                  label = ifelse(ipsc_de!='aaa_neither', gene, ''))) + 
      ggtitle("DEG in iPSC subtypes") + 
      ggrepel::geom_text_repel(max.overlaps = 30) +
      scale_color_manual('', values = c('downregulated' = down_regulated_blue, 'upregulated' = up_regulated_red, 'aaa_neither' = 'gray')) + 
      xlab("UMAP_1") + ylab("UMAP_2") + 
      geom_point(size = 1) + theme_bw() + theme(legend.position = 'none')
    print(p)
    
    ## pou5f1 lfc
    df4plot <- umap_coords %>%
      mutate(pou5f1_lfc = magpie_res_split[["POU5F1"]]$lfc[match(gene, magpie_res_split[["POU5F1"]]$downstream_gene_name)]) %>%
      mutate(pou5f1_lfc = ifelse(pou5f1_lfc  < -0.5, -0.5, ifelse(pou5f1_lfc > 0.5, 0.5, pou5f1_lfc )))
    p <- ggplot(df4plot %>% arrange(pou5f1_lfc), aes(x = umap_1, y = umap_2, 
                                                  col = pou5f1_lfc
                                                  )) + 
      ggtitle("DEG in iPSC subtypes") + 
      scale_color_gradient2(
        low = down_regulated_blue, mid = "white", high = up_regulated_red, 
        midpoint = 0
      ) + 
      #scale_alpha_continuous() + 
      xlab("UMAP_1") + ylab("UMAP_2") + 
      geom_point() + theme_bw() + theme(legend.position = 'none')
    print(p)
    
    dev.off()
  }
  
  
  ## expression level
  if (T){
    pdf(paste0(plotsdir, '/', date, "_downstream_umap_wt_expression_level.pdf"),
        width = 8, height = 6)
    df4plot <- umap_coords
    df4plot$expr_de = magpie_res_split[[1]]$control_norm_expr[match(df4plot$gene, magpie_res_split[[1]]$downstream_gene_name)]
    p <- ggplot(df4plot, aes(x = umap_1, y = umap_2, 
                             col = expr_de)) + 
      ggtitle("Wild-type Expression") + 
      xlab("UMAP_1") + ylab("UMAP_2") + 
      geom_point(size = 1) + theme_bw() + theme(legend.position = 'none')
    print(p)
    df4plot <- umap_coords
    df4plot$ipsc_de = ifelse(iPSC_subtype_deg$avg_log2FC[match(umap_coords$gene, iPSC_subtype_deg$downstream_gene_name)] > 0, "upregulated", "downregulated")
    df4plot$ipsc_de[is.na(df4plot$ipsc_de)] <- 'aaa_neither'
    p <- ggplot(df4plot %>% arrange(ipsc_de), aes(x = umap_1, y = umap_2, 
                                                  col = ipsc_de,
                                                  label = ifelse(ipsc_de!='aaa_neither', gene, ''))) + 
      ggtitle("DEG in iPSC subtypes") + 
      ggrepel::geom_text_repel(max.overlaps = 30) +
      scale_color_manual('', values = c('downregulated' = down_regulated_blue, 'upregulated' = up_regulated_red, 'aaa_neither' = 'gray')) + 
      xlab("UMAP_1") + ylab("UMAP_2") + 
      geom_point(size = 1) + theme_bw() + theme(legend.position = 'none')
    print(p)
    dev.off()
  }
  
  
}

## ---- Run Rmarkdown


# run Rmd and create htmls
rmarkdown::render(input = paste0(CodeFolder, ExperimentName, "/pipeline/", section_name, ".Rmd" ),
                  output_file = file.path(plotsdir, paste0(date, "_", section_name, ".html")),
                  params =  list(experiment_name = ExperimentName,
                                 date = date,
                                 section_name = section_name))
 
## ---- What's driving the singal?



## ---- scratch

p <- ggplot(umap_coords, aes(x = umap_1, y = umap_2, 
                            label = gene)) + 
  ggtitle("Downstream Genes") + 
  ggrepel::geom_text_repel(max.overlaps = 50) + 
  xlab("UMAP_1") + ylab("UMAP_2") + 
  geom_point(size = 0.4) + theme_bw() + theme(legend.position = 'none')

  
plotly::ggplotly(p)

genes2annotate <- gprofiler_annotations[['GO:0007259']]

umap_coords$gene[!(umap_coords$gene %in% unlist(clusters2annotate))]
ggplot(umap_coords, aes(x = umap_1, y = umap_2, 
                                label = ifelse(gene %in% genes2annotate, gene, ''))) + 
  ggtitle("What's left") + 
  #scale_color_manual('', values = c('downregulated' = down_regulated_blue, 'upregulated' = up_regulated_red, 'neither' = 'gray')) + 
  ggrepel::geom_text_repel(max.overlaps = 50) + 
  xlab("UMAP_1") + ylab("UMAP_2") + 
  geom_point(size = 0.4) + theme_bw() + theme(legend.position = 'none')

claudia <- gprofiler2::gost(genes_left)$result


df4plot <- data.frame(
  thy1 = magpie_res %>% dplyr::filter(downstream_gene_name == "TAGLN") %>% .$lfc,
  id3 = magpie_res %>% dplyr::filter(downstream_gene_name == "MYL9") %>% .$lfc,
  target = magpie_res %>% dplyr::filter(downstream_gene_name == "THY1") %>% .$target
)
p <- ggplot(df4plot, aes(x =thy1, y = id3,
                    label = target)) + 
  geom_point(size = 0.4) + theme_bw()
plotly::ggplotly(p)

genes_of_interset <- msig_db[[pathway_ind]]
df4plot <- umap_coords
df4plot$in_pathway <- ifelse(df4plot$gene %in% genes_of_interset, 'y', 'n')
ggplot(df4plot %>% arrange(in_pathway), aes(x = umap_1, y = umap_2, 
                                            col = in_pathway,
                                            label = ifelse(in_pathway == 'y', gene, ''))) + 
  ggtitle(names(msig_db)[pathway_ind]) + 
  ggrepel::geom_text_repel(max.overlaps = Inf) + 
  scale_color_manual('', values = c('y' = up_regulated_red, 'n' = 'gray')) + 
  xlab("UMAP_1") + ylab("UMAP_2") + 
  geom_point(size = 1) + theme_bw() + theme(legend.position = 'none')



