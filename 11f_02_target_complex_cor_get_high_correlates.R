suppressWarnings(suppressMessages(library(data.table, warn.conflicts = F)))
suppressWarnings(suppressMessages(library(dplyr, warn.conflicts = F)))
suppressWarnings(suppressMessages(library(ggplot2, warn.conflicts = F)))
suppressWarnings(suppressMessages(library(dplyr, warn.conflicts = F)))
suppressWarnings(suppressMessages(library(ggpubr, warn.conflicts = F)))
suppressWarnings(suppressMessages(library(ggExtra, warn.conflicts = F)))
suppressWarnings(suppressMessages(library(doParallel, warn.conflicts = F)))

# set relevent i/o paths
section_name <- "11f_target_complex_cor"

# set relevent i/o paths
date <- "2022-08-15"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
ExperimentName <- "Magpie"
REDO <- T

args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  ##needed arguments
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--experiment_name"){ ExperimentName <- args[[ind + 1]] }
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
  
  ### not needed
  if (arg == "--utils_path"){ utils_path <- args[[ind + 1]] }
  if (arg == "--io_path"){ io_path <- args[[ind + 1]] }
  
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
#config <- paste0(date, "_min_degs-", min_deg_per_target)
outdir <- file.path(OutFolder, section_name)
if(!dir.exists(outdir)) dir.create(outdir, recursive = T)
plotsdir <- file.path(HTMLFolder, "pipeline", section_name)
if(!dir.exists(plotsdir)) dir.create(plotsdir, recursive = T)
#if(!dir.exists(file.path(plotsdir, "/by_pair/"))) dir.create(file.path(plotsdir, "/by_pair/"), recursive = T)

#if (min_deg_per_target <= 10){min_deg_per_target <- 10}


## ---- LoadData

complex_meta <- fread(paste0(OutFolder, "11d_get_clusters_to_aggregate/", date, "_complex_meta.tsv"))#fread(complex_meta_path)
genes_in_complexes <- unique(unlist(strsplit(complex_meta$genes_in_complex, "_")))

off_target <- fread(paste0(OutFolder, "/Preprocess_External_Datasets/off_target/off_target_pairs_magpie-", date, "-version.tsv.gz"), header = T) %>%
  dplyr::filter(target_gene != off_target_gene)
guide_consistency <- fread(paste0(OutFolder, "/8a_reproducibility_across_guides/", date, "_correlation_across_guides_by_gene.tsv.gz"))

high_cor_fnm <- paste0(outdir, "/", date, "_high_target_complex_cor.tsv")
if (file.exists(high_cor_fnm) & REDO == F){
  high_cor <- fread(high_cor_fnm)
} else {
  target_complex_cor <- mapply(1:dim(complex_meta)[1], FUN = function(i){
    fnm <- paste0(outdir, "/by_complex/", date, "_", i, "_target_complex_cor.tsv")
    if(file.exists(fnm)){
      rtn <- fread(paste0(outdir, "/by_complex/", date, "_", i, "_target_complex_cor.tsv"))
    } else {rtn <- NULL}
    return(rtn)
  }, SIMPLIFY = F) %>% bind_rows() %>% as.data.frame()
  fwrite(target_complex_cor, paste0(outdir, "/", date, "_target_complex_cor.tsv"), sep = "\t")
  high_cor <- target_complex_cor %>%
    dplyr::filter(cor_all > sig_cor_thresh_target & !(gene_1 %in% genes_in_complexes) & !(gene_1 %in% off_target$target_gene) & !(gene_1 %in% off_target$off_target_gene)) %>%
    mutate(guide_consistency = guide_consistency$max_cor_all[match(gene_1, guide_consistency$gene)]) %>%
    mutate(complex_name = complex_meta$complex_name[match(gene_2, complex_meta$genes_in_complex)]) %>%
    dplyr::select( "gene_1", "n_cells_1", "gene_2","complex_name", "cor_all", "cor_all_pval", "guide_consistency")
  fwrite(high_cor, paste0(outdir, "/", date, "_high_target_complex_cor.tsv"), sep = "\t")
}

lfc_by_complex <- lapply(1:dim(complex_meta)[1], FUN = function(complex_ind){
  fnm <- paste0(OutFolder, "13a_calc_lfcs_transcriptome_wide_by_complex/by_gene/", complex_ind, "_", date, ".tsv.gz")
  if(file.exists(fnm)){
    complex_lfc <- fread(fnm)
    complex_lfc <- mutate(complex_lfc, complex_name = complex_meta$complex_name[complex_ind],
                          downstream_gene_name = unlist(lapply(strsplit(downstream_gene, ":"), "[[", 2)))
  } else {complex_lfc <- NULL}
})


lfcs_by_target_split <- fread(paste0(OutFolder, "6b_calc_lfcs_transcriptome_wide_by_gene", "/", date, "_combined_with_adj_pval.tsv.gz"))
lfcs_by_target_split <- split.data.frame(lfcs_by_target_split, f = lfcs_by_target_split$target)


lfcs_by_guide_split <- fread(paste0(OutFolder, "6c_calc_lfcs_transcriptome_wide_by_guide", "/", date, "_combined_with_adj_pval.tsv.gz"))
lfcs_by_guide_split <- split.data.frame(lfcs_by_guide_split, f = lfcs_by_guide_split$guide)

## ---- Target_Complex_Cor


pdf(paste0(plotsdir, "/", date, "_high_cor_scatters.pdf"), width = 4.5, height = 4.5)
for (i in 1:dim(high_cor)[1]){
  #print(i)
  tg <- high_cor$gene_1[i]; dg <- high_cor$gene_2[i]
  complex_ind <- which(complex_meta$genes_in_complex == dg)
  
  
  df4plot <- data.frame(
    target_lfc = lfcs_by_target_split[[tg]]$lfc,
    complex_lfc = lfc_by_complex[[complex_ind]]$lfc
  )
  p1 <- ggplot(df4plot, aes(x = target_lfc, y = complex_lfc)) + 
    xlab(tg) + ylab(complex_meta$complex_name[complex_ind])  + 
    annotate("text", x = min(df4plot$target_lfc)*0.8, y = max(df4plot$complex_lfc)*0.8, 
             label = paste0("r = ", round(high_cor$cor_all[i], 4), ",\np = ", formatC(high_cor$cor_all_pval[i])), col = 'red') + 
    geom_point(size = 0.5) + theme_bw()
  
  #pdf(paste0(plotsdir, '/',date, "_target_complex_cor_", i, ".pdf"), width = 4.5, height = 4.5)
  print(p1)
  #dev.off()
  
}
dev.off()


pdf(paste0(plotsdir, "/", date, "_high_cor_scatters_with_guide_consistency.pdf"), width = 10, height = 4.5)
for (i in 1:dim(high_cor)[1]){
  print(i)
  tg <- high_cor$gene_1[i]; dg <- high_cor$gene_2[i]
  complex_ind <- which(complex_meta$genes_in_complex == dg)
  df4plot <- data.frame(
    target_lfc = lfcs_by_target_split[[tg]]$lfc,
    complex_lfc = lfc_by_complex[[complex_ind]]$lfc
  )
  p1 <- ggplot(df4plot, aes(x = target_lfc, y = complex_lfc)) + 
    xlab(tg) + ylab(complex_meta$complex_name[complex_ind])  + 
    annotate("text", x = min(df4plot$target_lfc)*0.8, y = max(df4plot$complex_lfc)*0.8, 
               label = paste0("r = ", round(high_cor$cor_all[i], 4), ",\np = ", formatC(high_cor$cor_all_pval[i])), col = 'red') + 
    geom_point(size = 0.5) + theme_bw()
  
  guide_inds <- which(unlist(lapply(strsplit(names(lfcs_by_guide_split), "-"), "[[", 1)) == tg)
  
  if (length(guide_inds) > 1){
    df4plot <- lapply(lfcs_by_guide_split[guide_inds], FUN = function(x){x$lfc}) %>% bind_cols() %>% as.data.frame()
    #print(dim(df4plot)[2])
    names(df4plot) <- paste0("guide_", 1:dim(df4plot)[2])
    highest_cor <- cor(df4plot)
    highest_cor <-  highest_cor[upper.tri( highest_cor, diag = F)] 
    highest_cor <- which.max(highest_cor)
    if (highest_cor == 1){
      g1 <- 1; g2 <- 2
    } else if (highest_cor == 2){
      g1 <- 1; g2 <- 3
    } else if (highest_cor == 3){
      g1 <- 2; g2 <- 3
    }
    df4plot <- dplyr::select(df4plot, c(paste0("guide_", g1), paste0("guide_", g2)))
    names(df4plot) <- c("guide_1", "guide_2")
    p2 <- ggplot(df4plot, aes(guide_1, y = guide_2)) + 
      xlab("Guide 1") + ylab("Guide 2") + 
      geom_point(size = 0.5) + theme_bw()
    #pdf(paste0(plotsdir,'/', date, "_target_complex_cor_", i, ".pdf"), width = 10, height = 4.5)
    print(ggarrange(p1, p2))
    #dev.off()
  } else {
   #pdf(paste0(plotsdir, '/',date, "_target_complex_cor_", i, ".pdf"), width = 4.5, height = 4.5)
    #print(p1)
    #dev.off()
  }
  
}
dev.off()

## ---- Case Studies

## spliceosome
#complex_of_interest <- c("Spliceosome, tri-SNP complex")
#genes_in_complex <- complex_meta %>% filter(complex_name == complex_of_interest) %>% .$genes_in_complex
#correlated_genes <- filter(high_cor, gene_2 == genes_in_complex) %>% .$gene_1

#genes_of_interest <- c("DDX41", "NKAP")

#for (i in 1:length(genes_of_interest)){
#  tg <- genes_of_interest[i]; dg <- genes_in_complex
#  complex_ind <- which(complex_meta$genes_in_complex == dg)
#  df4plot <- data.frame(
#    target_lfc = lfcs_by_target_split[[tg]]$lfc,
#    complex_lfc = lfc_by_complex[[complex_ind]]$lfc
#  )
#  p <- ggplot(df4plot, aes(x = target_lfc, y = complex_lfc)) + 
#    xlab(tg) + ylab(complex_meta$complex_name[complex_ind])  + 
#    annotate("text", x = min(df4plot$target_lfc)*0.8, y = max(df4plot$complex_lfc)*0.8, 
#             label = paste0("r = ", round(high_cor$cor_all[i], 4), ",\np = ", formatC(high_cor$cor_all_pval[i])), col = 'red') + 
#    geom_point(size = 0.5) + theme_bw()
#  print(p)
#}

## cholesterol biosynthesis pathway




