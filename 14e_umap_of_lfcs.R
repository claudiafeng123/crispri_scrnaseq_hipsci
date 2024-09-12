#xvfb-run -a /software/R-4.1.3/bin/R
date <- "2022-10-05"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
section_name <- "14e_umap_of_lfcs"
ExperimentName <- "Pica" # options are: "Magpie", "Pica" or "both"
## assigned - do for "Magpie", "Pica" and "both"

library(ggpubr)
library(data.table)
library(Seurat)
library(dplyr)
library(reshape2)
library(tidyverse)
library(umap)
library(plotly)
#library(viridis)

# options
args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--correction"){ correction <- args[[ind + 1]] }
  if (arg == "--experiment_name"){ ExperimentName <- args[[ind + 1]] }
  if (arg == "--perturbation_status"){ perturbation_status <- args[[ind + 1]] }
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
}

print(paste0("Experiment name: ", ExperimentName))

# i/o
print(paste0("date: ", date))
source(file.path(HomeFolder, "scripts/io/", paste0(ExperimentName, "_io.R")))
source(file.path(HomeFolder, "scripts/io/", paste0("Magpie", "_Utils.R")))
outdir <- file.path(OutFolder, section_name)
plotsdir <- file.path(HTMLFolder, "pipeline/", section_name)
if(!dir.exists(outdir)) dir.create(outdir)
if(!dir.exists(plotsdir)) dir.create(plotsdir)

## ---- LoadData
# 

lfcs <- fread("/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Pica/6f_calc_lfcs_transcriptome_wide_by_gene_per_line/2022-10-05_all_lines_with_adj_pval.tsv.gz")


## ---- LFCs (no adjustment)

deg <- lfcs %>%
  dplyr::filter(target != downstream_gene_name & abs(lfc) > 0.25 & abs(pval_adj) < sig_pval_thresh) %>%
  .$downstream_gene_name %>% table() %>% as.data.frame() %>%
  dplyr::filter(Freq > 5) %>% .$`.`
df4umap <- lfcs %>%
  dplyr::filter(downstream_gene_name %in% as.character(deg)) %>%
  mutate(to_match = paste0(cell_line, ':', target)) %>%
  dplyr::select(c("to_match", "downstream_gene_name", 'lfc')) %>%
  dcast(to_match ~ downstream_gene_name , value.var = 'lfc') %>%
  column_to_rownames('to_match')
df4umap['no_change:no_change',] <- 0

umap_res <- umap(df4umap)
umap_coords <- data.frame(
  to_match = row.names(umap_res$layout),
  target = unlist(lapply(strsplit(row.names(umap_res$layout), ":"), "[[", 2)),
  cell_line = unlist(lapply(strsplit(row.names(umap_res$layout), ":"), "[[", 1)),
  donor = unlist(lapply(strsplit(row.names(umap_res$layout), "_"), "[[", 1)),
  UMAP_1 = umap_res$layout[,1],
  UMAP_2 = umap_res$layout[,2]
)

umap_coords <- fread(paste0(outdir, '/', date, "_per_target_per_line_umap.tsv"))
umap_coords$to_match <- paste0(umap_coords$target, ":", umap_coords$cell_line)
p <- ggplot(umap_coords, aes(x = UMAP_1, y = UMAP_2, 
                             col = cell_line, label = to_match)) + 
  geom_point(size = 0.4) + 
  #geom_text(umap_coords %>% dplyr::filter(cell_line == 'no_change')) + 
  theme_bw() + theme(legend.position = 'none')
ggplotly(p)

fwrite(umap_coords, paste0(outdir, '/', date, "_per_target_per_line_umap.tsv"), sep='\t')

## ---- LFCs (no adjustment)

df4umap <- lfcs %>%
  dplyr::filter(downstream_gene_name %in% as.character(deg)) %>%
  mutate(to_match = paste0(cell_line, ':', target)) %>%
  dplyr::select(c("to_match", "downstream_gene_name", 'lfc')) %>%
  dcast(to_match ~ downstream_gene_name , value.var = 'lfc') %>%
  column_to_rownames('to_match')
df4umap['no_change:no_change',] <- 0

umap_res <- umap(df4umap)
umap_coords <- data.frame(
  to_match = row.names(umap_res$layout),
  target = unlist(lapply(strsplit(row.names(umap_res$layout), ":"), "[[", 2)),
  cell_line = unlist(lapply(strsplit(row.names(umap_res$layout), ":"), "[[", 1)),
  donor = unlist(lapply(strsplit(row.names(umap_res$layout), "_"), "[[", 1)),
  UMAP_1 = umap_res$layout[,1],
  UMAP_2 = umap_res$layout[,2]
)

p <- ggplot(umap_coords, aes(x = UMAP_1, y = UMAP_2, 
                             label = to_match,col = cell_line)) + 
  geom_point(size = 0.4) + 
  geom_text(umap_coords %>% dplyr::filter(cell_line == 'no_change'))
theme_bw() + theme(legend.position = 'none')
ggplotly(p)

fwrite(umap_coords, paste0(outdir, '/', date, "_per_target_per_line_umap.tsv"), sep='\t')

## ---- LFCs (subtract mean)

var_expl_df <- fread("/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Pica/15c_calc_var_expl_by_target/15c_05_calc_var_expl_by_target_check_pvals/2023-10-18_heritability_with_pval_adj.tsv.gz")

pairs2consider <- var_expl_df %>%
  dplyr::filter(lfc_dLL_donor_pval_adj < 0.1) 
pairs2consider <-   paste0(pairs2consider$target, ":", pairs2consider$downstream_gene_name)

norm_lfcs <- lfcs %>%
  mutate(to_match=paste0(target, ":", downstream_gene_name)) %>%
  dplyr::filter(to_match %in% pairs2consider) %>%
  group_by(to_match) %>%
  reframe(
    target = unique(target),
    downstream_gene_name = unique(downstream_gene_name),
    cell_line = cell_line,
    norm_lfc = (lfc - mean(lfc))/var(lfc)
  ) %>%
  dplyr::filter(!(is.na(norm_lfc)))
df4umap <- norm_lfcs %>%
  mutate(to_match = paste0(target, ":", downstream_gene_name)) %>%
  dplyr::filter(to_match %in% pairs2consider) %>%
  dplyr::select(c("to_match", 'cell_line', 'norm_lfc')) %>%
  reshape2::dcast(to_match ~ cell_line, value.var = 'norm_lfc') %>%
  column_to_rownames('to_match')
df4umap['no_change:no_change',] <- 0
df4umap[is.na(df4umap)] <- 0

umap_res <- umap(df4umap)
umap_coords <- data.frame(
  to_match = row.names(umap_res$layout),
  target = unlist(lapply(strsplit(row.names(umap_res$layout), ":"), "[[", 2)),
  downstream_gene_name = unlist(lapply(strsplit(row.names(umap_res$layout), ":"), "[[", 1)),
  donor = unlist(lapply(strsplit(row.names(umap_res$layout), "_"), "[[", 1)),
  UMAP_1 = umap_res$layout[,1],
  UMAP_2 = umap_res$layout[,2]
)

p <- ggplot(umap_coords, aes(x = UMAP_1, y = UMAP_2, 
                             label = to_match, col = target)) + 
  geom_point(size = 0.4) + 
  #geom_text(umap_coords %>% dplyr::filter(cell_line == 'no_change'))
  theme_bw() + theme(legend.position = 'none')
ggplotly(p)


## do per pair
lfcs_subset <- lfcs %>% 
  mutate(to_match = paste0(target, ":", downstream_gene_name)) %>%
  dplyr::filter(to_match %in% pairs2consider) %>%
  mutate(donor = unlist(lapply(strsplit(cell_line, "_"), "[[", 1)))
lfcs_subset <- split.data.frame(lfcs_subset, f = lfcs_subset$to_match)

pdf('scratch/2024-03-22_cell_line_embedding.pdf', width = 10, height = 3.5)
lapply(lfcs_subset, FUN = function(df){
  tm <- paste0(df$target[1], ":", df$downstream_gene_name[1])
  df4plot <- umap_coords %>%
    mutate(of_interest = ifelse(to_match == tm, 'y', 'n'))
  p_umap <- ggplot(df4plot, aes(x = UMAP_1, y = UMAP_2, 
                                    col = of_interest, size = of_interest)) + 
    scale_color_manual('', values = c('y' = 'red', 'n' = 'gray')) + 
    scale_size_manual('', values = c('y' = 2, 'n' = 0.4)) + 
    geom_point() + 
    theme_bw() + theme(legend.position = 'none')
  p_bar <- ggplot(df, aes(x = cell_line, y = lfc, fill = donor)) + 
    ggtitle(paste0(tm)) + 
    geom_bar(stat = 'identity') + theme_bw() + theme(legend.position = 'none', axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  print(ggarrange(p_umap, p_bar, widths = c(0.5, 1)))
  
})
dev.off()


fwrite(umap_coords, paste0(outdir, '/', date, "_cell_line_embedding.tsv"), sep='\t')

## ---- Save

## ----

sessionInfo()



pca_res <- prcomp(x = t(df4umap), center = F, scale. = F)
pca_coords <- data.frame(
  to_match = row.names(pca_res$rotation),
  target = unlist(lapply(strsplit(row.names(pca_res$rotation), ":"), "[[", 2)),
  cell_line = unlist(lapply(strsplit(row.names(pca_res$rotation), ":"), "[[", 1)),
  donor = unlist(lapply(strsplit(row.names(pca_res$rotation), ":"), "[[", 1)),
  PCA_1 = pca_res$rotation[,1],
  PCA_2 = pca_res$rotation[,2]
)
p <- ggplot(pca_coords, aes(x = PCA_1, y = PCA_2, 
                            label = to_match, col = cell_line)) + 
  geom_point(size = 0.4) + 
  #geom_text(umap_coords %>% dplyr::filter(cell_line == 'no_change'), aes(x = UMAP_1, y = UMAP_2)) + 
  theme_bw() + theme(legend.position = 'none')
ggplotly(p)
