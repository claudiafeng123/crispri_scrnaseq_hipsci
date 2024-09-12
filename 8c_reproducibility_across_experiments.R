
suppressMessages(library(data.table, warn.conflicts = F))
suppressMessages(library(ggplot2, warn.conflicts = F))
suppressMessages(library(dplyr, warn.conflicts = F))
suppressMessages(library(ggrastr, warn.conflicts = F))
suppressMessages(library(tidyverse, warn.conflicts = F))

section_name <- "8c_reproducibility_across_experiments"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
magpie_date <- "2022-08-15"
pica_date <- "2022-10-05"
replogle_date <- "2022-09-12"


args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  ##needed arguments
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  ## not needed
  if (arg == "--utils_path"){ utils_path <- args[[ind + 1]] }
  if (arg == "--io_path"){ io_path <- args[[ind + 1]] }
}



# set relevent i/o paths
source(paste0(HomeFolder, "scripts/io/Magpie_io.R"))
source(paste0(HomeFolder, "scripts/io/Magpie_Utils.R"))
outdir <- file.path(OutFolder, section_name)
plotsdir <- file.path(HTMLFolder, "pipeline/", section_name)
if(!dir.exists(outdir)) dir.create(outdir)
if(!dir.exists(plotsdir)) dir.create(plotsdir)
#if(!dir.exists(paste0(plotsdir, "/cor_plots/"))) dir.create(paste0(plotsdir, "/cor_plots/"))

## ---- LoadData

#GuideMetadata <- fread(GuideMetadataPath)
#target_genes <- unique(GuideMetadata$gene)

# by gene
magpie_res <- fread(file.path(OutFolder, "6b_calc_lfcs_transcriptome_wide_by_gene", paste0(magpie_date, "_combined_with_adj_pval.tsv.gz")))
pica_res <- fread(file.path(gsub(OutFolder, pattern = "Magpie", replacement = "Pica"), "6b_calc_lfcs_transcriptome_wide_by_gene", paste0(pica_date, "_combined_with_adj_pval.tsv.gz")))
replogle_rpe1 <- fread(file.path(ProjectFolder, "outs/Other/Replogle/", "3b_calc_lfcs_transcriptome_wide_by_gene", paste0(replogle_date, "_RPE1_raw_combined_with_adj_pval.tsv.gz") ))
replogle_k562_gwps <- fread(file.path(ProjectFolder, "outs/Other/Replogle/", "3b_calc_lfcs_transcriptome_wide_by_gene", paste0(replogle_date, "_K562_gwps_combined_with_adj_pval.tsv.gz") ))


## sort out p-values
adj_pval <- p.adjust(c(magpie_res$pval_lm, pica_res$pval_lm, replogle_rpe1$pval_lm, replogle_k562_gwps$pval_lm), method = "BH")
magpie_res$pval_adj <- adj_pval[1:dim(magpie_res)[1]]
pica_res$pval_adj <- adj_pval[(dim(magpie_res)[1] + 1):(dim(magpie_res)[1] + dim(pica_res)[1])]
replogle_rpe1$pval_adj <- adj_pval[(dim(magpie_res)[1] + dim(pica_res)[1] + 1):(dim(magpie_res)[1] + dim(pica_res)[1] + dim(replogle_rpe1)[1])]
replogle_k562_gwps$pval_adj <- adj_pval[(dim(magpie_res)[1] + dim(pica_res)[1] + dim(replogle_rpe1)[1] + 1):(dim(magpie_res)[1] + dim(pica_res)[1] + dim(replogle_rpe1)[1] + dim(replogle_k562_gwps)[1])]

# by complex

## ---- Compare On-Target LFC

n_deg <- filter(replogle_rpe1, pval_adj < sig_pval_thresh) %>% .$target %>% table()
df4plot1 <- filter(replogle_rpe1, downstream_gene  == target) %>%
  mutate(cell_type = "RPE1",
         n_deg = as.numeric(n_deg)[match(target, names(n_deg))]) %>%
  select(c("target", "downstream_gene", "cell_type", "lfc", "pval_adj", "n_deg")) %>%
  mutate(n_deg = ifelse(is.na(n_deg), 0, n_deg))
n_deg <- filter(replogle_k562_gwps, pval_adj < sig_pval_thresh) %>% .$target %>% table()
df4plot2 <- filter(replogle_k562_gwps, downstream_gene  == target) %>%
  mutate(cell_type = "K562",
         n_deg = as.numeric(n_deg)[match(target, names(n_deg))]) %>%
  select(c("target", "downstream_gene", "cell_type", "lfc", "pval_adj", "n_deg")) %>%
  mutate(n_deg = ifelse(is.na(n_deg), 0, n_deg))
n_deg <- filter(magpie_res, pval_adj < sig_pval_thresh) %>% .$target %>% table()
df4plot3 <- filter(magpie_res, downstream_gene_name  == target) %>%
  mutate(cell_type = "iPSC",
         n_deg = as.numeric(n_deg)[match(target, names(n_deg))]) %>%
  select(c("target", "downstream_gene", "cell_type", "lfc", "pval_adj", "n_deg")) %>%
  mutate(n_deg = ifelse(is.na(n_deg), 0, n_deg))
names(df4plot3) <- c("target", "downstream_gene", "cell_type", "lfc", "pval_adj", "n_deg")
target_meta_combined <- as.data.frame(bind_rows(df4plot1, df4plot2, df4plot3)) %>%
  select(c("target", "cell_type", "lfc", "pval_adj", "n_deg"))
#df4plot$pval_adj <- p.adjust(df4plot$pval_adj, method = "BH")


pdf(paste0(plotsdir, "/", magpie_date, "_target_gene_lfc.pdf"), width = 5, height = 10)
ggplot(target_meta_combined, aes(x = lfc, y = -log10(pval_adj))) +
  facet_wrap(~cell_type, ncol = 1) + 
  xlab("LFC") + ylab("-Log10 Adjusted P-value") +
  geom_point() + theme_bw()
dev.off()

## ---- Number Of DEG

df4plot <- select(target_meta_combined, c("target", "n_deg", "cell_type")) %>%
  reshape2::dcast(target ~ cell_type, value.var = "n_deg")
pdf(paste0(plotsdir, "/", magpie_date, "_n_deg_correlation.pdf"), width = 10, height = 10)
GGally::ggpairs(df4plot[,-1]) + theme_bw()
ggplot(df4plot, aes(x = K562, y = iPSC, label = ifelse(iPSC > 200 | K562 > 1500, target, ""))) +
  ggrepel::geom_text_repel() + 
  scale_x_log10() + scale_y_log10() + 
  geom_point() + theme_bw()
ggplot(df4plot, aes(x = RPE1, y = iPSC, label = ifelse(iPSC > 200 | RPE1 > 1500, target, ""))) +
  ggrepel::geom_text_repel() + 
  scale_x_log10() + scale_y_log10() + 
  geom_point() + theme_bw()
dev.off()



## ---- PlotConsistency
#consistency of target down-regulation across guides (all target x downstream gene pairs)

## magpie vs. pica
magpie_deg <- unique(magpie_res$downstream_gene_name)
pica_deg <- unique(pica_res$downstream_gene_name)
common_deg <- intersect(magpie_deg, pica_deg)
magpie_target <- unique(magpie_res$target)
pica_target <- unique(pica_res$target)
common_target <- intersect(magpie_target, pica_target)

df1 <- magpie_res %>% select(c("target", "downstream_gene_name", "lfc")) %>% filter(target %in% common_target & downstream_gene_name %in% common_deg)
names(df1) <- c("target", "downstream_gene_name", "magpie")
df2 <- pica_res %>% select(c("target", "downstream_gene_name", "lfc")) %>% filter(target %in% common_target & downstream_gene_name %in% common_deg)
names(df2) <- c("target", "downstream_gene_name", "pica")
df4plot <- left_join(df1, df2) %>%
  filter(target != downstream_gene_name  )
pdf(paste0(plotsdir, "/", magpie_date, "_magpie_v_pica.pdf") , width = 5, height = 5)
p <- ggplot(df4plot, aes(x = magpie, y = pica)) + 
  xlab("Genome-wide Screen") + ylab("Targeted Screen") + 
  geom_point(alpha = 0.2, size = 0.2) + theme_bw()
rasterize(p, layers='Point', dpi=300)
dev.off()


## magpie vs. rpe1
magpie_deg <- unique(magpie_res$downstream_gene_name)
replogle_rpe1_deg <- unique(replogle_rpe1$downstream_gene)
common_deg <- intersect(magpie_deg, replogle_rpe1_deg)
magpie_target <- unique(magpie_res$target)
replogle_rpe1_target <- unique(replogle_rpe1$target)
common_target <- intersect(magpie_target, replogle_rpe1_target)

df1 <- magpie_res %>% select(c("target", "downstream_gene_name", "lfc")) %>% filter(target %in% common_target & downstream_gene_name %in% common_deg)
names(df1) <- c("target", "downstream_gene_name", "magpie")
df2 <- replogle_rpe1 %>% select(c("target", "downstream_gene", "lfc")) %>% filter(target %in% common_target & downstream_gene %in% common_deg)
names(df2) <- c("target", "downstream_gene_name", "replogle_rpe1")
df4plot <- left_join(df1, df2) %>%
  filter(target != downstream_gene_name  )
pdf(paste0(plotsdir, "/", magpie_date, "_magpie_v_rpe1.pdf"), width = 5, height = 5)
p <- ggplot(df4plot, aes(x = magpie, y = replogle_rpe1)) + 
  xlab("iPSCs") + ylab("RPE1 Cells") + 
  geom_point(alpha = 0.2, size = 0.2) + theme_bw()
rasterize(p, layers='Point', dpi=300)
dev.off()


## magpie vs. k562 (do genome-wide only)
magpie_deg <- unique(magpie_res$downstream_gene_name)
replogle_k562_gwps_deg <- unique(replogle_k562_gwps$downstream_gene)
common_deg <- intersect(magpie_deg, replogle_k562_gwps_deg)
magpie_target <- unique(magpie_res$target)
replogle_k562_gwps_target <- unique(replogle_k562_gwps$target)
common_target <- intersect(magpie_target, replogle_k562_gwps_target)

df1 <- magpie_res %>% select(c("target", "downstream_gene_name", "lfc")) %>% filter(target %in% common_target & downstream_gene_name %in% common_deg)
names(df1) <- c("target", "downstream_gene_name", "magpie")
df2 <- replogle_k562_gwps %>% select(c("target", "downstream_gene", "lfc")) %>% filter(target %in% common_target & downstream_gene %in% common_deg)
names(df2) <- c("target", "downstream_gene_name", "replogle_k562_gwps")
df4plot <- left_join(df1, df2) %>%
  filter(target != downstream_gene_name  )
pdf(paste0(plotsdir, "/", magpie_date, "_magpie_v_k562_gwps.pdf"), width = 5, height = 5)
p <- ggplot(df4plot, aes(x = magpie, y = replogle_k562_gwps)) + 
  xlab("iPSCs") + ylab("K562 Cells") + 
  geom_point(alpha = 0.2, size = 0.2) + theme_bw()
rasterize(p, layers='Point', dpi=300)
dev.off()

### k562 vs. rpe1
replogle_rpe1_deg <- unique(replogle_rpe1$downstream_gene)
replogle_k562_gwps_deg <- unique(replogle_k562_gwps$downstream_gene)
common_deg <- intersect(replogle_k562_gwps_deg, replogle_rpe1_deg)
replogle_k562_gwps_target <- unique(replogle_k562_gwps$target)
replogle_rpe1_target <- unique(replogle_rpe1$target)
common_target <- intersect(replogle_rpe1_target, replogle_k562_gwps_target)

df1 <- replogle_rpe1 %>% select(c("target", "downstream_gene", "lfc")) %>% filter(target %in% common_target & downstream_gene %in% common_deg)
names(df1) <- c("target", "downstream_gene_name", "replogle_rpe1")
df2 <- replogle_k562_gwps %>% select(c("target", "downstream_gene", "lfc")) %>% filter(target %in% common_target & downstream_gene %in% common_deg)
names(df2) <- c("target", "downstream_gene_name", "replogle_k562_gwps")
df4plot <- left_join(df1, df2) %>%
  filter(target != downstream_gene_name  )
pdf(paste0(plotsdir, "/", magpie_date, "_rpe1_v_k562_gwps.pdf"), width = 5, height = 5)
p <- ggplot(df4plot, aes(x = replogle_rpe1, y = replogle_k562_gwps)) + 
  xlab("RPE1 Cells") + ylab("K562 Cells") + 
  geom_point(alpha = 0.2, size = 0.2) + theme_bw()
rasterize(p, layers='Point', dpi=300)
dev.off()

## ---- Correlation Between Targets
## some examples

magpie_deg <- unique(magpie_res$downstream_gene_name)
replogle_rpe1_deg <- unique(replogle_rpe1$downstream_gene)
common_deg <- intersect(magpie_deg, replogle_rpe1_deg)
magpie_target <- unique(magpie_res$target)
replogle_rpe1_target <- unique(replogle_rpe1$target)
common_target <- intersect(magpie_target, replogle_rpe1_target)

rpe1_df_split <- filter(replogle_rpe1, target %in% common_target & downstream_gene %in% common_deg) %>%
  select(c("target", "downstream_gene", "pval_adj", "lfc", "n_perturbed"))
names(rpe1_df_split) <- c("target", "downstream_gene_name", "pval_adj", "lfc", "n_perturbed")
rpe1_df_split <- split.data.frame(rpe1_df_split, f = rpe1_df_split$target)
magpie_df_split <- filter(magpie_res, target %in% common_target & downstream_gene_name %in% common_deg) %>%
  select(c("target", "downstream_gene_name", "pval_adj", "lfc", "n_perturbed"))
names(magpie_df_split) <- c("target", "downstream_gene_name", "pval_adj", "lfc", "n_perturbed")
magpie_df_split <- split.data.frame(magpie_df_split, f = magpie_df_split$target)

silencer <- mapply(common_target, FUN = function(tg){
  df1 <- rpe1_df_split[[tg]]
  df2 <- magpie_df_split[[tg]]
  rtn <- calc_cor(df1= df1, df2 = df2,
           xlab1 = "RPE1 Cells", ylab1 = "iPSCs", plot_title = paste0("Downstream Effects of ", tg, " Knockdowns"))
  if (rtn$df$n_deg_2 > 10 & rtn$df$n_deg_1 > 10 & rtn$df$common_deg > 10){
    rtn <- rtn
  } else {rtn <- rtn$df}
  return(rtn)
}, SIMPLIFY = F)

## all correlations
correlation_by_target <- mapply(silencer, FUN = function(rtn){
  if (class(rtn) != "data.frame"){
    rtn <- rtn$df
  }
  return(rtn)
}, SIMPLIFY = F) %>% bind_rows() %>% as.data.frame()
fwrite(correlation_by_target, paste0(outdir, "/", magpie_date, "_rpe1_vs_magpie_correlation_by_target.tsv"))
pdf(paste0(plotsdir, "/", magpie_date, "_rpe1_vs_magpie_correlation_by_target.pdf" ), width = 5, height = 4)
ggplot(correlation_by_target, aes( x = cor_all_deg )) + 
  xlab("Correlation of Downstream Effects (DEG only)") + ylab("") + ggtitle("RPE1 vs. iPS Cells") + 
  geom_density() + theme_bw()
df4plot <- correlation_by_target %>%
  mutate(large_effect = as.character(n_deg_2 > 10 & n_deg_1 > 10))
ggplot(df4plot, aes( x = cor_all_deg, fill = large_effect)) + 
  xlab("Correlation of Downstream Effects (DEG only)") + ylab("") + ggtitle("RPE1 vs. iPS Cells") + 
  geom_density(alpha = 0.4) + theme_bw()
dev.off()

## plot highly correlated pairs
pdf(paste0(plotsdir, "/", magpie_date, "_rpe1_vs_magpie_scatters_by_target.pdf" ), width = 5.5, height = 6)
s2 <- mapply(silencer, FUN = function(rtn){
  if ("p" %in% names(rtn)){
    print(rtn$p)
  }
})
dev.off()


magpie_deg <- unique(magpie_res$downstream_gene_name)
replogle_k562_gwps_deg <- unique(replogle_k562_gwps$downstream_gene)
common_deg <- intersect(magpie_deg, replogle_k562_gwps_deg)
magpie_target <- unique(magpie_res$target)
replogle_k562_gwps_target <- unique(replogle_k562_gwps$target)
common_target <- intersect(magpie_target, replogle_k562_gwps_target)

k562_gwps_df_split <- filter(replogle_k562_gwps, target %in% common_target & downstream_gene %in% common_deg) %>%
  select(c("target", "downstream_gene", "pval_adj", "lfc", "n_perturbed"))
names(k562_gwps_df_split) <- c("target", "downstream_gene_name", "pval_adj", "lfc", "n_perturbed")
k562_gwps_df_split <- split.data.frame(k562_gwps_df_split, f = k562_gwps_df_split$target)
magpie_df_split <- filter(magpie_res, target %in% common_target & downstream_gene_name %in% common_deg) %>%
  select(c("target", "downstream_gene_name", "pval_adj", "lfc", "n_perturbed"))
names(magpie_df_split) <- c("target", "downstream_gene_name", "pval_adj", "lfc", "n_perturbed")
magpie_df_split <- split.data.frame(magpie_df_split, f = magpie_df_split$target)

silencer <- mapply(common_target, FUN = function(tg){
  df1 <- k562_gwps_df_split[[tg]]
  df2 <- magpie_df_split[[tg]]
  rtn <- calc_cor(df1= df1, df2 = df2,
                  xlab1 = "K562 Cells", ylab1 = "iPSCs", plot_title = paste0("Downstream Effects of ", tg, " Knockdowns"))
  if (rtn$df$n_deg_2 > 10 & rtn$df$n_deg_1 > 10 & rtn$df$common_deg > 10){
    rtn <- rtn
  } else {rtn <- rtn$df}
  return(rtn)
}, SIMPLIFY = F)

## all correlations
correlation_by_target <- mapply(silencer, FUN = function(rtn){
  if (class(rtn) != "data.frame"){
    rtn <- rtn$df
  }
  return(rtn)
}, SIMPLIFY = F) %>% bind_rows() %>% as.data.frame()
fwrite(correlation_by_target, paste0(outdir, "/", magpie_date, "_k562_gwps_vs_magpie_correlation_by_target.tsv"))

pdf(paste0(plotsdir, "/", magpie_date, "_k562_gwps_vs_magpie_correlation_by_target.pdf" ), width = 5, height = 4)
ggplot(correlation_by_target, aes( x = cor_all_deg )) + 
  xlab("Correlation of Downstream Effects (DEG only)") + ylab("") + ggtitle("K562 vs. iPS Cells") + 
  geom_density() + theme_bw()
df4plot <- correlation_by_target %>%
  mutate(large_effect = as.character(n_deg_2 > 10 & n_deg_1 > 10))
ggplot(df4plot, aes( x = cor_all_deg, fill = large_effect)) + 
  xlab("Correlation of Downstream Effects (DEG only)") + ylab("") + ggtitle("K562 vs. iPS Cells") + 
  geom_density(alpha = 0.4) + theme_bw()
dev.off()

## plot highly correlated pairs
pdf(paste0(plotsdir, "/", magpie_date, "_k562_gwps_vs_magpie_scatters_by_target.pdf" ), width = 5.5, height = 6)
s2 <- mapply(silencer, FUN = function(rtn){
  if ("p" %in% names(rtn)){
    print(rtn$p)
  }
})
dev.off()


correlation_by_target_k562_gwps <- fread( paste0(outdir, "/", magpie_date, "_k562_gwps_vs_magpie_correlation_by_target.tsv"))
correlation_by_target_rpe1 <- fread( paste0(outdir, "/", magpie_date, "_rpe1_vs_magpie_correlation_by_target.tsv"))

df4plot1 <- correlation_by_target_k562_gwps %>%
  select(gene_1, common_deg, cor_all_deg) %>%
  mutate(`Cell Type` = "K562")
df4plot2 <- correlation_by_target_rpe1 %>%
  select(gene_1, common_deg, cor_all_deg) %>%
  mutate(`Cell Type` = "RPE1")
df4plot <- as.data.frame(bind_rows(df4plot1, df4plot2))
pdf(paste0(plotsdir, "/", magpie_date, "_correlation_by_target.pdf" ), width = 5.5, height = 4)
ggplot(df4plot, aes(x = cor_all_deg, fill = `Cell Type`)) + 
  xlab("Correlation") + ggtitle("Correlation of Downstream Profiles (DEG only)") + ylab("") + 
  geom_density(alpha = 0.4) + geom_vline(xintercept = 0, col = 'gray', lty = 3) + theme_bw()
dev.off()

## correlation between magpie and pica screens
common_targets <- intersect(unique(magpie_res$target), unique(pica_res$target))
common_downstream <- intersect(unique(magpie_res$downstream_gene_name), unique(pica_res$downstream_gene_name))
magpie_res4comparison <- magpie_res %>%
  filter(target %in% common_targets & downstream_gene_name %in% common_downstream) 
magpie_res4comparison <- split.data.frame(magpie_res4comparison, f = magpie_res4comparison$target)
pica_res4comparison <- pica_res %>%
  filter(target %in% common_targets & downstream_gene_name %in% common_downstream)
pica_res4comparison <- split.data.frame(pica_res4comparison, f = pica_res4comparison$target)

cor_magpie_pica <- mapply(common_targets, FUN = function(tg){
  rtn <- calc_cor(df1 = magpie_res4comparison[[tg]], 
           df2 = pica_res4comparison[[tg]],
           by_pval = F,
           xlab1 = "Genome-wide Screen", ylab1 = "Targeted Screen")
  return(rtn$df)
}, SIMPLIFY = F) %>% bind_rows() %>% as.data.frame()
fwrite(cor_magpie_pica, paste0(outdir, "/", magpie_date, "_", pica_date, "_magpie_pica_comparison_by_abs_lfc.tsv"), sep = "\t")
cor_magpie_pica <- fread(paste0(outdir, "/", magpie_date, "_", pica_date, "_magpie_pica_comparison_by_abs_lfc.tsv"))

cor_mat <- cor_magpie_pica %>%
  dplyr::select(
    "Correlation" = "cor_all",
    "Correlation (DEG only)" = "cor_all_deg",
    "DEG (Union Across Both Screens)" = "total_deg",
    "Cells in Genome-wide Screen" = 'n_cells_1',
    "Cells in Targeted Screen" = "n_cells_2",
    'DEG (Genome-wide Screen)' = 'n_deg_2',
    "DEG (Targeted Screen)" = 'n_deg_1',
    "DEG (Intersection Across Both Screens)" = 'common_deg'
  ) %>%
  cor()
pdf(paste0(plotsdir, "/", magpie_date, "_", pica_date, "_magpie_pica_correlation_heatmap.pdf"), height = 4, width = 7)
print(my_heatmap(cor_mat, min_c = -1, max_c = 1, treeheight_row = 0, treeheight_col = 0, show_colnames = F))
dev.off()

order_vars <- order(cor_mat["Correlation",colnames(cor_mat) != "Correlation"])
names_var <- names(cor_mat["Correlation",colnames(cor_mat) != "Correlation"])[order_vars]
pvals_cor_test <- sapply(names_var,
                         function(nm) cor.test(cor_mat[,nm], cor_mat[,"Correlation"], method = "spearman")$p.value)
padj <- p.adjust(pvals_cor_test, method = "BH")
#stopifnot(names(padj) == names_var)


df4plot <- data.frame(
  Correlation = cor_mat["Correlation", -1] %>% as.vector(),
  Name = names(cor_mat["Correlation", ])[-1],
  pval_adj = padj[match(names_var, names(padj))]
)

ggplot(df4plot, aes(x = reorder(Name, Correlation), y = Correlation,
                    fill = ifelse(pval_adj < 0.05, "sig", "not sig"))) + 
  scale_fill_manual(values = c("sig" = 'navy', 'not sig' = 'gray')) +
  xlab("") + ylab("Correlation") + ggtitle("Predictors of Consistency in Downstream Profiles") + 
  geom_text(aes(label=ifelse(padj<0.001,"***", ifelse(padj<0.01,"**", ifelse(padj<0.05,"*","")))), 
            vjust=0.6 - sign(df4plot$Correlation)/2) + 
  geom_bar(stat = 'identity') + theme_bw()+ theme(legend.position = 'none',
                                                  axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                                                  plot.margin = margin(0.5,0.5,1,2, "cm"))

order_vars <- order(cor_mat["Correlation (DEG only)",colnames(cor_mat) != "Correlation (DEG only)"])
names_var <- names(cor_mat["Correlation (DEG only)",colnames(cor_mat) != "Correlation (DEG only)"])[order_vars]
pvals_cor_test <- sapply(names_var,
                         function(nm) cor.test(cor_mat[,nm], cor_mat[,"Correlation (DEG only)"], method = "spearman")$p.value)
padj <- p.adjust(pvals_cor_test, method = "BH")
#stopifnot(names(padj) == names_var)


df4plot <- data.frame(
  Correlation = cor_mat["Correlation (DEG only)", -1] %>% as.vector(),
  Name = names(cor_mat["Correlation (DEG only)", ])[-1],
  pval_adj = padj[match(names_var, names(padj))]
)

ggplot(df4plot, aes(x = reorder(Name, Correlation), y = Correlation,
                    fill = ifelse(pval_adj < 0.05, "sig", "not sig"))) + 
  scale_fill_manual(values = c("sig" = 'navy', 'not sig' = 'gray')) +
  xlab("") + ylab("Correlation (DEG only)") + ggtitle("Predictors of Consistency in Downstream Profiles") + 
  geom_text(aes(label=ifelse(padj<0.001,"***", ifelse(padj<0.01,"**", ifelse(padj<0.05,"*","")))), 
            vjust=0.6 - sign(df4plot$Correlation)/2) + 
  geom_bar(stat = 'identity') + theme_bw()+ theme(legend.position = 'none',
                                                  axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                                                  plot.margin = margin(0.5,0.5,1,2, "cm"))


pdf(paste0(plotsdir, "/", magpie_date, "_", pica_date, "_magpie_pica_comparison_cor_density.pdf"), width = 6, height = 4.5)
ggplot(cor_magpie_pica, aes(x = cor_all)) + 
  xlab("Correlation") + ylab("") + ggtitle("Transcriptional Profiles of Targeted vs. Genome-wide Screen") + 
  geom_density() + theme_bw()
ggplot(cor_magpie_pica, aes(x = cor_all_deg)) + 
  xlab("Correlation (DEG only)") + ylab("") + ggtitle("Transcriptional Profiles of Targeted vs. Genome-wide Screen") + 
  geom_density() + theme_bw()
dev.off()

pdf(paste0(plotsdir, "/", magpie_date, "_", pica_date, "_magpie_pica_comparison.pdf"), width = 4.5, height = 4.5)
ggplot(cor_magpie_pica, aes(x = n_cells_1, y = cor_all)) + 
  xlab("Cells in Genome-wide Screen") + ylab("Correlation") + 
  geom_point() + theme_bw()
ggplot(cor_magpie_pica, aes(x = n_deg_2, y = cor_all)) + 
  xlab("DEG (Targeted Screen)") + ylab("Correlation") + 
  geom_point() + theme_bw()
ggplot(cor_magpie_pica, aes(x = n_deg_2, y = cor_all_deg)) + 
  xlab("DEG (Targeted Screen)") + ylab("Correlation") + 
  geom_point() + theme_bw()
ggplot(cor_magpie_pica, aes(x = n_cells_2, y = cor_all_deg)) + 
  xlab("Cells in Targeted Screen") + ylab("Correlation") + 
  geom_point() + theme_bw()
ggplot(cor_magpie_pica, aes(x = n_cells_1, y = cor_all_deg)) + 
  xlab("Cells in Genome-wide Screen") + ylab("Correlation") + 
  geom_point() + theme_bw()
ggplot(cor_magpie_pica, aes(x = n_deg_1, y = cor_all_deg)) + 
  scale_x_log10() + 
  xlab("DEG (Genome-wide Screen)") + ylab("Correlation") + 
  geom_point() + theme_bw()
dev.off()

## ---- More Magpie vs. Pica

df4comparison <- left_join(
  pica_res %>% dplyr::select(c('target', 'downstream_gene_name', 
                               'pica_lfc' = 'lfc', 'pica_pval_adj' = 'pval_adj',
                               'pica_n_perturbed' = 'n_perturbed')),
  magpie_res%>% dplyr::select(c('target', 'downstream_gene_name', 
                                  'magpie_lfc' = 'lfc', 'magpie_pval_adj' = 'pval_adj',
                                'magpie_n_perturbed' = 'n_perturbed'))
) %>% dplyr::filter(!(is.na(pica_lfc) | is.na(magpie_lfc)))

abs_err_cutoffs <- c(0.01, 0.05, 0.1, 0.2, 0.25, 0.5)
err_df <- lapply(abs_err_cutoffs, FUN = function(e){
  sum(abs(df4comparison$pica_lfc - df4comparison$magpie_lfc) < e)
})
df4plot <- data.frame(
  err = abs_err_cutoffs,
  n_genes = unlist(err_df),
  frac_genes = unlist(err_df)/dim(df4comparison)[1]
)

ggplot(df4plot, aes(x = as.character(err), y = frac_genes)) + 
  xlab("Abs. Difference in LFC") + ylab("Fraction of Genes") + 
  geom_bar(stat = 'identity', fill = target_downstream_col) + theme_bw()

## ---- By Complex



