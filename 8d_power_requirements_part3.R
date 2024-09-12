
suppressMessages(library(data.table, warn.conflicts = F))
suppressMessages(library(ggplot2, warn.conflicts = F))
suppressMessages(library(dplyr, warn.conflicts = F))
suppressMessages(library(ggpubr, warn.conflicts = F))


section_name <- "8d_power_requirements"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
ExperimentName <- "Magpie"

#for testing
#target_gene_name <- "PAF1"
#n_cells <- 8
date <- "2022-08-15"
#num_cells <- 2:50
#n_iterations <- 100
REDO <- F

args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  ##needed arguments
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--experiment_name"){ ExperimentName <- args[[ind + 1]] }
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

GuideMetadata <- fread(GuideMetadataPath)
target_genes <- unique(GuideMetadata$gene)

target_meta <- fread(paste0(OutFolder, "/7b_target_summary/", date, "_target_meta_data.csv"))
target_gene_list <- list.files(paste0(OutFolder, "/8d_power_requirements/by_gene/"))

## ---- Subsampling stats

if (REDO == F){
  fnms <- grep(list.files(paste0(OutFolder, "/8d_power_requirements/"), pattern = "tsv"),
               pattern = "targets2consider", invert = T, value = T)
  dfs4plot <-  mapply(fnms, FUN = function(fnm){
    rtn <- fread(paste0(OutFolder, "/8d_power_requirements/", fnm))
  }, SIMPLIFY = F)
  target_gene_list <- gsub(fnms, pattern = paste0(date, "_|[.]tsv"), replacement = "")
  names(dfs4plot) <- target_gene_list
} else {
  
  lfcs_by_target <- fread(paste0(OutFolder, "/", "6b_calc_lfcs_transcriptome_wide_by_gene", "/", date, "_combined_with_adj_pval.tsv.gz"))
  lfcs_by_target <- split.data.frame(lfcs_by_target, f =lfcs_by_target$target)
  all_downstream_genes <- lfcs_by_target[[1]]$downstream_gene
  n_downstream_genes <- length(all_downstream_genes)
  
  
  dfs4plot <- mapply(target_gene_list, FUN = function(target_gene_name){
    
    fnm_out <- paste0(outdir, "/", date, "_", target_gene_name, ".tsv")
    if (REDO == T | !(file.exists(fnm_out))){
      
      lfcs <- lfcs_by_target[[target_gene_name]]
      
      lfcs$lfc_upper <- 1.1*lfcs$lfc; lfcs$lfc_lower <- 0.9*lfcs$lfc
      lfcs$pval_adj <- p.adjust(lfcs$pval_lm, method = "BH")
      degs_0.05 <- filter(lfcs, pval_adj < 0.05 ) %>% .$downstream_gene
      degs_0.1 <- filter(lfcs, pval_adj < 0.1 ) %>% .$downstream_gene
      degs <- filter(lfcs, pval_adj < sig_pval_thresh)
      sig_lfcs <- lfcs$lfc[match(degs$downstream_gene, lfcs$downstream_gene)]
      num_cells <- grep(gsub(list.files(paste0(OutFolder, "/8d_power_requirements/by_gene/", target_gene_name, "/"), pattern = date),
                             pattern = paste0(target_gene_name, "_", date, "_|-cells_.*"), replacement = ""),
                        pattern = "residuals", value = T, invert = T) %>% 
        unique() %>%
        as.numeric() %>% sort()
      
      
      
      df4plot <- mapply(num_cells, FUN = function(n_cells){
        
        fnms <- list.files(paste0(OutFolder, "/8d_power_requirements/by_gene/", target_gene_name, "/"), pattern = paste0(date, "_", n_cells, "-cells"))
        cor_per_cell_no <- mapply(fnms, FUN = function(fnm){
          
          df <- fread(paste0(OutFolder, "/8d_power_requirements/by_gene/", target_gene_name, "/", fnm))
          df$pval_adj <- p.adjust(df$pval_lm, method = "BH")
          
          cor_deg <- cor(sig_lfcs, lfcs4cor_sig)
          cor_all <- cor(df$lfc, lfcs$lfc)
          
          subsampled_degs_0.05 <- filter(df, pval_adj < 0.05) %>% .$downstream_gene
          subsampled_degs_0.1 <- filter(df, pval_adj < 0.1) %>% .$downstream_gene
          
          ## fraction of DEGs within 10% of true value
          n_genes_within_10pct <- length(which(df$lfc < lfcs$lfc_upper & df$lfc > lfcs$lfc_lower))
          frac_genes_within_10pct <- n_de_genes_within_10pct/length(df$lfc)
          n_de_genes_within_10pct <- length(which(df$lfc < lfcs$lfc_upper & df$lfc > lfcs$lfc_lower & abs(df$lfc) > 0.1))
          frac_de_genes_within_10pct <- n_de_genes_within_10pct/length(which(abs(df$lfc) > 0.1))
          
          #fraction of DEGs within 25% of true value
          n_genes_within_25pct <- length(which(df$lfc < lfcs$lfc*1.25 & df$lfc > lfcs$lfc*0.75))
          frac_genes_within_25pct <- n_de_genes_within_10pct/length(df$lfc)
          n_de_genes_within_25pct <- length(which(df$lfc < lfcs$lfc*1.25 & df$lfc > lfcs$lfc*0.75 & abs(df$lfc) > 0.1))
          frac_de_genes_within_25pct <- n_de_genes_within_25pct/length(which(abs(df$lfc) > 0.1))
          
          ## compare recovery of DEGs
          
          rtn1 <- data.frame(
            gene = target_gene_name,
            n_cells = n_cells,
            it = gsub(fnm, pattern = paste0(".*cells_|[.]tsv[.]gz"), replacement = ""),
            cor_deg = cor_deg,
            cor_all = cor_all,
            n_deg_unsampled_0.05 = length(degs_0.05),
            n_deg_sampled_0.05 = length(subsampled_degs_0.05),
            tp_0.05 = length(which(subsampled_degs_0.05  %in% degs_0.05)),
            n_deg_unsampled_0.1 = length(degs_0.1),
            n_deg_sampled_0.1=length(subsampled_degs_0.1),
            tp_0.1 = length(which(subsampled_degs_0.1  %in% degs_0.1)),
            n_genes_within_10pct = n_genes_within_10pct,
            n_de_genes_within_10pct = n_de_genes_within_10pct,
            n_genes_within_25pct = n_genes_within_25pct,
            n_de_genes_within_25pct = n_de_genes_within_25pct
          )
          return(rtn1)
        }, SIMPLIFY = F) %>% bind_rows %>% as.data.frame()
        
        fdr_df <- cor_per_cell_no %>%
          mutate(fn_0.05 = n_deg_unsampled_0.05 - tp_0.05,
                 fp_0.05 = n_deg_sampled_0.05 - tp_0.05,
                 tn_0.05 = n_downstream_genes - n_deg_sampled_0.05 - n_deg_unsampled_0.05 + tp_0.05,
                 fn_0.1 = n_deg_unsampled_0.1 - tp_0.1,
                 fp_0.1 = n_deg_sampled_0.1 - tp_0.1,
                 tn_0.1 = n_downstream_genes - n_deg_sampled_0.1 - n_deg_unsampled_0.1 + tp_0.1) %>%
          mutate(
            frac_deg_recovered_0.05 = tp_0.05/(tp_0.05 + fn_0.05),
            frac_false_deg_recovered_0.05 = fp_0.05/(fp_0.05 + tp_0.05),
            frac_deg_recovered_0.1 = tp_0.1/(tp_0.1 + fn_0.1),
            frac_false_deg_recovered_0.1 = fp_0.1/(fp_0.1 + tp_0.1)
          ) %>% select(c("gene", "n_cells", "it", 
                         "n_deg_sampled_0.05", "tp_0.05", 
                         "n_deg_sampled_0.1", "tp_0.1",
                         "n_deg_unsampled_0.05", "frac_deg_recovered_0.05", "frac_false_deg_recovered_0.05", 
                         "n_deg_unsampled_0.1", "frac_deg_recovered_0.1", "frac_false_deg_recovered_0.1"))
        
        
        rtn <- data.frame(
          target = target_gene_name,
          n_cells = n_cells,
          n_deg = length(degs$downstream_gene),
          
          median_sensitivity_0.05 = median(fdr_df$frac_deg_recovered_0.05),
          min_sensitivity_0.05 = quantile(fdr_df$frac_deg_recovered_0.05, 0.1),
          mean_sensitivity_0.05 = median(fdr_df$frac_deg_recovered_0.05),
          max_sensitivity_0.05 = quantile(fdr_df$frac_deg_recovered_0.05, 0.9),
          median_sensitivity_0.1 = median(fdr_df$frac_deg_recovered_0.1),
          min_sensitivity_0.1 = quantile(fdr_df$frac_deg_recovered_0.1, 0.1),
          mean_sensitivity_0.1 = median(fdr_df$frac_deg_recovered_0.1),
          max_sensitivity_0.1 = quantile(fdr_df$frac_deg_recovered_0.1, 0.9),
          
          median_fdr_0.05 = median(fdr_df$frac_false_deg_recovered_0.05),
          min_fdr_0.05 = quantile(fdr_df$frac_false_deg_recovered_0.05, 0.1),
          mean_fdr_0.05 = median(fdr_df$frac_false_deg_recovered_0.05),
          max_fdr_0.05 = quantile(fdr_df$frac_false_deg_recovered_0.05, 0.9),
          median_fdr_0.1 = median(fdr_df$frac_false_deg_recovered_0.1),
          min_fdr_0.1 = quantile(fdr_df$frac_false_deg_recovered_0.1, 0.1),
          mean_fdr_0.1 = median(fdr_df$frac_false_deg_recovered_0.1),
          max_fdr_0.1 = quantile(fdr_df$frac_false_deg_recovered_0.1, 0.9),
          
          min_conf_deg = quantile(cor_per_cell_no$cor_deg, 0.1),
          average_deg = mean(cor_per_cell_no$cor_deg),
          med_deg = median(cor_per_cell_no$cor_deg),
          max_conf_deg = quantile(cor_per_cell_no$cor_deg, 0.9),
          min_conf_all = quantile(cor_per_cell_no$cor_all, 0.1),
          average_all = mean(cor_per_cell_no$cor_all),
          med_all = median(cor_per_cell_no$cor_all),
          max_conf_all = quantile(cor_per_cell_no$cor_all, 0.9)
        )
        return(rtn)
      }, SIMPLIFY = F)
      df4plot <- as.data.frame(bind_rows(df4plot))
      fwrite(df4plot, fnm_out, quote = F, sep = "\t")
    } else{
      df4plot <- fread(fnm_out)
    }
    rtn <- df4plot
    return(rtn)
  }, SIMPLIFY = F)
  
}

## ---- Downsampled_Gene_Metadata
#info on the genes that were subsampled

n_deg <- data.frame(
  target = names(dfs4plot),
  total_cells = target_meta$n_cells[match(names(dfs4plot), target_meta$gene)],
  n_deg = target_meta$n_downstream_excl_target[match(names(dfs4plot), target_meta$gene)]
) %>% mutate(
  n_deg_bin = ifelse(n_deg > 250, "N > 250",
                     ifelse(
                       n_deg > 75, "N > 75",
                       ifelse(
                         n_deg > 50, "N > 50", 
                         ifelse (
                           n_deg > 25, "25 < N < 50", "10 < N < 25"
                         )
                       )
                     ))
)

## ---- Correlation

cor_df <- mapply(1:length(dfs4plot), FUN = function(j){
  df <- dfs4plot[[j]]$med_all
}) %>% as.data.frame()
colnames(cor_df) <- names(dfs4plot)
cor_df <- cor_df %>% mutate(n_cells = 2:50)
cor_df <- melt(cor_df, id.vars = c("n_cells"))
cor_df <- left_join(cor_df, n_deg,
                    by = c("variable" = "target"))

ggplot(cor_df, aes(x = n_cells, y = value,
                   col = n_deg_bin)) + 
  ggtitle("Correlation of Downstream Effect due to Knockdown") + ylab("Correlation") + xlab("Number of Cells") + 
  geom_point() + theme_bw()

## ---- FDR

## large effect 

fdr_df <- mapply(1:length(dfs4plot), FUN = function(j){
  df <- dfs4plot[[j]]$median_fdr_0.05
}) %>% as.data.frame()
colnames(fdr_df) <- names(dfs4plot)
fdr_df <- fdr_df %>% mutate(n_cells = 2:50)
fdr_df <- melt(fdr_df, id.vars = c("n_cells"))
fdr_df$n_total <- target_meta$n_downstream_excl_target[match(fdr_df$variable, target_meta$gene)]

ggplot(fdr_df, aes(x = n_cells, y = value,
                   col = n_total)) + 
  ylab("FDR") + 
  geom_point()

mapply(n_cells = 2:50, FUN = function(i){
  
})






#ggplot(df4plot, aes(x = n_cells, y = med)) + 
#  geom_point()

## Sensitivity (fraction of deg correctly recovered)
pdf(paste0(plotsdir, "/", date, "_fdr.pdf"), width = 9, height = 4)
mapply(target_gene_list, FUN = function(target_gene_name){
  df4plot <- dfs4plot[[target_gene_name]]
  
  p1 <- ggplot(df4plot, aes(x = n_cells, y = median_sensitivity_0.1, ymin = min_sensitivity_0.1, ymax = max_sensitivity_0.1)) + 
    geom_line(color = "navy", size = 2) +
    geom_ribbon(fill = "gray", alpha = 0.5) + 
    geom_vline(xintercept = 10, col = "red", lty = 3) + 
    ylab("Sensitivity (DEG, p-value < 0.1)") + xlab("Number of Cells") + ggtitle(paste0("Subsampling of ", target_gene_name)) + 
    theme_classic()
  
  p2 <- ggplot(df4plot, aes(x = n_cells, y = median_fdr_0.1, ymin = min_fdr_0.1, ymax = max_fdr_0.1)) + 
    geom_line(color = "navy", size = 2) +
    geom_ribbon(fill = "gray", alpha = 0.5) + 
    geom_vline(xintercept = 10, col = "red", lty = 3) + 
    ylab("FDR (DEG, p-value < 0.1)") + xlab("Number of Cells") + ggtitle(paste0("Subsampling of ", target_gene_name)) + 
    theme_classic()
  
  return(ggarrange(p1, p2, nrow = 1))
}, SIMPLIFY = F)
dev.off()

## Sensitivity (fraction of deg correctly recovered)
pdf(paste0(plotsdir, "/", date, "_fdr.pdf"), width = 9, height = 4)
mapply(target_gene_list, FUN = function(target_gene_name){
  df4plot <- dfs4plot[[target_gene_name]]
  
  p1 <- ggplot(df4plot, aes(x = n_cells, y = median_sensitivity_0.05, ymin = min_sensitivity_0.05, ymax = max_sensitivity_0.05)) + 
    geom_line(color = "navy", size = 2) +
    geom_ribbon(fill = "gray", alpha = 0.5) + 
    geom_vline(xintercept = 10, col = "red", lty = 3) + 
    ylab("Sensitivity (DEG, p-value < 0.1)") + xlab("Number of Cells") + ggtitle(paste0("Subsampling of ", target_gene_name)) + 
    theme_classic()
  
  p2 <- ggplot(df4plot, aes(x = n_cells, y = median_fdr_0.05, ymin = min_fdr_0.05, ymax = max_fdr_0.05)) + 
    geom_line(color = "navy", size = 2) +
    geom_ribbon(fill = "gray", alpha = 0.5) + 
    geom_vline(xintercept = 10, col = "red", lty = 3) + 
    ylab("FDR (DEG, p-value < 0.1)") + xlab("Number of Cells") + ggtitle(paste0("Subsampling of ", target_gene_name)) + 
    theme_classic()
  
  return(ggarrange(p1, p2, nrow = 1))
}, SIMPLIFY = F)
dev.off()

df1 <- fread("/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Magpie/8d_power_requirements/by_gene/CNOT1/CNOT1_2022-08-15_17-cells_60.tsv.gz")
df2 <- magpie_lfcs_split %>%
  bind_rows() %>% as.data.frame() %>%
  filter(target == "CNOT1") 
df4plot <- data.frame(
  gene = df1$downstream_gene,
  subsampled = df1$lfc,
  subsampled_pval = df1$pval_lm,
  all = df2$lfc[match(df1$downstream_gene, df2$downstream_gene)],
  all_pval = df2$pval_lm[match(df1$downstream_gene, df2$downstream_gene)]
) %>%
  mutate(
    deg_status = ifelse(all_pval < 0.05 & subsampled_pval < 0.05, "TP", 
                        ifelse(all_pval > 0.05 & subsampled_pval < 0.05, "FP", 
                               ifelse(all_pval < 0.05 & subsampled_pval > 0.05, "FN", "TN")))
  )
ggplot(df4plot, aes(x = all, y = subsampled, col = deg_status)) + 
  ggtitle("CNOT1") + ylab("Subsampled (17 cells)") + 
  geom_point(alpha = 0.4) + theme_bw()
ggplot(df4plot, aes(x = all, y = subsampled)) + 
  ggtitle("CNOT1") + ylab("Subsampled (17 cells)") + 
  geom_point(alpha = 0.4) + theme_bw()


p <- ggplot(df4plot, aes(x = n_cells, y = median_sensitivity_0.05, ymin = min_sensitivity_0.05, ymax = max_sensitivity_0.05)) + 
  geom_line(color = "navy", size = 2) +
  geom_ribbon(fill = "gray", alpha = 0.5) + 
  geom_vline(xintercept = 10, col = "red", lty = 3) + 
  ylab("Sensitivity (DEG, p-value < 0.05)") + xlab("Number of Cells") + ggtitle(paste0("Subsampling of ", target_gene_name)) + 
  theme_classic()
print(p)

p <- ggplot(df4plot, aes(x = n_cells, y = median_fdr_0.05, ymin = min_fdr_0.05, ymax = max_fdr_0.05)) + 
  geom_line(color = "navy", size = 2) +
  geom_ribbon(fill = "gray", alpha = 0.5) + 
  geom_vline(xintercept = 10, col = "red", lty = 3) + 
  ylab("FDR (DEG, p-value < 0.05)") + xlab("Number of Cells") + ggtitle(paste0("Subsampling of ", target_gene_name)) + 
  theme_classic()
print(p)

print(p2)




pdf(paste0(plotsdir, "/", date, "_correlation_simulation_deg_only.pdf"), width = 6, height = 4)
mapply(target_gene_list, FUN = function(target_gene_name){
  df4plot <- dfs4plot[[target_gene_name]]
  p <- ggplot(df4plot, aes(x = n_cells, y = med_deg, ymin = min_conf_deg, ymax = max_conf_deg)) + 
    geom_line(color = "navy", size = 2) +
    geom_ribbon(fill = "gray", alpha = 0.5) + 
    geom_vline(xintercept = 10, col = "red", lty = 3) + 
    ylab("Median Correlation (DEG Only)") + xlab("Number of Cells") + ggtitle(paste0("Subsampling of ", target_gene_name)) + 
    theme_classic()
  print(p)
})
dev.off()

pdf(paste0(plotsdir, "/", date, "_correlation_simulation_all_genes.pdf"), width = 6, height = 4)
mapply(target_gene_list, FUN = function(target_gene_name){
  df4plot <- dfs4plot[[target_gene_name]]
  p <- ggplot(df4plot, aes(x = n_cells, y = med_all, ymin = min_conf_all, ymax = max_conf_all)) + 
    geom_line(color = "navy", size = 2) +
    geom_ribbon(fill = "gray", alpha = 0.5) + 
    geom_vline(xintercept = 10, col = "red", lty = 3) + 
    ylab("Median Correlation (DEG Only)") + xlab("Number of Cells") + ggtitle(paste0("Subsampling of ", target_gene_name)) + 
    theme_classic()
  print(p)
})
dev.off()

pdf(paste0(plotsdir, "/", date, "_correlation_simulation.pdf"), width = 15, height = 12)
df4plot <- as.data.frame(bind_rows(dfs4plot))
ggplot(df4plot, aes(x = n_cells, y = med_all, ymin = min_conf_all, ymax = max_conf_all)) + 
  geom_line(color = "navy", size = 2) +
  geom_ribbon(fill = "gray", alpha = 0.5) + 
  geom_vline(xintercept = 10, col = "red", lty = 3) + 
  facet_wrap(~target) + 
  ylab("Median Correlation (All Genes)") + xlab("Number of Cells") + 
  theme_classic()
ggplot(df4plot, aes(x = n_cells, y = med_deg, ymin = min_conf_deg, ymax = max_conf_deg)) + 
  geom_line(color = "navy", size = 2) +
  geom_ribbon(fill = "gray", alpha = 0.5) + 
  geom_vline(xintercept = 10, col = "red", lty = 3) + 
  facet_wrap(~target) + 
  ylab("Median Correlation (DEG Only)") + xlab("Number of Cells") + 
  theme_classic()
dev.off()

## median correlation for 10 and 25 cells, w.r.t. number of DEG
pdf(paste0(plotsdir, "/", date, "_correlation_vs_n_deg.pdf"), width = 8, height = 8)
for (n in c(10, 25)){
  df4plot <- filter(as.data.frame(bind_rows(dfs4plot)), n_cells == n)
  p1 <- ggplot(df4plot, aes(x = n_deg, y = med_deg)) + 
    xlab("Number of Downstream DE Genes") + ylab("Median Correlation (DEG Only)") + 
    ggtitle(paste0(n, ' cells')) + 
    geom_point() + theme_bw()
  p2 <- ggplot(df4plot, aes(x = n_deg, y = med_all)) + 
    xlab("Number of Downstream DE Genes") + ylab("Median Correlation (All Genes)") + 
    ggtitle('') + 
    geom_point() + theme_bw()
  p3 <- ggplot(df4plot, aes(x = n_deg, y = min_conf_deg)) + 
    xlab("Number of Downstream DE Genes") + ylab("Correlation (10%) (DEG Only)") + 
    geom_point() + theme_bw()
  p4 <- ggplot(df4plot, aes(x = n_deg, y = min_conf_all)) + 
    xlab("Number of Downstream DE Genes") + ylab("Correlation (10%) (All Genes)") + 
    geom_point() + theme_bw()
  p <- ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
  print(p)
}
dev.off()




#### Scratch

lfcs_by_target <- fread(paste0(OutFolder, "/", "6b_calc_lfcs_transcriptome_wide_by_gene", "/", date, "_combined_with_adj_pval.tsv.gz"))
lfcs_by_target <- split.data.frame(lfcs_by_target, f =lfcs_by_target$target)
all_downstream_genes <- lfcs_by_target[[1]]$downstream_gene
n_downstream_genes <- length(all_downstream_genes)

iter <- 1; n_cells <- 40

fnms <- list.files(paste0(outdir, "/by_gene/"), pattern = paste0(n_cells, "-cells_", iter, '[.]tsv[.]gz'), recursive = T)
dfs <- mapply(paste0(outdir, "/by_gene/", fnms), FUN = fread, SIMPLIFY = F)
names(dfs) <- unlist(lapply(strsplit(fnms, "/"), "[[", 1))

plotlist <- mapply(1:length(dfs), FUN = function(i){
  print(i)
  tg <- names(dfs)[i]
  #unsampled <- lfcs_by_target[[tg]]
  #sampled <- dfs[[tg]]
  
  rtn <- data.frame(
    downstream_gene_name = lfcs_by_target[[tg]]$downstream_gene_name,
    lfc_unsampled =lfcs_by_target[[tg]]$lfc,
    lfc_unsampled_pval_adj = p.adjust(lfcs_by_target[[tg]]$pval_lm, method = "BH"),
    lfc_sampled =dfs[[tg]]$lfc,
    lfc_sampled_pval_adj = p.adjust(dfs[[tg]]$pval_lm, method = "BH")
  ) %>%
    mutate(
      err = abs(lfc_sampled - lfc_unsampled),
      pct_err = abs((lfc_sampled - lfc_unsampled)/lfc_unsampled),
      deg_stat = ifelse(lfc_unsampled_pval_adj < 0.05 & lfc_sampled_pval_adj < 0.05, "TP",
                        ifelse(lfc_unsampled_pval_adj < 0.05 & lfc_sampled_pval_adj > 0.05, "FN", 
                               ifelse( lfc_unsampled_pval_adj > 0.05 & lfc_sampled_pval_adj < 0.05, "FP", "TN")))
    )
  p <- ggplot(rtn, aes(x = lfc_sampled, y = lfc_unsampled)) + geom_point(size = 0.4) + ggtitle(tg) + xlab("All") + ylab("Sampled") + theme_bw()
 return(p)
  
}, SIMPLIFY = F)
ggarrange(plotlist = plotlist)












