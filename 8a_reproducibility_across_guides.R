
suppressMessages(library(data.table, warn.conflicts = F))
suppressMessages(library(ggplot2, warn.conflicts = F))
suppressMessages(library(dplyr, warn.conflicts = F))
suppressMessages(library(doParallel, warn.conflicts = F))
suppressMessages(library(ggpubr, warn.conflicts = F))

section_name <- "8a_reproducibility_across_guides"
analysis_name <- ""
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
date <- "2024-05-17"
min_num_de_genes <- 3
min_observed_timepoints <- 2
sig_pval_thresh <- 0.1
REDO <- F
ExperimentName <- "Pica"


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
source(paste0(HomeFolder, "scripts/io/", ExperimentName, "_io.R"))
source(paste0(HomeFolder, "scripts/io/Magpie_Utils.R"))
outdir <- file.path(OutFolder, section_name)
plotsdir <- file.path(HTMLFolder, "pipeline/", section_name)
if(!dir.exists(outdir)) dir.create(outdir)
if(!dir.exists(plotsdir)) dir.create(plotsdir)
if(!dir.exists(paste0(plotsdir, "/cor_plots/"))) dir.create(paste0(plotsdir, "/cor_plots/"))

## ---- LoadData

gene_res <- fread(file.path(OutFolder, "6b_calc_lfcs_transcriptome_wide_by_gene", paste0(date, "_combined_with_adj_pval.tsv.gz")))
guide_res <- fread(file.path(OutFolder, "6c_calc_lfcs_transcriptome_wide_by_guide", paste0(date, "_combined_with_adj_pval.tsv.gz")))
GuideMetadata <- fread(GuideMetadataPath)
target_genes <- unique(GuideMetadata$gene)
fnm <- sort(grep(list.files(file.path(OutFolder, "7b_target_summary"), pattern = 'target_meta_data.csv', full.name =T),
                 pattern = 'params', invert = T, value = T),
            decreasing = T)[1]
target_meta <- fread(fnm)
gene_res_split <- split.data.frame(gene_res, f = gene_res$target)
guide_res_split <- split.data.frame(guide_res, f = guide_res$target)

## ---- PlotConsistency
#consistency of target down-regulation across guides

target_lfc_by_guide <- guide_res %>% filter(target == downstream_gene_name)
target_lfc_by_gene <- gene_res %>% filter(target == downstream_gene_name)

df4plot <- target_lfc_by_guide %>% select(c("target", "guide", "lfc"))
df4plot$mean_target_lfc <- target_lfc_by_gene$lfc[match(df4plot$target, target_lfc_by_gene$target)]
df4plot$pool <- GuideMetadata$group[match(df4plot$target, GuideMetadata$gene)]

p1 <- ggplot(df4plot, aes(x = mean_target_lfc, y = lfc, col = pool)) + 
  xlab("Mean LFC of Target Gene") + ylab("LFC of Target Gene") + ggtitle("Consistency of Target Gene LFC Across Guides") + 
  geom_point() + scale_color_manual(values = cols4Pools) + theme_bw()
ggsave(p1, filename = paste0(plotsdir, "/", date, "_target_gene_de_consistency_across_guides.pdf"), width = 7, height = 6)


## ---- Correlation of Downstream Effects

registerDoParallel(20)
cor_across_guides <- foreach(tg = names(guide_res_split)) %dopar% {
  res_by_guide <- split.data.frame(guide_res_split[[tg]], f= guide_res_split[[tg]]$guide)
  if (length(res_by_guide) == 1){
    rtn <- NULL
  } else if (length(res_by_guide) == 2){
    rtn <- calc_cor(df1 = res_by_guide[[1]], df2 = res_by_guide[[2]])$df
    rtn$guide_1 <- c(res_by_guide[[1]]$guide[1])
    rtn$guide_2 <- c(res_by_guide[[2]]$guide[1])
  } else if (length(res_by_guide) == 3){
    rtn <- bind_rows(
      calc_cor(df1 = res_by_guide[[1]], df2 = res_by_guide[[2]])$df,
      calc_cor(df1 = res_by_guide[[2]], df2 = res_by_guide[[3]])$df,
      calc_cor(df1 = res_by_guide[[3]], df2 = res_by_guide[[1]])$df
    ) %>% as.data.frame()
    rtn$guide_1 <- c(res_by_guide[[1]]$guide[1], res_by_guide[[2]]$guide[1], res_by_guide[[3]]$guide[1])
    rtn$guide_2 <- c(res_by_guide[[2]]$guide[1], res_by_guide[[3]]$guide[1], res_by_guide[[1]]$guide[1])
  } else {
    rtn <- NULL
  }
  return(rtn)
} %>% bind_rows() %>% as.data.frame()
fwrite(cor_across_guides, file.path(outdir, paste0(date, "_correlation_across_guides.tsv.gz")),
       sep = "\t", compress = "gzip")

cor_across_guides <- mapply(target_genes, FUN = function(target_gene){
  #print(target_gene)
  
  lfc_by_target <- gene_res %>% filter(target == target_gene) %>% select(c("target", 
                                                                          "downstream_gene_name", "n_perturbed",
                                                                          "lfc", "pval_adj"))
  lfc_by_guide <- guide_res %>% filter(target == target_gene) %>% select(c("target", "guide", 
                                                                           "downstream_gene_name", "n_perturbed",
                                                                           "lfc", "pval_adj"))
  guides_present <- unique(lfc_by_guide$guide)
  sig_de_genes <- unique(lfc_by_target %>% filter(pval_adj < sig_pval_thresh) %>% .$downstream_gene_name)
  sig_de_genes_excl_target <- unique(lfc_by_target %>% filter(pval_adj < sig_pval_thresh & downstream_gene_name != target_gene) %>% .$downstream_gene_name)
  if (length(guides_present) > 1 & dim(lfc_by_target)[1] > 0){
    
    if (length(sig_de_genes) > 2){
      df4cor <- lfc_by_guide %>% 
        filter(downstream_gene_name %in% sig_de_genes) %>%
        select(c("downstream_gene_name", "guide", "lfc")) %>% 
        reshape2::dcast(downstream_gene_name ~ guide, value.var = "lfc")
      cors2calc <- combn(1:length(guides_present), 2)
      cors <- mapply(1:dim(cors2calc)[2], FUN = function(i){
        cor(df4cor[, cors2calc[1,i]+1], df4cor[, cors2calc[2,i]+1])
      })
      
      rtn <- data.frame(
        target_gene = target_gene,
        num_guides = length(guides_present),
        num_sig_de_genes = length(sig_de_genes),
        num_sig_de_genes_excl_target = length(sig_de_genes_excl_target),
        min_cor = min(cors),
        med_cor = median(cors),
        max_cor = max(cors)
      )
      
      
      # plot
      if (length(sig_de_genes) > 25){
        
        df4plot <- lfc_by_guide %>% 
          select(c("downstream_gene_name", "guide", "lfc")) %>% 
          reshape2::dcast(downstream_gene_name ~ guide, value.var = "lfc") %>% 
          mutate(is_de = as.character(downstream_gene_name %in% sig_de_genes)) %>%
          arrange(by_group = is_de)
        
        pdf(paste0(plotsdir, "/cor_plots/", target_gene, "_cor.pdf"), width = 8, height = 8)
        print(GGally::ggpairs(df4plot, columns = grep(colnames(df4plot), pattern = target_gene, value = T), 
                        mapping = ggplot2::aes(color = is_de, alpha = 0.3)) + theme_bw() + 
          scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) + scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "grey")))
        dev.off()
        
      }
      
    } else {
      rtn <- data.frame(
        target_gene = target_gene,
        num_guides = length(guides_present),
        num_sig_de_genes = length(sig_de_genes),
        num_sig_de_genes_excl_target = length(sig_de_genes_excl_target),
        min_cor = 0,
        med_cor = 0,
        max_cor = 0
      )
    }
    
    
  }  else {
    rtn <- data.frame(
      target_gene = target_gene,
      num_guides = length(guides_present),
      num_sig_de_genes = length(sig_de_genes),
      num_sig_de_genes_excl_target = length(sig_de_genes_excl_target),
      min_cor = 0,
      med_cor = 0,
      max_cor = 0
    )
  }
  
  return(rtn)
  
}, SIMPLIFY = F)
cor_across_guides <- as.data.frame(bind_rows(cor_across_guides))
fwrite(cor_across_guides, file.path(outdir, paste0(date, "_correlation_across_guides.tsv.gz")),
       sep = "\t", compress = "gzip")



## reformat by gene

cor_across_guides_by_gene <- cor_across_guides %>%
  group_by(gene = gene_1) %>%
  summarize(
    n_guides = length(unique(c(guide_1, guide_2))),
    n_cells = sum(n_cells_1, n_cells_2)/(length(unique(c(guide_1, guide_2)))-1),
    min_cells = min(c(n_cells_1, n_cells_2)),
    max_cells = max(c(n_cells_1, n_cells_2)),
    min_cor_all = min(cor_all, na.rm = T),
    median_cor_all = median(cor_all, na.rm = T),
    max_cor_all = max(cor_all, na.rm = T),
    min_cor_deg = min(cor_all_deg, na.rm = T),
    median_cor_deg = median(cor_all_deg, na.rm = T),
    max_cor_deg = max(cor_all_deg, na.rm = T),
    min_cor_all_deg_excl_target = min(cor_all_deg_excl_target, na.rm = T),
    median_cor_all_deg_excl_target = median(cor_all_deg_excl_target, na.rm = T),
    max_cor_all_deg_excl_target = max(cor_all_deg_excl_target, na.rm = T),
    min_cor_common_deg = min(cor_common_deg, na.rm = T),
    median_cor_common_deg = median(cor_common_deg, na.rm = T),
    max_cor_common_deg = max(cor_common_deg, na.rm = T),
    min_cor_common_deg_excl_target = min(cor_common_deg_excl_target, na.rm = T),
    median_cor_common_deg_excl_target = median(cor_common_deg_excl_target, na.rm = T),
    max_cor_common_deg_excl_target = max(cor_common_deg_excl_target, na.rm = T)
  ) %>%
  mutate(
    n_deg_excl_target = target_meta$n_downstream_excl_target[match(gene, target_meta$gene)]
  )

if (ExperimentName == "Magpie"){
  cor_across_guides_by_gene <- cor_across_guides_by_gene %>%
    mutate(guide_strength = GuideMetadata$group[match(gene, GuideMetadata$gene)]) %>%
  mutate(guide_strength = ifelse(guide_strength == "Strong", "Fitness Genes", "Non-Fitness Genes"))
}
fwrite(cor_across_guides_by_gene, file.path(outdir, paste0(date, "_correlation_across_guides_by_gene.tsv.gz")),
       sep = "\t", compress = "gzip")


## ---- Plots

cor_across_guides_by_gene <- fread(file.path(outdir, paste0(date, "_correlation_across_guides_by_gene.tsv.gz")))
cor_across_guides <- fread(file.path(outdir, paste0(date, "_correlation_across_guides.tsv.gz")))

plot_guide_cor <- function(target, guide_1, guide_2){
  n_cells_1 <- guide_res_split[[target]] %>%
    filter(guide == guide_1) %>% .$n_perturbed %>% unique()
  n_cells_2 <- guide_res_split[[target]] %>%
    filter(guide == guide_2) %>% .$n_perturbed %>% unique()
  degs <- unique(c(guide_res_split[[target]] %>%
                      filter(guide == guide_1 & pval_adj < sig_pval_thresh) %>% .$downstream_gene_name,
                    guide_res_split[[target]] %>%
                      filter(guide == guide_2 & pval_adj < sig_pval_thresh) %>% .$downstream_gene_name
                    ))
  df4plot <- data.frame(
    downstream_gene_name = guide_res_split[[1]]$downstream_gene_name,
    lfc_1 = guide_res_split[[target]] %>%
      filter(guide == guide_1) %>% .$lfc,
    lfc_2 = guide_res_split[[target]] %>%
      filter(guide == guide_2) %>% .$lfc,
    is_deg = ifelse(guide_res_split[[1]]$downstream_gene_name %in% degs, "DEG", "not DEG")
  )
  
  rtn <- ggplot(df4plot, aes(x = lfc_1, y = lfc_2,
                      col = is_deg)) + 
    scale_color_manual(values = c("DEG" = 'red', 'not DEG' = 'gray')) + 
    xlab(paste0("Guide #1, ", n_cells_1, " cells")) + ylab(paste0("Guide #2, ", n_cells_2, " cells")) + ggtitle(target) + 
    xlim(c(min(c(df4plot$lfc_1, df4plot$lfc_2)), max(c(df4plot$lfc_1, df4plot$lfc_2)) )) + ylim(c(min(c(df4plot$lfc_1, df4plot$lfc_2)), max(c(df4plot$lfc_1, df4plot$lfc_2)) )) + 
    geom_point(size = 0.4, alpha = 0.4) + theme_bw() + theme(legend.position = 'none', aspect.ratio=1)
  return(rtn)
}

rmarkdown::render(input = file.path(CodeFolder, ExperimentName, "pipeline", paste0(section_name, ".Rmd")), 
                  output_file = file.path(plotsdir, paste0(date, '_', section_name, '.html')), 
                  params = list(
                    home_folder = HomeFolder, 
                    section_name = section_name,
                    date =date,
                    experiment_name = ExperimentName
                  ))


## ---- SessionInfo()

sessionInfo()

