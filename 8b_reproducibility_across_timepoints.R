
suppressMessages(library(data.table, warn.conflicts = F))
suppressMessages(library(ggplot2, warn.conflicts = F))
suppressMessages(library(dplyr, warn.conflicts = F))
suppressMessages(library(ggpubr, warn.conflicts = F))

section_name <- "8b_reproducibility_across_timepoints"
analysis_name <- ""
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
date <- "2022-08-15"
min_num_de_genes <- 3
min_observed_timepoints <- 2
sig_pval_thresh <- 0.1


args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  ##needed arguments
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  ## not needed
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
  if (arg == "--utils_path"){ utils_path <- args[[ind + 1]] }
  if (arg == "--io_path"){ io_path <- args[[ind + 1]] }
  if (arg == "--min_num_de_genes"){ min_num_de_genes <- as.numeric(args[[ind + 1]]) }
  if (arg == "--min_observed_timepoints"){ min_observed_timepoints <- as.numeric(args[[ind + 1]]) }
  if (arg == "--sig_pval_thresh"){ sig_pval_thresh <- as.numeric(args[[ind + 1]]) }
}



# set relevent i/o paths
source(paste0(HomeFolder, "scripts/io/Magpie_io.R"))
source(paste0(HomeFolder, "scripts/io/Magpie_Utils.R"))
outdir <- file.path(OutFolder, section_name)
plotsdir <- file.path(HTMLFolder, "pipeline/", section_name)
if(!dir.exists(outdir)) dir.create(outdir)
if(!dir.exists(plotsdir)) dir.create(plotsdir)


## ---- LoadData

GuideMetadata <- fread(GuideMetadataPath)
guide_strengths <- c("Moderate", "Strong")
timepoints <- list("Moderate" = c("D3", "D6"),
                   "Strong" = c("D3", "D4", "D5"))
gene_res <- fread(file.path(OutFolder, "6b_calc_lfcs_transcriptome_wide_by_gene", paste0(date, "_combined_with_adj_pval.tsv.gz")))
res_by_timepoint <- mapply(guide_strengths, FUN = function(guide_strength){
  res_by_timepoint <- mapply(timepoints[[guide_strength]], FUN = function(timepoint){
    fnm <- file.path(OutFolder, "6d_calc_lfcs_transcriptome_wide_by_timepoint", paste0(date, "_", guide_strength, "_", timepoint, "_combined_with_adj_pval.tsv.gz"))
    rtn <- fread(fnm)
    rtn$timepoint <- timepoint
    rtn$guide_strength <- guide_strength
    rtn <- as.data.frame(rtn)
  }, SIMPLIFY = F)
}, SIMPLIFY = F)
print("data loaded")

## ---- CheckLFC

print("lfc of target gene...")
target_gene_lfc <- mapply(guide_strengths, FUN = function(guide_strength){
  res_by_timepoint_filtered <- mapply(timepoints[[guide_strength]], FUN = function(timepoint){
    target_gene_lfc <- res_by_timepoint[[guide_strength]][[timepoint]] %>% dplyr::filter(target == downstream_gene_name)
  }, SIMPLIFY = F)
  res_by_timepoint_filtered <- as.data.frame(bind_rows(res_by_timepoint_filtered))
  return(res_by_timepoint_filtered)
}, SIMPLIFY = F)
target_gene_lfc <- as.data.frame(bind_rows(target_gene_lfc))
fwrite(target_gene_lfc, paste0(outdir, "/", date, "_target_gene_lfc_by_dat.tsv"), sep = '\t')

pdf(paste0(plotsdir, "/", date, "_target_gene_lfc_by_day_scatter.pdf"), width = 8, height = 4)
silencer <- mapply(guide_strengths, FUN = function(g){
  p <- ggplot(target_gene_lfc %>% dplyr::filter(guide_strength == g), aes(x = lfc, y = -log10(pval_adj), col = guide_strength)) + 
    facet_wrap(~timepoint) + 
    ylab("Log10 Adjusted P-value") + xlab("LFC (Base 10)") + 
    scale_color_manual(values = cols4Pools) + 
    geom_point() + theme_bw()
  print(p)
})
dev.off()

pdf(paste0(plotsdir, "/", date, "_target_gene_lfc_by_day.pdf"), width = 8, height = 8)
for (gs in guide_strengths){
  df4plot <- target_gene_lfc %>% dplyr::filter(guide_strength == gs)
  tps <- sort(as.numeric(gsub(unique(df4plot$timepoint), pattern = "D", replacement = "")))
  timepoints_combn <- combn(tps,2)
  if (dim(timepoints_combn)[2] == 1){
    t <- 1
    p <- ggboxplot(df4plot %>% dplyr::filter(timepoint %in% paste0("D", c(timepoints_combn[1,t], timepoints_combn[2,t]))), 
                   x = "timepoint", y = "lfc", fill = "guide_strength",
                   line.color = "gray", line.size = 0.4) + 
      xlab("Sequencing Timepoint") + ylab("LFC (Base 10) of Target Gene") + 
      scale_fill_manual(values = cols4Pools) + theme_bw()
    p <- p + stat_compare_means()
    print(p)
  } else {
    for (t in 1:(dim(timepoints_combn)[2])){
      p <- ggboxplot(df4plot %>% dplyr::filter(timepoint %in% paste0("D", c(timepoints_combn[1,2], timepoints_combn[2,t]))), 
                     x = "timepoint", y = "lfc", fill = "guide_strength",
                     line.color = "gray", line.size = 0.4) + 
        xlab("Sequencing Timepoint") + ylab("LFC (Base 10) of Target Gene") + 
        scale_fill_manual(values = cols4Pools) + theme_bw()
      p <- p + stat_compare_means()
      print(p)
    }
  }
  
  
}
dev.off()


## ---- PlotConsistency

#df4plot <- mapply(res_by_timepoint)

#df4plot <- as.data.frame(bind_rows(bind_rows(res_by_timepoint)))
#df4plot <- target_lfc_by_guide %>% select(c("target", "guide", "lfc"))
#df4plot$mean_target_lfc <- target_lfc_by_gene$lfc[match(df4plot$target, target_lfc_by_gene$target)]
#df4plot$pool <- GuideMetadata$group[match(df4plot$target, GuideMetadata$gene)]

#p1 <- ggplot(df4plot, aes(x = mean_target_lfc, y = lfc, col = pool)) + 
#  xlab("Mean LFC of Target Gene") + ylab("LFC of Target Gene") + ggtitle("Consistency of Target Gene LFC Across Guides") + 
#  geom_point() + scale_color_manual(values = cols4Pools) + theme_bw()
#ggsave(p1, filename = paste0(plotsdir, "/", date, "_consistency_across_timepoints.pdf"), width = 7, height = 6)


## ---- Correlation of Downstream Effects

print("checking correlation across timepoints")
cor_across_timepoints <- mapply(guide_strengths, FUN = function(guide_strength){
  target_genes <- GuideMetadata %>% filter(group == guide_strength) %>% .$gene %>% unique()
  lfc_res <- as.data.frame(bind_rows(res_by_timepoint[[guide_strength]]))
  cor_across_timepoints <- mapply(target_genes, FUN = function(target_gene){
    print(target_gene)
    lfc_by_target <- lfc_res %>% filter(target == target_gene) %>% select(c("target", 
                                                                            "guide_strength", "timepoint",
                                                                            "downstream_gene_name", "n_perturbed",
                                                                            "lfc", "pval_adj"))
    timepoints_present <- unique(lfc_by_target$timepoint)
    sig_de_genes <- unique(lfc_by_target %>% filter(pval_adj < sig_pval_thresh) %>% .$downstream_gene_name)
    if (length(timepoints_present) > 1){
      
      if (length(sig_de_genes) > 2){
        df4cor <- lfc_by_target %>% 
          filter(downstream_gene_name %in% sig_de_genes) %>%
          select(c("downstream_gene_name", "timepoint", "lfc")) %>% 
          reshape2::dcast(downstream_gene_name ~ timepoint, value.var = "lfc")
        cors2calc <- combn(1:length(timepoints_present), 2)
        cors <- mapply(1:dim(cors2calc)[2], FUN = function(i){
          cor(df4cor[, cors2calc[1,i]+1], df4cor[, cors2calc[2,i]+1])
        })
        
        rtn <- data.frame(
          target_gene = target_gene,
          num_timepoints = length(timepoints_present),
          num_sig_de_genes = length(sig_de_genes),
          min_cor = min(cors),
          med_cor = median(cors),
          max_cor = max(cors)
        )
      } else {
        
        rtn <- data.frame(
          target_gene = target_gene,
          num_timepoints = length(timepoints_present),
          num_sig_de_genes = length(sig_de_genes),
          min_cor = 0,
          med_cor = 0,
          max_cor = 0
        )
        
      }
      
       
    }  else {
      rtn <- data.frame(
        target_gene = target_gene,
        num_timepoints = length(timepoints_present),
        num_sig_de_genes = length(sig_de_genes),
        min_cor = 0,
        med_cor = 0,
        max_cor = 0
      )
    }
    
    return(rtn)
    
  }, SIMPLIFY = F)
  
  cor_across_timepoints <- as.data.frame(bind_rows(cor_across_timepoints))
  
  return(cor_across_timepoints)
  
}, SIMPLIFY = F)

pdf(paste0(plotsdir, "/", date, "_timepoint-cor-density.pdf"), width = 12, height = 5)
silencer <- mapply(guide_strengths, FUN = function(guide_strength){
  df4plot <- cor_across_timepoints[[guide_strength]]
  df4plot <- as.data.frame(bind_rows(df4plot))
  df4plot <- df4plot %>% filter(min_cor != 0)
  df4plot$guide_strength <- guide_strength
  p1 <- ggplot(df4plot, aes(x = num_sig_de_genes, y = med_cor, col = guide_strength)) + 
    scale_x_log10() +
    xlab("Number of Downstream DE Genes") + ylab("Correlation Across Timepoints") + 
    geom_point() + scale_color_manual(values = cols4Pools) + theme_bw()
  df4plot <- df4plot %>% select(c("target_gene", "min_cor", "max_cor"))
  df4plot <- reshape2::melt(df4plot, id = "target_gene")
  p2 <- ggplot(df4plot, aes(x = value, fill = variable)) + 
    xlab("Correlation Between Timepoints") + ylab("Density") + 
    geom_density(alpha = 0.4) + scale_color_manual(values = cols4Pools) + theme_bw()
  print(ggarrange(p1, p2))
})
dev.off()

## ---- Number of DE Genes

pdf(paste0(plotsdir, "/", date, "_num_de_genes_cor.pdf"), width = 12, height = 5)
target_meta <- mapply(guide_strengths, FUN = function(guide_strength){
  ts <- timepoints[[guide_strength]]
  target_genes <- sort(unique(GuideMetadata %>% filter(group == guide_strength) %>% .$gene))
  num_de_genes <- mapply(ts, FUN = function(t){
    res_df <- res_by_timepoint[[guide_strength]][[t]]
    res_df.split <- split.data.frame(res_df, f = res_df$target)
    num_de_genes <- mapply(res_df.split, FUN = function(x){
      length(which(x$pval_adj < 0.05))
    })
    rtn <- data.frame(
      gene = target_genes,
      num_de_genes = num_de_genes[match(target_genes, names(num_de_genes))],
      timepoint = t,
      guide_strength = guide_strength,
      guide_subgroup = GuideMetadata$subgroup[match(target_genes, GuideMetadata$gene)]
    )
  }, SIMPLIFY = F)
  num_de_genes <- as.data.frame(bind_rows(num_de_genes))
  
  if (guide_strength == "Strong"){
    df4plot <- reshape2::dcast(num_de_genes, gene + guide_strength + guide_subgroup ~  timepoint, value.var = "num_de_genes")
    p <- ggplot(df4plot, aes(x = D3, y = D5, col = guide_subgroup, 
                             label = ifelse(D3 > 10 | D5 > 10, gene, ""))) + 
      ggrepel::geom_label_repel() + 
      ggtitle("Number of DE Genes") + 
      geom_point() + theme_bw()
    print(p)
  }
  
  return(num_de_genes)
}, SIMPLIFY = F)
dev.off()
saveRDS(target_meta, paste0(outdir, "/", date, "_target_meta.RDS"))

## ---- SessionInfo()

sessionInfo()



## ---- DataPreprocessing

#target_lfc_by_guide <- guide_res %>% filter(target == downstream_gene_name)
#target_lfc_by_gene <- gene_res %>% filter(target == downstream_gene_name)

#df4plot <- target_lfc_by_guide %>% select(c("target", "guide", "lfc"))
#df4plot$mean_target_lfc <- target_lfc_by_gene$lfc[match(df4plot$target, target_lfc_by_gene$target)]
#df4plot$pool <- GuideMetadata$group[match(df4plot$target, GuideMetadata$gene)]

