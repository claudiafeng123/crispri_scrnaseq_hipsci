suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(doParallel))
suppressMessages(library(umap))

# set relevent i/o paths
section_name <- "14b_variance_due_to_guide"


# set relevent i/o paths
date <- "2022-10-05"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
ExperimentName <- "Pica"

# options for testing
io_path <- paste0("scripts/io/", ExperimentName, "_io.R")
utils_path <- "scripts/io/Magpie_Utils.R"
REDO <- F


args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  ##needed arguments
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--experiment_name"){ ExperimentName <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--utils_path"){ utils_path <- args[[ind + 1]] }
  if (arg == "--io_path"){ io_path <- args[[ind + 1]] }
  if (arg == "--rnd_seed"){ rnd_seed <- as.numeric(args[[ind + 1]]) }
  if (arg == "--REDO"){ REDO <- as.logical(args[[ind + 1]]) }
}

print(paste0(section_name))
print(paste("date:", date))
print(paste("home folder:", HomeFolder))
if(exists("donors2include")){print('including donors: '); print(donors2include)}



setwd(HomeFolder)
source(io_path)
source(utils_path)

GuideMetadata <- fread(GuideMetadataPath)
LineMetadata <- fread(LineMetadataPath)

outdir <- file.path(file.path(OutFolder, section_name))
if(!dir.exists(outdir)) dir.create(outdir, recursive = T)
plotsdir <- file.path(file.path(HTMLFolder, 'pipeline', section_name))
if(!dir.exists(file.path(plotsdir))) dir.create(file.path(plotsdir), recursive = T)


## ---- LoadData

## cell metadata (for guides)
most_recent <- sort(gsub(
  list.files(file.path(OutFolder, "5a_combine"), pattern = "perturbed_cells_metadata.tsv", full.names = T),
  pattern = paste0(".*5a_combine[/]|_.*"),
  replacement = ""),
  decreasing = T)[1]
cell_metadata_fnms <- grep(list.files(file.path(OutFolder, "5a_combine"), pattern = "perturbed_cells_metadata.tsv", full.names = T),
                           pattern = most_recent, value = T)
cell_metadata <- lapply(cell_metadata_fnms, FUN = fread) %>%
  bind_rows() %>% as.data.frame() %>%
  dplyr::select(c("V1", "guide" = 'feature_call'))

mean_lfcs <- fread("/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Pica/6b_calc_lfcs_transcriptome_wide_by_gene/2022-10-05_combined_with_adj_pval.tsv.gz")
targets_with_effect <- mean_lfcs %>%
  dplyr::filter(target != downstream_gene_name & pval_adj < sig_pval_thresh & abs(lfc) > 0.25) %>%
  .$target %>% table() %>% as.data.frame()

## ---- Choose Genes

## gonna be slightly naughty here and use section 15
heritability_df <- fread("/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Pica/15c_calc_var_expl_by_target/15c_05_calc_var_expl_by_target_check_pvals/2024-01-10_heritability_with_pval_adj.tsv.gz")
var_expl_df <- heritability_df %>%
  dplyr::select(c("target", 'downstream_gene_name', contains('varExpl')))
pairs2plot <- var_expl_df %>%
  dplyr::filter(lfc_varExpl_guide > 0.3)
## subset this to targets with lots of downstream effect


pairs2plot <- split.data.frame(pairs2plot, f = pairs2plot$target)

## ---- MakeJitters


targets2plot <- sort(unique(names(pairs2plot)))

registerDoParallel(5)
pdf(paste0(plotsdir, '/', date, '_all_guide_effects.pdf'), width = 4, height = 5)
foreach (tg = targets2plot) %dopar% {
  
  guide_dependent_downstream <- pairs2plot[[tg]]
  kd_expr <- get_knockdown_expression(target_gene = tg)
  kd_expr <- left_join(kd_expr,cell_metadata, by = c("id" = "V1"))
  
  #tg <- targets2plot[1]
  #d <- guide_dependent_downstream$downstream_gene_name[1]
  plotlist <- mapply(d = guide_dependent_downstream$downstream_gene_name, FUN = function(d){
    dgs <- c(d, tg)
    df4lmm <- mapply(dg = dgs, FUN = function(dg){
      df4lmm1 <- get_df4lmm(kd_expr = kd_expr, downstream_gene_name = dg, include_guide = T)
      df4lmm1 <- left_join(df4lmm1, cell_metadata %>% dplyr::select(c("id" = "V1", 'Guide' = 'guide')))
      return(df4lmm1)
    }, SIMPLIFY = F) %>% bind_rows() %>% as.data.frame()
    
    df4lmm$Guide <- factor(df4lmm$Guide)
    levels(df4lmm$Guide) <- paste0("Guide ", 1:length(levels(df4lmm$Guide)))
    p <- ggplot(df4lmm %>% dplyr::filter(downstream_gene_name == tg), aes(x = Guide, y = lfc, col = Guide)) + 
      xlab('') +ylab("LFC") + ggtitle(paste0(d, ' expression in ', tg, ' knockdowns')) + 
      geom_jitter(size = 1, col = 'gray', alpha = 1) + 
      geom_jitter(data = df4lmm %>% dplyr::filter(downstream_gene_name == d), 
                  aes(col = Guide), size = 0.6, alpha = 0.4) + 
      geom_hline(yintercept = 0, col = 'red', lty = 3) + 
      theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = 'top') +
      guides(col = guide_legend(override.aes = list(alpha = 1, size =4),
                                nrow =1, title = ''))
    return(p)
  }, SIMPLIFY = F)
  
  
  #pdf(paste0(plotsdir, '/', date, '_', tg, '_guide_effects.pdf'), width = 4, height = 5)
  print(plotlist)
  #dev.off()
  
  return(plotlist)
}
dev.off()


## ---- HighEffects
# do more detailed plots

targets2plot <- sort(unique(names(pairs2plot)))

registerDoParallel(5)
#pdf(paste0(plotsdir, '/', date, '_all_guide_effects_vs_high_effects.pdf'), width = 4, height = 5)
foreach (tg = targets2plot) %dopar% {
  
  guide_dependent_downstream <- pairs2plot[[tg]]
  kd_expr <- get_knockdown_expression(target_gene = tg)
  kd_expr <- left_join(kd_expr,cell_metadata, by = c("id" = "V1"))
  
  #tg <- targets2plot[1]
  #d <- guide_dependent_downstream$downstream_gene_name[1]
  high_effects <- mean_lfcs %>% dplyr::filter(target == tg & pval_adj < sig_pval_thresh & target != downstream_gene_name & abs(lfc) > 0.1) %>% arrange(desc(abs(lfc)))  %>% .$downstream_gene_name
  downstreams2consider <- unique(c(guide_dependent_downstream$downstream_gene_name, high_effects))
  plotlist <- mapply(d = downstreams2consider, FUN = function(d){
    dgs <- c(d, tg)
    df4lmm <- mapply(dg = dgs, FUN = function(dg){
      df4lmm1 <- get_df4lmm(kd_expr = kd_expr, downstream_gene_name = dg, include_guide = T)
      df4lmm1 <- left_join(df4lmm1, cell_metadata %>% dplyr::select(c("id" = "V1", 'Guide' = 'guide')))
      return(df4lmm1)
    }, SIMPLIFY = F) %>% bind_rows() %>% as.data.frame()
    
    df4lmm$Guide <- factor(df4lmm$Guide)
    levels(df4lmm$Guide) <- paste0("Guide ", 1:length(levels(df4lmm$Guide)))
    
    ## expression
    p_jitter <- ggplot(df4lmm, aes(x = guide, y = lfc, col = cell_line)) + 
      xlab('') +ylab("LFC") + ggtitle(paste0(d, ' expression in ', tg, ' knockdowns')) + 
      facet_wrap(~downstream_gene_name) + 
      ggtitle(ifelse(d %in% guide_dependent_downstream$downstream_gene_name, paste0(d, "*"), d)) + 
      geom_jitter(size = 0.6) + 
      geom_hline(yintercept = 0, col = 'red', lty = 3) + 
      theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = 'none')# +
      
    ## breakdown variance explained
    df4plot <- var_expl_df %>%
      dplyr::filter(target == tg & downstream_gene_name == d) %>%
      reshape2::melt(id = c("target", "downstream_gene_name")) %>%
      dplyr::filter(variable != 'lfc_varExpl_Residual') %>%
      mutate(variable= gsub(variable, pattern = 'lfc_varExpl_', replacement = ''))
    p_var_expl <- ggplot(df4plot, aes(y = variable, x = value)) + 
      xlab("Variance explained") + ylab("") + 
      geom_bar(stat = 'identity', fill = 'navy') + theme_bw()
    rtn <- cowplot::plot_grid(p_jitter, p_var_expl, ncol = 1)
    #rtn
    return(rtn)
  }, SIMPLIFY = F)
  
  pdf(paste0(plotsdir, '/', date, '_', tg, '_guide_effects_with_high_effects_as_control.pdf'), width = 6, height = 9)
  print(plotlist)
  dev.off()
  
  return(plotlist)
}
#dev.off()


## ---- Find Repeated ones

pairs_of_interest <- table(as.data.frame(bind_rows(pairs2plot)) %>% .$target)
pairs_of_interest <- names(pairs_of_interest)[which(pairs_of_interest > 2)]
