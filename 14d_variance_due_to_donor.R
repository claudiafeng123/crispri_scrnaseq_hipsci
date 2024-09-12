suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(doParallel))
suppressMessages(library(umap))
suppressMessages(library(ggpubr))

# set relevent i/o paths
section_name <- "14c_variance_due_to_line"


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

heritability_df <- fread("/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Pica/15c_calc_var_expl_by_target/15c_05_calc_var_expl_by_target_check_pvals/2024-01-10_heritability_with_pval_adj.tsv.gz")
var_expl_df <- heritability_df %>%
  dplyr::select(c("target", 'downstream_gene_name', contains('varExpl')))
pairs2plot <- var_expl_df %>%
  dplyr::filter(lfc_varExpl_cell_line > 0.05) %>%
  arrange(desc(lfc_varExpl_cell_line))

pica_lfcs <- fread("/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Pica/6f_calc_lfcs_transcriptome_wide_by_gene_per_line/2022-10-05_all_lines_with_adj_pval.tsv.gz")
#pica_lfcs_split <- pica_lfcs %>% mutate(to_match = paste0(target, "_", downstream_gene_name) )
#pica_lfcs_split <- split.data.frame(pica_lfcs_split, f = pica_lfcs_split$to_match)
pica_lfcs_by_target <- fread("/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Pica/6b_calc_lfcs_transcriptome_wide_by_gene/2022-10-05_combined_with_adj_pval.tsv.gz")

## ---- The Rough Estimates

pica_lfcs_subset <- pica_lfcs %>% 
  dplyr::filter(paste0(target, "_", downstream_gene_name) %in% paste0(heritability_df$target, "_", heritability_df$downstream_gene_name)) %>%
  mutate(donor =unlist(lapply(strsplit(cell_line, "_"), "[[", 1))) %>%
  mutate(to_match = paste0(target, "_", downstream_gene_name) ) %>%
  dplyr::select(c("to_match", "donor", 'cell_line', 'lfc'))
pica_lfcs_subset <- split.data.frame(pica_lfcs_subset, f = pica_lfcs_subset$to_match )

var_expl_donor <- lapply(1:length(pica_lfcs_subset), FUN = function(i){
  print(i)
  df <- pica_lfcs_subset[[i]]
  fit <- tryCatch(expr = lme4::lmer(lfc ~ (1|donor), df), e = NA)
  rtn <- ifelse(!(is.na(fit)), variancePartition::calcVarPart(fit, showWarnings = ), NA)
  return(rtn)
})


## ---- MakePlots

pica_lfcs_subset <- pica_lfcs %>% dplyr::filter(paste0(target, "_", downstream_gene_name) %in% paste0(pairs2plot$target, "_", pairs2plot$downstream_gene_name) | (target == downstream_gene_name & target %in% pairs2plot$target)) %>%
  mutate(donor =unlist(lapply(strsplit(cell_line, "_"), "[[", 1)))

## Take a look at why these are a thing
tg <- 'GON4L'; dg <- 'H1F0'
plotlist <- mapply(1:dim(pairs2plot)[1], FUN = function(i){
  tg <- pairs2plot$target[i]; dg <- pairs2plot$downstream_gene_name[i]
  df4plot <- pica_lfcs_subset %>%
    dplyr::filter(target == tg & downstream_gene_name == dg)
  p1 <- ggplot(df4plot, aes(x=control_norm_expr + control_line_coef, y = lfc, 
                      col = donor, label = cell_line)) +
    ggrepel::geom_text_repel() + 
    xlab("Control Expression") + ylab("LFC") + ggtitle(paste0(dg, " Expression in ", tg," Knockdowns")) +
    geom_point() + theme_bw() + theme(legend.position = 'none')
  p2 <- ggplot(df4plot, aes(x=control_norm_expr + control_line_coef, y = control_norm_expr + control_line_coef+lfc, 
                            col = donor, label = cell_line)) +
    ggrepel::geom_text_repel() + 
    xlab("Control Expression") + ylab("Post-knockdown Expression") + ggtitle(paste0(dg, " Expression in ", tg," Knockdowns")) +
    geom_point() + theme_bw() + theme(legend.position = 'none')
  df4plot <- pica_lfcs_subset %>%
    dplyr::filter(target == tg & downstream_gene_name == tg)
  p3 <- ggplot(df4plot, aes(x = n_perturbed, y = lfc, col = donor, label = cell_line)) + 
    xlab("# of Cells") + ylab("On-target")+ 
    ggrepel::geom_text_repel() + 
    geom_point()+ theme_bw() + theme(legend.position = 'none')
  return(ggarrange(p1, p2, p3, ncol = 3))
}, SIMPLIFY = F)
pdf(paste0(plotsdir, '/', date, "_var_due_to_line.pdf"), width = 13.5, height = 5)
plotlist
dev.off()

## plot oct4 expression of each line (wild-type)

plot_jitter(df4lmm = control_expression %>% mutate(target = NonTargetGeneName), downstream_gene_name = 'POU5F1')

## ---- On-target vs. Line Effect

## Take a look at why these are a thing
tg <- 'GON4L'; dg <- 'H1F0'
plotlist <- mapply(1:dim(pairs2plot)[1], FUN = function(i){
  tg <- pairs2plot$target[i]; dg <- pairs2plot$downstream_gene_name[i]
  df4plot <- left_join(
    pica_lfcs_subset %>%
    dplyr::filter(target == tg & downstream_gene_name == dg) %>% 
    dplyr::select(c("cell_line", 'lfc')),
    pica_lfcs_subset %>%
      dplyr::filter(target == tg & downstream_gene_name == tg) %>% 
      dplyr::select(c("cell_line", 'on_target_expr' = 'lfc'))
  ) %>%
    mutate(donor = unlist(lapply(strsplit(cell_line, "_"), "[[", 1)))
  p1 <- ggplot(df4plot, aes(x=on_target_expr, y = lfc, 
                            col = donor, label = cell_line)) +
    ggrepel::geom_text_repel() + 
    xlab(tg) + ylab(dg) + ggtitle(paste0(dg, " Expression in ", tg," Knockdowns")) +
    geom_point() + theme_bw() + theme(legend.position = 'none')
 
  return(p1)
}, SIMPLIFY = F)
pdf(paste0(plotsdir, '/', date, "_trans_vs_cis.pdf"), width = 4.5, height = 5)
plotlist
dev.off()

## state of pluripotency?
## plot wild-type ipsc levels

## ---- Check if this a general pattern?


## ---- BSD
## BSD hits likely to CRISPR efficacy

top_variable_hits <- pica_lfcs_by_target %>% dplyr::filter(target != downstream_gene_name) %>%
  slice_max(v, n = 200)
tg <- 'GON4L'; dg <- 'H1F0'
plotlist <- mapply(1:dim(pairs2plot)[1], FUN = function(i){
  tg <- pairs2plot$target[i]; dg <- pairs2plot$downstream_gene_name[i]
  df4plot <- pica_lfcs_subset %>%
    dplyr::filter(target == tg & downstream_gene_name == dg)
  p1 <- ggplot(df4plot, aes(x=control_norm_expr + control_line_coef, y = lfc, 
                            col = donor, label = cell_line)) +
    ggrepel::geom_text_repel() + 
    xlab("Control Expression") + ylab("LFC") + ggtitle(paste0(dg, " Expression in ", tg," Knockdowns")) +
    geom_point() + theme_bw() + theme(legend.position = 'none')
  p2 <- ggplot(df4plot, aes(x=control_norm_expr + control_line_coef, y = control_norm_expr + control_line_coef+lfc, 
                            col = donor, label = cell_line)) +
    ggrepel::geom_text_repel() + 
    xlab("Control Expression") + ylab("Post-knockdown Expression") + ggtitle(paste0(dg, " Expression in ", tg," Knockdowns")) +
    geom_point() + theme_bw() + theme(legend.position = 'none')
  df4plot <- pica_lfcs_subset %>%
    dplyr::filter(target == tg & downstream_gene_name == tg)
  p3 <- ggplot(df4plot, aes(x = n_perturbed, y = lfc, col = donor, label = cell_line)) + 
    xlab("# of Cells") + ylab("On-target")+ 
    ggrepel::geom_text_repel() + 
    geom_point()+ theme_bw() + theme(legend.position = 'none')
  return(ggarrange(p1, p2, p3, ncol = 3))
}, SIMPLIFY = F)
pdf(paste0(plotsdir, '/', date, "_var_due_to_line.pdf"), width = 13.5, height = 5)
plotlist
dev.off()

## mean effect

## relative expression of bsd in control cells
kd_expr <- get_knockdown_expression(target_gene = 'POU5F1')
df4lmm <- get_df4lmm(kd_expr = control_expression %>% mutate(target = NonTargetGeneName), downstream_gene_name = "POU5F1", include_guide = F)
ggplot(df4lmm , aes(x = cell_line, y = post_kd_expr, col = unlist(lapply(strsplit(cell_line, "_"), "[[", 1)))) + 
  geom_jitter() + theme_bw() + theme(legend.position = 'none')




## ---- SNRPN

SNRPN <- get_knockdown_expression("SNRPN")
pairs2plot <- heritability_df %>%
  dplyr::filter(target == "SNRPN") %>%
  arrange(desc(bulk_mean_lfc))
## Take a look at why these are a thing
plotlist <- lapply(1:dim(pairs2plot)[1], FUN = function(i){
  dg <- pairs2plot$downstream_gene_name[i]
  plot_jitter(df4lmm = get_df4lmm(kd_expr = SNRPN, downstream_gene_name = dg, include_guide = F), mean_only = T) + ggtitle(dg)
})
pdf(paste0(plotsdir, '/', date, "_snrpn.pdf"), width = 5, height = 3)
plotlist
dev.off()
