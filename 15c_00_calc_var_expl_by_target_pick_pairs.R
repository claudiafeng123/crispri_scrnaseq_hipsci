
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(tidyverse))



# set relevent i/o paths
section_name <- "15c_calc_var_expl_by_target"
subsection_name <- "15c_00_calc_var_expl_by_target_pick_pairs"

# set relevent i/o paths
date <- "2024-06-03"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
ExperimentName <- "Pica"

# options for testing
io_path <- paste0("scripts/io/", ExperimentName, "_io.R")
utils_path <- "scripts/io/Magpie_Utils.R"
#donors2include <- 'eipl_iudw_jejf_kolf_paab'


args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  ##needed arguments
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  
  ## not needed
  if (arg == "--donors2include"){ donors2include <- args[[ind + 1]] }
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
  if (arg == "--utils_path"){ utils_path <- args[[ind + 1]] }
  if (arg == "--io_path"){ io_path <- args[[ind + 1]] }
  if (arg == "--subsection_name"){ subsection_name <- args[[ind + 1]] }
  
}

print(paste0(section_name))
print(paste("date:", date))
print(paste("home folder:", HomeFolder))


setwd(HomeFolder)
source(io_path)
source(utils_path)

outdir <- file.path(file.path(OutFolder, section_name, subsection_name))
if(!dir.exists(outdir)) dir.create(outdir, recursive = T)
plotsdir <- file.path(HTMLFolder, "pipeline", section_name, subsection_name) 
if(!dir.exists(plotsdir)) dir.create(plotsdir, recursive = T)

## ---- LoadData

lfc_by_line_fnm <- sort(list.files(file.path(OutFolder, "6f_calc_lfcs_transcriptome_wide_by_gene_per_line"), pattern = "_all_lines_with_adj_pval.tsv.gz", full.names = T), decreasing = T)[1]
lfc_by_line <- fread(lfc_by_line_fnm)
if (exists('donors2include')){
  d <- unlist(strsplit(donors2include, "_"))
  lfc_by_line <- lfc_by_line %>%
    dplyr::filter(unlist(lapply(strsplit(cell_line, "_"), "[[", 1)) %in% d)
}
lfc_pair_meta <- lfc_by_line %>%
  mutate(pair = paste0(target, "_", downstream_gene_name)) %>%
  group_by(pair) %>%
  summarize(max_lfc = max(lfc),
            min_lfc = min(lfc),
            min_pval_adj = min(pval_adj),
            n_significant = length(which(pval_adj < sig_pval_thresh)),
            opp_sign = max(lfc)*min(lfc)) %>%
  ungroup()
lfc_pair_meta$target <- unlist(lapply(strsplit(lfc_pair_meta$pair, "_"), "[[", 1))
lfc_pair_meta$downstream_gene_name <- unlist(lapply(strsplit(lfc_pair_meta$pair, "_"), "[[", 2))

pairs_to_permute <- lfc_pair_meta %>% 
  filter(min_pval_adj < 0.1 & target != downstream_gene_name & target != NonTargetGeneName)
pairs_to_permute <- pairs_to_permute %>%
  dplyr::select(c('target', 'downstream_gene_name', 'max_lfc', 'min_lfc', 'n_significant', 'opp_sign'))

## ----

pdf(paste0(plotsdir, '/', date, '_max_vs_min_lfc.pdf'), width = 7, height = 4.5)
scale_min <- min(lfc_pair_meta$min_lfc)
scale_max <- max(lfc_pair_meta$max_lfc)
ggplot(lfc_pair_meta, aes(x = min_lfc, y = max_lfc, 
                          col = ifelse(n_significant >= 1, "Permuted", "Not Permuted"))) + 
  xlim(c(scale_min, scale_max)) + ylim(c(scale_min, scale_max)) + 
  scale_color_manual('', values = c('Permuted' = 'black', "Not Permuted" = 'gray')) + 
  xlab("Minimum LFC") + ylab("Maximal LFC") + 
  geom_point(size = 0.4, alpha = 0.4) + geom_abline(slope = 1, lty = 3, col = 'red') + theme_bw()
dev.off()

## ---- Write File

fwrite(lfc_pair_meta, paste0(outdir, '/', date, "_pairs_lfc_meta.tsv"), sep = '\t')
fwrite(pairs_to_permute, paste0(outdir, '/', date, "_pairs_to_permute.tsv"), sep = '\t')

## ---- SessionInfo

sessionInfo()

