HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
section_name <- "11a_target_target_cor"
#date <- "2022-10-05"
ExperimentName <- "Magpie"
by_line <- F#; #l <- "tolg_4"

suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(doParallel))
suppressMessages(library(ggplot2))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
suppressMessages(library(ggExtra))
suppressMessages(library(ggpubr))

# options
args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
  if (arg == "--experiment_name"){ ExperimentName <- args[[ind + 1]] }
  if (arg == "--target_gene"){ target_gene <- args[[ind + 1]] }
  if (arg == "--by_line"){ by_line <- as.logical(args[[ind + 1]]) }
  if (arg == "--cell_line"){ l <- args[[ind + 1]] }
}


# i/o
source(file.path(HomeFolder, paste0("scripts/io/", ExperimentName, "_io.R")))
source(file.path(HomeFolder, "scripts/io/Magpie_Utils.R"))
if (by_line == T){
  outdir <- file.path(OutFolder, section_name, l, 'by_target/')
} else {
  outdir <- file.path(OutFolder, section_name, 'by_target/')
}
plotsdir <- file.path(HTMLFolder, "pipeline/", section_name, "cor_plots")
if(!dir.exists(outdir)) dir.create(outdir, recursive = T)
if(!dir.exists(plotsdir)) dir.create(plotsdir, recursive = T)
if (!exists("min_degs")){min_degs <- min_deg_per_target}
if (!exists("min_cells")){min_cells <- min_cells_per_gene}

print(paste0('results written to: ', outdir))
print(paste0('Experiment Name: ', ExperimentName))
print(paste0('cell_line: ', ifelse(by_line == T, l, 'mean across lines')))

## ---- LoadData

if (by_line == T){
  lfcs <- fread(file.path(OutFolder, "6f_calc_lfcs_transcriptome_wide_by_gene_per_line", paste0(date, "_all_lines_with_adj_pval.tsv.gz")))
  lfcs <- filter(lfcs, cell_line == l)
} else {
  lfcs <- fread(file.path(OutFolder, "6b_calc_lfcs_transcriptome_wide_by_gene", paste0(date, "_combined_with_adj_pval.tsv.gz")))
}


target_meta <- fread(file.path(OutFolder, "7b_target_summary", paste0(date, "_target_meta_data.csv")))

off_target <- fread(file.path(OutFolder, "Preprocess_External_Datasets", "off_target", paste0("off_target_pairs_", tolower(ExperimentName), "-", date, "-version.tsv.gz")),
                    header = T) %>% filter(target_gene != off_target_gene)
complex_pair_fnm <- sort(grep(list.files(file.path(OutFolder, "Preprocess_External_Datasets", "omnipath"), pattern = "protein_complex_pairs"),
                              pattern= date, value = T))[1]
complex_pairs <- fread(file.path(OutFolder, "Preprocess_External_Datasets", "omnipath", complex_pair_fnm ),
                       header = T)

## ---- CalculateCor

targets2consider <- unique(lfcs$target)
lfcs_split <- lfcs %>% filter(target %in% targets2consider)
lfcs_split <- split.data.frame(lfcs, f = lfcs$target)

if (target_gene %in% targets2consider){
  target_target_cor <- mapply(i = 1:length(targets2consider), FUN = function(i){
    print(i)
    df1 <- lfcs_split[[target_gene]]
    df2 <- lfcs_split[[targets2consider[i]]]
    c <- calc_cor(df1, df2, by_pval = T, sig_pval_thresh = sig_pval_thresh)
    c$df$is_protein_complex_pair <- (paste0(c$df$gene_1, "_", c$df$gene_2) %in% paste0(complex_pairs$gene_1, "_", complex_pairs$gene_2) | 
                                       paste0(c$df$gene_2, "_", c$df$gene_1) %in% paste0(complex_pairs$gene_1, "_", complex_pairs$gene_2))
    c$df$is_off_target <- (paste0(c$df$gene_1, "_", c$df$gene_2) %in% paste0(off_target$target_gene, "_", off_target$off_target_gene) | 
                             paste0(c$df$gene_2, "_", c$df$gene_1) %in% paste0(off_target$target_gene, "_", off_target$off_target_gene))
    
    return(c$df)
    
  }, SIMPLIFY = F)
  
  target_target_cor <- as.data.frame(bind_rows(target_target_cor))
  
  if (by_line == T){
    target_target_cor <- target_target_cor %>%
      mutate(cell_line = l)
    fwrite(target_target_cor, paste0(outdir, "/", date, "_", target_gene,"_", l, "_", "target_target_cor.tsv"), 
           sep = "\t", quote = F)
  } else {
    fwrite(target_target_cor, paste0(outdir, "/", date, "_", target_gene,"_", "target_target_cor.tsv"), 
           sep = "\t", quote = F)
  }
  
}

