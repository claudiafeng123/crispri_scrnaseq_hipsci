
suppressMessages(library(data.table, warn.conflicts = F))
suppressMessages(library(ggplot2, warn.conflicts = F))
suppressMessages(library(dplyr, warn.conflicts = F))

section_name <- "8d_power_requirements"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
ExperimentName <- "Pica"
date <- "2022-10-05"

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
source(paste0(HomeFolder, "scripts/io/", ExperimentName, "_io.R"))
source(paste0(HomeFolder, "scripts/io/Magpie_Utils.R"))
outdir <- file.path(OutFolder, section_name)
plotsdir <- file.path(HTMLFolder, "pipeline/", section_name)
if(!dir.exists(outdir)) dir.create(outdir)
if(!dir.exists(plotsdir)) dir.create(plotsdir)
#if(!dir.exists(paste0(plotsdir, "/cor_plots/"))) dir.create(paste0(plotsdir, "/cor_plots/"))

## ---- LoadData

GuideMetadata <- fread(GuideMetadataPath)
target_genes <- unique(GuideMetadata$gene)

#genes_done

#
target_df <- fread(paste0(OutFolder, "/", "7b_target_summary", "/", date, "_target_meta_data.csv"))

#genes2consider <- dplyr::filter(target_df, (n_downstream_excl_target > 10 & n_cells >= 1000) | (n_downstream_excl_target < 100 & n_downstream_excl_target > 25 & n_cells >= 85)) %>%
#  dplyr::select(c("gene", "n_downstream_excl_target", "n_cells"))

genes2consider <- dplyr::filter(target_df, (n_downstream_excl_target > 10 & n_cells >= 1000)) %>%
  dplyr::select(c("gene", "n_downstream_excl_target", "n_cells"))
  
## try to avoid overcount
#genes2consider <- genes2consider %>% filter(!(gene %in% c("CNOT1", "CTR9", "INTS8", "MAX", "NRF1", "PAF1", "RTF1", "U2AF1", "ZC3H18", "C1orf159", "CNOT3", "FAM32A", "LENG8", "NFRKB", "NUDT21", "PARD3B", "TAF8", "VHL")))
fwrite(genes2consider, paste0(outdir, "/", date, "_targets2consider.tsv"), sep = "\t", quote = F)
genes2consider %>% select(c("gene", "n_downstream_excl_target", "n_cells"))


## ---- ## power from lines

cell_metadata_fnms <- grep(list.files(paste0(OutFolder, "5a_combine"), pattern = "_perturbed_cells_metadata.tsv"), pattern = date, value = T)
cell_metadata <- lapply(cell_metadata_fnms, FUN = function(f){fread(file.path(OutFolder, "5a_combine", f))}) %>% bind_rows() %>% as.data.frame()
cells_per_gene_per_line <- as.data.frame(table(cell_metadata %>% 
                                                 
                                                          dplyr::select(c("cell_line", "target_gene_name"))))
genes_for_power_of_line <- cells_per_gene_per_line %>%
  filter(target_gene_name %in% genes2consider$gene) %>%
  filter(Freq > 100) %>%
  mutate(donor = unlist(lapply(strsplit(as.character(cell_line), split = "_"), "[[", 1))) %>%
  group_by(target_gene_name) %>%
  summarize(n_donor = length(unique(donor)),
            n_lines = length(unique(cell_line)),
            n_paired_lines = length(which(table(unlist(lapply(strsplit(as.character(unique(cell_line)), split = "_"), "[[", 1))) == 2)),
            lines = paste(sort(cell_line), collapse = '; '),
            paired_lines = paste(
              names(table(unlist(lapply(strsplit(as.character(unique(cell_line)), split = "_"), "[[", 1))))[which(table(unlist(lapply(strsplit(as.character(unique(cell_line)), split = "_"), "[[", 1))) == 2)], 
              collapse = '; ')
  ) 

## current power
barplot(table(genes_for_power_of_line$n_paired_lines),
        xlab = "# Paired Lines", ylab = "# Knockdowns")

genes_for_power_of_line <- left_join(genes_for_power_of_line, genes2consider,
                                     c('target_gene_name' = 'gene')) %>%
  filter(n_paired_lines == 5) 
fwrite(genes_for_power_of_line, paste0(outdir, '/', date, '_genes_for_line_power.tsv'), sep = '\t')


##
n_paired_lines <- 5
decimals <- 1:(2^n_paired_lines)
m <- sapply(decimals,function(x){ as.integer(intToBits(x))})
m <- t(m[1:(n_paired_lines),1:(2^n_paired_lines - 1)])
m <- m[which(apply(m, 1, sum) > 1),]
m <- apply(m, 1, FUN = function(x){which(x == 1)})


## 
iterations_to_run <- lapply(1:dim(genes_for_power_of_line)[1], FUN = function(i){
  g <- genes_for_power_of_line$target_gene_name[i]
  lines <- unlist(strsplit(genes_for_power_of_line$paired_lines[i], '; '))
  rtn <- data.frame(
    g = g,
    donors = unlist(lapply(m, FUN = function(x){paste(lines[x], collapse = '_')}))
)
  return(rtn)
}) %>% bind_rows() %>% as.data.frame()
fwrite(iterations_to_run, paste0(outdir, '/', date, '_line_requirements.tsv'), sep = '\t')







