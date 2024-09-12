suppressMessages(library(Matrix))
suppressMessages(library(Seurat))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(tidyverse))
suppressMessages(library(parallel))
suppressMessages(library(doParallel))
suppressMessages(library(lme4))
#suppressMessages(library(lme4qtl))
suppressMessages(library(lattice))
suppressMessages(library(variancePartition))
suppressMessages(library(gtools))
suppressMessages(library(lmtest))



# set relevent i/o paths
section_name <- "15c_calc_var_expl_by_target"

# set relevent i/o paths
date <- "2024-06-03"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
ExperimentName <- "Pica"

# options for testing
io_path <- paste0("scripts/io/", ExperimentName, "_io.R")
utils_path <- "scripts/io/Magpie_Utils.R"
control_tag <- c("unassigned", "NonTarget") # use all unassigned cells as control
target_gene <- "POU5F1"
keep_singleton_lines <- FALSE
include_line_quality <- T
include_sex <- T
num_perms_to_beat <- 10
max_num_perms_to_compute <- 10^4
start_perm_ind <- 1
perm_seed <- 0
subsection_name <- "15c_03_calc_var_expl_by_target_permute_lines_all_pairs"
val2calc <- "post_kd_expr" #'lfc' or 'post_kd_expr'
donors2include <- 'eipl_iudw_jejf_kolf_paab'

args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  ##needed arguments
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--target_gene"){ target_gene <- args[[ind + 1]] }
  if (arg == "--downstream_gene_name"){ downstream_gene_name <- args[[ind + 1]] }
  if (arg == "--val2calc"){ val2calc <- args[[ind + 1]] }
  if (arg == "--donors2include"){ donors2include <- args[[ind + 1]] }
  
  ## not needed
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
  if (arg == "--utils_path"){ utils_path <- args[[ind + 1]] }
  if (arg == "--io_path"){ io_path <- args[[ind + 1]] }
  
  ## linear model parameters
  if (arg == "--keep_singleton_lines"){ keep_singleton_lines <- as.logical(args[[ind + 1]]) }
  if (arg == "--include_line_quality"){ include_line_quality <- as.logical(args[[ind + 1]]) }
  if (arg == "--include_sex"){ include_sex <- as.logical(args[[ind + 1]]) }
  
  ## parameters
  if (arg == "--perm_seed"){ perm_seed <- as.numeric(args[[ind + 1]]) }
  if (arg == "--num_perms_to_beat"){ num_perms_to_beat <- as.numeric(args[[ind + 1]]) }
  if (arg == "--start_perm_ind"){ start_perm_ind <- as.numeric(args[[ind + 1]]) }
  if (arg == "--max_num_perms_to_compute"){ max_num_perms_to_compute <- as.numeric(args[[ind + 1]]) }
  
}

print(paste0(section_name))
print(paste("date:", date))
print(paste("home folder:", HomeFolder))
if(exists("donors2include")){print('including donors: '); print(donors2include)}
print(paste("val2calc:", val2calc))



setwd(HomeFolder)
source(io_path)
source(utils_path)

GuideMetadata <- fread(GuideMetadataPath)
LineMetadata <- fread(LineMetadataPath)
sexes <- LineMetadata$sex; names(sexes) <- LineMetadata$name

outdir <- file.path(file.path(OutFolder, section_name, subsection_name, val2calc, target_gene))
if(!dir.exists(outdir)) dir.create(outdir, recursive = T)
plotsdir <- file.path(HTMLFolder, "pipeline", section_name) 
if(!dir.exists(plotsdir)) dir.create(plotsdir, recursive = T)

fnm_out_suffix <- paste0('permute-cells-FALSE_permute-lines-TRUE_include-line-quality-', include_line_quality, '_include-sex-', include_sex)


## ---- load data

gene_expression <- get_knockdown_expression(target_gene = target_gene)
if (exists("donors2include")){
  d <- unlist(strsplit(donors2include, '_'))
  gene_expression <- filter(gene_expression, donor %in% d)
}
n_cells <- dim(gene_expression)[1]
n_lines <- length(unique(gene_expression$cell_line))
n_donors <- length(unique(gene_expression$donor))
print("cells per line:")
table(gene_expression$cell_line)


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
  dplyr::filter(target_gene_name == target_gene) %>%
  dplyr::select(c("V1", "guide" = 'feature_call'))
## add to gene expression df
gene_expression <- left_join(gene_expression,cell_metadata, by = c("id" = "V1"))
print('gene expression loaded!')



## is the target expressed?
target_ind <- grep(colnames(gene_expression), pattern = paste0(":", target_gene, ":"))
is_target_expressed <- (length(target_ind) == 1)
if (is_target_expressed){
  target_colname <- colnames(gene_expression)[target_ind]
  most_recent <- sort(grep(list.files(paste0(OutFolder, "/6a_run_control_lm/")), pattern = 'bsub_outs', invert = T, value = T), decreasing = T)[1]
  target_control_lm <- readRDS(paste0(OutFolder, "/6a_run_control_lm/", most_recent, '/', paste0("gene-", gsub(target_colname, pattern = ":", replacement = "-"), ".RDS")))
  target_control_lm_coefficients <- target_control_lm$coefficients[grep(names(target_control_lm$coefficients), pattern = "cell_line")]
  names(target_control_lm_coefficients) <- gsub(names(target_control_lm_coefficients), pattern= "cell_line", replacement = "")
}


## true values
dLL <- fread(file.path(OutFolder, "15c_calc_var_expl_by_target", "15c_02_calc_var_expl_by_target_variance_decomposition", paste0(date, "_var_expl_combined_no_pval.tsv.gz"))) %>%
  dplyr::select(c("target", "downstream_gene_name", 'lfc_dLL_donor', 'lfc_dLL_line')) %>%
  dplyr::filter(target == target_gene)

## ---- Throw out singleton lines

print("removing singleton lines...")
lines2keep <- data.frame(
  line = unique(gene_expression$cell_line),
  donor = unlist(lapply(strsplit(unique(gene_expression$cell_line), "_"), "[[", 1)),
  n_cells = as.numeric(table(gene_expression$cell_line))
) %>% filter(n_cells > min_cells_per_gene_per_donor)
donors2keep <- table(lines2keep$donor)
donors2keep  <- names(donors2keep)[which(donors2keep  == 2)]
lines2keep <- lines2keep %>% filter(donor %in% donors2keep)
print(paste0("keeping ", paste(lines2keep$line, collapse = ', ')))
gene_expression <- filter(gene_expression, cell_line %in% lines2keep$line)
print("cells per line:")
print(table(gene_expression$cell_line))
n_cells_used <- dim(gene_expression)[1]
n_lines_used <- length(unique(gene_expression$cell_line))


## ---- Permutations
## permute lines only here 
## cell permutation comes later

## permute lines
fnm <- sort(list.files(file.path(OutFolder, "15a_get_line_perms"), pattern = '_target_line_mapping', full.names = T), decreasing = T)[1]
target_line_mapping <- fread(fnm)
if (exists('donors2include')){
  lines_per_target <- fread(fnm)
  ds <- unlist(strsplit(donors2include, '_'))
  fnm <- sort(list.files(file.path(OutFolder, "15a_get_line_perms"), pattern = '_lines_per_target', full.names = T), decreasing = T)[1]
  line_ind_table <- fread(fnm)
  target_line_mapping$paired_lines <- unlist(lapply(1:length(target_line_mapping$paired_lines), FUN = function(i){
    well_powered_lines <- unlist(strsplit(target_line_mapping$paired_lines[i], ";"))
    well_powered_lines <- well_powered_lines[which(well_powered_lines %in% ds)]
    well_powered_lines <- paste(well_powered_lines, collapse = ';')
    return(well_powered_lines)
  }))
  target_line_mapping$n_paired_lines <- unlist(lapply(strsplit(target_line_mapping$paired_lines, ";"), length))
  target_line_mapping$line_set_name <- line_ind_table$line_set_name[match(target_line_mapping$paired_lines, line_ind_table$lines)]
}
target_line_mapping <- target_line_mapping %>% dplyr::filter(gene == target_gene)

fnm <- sort(list.files(file.path(OutFolder, "15a_get_line_perms"), pattern = paste0("_line_set-", target_line_mapping$line_set_name, "_", target_line_mapping$n_paired_lines, "_perms_rnd-seed-", perm_seed), full.names = T), decreasing = T)[1]
all_permutations <- readRDS(fnm)
all_permutations <- all_permutations[
  (start_perm_ind):dim(all_permutations)[1],]


## ---- Calculate Heritability

## decide which models need to be fun
models_to_test <- var_expl_models
fixed_effects = c("on_target_expr", "sex")
if (include_sex == F){
  models_to_test <- gsub(models_to_test, pattern = ' [+] sex', replacement = '')
  fixed_effects <- fixed_effects[which(fixed_effects != "sex")]
}
if (is_target_expressed == F | include_line_quality == F){
  models_to_test <- gsub(models_to_test, pattern = ' [+] on_target_expr', replacement = '')
  fixed_effects <- fixed_effects[which(fixed_effects != "on_target_expr")]
}
print(models_to_test)

fnm <- sort(list.files(file.path(OutFolder,  section_name,'15c_00_calc_var_expl_by_target_pick_pairs'), pattern = "_pairs_to_permute", full.names = T), decreasing = T)[1]
all_pairs2test <- fread(fnm)
downstream_genes2test <- all_pairs2test %>%
  dplyr::filter(target == target_gene) %>% .$downstream_gene_name

#already_done <- unlist(lapply(strsplit(gsub(list.files(outdir), pattern = paste0(date, '_', target_gene, '_'), replacement = ''), split= "_"), "[[", 1))
#downstream_genes2test <- downstream_genes2test[which(!(downstream_genes2test %in% already_done))]
print(paste0('testing ', length(downstream_genes2test), ' downstream genes'))

if (length(downstream_genes2test) > 1){
  registerDoParallel(min(50, parallel::detectCores(), length(downstream_genes2test)))
  compute_permuted_LLs <- foreach(dg = downstream_genes2test) %dopar% {
    df2write <- run_perm(df4lmm = get_df4lmm(kd_expr = gene_expression, 
                                             downstream_gene_name = dg), 
                         val2calc = val2calc,
                         models_to_test = models_to_test,
                         all_permutations2test = all_permutations)
    fwrite(df2write, paste0(outdir, '/', date, "_", target_gene, "_", dg, "_", start_perm_ind, "-", dim(df2write)[1] - 1, ".tsv.gz" ), sep = '\t')
  }
  
} else if (length(downstream_genes2test) == 1) {
  dg <- downstream_genes2test
  df2write <- run_perm(df4lmm = get_df4lmm(kd_expr = gene_expression, 
                                           downstream_gene_name = dg), 
                       h = dLL %>% filter(downstream_gene_name == dg) %>% .$lfc_dLL_donor,
                       all_permutations2test = all_permutations)
  fwrite(df2write, paste0(outdir, '/', date, "_", target_gene, "_", dg, "_", start_perm_ind, "-", dim(df2write)[1] - 1, ".tsv.gz" ), sep = '\t')
}


## ---- SessionInfo

sessionInfo()

