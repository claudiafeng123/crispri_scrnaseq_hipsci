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
subsection_name <- "15c_01_calc_var_expl_by_target_run_lme"

# set relevent i/o paths
date <- "2024-06-03"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
ExperimentName <- "Pica"

# options for testing
io_path <- paste0("scripts/io/", ExperimentName, "_io.R")
utils_path <- "scripts/io/Magpie_Utils.R"
control_tag <- c("unassigned", "NonTarget") # use all unassigned cells as control
target_gene <- 'AAMP'
keep_singleton_lines <- FALSE
include_line_quality <- T
include_sex <- F
donors2include <- 'eipl_iudw_jejf_kolf_paab'

args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  ##needed arguments
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--target_gene"){ target_gene <- args[[ind + 1]] }
  if (arg == "--donors2include"){ donors2include <- args[[ind + 1]] }
  
  ## not needed
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
  if (arg == "--subsection_name"){ subsection_name <- args[[ind + 1]] }
  if (arg == "--utils_path"){ utils_path <- args[[ind + 1]] }
  if (arg == "--io_path"){ io_path <- args[[ind + 1]] }
  
  ## permuting
  if (arg == "--keep_singleton_lines"){ keep_singleton_lines <- as.logical(args[[ind + 1]]) }
  if (arg == "--include_line_quality"){ include_line_quality <- as.logical(args[[ind + 1]]) }
  if (arg == "--include_sex"){ include_sex <- as.logical(args[[ind + 1]]) }
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
sexes <- LineMetadata$sex; names(sexes) <- LineMetadata$name

outdir <- file.path(file.path(OutFolder, section_name, subsection_name))
if(!dir.exists(outdir)) dir.create(outdir, recursive = T)
plotsdir <- file.path(HTMLFolder, "pipeline", section_name) 
if(!dir.exists(plotsdir)) dir.create(plotsdir, recursive = T)

fnm_out_suffix <- paste0('permute-cells-FALSE_permute-lines-FALSE_include-line-quality-', include_line_quality, '_include-sex-', include_sex)


## ---- load data


if (target_gene == NonTargetGeneName){
  gene_expression <- get_control_expression() 
  gene_expression$target <- target_gene
} else {
  gene_expression <- get_knockdown_expression(target_gene = target_gene)
}
if (exists("donors2include")){
  d <- unlist(strsplit(donors2include, '_'))
  gene_expression <- filter(gene_expression, donor %in% d)
}
n_cells <- dim(gene_expression)[1]
n_lines <- length(unique(gene_expression$cell_line))
n_donors <- length(unique(gene_expression$donor))
print("cells per line:")
table(gene_expression$cell_line)

## is the target expressed?
target_ind <- grep(colnames(gene_expression), pattern = paste0(":", target_gene, ":"))
is_target_expressed <- ((length(target_ind) == 1) & (include_line_quality == T))
if (is_target_expressed){
  target_colname <- colnames(gene_expression)[target_ind]
  most_recent <- sort(list.files(paste0(OutFolder, "/6a_run_control_lm/"), pattern = '-'), decreasing = T)[1]
  target_control_lm_fnm <- list.files(paste0(OutFolder, "/6a_run_control_lm/", most_recent), 
                                      full.names = T,
                                      pattern = paste0("gene-", gsub(target_colname, pattern = ":", replacement = "-")))
  target_control_lm <- readRDS(target_control_lm_fnm)
  target_control_lm_coefficients <- target_control_lm$coefficients[grep(names(target_control_lm$coefficients), pattern = "cell_line")]
  names(target_control_lm_coefficients) <- gsub(names(target_control_lm_coefficients), pattern= "cell_line", replacement = "")
}

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


## ---- Throw out singleton lines

if (keep_singleton_lines == F){
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
} else {
  lines2keep <- unique(gene_expression$cell_line)
}
n_cells_used <- dim(gene_expression)[1]
n_lines_used <- length(unique(gene_expression$cell_line))


## ---- Calculate Heritability

### some functions
most_recent <- sort(unique(gsub(grep(list.files(paste0(OutFolder, "6a_run_control_lm")), pattern = 'bsub_outs', invert = T, value = T), pattern = "_.*", replacement = "")), decreasing = T)[1]
downstream_gene_fnms <- list.files(paste0(OutFolder, "6a_run_control_lm/", most_recent))

## for testing a small number of downstream genes
#downstream_gene_fnms <- downstream_gene_fnms[840:850]
#
#downstream_gene_fnm <- grep(downstream_gene_fnms, pattern = "ZNF257", value = T)
#downstream_gene_fnm<- downstream_gene_fnms[1]
#downstream_gene_fnm<- downstream_gene_fnms[6406]

## decide which models need to be fun
models_to_test <- var_expl_models
if (is_target_expressed == F){
  models_to_test <- gsub(models_to_test, pattern = ' [+] on_target_expr', replacement = '')
  fe <- fe[which(fe != "on_target_expr")]
}

print(models_to_test)
#downstream_gene_fnm = downstream_gene_fnms[1]
downstream_gene_fnm = grep(downstream_gene_fnms, pattern = paste0("-", "FAM135A", "-"), value = T)

registerDoParallel(50)
res <- foreach(downstream_gene_fnm = downstream_gene_fnms) %dopar% {
  
  res <- list()
  downstream_gene_ind <- which(gsub(colnames(gene_expression), pattern = ":", replacement = '-') == gsub(downstream_gene_fnm, pattern = paste0("gene-|[.]RDS"), replacement = ""))
  downstream_gene_name <- unlist(lapply(strsplit(colnames(gene_expression)[downstream_gene_ind], ":"), "[[", 2))
  
  df4lmm_perturbed <- get_df4lmm(kd_expr = gene_expression,
                                 downstream_gene_name = downstream_gene_name,
                                 lines2include = lines2keep$line,
                                 include_line_quality = (is_target_expressed),
                                 include_sex = include_sex)
  
  
  
  ms <- paste0('lfc', models_to_test); names(ms) <- names(models_to_test)
  res_lfc <- suppressMessages(suppressWarnings(
    calc_fit_meta(df4lmm = df4lmm_perturbed, 
                  models_to_test = ms, compute_coefs = length(fe)>= 1, 
                  include_variance_explained = T)
  ))
  
  
  names(res_lfc$var_expl)[-(1:2)] <- paste0('lfc_', names(res_lfc$var_expl)[-(1:2)])
  ms <- paste0('post_kd_expr', models_to_test); names(ms) <- names(models_to_test)
  res_post_kd_expr <- suppressMessages(suppressWarnings(
    calc_fit_meta(df4lmm = df4lmm_perturbed, 
                  models_to_test = ms, compute_coefs = length(fe)>= 1, 
                  fixed_effects = fe)
  ))
  names(res_post_kd_expr$var_expl)[-(1:2)] <- paste0('post_kd_expr_', names(res_post_kd_expr$var_expl)[-(1:2)])
  
  res$var_expl <- data.frame(
    target = target_gene,
    downstream_gene_name = downstream_gene_name,
    res_lfc$var_expl[-(1:2)],
    res_post_kd_expr$var_expl[-(1:2)]
  )
  ## add some metadata
  ## also add some metadata for lfc
  ## bulk variance
  ## mean lfc (bulk)
  ## maximum lfc (abs)
  bulk_stats <- df4lmm_perturbed %>% group_by(cell_line) %>% 
    summarize(mean_lfc = mean(lfc))
  res$var_expl$bulk_mean_lfc = mean(bulk_stats$mean_lfc)
  res$var_expl$bulk_lfc_var = var(bulk_stats$mean_lfc)
  res$var_expl$max_abs_lfc = max(abs(bulk_stats$mean_lfc))
  res$var_expl$bulk_range_lfc = diff(range(bulk_stats$mean_lfc))
  
  
  return(res)
  
  
}




## ---- Write

## variance explained
#var_expl <- res$var_expl
var_expl <- lapply(res, FUN = function(l){l$var_expl}) %>% bind_rows() %>% as.data.frame()
## decide which models need to be fun
var_expl$lfc_dLL_line <- var_expl$lfc_LL_m1 - var_expl$lfc_LL_m0
var_expl$post_kd_expr_dLL_line <- var_expl$post_kd_expr_LL_m1 - var_expl$post_kd_expr_LL_m0
var_expl$lfc_dLL_donor <- var_expl$lfc_LL_m2 - var_expl$lfc_LL_m1
var_expl$post_kd_expr_dLL_donor <- var_expl$post_kd_expr_LL_m2 - var_expl$post_kd_expr_LL_m1
fwrite(var_expl, paste0(outdir, '/', date, "_", fnm_out_suffix, "_", target_gene, "_var_expl.tsv.gz"), compress = "gzip", sep = '\t')

## ---- SessionInfo

sessionInfo()

## AARS had error at 840
## EXOSC2 had error at 2690

#fit1 <- lmer(lfc ~ on_target_expr + (1|cell_line) + (1|donor), df4lmm)
#fit2 <- lmer(lfc ~ on_target_expr + (1|cell_line), df4lmm)
#lrtest(fit1, fit2)

