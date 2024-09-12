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
suppressMessages(library(lattice))
suppressMessages(library(variancePartition))
suppressMessages(library(gtools))
suppressMessages(library(lmtest))
suppressMessages(library(org.Hs.eg.db))

# set relevent i/o paths
section_name <- "15b_downstream_heritability"
subsection_long_name <- "15b_01_downstream_heritability_calc_downstream_meta"

# set relevent i/o paths
date <- "2024-05-17"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
ExperimentName <- "Pica"

# options for testing
io_path <- paste0("scripts/io/", ExperimentName, "_io.R")
utils_path <- "scripts/io/Magpie_Utils.R"
control_tag <- c("unassigned", "NonTarget") # use all unassigned cells as control
keep_singleton_lines <- FALSE
line_set_ind <- 7
perm_seed <- 0
start_perm_ind <- 1
donors2include='eipl_iudw_jejf_kolf_paab'

args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  ##needed arguments
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  
  ## donors to include
  if (arg == "--donors2include"){ donors2include <- args[[ind + 1]] }
  if (arg == "--line_set_ind"){ line_set_ind <- as.numeric(args[[ind + 1]]) }
  
  ## not needed
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
  if (arg == "--utils_path"){ utils_path <- args[[ind + 1]] }
  if (arg == "--io_path"){ io_path <- args[[ind + 1]] }
  
  ## permuting
  if (arg == "--permute_cells"){ permute_cells <- as.logical(args[[ind + 1]]) }
  if (arg == "--permute_lines"){ permute_lines <- as.logical(args[[ind + 1]]) }
  if (arg == "--rnd_seed_start"){ rnd_seed_start <- as.numeric(args[[ind + 1]]) }
  if (arg == "--max_seeds2test"){ max_seeds2test <- as.numeric(args[[ind + 1]]) }
  if (arg == "--n_perms"){ n_perms <- as.numeric(args[[ind + 1]]) }
  if (arg == "--keep_singleton_lines"){ keep_singleton_lines <- as.logical(args[[ind + 1]]) }
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

fnm <- sort(list.files(file.path(OutFolder, "15a_get_line_perms"), pattern = "_lines_per_target.tsv", full.names = T), decreasing = T)[1]
donor_sets <- fread(fnm )
if (exists("donors2include")){
  d <- sort(unlist(strsplit(donors2include, '_')))
} else {
  d <- sort(unlist(strsplit(donor_sets$lines[which(donor_sets$line_set_name == line_set_ind)], ";")))
}

outdir <- file.path(file.path(OutFolder, section_name, subsection_long_name))
if(!dir.exists(outdir)) dir.create(outdir, recursive = T)
plotsdir <- file.path(HTMLFolder, "pipeline", section_name, subsection_long_name) 
if(!dir.exists(plotsdir)) dir.create(plotsdir, recursive = T)


## ---- Calculate Heritability

### some functions
most_recent <- sort(unique(gsub(list.files(paste0(OutFolder, "6b_calc_lfcs_transcriptome_wide_by_gene/", "by_downstream")),
                                pattern = "_.*", replacement = "")), decreasing = T)[1]
downstream_gene_fnms <- list.files(paste0(OutFolder, "6b_calc_lfcs_transcriptome_wide_by_gene/", "by_downstream"), pattern = most_recent)
#downstream_gene_fnms <- downstream_gene_fnms[(start_ind+1):(start_ind + no_genes2test)]

n_donors <- length(d)
fnm <- sort(list.files(file.path(OutFolder, '15a_get_line_perms'), 
                       pattern = paste0("_line_set-", line_set_ind, "_", n_donors, "_perms_rnd-seed-", perm_seed, ".RDS"), full.names = T),
                       decreasing = T)[1]
all_permutations <- readRDS(fnm)
all_permutations <- all_permutations[
  (start_perm_ind):dim(all_permutations)[1],]


registerDoParallel(50)
res <- foreach(downstream_gene_fnm = downstream_gene_fnms) %dopar% {
  
  downstream_gene_name <- unlist(lapply(strsplit(gsub(downstream_gene_fnm, pattern = paste0(date, "_|-Gene-Expression.*"), replacement = ''), "-"),
                                        FUN = function(x){paste(x[-1], collapse = '-')}))
  df4lmm_control <- fread(paste0(OutFolder, "6b_calc_lfcs_transcriptome_wide_by_gene/", "by_downstream/", downstream_gene_fnm)) %>%
    dplyr::filter(target_gene_name %in% c(NonTargetGeneName, "unassigned") & 
                    donor %in% d) %>%
    mutate(post_kd_expr = regressed_lfc + cell_line_coef,
           downstream_gene_name = downstream_gene_name)
  
  ## with fixed effects
  models_to_test <- gsub(var_expl_models, pattern = '[+] on_target_expr |[+] \\(1[|]guide\\)', replacement = '')
  models_to_test <- paste0('post_kd_expr', models_to_test)
  names(models_to_test) <- names(var_expl_models)
  res <- suppressMessages(suppressWarnings(
    calc_fit_meta(df4lmm = df4lmm_control, 
                  fixed_effects = c(""),
                  random_effects = c("cell_line", "donor"),
                  models_to_test = models_to_test, compute_coefs = F)
  ))
  
  
  ## add some metadata
  bulk_stats <- df4lmm_control %>% group_by(cell_line) %>% 
    summarize(mean_expr = mean(post_kd_expr))
  res$var_expl$bulk_mean_control_expr = mean(bulk_stats$mean_expr)
  res$var_expl$bulk_control_expr_var = var(bulk_stats$mean_expr)
  
  
  return(res$var_expl)
  
  
}



## ---- Write

## variance explained
#var_expl <- res$var_expl
var_expl <- res %>% bind_rows() %>% as.data.frame()
var_expl <- var_expl %>%
  mutate(
    lfc_dLL_line = LL_m1 - LL_m0,
    lfc_dLL_donor = LL_m2 - LL_m1
  ) 

## add in chromosome

downstream_gene_entrez_ids <- select(org.Hs.eg.db, 
                                     keys = var_expl$downstream_gene_name,
                                     columns = c("ENTREZID", "SYMBOL"),
                                     keytype = "SYMBOL") %>%
  dplyr::filter(!(is.na(ENTREZID)))
downstream_gene_entrez_ids$CHROM <- unlist(lapply(as.list(org.Hs.egCHR[downstream_gene_entrez_ids$ENTREZID]), FUN = function(x){paste(x, collapse = '_')}))
var_expl$chr <- downstream_gene_entrez_ids$CHROM[match(var_expl$downstream_gene_name, downstream_gene_entrez_ids$SYMBOL)]

fwrite(var_expl, paste0(outdir, '/', date, "_", paste0(d, collapse = '_'), ".tsv.gz"), compress = "gzip", sep = '\t')

## ---- SessionInfo

sessionInfo()

## AARS had error at 840
## EXOSC2 had error at 2690

#fit1 <- lmer(lfc ~ on_target_expr + (1|cell_line) + (1|donor), df4lmm)
#fit2 <- lmer(lfc ~ on_target_expr + (1|cell_line), df4lmm)
#lrtest(fit1, fit2)

