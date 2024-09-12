library(tidyverse)
library(data.table)
library(parallel)

# options
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
section_name="9b_target_gene_downregulation_by_line"
#ExperimentName <- "Magpie"
#date <- "2022-08-15"
ExperimentName <- "Pica"
date <- "2022-10-05"
REDO <- FALSE
n_rnd_assignments <- 100


args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--experiment_name"){ ExperimentName <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
  if (arg == "--REDO"){ REDO <- as.logical(args[[ind + 1]]) }
}

source(paste0(HomeFolder, "scripts/io/", ExperimentName, "_io.R"))
source(paste0(HomeFolder, "scripts/io/Magpie_Utils.R"))
outdir <- file.path(OutFolder, section_name)
plotsdir <- file.path(HTMLFolder, "pipeline/", section_name, "/")
if(!dir.exists(outdir)) dir.create(outdir)
if(!dir.exists(plotsdir)) dir.create(plotsdir)

## ---- LoadData

if (!(exists("min_cells"))){min_cells <- min_cells_per_gene}
if (!(exists("min_degs"))){min_degs <- min_deg_per_target}

GuideMetadata <- fread(GuideMetadataPath)
experiment_summary <- readRDS(paste0(OutFolder, "/7a_experiment_summary/", date, "_7a_experiment_summary.RDS"))


## ---- Data Preprocessing 

## target meta
target_df_out <- file.path(OutFolder, section_name, paste0(date, "_target_meta_by_line.csv"))
if (REDO == T | !(file.exists(target_df_out))){
  lfc_res_fnm <- file.path(OutFolder, "6f_calc_lfcs_transcriptome_wide_by_gene_per_line", paste0(date, "_all_lines_with_adj_pval.tsv.gz"))
  lfc_res <- fread(lfc_res_fnm)
  #df4lmm <- fread(file.path(OutFolder, "6j_calc_lfc_per_gene_per_donor", paste0(date, "_df4lmm.tsv.gz")))
  #mean_downstream_expr_by_line <- filter(df4lmm, target == target[1]) %>% select(c("downstream_gene_name", "mean_downstream_expression"))
  n_deg <- lfc_res %>% 
    filter(pval_adj < sig_pval_thresh) %>% 
    select(c("cell_line", "target")) %>% 
    table() %>% as.data.frame() %>%
    mutate(to_match = paste0(cell_line, "_", target))
  n_deg_excl_target <- lfc_res %>% 
    filter(pval_adj < sig_pval_thresh & target  != downstream_gene_name) %>% 
    select(c("cell_line", "target")) %>% 
    table() %>% as.data.frame() %>%
    mutate(to_match = paste0(cell_line, "_", target))
  
  target_df <- filter(lfc_res, target == downstream_gene_name) %>% 
    select(c("target", "cell_line", "lfc", "pval_adj", "n_perturbed", "control_norm_expr"))
  all_lines <- unique(target_df$cell_line)
  ## add in downstream gene expression by line
  control_norm_expr_by_line <- mapply(tg = unique(target_df$target), FUN = function(tg){
    control_fnm <- grep(list.files(paste0(OutFolder, "/6a_run_control_lm/"), pattern = paste0("-", gsub(tg, pattern = "_", replacement = "-"), "-Gene-")), pattern = date, value = T)
    control_lm <- readRDS(file.path(OutFolder, "6a_run_control_lm", control_fnm))
    
    rtn <- data.frame(
      target = tg,
      cell_line = all_lines,
      line_control_expr = as.numeric(control_lm$coefficients['(Intercept)'])
    ) %>%
      mutate(line_control_expr = ifelse(cell_line %in% gsub(names(control_lm$coefficients), pattern = "cell_line", replacement = ""), 
                                        line_control_expr + control_lm$coefficients[paste0('cell_line', cell_line)],
                                line_control_expr))
    
  }, SIMPLIFY = F) %>% bind_rows() %>% as.data.frame()
  
  target_df <- left_join(target_df, control_norm_expr_by_line) %>%
    mutate(to_match = paste0(cell_line, "_", target))
               
  target_df <- target_df %>%
    mutate(n_deg = n_deg$Freq[match(to_match, n_deg$to_match)]) %>%
    mutate(n_deg_excl_target = n_deg_excl_target$Freq[match(to_match, n_deg_excl_target$to_match)]) %>%
    select(c("target", "cell_line", "lfc", "pval_adj", "n_perturbed", "n_deg", "n_deg_excl_target", "line_control_expr", "control_norm_expr"))
  
  fwrite(target_df, target_df_out)
  
} else {
  target_df <- fread(target_df_out)
}




## line metadata 

line_metadata_out <- paste0(outdir, "/", date, "_line_metadata.tsv")

if (REDO == T | !(file.exists(line_metadata_out))){
  line_metadata <- fread(LineMetadataPath)
  line_metadata <- split.data.frame(line_metadata, f = line_metadata$name)
  line_metadata <- mapply(line_metadata, FUN = function(df){
    if (dim(df)[1] > 1){
      dCas9 = mean(df$dCas9)
      mScarlet = mean(df$mScarlet)
      df <- df[1,]
      df$dCas9 = dCas9
      df$mScarlet = mScarlet
    }
    return(df)
  }, SIMPLIFY = F) %>% bind_rows() %>% as.data.frame()
  
  
  assigned_cells_per_donor <- experiment_summary$assigned_cells_per_donor
  frac_cells_assigned <- data.frame(
    name = assigned_cells_per_donor$cell_line,
    frac_cells_assigned = assigned_cells_per_donor$assigned/(apply(assigned_cells_per_donor[, -1], 1, sum))
  )
  
  line_metadata <- left_join(line_metadata, frac_cells_assigned) %>%
    mutate(low_target_downregulation = (name %in% c("fiaj_3", "tolg_4") )) %>%
    mutate(low_fraction_of_assigned_cells = (name %in% c("fiaj_3", "tolg_4", "pipw_5", "oikd_2"))) %>%
    mutate(dCas9_level = ifelse(is.na(dCas9), NA, ifelse(name %in% c("pipw_5", "oikd_2"), "low", ifelse(name %in% c("eipl_3", "kolf_3", "tolg_6"), "moderate", "high")) )) %>%
    mutate(mScarlet_level = ifelse(is.na(dCas9), NA, ifelse( name %in% c("fiaj_3", "kolf_3"), "low", ifelse(name %in% c("eipl_1", "kolf_2", "oikd_2", "pipw_5", "zapk_3"), "moderate", "high")  ))) %>%
    mutate(pluripotency_level = ifelse(name %in% c("pipw_4", "pipw_5", "eipl_1", "eipl_3", "fiaj_3"), "low", ifelse(name %in% c("fiaj_1", "iudw_1", "iudw_4", "oikd_2", "oikd_5"), "moderate", "high")) ) %>%
    mutate(novelty_level = ifelse(name %in% c("kolf_2", "kolf_3"), "low", "high") )
  
  ## add in line with small effeect size (use number of DEGs)
  names(line_metadata)[1] <- "cell_line"
  
  fwrite(line_metadata, line_metadata_out, sep = "\t", quote = F)
  
} else {
  line_metadata <- fread(line_metadata_out)
}

## crispr_stats

crispr_stats_out <- paste0(outdir, "/", date, "_crispr_sequence_stats.tsv")

if (REDO == T | !(file.exists(crispr_stats_out))){
  crispr_sequence_stats <- mapply(l = sort(unique(line_metadata$cell_line)), FUN = function(l){
    df <- cell_metadata %>% filter(cell_line == l)
    rtn <- data.frame(
      line = l,
      sequence = crispr_sequences,
      mean_sequence_counts = apply(select(df, all_of(paste0(crispr_sequences, "_counts"))), 2, mean),
      mean_sequence_counts_assigned = apply(select(df %>% filter(assignment_status != "unassigned"), all_of(paste0(crispr_sequences, "_counts"))), 2, mean),
      mean_sequence_counts_unassigned = apply(select(df %>% filter(assignment_status == "unassigned"), all_of(paste0(crispr_sequences, "_counts"))), 2, mean),
      frac_nonzeros = apply(select(df, all_of(paste0(crispr_sequences, "_counts"))), 2, FUN = function(x){length(which(x > 0))/length(x)}),
      frac_nonzeros_assigned = apply(select(df %>% filter(assignment_status != "unassigned"), all_of(paste0(crispr_sequences, "_counts"))), 2, FUN = function(x){length(which(x > 0))/length(x)}),
      frac_nonzeros_unassigned = apply(select(df %>% filter(assignment_status == "unassigned"), all_of(paste0(crispr_sequences, "_counts"))), 2, FUN = function(x){length(which(x > 0))/length(x)})
    )
    return(rtn)
  }, SIMPLIFY = F)
  crispr_sequence_stats <- as.data.frame(bind_rows(crispr_sequence_stats))
  
  fwrite(crispr_sequence_stats, crispr_stats_out,
         quote = F, sep = "\t")
  #crispr_sequence_stats <- fread(paste0(outdir, "/", date, "_crispr_sequence_stats.tsv"))
  
} else {
  crispr_sequence_stats <- fread(crispr_stats_out)
}



## ---- MakePlots

fnm_html <- paste0(section_name,".html")
rmarkdown::render(file.path(CodeFolder, "Magpie", "pipeline", paste0(section_name, ".Rmd")),
                  output_file = file.path(plotsdir, fnm_html))


## ---- SessionInfo

sessionInfo()

## ---- Scratch

