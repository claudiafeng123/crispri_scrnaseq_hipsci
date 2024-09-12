
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))
suppressMessages(library(doParallel))
suppressMessages(library(ggpubr))
suppressMessages(library(tidyverse))

section_name <- "10f_trans_eQTL_overlap"
subsection_name <- "10f_02_trans_eQTL_overlap_compute_cis_effect_size"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
ExperimentName <- "Magpie"
date <- "2022-08-15" 
start_ind <-1700
no_genes2test <- 100

args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  ##needed arguments
  if (arg == "--date"){ magpie_version <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--start_ind"){ start_ind <- as.numeric(args[[ind + 1]]) }
  if (arg == "--no_genes2test"){ no_genes2test <- as.numeric(args[[ind + 1]]) }
  
}

# set relevent i/o paths
source(paste0(HomeFolder, "scripts/io/", ExperimentName, "_io.R"))
source(paste0(HomeFolder, "scripts/io/Magpie_Utils.R"))
outdir <- file.path(OutFolder, section_name, subsection_name)
plotsdir <- file.path(HTMLFolder, 'pipeline/', section_name, subsection_name)
if(!dir.exists(outdir)) dir.create(outdir, recursive = T)
if(!dir.exists(plotsdir)) dir.create(plotsdir, recursive = T)
print(section_name)
control_tag <- c("unassigned", NonTargetGeneName) # use all unassigned cells as control


## ---- LoadData

magpie_gts <- fread(paste0(OutFolder, 'Preprocess_External_Datasets/eQTLs/', date, '_gts_at_cis-eQTLs.tsv'))
target_meta <- fread(paste0(OutFolder, "7b_target_summary/", date, "_target_meta_data.csv"))

inlet_strengths <- c("Moderate", "Strong")
inlets_combined <- mapply(inlet_strengths, FUN = function(inlet_strength){
  rtn1 <- readRDS(file.path(OutFolder, "5a_combine", paste0(date, "_", inlet_strength, "_unassigned_cells_combined.RDS")))
  rtn2 <- readRDS(file.path(OutFolder, "5a_combine", paste0(date, "_", inlet_strength, "_nontarget_cells_combined.RDS")))
  rtn <-  merge(x = rtn1,
                y = rtn2,
                add.cell.ids = c("unassigned", NonTargetGeneName),
                project = paste(ExperimentName, "combined", sep = "_"),
                merge.data = TRUE)
  return(rtn)
}, SIMPLIFY = F)
magpie_gts <- magpie_gts %>% 
  dplyr::filter(gene_name %in% unlist(lapply(strsplit(row.names(inlets_combined$Moderate[["RNA"]]), ":"), "[[", 2)))

## number of cells infected with the guide
is_control_cell <- c(inlets_combined$Moderate$target_gene_name %in% control_tag, inlets_combined$Strong$target_gene_name %in% control_tag)
n_control <- length(which(is_control_cell))


cell_meta <- data.frame(cell_line = c(inlets_combined$Moderate$cell_line, inlets_combined$Strong$cell_line),
                        inlet = c(inlets_combined$Moderate$orig.ident, inlets_combined$Strong$orig.ident),
                        nCount_RNA = c(inlets_combined$Moderate$nCount_RNA, inlets_combined$Strong$nCount_RNA),
                        percent_MT = c(inlets_combined$Moderate$percent_MT, inlets_combined$Strong$percent_MT),
                        s.score = c(inlets_combined$Moderate$S.Score, inlets_combined$Strong$S.Score),
                        g2m.score = c(inlets_combined$Moderate$G2M.Score, inlets_combined$Strong$G2M.Score)
)
print(paste("Nontargeting/unassigned cells:", n_control))




## ---- RunLM

snps2test <- unique(magpie_gts$snp_id)[(start_ind + 1):(min(length(rownames(inlets_combined[[1]])), (start_ind+no_genes2test)))]
genes2test <- unique(magpie_gts %>%
                       dplyr::filter(snp_id %in% snps2test) %>%
                       .$gene_name)
print(paste("Testing", length(snps2test), "SNPs for perturbation effects across", length(genes2test), "genes"))
normalized_expression <- lapply(inlets_combined, FUN = function(x){
  x <- x@assays$RNA@data
  row.names(x) <- unlist(lapply(strsplit(row.names(x), ":"), "[[", 2))
  return(x[genes2test,])
}); names(normalized_expression) <- c("Moderate", "Strong")
#run lm

n_core <- parallel::detectCores()
print(paste("cores:", n_core))
registerDoParallel(min(c(n_core, 20, length(no_genes2test))))
control_lm <- foreach (snp_ind = 1:length(snps2test)) %dopar% {
  
  snp_name <- snps2test[snp_ind]
  gene_names <- gsub(magpie_gts %>%
                      dplyr::filter(snp_id == snp_name) %>%
                      .$gene_name, pattern = "_", replacement = "-")
  
  rtn <- lapply(gene_names, FUN = function(gene_name){
    y4lm <- c(normalized_expression$Moderate[gene_name,], normalized_expression$Strong[gene_name,])
    line_n_alt <- unlist(lapply(as.vector(t(magpie_gts[snp_ind,-c(1:9)])), FUN = function(x){length(which(unlist(strsplit(x, "[|]")) == 1))}))
    line_gts <- as.vector(t(magpie_gts[snp_ind,-c(1:9)]))
    line_gts <- ifelse(line_gts == '0|1', '1|0', line_gts)
    names(line_gts) <- names(line_n_alt) <- colnames(magpie_gts)[-c(1:9)]
    df4lm <- cell_meta %>%
      mutate(n_alt = line_n_alt[match(cell_line, names(line_gts))],
             gt = line_gts[match(cell_line, names(line_gts))]) %>%
      mutate(y4lm = y4lm)
    
    fit <- lm(y4lm ~ n_alt + inlet + nCount_RNA + percent_MT + s.score + g2m.score, data = df4lm)
    
    ## make a plot just for fun
    p <- ggplot(df4lm, aes(x = cell_line, y = y4lm, fill = gt)) + 
      xlab("Cell line") + ylab("Normalized expression") + ggtitle(paste0(snps2test[snp_ind], ' - ', gene_name)) + 
      geom_boxplot(outlier.shape = NA) + 
      annotate('text', label = paste0('cis_beta = ', formatC(magpie_gts$beta[snp_ind], 4)),
               x = -Inf, y = Inf, hjust = -0.25, vjust = 2) + 
      theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    pdf(paste0(plotsdir, '/', snps2test[snp_ind], '.pdf'), width = 6, height = 4)
    print(p)
    dev.off()
    
    rtn <- data.frame(
      snp_id = snps2test[snp_ind],
      gene_name = gene_name,
      control_beta = summary(fit)$coefficients['n_alt', 1],
      control_p_value = summary(fit)$coefficients['n_alt', 4],
      cis_eqtl_beta = magpie_gts$beta[snp_ind],
      cis_eqtl_p_value = magpie_gts$p_value[snp_ind],
      n_ref_hom = length(which(df4lm$n_alt == 0)),
      n_het = length(which(df4lm$n_alt == 1)),
      n_alt_hom = length(which(df4lm$n_alt == 2))
    )
    return(rtn)
  }) %>% bind_rows() %>% as.data.frame()
  
  #write 
  fnm = paste0(snps2test[snp_ind], ".RDS")
  saveRDS(rtn, file = file.path(outdir, fnm))
  #print(paste0(gene_name, "done!"))
}


sessionInfo()


