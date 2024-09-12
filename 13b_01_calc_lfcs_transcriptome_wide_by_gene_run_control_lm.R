suppressMessages(library(Matrix))
suppressMessages(library(Seurat))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyverse))
suppressMessages(library(parallel))
suppressMessages(library(doParallel))

# set relevent i/o paths
section_name <- "13b_calc_lfcs_transcriptome_wide_by_gene"

# set relevent i/o paths
date <- "2022-08-15"
target_gene <- "POU5F1"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"

# options for testing
control_tag <- c("unassigned", "NonTarget") # use all unassigned cells as control
no_genes2test <- 10
start_ind <- 0
PLOT <- FALSE

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
  if (arg == "--scale_factor"){ scale_factor <- as.numeric(args[[ind + 1]]) }
  if (arg == "--start_ind"){ start_ind <- as.numeric(args[[ind + 1]]) }
  if (arg == "--no_genes2test"){ no_genes2test <- as.numeric(args[[ind + 1]]) }
  if (arg == "--lfc_base"){ lfc_base <- as.numeric(args[[ind + 1]]) }
  if (arg == "--PLOT"){ PLOT <- as.logical(args[[ind + 1]]) }
  if (arg == "--cell_type"){ cell_type <- args[[ind + 1]] }
  if (arg == "--screen"){ screen <- args[[ind + 1]] }
  if (arg == "--section_name"){section_name <- args[[ind + 1]] }
}

print(paste("date:", date))
print(paste("home folder:", HomeFolder))


setwd(HomeFolder)
io_path <- paste0("scripts/io/", ExperimentName, "_io.R")
utils_path <- "scripts/io/Magpie_Utils.R"
source(io_path)
source(utils_path)



if (ExperimentName == "Replogle"){
  outdir <- file.path(OutFolder, section_name, paste0(cell_type, "_", screen))
  plotsdir <- file.path(HTMLFolder, section_name)
} else {
  outdir <- file.path(OutFolder, section_name)
  plotsdir <- file.path(HTMLFolder, "pipeline", section_name)
}
if(!dir.exists(outdir)) dir.create(outdir, recursive = T)
if(!dir.exists(plotsdir)) dir.create(plotsdir)
print(outdir)

if (ExperimentName == "Magpie"){
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
  ## number of cells infected with the guide
  is_control_cell <- c(inlets_combined$Moderate$target_gene_name %in% control_tag, inlets_combined$Strong$target_gene_name %in% control_tag)
  n_control <- length(which(is_control_cell))
  genes2test <- rownames(inlets_combined[[1]])[(start_ind + 1):(min(length(rownames(inlets_combined[[1]])), (start_ind+no_genes2test)))]
  print(paste("Nontargeting/unassigned cells:", n_control))
  normalized_expression <- lapply(inlets_combined, FUN = function(x){
    x@assays$RNA@data[genes2test,]
  }); names(normalized_expression) <- c("Moderate", "Strong")
  #run lm
  data4lm <- data.frame(cell_line = c(inlets_combined$Moderate$cell_line, inlets_combined$Strong$cell_line),
                        inlet = c(inlets_combined$Moderate$orig.ident, inlets_combined$Strong$orig.ident),
                        nCount_RNA = c(inlets_combined$Moderate$nCount_RNA, inlets_combined$Strong$nCount_RNA),
                        percent_MT = c(inlets_combined$Moderate$percent_MT, inlets_combined$Strong$percent_MT),
                        s.score = c(inlets_combined$Moderate$S.Score, inlets_combined$Strong$S.Score),
                        g2m.score = c(inlets_combined$Moderate$G2M.Score, inlets_combined$Strong$G2M.Score)
  )
  
} else if (ExperimentName == "Pica"){
  
  #fnms <- grep(list.files(file.path(OutFolder, "5a_combine"), pattern = date), pattern = "control", value = T)
  inlets_combined <- readRDS(file.path(OutFolder, "5a_combine", paste0(date, "_ALL_control_cells.RDS")))
  print("seurat object loaded")
  genes2test <- rownames(inlets_combined)[(start_ind + 1):(min(length(rownames(inlets_combined)), (start_ind+no_genes2test)))]
  normalized_expression <- inlets_combined@assays$RNA@data[genes2test,]
  #run lm
  print("cell metadata obtained!")
  data4lm <- data.frame(cell_line = inlets_combined$cell_line,
                        inlet = c(inlets_combined$orig.ident),
                        nCount_RNA = c(inlets_combined$nCount_RNA),
                        percent_MT = c(inlets_combined$percent_MT),
                        s.score = c(inlets_combined$S.Score),
                        g2m.score = c(inlets_combined$G2M.Score)
  )
  
} 

#only test a subset of genes

print(paste("Testing", length(genes2test), "genes for perturbation effects"))

# renormalize
#log(t(t(counts_matrix)*10000/inlets_combined$nCount_RNA) + ps)
##not fixed since i don't think we're using
if (scale_factor != 1) {
  inlets_combined <- NormalizeData(inlets_combined, scale.factor = scale_factor*10^5)
  print("seurat object renormalized!")
}

n_core <- parallel::detectCores()
print(paste("cores:", n_core))
registerDoParallel(min(c(n_core, 20, length(no_genes2test))))
control_lm <- foreach (i = 1:length(genes2test)) %dopar% {
  #print(i)
  gene_name <- gsub(genes2test[i], pattern = ":", replacement = "-")
  #print(gene_name)
  if (ExperimentName == "Magpie"){
    y4lm <- c(normalized_expression$Moderate[i,], normalized_expression$Strong[i,])
    fit <- lm(y4lm ~ cell_line + inlet + nCount_RNA + percent_MT + s.score + g2m.score, data = data4lm)
  } else if (ExperimentName == "Pica"){
    y4lm <- normalized_expression[i,]
    fit <- lm(y4lm ~ cell_line + inlet + nCount_RNA + percent_MT + s.score + g2m.score, data = data4lm)
  } else if (ExperimentName == "Replogle"){
    y4lm <- normalized_expression[i,]
    fit <- lm(y4lm ~ nCount_RNA + percent_MT + s.score + g2m.score, data = data4lm)
  }
  #return(fit)
  
  #some betas are NA [change so that they're the mean effect across other values with the same variable]
  rtn <- list(residuals_var = var(fit$residuals), coefficients = fit$coefficients)
  rtn$coefficients[which(is.na(rtn$coefficients))] <- 0
  
  #n_perturbed
  if (ExperimentName == "Magpie"){
    n_control = length(inlets_combined$Moderate$orig.ident) + length(inlets_combined$Strong$orig.ident)
    rtn$n_control <- n_control
    #control_norm_expr
    rtn$control_num_reads <- sum(inlets_combined$Moderate[["RNA"]]@counts[genes2test[i],]) + sum(inlets_combined$Strong[["RNA"]]@counts[genes2test[i],])
    rtn$control_norm_expr <- (sum(inlets_combined$Moderate[["RNA"]]@data[genes2test[i],]) + sum(inlets_combined$Strong[["RNA"]]@data[genes2test[i],]))/n_control
    
    rtn$av_log_expression_control <- log((sum(expm1(as.matrix(inlets_combined$Moderate[["RNA"]]@data[genes2test[i], ]))) +
                                            sum(expm1(as.matrix(inlets_combined$Strong[["RNA"]]@data[genes2test[i], ]))))/n_control + 1/scale_factor, base = lfc_base) 
    
  } else {
    n_control = length(inlets_combined$orig.ident)
    rtn$n_control <- n_control
    #control_norm_expr
    rtn$control_num_reads <- sum(inlets_combined[["RNA"]]@counts[genes2test[i],])
    rtn$control_norm_expr <- sum(inlets_combined[["RNA"]]@data[genes2test[i],])/n_control
    
    rtn$av_log_expression_control <- log(sum(expm1(as.matrix(inlets_combined[["RNA"]]@data[genes2test[i], ]))) /n_control + 1/scale_factor, base = lfc_base) 
    
  }
  
  #control_num_reads
  
  #write 
  if (scale_factor != 1){
    fnm = paste0("sf-", scale_factor, "_", date, "_gene-", gsub(gene_name, pattern = "/", replacement = "-"), ".RDS")
  } else {
    if (ExperimentName == "Replogle"){
      fnm = paste0(date, "_gene-", gsub(gene_name, pattern = "/", replacement = "-"), "_", cell_type, ".RDS")
    } else {
      fnm = paste0(date, "_gene-", gsub(gene_name, pattern = "/", replacement = "-"), ".RDS")
    }
    
  }
  saveRDS(rtn, file = file.path(outdir, fnm))
  #print(paste0(gene_name, "done!"))
}
#names(control_lm) <- genes2test

sessionInfo()


