suppressWarnings(suppressMessages(library(Matrix, warn.conflicts = F)))
suppressWarnings(suppressMessages(library(Seurat, warn.conflicts = F)))
suppressWarnings(suppressMessages(library(data.table, warn.conflicts = F)))
suppressWarnings(suppressMessages(library(dplyr, warn.conflicts = F)))
suppressWarnings(suppressMessages(library(ggplot2, warn.conflicts = F)))
suppressWarnings(suppressMessages(library(tidyverse, warn.conflicts = F)))
suppressWarnings(suppressMessages(library(parallel, warn.conflicts = F)))
suppressWarnings(suppressMessages(library(doParallel, warn.conflicts = F)))




## for magpie
#date <- "2022-08-15"
#HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
#ExperimentName <- "Magpie"
#magpie_version <- "2022-08-15"
#perturbed_cell_infiles <- "/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Magpie/5a_combine/2022-08-15_Strong_perturbed_cells_combined.RDS"
#perturbed_cell_infiles <- "/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Magpie/5a_combine/"
#grep_pattern <- "unassigned_cells_combined[.]RDS"
#control_lm_folder <- "/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Magpie/6a_run_control_lm/"
#outdir <- "/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Magpie/13a_calc_lfcs_transcriptome_wide_by_complex/by_gene/"
#subsample <- T
#num_cells <- 25
#rnd.seed <- 0
#by_gene <- T
#by_guide <- F
#by_timepoint <- F



## for pica
#date <- "2022-10-05"
#date_out <- "2024-05-17"
#HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
#ExperimentName <- "Pica"
## control cells
#control_lm_folder <- "/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Pica/6a_run_control_lm/2024-05-17/"
#section_name <- "6f_calc_lfcs_transcriptome_wide_by_gene_per_line"
#outdir <- "/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Pica/6f_calc_lfcs_transcriptome_wide_by_gene_per_line/fiaj_3/by_target_by_line/"
#target_gene <- "unassigned"
#guide_name <- 'MBD3_AGCGGAGGCAGCTGAGGCGG'
#cell_line <- "tolg_4"; excl_line <- F
#perturbed_cell_infiles <- paste0("/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Pica/5a_combine/2022-10-05_", cell_line, "_perturbed_cells.RDS")
#by_gene <- F
#by_guide <- T
#by_timepoint <- F

#for replogle
#date <- "2022-09-12"
#HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
#ExperimentName <- "Pica"
#magpie_version <- "2022-08-15"
#perturbed_cell_infiles <- "/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Other/Replogle/2a_make_sc_seurat_object/RPE1_raw_sc.RDS"
#control_lm_folder <- "/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Other/Replogle/2b_run_control_lm/"
#section_name <- "2d_calc_lfcs_transcriptome_wide_by_complex"
#outdir <- "/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Other/Replogle/2d_calc_lfcs_transcriptome_wide_by_complex/by_gene/"
#by_complex <- T
#complex_ind <- 1
#by_gene <- F
#by_guide <- F
#by_timepoint <- F

## for pica and magpie
## this will typically be for more downstream genes

by_gene <- T
by_line <- F
excl_line <- F
by_donor <- F
by_complex <- F
by_timepoint <- F
by_guide <- F

scale_factor <- 1
lfc_base <- 10
subsample <- F
grep_pattern <- "perturbed_cells"
subsample_by_target <- F
max_cells <- 2000


args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  ##needed arguments
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--date_out"){ date_out <- args[[ind + 1]] }
  if (arg == "--magpie_version"){ magpie_version <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--experiment_name"){ ExperimentName <- args[[ind + 1]] }
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
  if (arg == "--target_gene"){ target_gene <- args[[ind + 1]] } #target not required if you're looking at complexes
  if (arg == "--guide_name"){ guide_name <- args[[ind + 1]] }
  
  ## in/out, outfile optional
  if (arg == "--perturbed_cell_infiles"){  perturbed_cell_infiles <- args[[ind + 1]]}
  if (arg == "--control_lm_folder"){ control_lm_folder <- args[[ind + 1]] }
  if (arg == "--outdir"){ outdir <- args[[ind + 1]] }
  if (arg == "--inlet_strength"){ inlet_strength <- args[[ind + 1]] }
  if (arg == "--grep_pattern"){ grep_pattern <- args[[ind + 1]] } #target not required if you're looking at complexes
  
  
  ### settings (on what you're testing against)
  ### optionality is dependent on the settings themselves
  if (arg == "--by_gene"){ by_gene <- as.logical(args[[ind + 1]]) }
  if (arg == "--by_guide"){ by_guide <- as.logical(args[[ind + 1]]) }
  if (arg == "--by_line"){ by_line <- as.logical(args[[ind + 1]]) }
  if (arg == "--excl_line"){ excl_line <- as.logical(args[[ind + 1]]) }
  if (arg == "--cell_line"){ cell_line <- args[[ind + 1]] }
  if (arg == "--by_donor"){ by_donor <- as.logical(args[[ind + 1]]) }
  if (arg == "--donor"){ donor <- args[[ind + 1]] }
  if (arg == "--by_timepoint"){ by_timepoint <- as.logical(args[[ind + 1]]) }
  if (arg == "--timepoint"){ timepoint <- as.numeric(args[[ind + 1]]) }
  if (arg == "--by_complex"){ by_complex <- as.logical(args[[ind + 1]]) }
  if (arg == "--complex_ind"){ complex_ind <- as.numeric(args[[ind + 1]]) }
  if (arg == "--complex_meta_path"){ complex_meta_path <- args[[ind + 1]] }
  if (arg == "--subsample"){ subsample <- as.logical(args[[ind + 1]]) }
  if (arg == "--num_cells"){ num_cells <- as.numeric(args[[ind + 1]]) }
  if (arg == "--max_cells"){ max_cells <- as.numeric(args[[ind + 1]]) }
  if (arg == "--rnd_seed"){ rnd.seed <- as.numeric(args[[ind + 1]]) }
  if (arg == "--subsample_by_target"){ subsample_by_target <- as.logical(args[[ind + 1]]) }
  if (arg == "--subsample_high"){ subsample_high <- as.logical(args[[ind + 1]]) }
  
  
  ### not needed
  if (arg == "--utils_path"){ utils_path <- args[[ind + 1]] }
  if (arg == "--io_path"){ io_path <- args[[ind + 1]] }
  if (arg == "--lfc_base"){ lfc_base <- as.numeric(args[[ind + 1]]) }
  if (arg == "--scale_factor"){ scale_factor <- as.numeric(args[[ind + 1]]) }
  if (arg == "--max_cells"){ max_cells <- as.numeric(args[[ind + 1]]) }
  if (arg == "--min_cells_per_gene"){ min_cells_per_gene1 <- as.numeric(args[[ind + 1]]) }
  if (arg == "--min_genes_per_complex"){ min_genes_per_complex1 <- as.numeric(args[[ind + 1]]) }
  if (arg == "--excl_line"){ excl_line <- as.logical(args[[ind + 1]]) }
  if (arg == "--n_core"){ n_core <- as.numeric(args[[ind + 1]]) }
  if (arg == "--redo"){ redo <- args[[ind + 1]] }
}

#scale_factor <- scale_factor*10^5
setwd(HomeFolder)
if (!(exists("io_path"))){io_path <- paste0(HomeFolder, "scripts/io/", ExperimentName, "_io.R")}
if (!(exists("utils_path"))){ utils_path <- paste0(HomeFolder, "/scripts/io/Magpie_Utils.R")}
source(io_path)
source(utils_path)
control_tag <- c("unassigned", "NonTarget") # use all unassigned cells as control

if (exists("min_cells_per_gene1")){ min_cells_per_gene <- min_cells_per_gene1}
if (exists("min_genes_per_complex1")){ min_genes_per_complex <- min_genes_per_complex1}

print(paste("date:", date))
print(paste("home folder:", HomeFolder))
print(paste("Experiment Name:", ExperimentName))
print(paste("min_cells_per_gene:", min_cells_per_gene))
print(paste("section_name:", section_name))


## set io
if (!(exists("outdir"))){paste0("no outdir specified!")} else {print(paste("writing output to:", outdir))}
if(!dir.exists(outdir)) dir.create(outdir, recursive = T)
#plotsdir <- file.path(HTMLFolder, "pipeline", section_name)
#if(!dir.exists(plotsdir)) dir.create(plotsdir, recursive = T)




## protein complex in question
if (by_complex == TRUE){
  if (!(exists("complex_meta_path"))){
    complex_meta_path <- sort(grep(list.files(file.path(ProjectFolder, "outs", "Magpie", "Preprocess_External_Datasets", "omnipath"), pattern = "complex_meta"),
                                   pattern = magpie_version, value = T), decreasing = T)[1]
    complex_meta <- fread(file.path(ProjectFolder, "outs", "Magpie", "Preprocess_External_Datasets", "omnipath", complex_meta_path))
  } else{
    complex_meta <- fread(complex_meta_path)
  }
  complex_genes <- unlist(strsplit(complex_meta$genes_in_complex[complex_ind], "_"))
  print(complex_genes)
} else{
  complex_genes <- 1:(min_genes_per_complex+1)
}
num_genes_in_complex <- length(complex_genes)


## if by guide
if (by_guide == T){
  print(guide_name)
  guide_name <- gsub(guide_name, pattern = "-", replacement = "_")
  target_gene <- get_target_gene(guide_name)
} else if (by_complex == T){
  guide_name <- target_gene <- paste0(complex_genes, collapse = "_")
} else {
  target_gene <- gsub(target_gene, pattern = "-", replacement = "_")
  guide_name <- target_gene
  print(target_gene)
}



# load data
## check if it's a folder or a single file
if ( (by_complex == T & num_genes_in_complex < min_genes_per_complex) ){
  assigned_cells <- NULL
  n_perturbed <- 0
} else 
{
  ## if you were given a file name rather than a folder
  if (substr(perturbed_cell_infiles, nchar(perturbed_cell_infiles) - 3, nchar(perturbed_cell_infiles)) == ".RDS"){
    print('reading data from')
    print(perturbed_cell_infiles)
    
    assigned_cells <- readRDS(file = perturbed_cell_infiles)
    
    if (by_timepoint == T){
      if (length(which(assigned_cells$Timepoint == paste0("D", timepoint))) == 0){
        assigned_cells <- NULL
      } else {
        assigned_cells <- subset(assigned_cells, Timepoint == paste0("D", timepoint))
      }
    }
    if (by_guide == T){
      if (length(which(assigned_cells$feature_call == guide_name)) == 0){
        assigned_cells <- NULL
      } else {
        assigned_cells <- subset(assigned_cells, feature_call == guide_name)
      }
    } else if (by_complex == T){
      if (length(which(assigned_cells$target_gene_name %in% complex_genes)) == 0){
        assigned_cells <- NULL
      } else {
        assigned_cells <- subset(assigned_cells, target_gene_name %in% complex_genes)
      }
    } else { #then it'll be by target gene
      if (length(which(assigned_cells$target_gene_name == target_gene)) == 0){
        assigned_cells <- NULL
      } else {
        assigned_cells <- subset(assigned_cells, target_gene_name == target_gene)
      }
    }
  } else {
    
    assigned_cells_paths <- grep(grep(file.path(perturbed_cell_infiles, list.files(perturbed_cell_infiles, pattern = grep_pattern)),
                                      pattern = date, value = T),
                                 pattern = "[.]RDS", value = T)
    print('reading data from')
    print(assigned_cells_paths)
    if (excl_line == T){
      assigned_cells_paths <- grep(assigned_cells_paths, pattern = cell_line, invert = T, value = T )
    }
    assigned_cells <- mapply(assigned_cells_paths, FUN = function(f){
      #print(f)
      rtn <- readRDS(file.path(f))
      
      if (by_timepoint == T){
        if (length(which(rtn$Timepoint == paste0("D", timepoint))) > 1){
          rtn <- subset(rtn, Timepoint == paste0("D", timepoint))
        } else {
          rtn <- NULL
        }
      }
      if (by_guide == T){
        if (length(which(rtn$feature_call == guide_name)) > 1){
          rtn <- subset(rtn, feature_call == guide_name)
        } else {
          rtn <- NULL
        }
      } else if (by_complex == T){
        if (length(which( rtn$target_gene_name %in% complex_genes)) > 1){
          rtn <- subset(rtn, target_gene_name %in% complex_genes)
        } else {
          rtn <- NULL
        }
      } else { #then it'll be by target gene
        if (length(which( rtn$target_gene_name == target_gene )) > 1){
          rtn <- subset(rtn, target_gene_name == target_gene)
        } else {
          rtn <- NULL
        }
      }
      
      return(rtn)
      
    }, SIMPLIFY = F)
    
    names(assigned_cells) <- gsub(assigned_cells_paths, pattern = paste0(".*", date, "_|_perturbed_cells.*"), replacement = "")
    assigned_cells[sapply(assigned_cells, is.null)] <- NULL
    if (length(assigned_cells) > 1){
      assigned_cells <- merge(x = assigned_cells[[1]],
                              y = assigned_cells[-1],
                              add.cell.ids = names(assigned_cells),
                              project = paste(ExperimentName, "combined", sep = "_"),
                              merge.data = TRUE)
    } else if (length(assigned_cells) == 1){
      assigned_cells <- assigned_cells[[1]]
    } else {
      assigned_cells <- NULL
    }
  }
  
  print('data loaded!')
  ifelse(!(is.null(assigned_cells)), print(paste0(dim(assigned_cells@meta.data)[1], " cells")), print('0 cells') )
  
  if (is.null(assigned_cells)){n_perturbed <- 0} else {
    n_perturbed <- dim(assigned_cells@meta.data)[1]
    if (by_complex == T){
      cells_per_knockdown <- table(assigned_cells$target_gene_name)
      num_genes_in_complex <- length(which(cells_per_knockdown > min_cells_per_guide))
    }
  }
}

if (n_perturbed > 0){
  if (dim(assigned_cells@meta.data)[1] > max_cells){
    subsample <- T; num_cells <- max_cells; subsample_by_target <- F; rnd.seed <- 0
  }
}
if (subsample == T){
  all_barcodes <- row.names(assigned_cells@meta.data)
  if(exists("rnd.seed")){set.seed(rnd.seed)}
  if (subsample_by_target == F){
    barcodes2keep <- sample(all_barcodes, size = num_cells, replace = T)
  } else {
    target_gene_ind <- which(unlist(lapply(strsplit(row.names( assigned_cells), ":" ), "[[", 2)) == target_gene)
    target_gene_expr <- as.vector(assigned_cells[["RNA"]][target_gene_ind,])
    guide_inds <- which(unlist(lapply(strsplit(row.names( assigned_cells[["guides"]]), "-" ), "[[", 1)) == target_gene)
    n_non_zero <- length(which(target_gene_expr > 0))
    barcodes2keep <- colnames(assigned_cells[["RNA"]])[order(target_gene_expr, decreasing = subsample_high)[1:min(num_cells, n_non_zero)]]
    expr <- data.frame(
      cell_barcode = barcodes2keep,
      target_expr = as.vector(assigned_cells[["RNA"]][target_gene_ind, barcodes2keep]),
      t(assigned_cells[["RNA"]]@counts[c("BSD:BSD:Gene-Expression", "dCas9-KRAB-MeCP2:dCas9-KRAB-MeCP2:Gene-Expression", "mScarlet:mScarlet:Gene-Expression"),
                              barcodes2keep]),
      total_guides = assigned_cells$nCount_guides[barcodes2keep],
      t(assigned_cells[["guides"]][guide_inds, barcodes2keep])
    )
    saveRDS(expr, paste0(outdir, "/", date, "_", target_gene, "_n-cells-", num_cells, "_subsample-", ifelse(subsample_high == T, "high", "low"), ".RDS"))
  }
  assigned_cells <- subset(assigned_cells, cells = barcodes2keep)
}
#

## ---- execute lm


#run lm
if ( ( n_perturbed < min_cells_per_gene) | (by_complex == T & num_genes_in_complex < min_genes_per_complex)){
  print("minimum cells per gene/guide:")
  print(min_cells_per_gene)
  print(paste0("number of cells: ", n_perturbed))
  print("minimum genes per complex:")
  print(min_genes_per_complex)
  print(paste0("number of genes: ", length(num_genes_in_complex)))
  print("Not enough knockdown cells! (or protein complex too small)")
  #file.create(file = file.path(outdir, fnm))
  ## write file with table on why target wasn't included
} else {
  
  #Number of genes expressed
  print(paste("Number of genes expressed:", dim(assigned_cells[["RNA"]]@data)[1]))
  
  # calc LFS 
  print(paste("Computing downstream LFCs for", guide_name))
  
  
  #run lm
  print("running LM!")
  
  genes2test <- gsub(list.files(control_lm_folder), pattern = "^gene-|[.]RDS", replacement = "")
  adjusted_row_names <- gsub(rownames(assigned_cells), pattern = ":|/", replacement = "-")
  genes2test <- rownames(assigned_cells)[match(genes2test, adjusted_row_names)]
  
  
  #also keep only genes where there is a lm run for the control cells
  print(paste("Testing", length(genes2test), "genes for perturbation effects"))
  
  # renormalize
  #log(t(t(counts_matrix)*10000/inlets_combined$nCount_RNA) + ps)
  if (scale_factor != 1) {
    assigned_cells <- NormalizeData(assigned_cells, scale.factor = scale_factor*10^5)
    print("seurat object renormalized!")
    normalized_expression <- assigned_cells@assays$RNA@data[genes2test,]
  } else {
    normalized_expression <- assigned_cells@assays$RNA@data[genes2test,]
  }
  
  #cell metadata for perturbed cells
  perturbed_cell_metadata <- data.frame(cell_line = c(assigned_cells$cell_line),
                                        inlet = c(assigned_cells$orig.ident),
                                        nCount_RNA = c(assigned_cells$nCount_RNA),
                                        percent_MT = c(assigned_cells$percent_MT),
                                        s.score = c(assigned_cells$S.Score),
                                        g2m.score = c(assigned_cells$G2M.Score)
  )
  if (exists("n_core")){ n_core <- min(n_core, parallel::detectCores(), length(genes2test)) } else {n_core <- min(parallel::detectCores(), length(genes2test))}
  
  registerDoParallel(n_core)
  lm_out <- foreach (i = 1:length(genes2test)) %dopar% {
    #print(i)
    gene2test <- genes2test[i]
    
    control_fit_path <- paste0("gene-", gsub(gene2test, pattern = ":|/", replacement = "-"), ".RDS")
    #print(paste0('reading control model from ', file.path(control_lm_folder, control_fit_path)))
    control_lm <- readRDS(file.path(control_lm_folder, control_fit_path))
    reformatted_perturbed_cell_metadata <- reformat_metadata_matrix(in_matrix = perturbed_cell_metadata, coefficients = control_lm$coefficients)
    lm_res <- run_lm(perturbed_normalized_expression = normalized_expression[gene2test,], 
                     perturbed_cell_metadata = reformatted_perturbed_cell_metadata,
                     control_fit = control_lm)
    lm_res$n_control = control_lm$n_control
    lm_res$control_num_reads = control_lm$control_num_reads
    lm_res$control_norm_expr = control_lm$control_norm_expr
    lm_res$av_log_expression_control = control_lm$av_log_expression_control
    
    if (by_line == T){
      lm_res$control_line_coef <- ifelse(paste0('cell_line', cell_line) %in% names(control_lm$coefficients), control_lm$coefficients[paste0('cell_line', cell_line)], 0)
    }
    
    rm(control_lm); gc()
    return(lm_res)
  }
  
  names(lm_out) <- genes2test
  if (by_line == T){
    lm_fit <- as.data.frame(t(bind_rows(lapply(lm_out, 
                                               FUN = function(x){
                                                 return(c(x$lfc, x$v, x$pval, x$n_control, x$control_num_reads, x$control_norm_expr, x$av_log_expression_control, x$control_line_coef))
                                               }))))
    names(lm_fit) <- c("lfc", "v", "pval", "n_control", "control_num_reads", "control_norm_expr", "av_log_expression_control", "control_line_coef")
    
  } else {
    lm_fit <- as.data.frame(t(bind_rows(lapply(lm_out, 
                                               FUN = function(x){
                                                 return(c(x$lfc, x$v, x$pval, x$n_control, x$control_num_reads, x$control_norm_expr, x$av_log_expression_control))
                                               }))))
    names(lm_fit) <- c("lfc", "v", "pval", "n_control", "control_num_reads", "control_norm_expr", "av_log_expression_control")
  }
  
  #print residuals
  cnms <- c("id", names(lm_out))
  lm_res <- data.frame(
    id = row.names(perturbed_cell_metadata),
    as.data.frame(bind_cols(mapply(lm_out, FUN = function(x){return(x$res)})))
  )
  colnames(lm_res) <- cnms
  
  if (exists("inlet_strength")){
    fwrite(lm_res, file = file.path(outdir, paste0(date_out, "_", inlet_strength, "_", gsub(guide_name, pattern = "/", replacement = "-"), "_knockdowns_residuals.tsv.gz")))
  } else if (by_complex == F & subsample == F){
    fwrite(lm_res, file = file.path(outdir, paste0(date_out, "_", gsub(guide_name, pattern = "/", replacement = "-"), "_knockdowns_residuals.tsv.gz")))
  } else if (subsample == T) {
    if (subsample_by_target == T){
      fnm <- paste0(gsub(target_gene, pattern = "/", replacement = "_"), "_", date_out, "_", num_cells, "-cells_subsample-", ifelse(subsample_high == T, "high", "low"), "_knockdowns_residuals.tsv.gz")
    } else {
      fnm <- paste0(gsub(target_gene, pattern = "/", replacement = "_"), "_", date_out, "_", num_cells, "-cells_", rnd.seed, "_knockdowns_residuals.tsv.gz")
    }
    fwrite(lm_res, row.names = T, file = file.path(outdir, fnm))
  } else{
    fwrite(lm_res, row.names = T, file = file.path(outdir, paste0(date_out, "_", complex_ind, "_knockdowns_residuals.tsv.gz")))
  }
  
  
  
  #row.names(lm_res) <- genes2test
  print("LM complete!")
  
  ## other stuff
  lm_fit$perturbed_norm_expr <- rowMeans(normalized_expression[genes2test,])
  lm_fit$av_log_expression_perturbed <- log(rowMeans(expm1(normalized_expression[genes2test,])) + 1/scale_factor, base = lfc_base)
  lm_fit$lfc_unadjusted <- lm_fit$av_log_expression_perturbed - lm_fit$av_log_expression_control
  lm_fit$perturbed_num_reads <- rowSums(assigned_cells@assays$RNA@counts[genes2test,])
  
  print("writing df!")
  rtn <- data.frame(downstream_gene = genes2test,
                    n_control = lm_fit$n_control,
                    n_perturbed = n_perturbed,
                    lfc = lm_fit$lfc,
                    lfc_unadjusted = lm_fit$lfc_unadjusted,
                    v = lm_fit$v,
                    pval_lm = lm_fit$pval,
                    control_norm_expr = lm_fit$control_norm_expr,
                    perturbed_norm_expr = lm_fit$perturbed_norm_expr,
                    control_num_reads = lm_fit$control_num_reads,
                    perturbed_num_reads = lm_fit$perturbed_num_reads)
  
  if (by_guide == T){
    print('adding guide and target')
    rtn <- mutate(rtn, guide = guide_name)
    rtn <- mutate(rtn, target = target_gene)
  } 
  if (by_complex == T){
    print('adding complex')
    rtn <- mutate(rtn, genes_in_complex = paste(complex_genes, collapse = "_"))
  } 
  if (by_timepoint == T){
    print('adding timepoint')
    rtn <- mutate(rtn, timepoint = paste0("D", timepoint))
  }
  if (by_gene == T){
    print('adding target')
    rtn <- mutate(rtn, target = target_gene)
  }
  if (by_line == T){
    print('adding line')
    rtn <- mutate(rtn, control_line_coef = lm_fit$control_line_coef, cell_line = cell_line)
  }
  
  if (exists("inlet_strength")){
    guide_name <- paste0(inlet_strength, "_", guide_name)
    target_gene <- paste0(inlet_strength, "_", target_gene)
  }
  if (by_guide == T) {
    fnm <- paste0(gsub(guide_name, pattern = "/", replacement = "_"), "_", date_out, ".tsv.gz")
  } else if (by_complex == T){
    fnm <- paste0(complex_ind, "_", date_out, ".tsv.gz")
  } else {
    if (subsample == T){
      if (subsample_by_target == T){
        fnm <- paste0(gsub(target_gene, pattern = "/", replacement = "_"), "_", date_out, "_", num_cells, "-cells_subsample-", ifelse(subsample_high == T, "high", "low"), ".tsv.gz")
      } else {
        fnm <- paste0(gsub(target_gene, pattern = "/", replacement = "_"), "_", date_out, "_", num_cells, "-cells_", rnd.seed, ".tsv.gz")
      }
    } else {
      fnm <- paste0(gsub(target_gene, pattern = "/", replacement = "_"), "_", date_out, ".tsv.gz")
    }
  }
  #write
  print(paste("writing df to:", file.path(outdir, fnm)))
  fwrite(rtn, file = file.path(paste0(outdir, fnm)), row.names = F, quote = F, sep = "\t", compress = "gzip")
  print("df written!")
}

#sessionInfo()

