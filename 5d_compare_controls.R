library(Matrix)
library(data.table)
library(dplyr)
library(tidyverse)
library(Seurat)
library(parallel)
# print(paste("Number of cores:", parallel::detectCores()))

date <- "2022-10-05"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
guide_group <- NA # "Moderate" #"Strong", "All_Inlets"
ExperimentName <- "Pica"
section_name <- "5d_compare_controls"




# options
args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--guide_group"){ guide_group <- args[[ind + 1]] }
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
  if (arg == "--ncores"){ ncores <- args[[ind + 1]] } 
  if (arg == "--experiment_name"){ ExperimentName <- args[[ind + 1]] } 
}

# i/o
source(file.path(HomeFolder, paste0("scripts/io/", ExperimentName, "_io.R")))
outdir <- file.path(OutFolder, section_name)
plotsdir <- file.path(HTMLFolder, "pipeline/", section_name)
if(!dir.exists(outdir)) dir.create(outdir)
if(!dir.exists(plotsdir)) dir.create(plotsdir)
guide_groups <- NA #c("Moderate", "Strong")

GuideMetadata <- fread(GuideMetadataPath)
outfile <- file.path(outdir, paste0(date, "_", guide_group,
                                    "_comparison_controls.csv"))
if (file.exists(outfile)){
  all_res <- fread(outfile)
} else {
  
  # load data
  if (ExperimentName == "Magpie"){
    if (guide_group == "All_Inlets"){
      unassigned <- mapply(guide_groups, FUN = function(guide_group){
        readRDS(file.path(OutFolder, "5a_combine",
                          paste0(date, "_", guide_group, "_unassigned_cells_combined.RDS")))
      })
      nontarget <- mapply(guide_groups, FUN = function(guide_group){
        readRDS(file.path(OutFolder, "5a_combine",
                          paste0(date, "_", guide_group, "_nontarget_cells_combined.RDS")))
      })
      perturbed <- mapply(guide_groups, FUN = function(guide_group){
        readRDS(file.path(OutFolder, "5a_combine",
                          paste0(date, "_", guide_group, "_perturbed_cells_combined.RDS")))
      })
      
      cell_metadata <- as.data.frame(bind_rows(
        unassigned[[1]]@meta.data,
        unassigned[[2]]@meta.data,
        nontarget[[1]]@meta.data,
        nontarget[[2]]@meta.data,
        perturbed[[1]]@meta.data,
        perturbed[[2]]@meta.data
      )) %>% 
        mutate(ident = c(rep("unassigned", dim(unassigned[[1]]@assays$RNA@data)[2] + dim(unassigned[[2]]@assays$RNA@data)[2]),
                         rep("nontarget", dim(nontarget[[1]]@assays$RNA@data)[2] + dim(nontarget[[2]]@assays$RNA@data)[2]),
                         rep("perturbed", dim(perturbed[[1]]@assays$RNA@data)[2] + dim(perturbed[[2]]@assays$RNA@data)[2])
        ))
      
      # unassigned cells
      print(paste("Number of unassigned cells:", dim(unassigned[[1]]@meta.data)[1] + dim(unassigned[[2]]@meta.data)[1]))
      
      # control cells
      print(paste("Number of cells with NonTargeting guides:", dim(nontarget[[1]]@meta.data)[1] + dim(nontarget[[2]]@meta.data)[1]))
      #stopifnot(sum(NTs) > 3)
      
      # test all genes (min. expression)
      genes_to_test <- rownames(perturbed[[1]])
      print(paste("Number of genes tested:", length(genes_to_test)))
      
      
    } else {
      unassigned <- readRDS(file.path(OutFolder, "5a_combine",
                                      paste0(date, "_", guide_group, "_unassigned_cells_combined.RDS")))
      nontarget <- readRDS(file.path(OutFolder, "5a_combine",
                                     paste0(date, "_", guide_group, "_nontarget_cells_combined.RDS")))
      perturbed <- readRDS(file.path(OutFolder, "5a_combine",
                                     paste0(date, "_", guide_group, "_perturbed_cells_combined.RDS")))
      
      cell_metadata <- as.data.frame(bind_rows(unassigned@meta.data,
                                               nontarget@meta.data,
                                               perturbed@meta.data)) %>% 
        mutate(ident = c(rep("unassigned", dim(unassigned@assays$RNA@data)[2]),
                         rep("nontarget", dim(nontarget@assays$RNA@data)[2]),
                         rep("perturbed", dim(perturbed@assays$RNA@data)[2])))
      
      # unassigned cells
      print(paste("Number of unassigned cells:", dim(unassigned@meta.data)[2]))
      
      # control cells
      print(paste("Number of cells with NonTargeting guides:", dim(nontarget@meta.data)[2]))
      #stopifnot(sum(NTs) > 3)
      
      # test all genes (min. expression)
      genes_to_test <- rownames(perturbed)
      print(paste("Number of genes tested:", length(genes_to_test)))
      
      
    }
    
    
  }
  # load data
  if (ExperimentName == "Pica"){
    
    unassigned <- mapply(guide_groups, FUN = function(guide_group){
      readRDS(file.path(OutFolder, "5a_combine",
                        paste0(date, "_", guide_group, "_unassigned_cells_combined.RDS")))
    })
    nontarget <- mapply(guide_groups, FUN = function(guide_group){
      readRDS(file.path(OutFolder, "5a_combine",
                        paste0(date, "_", guide_group, "_nontarget_cells_combined.RDS")))
    })
    perturbed <- mapply(guide_groups, FUN = function(guide_group){
      readRDS(file.path(OutFolder, "5a_combine",
                        paste0(date, "_", guide_group, "_perturbed_cells_combined.RDS")))
    })
    
    cell_metadata <- as.data.frame(bind_rows(
      unassigned[[1]]@meta.data,
      unassigned[[2]]@meta.data,
      nontarget[[1]]@meta.data,
      nontarget[[2]]@meta.data,
      perturbed[[1]]@meta.data,
      perturbed[[2]]@meta.data
    )) %>% 
      mutate(ident = c(rep("unassigned", dim(unassigned[[1]]@assays$RNA@data)[2] + dim(unassigned[[2]]@assays$RNA@data)[2]),
                       rep("nontarget", dim(nontarget[[1]]@assays$RNA@data)[2] + dim(nontarget[[2]]@assays$RNA@data)[2]),
                       rep("perturbed", dim(perturbed[[1]]@assays$RNA@data)[2] + dim(perturbed[[2]]@assays$RNA@data)[2])
      ))
    
    # unassigned cells
    print(paste("Number of unassigned cells:", dim(unassigned[[1]]@meta.data)[1] + dim(unassigned[[2]]@meta.data)[1]))
    
    # control cells
    print(paste("Number of cells with NonTargeting guides:", dim(nontarget[[1]]@meta.data)[1] + dim(nontarget[[2]]@meta.data)[1]))
    #stopifnot(sum(NTs) > 3)
    
    # test all genes (min. expression)
    genes_to_test <- rownames(perturbed[[1]])
    print(paste("Number of genes tested:", length(genes_to_test)))
    
    
    
  }
  
  
  # test for differences
  all_res <- lapply(c("NTs_vs_unassigned", "assigned_vs_unassigned", "NTs_vs_perturbed"), function(comp) {
    print(comp)
    res <- mclapply(genes_to_test, function(gene2test){
      print(gene2test)
      #y is the expression vector for the target gene in question
      data4lm <- cell_metadata
      if (guide_group == "All_Inlets"){
        data4lm$expr <- c(unassigned[[1]]@assays$RNA@data[gene2test,],unassigned[[2]]@assays$RNA@data[gene2test,],
                          nontarget[[1]]@assays$RNA@data[gene2test,], nontarget[[2]]@assays$RNA@data[gene2test,],
                          perturbed[[1]]@assays$RNA@data[gene2test,], perturbed[[2]]@assays$RNA@data[gene2test,])
      } else {
        data4lm$expr <- c(unassigned@assays$RNA@data[gene2test,],
                          nontarget@assays$RNA@data[gene2test,],
                          perturbed@assays$RNA@data[gene2test,])
      }
      
      
      if(comp == "NTs_vs_unassigned"){
        data4lm <- data4lm %>% filter(ident %in% c("nontarget", "unassigned"))
        data4lm$to_test <- 0
        data4lm$to_test[which(data4lm$ident == "unassigned")] <- 1
        ident1_name <- "nontarget"
        ident2_name <- "unassigned"
      } else if (comp == "assigned_vs_unassigned"){
        data4lm <- data4lm
        data4lm$to_test <- 0
        data4lm$to_test[which(data4lm$ident == "unassigned")] <- 1
        ident1_name <- "assigned"
        ident2_name <- "unassigned"
      } else if (comp == "NTs_vs_perturbed"){
        data4lm <- data4lm %>% filter(ident %in% c("nontarget", "perturbed"))
        data4lm$to_test <- 0
        data4lm$to_test[which(data4lm$ident == "perturbed")] <- 1
        ident1_name <- "nontarget"
        ident2_name <- "perturbed"
      } else {
        stop("Comparison type not defined.")
      }
      
      
      fit <- lm(expr ~ cell_line + S.Score + G2M.Score + orig.ident + to_test + nCount_RNA + percent_MT, data = data4lm)
      beta <- coef(fit)["to_test"]
      pval <- summary(fit)$coefficients["to_test",4]
      data.frame(lfc = beta,
                 pval_lm = pval,
                 n_ident1 = length(which(data4lm$to_test == 0)),
                 n_ident2 = length(which(data4lm$to_test == 1)),
                 gene = gene2test,
                 mean_expr = mean(data4lm$expr),
                 inlets = length(unique(data4lm$orig.ident)),
                 ident2_name = ident2_name,
                 ident1_name = ident1_name)
    }, mc.cores = ncores) %>% bind_rows()
    
    res$pval_adj <- p.adjust(res$pval_lm, method = "BH")
    res$comp <- comp
    res
  }) %>% bind_rows()
  
  # save results
  write.csv(all_res, file.path(outdir, paste0(date, "_", guide_group,
                                              "_comparison_controls.csv")))
  
  
}


all_res$comp <- factor(all_res$comp)
all_res  <- all_res %>% mutate(downstream_gene_name = unlist(lapply(strsplit(all_res$gene, ":"), "[[", 2)))
levels(all_res$comp) <- c("Assigned vs. Unassigned", "Non-target vs. Perturbed", "Non-Target vs. Unassigned")

# make plot
pdf(file.path(plotsdir,
              paste0(date, "_", guide_group,
                    "_comparison_controls.pdf")), width = 10, height = 4)
print(ggplot(all_res,  aes(x = lfc, y= -log10(pval_lm),
                           label = ifelse(pval_adj < sig_pval_thresh, downstream_gene_name, ""))) +
        ggrepel::geom_label_repel() +
        xlab("LFC (Natural Log)") + ylab("p-value (-log 10)") + ggtitle(paste0(guide_group, " Cells")) + 
        facet_wrap(~comp) + geom_point() + theme_bw() + geom_hline(yintercept = -log10(0.05), col = "red", lty = 3))
print(ggplot(all_res,  aes(x = lfc, y= -log10(pval_adj))) +
        xlab("LFC (Natural Log)") + ylab("Adjusted p-value (-log 10)") + ggtitle(paste0(guide_group, " Cells")) + 
        facet_wrap(~comp) + geom_point() + theme_bw() + geom_hline(yintercept = -log10(0.05), col = "red", lty = 3))
print(ggplot(all_res %>% filter(downstream_gene_name %in% GuideMetadata$gene),  aes(x = lfc, y= -log10(pval_adj))) +
        xlab("LFC (Natural Log)") + ylab("Adjusted p-value (-log 10)") + ggtitle(paste0(guide_group, " Cells (Targets only)")) + 
        facet_wrap(~comp) + geom_point() + theme_bw() + geom_hline(yintercept = -log10(0.05), col = "red", lty = 3))
dev.off()

# print out
print("Differentially expressed genes (non-target vs. unassigned):")
all_res %>% filter(pval_adj < sig_pval_thresh) %>%
  filter(comp == "Non-target vs. Unassigned") %>%
  select(lfc, gene, pval_lm, pval_adj,n_ident1, n_ident2) %>%
  arrange(-abs(lfc))

print("Differentially expressed genes (assigned vs. unassigned):")
all_res %>% filter(pval_adj < sig_pval_thresh) %>%
  filter(comp == "Assigned vs. Unassigned") %>%
  select(lfc, gene, pval_lm, pval_adj,n_ident1, n_ident2) %>%
  arrange(-abs(lfc))

