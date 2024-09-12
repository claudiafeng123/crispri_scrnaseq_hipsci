library(Matrix)
library(data.table)
library(dplyr)
library(tidyverse)
library(Seurat)
library(parallel)
# print(paste("Number of cores:", parallel::detectCores()))


date <- "2022-08-15"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
guide_group <- "Strong" #"Strong", "All_Inlets"
ExperimentName <- "Magpie"
section_name <- "5e_compare_timepoints"
cell_type <- "unperturbed" #perturbed, all



# options
args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--guide_group"){ guide_group <- args[[ind + 1]] }
  if (arg == "--experiment_name"){ ExperimentName <- args[[ind + 1]] }
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
  if (arg == "--cell_type"){ cell_type <- args[[ind + 1]] }
  if (arg == "--ncores"){ ncores <- args[[ind + 1]] } 
  if (arg == "--REDO"){ REDO <- as.logical(args[[ind + 1]]) } 
}

# i/o
source(file.path(HomeFolder, "scripts/io/Magpie_io.R"))
outdir <- file.path(OutFolder, section_name)
plotsdir <- file.path(HTMLFolder, "pipeline/", section_name)
if(!dir.exists(outdir)) dir.create(outdir)
if(!dir.exists(plotsdir)) dir.create(plotsdir)

fnm <- file.path(outdir, paste0(date, "_", guide_group, "_", cell_type, "_cells_comparison_timepoints.csv"))
print(fnm)

if (REDO == F & file.exists(fnm)){
  all_res <- fread(fnm)
} else {
  # load data
  if(guide_group == "Strong") {
    inlets_combined <- readRDS(file.path(OutFolder, "5a_combine",
                                         paste0(date, "_Strong_combined.RDS")))
  } else if (guide_group == "Moderate") {
    inlets_combined <- readRDS(file.path(OutFolder, "5a_combine",
                                         paste0(date, "_Moderate_combined.RDS")))
    
  } else {
    stop("Invalid guide group specified.")
  }
  
  if (cell_type == "unperturbed"){
    inlets_combined <- subset(inlets_combined, target_gene_name %in% c(NonTargetGeneName, "unassigned"))
  } else if (cell_type == "perturbed"){
    inlets_combined <- subset(inlets_combined, target_gene_name != NonTargetGeneName & target_gene_name != "unassigned")
  }
  
  # timepoints
  time_points <- unique(inlets_combined$Timepoint)
  print(table(inlets_combined$Timepoint))
  
  # test all genes (min. expression)
  genes_to_test <- rownames(inlets_combined)
  print(paste("Number of genes tested:", length(genes_to_test)))
  
  stopifnot(length(time_points) > 1)
  comps <- combn(time_points, 2)
  
  # test for differences
  all_res <- lapply(seq_len(ncol(comps)), function(compidx) {
    t1 <- comps[1,compidx]
    t2 <- comps[2,compidx]
    print(paste("Testing", t1, "vs.", t2))
    res <- mclapply(genes_to_test, function(gene2test){
      print(gene2test)
      #y is the expression vector for the target gene in question
      y <- inlets_combined@assays$RNA@data[gene2test,]
      ident.1 <- inlets_combined$Timepoint == t1
      ident.2 <- inlets_combined$Timepoint == t2
      data4lm <- data.frame(expr = c(y[ident.1],y[ident.2]),
                            time_point = c(rep(0, sum(ident.1)), rep(1, sum(ident.2))),
                            cell_line = c(inlets_combined$cell_line[ident.1], inlets_combined$cell_line[ident.2]),
                            s.score = c(inlets_combined$S.Score[ident.1], inlets_combined$S.Score[ident.2]),
                            g2m.score = c(inlets_combined$G2M.Score[ident.1], inlets_combined$G2M.Score[ident.2]),
                            # inlet = c(inlets_combined$orig.ident[ident.1], inlets_combined$orig.ident[ident.2]),
                            nCount_RNA = c(inlets_combined$nCount_RNA[ident.1], inlets_combined$nCount_RNA[ident.2]),
                            percent_MT = c(inlets_combined$percent_MT[ident.1], inlets_combined$percent_MT[ident.2])
      )
      fit <- lm(expr ~ cell_line + s.score + g2m.score + time_point + nCount_RNA + percent_MT, data = data4lm)
      beta <- coef(fit)["time_point"]
      pval <- summary(fit)$coefficients["time_point",4]
      data.frame(lfc = beta,
                 pval_lm = pval,
                 n_ident1 = sum(ident.1),
                 n_ident2 = sum(ident.2),
                 gene = gene2test,
                 mean_expr = mean(y),
                 inlets = length(unique(inlets_combined$orig.ident)),
                 ident2_name = t2,
                 ident1_name = t1)
    }, mc.cores = ncores) 
    res <- res %>% bind_rows()
    res$pval_adj <- p.adjust(res$pval_lm, method = "BH")
    res$comp <- paste0(t1,"_vs_", t2)
    res
  }) %>% bind_rows()
  
  # save results
  write.csv(all_res, fnm)
  
}

# make plot
pdf(file.path(plotsdir,
              paste0(date, "_", guide_group, "_", cell_type, "_cells_",
                     "_comparison_timepoints.pdf")), width = 12, height = 5)
all_res$comp <- gsub(gsub(all_res$comp, pattern = "_", replacement = " "), pattern = "vs", replacement = "vs.")
print(ggplot(all_res,  aes(x = lfc, y= -log10(pval_lm), col = mean_expr)) +
        ggtitle(paste0(cell_type, ", ", guide_group)) + xlab("LFC (Natural Log)") + ylab("-log10(p-value)") + 
        facet_wrap(~comp) + geom_point() + geom_hline(yintercept = -log10(0.05), col = "red", lty = 3) + theme_bw() )
dev.off()

# print out
all_res %>% filter(pval_adj < 0.1) %>%
  select(lfc, gene, pval_lm, pval_adj, n_ident1, n_ident2, ident1_name, ident2_name) %>%
  arrange(-abs(lfc))


## ---- Number of Cells
# do strong guides do well on d5?

#pdf(paste0(plotsdir, "/", date, "_target_gene_lfc_by_day_and_strength.pdf"), width = 8, height = 8)
#df4plot <- target_gene_lfc %>% dplyr::filter(guide_strength == "Strong")
#df4plot$day <- GuideMetadata$subgroup[match(df4plot$target, GuideMetadata$gene)]
#p <- ggplot(df4plot, aes(x = lfc, col = day)) + 
#  facet_wrap(~timepoint, ncol = 1) + 
#  xlab("LFC (Base 10)") + 
#  geom_density() + theme_bw()
#print(p)
#dev.off()

if (guide_group == "Strong" & cell_type == "all"){
  print("check number of cells for strong inlets")
  ## fraction of cells
  pdf(paste0(plotsdir, "/", date, "_cell_strong_timepoint_composition.pdf"), width = 5, height = 5)
  guide_strength <- "Strong"
  
  GuideMetadata <- fread(GuideMetadataPath)
  LaneMetadata <- fread(LaneMetadataPath)
  timepoints <- LaneMetadata %>% filter(Guide_Strength == guide_strength) %>% .$Timepoint %>% unique()
  
  perturbed_cell_metadata <- fread(paste0(OutFolder, "5a_combine", "/", date, "_", guide_strength, "_perturbed_cells_combined_meta.csv"))
  unassigned_cell_metadata <- fread(paste0(OutFolder, "5a_combine", "/", date, "_", guide_strength, "_unassigned_cells_combined_meta.csv"))
  nontarget_cell_metadata <- fread(paste0(OutFolder, "5a_combine", "/", date, "_", guide_strength, "_nontarget_cells_combined_meta.csv"))
  
  #num_unassigned <- unique(x$n_control)
  assigned_cell_counts <- as.data.frame(table(perturbed_cell_metadata %>% select(c("target_gene_name", "Timepoint"))))
  assigned_cell_counts <- reshape2::dcast(assigned_cell_counts, target_gene_name ~ Timepoint, value.var = "Freq")
  assigned_cell_counts$subgroup <- paste0("Group_", GuideMetadata$subgroup[match(gsub(assigned_cell_counts$target_gene_name, pattern = "_", replacement = "-"), GuideMetadata$gene)])
  assigned_cell_counts.split <- split.data.frame(assigned_cell_counts, f = assigned_cell_counts$subgroup)
  cell_counts_by_group <- as.data.frame(t(mapply(assigned_cell_counts.split, FUN = function(x){
    apply(x %>% select(contains("D")), 2, sum)
  })))
  unassigned_cells <- table(unassigned_cell_metadata$Timepoint)
  nontarget_cells <- table(nontarget_cell_metadata$Timepoint)
  cell_counts_by_group <- rbind(cell_counts_by_group, unassigned_cells, nontarget_cells)
  row.names(cell_counts_by_group)[c(dim(cell_counts_by_group)[1] - 1, dim(cell_counts_by_group)[1])] <- c("Unassigned", "Non-Target")
  cell_counts_by_group$guide_strength <- row.names(cell_counts_by_group)
  
  df4plot <- reshape2::melt(cell_counts_by_group, id = "guide_strength")
  names(df4plot) <- c("Guide_Strength", "Timepoint", "Num_Cells")
  p <- ggplot(df4plot, aes(x = Timepoint, y =Num_Cells, fill = Guide_Strength)) + 
    scale_fill_manual(values = c("Unassigned" = 'lightgray', 
                                 "Control" = "darkgray",
                                 "Group_D3" = "cornflowerblue",
                                 "Group_D4" = "lightblue",
                                 "Group_D5" = "darkblue")) + 
    ggtitle(guide_strength) + ylab("Fraction of Cells") + 
    geom_bar(stat = 'identity') + theme_bw()
  print(p)
  
  p <- ggplot(df4plot %>% filter(Guide_Strength != "Unassigned"), aes(x = Timepoint, y =Num_Cells, fill = Guide_Strength)) + 
    scale_fill_manual(values = c("Unassigned" = 'lightgray', 
                                 "Control" = "darkgray",
                                 "Group_D3" = "cornflowerblue",
                                 "Group_D4" = "lightblue",
                                 "Group_D5" = "darkblue")) + 
    ggtitle(guide_strength) + ylab("Fraction of Cells") + 
    geom_bar(stat = 'identity') + theme_bw()
  print(p)
  
  
  cell_counts_by_group[, timepoints] <- sweep(cell_counts_by_group[, timepoints],2,colSums(cell_counts_by_group[, timepoints]),`/`)
  df4plot <- reshape2::melt(cell_counts_by_group, id = "guide_strength")
  names(df4plot) <- c("Guide_Strength", "Timepoint", "Frac_Cells")
  
  
  
  p <- ggplot(df4plot, aes(x = Timepoint, y =Frac_Cells, fill = Guide_Strength)) + 
    scale_fill_manual(values = c("Unassigned" = 'lightgray', 
                                 "Control" = "darkgray",
                                 "Group_D3" = "cornflowerblue",
                                 "Group_D4" = "lightblue",
                                 "Group_D5" = "darkblue")) + 
    ggtitle(guide_strength) + ylab("Fraction of Cells") + 
    geom_bar(stat = 'identity') + theme_bw()
  print(p)
  
  cell_counts_by_group <- cell_counts_by_group %>% filter(guide_strength != "Unassigned")
  cell_counts_by_group[, timepoints] <- sweep(cell_counts_by_group[, timepoints],2,colSums(cell_counts_by_group[, timepoints]),`/`)
  df4plot <- reshape2::melt(cell_counts_by_group, id = "guide_strength")
  names(df4plot) <- c("Guide_Strength", "Timepoint", "Frac_Cells")
  
  p <- ggplot(df4plot, aes(x = Timepoint, y =Frac_Cells, fill = Guide_Strength)) + 
    scale_fill_manual(values = c("Control" = "darkgray",
                                 "Group_D3" = "cornflowerblue",
                                 "Group_D4" = "lightblue",
                                 "Group_D5" = "darkblue")) + 
    ggtitle(guide_strength) + ylab("Fraction of Assigned Cells") + 
    geom_bar(stat = 'identity') + theme_bw()
  print(p)
  
  dev.off()
  
}


## Target gene

