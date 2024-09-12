library(Matrix)
library(Seurat)
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(doParallel)

date <- "2022-10-05"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
ExperimentName <- "Pica"
ID <- "PC-P4-D3_I73"
section_name <- "4a_add_guides"
calc_lfcs <- F

# get options for running the script
args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  ## needed arguments
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--ID"){ ID <- args[[ind + 1]] }
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
  if (arg == "--experiment_name"){ ExperimentName <- args[[ind + 1]] }
  if (arg == "--calc_lfcs"){ calc_lfcs <- as.logical(args[[ind + 1]]) }
}
print(paste("ID:", ID))
print(paste("HomeFolder:", HomeFolder))
print(paste("Date:", date))

# set relevent i/o paths
source(file.path(HomeFolder, "scripts", "io", paste0(ExperimentName, "_io.R")))
#source(file.path(HomeFolder, "scripts/io/Magpie_Utils.R"))
outdir <- file.path(OutFolder, section_name, "/")
plotsdir <- file.path(HTMLFolder, "pipeline/", section_name, "/")
if(!dir.exists(outdir)) dir.create(outdir)
if(!dir.exists(plotsdir)) dir.create(plotsdir)
#date_new <- format(Sys.Date(), '%Y-%m-%d')

# load post-QC Seurat object for the inlet
SeuratObject <- readRDS(file.path(OutFolder, "3_Load_Data_QC", paste0(date, "_", ID, ".RDS")))
GuideMetadata <- fread(GuideMetadataPath)
LaneMetadata <- fread(LaneMetadataPath)
n_cells <- ncol(SeuratObject)

## calculate frequency of guides
nUMIS <- colSums(SeuratObject@assays$guides@counts)
max_UMI_across_guides <- apply(SeuratObject@assays$guides@counts, 2, max)
second_UMI_across_guides <- apply(SeuratObject@assays$guides@counts, 2,
                                  FUN=function(x){return(sort(x, decreasing = T)[2])})
freq_max_guide <- max_UMI_across_guides / nUMIS
df <- data.frame(cell_barcode = row.names(SeuratObject@meta.data),
                 n_umi_top_guide = max_UMI_across_guides,
                 nUMIs_total = nUMIS,
                 top_guide = max_UMI_across_guides / nUMIS,
                 second_guide = second_UMI_across_guides/nUMIS)
fwrite(df, sep = "\t", file = file.path(outdir, paste0(date, "_", ID, "_guide_freqs.txt")), row.names = F)
p <- qplot(df$top_guide, df$second_guide) +
  ggpointdensity::geom_pointdensity() +
  xlab("Frequency of most abundant guide in total guide UMIs") +
  ylab("Frequency of second most abundant guide in total guide UMIs")
ggsave(filename = paste0(plotsdir, ID, "_top-v-second-guide.pdf"), plot = p,
       width = 6, height = 5)


# set thresholds for guide assignment baed on freq. of most abundant guide
dens <- density(freq_max_guide, na.rm = T)
cutoff_windows = list(min = 0.5, max = 1)
th_freq <- optimize(approxfun(dens$x,dens$y),interval = c(cutoff_windows$min, cutoff_windows$max))$min
pdf(paste0(plotsdir, ID, "_top-guide-freq.pdf"), width = 7, height = 5)
plot(density(freq_max_guide, na.rm = T),
     main = "Frequency of the Top Guide")
abline(v = th_freq, col = "red", lty =3 )
dev.off()
print(th_freq)
write.table(th_freq, paste0(outdir, date, "_", ID, "_th_freq.txt"), row.names = F, col.names = F, quote = F)

# number of UMIs
th_UMIs <- 2
umi_counts <- data.frame(nUMIs = as.numeric(unlist(strsplit(SeuratObject$num_umis, split = "[|]"))))
ggplot(umi_counts, aes(x = nUMIs)) + 
  scale_x_log10() + scale_y_log10() + 
  ggtitle("nUMIs") + 
  geom_histogram(bins = 50) + 
  geom_vline(xintercept = th_UMIs, col = "red")
ggsave(filename = paste0(plotsdir, ID, "_nUMIs.pdf"), 
       width = 6, height = 5)


### store guide assignments in Seurat object
SeuratObject$CellRanger_num_features <- SeuratObject$num_features
SeuratObject$CellRanger_feature_call <- SeuratObject$feature_call
SeuratObject$CellRanger_num_umis <- SeuratObject$num_umis
SeuratObject$CellRanger_num_features <- SeuratObject$num_features

idx_max_guide <- apply(SeuratObject@assays$guides@counts, 2, which.max)
SeuratObject$guide_call <- sapply(seq_len(n_cells), function(i){
  c <- SeuratObject
  if(freq_max_guide[i] > th_freq & max_UMI_across_guides[i] > th_UMIs) {
    rownames(SeuratObject@assays$guides@counts)[idx_max_guide[i]]
  } else
  {
    "unassigned"
  }
})
SeuratObject$feature_call <- gsub(unlist(lapply(strsplit(SeuratObject$guide_call, ":"), "[[", 1)),
                                  pattern = "-", replacement = "_")
SeuratObject$target_gene_name <- unlist(lapply(SeuratObject$feature_call, FUN = function(x){
  if ( x == "unassigned"){
    return("unassigned")
  } else{
    paste(rev(rev(unlist(strsplit(x, split = "_")))[-1]), collapse = "_")
  }
}), "[[", 1)
#add num umis
SeuratObject$num_features <- mapply(SeuratObject$guide_call, FUN = function(x){
  if (x == "unassigned"){ return(0) }
  else {
    return(length(unlist(strsplit(x, "[|]"))))
  }
})
guide_ind <- match(SeuratObject$guide_call, row.names(SeuratObject$guides))
SeuratObject$num_umis <- as.vector(SeuratObject$guides@counts)[ 0:(ncol(SeuratObject$guides@counts)-1)*nrow(SeuratObject$guides@counts) + guide_ind ]


## check that the assigned guides are actually in the right group
inlet_guide_strength <- LaneMetadata$Guide_Strength[match(ID, LaneMetadata$ID)]
observed_guide_strength <- unique(GuideMetadata$group[match(SeuratObject$feature_call, GuideMetadata$guide_name)])

#freak out if the guide group doesn't match
#means that cellranger was run with the wrong guide library list
stopifnot(all(sort(observed_guide_strength) %in% sort(c(inlet_guide_strength, "Control")) ))
if ( ! all(sort(c(inlet_guide_strength, "Control")) == sort(observed_guide_strength) ) ){
  print('no control cells!')
}

#number of assigned cells
sort(table(SeuratObject$guide_call), decreasing = T)[1:10]
paste0("Number of Perturbed Cells: ", length(which(!(SeuratObject$guide_call %in% c("unassigned", NonTargetGeneName)))))
paste0("Number of ", NonTargetGeneName, " Cells: ", length(which(SeuratObject$target_gene_name == NonTargetGeneName)))
paste0("Number of Unassigned Cells: ", length(which(SeuratObject$target_gene_name == "unassigned")))
df4plot <- data.frame(
  status = c('Perturbed', 'Non-Target', 'Unassigned'),
  num_cells = c( length(which(!(SeuratObject$guide_call %in% c("unassigned", NonTargetGeneName)))),
                 length(which(SeuratObject$target_gene_name == NonTargetGeneName)),
                 length(which(SeuratObject$target_gene_name == "unassigned")) 
  )
)
pdf(paste0(plotsdir, ID, "_num_cells.pdf"), width = 6, height = 5)
ggplot(df4plot, aes(x = status, y = num_cells, fill = status)) + 
  ggtitle(ID) + xlab('') + ylab("Number of Cells") + 
  geom_bar(stat = 'identity') + theme_bw()
dev.off()

#caluclate lfcs
if (calc_lfcs == T){
  print("Calculating LFCs...")
  cells_per_gene <- table(SeuratObject$target_gene_name)
  genes2consider <- names(cells_per_gene[which(cells_per_gene > min_cells_per_gene)])
  highly_expressed_genes <- rowMeans(SeuratObject$RNA@data)
  highly_expressed_genes <- names(highly_expressed_genes)[which(highly_expressed_genes > min_expr_thresh)]
  highly_expressed_genes <- unlist(lapply(strsplit(highly_expressed_genes, ":"), "[[", 2))
  genes2consider <- genes2consider[which(genes2consider %in% highly_expressed_genes)]
  
  control_cells <- row.names(SeuratObject@meta.data)[SeuratObject$target_gene_name == NonTargetGeneName]
  is_control_cell = SeuratObject$target_gene_name == NonTargetGeneName
  registerDoParallel(20)
  target_lfc <- foreach(tg = genes2consider) %dopar% {
    
    gene_ind <- which(unlist(lapply(strsplit(row.names(SeuratObject$RNA@data), ":"), "[[", 2)) == tg)
    perturbed_cells <- row.names(SeuratObject@meta.data)[SeuratObject$target_gene_name == tg]
    is_perturbed_cell <- SeuratObject$target_gene_name == tg
    df4lm <- data.frame(
      expression = c(SeuratObject$RNA@data[gene_ind, is_perturbed_cell], SeuratObject$RNA@data[gene_ind, is_control_cell]),
      status = c(rep("Perturbed", length(which(is_perturbed_cell))), rep("Control", length(which(is_control_cell)))),
      cell_cycle = c(SeuratObject$Phase[is_perturbed_cell], SeuratObject$Phase[is_control_cell]),
      nCount_RNA = c(SeuratObject$nCount_RNA[is_perturbed_cell], SeuratObject$nCount_RNA[is_control_cell]),
      percent_MT = c(SeuratObject$percent_MT[is_perturbed_cell], SeuratObject$percent_MT[is_control_cell])
    )
    fit <- lm(expression ~ status + cell_cycle + nCount_RNA + percent_MT,
              data = df4lm)
    rtn <- list(target_lfc = coefficients(fit)['statusPerturbed'],
                pval = summary(fit)$coefficients['statusPerturbed', 4])
    return(rtn)
  }
  names(target_lfc ) <- genes2consider
  target_lfc  <- as.data.frame(bind_rows(target_lfc ))
  target_lfc$target <- genes2consider
  target_lfc$pval_adj <- p.adjust(target_lfc$pval, method = "BH")
  pdf(paste0(plotsdir, ID, "_target_lfc.pdf"), width = 6, height = 5)
  p <- ggplot(target_lfc, aes(x = target_lfc, y = -log10(pval_adj), 
                         label = ifelse( pval_adj < 0.05, target, ""))) + 
    ggrepel::geom_label_repel(max.overlaps = 20) + 
    xlab("Target LFC (Natural Log)") + ylab("-Log10 (P-value)") + ggtitle(ID) + 
    geom_point() + theme_bw()
  print(p)
  dev.off()
}


# save object
saveRDS(SeuratObject, file = file.path(outdir, paste0(ID, "_guides_assigned.RDS")))

sessionInfo()
