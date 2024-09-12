library(Matrix)
library(Seurat)
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggpubr)
#library(grid)

# ID <- "MP-MD-P5-D6_I4" #"MP-MD-P5-D6_I6 #MP-MD-P2-D6_I7
date <- "2022-08-15"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
ExperimentName <- "Pica"
ID <- "PC-P2-D3_I30"
calc_lfcs <- F

section_name <- "3_Load_Data_QC"

# get options for QC
args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  ## needed arguments
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--ID"){ ID <- args[[ind + 1]] }
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
  #if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
  if (arg == "--experiment_name"){ ExperimentName <- args[[ind + 1]] }
}
print(paste("Experiment:", ExperimentName))
print(paste("ID:", ID))
print(paste("HomeFolder:", HomeFolder))
print(paste("Date:", date))

## set colors for markers
vec1 <- c('dCas9-KRAB-MeCP2', 'BSD', 'mScarlet')
col1 <- c('#F8766D')
vec2 <- c('BFP', 'PURO')
col2 <- c('#00BFC4')

# set relevent i/o paths
source(file.path(HomeFolder, "scripts", "io", paste0(ExperimentName, "_io.R")))
source(file.path(HomeFolder, "scripts", "io", paste0("Magpie", "_Utils.R")))
CellRangerOutPath <- paste0(OutFolder, "1_RunCellRanger/",ID,"/")
outdir <- file.path(OutFolder, section_name, "/")
plotsdir <- file.path(HTMLFolder, "pipeline/", section_name, "/")
if(!dir.exists(outdir)) dir.create(outdir)
if(!dir.exists(plotsdir)) dir.create(plotsdir)
#date_new <- format(Sys.Date(), '%Y-%m-%d')

# create Seurat object from CellRanger output
counts <- readMM(paste0(CellRangerOutPath, "filtered_feature_bc_matrix/matrix.mtx.gz"))
features <- read.delim(paste0(CellRangerOutPath, "filtered_feature_bc_matrix/features.tsv.gz"),
                       header = FALSE, stringsAsFactors = FALSE)
features <- apply(features, 1, FUN = function(x){paste(x, collapse = ":")})
features <- gsub(features, pattern = " ", replacement = "-")
rownames(counts) <- features
colnames(counts) <- read.delim(paste0(CellRangerOutPath, "filtered_feature_bc_matrix/barcodes.tsv.gz"),
                               header = FALSE, stringsAsFactors = FALSE)$V1
SeuratObject <- CreateSeuratObject(counts[grepl(":Gene-Expression$", rownames(counts)),], project = ID)
SeuratObject[["guides"]] <- CreateAssayObject(counts = counts[grepl(":CRISPR-Guide-Capture$", rownames(counts)),])
SeuratObject[["percent_MT"]] <- PercentageFeatureSet(SeuratObject, pattern = ":MT-")

# add guide info
protospacers_calls <- read.csv(file.path(CellRangerOutPath, "crispr_analysis/protospacer_calls_per_cell.csv"))
SeuratObject@meta.data <- SeuratObject@meta.data %>%
  rownames_to_column("cell_barcode") %>%
  left_join(protospacers_calls, by = "cell_barcode") %>%
  column_to_rownames("cell_barcode")

# choose cell filters based on local minima
# find first local maxima to determine the minimum in-between
dens_nCount_RNA <- density(SeuratObject$nCount_RNA)
dens_nFeature_RNA <- density(SeuratObject$nFeature_RNA)

nCount_RNA_min_window <- optimize(approxfun(dens_nCount_RNA$x, -dens_nCount_RNA$y),
                                  interval = c(0, max(dens_nCount_RNA$x)/8))$min
nCount_RNA_max_window <- optimize(approxfun(dens_nCount_RNA$x, -dens_nCount_RNA$y),
                                  interval = c(nCount_RNA_min_window + (nCount_RNA_min_window - min(dens_nCount_RNA$x))*0.5, max(dens_nCount_RNA$x)))$min
nFeature_RNA_min_window <- optimize(approxfun(dens_nFeature_RNA$x,-dens_nFeature_RNA$y),
                                    interval = c(0, max(dens_nFeature_RNA$x)/8))$min
nFeature_RNA_max_window <- optimize(approxfun(dens_nFeature_RNA$x,-dens_nFeature_RNA$y),
                                  interval = c(nFeature_RNA_min_window+(nFeature_RNA_min_window - min(dens_nFeature_RNA$x))*0.5, max(dens_nFeature_RNA$x)/2))$min

cutoff_windows = list(nCount_RNA = list(min = nCount_RNA_min_window, max = nCount_RNA_max_window), 
                      nFeature_RNA = list(min = nFeature_RNA_min_window, max = nFeature_RNA_max_window))
#cutoff_windows = list(nCount_RNA = list(min = 2000, max = 10000), 
#                      nFeature_RNA = list(min = 50, max = 2000))

cutoff_features <- c("nCount_RNA", "nFeature_RNA", "mt")

cutoffs <- lapply(cutoff_features, function(cutoff_feature) {
  if(cutoff_feature == "mt"){
    dens <- density(SeuratObject$percent_MT)
    cutoff <- 10
  } else if (cutoff_feature == "nFeature_RNA"){
    dens <- density(SeuratObject$nFeature_RNA)
    cutoff <- 2000
    } else {
    dens <- density(SeuratObject[[cutoff_feature]][,1])
    cutoff <- optimize(approxfun(dens$x,dens$y),
                       interval = c(cutoff_windows[[cutoff_feature]]$min,
                                    cutoff_windows[[cutoff_feature]]$max))$min
  }
  # plot density and cutoffs
  pdf(file.path(plotsdir, paste0(ID,"_", cutoff_feature,".pdf")))
    if (cutoff_feature == "nCount_RNA"){
      plot(dens, main = cutoff_feature) + abline(v = cutoff, col = "red") +
        abline(v = nCount_RNA_min_window, lty = 3) + abline(v = nCount_RNA_max_window, lty = 3)
    } else{
      #else if (cutoff_feature == "nCount_RNA"){
      #plot(dens, main = cutoff_feature) + abline(v = cutoff, col = "red") +
    #    abline(v = nFeature_RNA_min_window, lty = 3) + abline(v = nFeature_RNA_max_window, lty = 3)#} 
      plot(dens, main = cutoff_feature) + abline(v = cutoff, col = "red")
    }
  dev.off()
  print(paste("Using cutoff of", cutoff, "for", cutoff_feature))
  return(cutoff)
})
names(cutoffs) <- cutoff_features
saveRDS(cutoffs, paste0(outdir, "/", date, "_", ID, "_cutoffs.RDS"))

p1 <- qplot(y =SeuratObject$nFeature_RNA, x = SeuratObject$nCount_RNA) +
  xlab("nCount_RNA") + ylab("nFeature_RNA") + 
  ggpointdensity::geom_pointdensity(size=0.3) + scale_color_viridis_c() 
ggsave(p1, filename = paste0(plotsdir, ID, "_nFeature-nCount.pdf"), width = 7, height = 6)

# Filter cells based on cutoffs and normalize data
SeuratObject <- subset(SeuratObject, subset = (nCount_RNA > cutoffs$nCount_RNA & percent_MT < cutoffs$mt & nFeature_RNA >= cutoffs$nFeature_RNA))
SeuratObject <- NormalizeData(SeuratObject)
SeuratObject <- FindVariableFeatures(SeuratObject)
SeuratObject <- ScaleData(SeuratObject)

# Calculate embeddings and clusters per inlet
SeuratObject <- RunPCA(SeuratObject)
SeuratObject <- RunUMAP(SeuratObject, dims = 1:20)
SeuratObject <- FindNeighbors(SeuratObject, dims = 1:20)
SeuratObject <- FindClusters(SeuratObject)

# Add cell cycle info
data("cc.genes.updated.2019")
s.genes <- cc.genes.updated.2019$s.genes
s.genes <- row.names(SeuratObject)[match(s.genes, unlist(lapply(strsplit(row.names(SeuratObject), ":"), "[[", 2)))]
g2m.genes <- cc.genes.updated.2019$g2m.genes
g2m.genes <- row.names(SeuratObject)[match(g2m.genes, unlist(lapply(strsplit(row.names(SeuratObject), ":"), "[[", 2)))]
SeuratObject <- Seurat::CellCycleScoring(SeuratObject,
                                         s.features = s.genes,
                                         g2m.features = g2m.genes)

# get expressions of guide and dCas9 markers
#get_expression_vector <- function(name, type = "normalized"){
#  print(rownames(SeuratObject@assays$RNA@data)[grep(name, rownames(SeuratObject@assays$RNA@data))])
#  if(type == "normalized")
#    x <-  SeuratObject@assays$RNA@data[grep(paste0(name, ":"), rownames(SeuratObject@assays$RNA@data)),]
#  else if(type == "counts")
#    x <- SeuratObject@assays$RNA@counts[grep(paste0(name, ":"), rownames(SeuratObject@assays$RNA@data)),]
#  return(x)
#}

extra.sequences <- gsub(list.files(paste0(ResourcesFolder, "/reference_sequences/Magpie_added_sequences/"),
                                   pattern = "[.]reference[.]fa"),
                        pattern = ".reference.fa", replacement = "")
df <- as.data.frame(mapply(extra.sequences, FUN = function(seq){
  expr <- get_expression_vector(seq, type = "counts", SeuratObject = SeuratObject)
}))
fwrite(df, paste0(outdir, "/", date, "_", ID, "_crispr_marker_counts.txt"), sep = "\t", quote = F, row.names = T)
#for (ID in LaneMetadata$ID){
#  SeuratObject <- readRDS(file.path(outdir,  paste0(date, "_", ID, ".RDS")))
#  df <- as.data.frame(mapply(extra.sequences, FUN = function(seq){
#    expr <- get_expression_vector(seq, type = "counts", SeuratObject = SeuratObject)
#  }))
#  fwrite(df, paste0(outdir, "/", date, "_", ID, "_crispr_marker_counts.txt"), sep = "\t", quote = F, row.names = T)
#}

# plot histogram of marker expression
marker_counts <- mapply(extra.sequences, FUN = function(seq){
  col <- col1
  if (seq %in% vec2){col <- col2}
  p <- ggplot(df, aes(get(seq))) + 
    ggtitle(seq) + xlab("Counts") + ylab("nCells") + 
    geom_histogram(binwidth = 1, col = "black", fill = col)
  return(p)
}, SIMPLIFY = F)
ggarrange(plotlist = marker_counts, ncol = length(extra.sequences))
ggsave(paste0(plotsdir, ID, "_marker-counts.pdf"), width = 3*length(extra.sequences), height = 3)


## plot correlation of markers' expression
df <- as.data.frame(mapply(extra.sequences, FUN = function(seq){
  expr <- get_expression_vector(seq, SeuratObject = SeuratObject)
}))
pdf(paste0(plotsdir, ID, "_marker-expr-correlation.pdf"), width = 8, height = 7)
GGally::ggpairs(df, 
                columns = colnames(df), 
                title = "Correlation of CRISPR Sequence Expressions",
                columnLabels = extra.sequences)
dev.off()
#ggsave(paste0(plotsdir, ID, "_marker-expr-correlation.pdf"), width = 8, height = 7)

## plot correlation of markers' expression for guide and dCas9 separately
#GGally::ggpairs(df, 
#                columns = vec1, 
#                title = "Correlation of dCas9 markers expresssion",
#                columnLabels = vec1)
#ggsave(paste0(plotsdir, ID, "_marker-expr-correlation_dCas9.pdf"), width = 8, height = 7)
#GGally::ggpairs(df, 
#                columns = vec2, 
#                title = "Correlation of guide markers expresssion",
#               columnLabels = vec2)
#ggsave(paste0(plotsdir, ID, "_marker-expr-correlation_guide.pdf"), width = 8, height = 7)

# save Seurat object per inlet
saveRDS(SeuratObject, file = file.path(outdir,  paste0(date, "_", ID, ".RDS")))

sessionInfo()
