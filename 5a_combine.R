library(Matrix)
library(Seurat)
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(R.utils)
library(ggpubr)
library(SeuratObject)

HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
date <- "2022-10-05"
min_expr_thresh <- 0.1
section_name <- "5a_combine"
ExperimentName <- "Pica"
#date <- "2022-08-15"
#ExperimentName <- "Magpie"
#final_object <- "Pica"

#separate <- T
include_target_expr <- T

args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  ##needed arguments
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--final_object"){ final_object <- args[[ind + 1]] }
  if (arg == "--experiment_name"){ ExperimentName <- args[[ind + 1]] }
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
  if (arg == "--min_expr_thresh"){ min_expr_thresh <- as.numeric(args[[ind + 1]]) }
  if (arg == "--include_target_expr"){ include_target_expr <- as.logical(args[[ind + 1]]) }
}


# set relevent i/o paths
source(file.path(HomeFolder, "scripts/io/", paste0(ExperimentName, "_io.R")))
source(paste0(HomeFolder, "scripts/io/Magpie_Utils.R"))
outdir <- file.path(OutFolder, section_name)
plotsdir <- file.path(HTMLFolder, "pipeline/", section_name)
if(!dir.exists(outdir)) dir.create(outdir)
if(!dir.exists(plotsdir)) dir.create(plotsdir)
date_new <- date#format(Sys.Date(), '%Y-%m-%d')

LaneMetadata <- fread(LaneMetadataPath)
GuideMetadata <- fread(GuideMetadataPath)
IDs <- LaneMetadata$ID
control_tag <- c("unassigned", NonTargetGeneName)
target_genes <- unique(GuideMetadata$gene)

# read all data
files <- list.files(paste0(OutFolder, "4b_add_donors/"))
files <- files[grepl("_with_donor.RDS", files) & !grepl("combined", files)]
files <- files[gsub(files, pattern = paste0(date, "_|_with_donor.RDS"), replacement = "") %in% IDs]
#use only IDs where data already exists
IDs_available <- gsub(files, pattern = paste0(date, "_|_with_donor.RDS"), replacement = "")

print("IDs without data:")
print(IDs[which(!(IDs %in% IDs_available))])
print("IDs going into combined object:")
IDs <- IDs[which(IDs %in% IDs_available)]
print(paste(length(IDs), "inlets total"))

#load data
print("Loading data...")
inlet_list <- lapply(files, function(fnm){
    single_inlet <- readRDS(file.path(OutFolder, "4b_add_donors", fnm))
    ID <- gsub(fnm, pattern = paste0(date, "_|_with_donor.RDS"), replacement = "")
    
    mean_expr_per_gene <- rowMeans(single_inlet[["RNA"]]@data)
    
    # lowly-expressed genes
    target_gene_ind <- which(unlist(lapply(strsplit(names(mean_expr_per_gene), ":"), "[[", 2)) %in% target_genes)
    if (fnm == files[1]){
      ## whcih target genes aren't in the expression matrix
      in_expr_mat <- unlist(lapply(strsplit(names(mean_expr_per_gene), ":"), "[[", 2))[target_gene_ind]
      not_in_expr_mat <- target_genes[which(!(target_genes %in% in_expr_mat))]
      print('targets without expression:')
      print(not_in_expr_mat)
      
      #plot mean expression
      pdf(file.path(plotsdir, paste0(date_new, "_", ID, "_reads_per_gene.pdf")), width = 7, height = 6)
      hist(mean_expr_per_gene, breaks = 50,
           xlab = "Log-Normalized Mean Expression", main = "Histogram of Mean Expression")
      abline(v = min_expr_thresh, col = "red", lty = 2, lwd = 3)
      dev.off()
    }
    gene_inds <- sort(unique(c(which(mean_expr_per_gene > min_expr_thresh), 
                                     target_gene_ind)))
    gene_list <- names(mean_expr_per_gene)[gene_inds]
                   
    
    return(list(seurat_object = single_inlet, gene_list = gene_list))
  })
names(inlet_list) <- IDs_available

## saving crispr-sequence counts
crispr_sequences <- gsub(list.files(file.path(ResourcesFolder, "reference_sequences", "Magpie_added_sequences"), pattern = "[.]gtf"),
                         pattern = "[.]gtf", replacement = "")
seq_names <- paste(crispr_sequences, crispr_sequences, "Gene-Expression", sep = ":")
crispr_counts <- lapply(inlet_list, FUN = function(i){
  rtn <- i$seurat_object[["RNA"]]@counts[seq_names,]
  rtn <- as.data.frame(t(rtn))
  guide_strength <- ifelse(substr(i$seurat_object@meta.data$orig.ident, 4, 5) == "MD", "Moderate", "Strong")
  row.names(rtn) <- paste0(guide_strength, "_", as.character(i$seurat_object@meta.data$orig.ident), "_", row.names(rtn))
  return(rtn)
}) %>% bind_rows() %>% as.data.frame() %>%
  rownames_to_column("id")
colnames(crispr_counts) <- unlist(lapply(strsplit(colnames(crispr_counts), ":"), "[[", 1))
fwrite(crispr_counts, paste0(outdir, "/", date, "_crispr_sequence_counts.tsv"), sep = "\t")





# list of highly expressed genes, use intersection of genes expressed above threshold in all inlets
highly_expressed_genes <- c()
#gene_list <- 
for (x in inlet_list){
  if (length(highly_expressed_genes) == 0){
    highly_expressed_genes <- x$gene_list
    print(paste("Start with ID:", x$seurat_object$orig.ident[1]))
    print(paste(length(highly_expressed_genes), "genes"))
  } else {
    prev_num_genes <- length(highly_expressed_genes)
    highly_expressed_genes <- intersect(highly_expressed_genes, x$gene_list)
    print(paste("Add ID:", x$seurat_object$orig.ident[1]))
    print(paste(length(highly_expressed_genes), "genes (removed", prev_num_genes - length(highly_expressed_genes), "genes)"))
  }
}
saveRDS(highly_expressed_genes, paste0(outdir, "/", date, "_gene_list.RDS"))



#crispr_sequence_counts <- inlets_combined[["RNA"]]@counts[ which(unlist(lapply(strsplit(row.names(inlets_combined[["RNA"]]@counts), ":"), "[[", 2)) %in% crispr_sequences),]



# remove lowly expressed genes
print("Removing lowly expressed genes...")
inlet_list <- mapply(inlet_list, FUN = function(x){
  seurat_object <- x$seurat_object
  seurat_object[["RNA"]]@data <- seurat_object[["RNA"]]@data[highly_expressed_genes, ]
  seurat_object[["RNA"]]@counts <- seurat_object[["RNA"]]@counts[highly_expressed_genes, ]
  return(seurat_object)
})
names(inlet_list) <- IDs_available


## combine all inlets (strong, moderate or all)
if (final_object %in% c("Strong", "Moderate")){
  relevant_IDs <- LaneMetadata$ID[which(LaneMetadata$Guide_Strength == final_object)]
} else {
  relevant_IDs <- LaneMetadata$ID
}
relevant_IDs <- intersect(relevant_IDs, IDs_available)
relevant_inlets <- inlet_list[relevant_IDs]

## combine assigned and unassigned cells for strong and moderate
if (final_object %in% c("Strong", "Moderate")){
  inlets_combined <- merge(x = relevant_inlets[[1]],
                           y = relevant_inlets[-1],
                           add.cell.ids = names(relevant_inlets),
                           project = paste(ExperimentName, "combined", sep = "_"),
                           merge.data = TRUE)
  inlets_combined <- FindVariableFeatures(inlets_combined)
  inlets_combined <- ScaleData(inlets_combined)
  inlets_combined <- RunPCA(inlets_combined)
  pdf(file.path(plotsdir, paste0(final_object, "_ELBO_combined.pdf")))
  print(ElbowPlot(inlets_combined, ndims = 50) + geom_hline(yintercept = 2) + ggtitle(final_object))
  dev.off()
  inlets_combined <- RunUMAP(inlets_combined, dims = 1:20)
  #nms <- rownames(inlets_combined@meta.data)
  nms <- rownames(inlets_combined@meta.data)
  inlets_combined@meta.data <- inlets_combined@meta.data %>% select(c("orig.ident", "nCount_RNA", "nFeature_RNA", "nCount_guides", "nFeature_guides", "percent_MT", "num_umis",
                                                                      "feature_call", "guide_call", "target_gene_name", "num_features",
                                                                      "donor", "cell_line", "hipsci_line_name",
                                                                      "S.Score", "G2M.Score", "Phase"
                                                                      ))
  inlets_combined@meta.data <- left_join(inlets_combined@meta.data,
                                         select(LaneMetadata, ID, Guide_Group, Guide_Strength, Pool, Timepoint, Batch),
                                         by = c("orig.ident" = "ID"))
  rownames(inlets_combined@meta.data) <- nms
  
  print("seurat object combined successfully!")
  print(inlets_combined)
  saveRDS(inlets_combined, file = file.path(OutFolder, section_name, paste0(date_new, "_", final_object, "_combined.RDS")))
  write.csv(inlets_combined@meta.data,
            file = file.path(OutFolder, section_name, paste0(date_new, "_", final_object, "_combined_meta.csv")))
  print("seurat object written")
  
  ## create separate objects for assigned and unassigned cells for strong, moderate and combined
  print("perturbed cells")
  perturbed_cells <- subset(inlets_combined, target_gene_name != "unassigned" & target_gene_name != NonTargetGeneName)
  print("Number of perturbed cells:")
  print(dim(perturbed_cells@meta.data)[1])
  
  print("unassigned cells")
  unassigned_cells <- subset(inlets_combined, target_gene_name == "unassigned")
  print("Number of unassigned cells:")
  print(dim(unassigned_cells@meta.data)[1])
  
  print("nontargeting cells")
  nontarget_cells <- subset(inlets_combined, target_gene_name == NonTargetGeneName)
  print("Number of non-target cells:")
  print(dim(nontarget_cells@meta.data)[1])
  
  for (seurat_object in c("perturbed_cells", "unassigned_cells", "nontarget_cells")){
    print(seurat_object)
    inlets_combined <- get(seurat_object)
    if (seurat_object != "nontarget_cells"){
      inlets_combined <- FindVariableFeatures(inlets_combined)
      inlets_combined <- ScaleData(inlets_combined)
      inlets_combined <- RunPCA(inlets_combined)
      pdf(file.path(plotsdir, paste0(final_object, "_ELBO_combined.pdf")))
      print(ElbowPlot(inlets_combined, ndims = 50) + geom_hline(yintercept = 2) + ggtitle(final_object))
      dev.off()
      inlets_combined <- RunUMAP(inlets_combined, dims = 1:20)
      print("seurat object combined successfully!")
      print(inlets_combined)
    } 
    saveRDS(inlets_combined,
            file = file.path(OutFolder, section_name, paste0(date_new, "_", final_object, "_", seurat_object, "_combined.RDS")))
    write.csv(inlets_combined@meta.data,
              file = file.path(OutFolder, section_name, paste0(date_new, "_", final_object, "_", seurat_object, "_combined_meta.csv")))
    print(paste0(seurat_object, " object saved!"))
  }
  
  
} else if (ExperimentName == "Magpie"){
  #all inlets
  #will need to separate
  
  for (seurat_object in c("perturbed_cells", "unassigned_cells", "nontarget_cells")){
    print(paste0("merging ", seurat_object))
    
    data2merge <- mapply(relevant_inlets, FUN = function(x){
      if (seurat_object == "unassigned_cells"){
        x <- subset(x, target_gene_name == 'unassigned')
      } else {
        x <- subset(x, target_gene_name != 'unassigned')
        if (seurat_object == "nontarget_cells"){
          x <- subset(x, target_gene_name == NonTargetGeneName)
        } else if (seurat_object == "perturbed_cells"){
          x <- subset(x, target_gene_name != NonTargetGeneName)
        }
      }
      return(x)
    })
    inlets_combined <- merge(x = data2merge[[1]],
                             y = data2merge[-1],
                             add.cell.ids = names(relevant_inlets),
                             project = paste(ExperimentName, "combined", sep = "_"),
                             merge.data = TRUE)
    inlets_combined <- FindVariableFeatures(inlets_combined)
    inlets_combined <- ScaleData(inlets_combined)
    inlets_combined <- RunPCA(inlets_combined)
    pdf(file.path(plotsdir, paste0(final_object, "_ELBO_combined.pdf")))
    print(ElbowPlot(inlets_combined, ndims = 50) + geom_hline(yintercept = 2) + ggtitle(final_object))
    dev.off()
    inlets_combined <- RunUMAP(inlets_combined, dims = 1:20)
    #nms <- rownames(inlets_combined@meta.data)
    nms <- rownames(inlets_combined@meta.data)
    
    print("seurat object combined successfully!")
    print(inlets_combined)
    saveRDS(inlets_combined,
            file = file.path(OutFolder, section_name, paste0(date_new, "_", final_object, "_", seurat_object, "_combined.RDS")))
    print("saved to ")
    print(file.path(OutFolder, section_name, paste0(date_new, "_", final_object, "_", seurat_object, "_combined.RDS")))
    write.csv(inlets_combined@meta.data,
              file = file.path(OutFolder, section_name, paste0(date_new, "_", final_object, "_", seurat_object, "_combined_meta.csv")))
    print(paste0(seurat_object, " object saved!"))
  }
  
} else if (ExperimentName == "Pica"){
  
  #split by line
  inlet_list <- mapply(names(inlet_list), FUN = function(x){
    seurat_object <- inlet_list[[x]]
    seurat_object$inlet <- x
    SplitObject(seurat_object, split.by = "cell_line")
  })

  cell_lines <- sort(gsub(unique(unlist(strsplit(LaneMetadata$Donors, ";"))), pattern = "HPS.*-", replacement = ""))
  control_seurats <- mapply(cell_line =  cell_lines, FUN = function(cell_line){
    print(cell_line)
    data2merge <- mapply(inlet_list, FUN = function(x){
      #print(x)
      if (length(which(names(x) %in% cell_line)) == 1){
        if (dim(x[[cell_line]]@meta.data)[1] > 1){
          return(x[[cell_line]])
        } else {
          return(NULL)
        }
        
      #  x <- x[[cell_line]]
       # row.names(x@meta.data) = paste0(x$inlet, "_", row.names(x@meta.data))
       #return(x)
      } else {
        return(NULL)
      }
    }, SIMPLIFY = F)
    
    data2merge[sapply(data2merge, is.null)] <- NULL
    
    if (length(data2merge) > 1){
       data2merge <- merge(x = data2merge[[1]],
                        y = data2merge[-1],
                        add.cell.ids = names(data2merge),
                        project = paste(ExperimentName, "combined", sep = "_"),
                        merge.data = TRUE)
    data2merge <- FindVariableFeatures(data2merge)
    data2merge <- ScaleData(data2merge)
    data2merge <- RunPCA(data2merge)
    
    perturbed_cells <- subset(data2merge, target_gene_name != NonTargetGeneName & target_gene_name != "unassigned")
    fwrite(perturbed_cells@meta.data, quote = F, row.names = T, sep = "\t",
           file.path(outdir, paste0(date, "_", cell_line, "_perturbed_cells_metadata.tsv")))
    saveRDS(perturbed_cells, file.path(outdir, paste0(date, "_", cell_line, "_perturbed_cells.RDS")))
    control_cells <- subset(data2merge, target_gene_name == NonTargetGeneName )
    fwrite(control_cells@meta.data, quote = F, row.names = T, sep = "\t",
           file.path(outdir, paste0(date, "_", cell_line, "_control_cells_metadata.tsv")))
    saveRDS(control_cells, file.path(outdir, paste0(date, "_", cell_line, "_control_cells.RDS")))
    unassigned_cells <- subset(data2merge, target_gene_name == "unassigned")
    fwrite(unassigned_cells@meta.data, quote = F, row.names = T, sep = "\t",
           file.path(outdir, paste0(date, "_", cell_line, "_unassigned_cells_metadata.tsv")))
    saveRDS(unassigned_cells, file.path(outdir, paste0(date, "_", cell_line, "_unassigned_cells.RDS")))
    return(control_cells)
    } else{
      return(NULL)
    }
   
  }, SIMPLIFY = F)

  #names(control_seurats) <- cell_lines
  #control_seurats[sapply(control_seurats, is.null)] <- NULL

  control_fnms <- grep(grep(list.files(file.path(OutFolder, section_name), pattern = "control_cells.RDS"), 
    pattern = date, value = T), pattern = "ALL", invert = T, value = T)
  cell_lines <- gsub(control_fnms, pattern = paste0(date, "_|_control_cells.RDS"), replacement = "")

  control_seurats <- mapply(control_fnms, FUN = function(x){
    readRDS(file.path(OutFolder, "5a_combine", x))
    }, SIMPLIFY = F)
  names(control_seurats) <- cell_lines
  
  #merge control cells
  control_seurats <- merge(x = control_seurats[[1]],
                           y = control_seurats[-1],
                           add.cell.ids = names(cell_lines),
                           project = paste(ExperimentName, "combined", sep = "_"),
                           merge.data = TRUE)
  control_seurats <- FindVariableFeatures(control_seurats)
  control_seurats <- ScaleData(control_seurats)
  control_seurats<- RunPCA(control_seurats)
  saveRDS(control_seurats, file.path(outdir, paste0(date, "_ALL_control_cells.RDS")))
  
}




#number of genes
if (final_object == "All Inlets"){
  write.table(length(highly_expressed_genes), file = file.path(outdir, paste0(date_new, "_num_highly_expressed_genes.txt")), 
              sep = "\t", quote = F, row.names = F, col.names = F)
}

sessionInfo()
      