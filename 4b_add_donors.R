library(data.table)
library(Seurat)
library(dplyr)

ID <- "MP-ST-P4-D5_I20"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
date <- "2021-11-12"


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
}

print(paste("ID:", ID))
print(paste("HomeFolder:", HomeFolder))
print(paste("date:", date))

setwd(HomeFolder)
# set relevent i/o paths
source(paste0("scripts/io/", ExperimentName, "_io.R"))
section_name <- "4b_add_donors"
outdir <- file.path(OutFolder, section_name, "/")
plotsdir <- file.path(HTMLFolder, "pipeline/", section_name, "/")
if(!dir.exists(outdir)) dir.create(outdir)
if(!dir.exists(plotsdir)) dir.create(plotsdir)


if (!file.exists(paste0(OutFolder, "2_Genotyping/", ID, "/donor_ids.tsv"))){
  print('No vireo output!!!!!')
}

donor_ids <- fread(paste0(OutFolder, "2_Genotyping/", ID, "/donor_ids.tsv"))
row.names(donor_ids) <- donor_ids$cell
donor_ids$cell_line <- gsub(donor_ids$donor_id, pattern = ".*-", replacement = "")
donor_ids$hipsci_line_name <- donor_ids$donor_id
donor_ids$donor <- unlist(mapply(donor_ids$donor_id, FUN = function(x){
  x <- gsub(x, pattern = ".*-", replacement = "")
  x <- gsub(x, pattern = "_.*", replacement = "")
  return(x)
}))

SeuratObject <- readRDS(paste0(OutFolder, "4a_add_guides/", ID, "_guides_assigned.RDS"))

SeuratObject <- AddMetaData(SeuratObject, metadata = donor_ids %>% dplyr::select(c("donor", "cell_line", "hipsci_line_name")), col.name = NULL)
SeuratObject <- subset(SeuratObject, donor != "doublet" & donor != "unassigned")

#UMAP

# cells per donor
cells_per_donor <- table(SeuratObject$donor)
print(cells_per_donor)
pdf(paste0(plotsdir, ID, "_cells-per-donor.pdf"), width = 7, height = 5)
par(xpd = T)
bp <- barplot(cells_per_donor,
              las = 2, xlab = "Donor", ylab = "Number of Cells", main = ID)
text(x = bp, y = cells_per_donor + max(cells_per_donor*0.05),
     labels = cells_per_donor)
par(xpd = F)
dev.off()

saveRDS(SeuratObject, paste0(outdir, date, "_", ID, "_with_donor.RDS"))




# CombinedSeuratObject <- readRDS("scripts/Magpie_Analysis/britta/pipeline/out/2021_10_25_combined.RDS")
# 
# df <- data.frame(donor = character(0),
#                  cell_line = character(0),
#                  hipsci_line_name = character(0))
# for (ID in FinishedIDs){
#   
#   if (file.exists((paste0("outs/Magpie/2_Genotyping/", ID, "/donor_ids.tsv")))){
#     print(ID)
#     donor_ids <- as.data.frame(fread(paste0("outs/Magpie/2_Genotyping/", ID, "/donor_ids.tsv")))
#     row.names(donor_ids) <- paste(date, ID, "guides_assigned.RDS", donor_ids$cell, sep = "_")
#     donor_ids$cell_line <- gsub(donor_ids$donor_id, pattern = ".*-", replacement = "")
#     donor_ids$hipsci_line_name <- donor_ids$donor_id
#     donor_ids$donor <- unlist(mapply(donor_ids$donor_id, FUN = function(x){
#       x <- gsub(x, pattern = ".*-", replacement = "")
#       x <- gsub(x, pattern = "_.*", replacement = "")
#       return(x)
#     }))
#     df <- rbind(df, donor_ids %>% select(c("donor", "cell_line", "hipsci_line_name")))
#     }
#   
# }
# CombinedSeuratObject <- AddMetaData(CombinedSeuratObject, metadata = df, col.name = NULL)
# 
# saveRDS(CombinedSeuratObject, paste0("scripts/Magpie_Analysis/britta/pipeline/out/2021_10_25_combined.RDS"))


