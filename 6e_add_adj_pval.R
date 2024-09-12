

suppressMessages(library(data.table))

#date <- "2021-12-31"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
DOWNSTREAMGENENAME = T

args <- commandArgs(trailingOnly = TRUE)
infolder <- args[[1]]
infile <- args[[2]]
outfile <- args[[3]]
DOWNSTREAMGENENAME <- as.logical(args[[4]])

setwd(HomeFolder)
#source(io_path)
#source(utils_path)


# load file
res <- fread(file.path(infolder, infile))

# add gene name
if (DOWNSTREAMGENENAME == T){
  res$downstream_gene_name <- unlist(lapply(strsplit(res$downstream_gene, ":"), "[[", 2))
}
  
#add adjusted p-value
res$pval_adj <- p.adjust(res$pval_lm, method = "BH")

#write
fwrite(res, file = file.path(infolder, outfile), row.names = F, quote = F, sep = "\t", compress = "gzip")




