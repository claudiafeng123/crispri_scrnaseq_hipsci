
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(gtools))
suppressMessages(require(R.utils))
suppressMessages(library(Rfast))
suppressMessages(library(doParallel))

# set relevent i/o paths
section_name <- "15a_get_line_perms"


# set relevent i/o paths
date <- "2022-10-05"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
ExperimentName <- "Pica"

# options for testing
io_path <- paste0("scripts/io/", ExperimentName, "_io.R")
utils_path <- "scripts/io/Magpie_Utils.R"
rnd_seed <- 0


args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  ##needed arguments
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--experiment_name"){ ExperimentName <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--utils_path"){ utils_path <- args[[ind + 1]] }
  if (arg == "--io_path"){ io_path <- args[[ind + 1]] }
  if (arg == "--rnd_seed"){ rnd_seed <- as.numeric(args[[ind + 1]]) }
}

print(paste0(section_name))
print(paste("date:", date))
print(paste("home folder:", HomeFolder))



setwd(HomeFolder)
source(io_path)
source(utils_path)

outdir <- file.path(file.path(OutFolder, section_name))
if(!dir.exists(outdir)) dir.create(outdir, recursive = T)

## ---- LoadData

lines_per_target <- fread(paste0(outdir, '/', date, '_lines_per_target.tsv'))
num_lines <- unique(lines_per_target$n_lines)

## ---- Permute

total_perms <- mapply(i = num_lines, FUN = function(i){
  fnm <- list.files(outdir, 
                    pattern = paste0( "_", i, "_perms"), full.names = T)[1]
  all_permutations <- readRDS(fnm)
  rtn <- dim(all_permutations)[1]
  return(rtn)
})
total_perms <- data.frame(
  n_lines = num_lines,
  n_perms = total_perms)

fwrite(total_perms, paste0(outdir, '/', date, '_total_perms.tsv'), sep = '\t')



