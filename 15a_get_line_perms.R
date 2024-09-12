
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
if(exists("donors2include")){print('including donors: '); print(donors2include)}



setwd(HomeFolder)
source(io_path)
source(utils_path)

GuideMetadata <- fread(GuideMetadataPath)
LineMetadata <- fread(LineMetadataPath)
sexes <- LineMetadata$sex; names(sexes) <- LineMetadata$name

#outdir <- file.path(file.path(AnalysisOutFolder, section_name))
outdir <- file.path(file.path(OutFolder, section_name))
if(!dir.exists(outdir)) dir.create(outdir, recursive = T)
#plotsdir <- file.path(HTMLFolder, "pipeline", section_name) 
#if(!dir.exists(plotsdir)) dir.create(plotsdir, recursive = T)

## ---- LoadData

LineMetadata <- fread(LineMetadataPath)

target_df <- fread(paste0(OutFolder, "/7b_target_summary/", date, "_target_meta_data.csv"))
targets_per_line_set <- table(target_df$paired_lines)
lines_per_target <- names(targets_per_line_set)
lines_per_target <- data.frame(
  lines = lines_per_target,
  n_targets = as.numeric(targets_per_line_set),
  n_lines = unlist(lapply(strsplit(lines_per_target, ";"), length))
) %>% arrange(desc(n_lines)) %>%
  mutate(line_set_name = 1:length(lines_per_target))
fwrite(lines_per_target, paste0(outdir, '/', date, '_lines_per_target.tsv'))

n_lines_to_permute <- unique(lines_per_target$n_lines)

## target to line mapping
target_line_mapping <- target_df %>%
  mutate(line_set_name = lines_per_target$line_set_name[match(paired_lines, lines_per_target$lines)]) %>%
  dplyr::select(c('gene', 'paired_lines', 'n_paired_lines', 'line_set_name'))
fwrite(target_line_mapping, paste0(outdir, '/', date, '_target_line_mapping.tsv'))

## ---- Permute

for (n in unique(lines_per_target$n_lines)){
  ## get all possible permutations
  all_line_perms <- permutations(n = n, r = n)
  to_excl <- sweep(all_line_perms, 2, 1:dim(all_line_perms)[2], `-`)
  to_excl <- apply(to_excl, 1, FUN = function(x){length(which(x == 0))})
  all_line_perms <- all_line_perms[to_excl == 0,]
  line_selection_from_donor <- t(as.matrix(as.data.frame(lapply(strsplit(intToBin(0:(2 ^ n - 1)), split = ""), as.numeric))))
  line_selection_from_donor <- line_selection_from_donor[which(line_selection_from_donor[,1] == 1),]
  
  ### keep the selected line fixed
  ## don't include donor names yet
  line_meta <- data.frame(
    cell_line_ind = 1:(n*2),
    line_ind = rep(c(0,1), n),
    donor_ind = rep(1:n, rep(2, n))
  ) %>%
    mutate(to_match = paste0(donor_ind, "_", line_ind)) %>%
    mutate(to_match_pretty = paste0('donor_', donor_ind, '_line_', line_ind))
  
  
  all_perms <- lapply(1:dim(line_selection_from_donor)[1], FUN = function(selection_ind){
    #print(selection_ind)
    perm_set <- line_meta %>%
      mutate(to_permute = to_match %in% paste0(1:n, "_", line_selection_from_donor[selection_ind,])) %>%
      arrange(desc(to_permute))
    lines2permute <- filter(perm_set, to_permute)
    permuted_lines <- matrix(as.vector(lines2permute$cell_line_ind[all_line_perms]),
                             ncol = dim(all_line_perms)[2],
                             nrow = dim(all_line_perms)[1])
    unpermuted_lines <- matrix(perm_set %>% filter(to_permute == F) %>% .$cell_line_ind,byrow = T,
                               ncol = dim(all_line_perms)[2],
                               nrow = dim(all_line_perms)[1])
    donor_perms <- as.data.frame(suppressMessages(bind_cols(permuted_lines, 
                                                            unpermuted_lines)))
    colnames(donor_perms) <- paste0('donor_', perm_set$donor_ind, '_line_', perm_set$line_ind)
    donor_perms <- donor_perms[, order(colnames(donor_perms))]
    return(donor_perms)
  })
  all_perms <- as.data.frame(bind_rows(all_perms))
  
  #sort the new donor pairs
  possible_pairs <- combn(1:(2*n), 2) %>% t()
  possible_pairs <- paste0(possible_pairs[,1], ';', possible_pairs[,2])
  all_perms <- mapply(i = 1:n, FUN = function(i){
    subsetted_df <- all_perms %>% dplyr::select(contains(paste0("donor_", i)))
    subsetted_df <- t(Rfast::colSort(t(subsetted_df)))
    subsetted_df <- paste0(subsetted_df[,1], ";", subsetted_df[,2])
    subsetted_df <- match(subsetted_df, possible_pairs)
    return(subsetted_df)
  }, SIMPLIFY = F) %>% bind_cols() %>% as.data.frame()
  all_perms <- t(Rfast::colSort(t(all_perms)))
  ## return to the cell_line inds
  all_perms <- matrix(possible_pairs[all_perms], 
                      nrow = dim(all_perms)[1], ncol = dim(all_perms)[2])
  all_perms <- unique(all_perms)
  all_perms <- lapply(1:n, FUN = function(i){
    rtn <- data.frame(
      line_1 = unlist(lapply(strsplit(all_perms[,i], ";"), "[[", 1)),
      line_2 = unlist(lapply(strsplit(all_perms[,i], ";"), "[[", 2))
    )
    return(rtn)
  }) %>% bind_cols() %>% as.data.frame()
  colnames(all_perms) <- paste0('pair_', rep(1:n, rep(2,n)), '_line_', rep(1:2, n))
  saveRDS(all_perms, paste0(outdir, '/', date, '_', n, '_perms.RDS'))
  
  
  ## start dealing with different sets of lines
  line_sets2consider <- filter(lines_per_target, n_lines == n) %>% .$line_set_name
  for (line_set in line_sets2consider){
    donors <- lines_per_target %>% filter(line_set_name == as.character(line_set)) %>% .$lines
    donors <- unlist(strsplit(donors, ";"))
    lines <- LineMetadata %>% filter(donor %in% donors) %>% .$name %>% unique() %>% sort()
    line_pairings <- matrix(lines[as.numeric(as.matrix(all_perms))], 
                            nrow = dim(all_perms)[1],
                            ncol = dim(all_perms)[2])
    colnames(line_pairings) <- paste0('pair_', rep(1:n, rep(2,n)), '_line_', rep(1:2, n))
    saveRDS(line_pairings, paste0(outdir, '/', date, '_line_set-', line_set, '_', n, '_perms.RDS'))
    
    ## randomize
    set.seed(rnd_seed)
    line_pairings <- line_pairings[sample(x = 1:dim(line_pairings)[1], size = dim(line_pairings)[1], replace = F),]
    saveRDS(line_pairings, paste0(outdir, '/', date, '_line_set-', line_set, '_', n, '_perms_rnd-seed-', rnd_seed, '.RDS'))
    
  }
  
  
}








