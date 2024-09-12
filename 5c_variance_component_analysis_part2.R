
library(tidyverse)
library(reshape2)
library(data.table)

HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
section_name <- "5c_variance_component_analysis"
date <- "2024-05-23"
#date_out <- '2024-05-23'

#settings:
#n_genes_total <- 10

# options
args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
  if (arg == "--n_genes"){ n_genes_total <- as.numeric(args[[ind + 1]]) }  
  if (arg == "--redo"){ redo <- as.logical(args[[ind + 1]])}
}

# i/o
source(file.path(HomeFolder, "scripts/io/Magpie_io.R"))
outdir <- file.path(OutFolder, section_name)
plotsdir <- file.path(HTMLFolder, "pipeline/", section_name)
if(!dir.exists(outdir)) dir.create(outdir, recursive = T)
if(!dir.exists(plotsdir)) dir.create(plotsdir)
inlet_strengths <- c("Moderate", "Strong")

## ---- LoadData

latest <- gsub(sort(list.files(file.path(OutFolder, "5a_combine"), pattern = 'gene_list'), decreasing = T)[1], pattern = '_gene_list.txt', replacement = '')
genes2test <- fread(file.path(OutFolder, "5a_combine", paste0(latest, "_gene_list.txt")))$claudia

if (date == "2022-08-15"){
  var_explained_fnms <- grep(list.files(file.path(outdir, 'by_ind'), recursive = T, pattern = paste0("_n-genes-", n_genes_total)), pattern = date, value = T)
  var_explained_df <- lapply(file.path(outdir, 'by_ind', var_explained_fnms), fread) %>% bind_rows() %>% as.data.frame()
  var_explained_df <- reshape2::melt(var_explained_df, id.vars = c("downstream_gene", "downstream_gene_name"))
  colnames(var_explained_df) <- c("downstream_gene", "downstream_gene_name", "variable", "variance_explained")
  var_explained_df$variable <- as.character(var_explained_df$variable)
  fwrite(var_explained_df, paste0(outdir, '/', date, '_var_expl_df.tsv'), sep = '\t')
} else {
  var_explained_df <- fread(paste0(outdir, '/', date, '_variance_explained.tsv.gz'))
  var_explained_df <- reshape2::melt(var_explained_df, id.vars = c("downstream_gene_name"))
  colnames(var_explained_df) <- c( "downstream_gene_name", "variable", "variance_explained")
  var_explained_df$variable <- as.character(var_explained_df$variable)
  fwrite(var_explained_df, paste0(outdir, '/', date, '_var_expl_df.tsv'), sep = '\t')
}

## ----

# summarize
res <- lapply(c(0.01, 0.05, 0.1,0.2, 0.25,0.3, 0.4), function(th) {
  var_explained_df %>% group_by(variable) %>%
    summarise(n_genes = sum(variance_explained > th)) %>%
    mutate(`Variance Explained` = paste(">", th)) %>% 
    filter(variable != "Residual" & variable != "Residuals") %>%
    mutate(variable = ifelse(variable == "Batch", "5_Batch", variable)) %>%
    mutate(variable = ifelse(variable == "Batch:orig.ident", "6_Batch:orig.ident", variable)) %>%
    mutate(variable = ifelse(variable == "donor", "2_donor", variable)) %>%
    mutate(variable = ifelse(variable == "donor:cell_line", "3_donor:cell_line", variable)) %>%
    mutate(variable = ifelse(variable == "nCount_RNA", "7_nCount_RNA", variable)) %>%
    mutate(variable = ifelse(variable == "percent_MT", "8_percent_MT", variable)) %>%
    mutate(variable = ifelse(variable == "Phase", "4_Phase", variable)) %>%
    mutate(variable = ifelse(variable == "status", "1_status", variable)) %>%
    mutate(`Source of Variance` = ifelse(variable %in% c("2_donor", "3_donor:cell_line"), "Genetics",
                                         ifelse(variable %in% c("4_Phase"), "Cellular", 
                                                ifelse(variable %in% c("1_status"), "CRISPR", "Technical"))))
}) %>% bind_rows()
res$variable <- factor(res$variable)
levels(res$variable) <- c("Perturbation status", 
                          "Donor", "Cell line", 
                          "Cell cycle phase",
                          "Batch", "Inlet", "Number of UMIs", "% of Mitochondrial RNA")
write.csv(res,
          file = file.path(outdir,
                           paste0(date,"_", length(genes2test) ,"_genes-total.csv")))

#write.csv(res,
#          file = file.path(outdir,
#                           paste0(date_new, "_res_", feature_choice, "_",target_choice, "_",
#                                  inlet_group,"_",cells,"_", n_genes_total ,".csv")))


# make plots
#pdf(file.path(plotsdir,
#              paste0(date_new, "_var_comps_", feature_choice, "_",target_choice, "_",
#                     inlet_group,"_",cells,"_", n_genes_total ,".pdf")))
pdf(file.path(plotsdir,
              paste0(date_new, "_", length(genes2test) ,"_genes-total.pdf")), width = 6, height = 5)
ggplot(res, aes(x = variable, y = log10(n_genes+1) , alpha = `Variance Explained`, fill = `Source of Variance`)) +
  geom_bar(stat = "identity", width = 0.9, position = "dodge") +
  scale_fill_manual(values = c("Genetics" = "darkorchid4", "Cellular" = "cornflower blue", "Technical" = 'dark green', 'CRISPR' = 'darkred')) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
        axis.title.x = element_blank(), legend.position = "top",
        strip.background = element_rect(fill = "white"), legend.box="vertical", legend.margin=margin()) +
  #ggtitle(paste0("Var. expl.: ", target_choice, ",",inlet_group,",",cells)) +
  ylab(paste0("Log10 (Number of genes + 1)\n", length(genes2test), " highly expressed genes total"))
ggplot(res, aes(x = variable, y = n_genes , alpha = `Variance Explained`, fill = `Source of Variance`)) +
  geom_bar(stat = "identity", width = 0.9, position = "dodge") +
  scale_fill_manual(values = c("Genetics" = "darkorchid4", "Cellular" = "cornflower blue", "Technical" = 'dark green', 'CRISPR' = 'darkred')) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
        axis.title.x = element_blank(), legend.position = "top",
        strip.background = element_rect(fill = "white"), legend.box="vertical", legend.margin=margin()) +
  #ggtitle(paste0("Var. expl.: ", target_choice, ",",inlet_group,",",cells)) +
  ylab(paste0("Number of genes\n", length(genes2test), " highly expressed genes total"))
dev.off()


## idk what this does
# do additional downsampling analysis if frac_cell < 1
if(frac_cells < 1) {
  print("Starting down-sampling experiment")
  var_comps_downsampled <- lapply(1:n_repeats, function(i) {
    seed <- sample(1:100000, 1)
    set.seed(seed)
    cells_to_include <- sample(seq_len(nrow(meta_df)), frac_cells * nrow(meta_df))
    dd_per_gene <- mclapply(features2fit, function(genenm) {
      print(genenm)
      tmp <- meta_df
      tmp$gene_expr <- data[genenm,]
      tmp <- tmp[cells_to_include,]
      dd <- fit_lmm(data4fit = tmp, gene_name = genenm, inlet_group = inlet_group, target_choice = target_choice)
      dd$seed <- seed
      dd$frac_cells <- frac_cells
      return(dd)
    }, mc.cores = ncores)
    return(bind_rows(dd_per_gene))
  }) 
  var_comps_downsampled <- bind_rows(var_comps_downsampled)
  
  # save results
  write.csv(var_comps_downsampled,
            file = file.path(outdir,
                             paste0(date_new, "_var_comps_downsampled_", feature_choice, "_" ,target_choice, "_",inlet_group,"_",cells,"_", n_genes_total ,".csv")))
  
  # summarise
  res_downsampled <- lapply(c(0.1,0.2,0.3, 0.4), function(th) {
    var_comps_downsampled %>% group_by(variable, seed, frac_cells) %>%
      summarise(n_genes = sum(var_explained > th)) %>%
      mutate(variance.explained = paste(">", th)) %>% 
      filter(variable != "Residual") %>%
      mutate(variable = ifelse(variable == "orig.ident:batch", "inlet:batch", variable)) %>%
      mutate(variable = ifelse(variable == "Phase", "CellCycle", variable))
  }) %>% bind_rows()
  
  write.csv(res_downsampled,
            file = file.path(outdir,
                             paste0(date_new, "_res_downsampled_", feature_choice, "_",
                                    target_choice, "_",inlet_group,"_",cells,"_", n_genes_total ,".csv")))
  
  # make plot
  cols <- RColorBrewer::brewer.pal(9, "YlGnBu")[6:9]
  names(cols) <- sort(unique(res_downsampled$variance.explained))
  
  pdf(file.path(plotsdir,
                paste0(date_new, "_var_comps_downsampled_", feature_choice, "_",
                       target_choice, "_",inlet_group,"_",cells,"_", n_genes_total ,".pdf")))
  gg <- ggplot(res_downsampled, aes(x = variable, y = n_genes ,
                                    col = variance.explained,
                                    fill = variance.explained)) +
    # stat_summary(fun.data = "mean_cl_normal", size = 0.1) +
    geom_boxplot() +
    geom_point() +
    geom_point(data = res, col = "darkred") +
    # geom_bar(stat = "identity", width = 0.5,
    # data = res, alpha = 0, lty = 3) +
    theme_bw() + guides(fill ="none")  +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
          axis.title.x = element_blank(), legend.position = "top",
          strip.background = element_rect(fill = "white")) +
    ggtitle(paste0("Var. expl.: ", target_choice, ",",inlet_group,",",cells)) +
    ylab(paste("Number of genes (out of", n_genes_total, "hvg)")) +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) + facet_wrap(~variance.explained)
  print(gg)
  dev.off()
}


# sessionInfo
sessionInfo()