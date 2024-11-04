#!/usr/bin/env Rscript
# Analyse complexes
source("src/config.R")
plots <- list()

meta <- read_tsv("data/2022-08-15_complex_meta.tsv") %>%
  mutate(genes_in_complex = str_split(genes_in_complex, "_")) %>%
  unnest(genes_in_complex) %>%
  select(gene = genes_in_complex, complex = complex_name, everything()) %>%
  drop_na(complex)

lfc <- read_tsv("data/2022-08-15_lfcs_complex_genes_only.tsv") %>%
  left_join(select(meta, gene, complex) %>% distinct(), by = c("target" = "gene"))

complex_cor <- function(tbl, key) {
  message(key)
  tbl <- group_by(tbl, downstream_gene_name) %>%
    filter(any(pval_adj < 0.1)) %>%
    ungroup() %>%
    select(target, downstream_gene_name, lfc) %>%
    distinct()
  
  if (nrow(tbl) == 0) {
    return(tibble(gene1 = NA, gene2 = NA, cor = NA))
  }
  
  mat <- pivot_wider(tbl, names_from = target, values_from = lfc) %>%
    tblhelpr::tibble_to_matrix(-downstream_gene_name, row_names = .$downstream_gene_name)
  
  # P val calculation from https://stackoverflow.com/a/46209337
  r <- cor(mat)
  n <- t(!is.na(mat)) %*% !is.na(mat)
  t <- (r*sqrt(n-2))/sqrt(1-r^2)
  p <- 2*(1 - pt(abs(t),(n-2)))
  se <- sqrt((1-r*r)/(n-2))
  
  tibble(gene1 = rep(rownames(r), each = length(rownames(r))),
         gene2 = rep(rownames(r), times = length(rownames(r))),
         idx = matrix(c(gene1, gene2), ncol = 2),
         n = n[idx],
         cor = r[idx],
         t = t[idx],
         p = p[idx],
         se = se[idx]) %>%
    select(-idx) %>%
    mutate(p_adj = p.adjust(p))
}

cors <- drop_na(lfc, complex, lfc, pval_adj) %>%
  group_by(complex) %>%
  group_modify(complex_cor) %>%
  ungroup()

affected_genes <- drop_na(lfc, complex, lfc, pval_adj) %>%
  group_by(complex, downstream_gene_name) %>%
  mutate(shared_downstream_gene = sum(pval_adj < 0.1) > 2) %>%
  group_by(complex) %>%
  mutate(total_affected = n_distinct(downstream_gene_name[pval_adj < 0.1]),
         total_shared_affected = n_distinct(downstream_gene_name[shared_downstream_gene])) %>%
  group_by(complex, target, total_affected, total_shared_affected) %>%
  summarise(affected_genes = sum(pval_adj < 0.1),
            shared_affected_genes = sum(pval_adj[shared_downstream_gene] < 0.1),
            .groups = "drop_last") %>%
  mutate(prop_affected_genes = affected_genes / total_affected,
         prop_shared_affected_genes = shared_affected_genes / total_shared_affected) %>%
  ungroup()

complexes <- rename(cors, gene = gene1) %>%
  filter(!gene == gene2) %>%
  group_by(complex, gene) %>%
  summarise(n_cor = sum(p_adj < 0.05, na.rm = TRUE),
            prop_cor = n_cor / n(),
            mean_cor = mean(cor, na.rm = TRUE),
            .groups = "drop") %>%
  left_join(affected_genes, by = c("complex", "gene"="target")) %>%
  left_join(distinct(meta, complex, .keep_all=TRUE) %>% select(-gene), by = c("complex")) %>%
  left_join(filter(lfc, target == downstream_gene_name) %>% select(complex, gene = target, lfc, pval_adj),
            by = c("complex", "gene"))

# How does number of downstream effects overall relate to number shared by complex
plots$shared_vs_total_affected <- ggplot(complexes, aes(x = prop_affected_genes, y = prop_shared_affected_genes,
                                                        colour = affected_genes)) +
  geom_point() +
  scale_colour_distiller(trans = "pseudo_log", direction = 1)

# Mean correlation relates to number of genes affected?
plots$mean_cor_vs_affected <- filter(complexes, total_shared_affected > 10, complex_size > 4,
                                     max(prop_cor, na.rm = TRUE) > 0.1) %>%
  ggplot(aes(x = mean_cor, y = prop_affected_genes, colour = prop_cor)) +
  geom_point()

plots$dist_cor_groups <- filter(complexes, total_shared_affected > 10, complex_size > 4,
                                max(prop_cor, na.rm = TRUE) > 0.1) %>%
  ggplot(aes(colour = n_cor > 4, x = prop_affected_genes, y = after_stat(scaled))) +
  stat_density(geom = "line", position = "identity")

# Similar for shared genes?
plots$mean_cor_vs_shared_affected <- filter(complexes, total_shared_affected > 10, complex_size > 4,
                                     max(prop_cor, na.rm = TRUE) > 0.1) %>%
  ggplot(aes(x = mean_cor, y = prop_shared_affected_genes, colour = prop_cor)) +
  geom_point()

plots$dist_cor_groups_shared <- filter(complexes, total_shared_affected > 10, complex_size > 4,
                                max(prop_cor, na.rm = TRUE) > 0.1) %>%
  ggplot(aes(colour = mean_cor > 0.25, x = prop_shared_affected_genes, y = after_stat(scaled))) +
  stat_density(geom = "line", position = "identity")

# Do only KO'd genes correlate?
plots$ko_vs_cor <- filter(complexes, total_shared_affected > 10, complex_size > 4,
                           max(prop_cor, na.rm = TRUE) > 0.1, !is.na(pval_adj)) %>%
  ggplot(aes(x = pval_adj < 0.1, y = mean_cor)) + 
  geom_boxplot(notch = TRUE, varwidth = TRUE) +
  stat_compare_means(comparisons = list(c("TRUE", "FALSE")))

plots$lfc_vs_cor <- filter(complexes, total_shared_affected > 10, complex_size > 4,
                           max(prop_cor, na.rm = TRUE) > 0.1, !is.na(pval_adj)) %>%
  ggplot(aes(x = lfc, y = mean_cor, colour = log10(pval_adj))) + 
  geom_point() +
  scale_colour_gradientn(colours = c("darkred", "red","antiquewhite", "blue"),
                         values = scales::rescale(c(-150, -10, -1, 0), c(0,1)))

# Save plots
save_plotlist(plots, "figures/correlation", overwrite = "all")