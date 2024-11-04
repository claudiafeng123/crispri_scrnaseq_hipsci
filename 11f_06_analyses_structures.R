#!/usr/bin/env Rscript
# Analyse structures
source("src/config.R")
plots <- list()

### Load data
cor_complex <- readxl::read_xlsx("data/2022-08-15_high_target_complex_cor.xlsx", sheet = "2022-08-15_high_target_complex_")

uniprot_map <- read_tsv("meta/uniprot_mapping.tsv") %>% select(gene = From, uniprot = Entry)
cor_all <- read_tsv("data/Genome_Wide_Target_Cor_2022-08-15.tsv.gz") %>%
  left_join(uniprot_map, by = c("gene_1" = "gene"), relationship = "many-to-many") %>%
  rename(uniprot_1 = uniprot) %>%
  left_join(uniprot_map, by = c("gene_2" = "gene"), relationship = "many-to-many") %>%
  rename(uniprot_2 = uniprot)

meta <- readxl::read_xlsx("data/2022-08-15_high_target_complex_cor.xlsx", sheet = "modelling_targets", range = "A1:K27") %>%
  mutate(gene_2 = str_split(gene_2, "_"),
         uniprot_2 = str_split(uniprot_2, "_")) %>%
  left_join(select(cor_complex, gene_1, n_cells_1, complex_name, cor_all, cor_all_pval, guide_consistency), by = c("gene_1", "complex_name"))

dockq <- read_csv("data/pdockq.csv") %>%
  rename_with(~str_to_lower(str_replace_all(., "[ -+,]", "_"))) %>%
  select(-param) %>%
  extract(name, "complex", "./([^/]*)/.*") %>%
  separate(complex, str_c("gene", 1:7), "_", fill = "right", remove = FALSE) %>%
  slice_max(pdockq, n = 1, by = complex) %>%
  arrange(desc(pdockq)) %>%
  {bind_rows(
    ., 
    filter(., is.na(gene3)) %>% mutate(t = gene2, gene2 = gene1, gene1 = t) %>% select(-t)
  )}

foldx_summary <- read_tsv("data/foldx/interaction.tsv") %>%
  rename_with(~str_replace_all(str_to_lower(.), " ", "_")) %>%
  extract(pdb, c("uniprot_1", "uniprot_2", "extras"), "data/foldx/([A-Z0-9]*)_([A-Z0-9]*)_?([A-Z0-9_]*)_Repair.pdb") %>%
  filter(extras == "") %>%
  select(-extras)

hits <- readxl::read_xlsx("data/2022-08-15_high_target_complex_cor.xlsx", sheet = "hits") %>%
  left_join(filter(dockq, is.na(gene3)) %>% select(gene1, gene2, pdockq:len2), by = c("uniprot_1" = "gene1", "uniprot_2" = "gene2")) %>%
  left_join(foldx_summary, by = c("uniprot_1", "uniprot_2")) %T>%
  write_tsv("data/hits.tsv")

background_dockq <- bind_rows(
  `Random` = read_csv("data/random_pdockq.csv") %>%
    rename_with(~str_to_lower(str_replace_all(., "[ +,]", "_"))) %>%
    separate(name, c("gene1", "gene2"), "-"),
  `Negatome` = read_table("data/negatome.pdockq") %>%
    separate(pair, c("gene1", "gene2"), "_"),
  CORUM = read_csv("data/CORUM_pdockq.csv") %>%
    rename_with(~str_to_lower(str_replace_all(., "[ +,]", "_"))) %>%
    separate(name, c("gene1", "gene2"), "-"),
  .id = "dist"
) %T>%
  write_tsv("data/background_dockq.tsv")

#### Analysis
# Cor vs dockq
plots$dockq_vs_cor <- select(meta, uniprot_1, uniprot_2, cor_all, cor_all_pval, guide_consistency) %>%
  unnest(uniprot_2) %>%
  left_join(filter(dockq, is.na(gene3)), ., by = c("gene1" = "uniprot_1", "gene2" = "uniprot_2")) %>%
  drop_na(pdockq, cor_all, guide_consistency) %>%
  ggplot(aes(x = pdockq, y = cor_all, colour = guide_consistency)) +
  geom_point() +
  labs(x = "Predicted DockQ", y = expression("Pearson's"~rho)) +
  scale_colour_distiller(name = "Guide Consistency", palette = "Blues", direction = 1, limits = c(0,1))

plots$dockq_vs_cor_all <- inner_join(select(cor_all, uniprot_1, uniprot_2, gene_1, gene_2, cor_all),
           select(filter(dockq, is.na(gene3)), gene1, gene2, pdockq, complex),
           by = c("uniprot_1"="gene1", "uniprot_2"="gene2")) %T>%
  {print(cor.test(.$cor_all, .$pdockq))} %>%
  ggplot(aes(x = cor_all, y = pdockq)) +
  geom_point()

plots$dockq_vs_cor_deg <- inner_join(select(cor_all, uniprot_1, uniprot_2, gene_1, gene_2, cor_all_deg),
                                     select(filter(dockq, is.na(gene3)), gene1, gene2, pdockq, complex),
                                     by = c("uniprot_1"="gene1", "uniprot_2"="gene2")) %T>%
  {print(cor.test(.$cor_all_deg, .$pdockq))} %>%
  ggplot(aes(x = cor_all_deg, y = pdockq)) +
  geom_point()

# Distribution vs random
plots$dockq_distribution <- ggplot() +
  stat_density(mapping = aes(x = pdockq, y = after_stat(scaled), colour = dist), data = background_dockq, geom = "line", position = "identity") +
  stat_density(mapping = aes(x = pdockq, y = after_stat(scaled), colour = "Candidates"), data = hits, geom = "line", position = "identity") +
  geom_point(mapping = aes(x = pdockq, fill = "Miss"), data = filter(hits, !interface), shape = 21, stroke = 0, y = 0.0125, size = 2) +
  geom_point(mapping = aes(x = pdockq, fill = complex), data = filter(hits, interface), shape = 21, stroke = 0, y = 0.025, size = 2) +
  scale_colour_manual(values = c("Random" = "grey80", "Negatome" = "grey50", "CORUM" = "green", "Candidates" = "#377eb8"),
                      labels = c("Random" = "Random", "Negatome" = "Negatome", "CORUM" = "CORUM", "Candidates" = "Candidates")) +
  scale_fill_manual(values = c("Miss" = "black", "No" = "#E41A1C", "Clash" = "#FF7F00", "Good" = "cornflowerblue", "Known" = "purple"),
                      labels = c("Miss" = "No Binding", "No" = "Incoherent Complex", "Clash" = "Clashing Complex",
                                 "Good" = "Good Complex", "Known" = "Known Complex")) +
  labs(x = "Predicted DockQ", y = "Scaled Density") +
  lims(x = c(0,1)) +
  theme(legend.title = element_blank())
plots$dockq_distribution_png <- labeled_plot(plots$dockq_distribution, file_format = "png")

# Test distributions
# ks.test(hits$pdockq, filter(background_dockq, dist == "Random")$pdockq, alternative = "less")
# ks.test(hits$pdockq, filter(background_dockq, dist == "Negatome")$pdockq, alternative = "less")
# ks.test(hits$pdockq, filter(background_dockq, dist == "CORUM")$pdockq, alternative = "less")
# 
# ks.test(hits$pdockq, filter(background_dockq, dist == "Random")$pdockq)
# ks.test(hits$pdockq, filter(background_dockq, dist == "Negatome")$pdockq)
# ks.test(hits$pdockq, filter(background_dockq, dist == "CORUM")$pdockq)
# 
# wilcox.test(hits$pdockq, filter(background_dockq, dist == "Random")$pdockq, alternative = "greater")
# wilcox.test(hits$pdockq, filter(background_dockq, dist == "Negatome")$pdockq, alternative = "greater")
# wilcox.test(hits$pdockq, filter(background_dockq, dist == "CORUM")$pdockq, alternative = "greater")

# Binding energy distribution
plots$complex_interface_residues <- ggplot(hits, aes(x = interface, y = interface_residues, fill = complex)) +
  geom_boxplot() +
  labs(x = "Complex Formed", y = "Number of Interface Residues") +
  scale_fill_brewer(palette = "Dark2")

plots$interaction_energy <- ggplot(hits, aes(x = pdockq, y = clamp(interaction_energy, upper = 10), colour = complex)) +
  geom_point() +
  labs(x = "pDockQ", y = "Interaction Energy (clamped to <= 10)") +
  scale_fill_brewer(palette = "Dark2")


# Save plots
save_plotlist(plots, "figures/structure_analysis", overwrite = "all")  
