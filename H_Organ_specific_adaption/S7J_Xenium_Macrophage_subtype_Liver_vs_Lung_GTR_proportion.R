library(tidyverse)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(forcats)
library(scales)
library(rstatix)
library(Seurat)
library(qs)
library(RColorBrewer)

output_dir = getwd()

metadata_path = '/diskmnt/Projects/MetNet_analysis/Colorectal/data/Xenium_annotations/mCRC_N26_Xenium_banky_celltype_metadata_20250713.csv'
meta_df <- read_csv(metadata_path)
colnames(meta_df)[1] <- "barcode"
meta_df <- meta_df %>% column_to_rownames("barcode")

meta_df = meta_df %>% 
          mutate(GTR = if_else(tns_label > 0, 'GTR', 'NAT')) %>% 
          mutate(TMR = if_else(tn_label > 0, 'TMR', 'non-TMR')) %>% 
          mutate(GTR_cat = case_when(
              tn_label > 0 ~ 'TMR',
              tns_label == 0 ~ 'NAT',
              TRUE ~ 'TSR')) %>% 
          mutate(GTR_NB = case_when(
              GTR == 'NAT' ~ 'NAT',
              TRUE ~ neighborhoods)) %>% 
          mutate(tns_mask = case_when(
              tns_label > 0 ~ 'GTR',
              TRUE ~ 'APT')) %>% 
          mutate(GTR_cat2 = case_when(
              tns_label == 0 ~ 'APT',
              tn_label == 0 ~ 'PSIZ',
              t_label > 0 ~ 'Tumor',
              TRUE ~ 'Necrosis')) %>% 
          mutate(TMR_tumor_cat = case_when(
              tn_label > 0 & All_cell_type1 == 'Non-canonical' ~ 'Pro-invasive',
              tn_label > 0 & All_cell_type1 %in% c('Proliferative-like', 'Stem-like', 'Intestine-like') ~ 'Canonical',
              TRUE ~ 'Other'))

targets <- c("IFN_TAM", "SPP1_TAM", "M2_TAM", "RTM_KC", "RTM_alveolar", "TAM", "classical_monocyte", "PMN-like_monocyte" )
organ_col <- c(Liver = "brown", Lung = "steelblue1")

gtr_df <- meta_df %>% filter(GTR == "GTR", Organ %in% c("Liver","Lung"))
totals <- gtr_df %>% count(sample_id, Organ, name = "n_total")

sample_counts <- gtr_df %>%
  dplyr::filter(All_cell_type1 != 'Doublet') %>% 
  count(sample_id, Organ, All_cell_type1, name = "n") %>%
  filter(All_cell_type1 %in% targets) %>%
  complete(sample_id, Organ, All_cell_type1 = targets, fill = list(n = 0)) %>%
  left_join(totals, by = c("sample_id","Organ")) %>%
  mutate(prop = n / n_total)

# Wilcoxon (raw p) per subtype + label positions
p_anno <- sample_counts %>%
  group_by(All_cell_type1) %>%
  wilcox_test(prop ~ Organ) %>%
  mutate(p_label = paste0("p = ", signif(p, 3))) %>%
  left_join(
    sample_counts %>% group_by(All_cell_type1) %>% summarise(y = max(prop, na.rm = TRUE) * 1.05),
    by = "All_cell_type1"
  )

# Prepare per-sample proportions (denominator = all GTR cells in each sample)
gtr_df <- meta_df %>% filter(GTR == "GTR", Organ %in% c("Liver","Lung"))
totals <- gtr_df %>% count(sample_id, Organ, name = "n_total")

sample_counts <- gtr_df %>%
  count(sample_id, Organ, All_cell_type1, name = "n") %>%
  filter(All_cell_type1 %in% targets) %>%
  complete(sample_id, Organ, All_cell_type1 = targets, fill = list(n = 0)) %>%
  left_join(totals, by = c("sample_id","Organ")) %>%
  mutate(prop = n / n_total)

# Wilcoxon (raw p) per subtype + label positions
p_anno <- sample_counts %>%
  group_by(All_cell_type1) %>%
  wilcox_test(prop ~ Organ) %>%
  mutate(p_label = paste0("p = ", signif(p, 3))) %>%
  left_join(
    sample_counts %>% group_by(All_cell_type1) %>% summarise(y = max(prop, na.rm = TRUE) * 1.05),
    by = "All_cell_type1"
  )

# Plot
ggplot(sample_counts, aes(x = Organ, y = prop, fill = Organ)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(aes(color = Organ), width = 0.15, height = 0, size = 1.6, alpha = 0.85) +
  facet_wrap(~ All_cell_type1, nrow = 2) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = organ_col) +
  scale_color_manual(values = organ_col) +
  labs(x = NULL, y = "Proportion of GTR cells",
       title = "Macrophage/Monocyte subtype proportions among GTR cells by organ (per sample)") +
  geom_text(data = transform(p_anno, x = 1.5),
            aes(x = x, y = y, label = p_label),
            inherit.aes = FALSE, size = 3.2) +
  theme_bw(base_size = 12) +   # <-- brings grids back
  theme(strip.text = element_text(face = "bold"),
        strip.background = element_rect(fill = "grey92", colour = NA)) -> p


pdf(file.path(output_dir, "Xenium_GTR_Myeloid_subtype_boxplot.pdf"), width = 6, height = 4)
p
dev.off()
