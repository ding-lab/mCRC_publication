library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(scales) 

xenium_anno = '/diskmnt/Projects/MetNet_analysis/Colorectal/data/Xenium_annotations/mCRC_N26_Xenium_banky_celltype_metadata_20250713.csv'
xenium_metadata = read.csv(xenium_anno, header = TRUE)
rownames(xenium_metadata) <- xenium_metadata$barcode
xenium_metadata <- xenium_metadata %>% select(-c(barcode))

xenium_metadata <- xenium_metadata%>% 
                   mutate(GTR = if_else(tns_label > 0, 'GTR', 'NAT')) %>% 
                   mutate(TMR = if_else(tn_label > 0, 'TMR', 'non-TMR')) %>% 
                   mutate(GTR_cat = case_when(
                        tn_label > 0 ~ 'TMR',
                        tns_label == 0 ~ 'APT',
                        TRUE ~ 'PSI'
                    ))

# Define the fibroblast classes
fibroblast_classes <- c("WNT5A_BMP", "WNT5A_infl", "stromal_fibroblast", "iCAF", "mCAF")

# Filter the data, group, and calculate proportions
fibroblast_proportions <- xenium_metadata %>%
  # Filter out 'Breast' organ and non-fibroblast cell types
  filter(Organ != "Breast", All_cell_type1 %in% fibroblast_classes) %>%
  
  # Group by the variables needed for counting
  group_by(Organ, GTR_cat, All_cell_type1) %>%
  
  # Count the number of cells (n) for each unique combination
  summarise(n = n(), .groups = "drop_last") %>%
  
  # Calculate the proportion of each fibroblast type within its GTR_cat/Organ group
  mutate(proportion = n / sum(n)) %>%
  
  # Arrange for better viewing
  ungroup() %>%
  arrange(Organ, GTR_cat, desc(proportion))

# Define the custom color palette and its inherent order
fibroblast_colors <- c(
  "WNT5A_BMP" = "#8dd3c7",
  "WNT5A_infl" = "#ffffb3",
  "iCAF" = "#bebada",
  "stromal_fibroblast" = "#fb8072",
  "mCAF" = "#80b1d3"
)

# Extract the names to use for setting factor levels (this ensures the order is respected)
cell_type_order <- names(fibroblast_colors)
gtr_cat_order <- c("TMR", "PSI", "APT")


## =========================================================================
## Factor Reordering Step (Updated)
## =========================================================================

fibroblast_proportions <- fibroblast_proportions %>%
  mutate(
    # Reorder the GTR_cat factor levels to be TMR, PSI, APT (from previous request)
    GTR_cat = factor(GTR_cat, levels = gtr_cat_order),
    
    # Reorder the All_cell_type1 factor levels based on your color palette
    All_cell_type1 = factor(All_cell_type1, levels = cell_type_order)
  )


## =========================================================================
## Stacked Bar Plot (Final Version)
## =========================================================================

bar_plot <- ggplot(fibroblast_proportions, 
                   aes(x = GTR_cat, y = proportion, fill = All_cell_type1)) +
  
  # Create the stacked bars
  # The stacking order and legend order are now controlled by the factor levels
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.2) +
  
  # Apply the custom fill colors
  scale_fill_manual(values = fibroblast_colors) +
  
  # Add labels to the bars for readability (only for proportions > 5%)
  geom_text(aes(label = ifelse(proportion > 0.05, 
                               scales::percent(proportion, accuracy = 0.1), 
                               "")), 
            position = position_stack(vjust = 0.5), 
            color = "black", size = 3) +
  
  # Separate the plots by Organ
  facet_wrap(~ Organ, scales = "free_x") +
  
  # Set labels and theme
  labs(title = "Compositional Proportion of Fibroblast Subtypes",
       subtitle = "Grouped by GTR Category and Organ (Excluding Breast)",
       x = "GTR Category (TMR, PSI, APT)",
       y = "Proportion of Total Fibroblasts",
       fill = "Fibroblast Subtype") +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold"))

out_dir = getwd()
pdf(file.path(out_dir, "barplot_fibroblast_subtype_by_organ_tissue_domain.pdf"), width = 5, height = 4)
bar_plot
dev.off()