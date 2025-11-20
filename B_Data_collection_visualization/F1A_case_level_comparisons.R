# ---- Libraries & setup ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(forcats)
  library(googlesheets4)
  library(rlang)
  library(ggplot2)
  library(ggpubr)
})

gs4_deauth()

source('/diskmnt/Users2/epeng/Projects/mCRC/scripts/jupyter_support_functions.R')
source('/diskmnt/Users2/epeng/Projects/mCRC/scripts/mCRC_colors.R')

output_dir <- getwd()
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

save_pdf <- function(plot, filename, width = 6, height = 6) {
  ggsave(
    filename = file.path(output_dir, filename),
    plot = plot,
    width = width, height = height,
    device = cairo_pdf # better text rendering
  )
}

# ---- Constants ----
sheet_url <- "https://docs.google.com/spreadsheets/d/1d4Bn3guQTzqf4tWx5fZFc772JC2GhpKVvsg70V7EPbQ/edit?gid=83959954#gid=83959954"

organ_cols <- c(Colorectum = '#C2B280', Rectum = '#604E97', Liver = 'brown', Lung = 'steelblue1', Other = '#d485b2')
assay_color_map <- c(WXS = "#463397", Bulk_RNA = "#91218c", snRNA = "#92ae31", snATAC = "#2b3514",
                     Visium = "#4277b6", Visium_HD = "#96cde6", Xenium = "#d32b1e",
                     Xenium_5K = "#e1a11a", CODEX = "#6f340d")

tissue_col <-  c(Primary = 'gold3', Metastasis = 'maroon3', Normal ='lightblue1')
organ_col  <- c(Colon = '#C2B280', Rectum = '#604E97', Liver = 'brown', Lung = 'steelblue1',
                Brain = 'bisque1', Adrenal = 'cyan4', Breast = 'violet', Spleen = 'lightcoral')
age_col    <- c('>50' = 'lightblue3', '<50' = 'lightblue1', 'NA' = 'gray75')
sex_col    <- c('M' = 'royalblue', 'F' = 'plum', 'NA' = 'gray75')
race_col   <- c('White' = 'darkgreen', 'African American' = 'firebrick2', 'NA' = 'gray75', 'Unknown' = 'gray75')
site_col   <- c('Left' = 'darkorange', 'Right' = 'cornflowerblue', 'Left/Right' = 'olivedrab3',
                'Rectum' = 'red3', 'Unknown' = 'gray75', 'NA' = 'gray75')
met_col    <- c('Yes' = 'maroon1', 'No' = 'cornflowerblue')
tx_col     <- c('Yes' = 'salmon', 'No' = 'antiquewhite', 'Unknown' = 'gray75', 'NA' = 'gray75')

# ---- Helpers ----
recode_availability <- function(x) dplyr::case_when(
  x %in% c("tumor_only","available") ~ 1L,
  x == "unavailable"                 ~ 0L,
  TRUE                               ~ NA_integer_
)

recode_countish <- function(x, allow = c("2","3","5","9")) {
  x <- as.character(x)
  dplyr::recode(x,
    !!!setNames(as.list(as.integer(allow)), allow),
    "available"   = 1L,
    "unavailable" = 0L,
    .default      = NA_integer_
  )
}

# Generic categorical plotter (bar or pie) that returns both plot & counts
plot_cat <- function(df, var, type = c("pie","bar"), title = NULL, palette = NULL,
                     order_levels = NULL, drop_labels = c("Unknown","NA")) {
  type <- match.arg(type)
  v <- enquo(var)

  tmp <- df %>%
    mutate(.val = as.character(!!v)) %>%
    filter(!is.na(.val), !(.val %in% drop_labels)) %>%
    mutate(.val = factor(.val, levels = if (is.null(order_levels)) unique(.val) else order_levels))

  counts <- tmp %>%
    count(.val, name = "n") %>%
    mutate(perc = n / sum(n) * 100)

  # Color scale resolver
  scale_fill_resolved <- function(p) {
    lvls <- levels(counts$.val)
    if (is.null(palette)) {
      p + scale_fill_brewer(palette = "Set2")
    } else if (is.character(palette) && length(palette) == 1) {
      p + scale_fill_brewer(palette = palette)
    } else {
      pal_vec <- palette
      # pad if unnamed & too short
      if (is.null(names(pal_vec)) && length(pal_vec) < length(lvls)) {
        pal_vec <- c(pal_vec, grDevices::hcl.colors(length(lvls) - length(pal_vec)))
      }
      if (!is.null(names(pal_vec))) pal_vec <- pal_vec[lvls]
      p + scale_fill_manual(values = pal_vec)
    }
  }

  if (nrow(counts) == 0) {
    return(list(plot = ggplot() + theme_void() + labs(title = paste0(title, " (no data)")), counts = counts))
  }

  if (type == "pie") {
    p <- ggplot(counts, aes(x = "", y = perc, fill = .val)) +
      geom_col(color = "white") +
      coord_polar(theta = "y") +
      geom_text(aes(label = paste0(.val, "\n", n, " (", round(perc, 1), "%)")),
                position = position_stack(vjust = 0.5), size = 3.5, lineheight = 0.9) +
      theme_void() +
      labs(title = title, fill = NULL)
    p <- scale_fill_resolved(p)
  } else {
    p <- ggplot(counts, aes(x = .val, y = n, fill = .val)) +
      geom_col(color = "black", width = 0.6) +
      geom_text(aes(label = paste0(n, " (", round(perc, 1), "%)")), vjust = -0.5, size = 3.5) +
      theme_minimal(base_size = 13) +
      labs(title = title, x = NULL, y = "Number of cases", fill = NULL) +
      theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
    p <- scale_fill_resolved(p)
  }

  list(plot = p, counts = counts)
}

# ---- Read tables ----
sample_tbl <- read_sheet(sheet_url, sheet = 'Sample_collections') %>%
  select(-c(MSI, APC_WXS, TP53_WXS, KRAS_WXS, Note))

clinical_tbl <- read_sheet(sheet_url, sheet = 'Clinical') %>%
  mutate(
    Age = if_else(Age_at_diagnosis > 50, '>50', '<50'),
    CRC_site = case_when(
      Site_in_Colon %in% c('Descending', 'Sigmoid', 'Left') ~ 'Left',
      Site_in_Colon %in% c('Ascending', 'Cecum', 'Transverse', 'Right') ~ 'Right',
      Site_in_Colon %in% c('Both') ~ 'Left/Right',
      TRUE ~ Site_in_Colon
    )
  ) %>%
  select(-c(SP_ID, Age_at_diagnosis, Site_in_Colon,
            MSI, APC, TP53, KRAS,
            Survival, FU_days, Note)) %>%
  dplyr::rename(Prior_Tx = Neoadjuvant_therapy)

# ---- Merge & standardize ----
sample_tbl <- sample_tbl %>%
  left_join(clinical_tbl, by = 'Case_ID') %>%
  mutate(
    Organ = if_else(Organ %in% c("Colon","Rectum"), "Colorectum", Organ),
    Bulk_DNA = recode_availability(Bulk_DNA),
    Bulk_RNA = recode_availability(Bulk_RNA),
    snRNA    = recode_countish(snRNA,    allow = c("2")),
    snATAC   = recode_countish(snATAC,   allow = c("2")),
    Visium   = recode_countish(Visium,   allow = c("2","3","5")),
    Visium_HD= recode_countish(Visium_HD),
    Xenium_regular = recode_countish(Xenium_regular, allow = c("2","9")),
    Xenium_5K      = recode_countish(Xenium_5K),
    Phenocycler    = recode_countish(Phenocycler, allow = c("2","3")),
    Pri_Met_Pair   = dplyr::recode(as.character(Pri_Met_Pair), "No" = 0L, "Yes" = 1L, .default = NA_integer_)
  ) %>%
  rename(Xenium = Xenium_regular, Pair = Pri_Met_Pair) %>%
  filter(Organ %in% c('Colorectum','Liver','Lung','Adrenal','Brain','Breast','Spleen','LN','Ovary','Pancreas')) %>%
  mutate(Organ = if_else(Organ %in% c('Adrenal','Brain','Breast','Spleen','LN','Ovary','Pancreas'), 'Other', Organ))

collection_tbl <- sample_tbl %>%
  select(Tissue_ID, Organ, Pair, Bulk_DNA, Bulk_RNA, snRNA, snATAC, Visium, Xenium, Xenium_5K, Phenocycler, Visium_HD)

collection_tbl2 <- collection_tbl %>%
  group_by(Organ) %>%
  summarise(
    WXS        = sum(Bulk_DNA,   na.rm = TRUE),
    Bulk_RNA   = sum(Bulk_RNA,   na.rm = TRUE),
    snRNA      = sum(snRNA,      na.rm = TRUE),
    snATAC     = sum(snATAC,     na.rm = TRUE),
    Visium     = sum(Visium,     na.rm = TRUE),
    Visium_HD  = sum(Visium_HD,  na.rm = TRUE),
    Xenium     = sum(Xenium,     na.rm = TRUE),
    Xenium_5K  = sum(Xenium_5K,  na.rm = TRUE),
    CODEX      = sum(Phenocycler,na.rm = TRUE),
    .groups = 'drop'
  )

# ---- Circular (donut+spike) plot (more robust) ----
ordered_orgs <- c("Liver","Lung","Colorectum","Other")

org4 <- collection_tbl2 %>%
  mutate(Organ = ifelse(Organ %in% ordered_orgs[-4], Organ, "Other")) %>%
  group_by(Organ) %>%
  summarise(across(where(is.numeric), ~sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(Organ = factor(Organ, levels = ordered_orgs)) %>%
  arrange(Organ)

assay_order <- c("WXS","Bulk_RNA","snRNA","snATAC","Visium","Visium_HD","Xenium","Xenium_5K","CODEX")
assay_vars  <- intersect(assay_order, setdiff(names(org4), "Organ"))

# Layout parameters
gap_frac <- 0.04; r_inner <- 0.30; r_outer <- 0.40; n_org <- length(ordered_orgs)

ring <- org4 %>%
  mutate(
    width  = 1 / n_org,
    xmid   = (as.integer(Organ) - 0.5) / n_org,
    xstart = xmid - width * (1 - gap_frac) / 2,
    xend   = xmid + width * (1 - gap_frac) / 2,
    ymin   = r_inner,
    ymax   = r_outer
  )

bars <- org4 %>%
  pivot_longer(all_of(assay_vars), names_to = "assay", values_to = "count") %>%
  mutate(assay = factor(assay, levels = assay_order)) %>%
  left_join(ring %>% select(Organ, xstart, xend), by = "Organ") %>%
  group_by(Organ) %>%
  arrange(assay, .by_group = TRUE) %>%
  mutate(
    w = xend - xstart,
    n_a = n(), i = row_number(),
    edge_pad_frac = 0.01, bar_fill = 0.96,
    w_eff  = w * (1 - 2 * edge_pad_frac),
    slot_w = w_eff / n_a,
    bar_w  = slot_w * bar_fill,
    xmin = xstart + w * edge_pad_frac + (i - 1) * slot_w + (slot_w - bar_w) / 2,
    xmax = xmin + bar_w
  ) %>%
  ungroup() %>%
  mutate(
    count = replace_na(count, 0L)
  )

max_count  <- max(bars$count, na.rm = TRUE)
max_count  <- ifelse(max_count == 0, 1, max_count) # avoid div-by-zero
radial_gap <- 0.01
height_mult <- 0.30

bars <- bars %>%
  mutate(
    ymin = r_outer + radial_gap,
    ymax = ymin + height_mult * (count / max_count)
  )

labels <- ring %>%
  transmute(
    Organ,
    x = (xstart + xend) / 2,
    y = (r_inner + r_outer) / 2,
    ang  = 90 - 360 * x,
    hjust = if_else(ang < -90, 1, 0),
    ang  = if_else(ang < -90, ang + 180, ang)
  )

p_circle <- ggplot() +
  geom_rect(data = ring,
            aes(xmin = xstart, xmax = xend, ymin = ymin, ymax = ymax, fill = Organ),
            colour = NA) +
  scale_fill_manual(values = organ_cols, name = "Organ") +
  ggnewscale::new_scale_fill() +
  geom_rect(data = bars,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = assay),
            alpha = 0.95, colour = NA) +
  scale_fill_manual(values = assay_color_map, name = "Assay") +
  geom_text(data = labels,
            aes(x = x, y = y, label = Organ, angle = ang, hjust = hjust),
            size = 4, fontface = "bold", colour = "white") +
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = r_inner,
           fill = "white", colour = NA) +
  coord_polar(theta = "x") +
  xlim(0, 1) +
  theme_void() +
  theme(legend.position = "bottom")

save_pdf(p_circle, "mCRC_assays_circle_plot.pdf", 6, 6)

# ---- Case-level tables & plots ----
case_tbl <- clinical_tbl %>%
  mutate(Met_status = case_when(
    Dx_with_Metastasis == 'No' ~ 'No_Met',
    Dx_with_Metastasis == 'Yes' & Lung_Met == 'No' ~ 'Met_wo_Lung',
    Lung_Met == 'Yes' ~ 'Lung_Met',
    TRUE ~ NA_character_
  ))

case_tbl3 <- case_tbl %>%
  mutate(CRC_site3 = case_when(
    CRC_site %in% c("Left","Left/Right") ~ "Left",
    CRC_site == "Right" ~ "Right",
    CRC_site == "Rectum" ~ "Rectum",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(CRC_site3), !is.na(Dx_with_Metastasis)) %>%
  mutate(
    CRC_site3 = factor(CRC_site3, levels = c("Left","Right","Rectum")),
    Dx_with_Metastasis = factor(Dx_with_Metastasis, levels = c("No","Yes"))
  )

# Pies
p_sex <- plot_cat(case_tbl, Sex, type = "pie", palette = sex_col,  order_levels = c("F","M"))$plot
p_age <- plot_cat(case_tbl, Age, type = "pie", palette = age_col,  order_levels = c("<50",">50"))$plot
p_met <- plot_cat(case_tbl, Dx_with_Metastasis, type = "pie", palette = met_col, order_levels = c("No","Yes"))$plot

save_pdf(p_met, "mCRC_sample_met_pieplot.pdf", 3, 3)
save_pdf(p_age, "mCRC_sample_age_pieplot.pdf", 6, 6)
save_pdf(p_sex, "mCRC_sample_sex_pieplot.pdf", 6, 6)

# ---- Pairwise Fisher across sites ----
site_levels <- levels(case_tbl3$CRC_site3)
pairs <- combn(site_levels, 2, simplify = FALSE)

pairwise_res <- map_dfr(pairs, function(pr) {
  sub <- filter(case_tbl3, CRC_site3 %in% pr)
  tab <- table(sub$CRC_site3, sub$Dx_with_Metastasis)
  ft  <- suppressWarnings(fisher.test(tab))
  tibble(
    group1 = pr[1],
    group2 = pr[2],
    p = ft$p.value,
    odds_ratio = unname(ft$estimate)
  )
})

# Bar plot with p-values
plot_df <- case_tbl3 %>%
  group_by(CRC_site3) %>%
  summarise(
    total_cases = n(),
    met_cases = sum(Dx_with_Metastasis == "Yes", na.rm = TRUE),
    prop_met = met_cases / total_cases,
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(met_cases, "/", total_cases),
    CRC_site3 = factor(CRC_site3, levels = c("Left","Right","Rectum"))
  )

ymax <- max(plot_df$prop_met, na.rm = TRUE)
step <- 0.06
upper <- ymax + 3 * step + 0.05

anno <- pairwise_res %>%
  mutate(
    group1 = factor(group1, levels = c("Left","Right","Rectum")),
    group2 = factor(group2, levels = c("Left","Right","Rectum")),
    y.position = c(ymax + step, ymax + 2 * step, ymax + 3 * step)[seq_len(n())],
    label = ifelse(is.na(p), "p = NA", paste0("p = ", signif(p, 3)))
  )

p1 <- ggplot(plot_df, aes(x = CRC_site3, y = prop_met, fill = CRC_site3)) +
  geom_col(width = 0.6, color = "black") +
  geom_text(aes(label = label), vjust = -0.5, size = 5) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent_format(),
    expand = expansion(mult = c(0, 0.18))
  ) +
  scale_fill_manual(values = site_col, guide = "none") +
  labs(
    x = "CRC Site",
    y = "Proportion with Metastasis at Diagnosis",
    title = "Metastasis at Diagnosis by CRC Site"
  ) +
  theme_mydefault() +
  theme(legend.position = "none") +
  stat_pvalue_manual(
    data = anno,
    label = "label",
    xmin = "group1", xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    inherit.aes = FALSE
  )

save_pdf(p1, "mCRC_primary_site_barplot.pdf", 8, 6)

# Prior Tx bar (using generic bar helper)
p_prior <- plot_cat(case_tbl, Prior_Tx, type = "bar", palette = tx_col, order_levels = c("Yes","No"))$plot
save_pdf(p_prior, "mCRC_sample_Tx_barplot.pdf", 3, 3)
