#!/usr/bin/env Rscript
# Local (neighborhood) KRAS VAF to evaluate intra-sample mutation heterogeneity.
#
# For each tumor cell in all 8 KRAS-mutant Xenium cases, sum KRAS ALT/WT transcripts from tumor cells within
# 100, 200, 300, 400, and 500 um, then assign that neighborhood VAF to the cell.
#
# Interpretation:
#   - Unimodal histogram tightly around the sample-wide VAF, at all radii
#     => spatially homogeneous KRAS (sampling noise, not clones)
#   - Spread that shrinks from 100 um to 500 um
#     => sparse transcript noise, not biological patches
#   - Bimodal histograms or spatial high/low patches that persist at 500 um
#     => real local heterogeneity
#
# Run in screen:
#   screen -S kras_local_vaf
#   cd /diskmnt/Projects/MetNet_analysis/Colorectal/mCRC_Manuscript_Script/F_Tumor_subtypes
#   mkdir -p logs
#   Rscript 15_Xenium_KRAS_local_VAF_neighborhood.R > logs/Xenium_KRAS_local_VAF_neighborhood.log 2>&1

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(ggplot2)
  library(qs)
  library(scales)
})

logmsg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ..., "\n")
  flush.console()
}

# -----------------------------------------------------------------------------
# Paths / parameters
# -----------------------------------------------------------------------------
mCRC_path <- Sys.getenv(
  "MCRC_XENIUM_SLIM_QS",
  unset = "/diskmnt/Projects/MetNet_analysis/Colorectal/epeng/Spatial_driver/mCRC/1_xenium_slim_spatial_dirver_crc_5K.qs"
)
metadata_path <- Sys.getenv(
  "MCRC_XENIUM_METADATA",
  unset = "/diskmnt/Projects/MetNet_analysis/Colorectal/data/Xenium_annotations/mCRC_N26_Xenium_banky_celltype_metadata_20250624.csv"
)
vaf_root <- Sys.getenv(
  "MCRC_KRAS_VAF_DIR",
  unset = "/diskmnt/Projects/MetNet_analysis/Colorectal/mCRC_Manuscript_Script/F_Tumor_subtypes/Xenium_local_VAF_output"
)
output_dir <- file.path(vaf_root, "Xenium_local_VAF")
spatial_dir <- file.path(output_dir, "spatial_by_sample")
cmets_bin_dir <- file.path(vaf_root, "local_VAF_CMETS")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(spatial_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cmets_bin_dir, recursive = TRUE, showWarnings = FALSE)

# All KRAS-mutant tumor cases in the slim object (8 after tumor subset).
# Confirmed from Layers() after Broad_cell_type2 == Tumor; CMD* samples have no tumor cells.
preferred_samples <- c(
  "CM1799-A1-Th1Fp1U1",
  "PE0319C1-A4-S1Fp1U1",
  "S15-1909-A5U1",
  "S15-1909-D2U1",
  "S16-38794-A4U2",
  "S16-38794-E2U1",
  "S17-16442-AFR1U1",
  "S20-46186-A16U1"
)

radii_um <- c(100, 200, 300, 400, 500)
min_depth <- 10
bin_um <- 20
kras_features <- c(
  "KRAS-p.G12D-ALT:T",
  "KRAS-p.G12V-ALT:A",
  "KRAS-p.G12V-WT"
)

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
sample_to_fov <- function(sample_id, fov_names) {
  candidates <- c(
    sample_id,
    gsub("-", ".", sample_id),
    gsub("-", "_", sample_id)
  )
  hit <- candidates[candidates %in% fov_names]
  if (length(hit) == 0) {
    hit <- fov_names[gsub("[._]", "-", fov_names) == sample_id]
  }
  if (length(hit) == 0) NA_character_ else hit[[1]]
}

# Viridis on white, with larger points so the tissue fill is readable.
spatial_bg <- "white"
vaf_na_color <- "grey55"

vaf_color_scale <- function(name = "KRAS VAF") {
  scale_color_viridis_c(
    option = "viridis",
    begin = 0.15,
    end = 1,
    limits = c(0, 1),
    na.value = vaf_na_color,
    name = name
  )
}

vaf_fill_scale <- function(name = "KRAS VAF") {
  scale_fill_viridis_c(
    option = "viridis",
    begin = 0.15,
    end = 1,
    limits = c(0, 1),
    na.value = vaf_na_color,
    name = name
  )
}

spatial_map_theme <- function(aspect = FALSE) {
  th <- theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", color = "black"),
      plot.subtitle = element_text(hjust = 0.5, size = 8, color = "grey30"),
      legend.position = "right",
      legend.title = element_text(color = "black"),
      legend.text = element_text(color = "black"),
      legend.background = element_rect(fill = spatial_bg, color = NA),
      legend.key = element_rect(fill = spatial_bg, color = NA),
      strip.text = element_text(size = 11, face = "bold", color = "black"),
      strip.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = spatial_bg, color = NA),
      panel.background = element_rect(fill = spatial_bg, color = NA)
    )
  if (aspect) th <- th + theme(aspect.ratio = 1)
  th
}

get_fov_coords <- function(obj, fov) {
  coords <- tryCatch(
    GetTissueCoordinates(obj, image = fov),
    error = function(e) {
      tryCatch(GetTissueCoordinates(obj[[fov]]), error = function(e2) NULL)
    }
  )
  if (is.null(coords) || nrow(coords) == 0) return(NULL)
  coords <- as.data.frame(coords)
  if (!"cell" %in% colnames(coords) && !is.null(rownames(coords))) {
    coords$cell <- rownames(coords)
  }
  if (!"x" %in% colnames(coords) && "imagecol" %in% colnames(coords)) {
    coords$x <- coords$imagecol
  }
  if (!"y" %in% colnames(coords) && "imagerow" %in% colnames(coords)) {
    coords$y <- coords$imagerow
  }
  coords
}

align_coords_to_cells <- function(coords, cell_ids) {
  if (any(coords$cell %in% cell_ids)) {
    return(coords %>% filter(cell %in% cell_ids))
  }
  coords$cell_short <- sub(".*_", "", coords$cell)
  map <- tibble(cell = cell_ids, cell_short = sub(".*_", "", cell_ids))
  coords %>%
    inner_join(map, by = "cell_short") %>%
    transmute(x = x, y = y, cell = cell.y)
}

# Disk-sum on a spatial grid. Needed for 100-500 um: a k-NN cap would truncate
# dense tumor (at 100 um, median neighbors were already ~200-300).
neighborhood_vaf <- function(xy, alt, wt, radii, bin_um, min_depth) {
  n <- nrow(xy)
  x <- xy[, 1]
  y <- xy[, 2]
  x0 <- min(x)
  y0 <- min(y)
  col <- pmax(1L, as.integer(floor((x - x0) / bin_um) + 1L))
  row <- pmax(1L, as.integer(floor((y - y0) / bin_um) + 1L))
  nc <- max(col)
  nr <- max(row)
  logmsg("    grid:", nr, "x", nc, "bins of", bin_um, "um")

  agg <- tibble(row = row, col = col, alt = alt, wt = wt) %>%
    group_by(row, col) %>%
    summarise(alt = sum(alt), wt = sum(wt), n = n(), .groups = "drop")

  grid_alt <- matrix(0, nr, nc)
  grid_wt <- matrix(0, nr, nc)
  grid_n <- matrix(0, nr, nc)
  grid_alt[cbind(agg$row, agg$col)] <- agg$alt
  grid_wt[cbind(agg$row, agg$col)] <- agg$wt
  grid_n[cbind(agg$row, agg$col)] <- agg$n

  add_shift <- function(dst, src, dy, dx) {
    i1 <- max(1L, 1L + dy)
    i2 <- min(nr, nr + dy)
    j1 <- max(1L, 1L + dx)
    j2 <- min(nc, nc + dx)
    if (i1 > i2 || j1 > j2) return(dst)
    dst[i1:i2, j1:j2] <- dst[i1:i2, j1:j2] + src[(i1 - dy):(i2 - dy), (j1 - dx):(j2 - dx)]
    dst
  }

  out <- vector("list", length(radii))
  names(out) <- as.character(radii)

  for (r in radii) {
    logmsg("    disk sum radius", r, "um")
    rp <- as.integer(ceiling(r / bin_um))
    r2 <- r * r
    sum_alt <- matrix(0, nr, nc)
    sum_wt <- matrix(0, nr, nc)
    sum_n <- matrix(0, nr, nc)
    for (dy in (-rp):rp) {
      for (dx in (-rp):rp) {
        if ((dy * bin_um)^2 + (dx * bin_um)^2 > r2) next
        sum_alt <- add_shift(sum_alt, grid_alt, dy, dx)
        sum_wt <- add_shift(sum_wt, grid_wt, dy, dx)
        sum_n <- add_shift(sum_n, grid_n, dy, dx)
      }
    }
    alt_sum <- sum_alt[cbind(row, col)]
    wt_sum <- sum_wt[cbind(row, col)]
    depth <- alt_sum + wt_sum
    out[[as.character(r)]] <- tibble(
      radius_um = r,
      n_neighbors = sum_n[cbind(row, col)],
      local_ALT = alt_sum,
      local_WT = wt_sum,
      local_depth = depth,
      local_VAF = ifelse(depth >= min_depth, alt_sum / depth, NA_real_)
    )
  }
  bind_rows(out)
}

# -----------------------------------------------------------------------------
# Load object
# -----------------------------------------------------------------------------
logmsg("Loading slim Xenium object")
mCRC <- qread(mCRC_path)

logmsg("Reading metadata")
meta_df <- read_csv(metadata_path, show_col_types = FALSE)
colnames(meta_df)[1] <- "barcode"
meta_df <- meta_df %>% column_to_rownames("barcode")
rownames(meta_df) <- gsub("__", "_", rownames(meta_df))

mCRC <- AddMetaData(
  mCRC,
  metadata = meta_df,
  col.name = c("tn_label", "Broad_cell_type2", "All_cell_type2")
)
mCRC$FolderPath <- NULL
rm(meta_df)
gc()

logmsg("Subsetting tumor cells in TMRs (tn_label > 0)")
mCRC_tumor <- subset(mCRC, subset = Broad_cell_type2 == "Tumor" & tn_label > 0)
rm(mCRC)
gc()

mCRC_tumor$CMETS_group <- case_when(
  mCRC_tumor$All_cell_type2 == "Non-canonical" ~ "CMETS",
  mCRC_tumor$All_cell_type2 %in% c("Stem-like", "Proliferative-like", "Intestine-like") ~ "Non-CMETS",
  TRUE ~ NA_character_
)

fov_names <- Images(mCRC_tumor)
logmsg("FOVs:", paste(fov_names, collapse = ", "))

avail <- gsub("^counts\\.", "", grep("^counts\\.", Layers(mCRC_tumor), value = TRUE))
spatial_samples <- c(
  preferred_samples[preferred_samples %in% avail],
  setdiff(avail, preferred_samples)
)
logmsg("Tumor samples to analyze (n =", length(spatial_samples), "):", paste(spatial_samples, collapse = ", "))

# -----------------------------------------------------------------------------
# Per-sample neighborhood VAF
# -----------------------------------------------------------------------------
local_list <- list()
sample_vaf_tbl <- list()

for (sid in spatial_samples) {
  logmsg("===", sid, "===")
  layer_name <- paste0("counts.", sid)
  if (!layer_name %in% Layers(mCRC_tumor)) {
    logmsg("WARNING: missing layer", layer_name)
    next
  }

  fov <- sample_to_fov(sid, fov_names)
  if (is.na(fov)) {
    logmsg("WARNING: no FOV for", sid)
    next
  }

  coords <- get_fov_coords(mCRC_tumor, fov)
  if (is.null(coords)) {
    logmsg("WARNING: no coordinates for", fov)
    next
  }

  counts_mat <- GetAssayData(mCRC_tumor, assay = "Xenium", layer = layer_name)
  coords <- align_coords_to_cells(coords, colnames(counts_mat))
  cells <- coords$cell
  if (length(cells) < 50) {
    logmsg("WARNING: too few aligned cells:", length(cells))
    next
  }

  feat <- kras_features[kras_features %in% rownames(counts_mat)]
  get_vec <- function(f) {
    if (f %in% rownames(counts_mat)) as.numeric(counts_mat[f, cells]) else rep(0, length(cells))
  }
  alt_g12d <- get_vec("KRAS-p.G12D-ALT:T")
  alt_g12v <- get_vec("KRAS-p.G12V-ALT:A")
  wt_vec <- get_vec("KRAS-p.G12V-WT")

  sum_g12d <- sum(alt_g12d)
  sum_g12v <- sum(alt_g12v)
  sum_wt <- sum(wt_vec)
  vaf_g12d <- ifelse(sum_g12d + sum_wt > 0, sum_g12d / (sum_g12d + sum_wt), NA_real_)
  vaf_g12v <- ifelse(sum_g12v + sum_wt > 0, sum_g12v / (sum_g12v + sum_wt), NA_real_)
  kras_call <- case_when(
    !is.na(vaf_g12d) & vaf_g12d >= 0.5 & (is.na(vaf_g12v) | vaf_g12d >= vaf_g12v) ~ "G12D",
    !is.na(vaf_g12v) & vaf_g12v >= 0.5 ~ "G12V",
    TRUE ~ NA_character_
  )
  if (is.na(kras_call)) {
    logmsg("WARNING: no KRAS call for", sid)
    next
  }

  alt_vec <- if (kras_call == "G12D") alt_g12d else alt_g12v
  sample_vaf <- if (kras_call == "G12D") vaf_g12d else vaf_g12v
  logmsg("  KRAS_call =", kras_call, "sample VAF =", sprintf("%.3f", sample_vaf), "n_cells =", length(cells))

  sample_vaf_tbl[[sid]] <- tibble(
    Sample = sid,
    KRAS_call = kras_call,
    sample_VAF = sample_vaf,
    sample_ALT = sum(alt_vec),
    sample_WT = sum_wt,
    n_tumor_cells = length(cells)
  )

  xy <- as.matrix(coords[, c("x", "y")])
  nb <- neighborhood_vaf(xy, alt_vec, wt_vec, radii_um, bin_um, min_depth)

  local_list[[sid]] <- bind_cols(
    tibble(
      Sample = sid,
      KRAS_call = kras_call,
      sample_VAF = sample_vaf,
      cell = cells,
      x = coords$x,
      y = coords$y,
      tn_label = as.numeric(mCRC_tumor@meta.data[cells, "tn_label"]),
      CMETS_group = mCRC_tumor@meta.data[cells, "CMETS_group"],
      cell_ALT = alt_vec,
      cell_WT = wt_vec
    )[rep(seq_len(nrow(coords)), times = length(radii_um)), ],
    nb
  )

  rm(counts_mat, xy)
  gc()
}

local_df <- bind_rows(local_list)
sample_vaf_df <- bind_rows(sample_vaf_tbl)
local_keep <- local_df %>% filter(!is.na(local_VAF))

logmsg("Cells with neighborhood depth >=", min_depth, ":", nrow(local_keep), "/", nrow(local_df))

write.csv(
  sample_vaf_df,
  file.path(output_dir, "Xenium_KRAS_local_VAF_sample_summary.csv"),
  row.names = FALSE
)

local_stats <- local_keep %>%
  group_by(Sample, KRAS_call, radius_um, sample_VAF) %>%
  summarise(
    n_cells = n(),
    VAF_mean = mean(local_VAF),
    VAF_sd = sd(local_VAF),
    VAF_IQR = IQR(local_VAF),
    VAF_min = min(local_VAF),
    VAF_max = max(local_VAF),
    frac_within_0.05 = mean(abs(local_VAF - sample_VAF) <= 0.05),
    frac_within_0.10 = mean(abs(local_VAF - sample_VAF) <= 0.10),
    median_depth = median(local_depth),
    median_neighbors = median(n_neighbors),
    .groups = "drop"
  ) %>%
  arrange(Sample, radius_um)

write.csv(
  local_stats,
  file.path(output_dir, "Xenium_KRAS_local_VAF_stats.csv"),
  row.names = FALSE
)
print(as.data.frame(local_stats), row.names = FALSE)

# Do not write the full per-cell table (can be tens of millions of rows).
# A downsampled table is enough to rebuild plots.
set.seed(1)
local_ds <- local_keep %>%
  group_by(Sample, radius_um) %>%
  group_modify(~ dplyr::slice_sample(.x, n = min(20000L, nrow(.x)))) %>%
  ungroup()
write.csv(
  local_ds,
  file.path(output_dir, "Xenium_KRAS_local_VAF_downsampled.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 1. Histogram: number of cells vs local VAF (the plot you asked for)
# -----------------------------------------------------------------------------
logmsg("Plotting histograms")
p_hist <- ggplot(local_keep, aes(x = local_VAF)) +
  geom_histogram(binwidth = 0.02, boundary = 0, fill = "steelblue", color = "white", linewidth = 0.15) +
  geom_vline(aes(xintercept = sample_VAF), linetype = "dashed", color = "red", linewidth = 0.4) +
  facet_grid(Sample ~ paste0(radius_um, " um"), scales = "free_y") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(
    x = "Neighborhood KRAS VAF",
    y = "Number of cells",
    title = "Distribution of local KRAS VAF around each tumor cell",
    subtitle = "Dashed red line = sample-wide VAF. Neighborhoods with < 10 KRAS transcripts are excluded."
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 8),
    plot.subtitle = element_text(size = 9)
  )

ggsave(file.path(output_dir, "Xenium_KRAS_local_VAF_histogram.pdf"), p_hist, width = 12, height = 18)
ggsave(file.path(output_dir, "Xenium_KRAS_local_VAF_histogram.png"), p_hist, width = 12, height = 18, dpi = 300)

# -----------------------------------------------------------------------------
# 2. Density overlay of the three radii (does spread shrink with scale?)
# -----------------------------------------------------------------------------
p_dens <- ggplot(local_keep, aes(x = local_VAF, color = factor(radius_um), fill = factor(radius_um))) +
  geom_density(alpha = 0.15, linewidth = 0.7) +
  geom_vline(aes(xintercept = sample_VAF), linetype = "dashed", color = "red", linewidth = 0.4) +
  facet_wrap(~ Sample, ncol = 4, scales = "free_y") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_color_brewer(palette = "Dark2", name = "Radius") +
  scale_fill_brewer(palette = "Dark2", name = "Radius") +
  labs(
    x = "Neighborhood KRAS VAF",
    y = "Density",
    title = "Local VAF spread vs neighborhood size",
    subtitle = "If spread shrinks from 100 um to 500 um, it is sparse-count noise rather than clones"
  ) +
  theme_bw()

ggsave(file.path(output_dir, "Xenium_KRAS_local_VAF_density.pdf"), p_dens, width = 14, height = 8)
ggsave(file.path(output_dir, "Xenium_KRAS_local_VAF_density.png"), p_dens, width = 14, height = 8, dpi = 300)

# -----------------------------------------------------------------------------
# 3. SD of local VAF vs radius (summary of heterogeneity by scale)
# -----------------------------------------------------------------------------
p_sd <- ggplot(local_stats, aes(x = factor(radius_um), y = VAF_sd, group = Sample, color = Sample)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2.5) +
  labs(
    x = "Neighborhood radius (um)",
    y = "SD of local KRAS VAF",
    title = "Local VAF heterogeneity as a function of spatial scale"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "Xenium_KRAS_local_VAF_sd_vs_radius.pdf"), p_sd, width = 9, height = 5.5)
ggsave(file.path(output_dir, "Xenium_KRAS_local_VAF_sd_vs_radius.png"), p_sd, width = 9, height = 5.5, dpi = 300)

# -----------------------------------------------------------------------------
# 4. Spatial maps of local VAF (shared 0-1 color scale)
#    Plot all tumor cells: depth < 10 is grey, matching TMR maps.
# -----------------------------------------------------------------------------
logmsg("Plotting spatial maps")
set.seed(1)
spatial_plot_df <- local_df %>%
  group_by(Sample, radius_um) %>%
  group_modify(~ dplyr::slice_sample(.x, n = min(80000L, nrow(.x)))) %>%
  ungroup() %>%
  mutate(
    Sample = factor(Sample, levels = spatial_samples),
    radius_lab = factor(paste0(radius_um, " um"), levels = paste0(radii_um, " um")),
    pt_size = ifelse(as.character(Sample) == "S15-1909-D2U1", 0.5, 0.25)
  )

tmr_centroids_local <- spatial_plot_df %>%
  filter(!is.na(tn_label), tn_label > 0) %>%
  group_by(Sample, radius_lab, tn_label) %>%
  summarise(
    x = median(x, na.rm = TRUE),
    y = median(y, na.rm = TRUE),
    .groups = "drop"
  )

p_spatial <- ggplot(spatial_plot_df, aes(x = x, y = y, color = local_VAF)) +
  geom_point(aes(size = pt_size), stroke = 0) +
  scale_size_identity() +
  geom_text(
    data = tmr_centroids_local,
    aes(x = x, y = y, label = tn_label),
    inherit.aes = FALSE,
    color = "black",
    fontface = "bold",
    size = 2.8
  ) +
  facet_grid(Sample ~ radius_lab, scales = "free") +
  vaf_color_scale("Local\nKRAS VAF") +
  labs(
    title = "Neighborhood KRAS VAF around each tumor cell",
    subtitle = "Shared 0-1 viridis scale. Grey = neighborhood with < 10 KRAS transcripts.",
    x = NULL,
    y = NULL
  ) +
  spatial_map_theme(aspect = TRUE)

ggsave(
  file.path(output_dir, "Xenium_KRAS_local_VAF_spatial.png"),
  p_spatial,
  width = 18,
  height = 26,
  dpi = 200,
  bg = spatial_bg
)
ggsave(
  file.path(output_dir, "Xenium_KRAS_local_VAF_spatial.pdf"),
  p_spatial,
  width = 18,
  height = 26,
  bg = spatial_bg
)

# Per-sample maps: 100-400 um, no TMR labels (one 2x2 figure per sample)
logmsg("Plotting per-sample 100-400 um maps (no TMR labels)")
spatial_100_400_dir <- file.path(output_dir, "spatial_100-400um")
dir.create(spatial_100_400_dir, recursive = TRUE, showWarnings = FALSE)
radii_plot <- c(100, 200, 300, 400)
big_dot_samples <- c(
  "CM1799-A1-Th1Fp1U1",
  "PE0319C1-A4-S1Fp1U1",
  "S15-1909-A5U1",
  "S16-38794-A4U2"
)

for (sid in spatial_samples) {
  dfr <- local_df %>%
    filter(Sample == sid, radius_um %in% radii_plot) %>%
    mutate(radius_lab = factor(paste0(radius_um, " um"), levels = paste0(radii_plot, " um")))
  if (nrow(dfr) == 0) next
  pt_size <- if (sid %in% big_dot_samples) 1.2 else if (sid == "S15-1909-D2U1") 0.55 else 0.28
  dfr <- dfr %>% arrange(desc(is.na(local_VAF)))
  p_s <- ggplot(dfr, aes(x = x, y = y, color = local_VAF)) +
    geom_point(size = pt_size, stroke = 0) +
    facet_wrap(~ radius_lab, ncol = 2) +
    vaf_color_scale("Local\nKRAS VAF") +
    labs(
      title = sid,
      subtitle = paste0(
        unique(dfr$KRAS_call),
        "  |  sample VAF = ",
        sprintf("%.3f", unique(dfr$sample_VAF)[1]),
        "  |  100-400 um  |  grey = depth < 10"
      )
    ) +
    spatial_map_theme(aspect = TRUE)
  ggsave(
    file.path(spatial_100_400_dir, paste0("Xenium_KRAS_local_VAF_spatial_100-400um_", sid, ".png")),
    p_s,
    width = 9.5,
    height = 8.8,
    dpi = 300,
    bg = spatial_bg
  )
  ggsave(
    file.path(spatial_100_400_dir, paste0("Xenium_KRAS_local_VAF_spatial_100-400um_", sid, ".pdf")),
    p_s,
    width = 9.5,
    height = 8.8,
    bg = spatial_bg
  )
}

# Per-sample maps at 100-500 um
for (sid in unique(local_df$Sample)) {
  pt_size <- if (sid == "S15-1909-D2U1") 0.7 else 0.35
  for (r in radii_um) {
    dfr <- local_df %>% filter(Sample == sid, radius_um == r)
    if (nrow(dfr) == 0) next
    tmr_cent <- dfr %>%
      filter(!is.na(tn_label), tn_label > 0) %>%
      group_by(tn_label) %>%
      summarise(
        x = median(x, na.rm = TRUE),
        y = median(y, na.rm = TRUE),
        .groups = "drop"
      )
    n_colored <- sum(!is.na(dfr$local_VAF))
    p_r <- ggplot(dfr, aes(x = x, y = y, color = local_VAF)) +
      geom_point(size = pt_size, stroke = 0) +
      geom_text(
        data = tmr_cent,
        aes(x = x, y = y, label = tn_label),
        inherit.aes = FALSE,
        color = "black",
        fontface = "bold",
        size = 3.5
      ) +
      coord_fixed() +
      vaf_color_scale("Local\nKRAS VAF") +
      labs(
        title = paste0(sid, "  |  ", r, " um neighborhood VAF"),
        subtitle = paste0(
          unique(dfr$KRAS_call),
          "  |  sample VAF = ",
          sprintf("%.3f", unique(dfr$sample_VAF)),
          "  |  colored cells = ", n_colored, "/", nrow(dfr),
          "  |  local SD = ",
          sprintf("%.3f", sd(dfr$local_VAF, na.rm = TRUE))
        ),
        x = NULL,
        y = NULL
      ) +
      spatial_map_theme()
    ggsave(
      file.path(spatial_dir, paste0("Xenium_KRAS_local_VAF_spatial_", r, "um_", sid, ".png")),
      p_r,
      width = 8,
      height = 7,
      dpi = 300,
      bg = spatial_bg
    )
    ggsave(
      file.path(spatial_dir, paste0("Xenium_KRAS_local_VAF_spatial_", r, "um_", sid, ".pdf")),
      p_r,
      width = 8,
      height = 7,
      bg = spatial_bg
    )
  }
}

# -----------------------------------------------------------------------------
# 5. Hex-binned mean local VAF at 200 um
# -----------------------------------------------------------------------------
p_hex <- ggplot(
  local_keep %>% filter(radius_um == 200),
  aes(x = x, y = y, z = local_VAF)
) +
  stat_summary_hex(fun = mean, bins = 80) +
  facet_wrap(~ Sample, scales = "free", ncol = 4) +
  vaf_fill_scale("Mean local\nKRAS VAF") +
  labs(
    title = "Hex-binned mean 200 um neighborhood VAF",
    x = NULL,
    y = NULL
  ) +
  spatial_map_theme(aspect = TRUE)

ggsave(
  file.path(output_dir, "Xenium_KRAS_local_VAF_hexbin_200um.png"),
  p_hex,
  width = 16,
  height = 10,
  dpi = 300,
  bg = spatial_bg
)
ggsave(
  file.path(output_dir, "Xenium_KRAS_local_VAF_hexbin_200um.pdf"),
  p_hex,
  width = 16,
  height = 10,
  bg = spatial_bg
)

# -----------------------------------------------------------------------------
# 6. CMETS vs Non-CMETS composition across local VAF categories
#
# Primary (absolute): low_depth, <75%, 75-90%, >90%
# Written to CMETS_KRAS/local_VAF_CMETS/.
# -----------------------------------------------------------------------------
logmsg("CMETS vs Non-CMETS local VAF category composition")

cmets_pal <- c("CMETS" = "firebrick2", "Non-CMETS" = "palegreen4")
abs_levels <- c("low_depth", "<75%", "75-90%", ">90%")
qtl_levels <- c("low_depth", "Q1 (lowest)", "Q2", "Q3", "Q4 (highest)")

wilson_ci <- function(x, n, z = 1.96) {
  prop <- ifelse(n > 0, x / n, NA_real_)
  lo <- hi <- rep(NA_real_, length(n))
  ok <- !is.na(n) & n > 0
  if (any(ok)) {
    p <- x[ok] / n[ok]
    nn <- n[ok]
    denom <- 1 + z^2 / nn
    center <- (p + z^2 / (2 * nn)) / denom
    half <- z * sqrt(p * (1 - p) / nn + z^2 / (4 * nn^2)) / denom
    lo[ok] <- pmax(0, center - half)
    hi[ok] <- pmin(1, center + half)
  }
  tibble(prop = prop, prop_lo = lo, prop_hi = hi)
}

assoc_one <- function(tab) {
  tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]
  n <- sum(tab)
  empty <- tibble(
    test = NA_character_,
    statistic = NA_real_,
    p_value = NA_real_,
    cramers_v = NA_real_,
    n = n,
    n_vaf_cats = nrow(tab)
  )
  if (nrow(tab) < 2 || ncol(tab) < 2 || n < 10) return(empty)
  chi <- suppressWarnings(chisq.test(tab))
  k <- min(nrow(tab), ncol(tab))
  v <- sqrt(as.numeric(chi$statistic) / (n * (k - 1)))
  tibble(
    test = "chi-square",
    statistic = as.numeric(chi$statistic),
    p_value = chi$p.value,
    cramers_v = v,
    n = n,
    n_vaf_cats = nrow(tab)
  )
}

bin_theme <- theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey95")
  )

analyze_cmets_bins <- function(df, cat_col, cat_levels, scheme, file_stub, xlab, title_extra, out_dir) {
  d <- df %>%
    mutate(vaf_cat = factor(.data[[cat_col]], levels = cat_levels))

  bin_counts <- d %>%
    count(Sample, radius_um, radius_lab, vaf_cat, CMETS_group, name = "n") %>%
    complete(
      nesting(Sample, radius_um, radius_lab),
      vaf_cat,
      CMETS_group,
      fill = list(n = 0)
    ) %>%
    group_by(Sample, radius_um, radius_lab, vaf_cat) %>%
    mutate(
      n_bin = sum(n),
      prop_in_bin = ifelse(n_bin > 0, n / n_bin, NA_real_)
    ) %>%
    group_by(Sample, radius_um, radius_lab, CMETS_group) %>%
    mutate(
      n_group = sum(n),
      prop_in_group = ifelse(n_group > 0, n / n_group, NA_real_)
    ) %>%
    ungroup() %>%
    mutate(scheme = scheme)

  overall_cmets <- d %>%
    count(Sample, radius_um, radius_lab, CMETS_group, name = "n") %>%
    group_by(Sample, radius_um, radius_lab) %>%
    summarise(
      n_cells = sum(n),
      n_cmets = sum(n[CMETS_group == "CMETS"]),
      overall_cmets_frac = n_cmets / n_cells,
      .groups = "drop"
    )

  cmets_in_bin <- bin_counts %>%
    filter(CMETS_group == "CMETS") %>%
    transmute(scheme, Sample, radius_um, radius_lab, vaf_cat, n_cmets = n, n_bin) %>%
    mutate(wilson_ci(n_cmets, n_bin)) %>%
    left_join(overall_cmets, by = c("Sample", "radius_um", "radius_lab")) %>%
    mutate(enrich = prop - overall_cmets_frac)

  tests <- bind_rows(
    d %>%
      group_by(Sample, radius_um, radius_lab) %>%
      group_modify(~ assoc_one(table(.x$vaf_cat, .x$CMETS_group))) %>%
      ungroup() %>%
      mutate(contrast = "all_bins"),
    d %>%
      filter(vaf_cat != "low_depth") %>%
      group_by(Sample, radius_um, radius_lab) %>%
      group_modify(~ assoc_one(table(.x$vaf_cat, .x$CMETS_group))) %>%
      ungroup() %>%
      mutate(contrast = "no_low_depth")
  )

  q14_levels <- c("Q1 (lowest)", "Q4 (highest)")
  if (all(q14_levels %in% cat_levels)) {
    tests <- bind_rows(
      tests,
      d %>%
        filter(vaf_cat %in% q14_levels) %>%
        mutate(vaf_cat = droplevels(vaf_cat)) %>%
        group_by(Sample, radius_um, radius_lab) %>%
        group_modify(~ assoc_one(table(.x$vaf_cat, .x$CMETS_group))) %>%
        ungroup() %>%
        mutate(contrast = "Q1_vs_Q4")
    )
  }

  tests <- tests %>%
    mutate(scheme = scheme) %>%
    group_by(scheme, contrast) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    ungroup() %>%
    arrange(contrast, Sample, radius_um)

  write.csv(bin_counts, file.path(out_dir, paste0(file_stub, "_counts.csv")), row.names = FALSE)
  write.csv(cmets_in_bin, file.path(out_dir, paste0(file_stub, "_prop_in_bin.csv")), row.names = FALSE)
  write.csv(tests, file.path(out_dir, paste0(file_stub, "_tests.csv")), row.names = FALSE)

  logmsg(scheme, "association tests:")
  print(
    tests %>% select(scheme, contrast, Sample, radius_um, n, n_vaf_cats, cramers_v, p_value, p_adj),
    n = 80
  )

  p_stack <- ggplot(bin_counts, aes(x = vaf_cat, y = n, fill = CMETS_group)) +
    geom_col(position = "fill", width = 0.85, color = NA) +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual(values = cmets_pal) +
    facet_grid(Sample ~ radius_lab) +
    labs(
      title = paste0("CMETS vs Non-CMETS composition in each local VAF category", title_extra),
      subtitle = "Each bar is 100% of tumor cells in that Sample x radius x VAF bin. low_depth = neighborhood KRAS transcripts < 10.",
      x = xlab,
      y = "Proportion of cells in bin",
      fill = NULL
    ) +
    bin_theme

  ggsave(file.path(out_dir, paste0(file_stub, "_prop_in_bin.pdf")), p_stack, width = 16, height = 20)
  ggsave(file.path(out_dir, paste0(file_stub, "_prop_in_bin.png")), p_stack, width = 16, height = 20, dpi = 300)

  p_frac <- ggplot(cmets_in_bin, aes(x = vaf_cat, y = prop)) +
    geom_hline(
      data = overall_cmets,
      aes(yintercept = overall_cmets_frac),
      linetype = "dashed",
      color = "grey30",
      linewidth = 0.4
    ) +
    geom_col(fill = "firebrick2", width = 0.8, alpha = 0.85) +
    geom_errorbar(aes(ymin = prop_lo, ymax = prop_hi), width = 0.2, linewidth = 0.35) +
    geom_text(
      data = ~ dplyr::filter(.x, n_bin > 0),
      aes(label = n_bin, y = pmax(dplyr::coalesce(prop_hi, prop), 0)),
      vjust = -0.35,
      size = 2.1,
      color = "grey20"
    ) +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.12))) +
    facet_grid(Sample ~ radius_lab) +
    labs(
      title = paste0("CMETS fraction in each local VAF category", title_extra),
      subtitle = "Dashed line = overall CMETS fraction in that sample. Error bars = Wilson 95% CI. Numbers = cells in bin.",
      x = xlab,
      y = "CMETS proportion"
    ) +
    bin_theme +
    theme(legend.position = "none")

  ggsave(file.path(out_dir, paste0(file_stub, "_fraction_by_bin.pdf")), p_frac, width = 16, height = 20)
  ggsave(file.path(out_dir, paste0(file_stub, "_fraction_by_bin.png")), p_frac, width = 16, height = 20, dpi = 300)

  p_group <- ggplot(bin_counts, aes(x = vaf_cat, y = prop_in_group, fill = CMETS_group)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.05))) +
    scale_fill_manual(values = cmets_pal) +
    facet_grid(Sample ~ radius_lab) +
    labs(
      title = paste0("Local VAF category distribution within CMETS vs Non-CMETS", title_extra),
      subtitle = "Each group sums to 100% within a Sample x radius.",
      x = xlab,
      y = "Proportion of cells in group",
      fill = NULL
    ) +
    bin_theme

  ggsave(file.path(out_dir, paste0(file_stub, "_group_distribution.pdf")), p_group, width = 16, height = 20)
  ggsave(file.path(out_dir, paste0(file_stub, "_group_distribution.png")), p_group, width = 16, height = 20, dpi = 300)

  invisible(list(counts = bin_counts, prop = cmets_in_bin, tests = tests))
}

cat_df <- local_df %>%
  filter(CMETS_group %in% c("CMETS", "Non-CMETS")) %>%
  mutate(
    Sample = factor(Sample, levels = spatial_samples),
    CMETS_group = factor(CMETS_group, levels = c("CMETS", "Non-CMETS")),
    radius_lab = factor(paste0(radius_um, " um"), levels = paste0(radii_um, " um")),
    vaf_abs = factor(
      case_when(
        is.na(local_VAF) | local_depth < min_depth ~ "low_depth",
        local_VAF <= 0.75 ~ "<75%",
        local_VAF <= 0.90 ~ "75-90%",
        TRUE ~ ">90%"
      ),
      levels = abs_levels
    )
  ) %>%
  group_by(Sample, radius_um) %>%
  mutate(
    vaf_q = {
      q <- rep(NA_integer_, n())
      ok <- !is.na(local_VAF)
      if (any(ok)) q[ok] <- ntile(local_VAF[ok], 4L)
      q
    },
    vaf_qtl = factor(
      case_when(
        is.na(local_VAF) | local_depth < min_depth ~ "low_depth",
        vaf_q == 1 ~ "Q1 (lowest)",
        vaf_q == 2 ~ "Q2",
        vaf_q == 3 ~ "Q3",
        vaf_q == 4 ~ "Q4 (highest)"
      ),
      levels = qtl_levels
    )
  ) %>%
  ungroup()

logmsg("Absolute VAF bin occupancy:")
print(
  cat_df %>% count(Sample, radius_um, vaf_abs) %>% arrange(Sample, radius_um, vaf_abs),
  n = 120
)

q_cuts <- cat_df %>%
  filter(!is.na(local_VAF)) %>%
  group_by(Sample, radius_um, radius_lab) %>%
  summarise(
    n_depth_ge10 = n(),
    q25 = quantile(local_VAF, 0.25, na.rm = TRUE),
    q50 = quantile(local_VAF, 0.50, na.rm = TRUE),
    q75 = quantile(local_VAF, 0.75, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(
  q_cuts,
  file.path(cmets_bin_dir, "Xenium_KRAS_local_VAF_CMETS_quartile_cutpoints.csv"),
  row.names = FALSE
)
logmsg("Within-sample local VAF quartile cutpoints:")
print(q_cuts, n = 40)

analyze_cmets_bins(
  cat_df,
  cat_col = "vaf_abs",
  cat_levels = abs_levels,
  scheme = "absolute",
  file_stub = "Xenium_KRAS_local_VAF_CMETS",
  xlab = "Local KRAS VAF category",
  title_extra = "",
  out_dir = cmets_bin_dir
)

analyze_cmets_bins(
  cat_df,
  cat_col = "vaf_qtl",
  cat_levels = qtl_levels,
  scheme = "quartile",
  file_stub = "Xenium_KRAS_local_VAF_CMETS_quartile",
  xlab = "Local KRAS VAF quartile (within sample x radius)",
  title_extra = " (within-sample quartiles)",
  out_dir = cmets_bin_dir
)

logmsg("Writing compact CMETS vs VAF summary figures")
summary_r <- Sys.getenv(
  "MCRC_KRAS_VAF_SUMMARY_R",
  unset = "/diskmnt/Projects/MetNet_analysis/Colorectal/mCRC_Manuscript_Script/F_Tumor_subtypes/Xenium_KRAS_local_VAF_CMETS_summary.R"
)
Sys.setenv(MCRC_KRAS_CMETS_VAF_DIR = cmets_bin_dir)
sys.source(summary_r, envir = new.env())

logmsg("Done. Local VAF maps in", output_dir)
logmsg("Per-sample spatial maps in", spatial_dir)
logmsg("CMETS vs VAF-bin plots in", cmets_bin_dir)
invisible(NULL)
