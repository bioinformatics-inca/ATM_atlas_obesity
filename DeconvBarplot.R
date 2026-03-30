# ==============================================================================
# BayesPrism visualization and clustering pipeline
# Description: Hierarchical clustering and composition plots of cell-state fractions
# Author: Gabriela Rapozo
# Last update: 2026-03-30
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Package setup
# ------------------------------------------------------------------------------

# Custom library path
lib <- "~/lbbc_members/scRNASeq/lib/"
.libPaths(lib)

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
  library(dplyr)
  library(tibble)
  library(dendextend)
  library(RColorBrewer)
  library(ggdendro)
  library(cowplot)
  library(patchwork)
  library(readr)
})

# ------------------------------------------------------------------------------
# 2. Plot settings and color palettes
# ------------------------------------------------------------------------------

size.axis.Y <- 12
dend.Top <- 10
dend.Bot <- 0
height <- FALSE

# Base theme
tema <- list(
  theme_bw(),
  theme(
    axis.title = element_blank(),
    axis.text.y = element_text(size = size.axis.Y, vjust = 1),
    axis.text.x = element_text(angle = 90),
    legend.position = "none",
    plot.margin = unit(c(0, 0, 0, 0), "mm")
  ),
  scale_x_discrete(expand = c(0, 0))
)

# Cell-state palette
pal <- c(
  "ABCA1+ PVM" = "#1b9e77",
  "CD83+ PVM" = "#33a02c",
  "FGFR1+ CXCL12+ PVM" = "#66c2a5",
  "MARCO+ PVM" = "#238b45",
  "LAM" = "#41ab5d",
  "Myeloid" = "#1f78b4",
  "Early macrophages" = "#4ea3d9",
  "Mast" = "#08519c",
  "NK/T" = "#9ecae1",
  "Endothelial" = "#8c510a",
  "Lymphatic endothelial" = "#d8b365",
  "Mesothelial" = "#cab2d6",
  "Pericyte" = "#6a3d9a",
  "Adipocyte" = "#ffd92f",
  "ASPC" = "#f16913"
)

# Annotation palettes
pal_depot <- c(
  "Visceral" = "#FDC086",
  "Subcutaneous" = "#FFFF99"
)

pal_t2d <- c(
  "T2D" = "#FBB4AE",
  "Non-diabetic" = "#B3CDE3"
)

# ------------------------------------------------------------------------------
# 3. Data loading
# ------------------------------------------------------------------------------

# Load BayesPrism output (cell-state fractions)
theta <- readRDS("~/czi_macroph_adip/data/deconv/theta_state_withoutmye.RDS")
colnames(theta) <- gsub("\\.", " ", colnames(theta))
estimate <- as.data.frame(theta)

# Load sample metadata
samples_rna_new <- read_csv(
  "~/czi_macroph_adip/data/deconv/BULK RNA SEQ CZI - INCA- SEQ 032026.xlsx - Página1.csv"
)

# Standardize sample IDs
samples_rna_new$ID <- paste0("mori-rna-", sprintf("%02d", samples_rna_new$ID))

# Recode depot labels
samples_rna_new$Depot <- ifelse(
  samples_rna_new$Depot == "SC", "Subcutaneous",
  ifelse(samples_rna_new$Depot == "VC", "Visceral", samples_rna_new$Depot)
)

# Recode diabetes status
samples_rna_new$T2D <- ifelse(
  samples_rna_new$T2D == "NO", "Non-diabetic",
  ifelse(samples_rna_new$T2D == "YES", "T2D", samples_rna_new$T2D)
)

# ------------------------------------------------------------------------------
# 4. Preprocessing and normalization
# ------------------------------------------------------------------------------

bayes.obj <- estimate
bayes.obj$sample <- rownames(theta)
sample.column <- "sample"

# Keep only columns defined in palette
bayes.obj <- bayes.obj %>%
  dplyr::select(all_of(sample.column), any_of(names(pal)))

# Identify numeric columns
value.columns <- colnames(bayes.obj)[!colnames(bayes.obj) %in% sample.column]

# Normalize to percentage
sum_ciber <- rowSums(bayes.obj[, value.columns])
bayes.obj[, value.columns] <- (bayes.obj[, value.columns] / sum_ciber) * 100

# Add metadata
bayes.obj <- bayes.obj %>%
  left_join(
    samples_rna_new %>% dplyr::select(ID, Depot, T2D),
    by = c("sample" = "ID")
  )

# ------------------------------------------------------------------------------
# 5. Hierarchical clustering
# ------------------------------------------------------------------------------

hc_complete <- hclust(
  as.dist(1 - cor(t(bayes.obj[, value.columns]), method = "pearson")),
  method = "ward.D2"
)

ord <- hc_complete$order

# Reorder samples
bayes.obj[[sample.column]] <- factor(
  bayes.obj[[sample.column]],
  levels = bayes.obj[[sample.column]][ord]
)

# Build dendrogram
dendrogram <- as.dendrogram(hc_complete)
ddata <- ggdendro::dendro_data(dendrogram, type = "rectangle")

scaleFUN.dend <- function(x) sprintf("%.2f", x)

dend <- ggplot(ggdendro::segment(ddata)) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = size.axis.Y, vjust = -0.2),
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    plot.margin = unit(c(dend.Top, 0.5, dend.Bot, 0.5), "mm")
  ) +
  scale_y_continuous(labels = scaleFUN.dend, expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0))

# Hide height labels if requested
if (!height) {
  dend <- dend + theme(
    axis.text.y = element_text(color = "white"),
    axis.ticks.y = element_blank()
  )
}

# ------------------------------------------------------------------------------
# 6. Annotation tracks
# ------------------------------------------------------------------------------

plot_ann_label <- function(data, ann_label, custom_pal) {
  data[[sample.column]] <- factor(
    data[[sample.column]],
    levels = levels(bayes.obj[[sample.column]])
  )
  
  ggplot(data, aes_string(x = sample.column, y = 1, fill = ann_label)) +
    geom_tile() +
    tema +
    theme(
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_blank()
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = custom_pal)
}

# Depot annotation
p.ann_depot <- plot_ann_label(bayes.obj, "Depot", pal_depot)

legend_depot <- cowplot::get_legend(
  p.ann_depot +
    theme(
      legend.text = element_text(size = 10),
      legend.position = "right",
      legend.key.size = unit(0.8, "line")
    ) +
    guides(fill = guide_legend(ncol = 1, title = "Depot"))
)

# Diabetes annotation
p.ann_t2d <- plot_ann_label(bayes.obj, "T2D", pal_t2d)

legend_t2d <- cowplot::get_legend(
  p.ann_t2d +
    theme(
      legend.text = element_text(size = 10),
      legend.position = "right",
      legend.key.size = unit(0.8, "line")
    ) +
    guides(fill = guide_legend(ncol = 1, title = "Diabetes"))
)

# ------------------------------------------------------------------------------
# 7. Stacked barplot
# ------------------------------------------------------------------------------

d <- bayes.obj %>%
  dplyr::select(-Depot, -T2D) %>%
  reshape2::melt(id.vars = sample.column)

d[[sample.column]] <- factor(
  d[[sample.column]],
  levels = levels(bayes.obj[[sample.column]])
)

p1 <- suppressWarnings(
  ggplot(d, aes_string(x = sample.column, y = "value", fill = "variable")) +
    geom_bar(stat = "identity", width = 1) +
    labs(fill = "Subpopulation", y = "Relative percent") +
    tema +
    scale_y_continuous(expand = c(0, 0)) +
    guides(fill = guide_legend(ncol = 2)) +
    scale_fill_manual(values = pal)
)

legend_main <- cowplot::get_legend(
  p1 +
    guides(fill = guide_legend(ncol = 1)) +
    theme(
      legend.text = element_text(size = 8),
      legend.position = "right"
    )
)

# ------------------------------------------------------------------------------
# 8. Final figure assembly
# ------------------------------------------------------------------------------

legends <- cowplot::plot_grid(
  legend_main,
  legend_depot,
  legend_t2d,
  ncol = 1,
  rel_heights = c(2, 1, 1)
)

plots <- cowplot::plot_grid(
  dend,
  p.ann_depot + theme(legend.position = "none"),
  p.ann_t2d + theme(legend.position = "none"),
  p1 + theme(legend.position = "none"),
  ncol = 1,
  align = "v",
  rel_heights = c(1.8, 0.45, 0.45, 8)
)

final_plot <- cowplot::plot_grid(
  plots,
  legends,
  ncol = 2,
  rel_widths = c(8, 2)
)

print(final_plot)