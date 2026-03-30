# ==============================================================================
# BayesPrism deconvolution pipeline
# Author: Gabriela Rapozo
# Last update: 2026-03-30
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Package setup
# ------------------------------------------------------------------------------

# Custom library path (HPC/local environment)
lib <- "~/lbbc_members/scRNASeq/lib/"
.libPaths(lib)

suppressPackageStartupMessages({
  library(BayesPrism)     # Deconvolution framework
  library(Seurat)         # Single-cell object handling
  library(dplyr)          # Data manipulation
  library(tibble)         # Data frame utilities
  library(EnsDb.Hsapiens.v86)  # Gene annotation (Ensembl v86)
})

# ------------------------------------------------------------------------------
# 2. Working directory
# ------------------------------------------------------------------------------

setwd("~/czi_macroph_adip/data/deconv/")

# ------------------------------------------------------------------------------
# 3. Bulk RNA-seq input (mixture)
# ------------------------------------------------------------------------------

# Load summarized experiment (Salmon output)
salmon.merged.gene_counts.SummarizedExperiment <- readRDS(
  "~/czi_macroph_adip/results/RNAseq_20260316/OUT_RNASEQ_20260316/star_salmon/salmon.merged.gene_counts.SummarizedExperiment.rds"
)

# Gene ID mapping file (tx2gene)
gene_id <- read.delim(
  "~/czi_macroph_adip/results/RNAseq_20260316/OUT_RNASEQ_20260316/star_salmon/tx2gene.tsv"
)

# Extract gene-level counts and filter lowly expressed genes
counts <- round(
  salmon.merged.gene_counts.SummarizedExperiment@assays@data@listData[["salmon.merged.gene_counts"]]
)

counts <- counts[rowSums(counts[, -1]) > 10, ]

# ------------------------------------------------------------------------------
# 4. Gene annotation and ID conversion (Ensembl → gene symbol)
# ------------------------------------------------------------------------------

anno <- genes(EnsDb.Hsapiens.v86, return.type = "data.frame")
anno <- anno[, c("gene_id", "gene_name", "gene_biotype")]

# Keep only protein-coding genes
anno_pc <- anno %>% dplyr::filter(gene_biotype == "protein_coding")

# Create mapping (Ensembl ID → gene symbol)
mapa <- setNames(anno_pc$gene_name, anno_pc$gene_id)

ensg <- rownames(counts)
new_rownames <- ifelse(ensg %in% names(mapa), mapa[ensg], ensg)

# Assign gene symbols and collapse duplicates by highest expression
counts$gene <- new_rownames
counts$sum_expr <- rowSums(counts[, -ncol(counts)])

counts <- counts %>%
  group_by(gene) %>%
  slice_max(order_by = sum_expr, n = 1, with_ties = FALSE) %>%
  ungroup()

counts <- counts %>% dplyr::select(-sum_expr)
counts <- counts %>% column_to_rownames("gene")

# Final bulk matrix (genes x samples)
bk.data <- counts

# Transpose to BayesPrism format (samples x genes)
bk.data <- t(bk.data)
colnames(bk.data) <- rownames(counts)
rownames(bk.data) <- colnames(counts)

# ------------------------------------------------------------------------------
# 5. Single-cell reference input
# ------------------------------------------------------------------------------

# Load Seurat object
sc.data <- qs::qread("srt_deconvolution.qs")
sc.data = subset(sc.data, subset = cell_state7 != "Myeloid")

# Extract cell type and state annotations
type_labels  <- sc.data$cell_state2
state_labels <- sc.data$cell_state7

# Extract raw counts (cells x genes)
mtx <- sc.data@assays$RNA_full@counts %>%
  as.matrix() %>%
  t()

# ------------------------------------------------------------------------------
# 6. Gene filtering and preprocessing
# ------------------------------------------------------------------------------

# Remove unwanted gene groups and lowly expressed genes
sc.data.subset <- cleanup.genes(
  input = mtx,
  input.type = "count.matrix",
  species = "hs",
  gene.group = c("Rb", "Mrp", "other_Rb", "chrM", "MALAT1", "chrX", "chrY"),
  exp.cells = 5
)

# Restrict to protein-coding genes
sc.data.subset <- select.gene.type(
  sc.data.subset,
  gene.type = "protein_coding"
)

# ------------------------------------------------------------------------------
# 7. Differential expression statistics (marker selection)
# ------------------------------------------------------------------------------

diff.exp.stat <- get.exp.stat(
  sc.dat = sc.data.subset[, colSums(sc.data.subset > 0) > 3],
  cell.type.labels = type_labels,
  cell.state.labels = state_labels,
  pseudo.count = 0.1,
  cell.count.cutoff = 50,
  n.cores = 50
)

# Select marker genes based on statistical thresholds
sc.data.subset <- select.marker(
  sc.dat = sc.data.subset,
  stat = diff.exp.stat,
  pval.max = 0.01,
  lfc.min = 0.1
)

# ------------------------------------------------------------------------------
# 8. BayesPrism model initialization
# ------------------------------------------------------------------------------

myPrism <- new.prism(
  reference = sc.data.subset,
  mixture = bk.data,
  input.type = "count.matrix",
  cell.type.labels = type_labels,
  cell.state.labels = state_labels,
  key = NULL,                  # No key cell type specified
  outlier.cut = 0.01,
  outlier.fraction = 0.1
)

# ------------------------------------------------------------------------------
# 9. Deconvolution
# ------------------------------------------------------------------------------

bp.res <- run.prism(
  prism = myPrism,
  n.cores = 50
)

# ------------------------------------------------------------------------------
# 10. Output: cell state fractions
# ------------------------------------------------------------------------------

theta_state <- get.fraction(
  bp = bp.res,
  which.theta = "first",
  state.or.type = "state"
) %>% as.data.frame()

# Save results
saveRDS(
  theta_state,
  "~/czi_macroph_adip/data/deconv/theta_state_withoutmye.RDS"
)
