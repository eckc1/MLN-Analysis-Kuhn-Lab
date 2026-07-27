library(harmony)
library(patchwork)
library(dplyr)
library(Seurat)
library(dplyr)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(purrr)
library(tibble)
library(Matrix)

# Ignore local RAM throttle 
#mem.maxVSize(vsize = Inf)

# Load merged data
#dat <- readRDS("dat_clusters_called.rds")
#dat <- readRDS("dat_clean_strict.rds")
#dat <- readRDS("dat_clean_marker_filtered_celltypes.rds")

dat <- readRDS("dat_clean_strict_adjusted_dc_cuttoff.rds")
head(dat)
Reductions(dat)
length(unique(dat$seurat_clusters))
table(dat@meta.data$marker_pass_strict)
table(dat$celltype_broad, dat$marker_pass_strict)
# Subset to only DC cells from your broad annotation
dc_dat <- subset(dat, subset = celltype_broad == "DC")
length(rownames(dc_dat@meta.data))
dc_dat <- subset(dat, subset = celltype_broad == "DC" & marker_pass_strict == TRUE)

# Plot DC cells only
dc_umap <- DimPlot(
  dc_dat,
  reduction = "wnn.umap",
  group.by = "pool_id",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Dendritic Cell clusters only")
ggsave(plot = dc_umap, filename="dc_umap.png")
table(dc_dat$pool_id)

table(dc_dat$pool_id)
table(dc_dat$orig.ident, dc_dat$pool_id)
table(dc_dat$demux_id, dc_dat$pool_id)


# If a WNN graph exists after subsetting, recluster directly on it
dc_dat <- FindClusters(
  dc_dat,
  graph.name = "wsnn",
  resolution = 0.2,
  cluster.name = "dc_wnn_clusters"
)

dc_umap_reclustered <- DimPlot(
  dc_dat,
  reduction = "wnn.umap",
  group.by = "dc_wnn_clusters",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("DC cells reclustered using WNN graph")

dc_umap_reclustered

table(dc_dat$dc_wnn_clusters)

# Make sure cluster IDs are character
dc_dat$dc_wnn_clusters <- as.character(dc_dat$dc_wnn_clusters)

# Remove dc_wnn_cluster 0
#dc_dat_no2 <- subset(
#  dc_dat,
 # subset = dc_wnn_clusters != "0"
#)

# Check result
#table(dc_dat_no2$dc_wnn_clusters)

#dc_umap_no2 <- DimPlot(
  #dc_dat_no2,
  #reduction = "wnn.umap",
  #group.by = "dc_wnn_clusters",
  #label = TRUE,
  #repel = TRUE
#) +
  #ggtitle("DC cells reclustered using WNN graph, cluster 2 removed")

#dc_umap_no2
#dc_dat <- dc_dat_no2


##################
# Run DEG Analysis

DefaultAssay(dc_dat) <- "RNA"

outdir <- "DC_Pseudobulk_Volcano_By_demux_id_unpaired"

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

# pool_id = treatment / condition
# demux_id  = mouse replicate within treatment
group_col <- "pool_id"
replicate_col <- "demux_id"

fdr_cutoff <- 0.05
logfc_cutoff <- 0.5
pvalue_cutoff <- 0.05


# Check metadata
if (!group_col %in% colnames(dc_dat@meta.data)) {
  stop("Column not found: ", group_col)
}

if (!replicate_col %in% colnames(dc_dat@meta.data)) {
  stop("Column not found: ", replicate_col)
}

cat("\nCells per treatment and mouse-within-treatment:\n")
print(table(dc_dat@meta.data[[group_col]], dc_dat@meta.data[[replicate_col]]))


# Create pseudobulk IDs, which represent unique mouse-level pseudobulk samples.
# Because pool_id repeats across treatments, the unique mouse sample identifier is: treatment + pool_id

dc_dat$treatment <- as.character(dc_dat@meta.data[[group_col]])
dc_dat$mouse_within_treatment <- as.character(dc_dat@meta.data[[replicate_col]])

dc_dat$pseudobulk_id <- paste0(
  "treat.", dc_dat$treatment,
  ".mouse.", dc_dat$mouse_within_treatment
)

# Build pseudobulk metadata
pb_meta <- dc_dat@meta.data %>%
  distinct(
    pseudobulk_id,
    treatment,
    mouse_within_treatment
  ) %>%
  arrange(treatment, mouse_within_treatment)

cat("\nExpected pseudobulk samples:\n")
print(pb_meta)

cat("\nExpected number of mouse pseudobulks per treatment:\n")
print(pb_meta %>% count(treatment, name = "n_mouse_replicates"))

# Should show 3 per treatment
cat("\nPseudobulk replicate table:\n")
print(table(pb_meta$treatment))


# Aggregate raw counts
counts <- GetAssayData(
  object = dc_dat,
  assay = "RNA",
  layer = "counts"
)

cell_meta <- dc_dat@meta.data[colnames(counts), , drop = FALSE]

# Make sure pseudobulk_id order matches cells
pseudobulk_factor <- factor(
  cell_meta$pseudobulk_id,
  levels = pb_meta$pseudobulk_id
)

# Sparse cell x pseudobulk design matrix
pb_design <- sparse.model.matrix(~ 0 + pseudobulk_factor)

colnames(pb_design) <- levels(pseudobulk_factor)

# genes x pseudobulk samples
pb_counts <- counts %*% pb_design

# Convert to sparse matrix
pb_counts <- as(pb_counts, "dgCMatrix")

# Match metadata to count columns
pb_meta <- as.data.frame(pb_meta)
rownames(pb_meta) <- pb_meta$pseudobulk_id
pb_meta <- pb_meta[colnames(pb_counts), , drop = FALSE]

cat("\nFinal pseudobulk treatment table:\n")
print(table(pb_meta$treatment))

cat("\nFinal pseudobulk treatment by mouse-within-treatment table:\n")
print(table(pb_meta$treatment, pb_meta$mouse_within_treatment))

# Save metadata and counts
write.csv(
  pb_meta,
  file = file.path(outdir, "DC_pseudobulk_metadata.csv"),
  row.names = TRUE
)

write.csv(
  as.data.frame(as.matrix(pb_counts)) %>% rownames_to_column("gene"),
  file = file.path(outdir, "DC_pseudobulk_counts.csv"),
  row.names = FALSE
)


# Run edgeR 
dge <- DGEList(counts = pb_counts)

keep <- filterByExpr(
  dge,
  group = pb_meta$treatment
)

dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

cat("\nGenes kept after filtering:\n")
print(nrow(dge))


# Pairwise unpaired edgeR comparisons

treatment_ids <- sort(unique(pb_meta$treatment))
comparisons <- combn(treatment_ids, 2, simplify = FALSE)

run_pairwise_edgeR <- function(comp) {
  
  group1 <- comp[1]
  group2 <- comp[2]
  
  comparison_name <- paste0(group2, "_vs_", group1)
  
  message("\nRunning comparison: ", comparison_name)
  
  samples_keep <- pb_meta$treatment %in% c(group1, group2)
  
  dge_sub <- dge[, samples_keep]
  meta_sub <- pb_meta[samples_keep, , drop = FALSE]
  
  meta_sub$treatment <- factor(meta_sub$treatment, levels = c(group1, group2))
  
  cat("\nSamples used for ", comparison_name, ":\n", sep = "")
  print(meta_sub)
  
  cat("\nReplicates per treatment:\n")
  print(table(meta_sub$treatment))
  
  if (any(table(meta_sub$treatment) < 2)) {
    stop("Fewer than 2 pseudobulk replicates in at least one group for ", comparison_name)
  }
  
  # Design for 9 total mice, 3 independent mice per treatment.
  design <- model.matrix(~ treatment, data = meta_sub)
  
  cat("\nDesign matrix:\n")
  print(design)
  
  dge_sub <- estimateDisp(dge_sub, design)
  fit <- glmQLFit(dge_sub, design)
  
  coef_name <- paste0("treatment", group2)
  
  if (!coef_name %in% colnames(design)) {
    stop(
      "Could not find coefficient: ", coef_name,
      "\nAvailable coefficients are:\n",
      paste(colnames(design), collapse = ", ")
    )
  }
  qlf <- glmQLFTest(fit, coef = coef_name)
  
  res <- topTags(qlf, n = Inf)$table %>%
    rownames_to_column("gene") %>%
    mutate(
      comparison = comparison_name,
      group1 = group1,
      group2 = group2
    )
  
  return(res)
}

all_results <- map_dfr(comparisons, run_pairwise_edgeR)

write.csv(
  all_results,
  file = file.path(outdir, "DC_pseudobulk_pairwise_demux_id_DE_results.csv"),
  row.names = FALSE
)


# Volcano plotting
plot_volcano <- function(res_df, comparison_name) {
  df <- res_df %>%
    filter(comparison == comparison_name) %>%
    mutate(
      neg_log10_FDR = -log10(FDR),
      sig = case_when(
        FDR < fdr_cutoff & logFC >= logfc_cutoff  ~ "Upregulated",
        FDR < fdr_cutoff & logFC <= -logfc_cutoff ~ "Downregulated",
        TRUE ~ "Not significant"
      )
    )
  
  top_up <- df %>%
    filter(sig == "Upregulated") %>%
    arrange(FDR, desc(logFC)) %>%
    slice_head(n = 10)
  
  top_down <- df %>%
    filter(sig == "Downregulated") %>%
    arrange(FDR, logFC) %>%
    slice_head(n = 10)
  
  top_labels <- bind_rows(top_up, top_down)
  
  p <- ggplot(df, aes(x = logFC, y = neg_log10_FDR, color = sig)) +
    geom_point(alpha = 0.75, size = 1.5) +
    geom_vline(
      xintercept = c(-logfc_cutoff, logfc_cutoff),
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = -log10(fdr_cutoff),
      linetype = "dashed"
    ) +
    ggrepel::geom_text_repel(
      data = top_labels,
      aes(label = gene),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.4,
      point.padding = 0.3
    ) +
    scale_color_manual(
      values = c(
        "Upregulated" = "red",
        "Downregulated" = "blue",
        "Not significant" = "gray70"
      )
    ) +
    theme_classic() +
    labs(
      title = paste0("Pseudobulk volcano: ", comparison_name),
      subtitle = paste0(
        "Significant: FDR < ", fdr_cutoff,
        " and |log2FC| >= ", logfc_cutoff
      ),
      x = "log2 fold-change",
      y = "-log10(FDR)",
      color = "Result"
    )
  
  return(p)
}

comparison_names <- unique(all_results$comparison)

for (comp_name in comparison_names) {
  
  p <- plot_volcano(all_results, comp_name)
  
  ggsave(
    filename = file.path(outdir, paste0("Volcano_", comp_name, ".png")),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
}


# Save significant and labeled genes
sig_results <- all_results %>%
  mutate(
    sig = case_when(
      FDR < fdr_cutoff & logFC >= logfc_cutoff  ~ "Upregulated",
      FDR < fdr_cutoff & logFC <= -logfc_cutoff ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  ) %>%
  filter(sig != "Not significant")

write.csv(
  sig_results,
  file = file.path(outdir, "DC_significant_DE_genes_FDR0.05_logFC1.csv"),
  row.names = FALSE
)

top_labeled_genes <- all_results %>%
  mutate(
    sig = case_when(
      FDR < fdr_cutoff & logFC >= logfc_cutoff  ~ "Upregulated",
      FDR < fdr_cutoff & logFC <= -logfc_cutoff ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  ) %>%
  filter(sig != "Not significant") %>%
  group_by(comparison, sig) %>%
  arrange(FDR, desc(abs(logFC)), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

write.csv(
  top_labeled_genes,
  file = file.path(outdir, "DC_top10_up_top10_down_labeled_genes.csv"),
  row.names = FALSE
)

# Save session info
sink(file.path(outdir, "sessionInfo.txt"))
sessionInfo()
sink()

cat("\nDone. Results saved to:\n")
cat(outdir, "\n")


pvalue_cutoff <- 0.05

# Raw p-value volcano plotting

plot_volcano_raw_p <- function(res_df, comparison_name) {
  
  df <- res_df %>%
    filter(comparison == comparison_name) %>%
    mutate(
      neg_log10_PValue = -log10(PValue),
      sig = case_when(
        PValue < pvalue_cutoff & logFC >= logfc_cutoff  ~ "Upregulated",
        PValue < pvalue_cutoff & logFC <= -logfc_cutoff ~ "Downregulated",
        TRUE ~ "Not significant"
      )
    )
  
  top_up <- df %>%
    filter(sig == "Upregulated") %>%
    arrange(PValue, desc(logFC)) %>%
    slice_head(n = 10)
  
  top_down <- df %>%
    filter(sig == "Downregulated") %>%
    arrange(PValue, logFC) %>%
    slice_head(n = 10)
  
  top_labels <- bind_rows(top_up, top_down)
  
  p <- ggplot(df, aes(x = logFC, y = neg_log10_PValue, color = sig)) +
    geom_point(alpha = 0.75, size = 1.5) +
    geom_vline(
      xintercept = c(-logfc_cutoff, logfc_cutoff),
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = -log10(pvalue_cutoff),
      linetype = "dashed"
    ) +
    ggrepel::geom_text_repel(
      data = top_labels,
      aes(label = gene),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.4,
      point.padding = 0.3
    ) +
    scale_color_manual(
      values = c(
        "Upregulated" = "red",
        "Downregulated" = "blue",
        "Not significant" = "gray70"
      )
    ) +
    scale_x_continuous(
      breaks = seq(
        floor(min(df$logFC, na.rm = TRUE)),
        ceiling(max(df$logFC, na.rm = TRUE)),
        by = 0.5
      )
    ) +
    theme_classic() +
    labs(
      title = paste0("Pseudobulk volcano using raw p-value: ", comparison_name),
      subtitle = paste0(
        "Significant: raw P < ", pvalue_cutoff,
        " and abs(log2FC) >= ", logfc_cutoff
      ),
      x = "log2 fold-change",
      y = "-log10(raw P-value)",
      color = "Result"
    )
  
  return(p)
}


exclude_genes <- c(
  "Trac", "Trbc1", "Trbc2", "Cd3d", "Cd3e", "Cd3g",
  "Cd79a", "Cd79b", "Ms4a1", "Ighm", "Ighd", "Jchain",
  "Iglc1", "Iglc2", "Iglc3"
)

all_results_no_ribo_no_contam <- all_results %>%
  filter(
    !grepl("^Rpl|^Rps", gene),
    !gene %in% exclude_genes
  )

for (comp_name in comparison_names) {
  
  p_raw_filtered <- plot_volcano_raw_p(all_results_no_ribo_no_contam, comp_name)
  
  ggsave(
    filename = file.path(outdir, paste0("Volcano_rawP_noRibo_noContam_", comp_name, ".png")),
    plot = p_raw_filtered,
    width = 8,
    height = 6,
    dpi = 300
  )
}



# DEG diagnostic summary

deg_summary <- all_results %>%
  group_by(comparison) %>%
  summarise(
    n_genes_tested = n(),
    min_PValue = min(PValue, na.rm = TRUE),
    min_FDR = min(FDR, na.rm = TRUE),
    n_rawP_0.05 = sum(PValue < 0.05, na.rm = TRUE),
    n_rawP_0.05_logFC_0.5 = sum(PValue < 0.05 & abs(logFC) >= 0.5, na.rm = TRUE),
    n_FDR_0.05 = sum(FDR < 0.05, na.rm = TRUE),
    n_FDR_0.10 = sum(FDR < 0.10, na.rm = TRUE),
    n_FDR_0.25 = sum(FDR < 0.25, na.rm = TRUE),
    .groups = "drop"
  )

print(deg_summary)

write.csv(
  deg_summary,
  file = file.path(outdir, "DC_DEG_diagnostic_summary.csv"),
  row.names = FALSE
)

# Save top-ranked genes for exploratory interpretation

top_genes_exploratory <- all_results %>%
  group_by(comparison) %>%
  arrange(PValue, .by_group = TRUE) %>%
  slice_head(n = 100) %>%
  ungroup()

write.csv(
  top_genes_exploratory,
  file = file.path(outdir, "DC_top100_genes_by_rawP_each_comparison.csv"),
  row.names = FALSE
)

# Plot number of DC cells per pseudobulk
cell_count_summary <- dc_dat@meta.data %>%
  count(pool_id, demux_id, pseudobulk_id, name = "n_cells")

p_cells <- ggplot(
  cell_count_summary,
  aes(x = pseudobulk_id, y = n_cells, fill = pool_id)
) +
  geom_col() +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Number of DC cells per pseudobulk sample",
    x = "Pseudobulk sample",
    y = "Number of DC cells",
    fill = "Treatment / pool_id"
  )

ggsave(
  filename = file.path(outdir, "DC_cells_per_pseudobulk_barplot.png"),
  plot = p_cells,
  width = 8,
  height = 5,
  dpi = 300
)

FeaturePlot(dc_dat, features = c("Itgax", "H2-Ab1", "Tcf4", "Nkg7", "Klrc1", "Cd3e", "Trac"), reduction = "wnn.umap")


