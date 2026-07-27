library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)
library(dplyr)
library(edgeR)
library(ggrepel)
library(purrr)
library(tibble)
library(Matrix)
library(tidyr)


# Ignore local RAM throttle 
mem.maxVSize(vsize = Inf)
# Load merged data
#dat <- readRDS("dat_clusters_called.rds")
dat <- readRDS("dat_clean_strict.rds")
head(dat)
Reductions(dat)
length(unique(dat$seurat_clusters))
colnames(dat@meta.data)
dat@meta.data$celltype_broad
DimPlot(dat, reduction = "wnn.umap", group.by = "celltype_broad")
# Subset to only CD4_T cells from your broad annotation
CD4_T_dat <- subset(dat, subset = celltype_broad == "CD4_T")
colnames(CD4_T_dat@meta.data)

# Plot CD4_T cells only
CD4_T_umap <- DimPlot(
  CD4_T_dat,
  reduction = "wnn.umap",
  group.by = "pool_id",
  #group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("CD4 T Cell clusters only")
#ggsave(plot = CD4_T_umap, filename="CD4_T_umap.png")
#ggsave(plot = CD4_T_umap, filename="CD4_T_umap_transcript.png")
class(CD4_T_dat)
dim(CD4_T_dat)

table(CD4_T_dat@meta.data$seurat_clusters)

Reductions(CD4_T_dat)

colnames(CD4_T_dat@meta.data)
table(CD4_T_dat$pool_id)
table(CD4_T_dat$orig.ident, CD4_T_dat$pool_id)
table(CD4_T_dat$demux_id, CD4_T_dat$pool_id)

# Reculster CD4 T cell object

DefaultAssay(CD4_T_dat) <- "RNA"
CD4_T_dat <- NormalizeData(CD4_T_dat)
CD4_T_dat <- FindVariableFeatures(CD4_T_dat)
CD4_T_dat <- ScaleData(CD4_T_dat)
CD4_T_dat <- RunPCA(CD4_T_dat, reduction.name = "pca")

DefaultAssay(CD4_T_dat) <- "ADT"
CD4_T_dat <- NormalizeData(CD4_T_dat, normalization.method = "CLR", margin = 2)
CD4_T_dat <- ScaleData(CD4_T_dat)
CD4_T_dat <- RunPCA(CD4_T_dat, reduction.name = "apca")

CD4_T_dat <- FindMultiModalNeighbors(
  CD4_T_dat,
  reduction.list = list("pca", "apca"),
  dims.list = list(1:20, 1:18)
)

CD4_T_dat <- FindClusters(
  CD4_T_dat,
  graph.name = "wsnn",
  resolution = 0.7
)

CD4_T_dat <- RunUMAP(
  CD4_T_dat,
  nn.name = "weighted.nn",
  reduction.name = "cd4.wnn.umap",
  reduction.key = "cd4wnnUMAP_"
)

CD4_T_umap <- DimPlot(
  CD4_T_dat,
  reduction = "wnn.umap",
  group.by = "pool_id",
  label = TRUE,repel = TRUE
) +
  ggtitle("CD4 T Cell clusters only")
ggsave(plot = CD4_T_umap, filename="reclustered_CD4_T_umap.png")


CD4_T_umap <- DimPlot(
  CD4_T_dat,
  reduction = "wnn.umap",
  group.by = "seurat_clusters",
  label = TRUE,repel = TRUE
) +
  ggtitle("CD4 T Cell clusters only")
ggsave(plot = CD4_T_umap, filename="reclustered_CD4_T_umap_transcript.png")




# Keep only TCRB / CD3 / CD4 positive cells

# Set ADT as default temporarily for subsetting
#DefaultAssay(dat) <- "ADT"
#adt_counts <- GetAssayData(dat, assay = "ADT", layer = "counts")

#cells_keep <- colnames(dat)[
 #   adt_counts["Ms.TCR.Bchain", ] > 1 &
#    adt_counts["Ms.CD3", ] > 1 &
#    adt_counts["Ms.CD4", ] > 1
#]

#CD4_T_dat <- subset(CD4_T_dat, cells = cells_keep)
#dim(CD4_T_dat)

# Run DEG Analysis

DefaultAssay(CD4_T_dat) <- "RNA"
#CD4_T_dat <- subset(
  #CD4_T_dat,
  #subset =
    #Cd3d > 0 &
   # Cd3e > 0 &
    #Trac > 0 &
    #Cd4 > 0 &
    #Cd8a == 0 &
    #Cd8b1 == 0 &
    #Cd79a == 0 &
    #Ms4a1 == 0 &
    #Bank1 == 0 &
    #Lyz2 == 0 &
    #Itgax == 0 &
   # `H2-Aa` == 0 &
   # `H2-Ab1` == 0
#)
# remove obvious non-cd4 t cell markers
CD4_T_dat <- subset(
  CD4_T_dat,
  subset =
    Cd3d > 0 &
    Cd3e > 0 &
    Trac > 0 &
    Cd4 > 0 &
    Cd8a == 0 &
    Cd8b1 == 0 &
    Cd79a == 0 &
    Ms4a1 == 0 &
    Bank1 == 0 &
    Lyz2 == 0 &
    Itgax == 0
)

dim(CD4_T_dat)
outdir <- "CD4_T_Pseudobulk_Volcano_By_demux_id_unpaired"

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

# pool_id = treatment / condition
# demux_id  = mouse replicate within treatment
group_col <- "pool_id"
replicate_col <- "demux_id"

# Thresholds
fdr_cutoff <- 0.05
logfc_cutoff <- 0.5

cat("\nCells per treatment and mouse-within-treatment:\n")
print(table(CD4_T_dat@meta.data[[group_col]], CD4_T_dat@meta.data[[replicate_col]]))


# Create pseudobulk IDs, which represent unique mouse-level pseudobulk samples.
# pool_id repeats across treatments, so the true unique mouse sample is: treatment + pool_id

CD4_T_dat$treatment <- as.character(CD4_T_dat@meta.data[[group_col]])
CD4_T_dat$mouse_within_treatment <- as.character(CD4_T_dat@meta.data[[replicate_col]])

CD4_T_dat$pseudobulk_id <- paste0(
  "treat.", CD4_T_dat$treatment,
  ".mouse.", CD4_T_dat$mouse_within_treatment
)

# Build pseudobulk metadata
pb_meta <- CD4_T_dat@meta.data %>%
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


# Aggregate raw counts and rename pseudobulk IDs.

counts <- GetAssayData(
  object = CD4_T_dat,
  assay = "RNA",
  layer = "counts"
)

cell_meta <- CD4_T_dat@meta.data[colnames(counts), , drop = FALSE]

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
  file = file.path(outdir, "CD4_T_pseudobulk_metadata.csv"),
  row.names = TRUE
)

write.csv(
  as.data.frame(as.matrix(pb_counts)) %>% rownames_to_column("gene"),
  file = file.path(outdir, "CD4_T_pseudobulk_counts.csv"),
  row.names = FALSE
)


# edgeR setup
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
  
  # 3 independent mice per treatment.Do not include mouse as a blocking factor 
  # because pool_id repeats across treatments but does not represent the same mouse 
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
  file = file.path(outdir, "CD4_T_pseudobulk_pairwise_demux_id_DE_results.csv"),
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
  file = file.path(outdir, "CD4_T_significant_DE_genes_FDR0.05_logFC1.csv"),
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
  file = file.path(outdir, "CD4_T_top10_up_top10_down_labeled_genes.csv"),
  row.names = FALSE
)


# session info

sink(file.path(outdir, "sessionInfo.txt"))
sessionInfo()
sink()

cat("\nDone. Results saved to:\n")
cat(outdir, "\n")


test_features<- FeaturePlot(
  CD4_T_dat,
  features = c(
    "Cd3d", "Cd3e", "Trac", "Cd4",
    "Cd8a", "Cd8b1",
    "Cd79a", "Ms4a1", "Cd74", "H2-Aa", "H2-Ab1",
    "Nkg7", "Bank1", "Lyz2", "Itgax"
  ),
  reduction = "wnn.umap",
  ncol = 4,
  order=TRUE
)
ggsave(plot=test_features, "cd4_test_features.png")

marker_violin_cd4 <- VlnPlot(
  CD4_T_dat,
  features = c("Cd3d", "Cd3e", "Trac", "Cd4", "Cd79a", "Cd74", "H2-Aa", "H2-Ab1", "Bank1"),
  group.by = "seurat_clusters",
  pt.size = 0
)
ggsave(plot=marker_violin_cd4, "marker_violin_cd4.png")





table(CD4_T_dat$treatment)
table(CD4_T_dat$pool_id)



CD4_T_dat@meta.data %>%
  dplyr::count(pool_id, demux_id) %>%
  tidyr::pivot_wider(
    names_from = demux_id,
    values_from = n,
    values_fill = 0
  )







# Cell counts by CD4 cluster and treatment
cluster_counts <- CD4_T_dat@meta.data %>%
  dplyr::count(pool_id, demux_id, seurat_clusters)

# Proportions within each mouse/treatment pseudobulk
cluster_props <- cluster_counts %>%
  dplyr::group_by(pool_id, demux_id) %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  dplyr::ungroup()

cluster_props




plot <- ggplot(cluster_props, aes(x = pool_id, y = prop, fill = pool_id)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2) +
  facet_wrap(~ seurat_clusters, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Treatment / pool_id",
    y = "Proportion of CD4 T cells",
    title = "CD4 T-cell state composition by treatment"
  )
ggsave(plot=plot, filename="cluster_composition_boxplots.png")



### Cluster-based analysis


# CD4 T-cell cluster-based analysis using existing embeddings

# Output directory
outdir <- "CD4_T_existing_cluster_based_analysis"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)


# Use existing CD4_T object
obj <- CD4_T_dat

cat("\nAvailable reductions in CD4_T_dat:\n")
print(Reductions(obj))

cat("\nMetadata columns:\n")
print(colnames(obj@meta.data))


# use the embedding present in the object
if ("wnn.umap" %in% Reductions(obj)) {
  reduction_use <- "wnn.umap"
} else if ("umap" %in% Reductions(obj)) {
  reduction_use <- "umap"
} else if ("cd4.wnn.umap" %in% Reductions(obj)) {
  reduction_use <- "cd4.wnn.umap"
} else {
  stop(
    "No usable UMAP reduction found. Available reductions are: ",
    paste(Reductions(obj), collapse = ", ")
  )
}

cat("\nUsing existing reduction:\n")
print(reduction_use)


# Required metadata columns
cluster_col <- "seurat_clusters"
treatment_col <- "pool_id"
replicate_col <- "demux_id"

required_cols <- c(cluster_col, treatment_col, replicate_col)

missing_cols <- required_cols[!required_cols %in% colnames(obj@meta.data)]

if (length(missing_cols) > 0) {
  stop(
    "Missing required metadata columns: ",
    paste(missing_cols, collapse = ", ")
  )
}


# Standardize metadata columns
obj$cluster_id <- as.character(obj@meta.data[[cluster_col]])
obj$treatment <- as.character(obj@meta.data[[treatment_col]])
obj$replicate_id <- as.character(obj@meta.data[[replicate_col]])

obj$treatment <- factor(obj$treatment, levels = c("601", "602", "603"))

cat("\nExisting CD4 cluster counts:\n")
print(table(obj$cluster_id, useNA = "ifany"))

cat("\nTreatment counts:\n")
print(table(obj$treatment, useNA = "ifany"))

cat("\nTreatment by replicate counts:\n")
print(table(obj$treatment, obj$replicate_id, useNA = "ifany"))


# Curated cluster labels
cluster_labels <- c(
  "0"  = "Unclear resting/memory-like T cells",
  "1"  = "Ribosomal-high / low-information cells",
  "2"  = "IFN-responsive / Ly6c+ activated T cells",
  "3"  = "Cytotoxic-like / activated T cells",
  "4"  = "Naive / central memory-like CD4 T cells",
  "5"  = "Activated/metabolic T cells",
  "6"  = "Regulatory T cells / Treg cells",
  "7"  = "Ribosomal-high / low-information cells",
  "8"  = "PD-1+ activated / Tfh-like or exhausted-like T cells",
  "9"  = "Unclear IFN/stress-response T cells",
  "10" = "Strongly activated effector T cells",
  "11" = "Unclear stress-related T cell state",
  "12" = "Th17-like cells"
)

obj$curated_label <- unname(cluster_labels[obj$cluster_id])

obj$curated_label_numbered <- ifelse(
  is.na(obj$curated_label),
  paste0("cluster_", obj$cluster_id, ": Unlabeled"),
  paste0("cluster_", obj$cluster_id, ": ", obj$curated_label)
)

cat("\nCluster label check:\n")
print(table(obj$cluster_id, obj$curated_label_numbered, useNA = "ifany"))


# Save labeled object and metadata
CD4_T_dat <- obj

saveRDS(
  CD4_T_dat,
  file = file.path(outdir, "CD4_T_existing_clusters_with_curated_labels.rds")
)

write.csv(
  CD4_T_dat@meta.data,
  file = file.path(outdir, "CD4_T_existing_clusters_metadata_with_curated_labels.csv"),
  row.names = TRUE
)


# UMAP plots using existing embedding
p_existing_clusters <- DimPlot(
  obj,
  reduction = reduction_use,
  group.by = "cluster_id",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("CD4 T cells colored by existing clusters") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  filename = file.path(outdir, "UMAP_existing_CD4_clusters.png"),
  plot = p_existing_clusters,
  width = 9,
  height = 7,
  dpi = 300
)

ggsave(
  filename = file.path(outdir, "UMAP_existing_CD4_clusters.pdf"),
  plot = p_existing_clusters,
  width = 9,
  height = 7
)

print(p_existing_clusters)

p_curated <- DimPlot(
  obj,
  reduction = reduction_use,
  group.by = "curated_label_numbered",
  label = FALSE
) +
  ggtitle("CD4 T cells colored by curated cluster labels") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.text = element_text(size = 8)
  )

ggsave(
  filename = file.path(outdir, "UMAP_existing_embedding_curated_cluster_labels.png"),
  plot = p_curated,
  width = 13,
  height = 8,
  dpi = 300
)

ggsave(
  filename = file.path(outdir, "UMAP_existing_embedding_curated_cluster_labels.pdf"),
  plot = p_curated,
  width = 13,
  height = 8
)

print(p_curated)

p_treatment <- DimPlot(
  obj,
  reduction = reduction_use,
  group.by = "treatment",
  label = FALSE
) +
  ggtitle("CD4 T cells colored by treatment / pool_id") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  filename = file.path(outdir, "UMAP_existing_embedding_by_treatment.png"),
  plot = p_treatment,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  filename = file.path(outdir, "UMAP_existing_embedding_by_treatment.pdf"),
  plot = p_treatment,
  width = 8,
  height = 6
)

print(p_treatment)


# Focus clusters
clusters_of_interest <- c("2", "4", "8", "10", "12")

obj$focus_cluster <- ifelse(
  obj$cluster_id %in% clusters_of_interest,
  obj$curated_label_numbered,
  "Other CD4 T-cell clusters"
)

p_focus_umap <- DimPlot(
  obj,
  reduction = reduction_use,
  group.by = "focus_cluster",
  label = FALSE
) +
  ggtitle("Existing CD4 embedding highlighting clusters 2, 4, 8, 10, and 12") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.text = element_text(size = 8)
  )

ggsave(
  filename = file.path(outdir, "UMAP_existing_embedding_focus_clusters_2_4_8_10_12.png"),
  plot = p_focus_umap,
  width = 13,
  height = 8,
  dpi = 300
)

ggsave(
  filename = file.path(outdir, "UMAP_existing_embedding_focus_clusters_2_4_8_10_12.pdf"),
  plot = p_focus_umap,
  width = 13,
  height = 8
)

print(p_focus_umap)

# Save object with focus labels
CD4_T_dat <- obj

saveRDS(
  CD4_T_dat,
  file = file.path(outdir, "CD4_T_existing_clusters_with_focus_labels.rds")
)


# Marker validation for focus clusters
DefaultAssay(obj) <- "RNA"

marker_genes <- c(
  # cluster 2 IFN / Ly6c activation
  "Tgtp1", "Usp18", "Rtp4", "Ly6c1", "Cxcr4",
  
  # cluster 4 naive / central memory-like
  "Il7r", "Tsc22d3", "Trat1", "Txk", "Pik3ip1",
  
  # cluster 8 PD-1 / Tfh-like / exhaustion-like
  "Pdcd1", "Tigit", "Ctla4", "Maf", "Cd44", "Tox2",
  
  # cluster 10 effector activation
  "Tnfrsf9", "Nr4a1", "Nr4a3", "Relb", "Nfkb2", "Mir155hg", "Slc7a5",
  
  # cluster 12 Th17-like
  "Ccr6", "Il1r1", "Rorc", "Rora", "Il18rap", "Ccr4"
)

marker_genes_present <- marker_genes[marker_genes %in% rownames(obj)]

cat("\nMarker genes present:\n")
print(marker_genes_present)

cat("\nMarker genes missing:\n")
print(setdiff(marker_genes, marker_genes_present))

focus_dat <- subset(obj, subset = cluster_id %in% clusters_of_interest)

if (length(marker_genes_present) > 0) {
  
  p_dot <- DotPlot(
    focus_dat,
    features = marker_genes_present,
    group.by = "curated_label_numbered"
  ) +
    RotatedAxis() +
    ggtitle("Marker validation for focus CD4 T-cell clusters") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(
    filename = file.path(outdir, "DotPlot_focus_cluster_marker_validation.png"),
    plot = p_dot,
    width = 15,
    height = 7,
    dpi = 300
  )
  
  ggsave(
    filename = file.path(outdir, "DotPlot_focus_cluster_marker_validation.pdf"),
    plot = p_dot,
    width = 15,
    height = 7
  )
  
  print(p_dot)
  
  p_feature <- FeaturePlot(
    obj,
    features = marker_genes_present,
    reduction = reduction_use,
    ncol = 5,
    order = TRUE
  )
  
  ggsave(
    filename = file.path(outdir, "FeaturePlot_focus_cluster_markers.png"),
    plot = p_feature,
    width = 18,
    height = 14,
    dpi = 300
  )
  
  print(p_feature)
}

# Composition analysis across treatments
# clusters 2, 4, 8, 10, and 12
composition_counts <- obj@meta.data %>%
  count(
    treatment,
    replicate_id,
    cluster_id,
    curated_label_numbered,
    name = "n_cells"
  ) %>%
  arrange(treatment, replicate_id, cluster_id)

write.csv(
  composition_counts,
  file = file.path(outdir, "Existing_CD4_cluster_counts_by_treatment_replicate.csv"),
  row.names = FALSE
)

replicate_totals <- obj@meta.data %>%
  count(treatment, replicate_id, name = "total_cd4_cells")

composition_props <- composition_counts %>%
  left_join(replicate_totals, by = c("treatment", "replicate_id")) %>%
  mutate(prop_of_cd4 = n_cells / total_cd4_cells)

write.csv(
  composition_props,
  file = file.path(outdir, "Existing_CD4_cluster_proportions_by_treatment_replicate.csv"),
  row.names = FALSE
)

focus_composition_props <- composition_props %>%
  filter(cluster_id %in% clusters_of_interest)

write.csv(
  focus_composition_props,
  file = file.path(outdir, "Focus_existing_cluster_proportions_by_treatment_replicate.csv"),
  row.names = FALSE
)

p_focus_props <- ggplot(
  focus_composition_props,
  aes(x = treatment, y = prop_of_cd4, fill = treatment)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.12, size = 2.5) +
  facet_wrap(~ curated_label_numbered, scales = "free_y") +
  theme_classic() +
  labs(
    title = "Existing CD4 cluster composition across treatment groups",
    x = "Treatment / pool_id",
    y = "Proportion of total CD4 T cells"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

ggsave(
  filename = file.path(outdir, "Focus_existing_cluster_proportions_across_treatments.png"),
  plot = p_focus_props,
  width = 13,
  height = 8,
  dpi = 300
)

ggsave(
  filename = file.path(outdir, "Focus_existing_cluster_proportions_across_treatments.pdf"),
  plot = p_focus_props,
  width = 13,
  height = 8
)

print(p_focus_props)


# Kruskal-Wallis tests across treatments
composition_tests <- focus_composition_props %>%
  group_by(cluster_id, curated_label_numbered) %>%
  summarise(
    n_replicates = n(),
    n_treatments = n_distinct(treatment),
    kruskal_p = kruskal.test(prop_of_cd4 ~ treatment)$p.value,
    mean_601 = mean(prop_of_cd4[treatment == "601"], na.rm = TRUE),
    mean_602 = mean(prop_of_cd4[treatment == "602"], na.rm = TRUE),
    mean_603 = mean(prop_of_cd4[treatment == "603"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(kruskal_FDR = p.adjust(kruskal_p, method = "BH")) %>%
  arrange(kruskal_p)

write.csv(
  composition_tests,
  file = file.path(outdir, "Focus_existing_cluster_composition_KruskalWallis_tests.csv"),
  row.names = FALSE
)

print(composition_tests)


# Pairwise Wilcoxon tests
pairwise_composition_tests <- focus_composition_props %>%
  group_by(cluster_id, curated_label_numbered) %>%
  group_modify(~ {
    
    if (n_distinct(.x$treatment) < 2) {
      return(tibble())
    }
    
    pw <- pairwise.wilcox.test(
      x = .x$prop_of_cd4,
      g = .x$treatment,
      p.adjust.method = "BH",
      exact = FALSE
    )
    
    as.data.frame(as.table(pw$p.value)) %>%
      filter(!is.na(Freq)) %>%
      rename(
        treatment_1 = Var1,
        treatment_2 = Var2,
        pairwise_p_adj = Freq
      )
  }) %>%
  ungroup()

write.csv(
  pairwise_composition_tests,
  file = file.path(outdir, "Focus_existing_cluster_composition_pairwise_Wilcoxon_tests.csv"),
  row.names = FALSE
)

print(pairwise_composition_tests)


# Th17-like cell analysis using existing cluster 12
th17_cluster <- "12"

th17_counts <- obj@meta.data %>%
  mutate(is_Th17_cluster = cluster_id == th17_cluster) %>%
  group_by(treatment, replicate_id) %>%
  summarise(
    n_Th17 = sum(is_Th17_cluster),
    total_CD4 = n(),
    prop_Th17 = n_Th17 / total_CD4,
    .groups = "drop"
  )

write.csv(
  th17_counts,
  file = file.path(outdir, "Existing_cluster12_Th17_counts_proportions_by_treatment_replicate.csv"),
  row.names = FALSE
)

print(th17_counts)

p_th17_counts <- ggplot(
  th17_counts,
  aes(x = treatment, y = n_Th17, fill = treatment)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.12, size = 3) +
  theme_classic() +
  labs(
    title = "Number of Th17-like cells across treatment groups",
    subtitle = "Th17-like cells defined as existing CD4 cluster 12",
    x = "Treatment / pool_id",
    y = "Number of Th17-like cells"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )

ggsave(
  filename = file.path(outdir, "Existing_cluster12_Th17_cell_counts_across_treatments.png"),
  plot = p_th17_counts,
  width = 7,
  height = 6,
  dpi = 300
)

ggsave(
  filename = file.path(outdir, "Existing_cluster12_Th17_cell_counts_across_treatments.pdf"),
  plot = p_th17_counts,
  width = 7,
  height = 6
)

print(p_th17_counts)

p_th17_props <- ggplot(
  th17_counts,
  aes(x = treatment, y = prop_Th17, fill = treatment)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.12, size = 3) +
  theme_classic() +
  labs(
    title = "Proportion of Th17-like cells across treatment groups",
    subtitle = "Th17-like cells defined as existing CD4 cluster 12",
    x = "Treatment / pool_id",
    y = "Proportion of total CD4 T cells"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )

ggsave(
  filename = file.path(outdir, "Existing_cluster12_Th17_proportions_across_treatments.png"),
  plot = p_th17_props,
  width = 7,
  height = 6,
  dpi = 300
)

ggsave(
  filename = file.path(outdir, "Existing_cluster12_Th17_proportions_across_treatments.pdf"),
  plot = p_th17_props,
  width = 7,
  height = 6
)

print(p_th17_props)


# Th17 statistical tests
th17_count_test <- kruskal.test(n_Th17 ~ treatment, data = th17_counts)
th17_prop_test <- kruskal.test(prop_Th17 ~ treatment, data = th17_counts)

th17_test_results <- data.frame(
  test = c("Th17 cell count", "Th17 proportion of CD4"),
  p_value = c(th17_count_test$p.value, th17_prop_test$p.value)
) %>%
  mutate(FDR = p.adjust(p_value, method = "BH"))

write.csv(
  th17_test_results,
  file = file.path(outdir, "Existing_cluster12_Th17_KruskalWallis_tests.csv"),
  row.names = FALSE
)

print(th17_test_results)

th17_pairwise_count <- pairwise.wilcox.test(
  x = th17_counts$n_Th17,
  g = th17_counts$treatment,
  p.adjust.method = "BH",
  exact = FALSE
)

th17_pairwise_prop <- pairwise.wilcox.test(
  x = th17_counts$prop_Th17,
  g = th17_counts$treatment,
  p.adjust.method = "BH",
  exact = FALSE
)

th17_pairwise_tests <- bind_rows(
  as.data.frame(as.table(th17_pairwise_count$p.value)) %>%
    filter(!is.na(Freq)) %>%
    rename(treatment_1 = Var1, treatment_2 = Var2, p_adj = Freq) %>%
    mutate(test = "Th17 cell count"),
  
  as.data.frame(as.table(th17_pairwise_prop$p.value)) %>%
    filter(!is.na(Freq)) %>%
    rename(treatment_1 = Var1, treatment_2 = Var2, p_adj = Freq) %>%
    mutate(test = "Th17 proportion of CD4")
)

write.csv(
  th17_pairwise_tests,
  file = file.path(outdir, "Existing_cluster12_Th17_pairwise_Wilcoxon_tests.csv"),
  row.names = FALSE
)

print(th17_pairwise_tests)

# Th17 marker score across clusters and treatments
DefaultAssay(obj) <- "RNA"

th17_markers <- c("Ccr6", "Il1r1", "Rorc", "Rora", "Il18rap", "Ccr4")
th17_markers_present <- th17_markers[th17_markers %in% rownames(obj)]

cat("\nTh17 markers present:\n")
print(th17_markers_present)

cat("\nTh17 markers missing:\n")
print(setdiff(th17_markers, th17_markers_present))

if (length(th17_markers_present) > 0) {
  
  obj <- AddModuleScore(
    obj,
    features = list(th17_markers_present),
    name = "Th17_score"
  )
  
  score_col <- "Th17_score1"
  
  p_th17_score_cluster <- VlnPlot(
    obj,
    features = score_col,
    group.by = "cluster_id",
    pt.size = 0
  ) +
    ggtitle("Th17 marker score across existing CD4 clusters") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(
    filename = file.path(outdir, "Th17_marker_score_by_existing_cluster.png"),
    plot = p_th17_score_cluster,
    width = 10,
    height = 6,
    dpi = 300
  )
  
  print(p_th17_score_cluster)
  
  th17_score_summary <- obj@meta.data %>%
    group_by(treatment, replicate_id, cluster_id, curated_label_numbered) %>%
    summarise(
      mean_Th17_score = mean(.data[[score_col]], na.rm = TRUE),
      median_Th17_score = median(.data[[score_col]], na.rm = TRUE),
      n_cells = n(),
      .groups = "drop"
    )
  
  write.csv(
    th17_score_summary,
    file = file.path(outdir, "Th17_marker_score_by_treatment_replicate_cluster.csv"),
    row.names = FALSE
  )
  
  p_th17_score_treatment_cluster <- ggplot(
    th17_score_summary,
    aes(x = treatment, y = mean_Th17_score, fill = treatment)
  ) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.12, size = 2) +
    facet_wrap(~ curated_label_numbered, scales = "free_y") +
    theme_classic() +
    labs(
      title = "Mean Th17 marker score by existing cluster and treatment",
      x = "Treatment / pool_id",
      y = "Mean Th17 module score per replicate"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"
    )
  
  ggsave(
    filename = file.path(outdir, "Th17_marker_score_by_existing_cluster_and_treatment.png"),
    plot = p_th17_score_treatment_cluster,
    width = 14,
    height = 10,
    dpi = 300
  )
  
  print(p_th17_score_treatment_cluster)
}

# Save object with scores
CD4_T_dat <- obj

saveRDS(
  CD4_T_dat,
  file = file.path(outdir, "CD4_T_existing_clusters_with_labels_and_Th17_score.rds")
)


# Pseudobulk transcriptional DE within each focus cluster
# Test transcriptional composition within clusters across treatments 601, 602, and 603
DefaultAssay(obj) <- "RNA"

pb_outdir <- file.path(outdir, "Pseudobulk_DE_within_focus_clusters")
dir.create(pb_outdir, showWarnings = FALSE, recursive = TRUE)

run_cluster_pseudobulk_DE <- function(cluster_use) {
  
  message("\n==============================")
  message("Running pseudobulk DE for cluster ", cluster_use)
  message("==============================")
  
  cluster_obj <- subset(obj, subset = cluster_id == cluster_use)
  
  cluster_label <- unique(cluster_obj$curated_label_numbered)
  cluster_label <- cluster_label[1]
  cluster_label_clean <- gsub("[^A-Za-z0-9_]+", "_", cluster_label)
  
  cluster_dir <- file.path(
    pb_outdir,
    paste0("cluster_", cluster_use, "_", cluster_label_clean)
  )
  dir.create(cluster_dir, showWarnings = FALSE, recursive = TRUE)
  
  cluster_obj$pseudobulk_id <- paste0(
    "treat.", cluster_obj$treatment,
    ".mouse.", cluster_obj$replicate_id
  )
  
  pb_meta <- cluster_obj@meta.data %>%
    distinct(pseudobulk_id, treatment, replicate_id) %>%
    arrange(treatment, replicate_id)
  
  rep_table <- table(pb_meta$treatment)
  
  write.csv(
    as.data.frame(rep_table),
    file = file.path(cluster_dir, "replicates_per_treatment.csv"),
    row.names = FALSE
  )
  
  cat("\nReplicates per treatment for cluster ", cluster_use, ":\n", sep = "")
  print(rep_table)
  
  if (any(rep_table < 2)) {
    warning(
      "Skipping cluster ", cluster_use,
      ": fewer than 2 pseudobulk replicates in at least one treatment."
    )
    return(NULL)
  }
  
  counts <- GetAssayData(
    cluster_obj,
    assay = "RNA",
    layer = "counts"
  )
  
  cell_meta <- cluster_obj@meta.data[colnames(counts), , drop = FALSE]
  
  pseudobulk_factor <- factor(
    cell_meta$pseudobulk_id,
    levels = pb_meta$pseudobulk_id
  )
  
  pb_design <- sparse.model.matrix(~ 0 + pseudobulk_factor)
  colnames(pb_design) <- levels(pseudobulk_factor)
  
  pb_counts <- counts %*% pb_design
  pb_counts <- as(pb_counts, "dgCMatrix")
  
  pb_meta <- as.data.frame(pb_meta)
  rownames(pb_meta) <- pb_meta$pseudobulk_id
  pb_meta <- pb_meta[colnames(pb_counts), , drop = FALSE]
  pb_meta$treatment <- factor(pb_meta$treatment, levels = c("601", "602", "603"))
  
  write.csv(
    pb_meta,
    file = file.path(cluster_dir, "pseudobulk_metadata.csv"),
    row.names = TRUE
  )
  
  write.csv(
    as.data.frame(as.matrix(pb_counts)) %>%
      rownames_to_column("gene"),
    file = file.path(cluster_dir, "pseudobulk_counts.csv"),
    row.names = FALSE
  )
  
  dge <- DGEList(counts = pb_counts)
  
  keep <- filterByExpr(
    dge,
    group = pb_meta$treatment
  )
  
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  
  cat("\nGenes kept for cluster ", cluster_use, ": ", nrow(dge), "\n", sep = "")
  
  if (nrow(dge) < 50) {
    warning(
      "Skipping cluster ", cluster_use,
      ": fewer than 50 genes after filterByExpr."
    )
    return(NULL)
  }
  

# edgeR design

  pb_meta$treatment <- factor(
    pb_meta$treatment,
    levels = c("601", "602", "603")
  )
  
  pb_meta$treatment_safe <- factor(
    paste0("treat_", as.character(pb_meta$treatment)),
    levels = c("treat_601", "treat_602", "treat_603")
  )
  
  design <- model.matrix(~ 0 + treatment_safe, data = pb_meta)
  
  colnames(design) <- gsub("^treatment_safe", "", colnames(design))
  
  cat("\nDesign matrix for cluster ", cluster_use, ":\n", sep = "")
  print(design)
  
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  
  contrasts <- makeContrasts(
    `602_vs_601` = treat_602 - treat_601,
    `603_vs_601` = treat_603 - treat_601,
    `603_vs_602` = treat_603 - treat_602,
    levels = design
  )
  
  cluster_results <- list()
  
  for (comp in colnames(contrasts)) {
    
    qlf <- glmQLFTest(
      fit,
      contrast = contrasts[, comp]
    )
    
    res <- topTags(qlf, n = Inf)$table %>%
      rownames_to_column("gene") %>%
      mutate(
        cluster_id = cluster_use,
        cluster_label = cluster_label,
        comparison = comp
      )
    
    write.csv(
      res,
      file = file.path(cluster_dir, paste0("DE_", comp, ".csv")),
      row.names = FALSE
    )
    
    cluster_results[[comp]] <- res
  }
  
  bind_rows(cluster_results)
}

all_cluster_DE <- map_dfr(
  clusters_of_interest,
  run_cluster_pseudobulk_DE
)

write.csv(
  all_cluster_DE,
  file = file.path(pb_outdir, "All_focus_cluster_pseudobulk_DE_results.csv"),
  row.names = FALSE
)

sig_cluster_DE <- all_cluster_DE %>%
  filter(FDR < 0.05, abs(logFC) >= 0.5)

write.csv(
  sig_cluster_DE,
  file = file.path(pb_outdir, "Significant_focus_cluster_pseudobulk_DE_FDR0.05_logFC0.5.csv"),
  row.names = FALSE
)


# Volcano plots for pseudobulk DE
plot_cluster_volcano <- function(res_df, cluster_use, comparison_use) {
  
  df <- res_df %>%
    filter(
      cluster_id == cluster_use,
      comparison == comparison_use
    ) %>%
    mutate(
      neg_log10_FDR = -log10(FDR),
      sig = case_when(
        FDR < 0.05 & logFC >= 0.5  ~ "Upregulated",
        FDR < 0.05 & logFC <= -0.5 ~ "Downregulated",
        TRUE ~ "Not significant"
      )
    )
  
  top_labels <- df %>%
    filter(sig != "Not significant") %>%
    arrange(FDR, desc(abs(logFC))) %>%
    slice_head(n = 15)
  
  p <- ggplot(df, aes(x = logFC, y = neg_log10_FDR, color = sig)) +
    geom_point(alpha = 0.75, size = 1.4) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    ggrepel::geom_text_repel(
      data = top_labels,
      aes(label = gene),
      size = 3,
      max.overlaps = Inf
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
      title = paste0("Cluster ", cluster_use, ": ", comparison_use),
      subtitle = "Pseudobulk DE within existing CD4 cluster",
      x = "log2 fold-change",
      y = "-log10(FDR)",
      color = "Result"
    )
  
  return(p)
}

if (nrow(all_cluster_DE) > 0) {
  
  volcano_dir <- file.path(pb_outdir, "Volcano_plots")
  dir.create(volcano_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (cl in unique(all_cluster_DE$cluster_id)) {
    for (comp in unique(all_cluster_DE$comparison)) {
      
      p_volcano <- plot_cluster_volcano(
        all_cluster_DE,
        cluster_use = cl,
        comparison_use = comp
      )
      
      ggsave(
        filename = file.path(
          volcano_dir,
          paste0("Volcano_cluster_", cl, "_", comp, ".png")
        ),
        plot = p_volcano,
        width = 8,
        height = 6,
        dpi = 300
      )
    }
  }
}


# Save final object and session info
CD4_T_dat <- obj

saveRDS(
  CD4_T_dat,
  file = file.path(outdir, "CD4_T_final_existing_cluster_based_analysis_object.rds")
)

sink(file.path(outdir, "sessionInfo.txt"))
sessionInfo()
sink()

cat("\nDone. Cluster-based analysis saved to:\n")
cat(outdir, "\n")



plotMDS(dge, labels = pb_meta$treatment, col = as.numeric(factor(pb_meta$treatment)))


logCPM <- edgeR::cpm(dge, log = TRUE, prior.count = 2)

pca <- prcomp(t(logCPM), scale. = TRUE)

pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  treatment = pb_meta$treatment,
  mouse = pb_meta$mouse_within_treatment
)

plot <- ggplot(pca_df, aes(PC1, PC2, color = treatment, label = mouse)) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel() +
  theme_classic()
ggsave(plot=plot, filename='cd4_pca.png')
#saveRDS(CD4_T_dat, "cleanCD4Tcell.rds")

