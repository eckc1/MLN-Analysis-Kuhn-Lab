library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)
library(dplyr)

# Ignore local RAM throttle 
mem.maxVSize(vsize = Inf)
# Load merged data
dat <- readRDS("noHashtag_Reclustered.rds")
head(dat)
Reductions(dat)
length(unique(dat$seurat_clusters))

# DEG Heatmaps
outdir <- "Top20_DEG_Heatmaps_By_Cluster_rescaled"

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

# Set assay and identities
DefaultAssay(dat) <- "RNA"
Idents(dat) <- "seurat_clusters"

# Find markers for all clusters
markers <- FindAllMarkers(
  object = dat,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# Save full marker table
#write.csv(markers,file = file.path(outdir, "AllClusterMarkers.csv"),row.names = FALSE)


# Get top 20 genes per cluster-
top20_markers <- markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 20) %>%
  ungroup()

#write.csv(top20_markers,file = file.path(outdir, "Top20Markers_PerCluster.csv"),row.names = FALSE)


# Scale all genes that will be used in heatmaps
all_heatmap_genes <- unique(top20_markers$gene)
all_heatmap_genes <- all_heatmap_genes[all_heatmap_genes %in% rownames(dat[["RNA"]])]

dat <- ScaleData(
  object = dat,
  assay = "RNA",
  features = all_heatmap_genes
)


# Save one heatmap per cluster
clusters <- sort(unique(top20_markers$cluster))

for (clust in clusters) {
  
  genes_use <- top20_markers %>%
    filter(cluster == clust) %>%
    pull(gene) %>%
    unique()
  
  genes_use <- genes_use[genes_use %in% rownames(dat[["RNA"]])]
  
  if (length(genes_use) == 0) {
    message("Skipping cluster ", clust, ": no valid genes found.")
    next
  }
  
  p <- DoHeatmap(
    object = dat,
    features = genes_use,
    group.by = "seurat_clusters",
    assay = "RNA",
    layer = "scale.data"
  ) +
    ggtitle(paste0("Top 20 DE Genes - Cluster ", clust)) +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
  
  #ggsave(
    #filename = file.path(outdir, paste0("Cluster_", clust, "_Top20_DEG_Heatmap.png")),
    #plot = p,width = 10,height = 8,dpi = 300
  #)
}

# save the updated object with the rescaled heatmap genes
#saveRDS(dat, file.path(outdir, "dat_with_top20_heatmap_genes_scaled.rds"))



### Annotate Full UMAP ###

# Make sure cluster IDs are character
dat$seurat_clusters <- as.character(dat$seurat_clusters)

# Define cluster groups
#dc_clusters <- c(16, 17)

#b_cell_clusters <- c(1, 13, 15)

#nk_clusters <- c(11)

#cd4_t_clusters <- c(6, 5, 4, 7, 2, 3)

#cd8_t_clusters <- c(0, 14, 18)


dc_clusters <- c(16)

b_cell_clusters <- c(1, 13)

nk_clusters <- c(17)

cd4_t_clusters <- c(12,9,8,7,6,5,4,3,2)

cd8_t_clusters <- c(0, 14, 18)


# Convert cluster vectors to character
dc_clusters <- as.character(dc_clusters)
b_cell_clusters <- as.character(b_cell_clusters)
nk_clusters <- as.character(nk_clusters)
cd4_t_clusters <- as.character(cd4_t_clusters)
cd8_t_clusters <- as.character(cd8_t_clusters)

# Find unassigned clusters
assigned_clusters <- unique(c(
  dc_clusters,
  b_cell_clusters,
  nk_clusters,
  cd4_t_clusters,
  cd8_t_clusters
))

all_clusters <- sort(unique(dat$seurat_clusters))
unassigned_clusters <- setdiff(all_clusters, assigned_clusters)

cat("Unassigned clusters:\n")
print(unassigned_clusters)

# Create annotation column
dat$celltype_broad <- "Unassigned"

dat$celltype_broad[dat$seurat_clusters %in% cd4_t_clusters] <- "CD4_T"
dat$celltype_broad[dat$seurat_clusters %in% cd8_t_clusters] <- "CD8_T"
dat$celltype_broad[dat$seurat_clusters %in% nk_clusters] <- "NK_cell"
dat$celltype_broad[dat$seurat_clusters %in% b_cell_clusters] <- "B_cell"
dat$celltype_broad[dat$seurat_clusters %in% dc_clusters] <- "DC"

# Check labels
table(dat$celltype_broad)
table(dat$seurat_clusters, dat$celltype_broad)

# Save table of number of cells per cluster and annotation
cluster_count_table <- dat@meta.data %>%
  count(seurat_clusters, celltype_broad, name = "n_cells") %>%
  arrange(celltype_broad, as.numeric(seurat_clusters))

write.csv(
  cluster_count_table,
  file = "Initial_cell_counts_per_cluster.csv",
  row.names = FALSE
)

# Also save broad cell type counts
celltype_count_table <- dat@meta.data %>%
  count(celltype_broad, name = "n_cells") %>%
  arrange(desc(n_cells))

write.csv(
  celltype_count_table,
  file = "Initial_cell_counts_per_celltype.csv",
  row.names = FALSE
)

# Subset to remove unassigned cells from plotting
dat_assigned <- subset(dat, subset = celltype_broad != "Unassigned")

# Plot only assigned cells
full_umap <- DimPlot(
  dat_assigned,
  reduction = "wnn.umap",
  group.by = "celltype_broad",
  label = FALSE,
  repel = TRUE
) +
  ggtitle("Broad cell type annotations")

full_umap

ggsave(
  filename = "Initial_full_umap_annotated_no_unassigned.png",
  plot = full_umap,width = 9,height = 7,dpi = 300
)

#saveRDS(dat, "dat_clusters_called.rds")



######## Clean Cell Clusters ##############

dat <- readRDS("dat_clusters_called.rds")

# Marker-based cleanup of broad cell type annotations

DefaultAssay(dat) <- "RNA"

#  keep only markers present in object
present_features <- function(obj, genes) {
  genes[genes %in% rownames(obj)]
}

# count how many marker genes are expressed per cell
marker_count <- function(obj, genes) {
  genes_use <- present_features(obj, genes)
  
  if (length(genes_use) == 0) {
    warning("None of these markers were found: ", paste(genes, collapse = ", "))
    return(rep(0, ncol(obj)))
  }
  
  mat <- GetAssayData(obj, assay = "RNA", layer = "data")[genes_use, , drop = FALSE]
  Matrix::colSums(mat > 0)
}


# Canonical mouse markers

t_markers <- c("Cd3d", "Cd3e", "Cd3g", "Trac")

cd4_markers <- c("Cd4", "Il7r", "Tcf7", "Lef1", "Ccr7")
cd8_markers <- c("Cd8a", "Cd8b1", "Gzmk", "Nkg7", "Ccl5")

b_markers <- c("Cd79a", "Cd79b", "Ms4a1", "Cd19", "Cd74", "H2-Ab1")
nk_markers <- c("Nkg7", "Ncr1", "Klrb1c", "Klrk1", "Prf1", "Gzmb", "Gzma")

dc_markers <- c("Itgax", "H2-Ab1", "H2-Eb1", "Cd74", "Flt3", "Clec10a", "Xcr1", "Zbtb46")


# Compute marker scores

dat$t_marker_n   <- marker_count(dat, t_markers)
dat$cd4_marker_n <- marker_count(dat, cd4_markers)
dat$cd8_marker_n <- marker_count(dat, cd8_markers)
dat$b_marker_n   <- marker_count(dat, b_markers)
dat$nk_marker_n  <- marker_count(dat, nk_markers)
dat$dc_marker_n  <- marker_count(dat, dc_markers)

# Moderate filtering based on gene expression
dat$marker_pass <- FALSE

dat$marker_pass[
  dat$celltype_broad == "CD4_T" &
    dat$t_marker_n >= 1 &
    dat$cd4_marker_n >= 1
] <- TRUE

dat$marker_pass[
  dat$celltype_broad == "CD8_T" &
    dat$t_marker_n >= 1 &
    dat$cd8_marker_n >= 1
] <- TRUE

dat$marker_pass[
  dat$celltype_broad == "B_cell" &
    dat$b_marker_n >= 1
    #dat$t_marker_n == 0
] <- TRUE

dat$marker_pass[
  dat$celltype_broad == "NK_cell" &
    dat$nk_marker_n >= 1 
    #dat$t_marker_n == 0
] <- TRUE

dat$marker_pass[
  dat$celltype_broad == "DC" &
    dat$dc_marker_n >= 1 
    #dat$b_marker_n < 1 &
    #dat$t_marker_n == 0
] <- TRUE


# Save before/after cleanup tables

marker_cleanup_table <- dat@meta.data %>%
  count(celltype_broad, marker_pass, name = "n_cells") %>%
  group_by(celltype_broad) %>%
  mutate(
    total_before = sum(n_cells),
    percent = round(100 * n_cells / total_before, 2)
  ) %>%
  ungroup() %>%
  arrange(celltype_broad, desc(marker_pass))

#write.csv(
  #marker_cleanup_table,
  #file = "Moderate_marker_cleanup_counts_by_celltype.csv",
  #row.names = FALSE
#)

cluster_marker_cleanup_table <- dat@meta.data %>%
  count(seurat_clusters, celltype_broad, marker_pass, name = "n_cells") %>%
  group_by(seurat_clusters, celltype_broad) %>%
  mutate(
    total_before = sum(n_cells),
    percent = round(100 * n_cells / total_before, 2)
  ) %>%
  ungroup() %>%
  arrange(celltype_broad, as.numeric(seurat_clusters), desc(marker_pass))

#write.csv(
 # cluster_marker_cleanup_table,
 # file = "Moderate_marker_cleanup_counts_by_cluster.csv",
 # row.names = FALSE
#)

print(marker_cleanup_table)
print(cluster_marker_cleanup_table)
tail(marker_cleanup_table)

# Apply marker filtering

dat_clean_markers <- subset(
  dat,
  subset = celltype_broad != "Unassigned" & marker_pass == TRUE
)

#saveRDS(dat_clean_markers, "dat_clean_marker_filtered_celltypes.rds")


# Plot cleaned annotation UMAP
clean_umap <- DimPlot(
  dat_clean_markers,
  reduction = "wnn.umap",
  group.by = "celltype_broad",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Marker-filtered broad cell type annotations")

clean_umap

#ggsave(
  #filename = "Moderate_full_umap_marker_filtered_celltypes.png",
  #plot = clean_umap,
  #width = 9,height = 7,dpi = 300)

##################################################

# More strict Filtering
# T lineage
t_core_markers <- c("Cd3d", "Cd3e", "Cd3g", "Trac")

# CD4 identity
cd4_core_markers <- c("Cd4")

# CD8 identity
cd8_core_markers <- c("Cd8a", "Cd8b1")

# B identity
b_core_markers <- c("Cd79a", "Cd79b", "Ms4a1", "Cd19")

# NK identity
nk_core_markers <- c("Ncr1", "Klrb1c", "Klrk1")

# Cytotoxic markers: useful, but not NK-specific
cytotoxic_markers <- c("Nkg7", "Prf1", "Gzmb", "Gzma", "Ccl5")

# DC / APC
dc_core_markers <- c("Itgax", "Flt3", "Zbtb46", "Clec10a", "Xcr1")

# MHC-II/APC markers: shared by B and DC
apc_shared_markers <- c("Cd74", "H2-Ab1", "H2-Eb1")

dat$t_core_n      <- marker_count(dat, t_core_markers)
dat$cd4_core_n    <- marker_count(dat, cd4_core_markers)
dat$cd8_core_n    <- marker_count(dat, cd8_core_markers)
dat$b_core_n      <- marker_count(dat, b_core_markers)
dat$nk_core_n     <- marker_count(dat, nk_core_markers)
dat$cytotoxic_n   <- marker_count(dat, cytotoxic_markers)
dat$dc_core_n     <- marker_count(dat, dc_core_markers)
dat$apc_shared_n  <- marker_count(dat, apc_shared_markers)

dat$marker_pass_strict <- FALSE

# CD4 T cells: T lineage plus CD4, not CD8
dat$marker_pass_strict[
  dat$celltype_broad == "CD4_T" &
    dat$t_core_n >= 1 &
    dat$cd4_core_n >= 1 &
    dat$cd8_core_n == 0
] <- TRUE

# CD8 T cells: T lineage plus CD8, not CD4
dat$marker_pass_strict[
  dat$celltype_broad == "CD8_T" &
    dat$t_core_n >= 1 &
    dat$cd8_core_n >= 1 &
    dat$cd4_core_n == 0
] <- TRUE

# B cells: B identity markers, no T, no NK
dat$marker_pass_strict[
  dat$celltype_broad == "B_cell" &
    dat$b_core_n >= 2 &
    dat$t_core_n == 0 &
    dat$nk_core_n == 0
] <- TRUE

# NK cells: NK/cytotoxic, but no T lineage
dat$marker_pass_strict[
  dat$celltype_broad == "NK_cell" &
    dat$nk_core_n >= 1 &
    dat$cytotoxic_n >= 1 &
    dat$t_core_n == 0
] <- TRUE

# DC cells: DC/APC markers, not B, not T, not NK
dat$marker_pass_strict[
  dat$celltype_broad == "DC" &
    (
      dat$dc_core_n >= 1 |
        dat$apc_shared_n >= 1
    ) &
    dat$b_core_n == 0 
    #dat$t_core_n == 0 &
    #dat$nk_core_n == 0
] <- TRUE

dat_clean_strict <- subset(
  dat,
  subset = celltype_broad != "Unassigned" & marker_pass_strict == TRUE
)

table(dat_clean_strict$celltype_broad)
table(dat_clean_strict$seurat_clusters, dat_clean_strict$celltype_broad)

strict_cleanup_table <- dat@meta.data %>%
  count(celltype_broad, marker_pass_strict, name = "n_cells") %>%
  group_by(celltype_broad) %>%
  mutate(
    total_before = sum(n_cells),
    percent = round(100 * n_cells / total_before, 2)
  ) %>%
  ungroup() %>%
  arrange(celltype_broad, desc(marker_pass_strict))

write.csv(
  strict_cleanup_table,
  "strict_marker_cleanup_counts_by_celltype.csv",
  row.names = FALSE
)

strict_cluster_cleanup_table <- dat@meta.data %>%
  count(seurat_clusters, celltype_broad, marker_pass_strict, name = "n_cells") %>%
  group_by(seurat_clusters, celltype_broad) %>%
  mutate(
    total_before = sum(n_cells),
    percent = round(100 * n_cells / total_before, 2)
  ) %>%
  ungroup() %>%
  arrange(celltype_broad, as.numeric(seurat_clusters), desc(marker_pass_strict))

write.csv(
  strict_cluster_cleanup_table,
  "strict_marker_cleanup_counts_by_cluster.csv",
  row.names = FALSE
)
dc_table <- strict_cleanup_table[strict_cleanup_table$celltype_broad=="DC",]
dc_table
strict_umap <- DimPlot(
  dat_clean_strict,
  reduction = "wnn.umap",
  group.by = "celltype_broad",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Strict marker-filtered broad cell type annotations")

strict_umap

ggsave(
  filename = "full_umap_strict_marker_filtered_celltypes.png",
  plot = strict_umap,
  width = 9,
  height = 7,
  dpi = 300
)
strict_umap_pool <- DimPlot(
  dat_clean_strict,
  reduction = "wnn.umap",
  group.by = "pool_id",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Strict marker-filtered broad cell type annotations")

ggsave(
 filename = "full_umap_strict_pool.png",
  plot = strict_umap_pool,
  width = 9,
  height = 7,
  dpi = 300
)

#saveRDS(dat_clean_strict, "dat_clean_strict.rds")

#saveRDS(dat_clean_strict, "dat_clean_strict_adjusted_dc_cuttoff.rds")

############## Feature Plots



# Clean
# Feature plots for canonical markers
dat <- readRDS("dat_clean_marker_filtered_celltypes.rds")
DefaultAssay(dat) <- "RNA"

# Marker lists
marker_lists <- list(
  T_markers = c("Cd3d", "Cd3e", "Cd3g", "Trac"),
  CD4_markers = c("Cd4", "Il7r", "Tcf7", "Lef1", "Ccr7"),
  CD8_markers = c("Cd8a", "Cd8b1", "Gzmk", "Nkg7", "Ccl5"),
  B_markers = c("Cd79a", "Cd79b", "Ms4a1", "Cd19", "Cd74", "H2-Ab1"),
  NK_markers = c("Nkg7", "Ncr1", "Klrb1c", "Klrk1", "Prf1", "Gzmb", "Gzma"),
  DC_markers = c("Itgax", "H2-Ab1", "H2-Eb1", "Cd74", "Flt3", "Clec10a", "Xcr1", "Zbtb46")
)

# Choose reduction
# Use "wnn.umap" if this is your original WNN object.
# Use "umap" if this is your reclustered/recomputed object.
reduction_use <- if ("wnn.umap" %in% Reductions(dat)) {
  "wnn.umap"
} else if ("umap" %in% Reductions(dat)) {
  "umap"
} else {
  stop("No UMAP reduction found. Expected 'wnn.umap' or 'umap'.")
}

cat("Using reduction:", reduction_use, "\n")

# Output folder
outdir <- "FeaturePlots_Canonical_Markers"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Save individual and combined feature plots

for (marker_group in names(marker_lists)) {
  
  message("Processing: ", marker_group)
  
  group_dir <- file.path(outdir, marker_group)
  dir.create(group_dir, showWarnings = FALSE, recursive = TRUE)
  
  markers <- marker_lists[[marker_group]]
  markers_present <- markers[markers %in% rownames(dat)]
  markers_missing <- setdiff(markers, markers_present)
  
  if (length(markers_missing) > 0) {
    message("Missing markers in ", marker_group, ": ", paste(markers_missing, collapse = ", "))
  }
  
  if (length(markers_present) == 0) {
    warning("No markers found for ", marker_group)
    next
  }
  
  # Save individual plots
  for (gene in markers_present) {
    
    p <- FeaturePlot(
      dat,
      features = gene,
      reduction = reduction_use,
      order = TRUE,
      pt.size = 0.2
    ) +
      ggtitle(paste0(marker_group, ": ", gene))
    
    ggsave(
      filename = file.path(group_dir, paste0(gene, "_FeaturePlot.png")),
      plot = p,
      width = 6,
      height = 5,
      dpi = 300
    )
  }
  
  # Save combined grid
  p_list <- FeaturePlot(
    dat,
    features = markers_present,
    reduction = reduction_use,
    order = TRUE,
    pt.size = 0.2,
    combine = FALSE
  )
  
  p_combined <- wrap_plots(p_list, ncol = 3) +
    plot_annotation(title = paste0(marker_group, " FeaturePlots"))
  
  ggsave(
    filename = file.path(outdir, paste0(marker_group, "_combined_FeaturePlots.png")),
    plot = p_combined,
    width = 15,
    height = ceiling(length(markers_present) / 3) * 4,
    dpi = 300,
    limitsize = FALSE
  )
}

# Save marker presence/missing summary
marker_summary <- bind_rows(
  lapply(names(marker_lists), function(marker_group) {
    markers <- marker_lists[[marker_group]]
    
    data.frame(
      marker_group = marker_group,
      marker = markers,
      present_in_object = markers %in% rownames(dat)
    )
  })
)

write.csv(
  marker_summary,
  file = file.path(outdir, "canonical_marker_presence_summary.csv"),
  row.names = FALSE
)

marker_summary


# Strict

# Feature plots for canonical markers
dat <- readRDS("dat_clean_strict.rds")
DefaultAssay(dat) <- "RNA"

# Marker lists
marker_lists <- list(
  T_markers = c("Cd3d", "Cd3e", "Cd3g", "Trac"),
  CD4_markers = c("Cd4", "Il7r", "Tcf7", "Lef1", "Ccr7"),
  CD8_markers = c("Cd8a", "Cd8b1", "Gzmk", "Nkg7", "Ccl5"),
  B_markers = c("Cd79a", "Cd79b", "Ms4a1", "Cd19", "Cd74", "H2-Ab1"),
  NK_markers = c("Nkg7", "Ncr1", "Klrb1c", "Klrk1", "Prf1", "Gzmb", "Gzma"),
  DC_markers = c("Itgax", "H2-Ab1", "H2-Eb1", "Cd74", "Flt3", "Clec10a", "Xcr1", "Zbtb46")
)


reduction_use <- if ("wnn.umap" %in% Reductions(dat)) {
  "wnn.umap"
} else if ("umap" %in% Reductions(dat)) {
  "umap"
} else {
  stop("No UMAP reduction found. Expected 'wnn.umap' or 'umap'.")
}

cat("Using reduction:", reduction_use, "\n")

# Output folder
outdir <- "Strict_FeaturePlots_Canonical_Markers"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)


# Save individual and combined feature plots

for (marker_group in names(marker_lists)) {
  
  message("Processing: ", marker_group)
  
  group_dir <- file.path(outdir, marker_group)
  dir.create(group_dir, showWarnings = FALSE, recursive = TRUE)
  
  markers <- marker_lists[[marker_group]]
  markers_present <- markers[markers %in% rownames(dat)]
  markers_missing <- setdiff(markers, markers_present)
  
  if (length(markers_missing) > 0) {
    message("Missing markers in ", marker_group, ": ", paste(markers_missing, collapse = ", "))
  }
  
  if (length(markers_present) == 0) {
    warning("No markers found for ", marker_group)
    next
  }
  
  # Save individual plots
  for (gene in markers_present) {
    
    p <- FeaturePlot(
      dat,
      features = gene,
      reduction = reduction_use,
      order = TRUE,
      pt.size = 0.2
    ) +
      ggtitle(paste0(marker_group, ": ", gene))
    
    ggsave(
      filename = file.path(group_dir, paste0(gene, "_FeaturePlot.png")),
      plot = p,
      width = 6,
      height = 5,
      dpi = 300
    )
  }
  
  # Save combined grid
  p_list <- FeaturePlot(
    dat,
    features = markers_present,
    reduction = reduction_use,
    order = TRUE,
    pt.size = 0.2,
    combine = FALSE
  )
  
  p_combined <- wrap_plots(p_list, ncol = 3) +
    plot_annotation(title = paste0(marker_group, " FeaturePlots"))
  
  ggsave(
    filename = file.path(outdir, paste0(marker_group, "_combined_FeaturePlots.png")),
    plot = p_combined,
    width = 15,
    height = ceiling(length(markers_present) / 3) * 4,
    dpi = 300,
    limitsize = FALSE
  )
  
}


# Save marker presence/missing summary

marker_summary <- bind_rows(
  lapply(names(marker_lists), function(marker_group) {
    markers <- marker_lists[[marker_group]]
    
    data.frame(
      marker_group = marker_group,
      marker = markers,
      present_in_object = markers %in% rownames(dat)
    )
  })
)

write.csv(
  marker_summary,
  file = file.path(outdir, "canonical_marker_presence_summary.csv"),
  row.names = FALSE
)

marker_summary
