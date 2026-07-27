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
# Subset to only Bcell cells from your broad annotation
Bcell_dat <- subset(dat, subset = celltype_broad == "B_cell")
colnames(Bcell_dat@meta.data)
# Plot Bcell cells only
Bcell_umap <- DimPlot(
  Bcell_dat,
  reduction = "wnn.umap",
  group.by = "pool_id",
  #group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Bcell T Cell clusters only")
ggsave(plot = Bcell_umap, filename="Bcell_umap.png")
#ggsave(plot = Bcell_umap, filename="Bcell_umap_transcript.png")

class(Bcell_dat)
dim(Bcell_dat)

table(Bcell_dat@meta.data$seurat_clusters)

Reductions(Bcell_dat)

colnames(Bcell_dat@meta.data)
table(Bcell_dat$pool_id)
table(Bcell_dat$orig.ident, Bcell_dat$pool_id)
table(Bcell_dat$demux_id, Bcell_dat$pool_id)


# Reculster Bcell T cell object

DefaultAssay(Bcell_dat) <- "RNA"
Bcell_dat <- NormalizeData(Bcell_dat)
Bcell_dat <- FindVariableFeatures(Bcell_dat)
Bcell_dat <- ScaleData(Bcell_dat)
Bcell_dat <- RunPCA(Bcell_dat, reduction.name = "pca")

DefaultAssay(Bcell_dat) <- "ADT"
Bcell_dat <- NormalizeData(Bcell_dat, normalization.method = "CLR", margin = 2)
Bcell_dat <- ScaleData(Bcell_dat)
Bcell_dat <- RunPCA(Bcell_dat, reduction.name = "apca")

Bcell_dat <- FindMultiModalNeighbors(
  Bcell_dat,
  reduction.list = list("pca", "apca"),
  dims.list = list(1:20, 1:18)
)

Bcell_dat <- FindClusters(
  Bcell_dat,
  graph.name = "wsnn",
  resolution = 0.3
)

Bcell_dat <- RunUMAP(
  Bcell_dat,
  nn.name = "weighted.nn",
  reduction.name = "Bcell.wnn.umap",
  reduction.key = "BcellwnnUMAP_"
)

Bcell_umap <- DimPlot(
  Bcell_dat,
  reduction = "wnn.umap",
  group.by = "pool_id",
  label = TRUE,repel = TRUE
) +
  ggtitle("B Cell clusters only")
ggsave(plot = Bcell_umap, filename="reclustered_Bcell_umap.png")


Bcell_umap <- DimPlot(
  Bcell_dat,
  reduction = "wnn.umap",
  group.by = "seurat_clusters",
  label = TRUE,repel = TRUE
) +
  ggtitle("B Cell clusters only")
ggsave(plot = Bcell_umap, filename="reclustered_Bcell_umap_transcript.png")


outdir <- "/Users/eckco/Desktop/Kuhn_Lab/Jing_Data/MLN/Demultiplexed/Recluster_No_Hashtag/Bcell_Analysis/Bcell_Pseudobulk_Volcano_By_demux_id_unpaired"

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
print(table(Bcell_dat@meta.data[[group_col]], Bcell_dat@meta.data[[replicate_col]]))


# Create pseudobulk IDs, which represent unique mouse-level pseudobulk samples.
# pool_id repeats across treatments, so the true unique mouse sample is: treatment + pool_id

Bcell_dat$treatment <- as.character(Bcell_dat@meta.data[[group_col]])
Bcell_dat$mouse_within_treatment <- as.character(Bcell_dat@meta.data[[replicate_col]])

Bcell_dat$pseudobulk_id <- paste0(
  "treat.", Bcell_dat$treatment,
  ".mouse.", Bcell_dat$mouse_within_treatment
)

# Build pseudobulk metadata
pb_meta <- Bcell_dat@meta.data %>%
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
  object = Bcell_dat,
  assay = "RNA",
  layer = "counts"
)

cell_meta <- Bcell_dat@meta.data[colnames(counts), , drop = FALSE]

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
  file = file.path(outdir, "Bcell_pseudobulk_metadata.csv"),
  row.names = TRUE
)

write.csv(
  as.data.frame(as.matrix(pb_counts)) %>% rownames_to_column("gene"),
  file = file.path(outdir, "Bcell_pseudobulk_counts.csv"),
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
  file = file.path(outdir, "Bcell_pseudobulk_pairwise_demux_id_DE_results.csv"),
  row.names = FALSE
)


# Raw p-value DEG results and volcano plots. Exploratory results: PValue < 0.05

pvalue_cutoff <- 0.05

# Save all raw-p significant genes
rawp_sig_results <- all_results %>%
  mutate(
    sig_rawP = case_when(
      PValue < pvalue_cutoff & logFC >= logfc_cutoff  ~ "Upregulated",
      PValue < pvalue_cutoff & logFC <= -logfc_cutoff ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  ) %>%
  filter(sig_rawP != "Not significant")

write.csv(
  rawp_sig_results,
  file = file.path(outdir, "Bcell_significant_DE_genes_rawP0.05_logFC0.5.csv"),
  row.names = FALSE
)

# Save top 10 up/down genes by raw p-value per comparison
rawp_top_labeled_genes <- all_results %>%
  mutate(
    sig_rawP = case_when(
      PValue < pvalue_cutoff & logFC >= logfc_cutoff  ~ "Upregulated",
      PValue < pvalue_cutoff & logFC <= -logfc_cutoff ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  ) %>%
  filter(sig_rawP != "Not significant") %>%
  group_by(comparison, sig_rawP) %>%
  arrange(PValue, desc(abs(logFC)), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

write.csv(
  rawp_top_labeled_genes,
  file = file.path(outdir, "Bcell_top10_up_top10_down_labeled_genes_rawP0.05_logFC0.5.csv"),
  row.names = FALSE
)

# Raw-p volcano plotting function
plot_volcano_rawP <- function(res_df, comparison_name) {
  
  df <- res_df %>%
    filter(comparison == comparison_name) %>%
    mutate(
      neg_log10_PValue = -log10(PValue),
      sig_rawP = case_when(
        PValue < pvalue_cutoff & logFC >= logfc_cutoff  ~ "Upregulated",
        PValue < pvalue_cutoff & logFC <= -logfc_cutoff ~ "Downregulated",
        TRUE ~ "Not significant"
      )
    )
  
  top_up <- df %>%
    filter(sig_rawP == "Upregulated") %>%
    arrange(PValue, desc(logFC)) %>%
    slice_head(n = 10)
  
  top_down <- df %>%
    filter(sig_rawP == "Downregulated") %>%
    arrange(PValue, logFC) %>%
    slice_head(n = 10)
  
  top_labels <- bind_rows(top_up, top_down)
  
  p <- ggplot(df, aes(x = logFC, y = neg_log10_PValue, color = sig_rawP)) +
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
    theme_classic() +
    labs(
      title = paste0("Pseudobulk volcano, raw p-value: ", comparison_name),
      subtitle = paste0(
        "Exploratory threshold: raw PValue < ", pvalue_cutoff,
        " and |log2FC| >= ", logfc_cutoff
      ),
      x = "log2 fold-change",
      y = "-log10(raw PValue)",
      color = "Result"
    )
  
  return(p)
}

# Save raw-p volcano plots
for (comp_name in comparison_names) {
  
  p_raw <- plot_volcano_rawP(all_results, comp_name)
  
  ggsave(
    filename = file.path(outdir, paste0("Volcano_rawP0.05_", comp_name, ".png")),
    plot = p_raw,
    width = 8,
    height = 6,
    dpi = 300
  )
}


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
  file = file.path(outdir, "Bcell_significant_DE_genes_FDR0.05_logFC1.csv"),
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
  file = file.path(outdir, "Bcell_top10_up_top10_down_labeled_genes.csv"),
  row.names = FALSE
)


# session info

sink(file.path(outdir, "sessionInfo.txt"))
sessionInfo()
sink()

cat("\nDone. Results saved to:\n")
cat(outdir, "\n")


dim(all_results)
head(all_results)

table(all_results$comparison)
summary(all_results$FDR)
summary(all_results$logFC)

dim(sig_results)




all_results %>%
  group_by(comparison) %>%
  summarise(
    n_genes = n(),
    min_PValue = min(PValue, na.rm = TRUE),
    min_FDR = min(FDR, na.rm = TRUE),
    n_PValue_0.05 = sum(PValue < 0.05, na.rm = TRUE),
    n_PValue_0.01 = sum(PValue < 0.01, na.rm = TRUE),
    n_FDR_0.1 = sum(FDR < 0.1, na.rm = TRUE),
    n_FDR_0.05 = sum(FDR < 0.05, na.rm = TRUE),
    .groups = "drop"
  )


logCPM <- edgeR::cpm(dge, log = TRUE, prior.count = 2)

pca <- prcomp(t(logCPM), scale. = TRUE)

pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  treatment = pb_meta$treatment,
  mouse = pb_meta$mouse_within_treatment
)


exploratory_results <- all_results %>%
  filter(PValue < 0.05 & abs(logFC) >= 0.5)

write.csv(
  exploratory_results,
  file = file.path(outdir, "Bcell_exploratory_DE_genes_PValue0.05_logFC0.5.csv"),
  row.names = FALSE
)

logCPM <- edgeR::cpm(dge, log = TRUE, prior.count = 2)

pca <- prcomp(t(logCPM), scale. = TRUE)

pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  treatment = pb_meta$treatment,
  mouse = pb_meta$mouse_within_treatment
)

ggplot(pca_df, aes(PC1, PC2, color = treatment, label = mouse)) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel() +
  theme_classic() +
  labs(title = "B-cell pseudobulk PCA")

