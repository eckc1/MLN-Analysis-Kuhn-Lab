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
library(readr)

# Ignore local RAM throttle 
mem.maxVSize(vsize = Inf)
# Load merged data
#dat <- readRDS("dat_clusters_called.rds")
dat <- readRDS("cleanCD4Tcell.rds")
Reductions(dat)
DimPlot(dat, reduction = "wnn.umap", label=TRUE)

# Top DEGs per cluster

outdir <- "Reclustered_CD4_Top20_DEGs_Per_Cluster"

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

# Make sure RNA is active assay for gene expression markers
DefaultAssay(dat) <- "RNA"

# Use current Seurat clusters as identities
Idents(dat) <- "seurat_clusters"

# Check cluster sizes
cluster_sizes <- table(Idents(dat))
print(cluster_sizes)

# Find positive marker genes for each cluster vs all other clusters
cluster_markers <- FindAllMarkers(
  object = dat,
  assay = "RNA",
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.10,
  logfc.threshold = 0.25
)

# Save all marker results
#write.csv(
  #cluster_markers,
  #file = file.path(outdir, "All_DEGs_Per_Cluster.csv"),
  #row.names = FALSE
#)

# Get top 20 genes per cluster
top20_markers <- cluster_markers %>%
  group_by(cluster) %>%
  arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 20) %>%
  ungroup()

# Save top 20 per cluster
#write.csv(
  #top20_markers,
  #file = file.path(outdir, "Top20_DEGs_Per_Cluster.csv"),
  #row.names = FALSE
#)

# save just gene names in wide format: one column per cluster
top20_gene_list <- top20_markers %>%
  group_by(cluster) %>%
  mutate(rank = row_number()) %>%
  ungroup() %>%
  dplyr::select(cluster, rank, gene) %>%
  tidyr::pivot_wider(
    names_from = cluster,
    values_from = gene,
    names_prefix = "cluster_"
  )

#write.csv(
 # top20_gene_list,
  #file = file.path(outdir, "Top20_Gene_Names_Per_Cluster_Wide.csv"),
  #row.names = FALSE
#)

cat("\nDone. Files saved to:\n")
cat(outdir, "\n")

# Annotated UMAP
cluster_labels <- c(
  "0"  = "Resting/Cxcr4+ CD4",
  "1"  = "Ribosomal-high CD4",
  "2"  = "IFN/Ly6c CD4",
  "3"  = "Cytotoxic-like CD4",
  "4"  = "Naive/CM-like CD4",
  "5"  = "Activated/metabolic CD4",
  "6"  = "Treg",
  "7"  = "Ribosomal-high CD4",
  "8"  = "Activated inhibitory CD4",
  "9"  = "IFN/stress-like CD4",
  "10" = "Activated CD4",
  "11" = "Stress/low-confidence CD4",
  "12" = "Th17-like CD4"
)
# Create labels in the same order as cells
tcell_subset_labels <- unname(cluster_labels[as.character(dat$seurat_clusters)])


# Check for missing labels
table(tcell_subset_labels, useNA = "ifany")

# Add labels to metadata
dat <- AddMetaData(
  object = dat,
  metadata = tcell_subset_labels,
  col.name = "tcell_subset"
)

# Remove cells with NA tcell_subset labels
dat <- subset(
  dat,
  subset = !is.na(tcell_subset)
)

table(dat@meta.data$seurat_clusters, dat@meta.data$tcell_subset, useNA = "ifany")

annotate_cd4 <-DimPlot(
  dat,
  reduction = "wnn.umap",
  group.by = "tcell_subset",
  label = FALSE,
  repel = TRUE
)
#ggsave(plot=annotate_cd4, filename="/Users/eckco/Desktop/Kuhn_Lab/Jing_Data/MLN/Demultiplexed/Recluster_No_Hashtag/CD4_TCell_Analysis/Reclustered_CD4_Top20_DEGs_Per_Cluster/annotate_cd4.png")

annotate_cd4<- DimPlot(
  dat,
  reduction = "wnn.umap",
  label = TRUE,
  repel = TRUE
)
#ggsave(plot=annotate_cd4, filename="/Users/eckco/Desktop/Kuhn_Lab/Jing_Data/MLN/Demultiplexed/Recluster_No_Hashtag/CD4_TCell_Analysis/Reclustered_CD4_Top20_DEGs_Per_Cluster/transcript_annotate_cd4.png")



# Feature plot for CCR6 and Klrb1a
p <- FeaturePlot(
  dat,
  features = c("Ccr6", "Klrb1a"),
  reduction = "wnn.umap",
  min.cutoff = "q10",
  max.cutoff = "q99",
  cols = c("lightgrey", "blue"),
  order = TRUE,
  combine = TRUE
)
#ggsave(plot = p, filename = "/Users/eckco/Desktop/Kuhn_Lab/Jing_Data/MLN/Demultiplexed/Recluster_No_Hashtag/CD4_TCell_Analysis/Reclustered_CD4_Top20_DEGs_Per_Cluster/TH17_markers.png")


# Define cell subsets
print(rownames(dat[["ADT"]]))

DefaultAssay(dat) <- "ADT"
adt_counts <- GetAssayData(dat, assay = "ADT", layer = "counts")

cells_keep <- colnames(dat)[
  adt_counts["Ms.TCR.Bchain", ] > 1 &
    adt_counts["Ms.CD3", ] > 1 &
    adt_counts["Ms.CD4", ] > 1
]
length(cells_keep)

DefaultAssay(dat) <- "RNA"
RNA_counts <- GetAssayData(dat, assay = "RNA", layer = "counts")

write.csv(
  data.frame(RNA_genes = rownames(dat)),
  file = "RNA_genes.csv",
  row.names = FALSE
)
# Define TH17 Cells
cells_keep_Additional_subset <- colnames(dat)[
  (
    RNA_counts["Ccr6", ] > 0.5 |
      RNA_counts["Klrb1a", ] > 0.5 |
      RNA_counts["Rorc", ] > 0.5
  ) &
    RNA_counts["Il21r", ] > 0.5
]
length(cells_keep_Additional_subset)


cluster1_sub <- subset(dat, cells = cells_keep_Additional_subset)

## TH17 Numbers

# Check available metadata columns
colnames(dat$pool_id)
head(dat)
table(dat$treatment)
# treatment column 
treatment_col <- "treatment"

# Get Th17 cell names
th17_cells <- colnames(cluster1_sub)

# Build dataframe of all cells with Th17 status and treatment group
th17_bar_df <- data.frame(
  cell = colnames(dat),
  treatment = dat@meta.data[colnames(dat), treatment_col],
  is_th17 = colnames(dat) %in% th17_cells
)

# Count number of Th17 cells per treatment group
th17_counts <- th17_bar_df %>%
  filter(is_th17) %>%
  group_by(treatment) %>%
  summarise(
    Th17_cells = n(),
    .groups = "drop"
  )

# Bar plot: number of Th17 cells per treatment group
p_th17_bar <- ggplot(th17_counts, aes(x = treatment, y = Th17_cells, fill = treatment)) +
  geom_col(width = 0.75, color = "black") +
  geom_text(
    aes(label = Th17_cells),
    vjust = -0.4,
    size = 5
  ) +
  theme_bw() +
  labs(
    title = "Number of Th17 Cells per Treatment Group",
    x = "Treatment Group",
    y = "Number of Th17 Cells",
    fill = "Treatment Group"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_th17_bar)

ggsave(
  "Th17_cells_per_treatment_group.png",
  p_th17_bar,
  width = 8,
  height = 6,
  dpi = 300
)


# Count total cells and Th17 cells per treatment group
th17_percent_df <- th17_bar_df %>%
  group_by(treatment) %>%
  summarise(
    Total_cells = n(),
    Th17_cells = sum(is_th17),
    Percent_Th17 = 100 * Th17_cells / Total_cells,
    .groups = "drop"
  )

# Bar plot: percent Th17 cells per treatment group
p_th17_percent <- ggplot(th17_percent_df, aes(x = treatment, y = Percent_Th17, fill = treatment)) +
  geom_col(width = 0.75, color = "black") +
  geom_text(
    aes(label = paste0(round(Percent_Th17, 2), "%")),
    vjust = -0.4,
    size = 5
  ) +
  theme_bw() +
  labs(
    title = "Percent Th17 Cells per Treatment Group",
    x = "Treatment Group",
    y = "Percent Th17 Cells",
    fill = "Treatment Group"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_th17_percent)

ggsave(
  "Percent_Th17_cells_per_treatment_group.png",
  p_th17_percent,
  width = 8,
  height = 6,
  dpi = 300
)

DimPlot(cluster1_sub)
dim(cluster1_sub)
head(cluster1_sub)

# Make pie chart of cluster1_sub as a slice of dat
cells_in_sub <- colnames(cluster1_sub)
all_cells <- colnames(dat)

pie_df <- data.frame(
  Group = c("cluster1_sub", "Other cells in dat"),
  Count = c(
    length(cells_in_sub),
    length(setdiff(all_cells, cells_in_sub))
  )
)

pie_df <- pie_df %>%
  mutate(
    Group_label = case_when(
      Group == "cluster1_sub" ~ "Th17 Cells",
      Group == "Other cells in dat" ~ "Other",
      TRUE ~ Group
    ),
    Fraction = Count / sum(Count),
    Percent = Fraction * 100,
    Label = paste0(
      Group_label, "\n",
      Count, " cells\n",
      round(Percent, 2), "%"
    )
  )

p_pie <- ggplot(pie_df, aes(x = "", y = Count, fill = Group_label)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  theme_void() +
  labs(
    title = "Percent Th17 Cells relative to full CD4 T Cell Cluster",
    fill = "Group"
  ) +
  geom_text(
    aes(label = Label),
    position = position_stack(vjust = 0.5),
    size = 5
  )

print(p_pie)
#ggsave("cluster1_sub_pie_chart.png", p_pie, width = 7, height = 7, dpi = 300)

outdir <- getwd()

# Ridge plots and feature plots
DefaultAssay(cluster1_sub) <- "RNA"
cluster1_sub <- NormalizeData(cluster1_sub, verbose = FALSE)

genes_to_plot <- c("Rora", "Rorc")
genes_present <- genes_to_plot[genes_to_plot %in% rownames(cluster1_sub)]

if (length(genes_present) == 0) {
  stop("None of the requested genes were found in the object.")
}

p_ridge <- RidgePlot(
  object = cluster1_sub,
  features = genes_present,
  group.by = "seurat_clusters",
  assay = "RNA",
  ncol = 1
) &
  theme_bw() &
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold")
  )

#ggsave(
 # filename = file.path(outdir, "RORA_RORC_RidgePlot.png"),
 # plot = p_ridge,
  #width = 10,
  #height = 8,
  #dpi = 300
#)

p_feature <- FeaturePlot(
  object = cluster1_sub,
  features = genes_present,
  reduction = "wnn.umap",
  cols = c("lightgrey", "red"),
  ncol = 2,
  order = TRUE
) &
  theme_bw() &
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold")
  )

#ggsave(
 # filename = file.path(outdir, "RORA_RORC_FeaturePlot.png"),
 # plot = p_feature,
 # width = 12,
 # height = 5,
  #dpi = 300
#)



# =================================
# Violin plots: Ccr6 and Klrb1a
# Comparisons: 602 vs 603, and 602 vs 601
# =================================

DefaultAssay(cluster1_sub) <- "RNA"

# Normalize if needed
cluster1_sub <- NormalizeData(cluster1_sub, verbose = FALSE)

# Make sure grouping column exists
if (!"pool_id" %in% colnames(cluster1_sub@meta.data)) {
  stop("pool_id column not found in cluster1_sub@meta.data")
}

# Genes of interest (case-insensitive matching to rownames)
genes_requested <- c("Ccr6", "Klrb1a")

match_gene <- function(gene, gene_names) {
  hit <- gene_names[tolower(gene_names) == tolower(gene)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

genes_present <- sapply(genes_requested, match_gene, gene_names = rownames(cluster1_sub))
names(genes_present) <- genes_requested

if (any(is.na(genes_present))) {
  missing_genes <- names(genes_present)[is.na(genes_present)]
  stop("These genes were not found in the object: ", paste(missing_genes, collapse = ", "))
}

print(genes_present)

# Function to make one comparison plot with both genes
make_violin_compare_plot <- function(
    seu,
    group_col = "pool_id",
    group1,
    group2,
    genes_vec,
    outdir,
    filename
) {
  
  meta_col <- seu@meta.data[[group_col]]
  keep_cells <- rownames(seu@meta.data)[meta_col %in% c(group1, group2)]
  
  if (length(keep_cells) == 0) {
    stop("No cells found for comparison: ", group1, " vs ", group2)
  }
  
  sub <- subset(seu, cells = keep_cells)
  sub[[group_col]] <- factor(sub[[group_col]][,1], levels = c(group1, group2))
  
  plot_df <- FetchData(sub, vars = c(group_col, genes_vec))
  colnames(plot_df)[1] <- "group"
  
  long_df <- do.call(
    rbind,
    lapply(genes_vec, function(g) {
      data.frame(
        group = plot_df$group,
        gene = g,
        expression = plot_df[[g]],
        stringsAsFactors = FALSE
      )
    })
  )
  
  long_df$group <- factor(long_df$group, levels = c(group1, group2))
  long_df$gene <- factor(long_df$gene, levels = genes_vec)
  
  # Compute p-values per gene
  pval_df <- do.call(
    rbind,
    lapply(genes_vec, function(g) {
      sub_df <- long_df %>% filter(gene == g)
      
      wt <- wilcox.test(expression ~ group, data = sub_df)
      ymax <- max(sub_df$expression, na.rm = TRUE)
      
      data.frame(
        gene = g,
        p_value = wt$p.value,
        label = paste0("p = ", formatC(wt$p.value, format = "e", digits = 2)),
        y_pos = ymax * 1.08 + 0.05,
        x = 1.5,
        stringsAsFactors = FALSE
      )
    })
  )
  
  p <- ggplot(long_df, aes(x = group, y = expression, fill = group)) +
    geom_violin(trim = FALSE, scale = "width", alpha = 0.7) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9) +
    geom_jitter(width = 0.15, size = 0.3, alpha = 0.35) +
    facet_wrap(~ gene, scales = "free_y", ncol = 2) +
    geom_segment(
      data = pval_df,
      aes(x = 1, xend = 2, y = y_pos, yend = y_pos),
      inherit.aes = FALSE
    ) +
    geom_text(
      data = pval_df,
      aes(x = x, y = y_pos, label = label),
      inherit.aes = FALSE,
      vjust = -0.4,
      size = 4
    ) +
    labs(
      title = paste0("Expression in ", group1, " vs ", group2),
      x = "Group",
      y = "Normalized expression"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      strip.text = element_text(face = "bold", size = 11),
      legend.position = "none"
    )
  
  #ggsave(
  #  filename = file.path(outdir, filename),
  #  plot = p,
  #  width = 10,
   # height = 5.5,
   # dpi = 300
  #)
  
  return(list(plot = p, pvals = pval_df))
}

# 602 vs 603
res_602_603 <- make_violin_compare_plot(
  seu = cluster1_sub,
  group_col = "pool_id",
  group1 = "602",
  group2 = "603",
  genes_vec = unname(genes_present),
  outdir = outdir,
  filename = "Violin_Ccr6_Klrb1a_602_vs_603.png"
)

# 602 vs 601
res_602_601 <- make_violin_compare_plot(
  seu = cluster1_sub,
  group_col = "pool_id",
  group1 = "602",
  group2 = "601",
  genes_vec = unname(genes_present),
  outdir = outdir,
  filename = "Violin_Ccr6_Klrb1a_602_vs_601.png"
)

print(res_602_603$pvals)
print(res_602_601$pvals)








# =================================
# Percent expressing (>0) for Ccr6 and Klrb1a
# Comparisons: 602 vs 601, and 602 vs 603
# =================================

library(Seurat)
library(dplyr)
library(ggplot2)

DefaultAssay(cluster1_sub) <- "RNA"

# use normalized data for consistency with your violin plots
cluster1_sub <- NormalizeData(dat, verbose = FALSE)

genes_requested <- c("Ccr6", "Klrb1a")

match_gene <- function(gene, gene_names) {
  hit <- gene_names[tolower(gene_names) == tolower(gene)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

genes_present <- sapply(genes_requested, match_gene, gene_names = rownames(dat))
names(genes_present) <- genes_requested

if (any(is.na(genes_present))) {
  missing_genes <- names(genes_present)[is.na(genes_present)]
  stop("These genes were not found in the object: ", paste(missing_genes, collapse = ", "))
}

make_percent_expressing_plot <- function(
    seu,
    group_col = "pool_id",
    group1,
    group2,
    genes_vec,
    outdir,
    filename_prefix
) {
  meta_col <- seu@meta.data[[group_col]]
  keep_cells <- rownames(seu@meta.data)[meta_col %in% c(group1, group2)]
  
  if (length(keep_cells) == 0) {
    stop("No cells found for comparison: ", group1, " vs ", group2)
  }
  
  sub <- subset(seu, cells = keep_cells)
  sub[[group_col]] <- factor(sub[[group_col]][,1], levels = c(group1, group2))
  
  expr_df <- FetchData(sub, vars = c(group_col, genes_vec))
  colnames(expr_df)[1] <- "group"
  
  # long format
  long_df <- do.call(
    rbind,
    lapply(genes_vec, function(g) {
      data.frame(
        group = expr_df$group,
        gene = g,
        expression = expr_df[[g]],
        expressed = expr_df[[g]] > 0,
        stringsAsFactors = FALSE
      )
    })
  )
  
  long_df$group <- factor(long_df$group, levels = c(group1, group2))
  long_df$gene  <- factor(long_df$gene, levels = genes_vec)
  
  # percent expressing summary
  pct_df <- long_df %>%
    group_by(gene, group) %>%
    summarise(
      n_cells = n(),
      n_expressing = sum(expressed, na.rm = TRUE),
      pct_expressing = 100 * n_expressing / n_cells,
      .groups = "drop"
    )
  
  # Fisher's exact test per gene
  pval_df <- do.call(
    rbind,
    lapply(genes_vec, function(g) {
      tmp <- long_df %>% filter(gene == g)
      
      tab <- table(
        factor(tmp$group, levels = c(group1, group2)),
        factor(tmp$expressed, levels = c(FALSE, TRUE))
      )
      
      ft <- fisher.test(tab)
      
      ymax <- max((pct_df %>% filter(gene == g))$pct_expressing, na.rm = TRUE)
      
      data.frame(
        gene = g,
        p_value = ft$p.value,
        odds_ratio = unname(ft$estimate),
        label = paste0("p = ", formatC(ft$p.value, format = "e", digits = 2)),
        y_pos = ymax + 5,
        x = 1.5,
        stringsAsFactors = FALSE
      )
    })
  )
  
  p <- ggplot(pct_df, aes(x = group, y = pct_expressing, fill = group)) +
    geom_col(width = 0.7, alpha = 0.85) +
    facet_wrap(~ gene, ncol = 2, scales = "free_y") +
    geom_segment(
      data = pval_df,
      aes(x = 1, xend = 2, y = y_pos, yend = y_pos),
      inherit.aes = FALSE
    ) +
    geom_text(
      data = pval_df,
      aes(x = x, y = y_pos, label = label),
      inherit.aes = FALSE,
      vjust = -0.4,
      size = 4
    ) +
    labs(
      title = paste0("Percent of cells expressing (>0): ", group1, " vs ", group2),
      x = "Group",
      y = "Percent expressing"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      strip.text = element_text(face = "bold", size = 11),
      legend.position = "none"
    )
  
  ggsave(
    filename = file.path(outdir, paste0(filename_prefix, ".png")),
    plot = p,
    width = 10,
    height = 5.5,
    dpi = 300
  )
  
  write.csv(
    pct_df,
    file = file.path(outdir, paste0(filename_prefix, "_percent_expressing_table.csv")),
    row.names = FALSE
  )
  
  write.csv(
    pval_df,
    file = file.path(outdir, paste0(filename_prefix, "_fisher_pvalues.csv")),
    row.names = FALSE
  )
  
  return(list(
    summary = pct_df,
    stats = pval_df,
    plot = p
  ))
}

# 602 vs 601
pct_602_601 <- make_percent_expressing_plot(
  seu = dat,
  group_col = "pool_id",
  group1 = "602",
  group2 = "601",
  genes_vec = unname(genes_present),
  outdir = outdir,
  filename_prefix = "Cell_Level_Percent_Expressing_Ccr6_Klrb1a_602_vs_601"
)

# 602 vs 603
pct_602_603 <- make_percent_expressing_plot(
  seu = dat,
  group_col = "pool_id",
  group1 = "602",
  group2 = "603",
  genes_vec = unname(genes_present),
  outdir = outdir,
  filename_prefix = "Cell_level_Percent_Expressing_Ccr6_Klrb1a_602_vs_603"
)

print(pct_602_601$summary)
print(pct_602_601$stats)

print(pct_602_603$summary)
print(pct_602_603$stats)



# Donor level Fischers Exact test


# =================================
# Sample-level percent expressing (>0)
# Ccr6 and Klrb1a
# =================================

library(Seurat)
library(dplyr)
library(ggplot2)

DefaultAssay(dat) <- "RNA"
dat <- NormalizeData(dat, verbose = FALSE)

# required metadata
required_cols <- c("orig.ident", "pool_id")
missing_cols <- setdiff(required_cols, colnames(dat@meta.data))
if (length(missing_cols) > 0) {
  stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "))
}

genes_requested <- c("Ccr6", "Klrb1a")

match_gene <- function(gene, gene_names) {
  hit <- gene_names[tolower(gene_names) == tolower(gene)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

genes_present <- sapply(genes_requested, match_gene, gene_names = rownames(dat))
names(genes_present) <- genes_requested

if (any(is.na(genes_present))) {
  missing_genes <- names(genes_present)[is.na(genes_present)]
  stop("These genes were not found in the object: ", paste(missing_genes, collapse = ", "))
}

make_sample_level_percent_plot <- function(
    seu,
    group1,
    group2,
    genes_vec,
    sample_col = "orig.ident",
    group_col = "pool_id",
    outdir,
    filename_prefix
) {
  # subset cells for the two groups
  keep_cells <- rownames(seu@meta.data)[seu@meta.data[[group_col]] %in% c(group1, group2)]
  
  if (length(keep_cells) == 0) {
    stop("No cells found for comparison: ", group1, " vs ", group2)
  }
  
  sub <- subset(seu, cells = keep_cells)
  
  expr_df <- FetchData(sub, vars = c(sample_col, group_col, genes_vec))
  colnames(expr_df)[1:2] <- c("sample_id", "group")
  
  expr_df$group <- factor(expr_df$group, levels = c(group1, group2))
  
  # long format
  long_df <- do.call(
    rbind,
    lapply(genes_vec, function(g) {
      data.frame(
        sample_id = expr_df$sample_id,
        group = expr_df$group,
        gene = g,
        expression = expr_df[[g]],
        expressed = expr_df[[g]] > 0,
        stringsAsFactors = FALSE
      )
    })
  )
  
  long_df$group <- factor(long_df$group, levels = c(group1, group2))
  long_df$gene <- factor(long_df$gene, levels = genes_vec)
  
  # sample-level percent expressing
  sample_pct_df <- long_df %>%
    group_by(gene, group, sample_id) %>%
    summarise(
      n_cells = n(),
      n_expressing = sum(expressed, na.rm = TRUE),
      pct_expressing = 100 * n_expressing / n_cells,
      .groups = "drop"
    )
  
  # group summary table
  summary_df <- sample_pct_df %>%
    group_by(gene, group) %>%
    summarise(
      n_samples = n(),
      mean_pct = mean(pct_expressing, na.rm = TRUE),
      median_pct = median(pct_expressing, na.rm = TRUE),
      sd_pct = sd(pct_expressing, na.rm = TRUE),
      min_pct = min(pct_expressing, na.rm = TRUE),
      max_pct = max(pct_expressing, na.rm = TRUE),
      .groups = "drop"
    )
  
  # sample-level Wilcoxon test per gene
  pval_df <- do.call(
    rbind,
    lapply(genes_vec, function(g) {
      tmp <- sample_pct_df %>% filter(gene == g)
      
      # only test if both groups have at least 2 samples
      nsamp <- tmp %>%
        group_by(group) %>%
        summarise(n = n_distinct(sample_id), .groups = "drop")
      
      if (nrow(nsamp) < 2 || any(nsamp$n < 2)) {
        pval <- NA_real_
        label <- "p = NA"
      } else {
        wt <- wilcox.test(pct_expressing ~ group, data = tmp, exact = FALSE)
        pval <- wt$p.value
        label <- paste0("p = ", formatC(pval, format = "e", digits = 2))
      }
      
      ymax <- max(tmp$pct_expressing, na.rm = TRUE)
      
      data.frame(
        gene = g,
        p_value = pval,
        label = label,
        y_pos = ymax + max(3, 0.08 * ymax),
        x = 1.5,
        stringsAsFactors = FALSE
      )
    })
  )
  
  p <- ggplot(sample_pct_df, aes(x = group, y = pct_expressing, fill = group)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.12, size = 2.2, alpha = 0.8) +
    facet_wrap(~ gene, ncol = 2, scales = "free_y") +
    geom_segment(
      data = pval_df,
      aes(x = 1, xend = 2, y = y_pos, yend = y_pos),
      inherit.aes = FALSE
    ) +
    geom_text(
      data = pval_df,
      aes(x = x, y = y_pos, label = label),
      inherit.aes = FALSE,
      vjust = -0.35,
      size = 4
    ) +
    labs(
      title = paste0("Sample-level percent expressing (>0): ", group1, " vs ", group2),
      x = "Group",
      y = "Percent expressing per sample"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      strip.text = element_text(face = "bold", size = 11),
      legend.position = "none"
    )
  
  ggsave(
    filename = file.path(outdir, paste0(filename_prefix, ".png")),
    plot = p,
    width = 10,
    height = 5.5,
    dpi = 300
  )
  
  write.csv(
    sample_pct_df,
    file = file.path(outdir, paste0(filename_prefix, "_sample_level_percentages.csv")),
    row.names = FALSE
  )
  
  write.csv(
    summary_df,
    file = file.path(outdir, paste0(filename_prefix, "_summary.csv")),
    row.names = FALSE
  )
  
  write.csv(
    pval_df,
    file = file.path(outdir, paste0(filename_prefix, "_wilcox_pvalues.csv")),
    row.names = FALSE
  )
  
  return(list(
    sample_level = sample_pct_df,
    summary = summary_df,
    stats = pval_df,
    plot = p
  ))
}

# 602 vs 601
sample_pct_602_601 <- make_sample_level_percent_plot(
  seu = dat,
  group1 = "602",
  group2 = "601",
  genes_vec = unname(genes_present),
  sample_col = "orig.ident",
  group_col = "pool_id",
  outdir = outdir,
  filename_prefix = "SampleLevel_Percent_Expressing_Ccr6_Klrb1a_602_vs_601"
)

# 602 vs 603
sample_pct_602_603 <- make_sample_level_percent_plot(
  seu = dat,
  group1 = "602",
  group2 = "603",
  genes_vec = unname(genes_present),
  sample_col = "orig.ident",
  group_col = "pool_id",
  outdir = outdir,
  filename_prefix = "SampleLevel_Percent_Expressing_Ccr6_Klrb1a_602_vs_603"
)

print(sample_pct_602_601$summary)
print(sample_pct_602_601$stats)

print(sample_pct_602_603$summary)
print(sample_pct_602_603$stats)


print("Script Complete")




th17_dot <-DotPlot(
  dat,
  features = c(
    "Rorc", "Rora", "Ccr6", "Il17a", "Il17f",
    "Il23r", "Ahr", "Stat3", "Batf", "Maf",
    "Il22", "Csf2", "Cxcr6"
  ),
  group.by = "seurat_clusters"
) + RotatedAxis()
#ggsave(plot=th17_dot, filename="th17_dot.png")







# filter object to just Th17 cells

dat <- cluster1_sub

outdir <- "/Users/eckco/Desktop/Kuhn_Lab/Jing_Data/MLN/Demultiplexed/Recluster_No_Hashtag/CD4_TCell_Analysis/Reclustered_CD4_Top20_DEGs_Per_Cluster/CD4_T_Pseudobulk_Volcano_By_demux_id_unpaired"

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
print(table(dat@meta.data[[group_col]], dat@meta.data[[replicate_col]]))


# Create pseudobulk IDs, which represent unique mouse-level pseudobulk samples.
# pool_id repeats across treatments, so the true unique mouse sample is: treatment + pool_id

dat$treatment <- as.character(dat@meta.data[[group_col]])
dat$mouse_within_treatment <- as.character(dat@meta.data[[replicate_col]])

dat$pseudobulk_id <- paste0(
  "treat.", dat$treatment,
  ".mouse.", dat$mouse_within_treatment
)

# Build pseudobulk metadata
pb_meta <- dat@meta.data %>%
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
  object = dat,
  assay = "RNA",
  layer = "counts"
)

cell_meta <- dat@meta.data[colnames(counts), , drop = FALSE]

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
#write.csv(
 # pb_meta,
 # file = file.path(outdir, "CD4_T_pseudobulk_metadata.csv"),
 # row.names = TRUE
#)

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
comparison_names <- unique(all_results$comparison)

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

#write.csv(
  #sig_results,
  #file = file.path(outdir, "CD4_T_significant_DE_genes_FDR0.05_logFC1.csv"),
  #row.names = FALSE
#)

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

