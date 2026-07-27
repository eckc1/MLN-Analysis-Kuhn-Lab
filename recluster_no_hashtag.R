library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)

# Ignore local RAM throttle 
mem.maxVSize(vsize = Inf)
# Load merged data
dat <- readRDS("/Users/eckco/Desktop/Kuhn_Lab/Jing_Data/MLN/Demultiplexed/Original_dat_clustered.rds")
head(dat)
Reductions(dat)
transcript <- DimPlot(dat,label=TRUE, reduction = "wnn.umap")
#ggsave(plot=transcript, filename="transcript_wnn.png")

demux_id <- DimPlot(dat,label=TRUE, reduction = "wnn.umap", group.by = "demux_id")
#ggsave(plot=demux_id, filename="demux_id_wnn.png")

sample_id <- DimPlot(dat,label=FALSE, reduction = "wnn.umap", group.by = "sample_id")
#ggsave(plot=sample_id, filename="sample_id_wnn.png")

pool_id <- DimPlot(dat,label=TRUE, reduction = "wnn.umap", group.by = "pool_id")
#ggsave(plot=pool_id, filename="pool_id_wnn.png")


# Check whether demux_id is driving clustering
Reductions(dat)

table(dat$seurat_clusters, dat$demux_id)
cluster_hashtag_counts <- prop.table(table(dat$seurat_clusters, dat$demux_id), margin = 1)
write.csv(cluster_hashtag_counts ,"cluster_hashtag_counts.csv")
#pc_embeddings <- Embeddings(dat, "pca")
pc_embeddings <- Embeddings(dat, "apca")
pc_df <- data.frame(
  PC1 = pc_embeddings[, 1],
  demux_id = dat$demux_id
)
#png(filename = "PC1_demux_id.png",width = 1800,height = 1200,res = 200)
#boxplot(PC1 ~ demux_id, data = pc_df, las = 2)
#dev.off()

pc_df <- as.data.frame(pc_embeddings[, 1:10])
pc_df$demux_id <- dat$demux_id

pvals <- sapply(colnames(pc_df)[1:10], function(pc) {
  summary(aov(pc_df[[pc]] ~ pc_df$demux_id))[[1]][["Pr(>F)"]][1]
})

#write.csv(pvals, "PCA_pvalues.csv")


pc_df <- as.data.frame(Embeddings(dat, "apca")[, 1:10])
pc_df$demux_id <- dat$demux_id

r2 <- sapply(colnames(pc_df)[1:10], function(pc) {
  fit <- lm(pc_df[[pc]] ~ pc_df$demux_id)
  summary(fit)$r.squared
})
#write.csv(r2,"r2_values.csv")

#png(filename = "Association_between_demux_id_and_ADT_PCs.png",width = 1800,height = 1200,res = 200)
#barplot(r2,las = 2,ylab = "Variance explained by demux_id",main = "Association between demux_id and ADT PCs")
#dev.off()

# Check loadings for suspect PCs
adt_loadings <- Loadings(dat[["apca"]])
for (pc in c(3, 4, 5, 8, 9)) {
  cat("\nTop positive loadings for ADT PC", pc, "\n")
  print(head(sort(adt_loadings[, pc], decreasing = TRUE), 15))
  cat("\nTop negative loadings for ADT PC", pc, "\n")
  print(head(sort(adt_loadings[, pc], decreasing = FALSE), 15))
}

# Remove hashtag features
DefaultAssay(dat) <- "ADT"
hto_features <- grep(
  "hash|hto|hashtag|totalseq|sample|demux|601|602|603",
  rownames(dat[["ADT"]]),value = TRUE,ignore.case = TRUE)
hto_features

adt_features_clean <- setdiff(rownames(dat[["ADT"]]), hto_features)
dat_nohto <- dat
dat_nohto[["ADT"]] <- subset(
  dat[["ADT"]],features = adt_features_clean)
# check removal was successful
rownames(dat[["ADT"]])
rownames(dat_nohto[["ADT"]])
hto_features %in% rownames(dat_nohto[["ADT"]])
dat <- dat_nohto 
setdiff(hto_features, rownames(dat))
#saveRDS(dat, "NoHashtag_dat.rds")

#############
# recluster full dataset
dat <- readRDS("NoHashtag_dat.rds")
head(dat)
Layers(dat)

# Cluster RNA Assay
DefaultAssay(dat) <- 'RNA'
dat<- NormalizeData(dat) %>% FindVariableFeatures() %>% 
  ScaleData() %>% RunPCA(reduction.name = 'pca')

# Cluster ADT Assay
DefaultAssay(dat) <- 'ADT'

VariableFeatures(dat) <- rownames(dat[["ADT"]])
dat <- NormalizeData(dat, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')


set.seed(4538)

# Batch correction with harmony could be implemented here

# Integrate multi-modal clustering
dat <- FindMultiModalNeighbors(
  dat,
  reduction.list = list("pca", "apca"),
  dims.list = list(1:15, 1:18),
  modality.weight.name = "wnn.weight",
  k.nn = 20
)

# UMAP + clustering
dat <- RunUMAP(dat, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
dat <- FindClusters(dat, graph.name = "wsnn", algorithm = 3, resolution = 0.2, verbose = FALSE)



# Keep clusters with >= 200 cells
clust_counts <- table(Idents(dat))
keep_clusters <- names(clust_counts)[clust_counts >= 200]
dat <- subset(dat, idents = keep_clusters)
DimPlot(dat)

# Process RNA assay
transcript <- DimPlot(dat,label=TRUE, reduction = "wnn.umap")
ggsave(plot=transcript, filename="transcript_wnn_noHashtag.png")

demux_id <- DimPlot(dat,label=FALSE, reduction = "wnn.umap", group.by = "demux_id")
ggsave(plot=demux_id, filename="demux_id_wnn_noHashtag.png")

sample_id <- DimPlot(dat,label=FALSE, reduction = "wnn.umap", group.by = "sample_id")
ggsave(plot=sample_id, filename="sample_id_wnn_noHashtag.png")

pool_id <- DimPlot(dat,label=FALSE, reduction = "wnn.umap", group.by = "pool_id")
ggsave(plot=pool_id, filename="pool_id_wnn_noHashtag.png")

#saveRDS(dat, "noHashtag_Reclustered.rds")


# Check whether demux_id is driving clustering
Reductions(dat)

table(dat$seurat_clusters, dat$demux_id)
cluster_hashtag_counts <- prop.table(table(dat$seurat_clusters, dat$demux_id), margin = 1)
write.csv(cluster_hashtag_counts ,"cluster_hashtag_counts.csv")
#pc_embeddings <- Embeddings(dat, "pca")
pc_embeddings <- Embeddings(dat, "apca")
pc_df <- data.frame(
  PC1 = pc_embeddings[, 1],
  demux_id = dat$demux_id
)
#png(filename = "PC1_demux_id.png",width = 1800,height = 1200,res = 200)
#boxplot(PC1 ~ demux_id, data = pc_df, las = 2)
#dev.off()

pc_df <- as.data.frame(pc_embeddings[, 1:10])
pc_df$demux_id <- dat$demux_id

pvals <- sapply(colnames(pc_df)[1:10], function(pc) {
  summary(aov(pc_df[[pc]] ~ pc_df$demux_id))[[1]][["Pr(>F)"]][1]
})

#write.csv(pvals, "PCA_pvalues.csv")


pc_df <- as.data.frame(Embeddings(dat, "apca")[, 1:10])
pc_df$demux_id <- dat$demux_id

r2 <- sapply(colnames(pc_df)[1:10], function(pc) {
  fit <- lm(pc_df[[pc]] ~ pc_df$demux_id)
  summary(fit)$r.squared
})
#write.csv(r2,"r2_values.csv")

#png(filename = "Association_between_demux_id_and_ADT_PCs.png",width = 1800,height = 1200,res = 200)
#barplot(r2,las = 2,ylab = "Variance explained by demux_id",main = "Association between demux_id and ADT PCs")
#dev.off()

# Check loadings for suspect PCs
adt_loadings <- Loadings(dat[["apca"]])
for (pc in c(3, 4, 5, 8, 9)) {
  cat("\nTop positive loadings for ADT PC", pc, "\n")
  print(head(sort(adt_loadings[, pc], decreasing = TRUE), 15))
  cat("\nTop negative loadings for ADT PC", pc, "\n")
  print(head(sort(adt_loadings[, pc], decreasing = FALSE), 15))
}