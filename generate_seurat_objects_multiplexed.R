library(Seurat)
library(dplyr)
library(tibble)
library(stringr)

# Set parent path to pooled sample folders
base_dir <- "/scratch/alpine/ceck@xsede.org/Jan_5_2026_Backup/Kuhn_Data/Jing_Data/Backup/MLN_Data/Demultiplex/Output"
sample_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
print(sample_dirs)

# Lists to hold Seurat objects
seurat_list <- list()
pooled_list <- list()

# Function to process each demultiplexed sample within a pooled group
process_sample <- function(pool_path, inner_sample_name) {
  pool_name <- basename(pool_path)
  sample_name <- paste0("pool", pool_name, "_sample", inner_sample_name)
  message("Processing sample: ", sample_name)

  # Define output and matrix paths
  output_folder <- file.path(pool_path, paste0(pool_name, "_output"))
gex_path <- file.path(
  pool_path,
  paste0(pool_name, "_output"),
  "outs",
  "per_sample_outs",
  inner_sample_name,
  "count",
  "sample_filtered_feature_bc_matrix"
)
  vdj_path <- file.path(output_folder, "outs", "per_sample_outs", inner_sample_name, "vdj_t", "filtered_contig_annotations.csv")

  if (!dir.exists(gex_path)) stop("GEX path not found: ", gex_path)

  # Load GEX and ADT matrices
  gex <- Read10X(gex_path)

  if (is.list(gex) && "Gene Expression" %in% names(gex)) {
    rna <- gex[["Gene Expression"]]
  } else if (!is.list(gex)) {
    rna <- gex
  } else {
    stop("Gene Expression matrix not found")
  }

  adt <- if (is.list(gex) && "Antibody Capture" %in% names(gex)) gex[["Antibody Capture"]] else NULL

  # Create Seurat object
  seu <- CreateSeuratObject(counts = rna, project = sample_name)

  if (!is.null(adt)) {
    seu[["ADT"]] <- CreateAssayObject(counts = adt)
  }

  # Add metadata
  seu$pool_id <- pool_name
  seu$demux_id <- inner_sample_name
  seu$sample_id <- sample_name
  seu$orig.ident <- sample_name

  # Rename cells so barcodes stay unique after merge
  seu <- RenameCells(seu, add.cell.id = sample_name)

  # Load and merge VDJ if available
  if (file.exists(vdj_path)) {
    vdj <- read.csv(vdj_path)
    vdj <- vdj %>%
      filter(productive == "true") %>%
      group_by(barcode) %>%
      slice(1) %>%
      ungroup()

    # Rename VDJ barcodes to match renamed Seurat cell names
    vdj$barcode <- paste0(sample_name, "_", vdj$barcode)

    matching_barcodes <- vdj$barcode %in% colnames(seu)
    if (any(matching_barcodes)) {
      vdj <- vdj[matching_barcodes, ]
      rownames(vdj) <- vdj$barcode
      metadata <- seu@meta.data %>% rownames_to_column("barcode")
      metadata <- left_join(metadata, vdj, by = "barcode") %>% column_to_rownames("barcode")
      seu@meta.data <- metadata
    }
  }

  return(seu)
}

# Loop through all pooled groups
for (pool_path in sample_dirs) {
  pool_name <- basename(pool_path)

  output_folder <- file.path(pool_path, paste0(pool_name, "_output"))
  per_sample_dir <- file.path(output_folder, "outs", "per_sample_outs")

  if (!dir.exists(per_sample_dir)) {
    message("Skipping ", pool_path, ": per_sample_outs not found")
    next
  }

  inner_sample_dirs <- list.dirs(per_sample_dir, full.names = TRUE, recursive = FALSE)
  inner_sample_names <- basename(inner_sample_dirs)

  pool_seurat_objects <- list()

  # Loop through each demultiplexed sample within the pool
  for (inner_sample_name in inner_sample_names) {
    seurat_obj <- tryCatch({
      process_sample(pool_path, inner_sample_name)
    }, error = function(e) {
      message("Error processing ", pool_name, "/", inner_sample_name, ": ", e$message)
      return(NULL)
    })

    if (!is.null(seurat_obj)) {
      obj_name <- paste0("pool", pool_name, "_sample", inner_sample_name)
      seurat_list[[obj_name]] <- seurat_obj
      pool_seurat_objects[[obj_name]] <- seurat_obj
    }
  }

  # Merge all demultiplexed samples within this pooled group
  if (length(pool_seurat_objects) > 0) {
    pooled_merged <- pool_seurat_objects[[1]]
    if (length(pool_seurat_objects) > 1) {
      pooled_merged <- merge(
        x = pool_seurat_objects[[1]],
        y = pool_seurat_objects[-1],
        project = paste0("pool", pool_name, "_merged")
      )
    }
    pooled_list[[paste0("pool", pool_name, "_merged")]] <- pooled_merged
  }
}

# Save each mouse-level Seurat object
for (name in names(seurat_list)) {
  saveRDS(seurat_list[[name]], file = paste0("seurat_", name, ".rds"))
}

# Save each pooled-group merged Seurat object
for (name in names(pooled_list)) {
  saveRDS(pooled_list[[name]], file = paste0("seurat_", name, ".rds"))
}

message("Script complete.")