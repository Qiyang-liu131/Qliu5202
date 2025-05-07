# Load necessary libraries
library(Seurat)

# Load your mouse dataset, the original dataset is available at
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104556
mouse_full <- Read10X(data.dir = 
 "/change/to/your/local/file/location/")

# Create Seurat object
mouse_seurat <- CreateSeuratObject(counts = mouse_full)

# Step 1: Quality Control (QC)
# Calculate mitochondrial gene percentage 
#(adjust pattern for rat mitochondrial genes, e.g., "^mt-")
mouse_seurat[["percent.mt"]] <- 
  PercentageFeatureSet(mouse_seurat, pattern = "^mt-")

# Filter cells based on QC metrics
mouse_seurat <- subset(mouse_seurat, 
  subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

# Step 2: Normalization
mouse_seurat <- NormalizeData(mouse_seurat, 
  normalization.method = "LogNormalize", scale.factor = 10000)

# Step 3: Define the markers (biomarkers you provided)
markers_sertoli <- c("Clu")
markers_leydig <- c("Insl3")
markers_peritubular <- c("Acta2")
markers_macrophage <- c("Apoe")
markers_spermatogonia <- c("Gfra1")
markers_spermatocytes <- c("Sycp3")

# Spermatid stages
markers_early_spermatid <- c("Lrriq1")
markers_late_round_spermatid <- c("Acrv1")
markers_elongated_spermatid <- c("Prm1")

# Germ Cell
markers_all_germ_cells <- c("Ddx4", "Mybl1", "Lca5l", "Arrdc5", "Tex38")

# Combine all markers into one list
all_markers <- list(
  markers_sertoli, markers_leydig, markers_peritubular,
  markers_macrophage, markers_spermatogonia,
  markers_spermatocytes, markers_early_spermatid,
  markers_late_round_spermatid, markers_elongated_spermatid
)

# Now include the four markers to the list of all markers 
# without using them for clustering
all_markers_with_extra <- unlist(c(all_markers, markers_all_germ_cells))

# Step 4: Subset Seurat object based on the markers 
# (including the four extra markers)
mouse_seurat <- 
  mouse_seurat[rownames(mouse_seurat) %in% all_markers_with_extra, ]

# Step 5: Scaling the data, excluding the 4 genes from the list of
# features used in clustering
mouse_seurat <- ScaleData(mouse_seurat, 
  features = setdiff(rownames(mouse_seurat), markers_all_germ_cells))

# Step 6: PCA for dimensionality reduction
# excluding the 4 genes from the features
mouse_seurat <- RunPCA(mouse_seurat, features = setdiff(rownames(mouse_seurat), 
  markers_all_germ_cells))

# Step 7: Clustering
mouse_seurat <- FindNeighbors(mouse_seurat, dims = 1:8)
mouse_seurat <- FindClusters(mouse_seurat, resolution = 0.8)

# Step 8: UMAP for Visualization
mouse_seurat <- RunUMAP(mouse_seurat, dims = 1:8)

# Step 10: Find marker genes for each cluster
cluster_markers <- FindAllMarkers(mouse_seurat, only.pos = TRUE, 
  min.pct = 0.25, logfc.threshold = 0.25)

# Step 11: Annotate clusters based on marker genes
cell_type_markers <- list(
  Sertoli = markers_sertoli,
  Leydig = markers_leydig,
  Peritubular = markers_peritubular,
  Macrophage = markers_macrophage,
  Spermatogonia = markers_spermatogonia,
  Spermatocytes = markers_spermatocytes,
  Early_Spermatids = markers_early_spermatid,
  Late_Round_Spermatids = markers_late_round_spermatid,
  Elongated_Spermatids = markers_elongated_spermatid
)

# Function to annotate clusters
annotate_clusters <- function(cluster_markers, cell_type_markers) {
  cluster_annotations <- list()
  
  for (cluster in unique(cluster_markers$cluster)) {
    cluster_genes <- cluster_markers[cluster_markers$cluster == cluster, "gene"]
    
    # Calculate overlap between cluster genes and cell type markers
    overlaps <- sapply(cell_type_markers, function(markers) {
      length(intersect(cluster_genes, markers))
    })
    
    # Assign the cell type with the highest overlap
    cluster_annotations[[as.character(cluster)]] <- names(which.max(overlaps))
  }
  
  return(cluster_annotations)
}

# Annotate clusters
cluster_annotations <- annotate_clusters(cluster_markers, cell_type_markers)

# Ensure cluster_annotations is a named vector
cluster_annotations <- unlist(cluster_annotations)

# Create a metadata vector for each cell based on its cluster ID
annotations <- cluster_annotations[as.character(Idents(mouse_seurat))]

# Assign names to match the Seurat object cell barcodes
names(annotations) <- colnames(mouse_seurat)

# Add the metadata to the Seurat object
mouse_seurat <- AddMetaData(
  object = mouse_seurat,
  metadata = annotations,
  col.name = "cell_type"
)

# Step 12: Visualize the annotated clusters
DimPlot(mouse_seurat, group.by = "cell_type", label = TRUE)

# Step 13: Visualize expression of 'Mybl1' gene on annotated UMAP
FeaturePlot(mouse_seurat, features = "Mybl1", reduction = "umap", 
            pt.size = 1, label = TRUE)

# Step 14: Visualize expression of 'Lca5l' gene on annotated UMAP
FeaturePlot(mouse_seurat, features = "Lca5l", reduction = "umap", 
            pt.size = 1, label = TRUE)

# Step 15: Visualize expression of 'Arrdc5' gene on annotated UMAP
FeaturePlot(mouse_seurat, features = "Arrdc5", reduction = "umap", 
            pt.size = 1, label = TRUE)

# Step 16: Save the Seurat Object
saveRDS(mouse_seurat, file = "annotated_mouse_seurat_with_extra_genes.rds")

?saveRDS()