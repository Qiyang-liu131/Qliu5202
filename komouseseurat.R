library(Seurat)

# Load your mouse dataset, the original dataset is available at
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6128372
KOmouse_testis_ko <- Read10X(data.dir = "/change/to/your/local/file/location/")
KOmouse_seurat <- CreateSeuratObject(counts = KOmouse_testis_ko)

# Add species information
KOmouse_seurat$species <- "KO_mouse"

# Normalize, scale, and perform PCA on KO mouse data
KOmouse_seurat <- NormalizeData(KOmouse_seurat)
KOmouse_seurat <- FindVariableFeatures(KOmouse_seurat)
KOmouse_seurat <- ScaleData(KOmouse_seurat)
KOmouse_seurat <- RunPCA(KOmouse_seurat)

# Assign cell type based on expression of biomarkers
KOmouse_seurat$stage <- "Unknown"
threshold <- 0.5  # Adjust this threshold as necessary
KOmouse_seurat$stage[which(KOmouse_seurat[["RNA"]]$data["Sycp3", ] > threshold)] <- "Spermatocytes"
KOmouse_seurat$stage[which(KOmouse_seurat[["RNA"]]$data["Cd46", ] > threshold)] <- "Spermatids"
KOmouse_seurat$stage[which(KOmouse_seurat[["RNA"]]$data["Stra8", ] > threshold)] <- "Spermatogonia"

# Check the classification
print(table(KOmouse_seurat$stage))

# Split the KO mouse data into separate Seurat objects based on stage
spermatogonia_rat <- subset(KOmouse_seurat, stage == "Spermatogonia")
spermatocytes_rat <- subset(KOmouse_seurat, stage == "Spermatocytes")
spermatids_rat <- subset(KOmouse_seurat, stage == "Spermatids")

# Merge the separate objects into one
merged_KOmouse_seurat <- merge(spermatogonia_rat, y = c(spermatocytes_rat, spermatids_rat))

# Normalize, scale, and perform PCA on merged object
merged_KOmouse_seurat <- NormalizeData(merged_KOmouse_seurat)
merged_KOmouse_seurat <- FindVariableFeatures(merged_KOmouse_seurat)
merged_KOmouse_seurat <- ScaleData(merged_KOmouse_seurat)
merged_KOmouse_seurat <- RunPCA(merged_KOmouse_seurat)

# Run UMAP on merged data
merged_KOmouse_seurat <- RunUMAP(merged_KOmouse_seurat, dims = 1:10)

# Plot UMAP showing all stages (Spermatogonia, Spermatocytes, Spermatids)
DimPlot(merged_KOmouse_seurat, reduction = "umap", group.by = "stage")

saveRDS(merged_KOmouse_seurat, "/change/to/your/local/file/location/")
