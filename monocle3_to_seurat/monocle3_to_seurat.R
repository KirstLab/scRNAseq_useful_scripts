# Set a seed increase reproducibility
set.seed(1407)

# Load required packages
require(Seurat)
require(monocle3)
require(tidyverse)

################################
## Parameters set by the user ##

MONOCLE3_OBJ <- "your_clustered_monocle3_obj.rds" # name of the input monocle3 object
SEURAT_OBJ <- "your_seurat_obj.rds" # name of the output seurat object

################################

# Read the RDS file
my.cds <- readRDS(MONOCLE3_OBJ)

# Add the cluster information to the monocle3 object
colData(my.cds)$cell_type <- monocle3::clusters(my.cds)

# Plot the monocle3 dataset, to confirm all is working
p1 <- plot_cells(my.cds,
           label_cell_groups = T,
           graph_label_size = 1.5,
           cell_size = 1,
           group_label_size = 7)

## Save the plot
ggsave("monocle3_plot.png", p1, width = 10, height = 10, units = "in")

# transfer the counts from the monocle object
rna_counts <- monocle3::exprs(my.cds)

## Create a new Seurat Object using the raw counts
adata <- CreateSeuratObject(counts = rna_counts,
                            assay = "RNA",
                            project = "scMedicago" # change this to your sample
                            )

# Merge the data normalized counts to the object

## Gather the normalized counts
rna_data <- monocle3::normalized_counts(my.cds)

## Merge the normalized counts to the object
rna_counts_assay <- CreateAssayObject(counts = rna_counts)
rna_counts_assay@key <- "RNA_"
adata@assays$RNA <- rna_counts_assay

# Check for "_" since Seurat doesn't accept it.
rownames(rna_data) <- sub(pattern = "_", "-", rownames(rna_data))
adata@assays$RNA@data <- rna_data

# Add the cell embeddings for UMAP
mat_emb <- reducedDims(my.cds)$UMAP
UMA_assay <- CreateDimReducObject(embeddings = mat_emb, key = 'UMAP_' )
UMA_assay@key <- "UMAP_"
adata@reductions$UMAP <- UMA_assay

# Add the cell embeddings for PCA
mat_emb <- reducedDims(my.cds)$PCA
pca_assay <- CreateDimReducObject(embeddings = mat_emb, key = 'pca_' )
pca_assay@key <- "pca_"
adata@reductions$pca <- pca_assay

# Add clustering and other meta information to the object
colData( my.cds )$seurat_clusters <- colData( my.cds )$cell_type
adata <- AddMetaData(adata, as.data.frame( colData( my.cds ) ) )

# Creates a UMAP plot to test if conversion worked
p2 <- Seurat::DimPlot(adata, reduction = 'UMAP',
                group.by = "seurat_clusters",
                pt.size = 1 )

## Save the plot
ggsave("seurat_plot.png", p2, width = 10, height = 10, units = "in")

# Save the seurat formated dataset as a RDS file
saveRDS(adata, SEURAT_OBJ)