## Integration of multiple samples using seurat v5.
### Note: using Seurat package to integrate samples with methods such as harmony may not be the best approach. The wrap does not expose all parameters, so sometimes it will be preferrable to run the integration directly in harmony, for example. Also, scVI won't work without a GPU.

require(Seurat)
require(SeuratWrappers)
require(SeuratObject)
require(reticulate)
require(future)
require(batchelor)
require(harmony)

project_name <- "My_project" # Set the project name

# Set the number of cores that can be used for parallel processing
future::plan("multisession", workers = 6)

# Create a directory to save the output images
system("mkdir -p out_images")

folder_path <- "data" # Set the path to the folder containing the samples. Within this folder, you mush have one folder for each sample containing filtered dataset using the cell ranger output format (matrix.mtx.gz, barcodes.tsv.gz, and features.tsv.gz)

list_of_folders <- list.dirs(path = folder_path, full.names = FALSE, recursive = FALSE)

# The loop below will read each sample, create a Seurat object, and merge them into a single object.
rm(merged_obj) # Remove the merged object if it already exists
for( i in 1:length( list_of_folders ) ) {
    
    sample_mt <- Seurat::Read10X(
        data.dir = paste0(folder_path, "/", list_of_folders[i]) 
    )
    
    data_obj <- Seurat::CreateSeuratObject(
        counts = sample_mt,
        project = project_name,
        min.cells = 0, ## Adjust this according to your sample
        min.features = 0 ## Adjust this according to your sample
    )
    
    data_obj[["percent.mt"]] <- PercentageFeatureSet(data_obj, pattern = "MTg") # Change the pattern according to your sample. This is the tag to identify mitochondrial genes.
    
    data_obj_violin <- VlnPlot(data_obj,
                               features = c("nFeature_RNA",
                                            "nCount_RNA",
                                            "percent.mt"),
                               ncol = 3)
    
    ggplot2::ggsave(plot = data_obj_violin,
                    paste0("out_images/Violin_plot_", list_of_folders[1], ".png"),
                    width = 16, height = 12)
    
    data_obj$sample <- list_of_folders[i]
    
    if (i == 1) {
        merged_obj <- data_obj
        
    } else if (i > 1) {
        merged_obj <- merge(merged_obj, data_obj)
    }
}

# After merging the samples, we join the layers and normalize the data
merged_obj[['RNA']]  <- JoinLayers(object = merged_obj[['RNA']]  )
merged_obj[["RNA"]] <- split(merged_obj[["RNA"]], f = merged_obj$sample)

# Check the merged object
merged_obj

# Traditional Seurat pipeline for normalization and dimensionality reduction
merged_obj <- NormalizeData(merged_obj)
merged_obj <- FindVariableFeatures(merged_obj)
merged_obj <- ScaleData(merged_obj)
merged_obj <- RunPCA(merged_obj, npcs = 100)

# Visualize QC metrics as a violin plot
VlnPlot(merged_obj,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3)

# If needed, you can filter the merged dataset
merged_obj <- subset(merged_obj,
                     subset = nFeature_RNA > 100 &
                         nFeature_RNA < 15000 &
                         percent.mt < 2)

# Visualize QC metrics as a violin plot after filtering
violin_plot <- VlnPlot(merged_obj,
                       features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                       ncol = 3)

ggplot2::ggsave(plot = violin_plot,
                "out_images/Violin_plot_integrated_sample.png",
                width = 16, height = 12)

Elbow_plot <- ElbowPlot(merged_obj, ndims = 100)

ggplot2::ggsave(plot = Elbow_plot,
                "out_images/Elbow_plot_integrated_sample.png",
                width = 16, height = 12)

## Observing the data priour batch correction
merged_obj <- RunUMAP(merged_obj,
                      dims = 1:50,
                      reduction = "pca",
                      reduction.name = "umap.unintegrated")

batch_p1_comb <- DimPlot( merged_obj,
                                reduction = umap.unintegrated,
                                group.by = "treatment")

# Save the object before batch correction
saveRDS(merged_obj, ("seurat_processed_merged_obj_before_batch_correction.rds") )

###################################################################
## Integrating the data using different batch correction methods ##
###################################################################

## CCA
merged_obj <- IntegrateLayers(
    object = merged_obj,
    method = CCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.cca",
    verbose = FALSE
)

Seurat::DimPlot(merged_obj, reduction = "integrated.cca", group.by = "treatment")

## RPCA
merged_obj <- IntegrateLayers(
    object = merged_obj,
    method = RPCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.rpca",
    verbose = FALSE
)

Seurat::DimPlot(merged_obj, reduction = "integrated.rpca", group.by = "treatment")

## Harmony
merged_obj <- IntegrateLayers(
    object = merged_obj, method = HarmonyIntegration,
    orig.reduction = "pca", new.reduction = "integrated.harmony",
    verbose = FALSE
)

Seurat::DimPlot(merged_obj, reduction = "integrated.harmony", group.by = "treatment")

## FastMNN
merged_obj <- IntegrateLayers(
    object = merged_obj, method = FastMNNIntegration,
    new.reduction = "integrated.mnn",
    verbose = FALSE
)

Seurat::DimPlot(merged_obj, reduction = "integrated.mnn", group.by = "treatment")
