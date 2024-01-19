# single-cell analysis package
require(Seurat)
require(monocle3)
require(velocyto.R)
# plotting and data science packages
require(tidyverse)
require(cowplot)
require(patchwork)
require(SeuratWrappers)

#####################################
## Load the velocyto output - loom ##
ldat <- SeuratWrappers::ReadVelocity("data/velocyto/merged_samples.loom")
adata <- Seurat::as.Seurat(x = ldat)

## Read the rds file with the data from monocle 3
rds_name <- "batched_integrated_clustered_complete_dataset.rds"

my.cds <- readRDS( paste0("data/monocle3_files/", rds_name) )
colData(my.cds)$cell_type <- monocle3::clusters(my.cds)

plot_cells(my.cds,
           label_cell_groups = T,
           graph_label_size = 1.5,
           cell_size = 1,
           group_label_size = 7) + 
    facet_wrap(~timepoint)

########################
## Filtering velocyto ##

# keep only the cells in the monocle3 object. Some string modification is necessary to match the cell IDs in both datasets
colnames_v <- colnames(my.cds)
colnames_v <- colnames_v %>%
    gsub(x = ., pattern = "-1", "x") %>%
    gsub(x = ., pattern = "_1$", ":A17_sep_2022_0h_10k") %>%
    gsub(x = ., pattern = "_2$", ":A17_sep_2022_24h_10k") %>%
    gsub(x = ., pattern = "_3$", ":A17_sep_2022_48h_10k") %>%
    gsub(x = ., pattern = "_4$", ":A17_sep_2022_96h_10k") %>%
    gsub(x = ., pattern = "_5$", ":Sunn_sep_2022_0h_10k") %>%
    gsub(x = ., pattern = "_6$", ":Sunn_sep_2022_24h_10k") %>%
    gsub(x = ., pattern = "_7$", ":Sunn_sep_2022_48h_10k") %>%
    gsub(x = ., pattern = "_8$", ":Sunn_sep_2022_96h_10k")

cells_to_keep <- as_tibble(colnames_v) %>%
    tidyr::separate(col = value, into = c("sample", "barcode"), sep = ":") %>%
    dplyr::mutate( name = paste(barcode, sample, sep = ":") ) %>%
    pull("name")

# Select the cells in the same dataset
adata@meta.data$cells <- rownames(adata@meta.data)
adata <- subset(adata,
                cells %in% cells_to_keep) 

## transfer the counts from the monocle object
rna_counts <- monocle3::exprs(my.cds)

colnames(rna_counts) <- colnames(rna_counts) %>%
    gsub(x = ., pattern = "-1", "x") %>%
    gsub(x = ., pattern = "_1$", ":A17_sep_2022_0h_10k") %>%
    gsub(x = ., pattern = "_2$", ":A17_sep_2022_24h_10k") %>%
    gsub(x = ., pattern = "_3$", ":A17_sep_2022_48h_10k") %>%
    gsub(x = ., pattern = "_4$", ":A17_sep_2022_96h_10k") %>%
    gsub(x = ., pattern = "_5$", ":Sunn_sep_2022_0h_10k") %>%
    gsub(x = ., pattern = "_6$", ":Sunn_sep_2022_24h_10k") %>%
    gsub(x = ., pattern = "_7$", ":Sunn_sep_2022_48h_10k") %>%
    gsub(x = ., pattern = "_8$", ":Sunn_sep_2022_96h_10k") %>%
    as_tibble() %>%
    tidyr::separate(col = value, into = c("sample", "barcode"), sep = ":") %>%
    dplyr::mutate( name = paste(barcode, sample, sep = ":") ) %>%
    pull("name")

rna_assay <- CreateAssayObject(counts = rna_counts)
rna_assay@key <- "RNA_"
adata@assays$RNA <- rna_assay

# Add the cell embeddings for UMAP
mat_emb <- reducedDims(my.cds)$UMAP
rownames(mat_emb) <- rownames(mat_emb) %>%
    gsub(x = ., pattern = "-1", "x") %>%
    gsub(x = ., pattern = "_1$", ":A17_sep_2022_0h_10k") %>%
    gsub(x = ., pattern = "_2$", ":A17_sep_2022_24h_10k") %>%
    gsub(x = ., pattern = "_3$", ":A17_sep_2022_48h_10k") %>%
    gsub(x = ., pattern = "_4$", ":A17_sep_2022_96h_10k") %>%
    gsub(x = ., pattern = "_5$", ":Sunn_sep_2022_0h_10k") %>%
    gsub(x = ., pattern = "_6$", ":Sunn_sep_2022_24h_10k") %>%
    gsub(x = ., pattern = "_7$", ":Sunn_sep_2022_48h_10k") %>%
    gsub(x = ., pattern = "_8$", ":Sunn_sep_2022_96h_10k") %>%
    as_tibble() %>%
    tidyr::separate(col = value, into = c("sample", "barcode"), sep = ":") %>%
    dplyr::mutate( name = paste(barcode, sample, sep = ":") ) %>%
    pull("name")

mat_emb <- mat_emb[rownames(mat_emb) %in% colnames(adata),  ]
order_v <- colnames(adata@assays$RNA)
mat_emb <- mat_emb[order(match(rownames(mat_emb), order_v)), , drop = FALSE]

UMA_assay <- CreateDimReducObject(embeddings = mat_emb, key = 'UMAP_' )
UMA_assay@key <- "UMAP_"
adata@reductions$UMAP <- UMA_assay

## Adds embeding from PCA
# Add the cell embeddings for UMAP
mat_emb <- reducedDims(my.cds)$PCA
rownames(mat_emb) <- rownames(mat_emb) %>%
    gsub(x = ., pattern = "-1", "x") %>%
    gsub(x = ., pattern = "_1$", ":A17_sep_2022_0h_10k") %>%
    gsub(x = ., pattern = "_2$", ":A17_sep_2022_24h_10k") %>%
    gsub(x = ., pattern = "_3$", ":A17_sep_2022_48h_10k") %>%
    gsub(x = ., pattern = "_4$", ":A17_sep_2022_96h_10k") %>%
    gsub(x = ., pattern = "_5$", ":Sunn_sep_2022_0h_10k") %>%
    gsub(x = ., pattern = "_6$", ":Sunn_sep_2022_24h_10k") %>%
    gsub(x = ., pattern = "_7$", ":Sunn_sep_2022_48h_10k") %>%
    gsub(x = ., pattern = "_8$", ":Sunn_sep_2022_96h_10k") %>%
    as_tibble() %>%
    tidyr::separate(col = value, into = c("sample", "barcode"), sep = ":") %>%
    dplyr::mutate( name = paste(barcode, sample, sep = ":") ) %>%
    pull("name")

mat_emb <- mat_emb[rownames(mat_emb) %in% colnames(adata),  ]
order_v <- colnames(adata@assays$RNA)
mat_emb <- mat_emb[order(match(rownames(mat_emb), order_v)), , drop = FALSE]

pca_assay <- CreateDimReducObject(embeddings = mat_emb, key = 'pca_' )
pca_assay@key <- "pca_"
adata@reductions$pca <- pca_assay

## Add clustering from monocle3
cluster_mt <- as.data.frame( monocle3::clusters(my.cds) )

cluster_mt$cells <- rownames(cluster_mt) %>%
    gsub(x = ., pattern = "-1", "x") %>%
    gsub(x = ., pattern = "_1$", ":A17_sep_2022_0h_10k") %>%
    gsub(x = ., pattern = "_2$", ":A17_sep_2022_24h_10k") %>%
    gsub(x = ., pattern = "_3$", ":A17_sep_2022_48h_10k") %>%
    gsub(x = ., pattern = "_4$", ":A17_sep_2022_96h_10k") %>%
    gsub(x = ., pattern = "_5$", ":Sunn_sep_2022_0h_10k") %>%
    gsub(x = ., pattern = "_6$", ":Sunn_sep_2022_24h_10k") %>%
    gsub(x = ., pattern = "_7$", ":Sunn_sep_2022_48h_10k") %>%
    gsub(x = ., pattern = "_8$", ":Sunn_sep_2022_96h_10k") %>%
    as_tibble() %>%
    tidyr::separate(col = value, into = c("sample", "barcode"), sep = ":") %>%
    dplyr::mutate( name = paste(barcode, sample, sep = ":") ) %>%
    pull("name")

cluster_mt <- cluster_mt %>%
    dplyr::rename(cluster = "monocle3::clusters(my.cds)")

meta_mt <- plyr::join(adata@meta.data, cluster_mt)
adata@meta.data$seurat_clusters <- meta_mt$cluster

DimPlot(adata, reduction = 'UMAP', group.by = "seurat_clusters")

## Save the dataset in a format compatible to scVelo, to be used downstream.
rds_name_2 <- paste0("data/seurat_formated_", rds_name)
saveRDS(adata, rds_name_2)

# Run velocity 
adata <- RunVelocity(object = adata,
                     deltaT = 1,
                     kCells = 25,
                     fit.quantile = 0.02,
                     ncores = 23)

ident.colors <- (scales::hue_pal())(n = length(x = levels(x = adata@meta.data$seurat_clusters)))
names(x = ident.colors) <- levels(x = adata)

cell.colors <- ident.colors[adata@meta.data$seurat_clusters]
names(x = cell.colors) <- colnames(x = adata)

png(filename = "velocyto_outputs/velocityR_whole_dataset_plot.png",
    width = 12,
    height = 12,
    bg = "white",
    units = "cm",
    res = 300)

show.velocity.on.embedding.cor(
    emb = Embeddings(object = adata, reduction = "UMAP"),
    vel = Tool(object = adata, slot = "RunVelocity"),
    n = 200,
    scale = "sqrt",
    cell.colors = ac(x = cell.colors, alpha = 0.5), 
    cex = 0.8,
    arrow.scale = 3,
    show.grid.flow = TRUE,
    min.grid.cell.mass = 0.5,
    grid.n = 40,
    arrow.lwd = 1,
    do.par = FALSE,
    cell.border.alpha = 0.1,
    n.cores = 23)

dev.off()
