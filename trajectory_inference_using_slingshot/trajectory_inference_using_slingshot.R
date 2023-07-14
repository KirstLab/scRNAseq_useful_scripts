# Set a seed increase reproducibility
set.seed(1407)

# Load required packages
require(slingshot)
require(Seurat)
require(tidyverse)
require(paletteer)
require(magrittr)

# Help functions for plotting
theme_umap <- function(base.size = 14) {
    ggplot2::theme_classic(base_size = base.size) +
        ggplot2::theme(
            axis.ticks = ggplot2::element_blank(),
            axis.text = ggplot2::element_blank(),
            plot.subtitle = ggplot2::element_text(face = "italic", size = 11),
            plot.caption = ggplot2::element_text(face = "italic", size = 11)
        )
}
guide_umap <- function(key.size = 4) {
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = key.size, alpha = 1)))
}

palette_cluster <- paletteer::paletteer_d("ggsci::default_jama")
palette_celltype <- paletteer::paletteer_d("ggsci::category20_d3")
palette_heatmap <- paletteer::paletteer_d("wesanderson::Zissou1")

################################
## Parameters set by the user ##

SEURAT_OBJ = "your_seurat_obj.rds" # # name of the input seurat object
SLING_Embeddings = "PCA" # alternatively, UMA. This decides which embeddings to use for slingshot, the plot will always be in UMAP.
STARTING_CLUSTER = 2
OUT_PREFIX = "your_output_prefix" # prefix for the output files

################################

# Reads the RDS file 
seu <- readRDS(SEURAT_OBJ)
Idents(seu) <- seu$seurat_clusters

# Plot the seurat dataset, to confirm all is working
(p1 <- DimPlot(seu, group.by = "seurat_clusters") )

ggsave(
    paste0(OUT_PREFIX, "_input_clusters_UMAP.png"),
    p1,
    width = 10,
    height = 10,
    units = "in"
)

# trajectory inference
DefaultAssay(object = seu) <- "RNA"

sling_res <- slingshot::slingshot( Embeddings(seu, SLING_Embeddings), 
                                   clusterLabels = seu$seurat_clusters, 
                                   start.clus = STARTING_CLUSTER, 
                                   approx_points = 1000,
                                   extend = "n")

# Save the lineages
sink( paste0(OUT_PREFIX, "_lineages.txt") )
slingshot::slingLineages(sling_res)
sink()

# Save the RDS file with the lineages
write_rds(sling_res,
          paste0(OUT_PREFIX, "_slingshot.rds") )

# Check how many columns (lineages) exist and generate the colnames
sling_pt <- slingshot::slingPseudotime(sling_res) %>% 
    as.data.frame() 

col_n <- paste("PT", seq(1:ncol(sling_pt)), sep = "")
sling_pt <- sling_pt %>% 
    magrittr::set_colnames(col_n) 

col_n2 <- paste("sling_pt", col_n, sep = "_")
# Add the time as metadata to Seurat obj
seu <- AddMetaData(seu, 
                   metadata = sling_pt, 
                   col.name = col_n2)

# For each lineage, generate a new UMAP showing the pseudotime.
# Cell ordering based on slingshot
for (c in 1:ncol(sling_pt)) {
    
    pt_lin <- paste0("PT", c)
    
    ps_plot <- Embeddings(seu, "UMAP") %>% 
        as.data.frame() %>% 
        magrittr::set_colnames(c("UMAP_1", "UMAP_2")) %>% 
        mutate(PT = sling_pt[[pt_lin]]) %>% 
        ggplot(aes(x = UMAP_1, y = UMAP_2, color = PT)) + 
        geom_point(size = 1, alpha = 0.75) + 
        labs(x = "UMAP 1", 
             y = "UMAP 2", 
             color = paste("Pseudotime lineage",c )) + 
        scale_color_gradientn(colors = palette_heatmap, 
                              labels = scales::label_number(accuracy = 0.1),
                              na.value="white") + 
        theme_umap()
    
    ggsave(
        paste0( OUT_PREFIX,
                "_Pseudotime_lineage_",
                c, "_UMAP.png"),
        ps_plot, width = 12, height = 8, dpi = 300, bg = "white")
    
}

## Combining the pseudotime of all lineages in one plot
sling_pt[is.na(sling_pt)] <- 0
sling_pt <- sling_pt %>%
    keep(is.numeric) %>% 
    rowwise() %>%
    mutate(maxval = max(across()))

p3 <- Embeddings(seu, "UMAP") %>%
    as.data.frame() %>%
    magrittr::set_colnames(c("UMAP_1", "UMAP_2")) %>%
    mutate(PT = sling_pt$maxval) %>%
    ggplot(aes(x = UMAP_1, y = UMAP_2, color = PT)) +
    geom_point(size = 1, alpha = 0.75) +
    labs(x = "UMAP 1",
         y = "UMAP 2",
         color = "Pseudotime") +
    scale_color_gradientn(colors = palette_heatmap,
                          labels = scales::label_number(accuracy = 0.1)) +
    theme_umap()

ggsave(
    paste0(OUT_PREFIX, "_Pseudotime_combined_lineages_UMAP.png"),
    p3,
    width = 12,
    height = 8,
    dpi = 300,
    bg = "white")
