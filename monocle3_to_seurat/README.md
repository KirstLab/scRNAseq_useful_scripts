# Converting data from Monocle3 object to Seurat object

Many tools expect data in Seurat format. This script converts data from a Monocle3 object to a Seurat object. Note that it expects the Monocle3 object to be saved in a RDS file and the dataset to be fully processed (i.e. normalized, clustered, etc.).

The following information is transferred from the Monocle3 object to the Seurat object:

* Raw counts
* Normalized expression matrix
* Dimensionality reduction (PCA and UMAP)
* Clusters
* Other metadata (e.g. cell type, sample, etc.)

To run the script, first set the name of the Monocle3 object file and the name of the Seurat object file within the script. Then, run the script in a terminal (or execute line by line in RStudio).

```bash
Rscript monocle3_to_seurat.R
```

Note that the script requires the following R packages:

* Seurat - <https://satijalab.org/seurat/>
* monocle3 - <https://cole-trapnell-lab.github.io/monocle3/>
* tidyverse - <https://www.tidyverse.org/>
