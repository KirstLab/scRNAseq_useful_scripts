# Trajectory inference with Slingshot

Trajectory inference is a method to infer the developmental trajectory of cells from single-cell RNA-seq data. The trajectory is inferred from the expression of genes that are differentially expressed across the cells. There are many software packages that can perform trajectory inference and no software is adequated for all datasets. Mainly, the choice of the software depends on the topology of the trajectory (e.g. linear or branched) and the number of cells in the dataset. For a review of trajectory inference methods, see [Saelens et al. (2019)](https://doi.org/10.1016/j.cell.2019.05.031).

The script stored in this folder performs trajectory inference using the [Slingshot](https://bioconductor.org/packages/release/bioc/html/slingshot.html) package, which has been shown to perform well in many scenarios. The script expects a Seurat object, including clusters and dimensionality reduction, as input. 

Before running the script, the user should set some parameters in the scripts, including the start and/or end of the expected trajectory. Optionally, the user can alter many other parameters controlling the behavior of Slingshot. To run the script, execute the following command in a terminal (or execute line by line in RStudio):

```bash
Rscript trajectory_inference_using_slingshot.R
``` 

Note that the script requires the following R packages:

slingshot - <https://bioconductor.org/packages/release/bioc/html/slingshot.html>
Seurat - <https://satijalab.org/seurat/>
tidyverse - <https://www.tidyverse.org/>
paletteer - <https://cran.r-project.org/web/packages/paletteer/index.html>
magrittr - <https://cran.r-project.org/web/packages/magrittr/index.html>
