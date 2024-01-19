# COPILOT using CellRanger output

Quality control of single-cell RNA-seq is a challenging task, specially for samples other than mammalian cells. While Cell Ranger attempts to separate cells from background (empty droplets, dead cells, etc.), it is not always successful. In part, because the algorithm have implicit assumptions about the data that are not always met.

COPILOT (<https://github.com/Hsu-Che-Wei/COPILOT>) is a tool that can be used to perform quality control of single-cell RNA-seq data by using user-specified set of genes that represent the signal of noise or low quality cells for cell-filtering, such as the mitochondrial genes. 

In its original implementation, COPILOT does not work well with data from Cell Ranger. Therefore, we changed the code to make it compatible with Cell Ranger output. In addition, we used docopt to make the execution of the script easier. Moreover, the edit version of COPILOT is used only for removing the background cells, while the original version also allows the removal of doublets and the processing of the data using Seurat.

If using this script, please cite the original COPILOT papers:

* A single-cell Arabidopsis root atlas reveals developmental trajectories in wild-type and cell identity mutants: <https://doi.org/10.1016/j.devcel.2022.01.008>.

* Protocol for fast scRNA-seq raw data processing using scKB and non-arbitrary quality control with COPILOT:  <https://doi.org/10.1016/j.xpro.2022.101729>.


## Usage

First, download the script and the folder "COPILOT-master" from this repository. Then, run the script as follows:

```console
Rscript copilot_cell_ranger.R --sample_name <name of your sample> --species <name of the species> --mito_string <string to identify the mitochondrial genes> --mito_threshold <threshold to identify cells with high mitochondrial content> --cloro_string <string to identify the chloroplast genes> --iterative_filtering <TRUE or FALSE> <path to the folder with the RAW COUNTS from Cell Ranger output>
```

Note that the script was tested using R 4.2, and requires the following R packages:

* docopt: <https://cran.r-project.org/web/packages/docopt/index.html>
* Seurat: <https://satijalab.org/seurat/>
* Matrix: <https://cran.r-project.org/web/packages/Matrix/index.html>
* rjson: <https://cran.r-project.org/web/packages/rjson/index.html>
* DoubletFinder: <https://bioconductor.org/packages/release/bioc/html/DoubletFinder.html>
* DropletUtils: <https://bioconductor.org/packages/release/bioc/html/DropletUtils.html>
* ggplot2: <https://ggplot2.tidyverse.org/>
* scales: <https://cran.r-project.org/web/packages/scales/index.html>