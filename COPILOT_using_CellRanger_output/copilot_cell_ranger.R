"Usage: copilot_cell_ranger.R (--sample_name <sn>) (--species <sp>) (--mito_string <ms>) (--mito_threshold <mt>) (--cloro_string <cs>) (--iterative_filtering <ifilt>) <input1>

Arguments:
    --sample_name    Sample name.
    --species   Species name.
    --mitho_string  String used to identify mitochondrial genes.
    --mito_threshold    % of expression that mitochondrial genes can provide. Used to define low quality cells.
    --cloro_string  String used to identify cloroplast genes.
    input1  Path to the raw matrices exported by Cell Ranger.
copilot_cell_ranger.R -h
" -> doc

set.seed(1407)

suppressMessages(require(docopt) )
suppressMessages(require(Seurat))
suppressMessages(require(Matrix))
suppressMessages(require(rjson))
suppressMessages(require(DoubletFinder))
suppressMessages(require(DropletUtils))
suppressMessages(require(ggplot2))
suppressMessages(require(scales))

# Loading functions from the copilot package
file.sources = list.files("COPILOT-master/R/", pattern="*.R",full.names = T)
sapply(file.sources,source,.GlobalEnv)

# retrieve the command-line arguments
opts <- docopt(doc)

seurat.data <- Seurat::Read10X( opts$input1,
  gene.column = 1,
)

outdir <- paste0("copilot_outputs/",
                 opts$sn,
                 "_filtered_matrix")
output_meta <- paste0("copilot_outputs/",
                      opts$sn,
                      "_metadata_of_quality_control")


if (opts$ifilt == TRUE) {
    
    filtering_ratio <- 0
} else {
    filtering_ratio <- 1
}

copilot_CR(
    sample.name = opts$sn,
    total.mtx = seurat.data,
    output_meta = output_meta,
    filtered.mtx.output.dir = outdir,
    species.name = opts$sp,
    mt.pattern = opts$ms, 
    mt.threshold = opts$mt,
    cp.pattern = opts$cs,
    top.percent = 0,
    estimate.doublet.rate = F,
    remove.doublet = F,
    do.seurat = F,
    do.annotation = F,
    filtering.ratio = filtering_ratio)