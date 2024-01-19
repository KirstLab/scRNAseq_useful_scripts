#' copilot_CR
#' @param sample.name User defined sample name (character), which should be the same as the name of directory that contains spliced and unspliced matrices if you are following scKB pipeline to produce raw counts matrices.
#' @param spliced.mtx Gene by cell matrix of spliced counts, which should have column and row names, Default is NULL.
#' @param unspliced.mtx Gene by cell matrix of unspliced counts, which should have column and row names. Default is NULL.
#' @param total.mtx Gene by cell matrix of total counts, which should have column and row names. Default is NULL.
#' @param filtered.mtx.output.dir Output directory for quality filtered matrices. Default is NULL.
#' @param species.name Species name (character). Default is "Not Provided".
#' @param transcriptome.name Name of transcriptome annotation file. (e.g. TAIR10 for Arabidopsis). Default is "Not Provided".
#' @param sample.stats Meta data of the sample in data.frame format. Default is NULL.
#' @param mt.pattern Pattern of mitochondrial gene names/ids (character; e.g. "ATMG") or list of mitochondrial genes (character vector). This argument is required to run copilot.
#' @param mt.threshold Threshold of mitochondrial expression percentage. Cell would be treated as dying cell if it has mitochodrial expression percentage higher than this threshold (numeric). Default is 5.
#' @param cp.pattern Pattern of chloroplast gene names/ids (character; e.g. "ATCG") or list of chloroplast genes (character vector). Default is NULL.
#' @param top.percent Percentage of cells that contain high numer of UMIs filtered (numeric). Default is 1.
#' @param filtering.ratio Metric that controls the stringency of cell filtering (lenient: 1; strict:0; moderate: 0 < filtering.ratio < 1; numeric). Default is 1.
#' @param estimate.doublet.rate Whether or not to estimate doublet rate according to 10X Genomics' estimation (boolean). Default is TRUE.
#' @param doublet.rate User specified doublet rate (numeric). Default is NULL.
#' @param remove.doublet Whether or not to remove doublets after quality filtering of gene and cell (boolean). Default is TRUE.
#' @param do.seurat Whether or not to perform normalization, PCA, UMAP and clustering using Seurat and output a Seurat object (boolean). Default is TRUE.
#' @param do.annotation Whether or not to do annotation (boolean). COPILOT only supports annotation on root of Arabidopsis thaliana and Oryza sativa. Default is FALSE.
#' @param unwanted.genes Gene IDs/names of unwanted genes (character vector, e.g. cell cycle related genes, organelle genes ... etc). Default is NULL.
#' @param HVG Whether or not to select highly variable genes (boolean). Default is FALSE.
#' @param HVGN Number of highly variable genes selected (numeric). Defalut is 200.
#' @param dir_to_bulk Directory to reference expression profile for annotation. Default is NULL.
#' @param dir_to_color_scheme Directory to color scheme file for annotation. Default is NULL.
#' @param clustering_alg Algorithm for clustering (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm, which requires the leidenalg python; numeric). Default is 3.
#' @param res Resolution used for clustering (numeric). Default is 0.5.
#' @param min.UMI.low.quality Minimum UMIs for a barcode to be considered as cell (numeric). Default is 100.
#' @param min.UMI.high.quality Minimum UMIs for a cell to be considered as high quality cell (numeric). Default is 300.
#' @param legend.position x y position of the legend on UMI histogram plot (numeric vector of length 2). Default is c(0.8,0.8).
#' @export
#' @import Seurat
#' @import Matrix
#' @import scales
#' @import rjson
#' @import DoubletFinder
#' @import DropletUtils

copilot_CR <- function(sample.name, spliced.mtx = NULL, unspliced.mtx = NULL, total.mtx = NULL, filtered.mtx.output.dir = NULL,
                    output_meta = NULL, species.name = "Not Provided", transcriptome.name = "Not Provided", sample.stats = NULL, mt.pattern , mt.threshold = 5, cp.pattern = NULL, top.percent = 1, filtering.ratio = 1, estimate.doublet.rate = TRUE, doublet.rate = NULL, remove.doublet = TRUE,
                    do.seurat = TRUE, do.annotation = FALSE, unwanted.genes = NULL, HVG = FALSE, HVGN = 200, dir_to_bulk = NULL, dir_to_color_scheme = NULL, clustering_alg = 3, res = 0.5, min.UMI.low.quality = 100, min.UMI.high.quality = 300, legend.position = c(0.8,0.8)) {
    if (missing(sample.name)) {
        stop("Error! Necessary argument \"sample.name\" is missing.")
    }
    
    if (missing(mt.pattern)) {
        stop("Error! Necessary argument \"mt.pattern\" is missing.")
    }
    
    library(Seurat)
    library(Matrix)
    library(rjson)
    library(DoubletFinder)
    library(DropletUtils)
    library(ggplot2)
    library(scales)
    
    # Loading functions from the copilot package
    file.sources = list.files("COPILOT-master/R/", pattern="*.R", full.names = T)
    sapply(file.sources,source,.GlobalEnv)
    
    if ( is.null(spliced.mtx) && is.null(unspliced.mtx) && is.null(total.mtx) ) {
        
        spliced <- readMM(paste0("./", sample.name, "/spliced.mtx"))
        rownames(spliced) <- as.character(read.table(paste0("./", 
                                                            sample.name, "/spliced.barcodes.txt"), header = F)$V1)
        colnames(spliced) <- as.character(read.table(paste0("./", 
                                                            sample.name, "/spliced.genes.txt"), header = F)$V1)
        spliced <- t(spliced)
        unspliced <- readMM(paste0("./", sample.name, "/unspliced.mtx"))
        rownames(unspliced) <- as.character(read.table(paste0("./", 
                                                              sample.name, "/unspliced.barcodes.txt"), header = F)$V1)
        colnames(unspliced) <- as.character(read.table(paste0("./", 
                                                              sample.name, "/unspliced.genes.txt"), header = F)$V1)
        unspliced <- t(unspliced)
        bcs_use <- intersect(colnames(spliced), colnames(unspliced))
        tot_genes <- Matrix::rowSums(spliced)
        genes_use <- rownames(spliced)[tot_genes > 0]
        sf <- spliced[genes_use, bcs_use]
        uf <- unspliced[genes_use, bcs_use]
        afr <- sf + uf
        tot_gene <- Matrix::colSums(afr)
        bcs_use <- colnames(afr)[tot_gene > min.UMI.low.quality]
        sf <- spliced[genes_use, bcs_use]
        uf <- unspliced[genes_use, bcs_use]
        afr <- sf + uf
    } else if (!is.null(total.mtx)) {
        afr <- total.mtx
        tot_gene <- Matrix::colSums(afr)
        tot_genes <- Matrix::rowSums(afr)
        bcs_use <- colnames(afr)[tot_gene > min.UMI.low.quality]
        genes_use <- rownames(afr)[tot_genes > 0]
        afr <- afr[genes_use, bcs_use]
    } else {
        spliced <- spliced.mtx
        unspliced <- unspliced.mtx
        bcs_use <- intersect(colnames(spliced), colnames(unspliced))
        tot_genes <- Matrix::rowSums(spliced)
        genes_use <- rownames(spliced)[tot_genes > 0]
        sf <- spliced[genes_use, bcs_use]
        uf <- unspliced[genes_use, bcs_use]
        afr <- sf + uf
        tot_gene <- Matrix::colSums(afr)
        bcs_use <- colnames(afr)[tot_gene > min.UMI.low.quality]
        sf <- spliced[genes_use, bcs_use]
        uf <- unspliced[genes_use, bcs_use]
        afr <- sf + uf
    }
    
    nc <- colSums(afr)
    lnc <- log10(nc)
    pmt <- (colSums(afr[grep(paste(mt.pattern, collapse = "|"), 
                             rownames(afr)), ])/nc) * 100
    lnc.mt <- lnc[pmt >= mt.threshold]
    if (length(lnc.mt) == 0) {
        mmt <- log10(min.UMI.high.quality)
    } else if (length(lnc.mt) == 1) {
        mmt <- lnc.mt
    } else {
        dd <- density(lnc.mt)
        mmt <- max(dd$x[find_modes(dd$y)][dd$x[find_modes(dd$y)] < 
                                              max(boxplot.stats(lnc.mt)$stats)])
    }
    
    lcu <- as.numeric(log10(cumsum(rev(sort(nc)))))
    snc <- as.numeric(rev(sort(nc)))
    lin <- log10(seq(1, ncol(afr), 1))
    dYu <- diff(lcu)/diff(lin)
    dX <- rowMeans(embed(lin, 2))
    local_min <- min(dYu[1:(length(which(lnc > mmt)) - 1)])
    cnu <- which(dYu <= local_min)[1]
    if (snc[cnu] < min.UMI.high.quality) {
        cnu <- length(which(nc > min.UMI.high.quality))
    }
    ct <- sort(nc, decreasing = TRUE)[cnu]
    filter.UMI.thres.n = c(min.UMI.low.quality, ct)
    filter.UMI.thres = c(ct, 1e+06)
    nf <- afr[, (colSums(afr) > filter.UMI.thres.n[1]) & (colSums(afr) <= 
                                                              filter.UMI.thres.n[2])]
    af <- afr[, (colSums(afr) > filter.UMI.thres[1]) & (colSums(afr) < 
                                                            filter.UMI.thres[2])]
    nnf <- log1p((nf/colSums(nf)) * 10000)
    naf <- log1p((af/colSums(af)) * 10000)
    np <- rowSums(nnf)/ncol(nnf)
    sp <- rowSums(naf)/ncol(naf)
    cornp <- apply(naf, 2, cor, y = np)
    corsp <- apply(naf, 2, cor, y = sp)
    sidx <- c()
    nidx <- c()
    for (i in 1:ncol(af)) {
        if (corsp[i] > cornp[i]) {
            sidx <- c(sidx, i)
        } else {
            nidx <- c(nidx, i)
        }
    }
    nidx_thres <- length(nidx) * filtering.ratio
    print(paste0("threshold cell number: ", nidx_thres))
    n_iteration <- 0
    print(paste0("removed cells: ", length(nidx)))
    while (length(nidx) != 0 & length(nidx) >= nidx_thres) {
        n_iteration <- n_iteration + 1
        nf <- cbind(nf, afr[, colnames(af)[nidx]])
        af <- afr[, colnames(af)[sidx]]
        nnf <- log1p((nf/colSums(nf)) * 10000)
        naf <- log1p((af/colSums(af)) * 10000)
        np <- rowSums(nnf)/ncol(nnf)
        sp <- rowSums(naf)/ncol(naf)
        cornp <- apply(naf, 2, cor, y = np)
        corsp <- apply(naf, 2, cor, y = sp)
        sidx <- c()
        nidx <- c()
        for (i in 1:ncol(af)) {
            if (corsp[i] > cornp[i]) {
                sidx <- c(sidx, i)
            }
            else {
                nidx <- c(nidx, i)
            }
        }
        print(paste0("iteration: ", n_iteration))
        print(paste0("removed cells: ", length(nidx)))
    }
    
    # Avoid an error when user input a total counts matrix directly
    if ( is.null(spliced.mtx) && is.null(unspliced.mtx) && is.null(total.mtx) ) {
        sf <- sf[, colnames(af)]
        uf <- uf[, colnames(af)]
    }
    
    message("Iteration finished")
    nc <- colSums(af)
    ng <- colSums(af > 0)
    lnc <- log10(nc)
    lng <- log10(ng)
    pmts <- (colSums(af[grep(paste(mt.pattern, collapse = "|"), 
                             rownames(af)), ])/nc) * 100
    ncn <- colSums(nf)
    ngn <- colSums(nf > 0)
    lncn <- log10(ncn)
    lngn <- log10(ngn)
    pmtn <- (colSums(nf[grep(paste(mt.pattern, collapse = "|"), 
                             rownames(nf)), ])/ncn) * 100
    select <- rep(paste("high quality"), length(lnc))
    label_mt <- paste0("percent mt >= ", mt.threshold)
    select[which(pmts >= mt.threshold)] <- label_mt
    quantile_idx <- which(lnc < quantile(lnc[which(select == 
                                                       "high quality")], probs = (1 - top.percent/100)))
    ssidx <- intersect(quantile_idx, which(select == "high quality"))
    select[ssidx] <- paste("high quality (", length(ssidx), ")", 
                           sep = "")
    select[which(select == "high quality")] <- paste("top ", 
                                                     top.percent, "% (", length(which(select == "high quality")), 
                                                     ")", sep = "")
    nonselect <- rep("low quality", length(lncn))
    nonselect[which(pmtn >= mt.threshold)] <- label_mt
    select <- c(select, nonselect)
    select[which(select == label_mt)] <- paste("percent mt >= ", 
                                               mt.threshold, " (", length(which(select == label_mt)), 
                                               ")", sep = "")
    select[which(select == "low quality")] <- paste0("low quality (", 
                                                     length(which(select == "low quality")), ")")
    lncc <- c(lnc, lncn)
    lngc <- c(lng, lngn)
    lnc_gg <- data.frame(select = select, lnc = lncc)
    lnc_gg$rank <- match(lnc_gg$lnc, sort(lnc_gg$lnc, decreasing = TRUE))
    lnc_gg$lnc <- 10^lnc_gg$lnc
    lnc_gg$select <- factor(lnc_gg$select, levels = c(paste0("high quality (", 
                                                             length(grep("high", select)), ")"), paste0("low quality (", 
                                                                                                        length(grep("low", select)), ")"), paste0("percent mt >= ", 
                                                                                                                                                  mt.threshold, " (", length(grep("percent", select)), 
                                                                                                                                                  ")"), paste0("top ", top.percent, "% (", length(grep("top", 
                                                                                                                                                                                                       select)), ")")))
    lng_gg <- data.frame(select = select, lng = lngc)
    lng_gg$rank <- match(lng_gg$lng, sort(lng_gg$lng, decreasing = TRUE))
    lng_gg$lng <- 10^lng_gg$lng
    lng_gg$select <- factor(lng_gg$select, levels = c(paste0("high quality (", 
                                                             length(grep("high", select)), ")"), paste0("low quality (", 
                                                                                                        length(grep("low", select)), ")"), paste0("percent mt >= ", 
                                                                                                                                                  mt.threshold, " (", length(grep("percent", select)), 
                                                                                                                                                  ")"), paste0("top ", top.percent, "% (", length(grep("top", 
                                                                                                                                                                                                       select)), ")")))
    af <- af[, ssidx]
    # Avoid an error when user input a total counts matrix directly
    if ( is.null(spliced.mtx) && is.null(unspliced.mtx) && is.null(total.mtx) ) {
        sf <- sf[, ssidx]
        uf <- uf[, ssidx]
    }
    
    if (is.null(total.mtx)) {
        kb_stats <- c(fromJSON(file = paste0("./", sample.name, 
                                             "/inspect.json")), fromJSON(file = paste0("./", sample.name, 
                                                                                       "/run_info.json")))
        tech <- tech <- strsplit(kb_stats$call, "\\s")[[1]][8]
        seq_stats <- data.frame(stat = c("Number of Reads Processed", 
                                         "Reads Pseudoaligned", "Reads on Whitelist", "Total UMI Counts", 
                                         "Sequencing Technology", "Species", "Transcriptome"), 
                                value = prettyNum(c(kb_stats$n_processed, paste0(kb_stats$p_pseudoaligned, 
                                                                                 " %"), paste0(round(kb_stats$percentageReadsOnWhitelist, 
                                                                                                     2), " %"), sum(afr), tech, species.name, transcriptome.name), 
                                                  big.mark = ","))
    } else {
        seq_stats <- data.frame(stat = c("Number of Reads Processed", 
                                         "Reads Pseudoaligned", "Reads on Whitelist", "Total UMI Counts", 
                                         "Sequencing Technology", "Species", "Transcriptome"), 
                                value = prettyNum(c(NA, NA, NA, 
                                                    sum(afr), NA, species.name, transcriptome.name), 
                                                  big.mark = ","))
    }
    
    if (estimate.doublet.rate) {
        p_cnts_in_cells <- round((sum(af)/sum(afr)) * 100, 2)
        med_cnts_cell <- median(colSums(af))
        med_genes_cell <- median(apply(af, 2, function(x) sum(x >= 
                                                                  1)))
        tot_genes_detected <- sum(rowSums(af) >= 1)
        if (length(table(lnc_gg$select)[grep("mt", names(table(lnc_gg$select)))]) == 
            0) {
            p_cnts_dying <- 0
        } else {
            p_cnts_dying <- as.numeric(round(table(lnc_gg$select)[grep("mt", 
                                                                       names(table(lnc_gg$select)))]/nrow(lnc_gg) * 100, 2))
        }
        p_cnts_high <- as.numeric(round(table(lnc_gg$select)[grep("high", 
                                                                  names(table(lnc_gg$select)))]/nrow(lnc_gg) * 100, 
                                        2))
        doublet.rate <- est.doub.rate(ncol(af))
        cell_stats <- data.frame(stat = c("Estimated Number of High Quality Cell", 
                                          "High Quality Cell", "Total UMI Counts in High Quality Cell", 
                                          "UMI Counts in High Quality Cell", "Median UMI Counts per High Quality Cell", 
                                          "Median Genes per High Quality Cell", "Total Genes Detected in High Quality Cell", 
                                          "Cell above Mitochondrial Expression Threshold", 
                                          "Estimated Doublet Rate in High Quality Cell"), value = prettyNum(c(ncol(af), 
                                                                                                              paste0(p_cnts_high, " %"), sum(af), paste0(p_cnts_in_cells, 
                                                                                                                                                         " %"), med_cnts_cell, med_genes_cell, tot_genes_detected, 
                                                                                                              paste0(p_cnts_dying, " %"), paste0(doublet.rate, 
                                                                                                                                                 " %")), big.mark = ","))
    } else {
        p_cnts_in_cells <- round((sum(af)/sum(afr)) * 100, 2)
        med_cnts_cell <- median(colSums(af))
        med_genes_cell <- median(apply(af, 2, function(x) sum(x >= 
                                                                  1)))
        tot_genes_detected <- sum(rowSums(af) >= 1)
        if (length(table(lnc_gg$select)[grep("mt", names(table(lnc_gg$select)))]) == 
            0) {
            p_cnts_dying <- 0
        } else {
            p_cnts_dying <- as.numeric(round(table(lnc_gg$select)[grep("mt", 
                                                                       names(table(lnc_gg$select)))]/nrow(lnc_gg) * 
                                                 100, 2))
        }
        p_cnts_high <- as.numeric(round(table(lnc_gg$select)[grep("high", 
                                                                  names(table(lnc_gg$select)))]/nrow(lnc_gg) * 100, 
                                        2))
        cell_stats <- data.frame(stat = c("Estimated Number of High Quality Cell", 
                                          "High Quality Cell", "Total UMI Counts in High Quality Cell", 
                                          "UMI Counts in High Quality Cell", "Median UMI Counts per High Quality Cell", 
                                          "Median Genes per High Quality Cell", "Total Genes Detected in High Quality Cell", 
                                          "Cell above Mitochondrial Expression Threshold"), 
                                 value = prettyNum(c(ncol(af), paste0(p_cnts_high, 
                                                                      " %"), sum(af), paste0(p_cnts_in_cells, " %"), 
                                                     med_cnts_cell, med_genes_cell, tot_genes_detected, 
                                                     paste0(p_cnts_dying, " %")), big.mark = ","))
    }
    if (remove.doublet) {
        if (estimate.doublet.rate) {
            parameters <- data.frame(para = c("Iteration of Filtering", 
                                              "Mitochondrial Expression Threshold", "Top High Quality Cell Filtered", 
                                              "Doublet Removed"), value = c(n_iteration, paste0(mt.threshold, 
                                                                                                " %"), paste0(top.percent, " %"), "Yes"))
        }
        else {
            parameters <- data.frame(para = c("Iteration of Filtering", 
                                              "Mitochondrial Expression Threshold", "Top High Quality Cell Filtered", 
                                              "Expected Doublet Rate in High Quality Cell", 
                                              "Doublet Removed"), value = c(n_iteration, paste0(mt.threshold, 
                                                                                                " %"), paste0(top.percent, " %"), paste0(doublet.rate, 
                                                                                                                                         " %"), "Yes"))
        }
    } else {
        if (estimate.doublet.rate) {
            parameters <- data.frame(para = c("Iteration of Filtering", 
                                              "Mitochondrial Expression Threshold", "Top High Quality Cell Filtered", 
                                              "Doublet Removed"), value = c(n_iteration, paste0(mt.threshold, 
                                                                                                " %"), paste0(top.percent, " %"), "No"))
        }
        else {
            parameters <- data.frame(para = c("Iteration of Filtering", 
                                              "Mitochondrial Expression Threshold", "Top High Quality Cell Filtered", 
                                              "Expected Doublet Rate in High Quality Cell", 
                                              "Doublet Removed"), value = c(n_iteration, paste0(mt.threshold, 
                                                                                                " %"), paste0(top.percent, " %"), paste0(doublet.rate, 
                                                                                                                                         " %"), "No"))
        }
    }
    system( paste0("mkdir -p ",output_meta) )
    bc_rank_plot(lnc_gg = lnc_gg, save = paste0(output_meta, "/", 
                                                "/bc_rank_plot.png"))
    UMI_hist_plot(lnc_gg = lnc_gg, save = paste0(output_meta, "/",
                                                 "/UMI_hist_plot.png"), legend.position = legend.position)
    Gene_hist_plot(lng_gg = lng_gg, save = paste0(output_meta, "/",  
                                                  "/Gene_hist_plot.png"))
    if (med_cnts_cell < 1000) {
        do.annotation <- FALSE
    }
    
    system( paste0("mkdir -p ", filtered.mtx.output.dir) )
    if (!do.seurat) {
        print_HTML_no_analysis(parameters = parameters, cell_stats = cell_stats, 
                               seq_stats = seq_stats, sample_stats = sample.stats, 
                               dir = paste0(output_meta, "/"), sample.name = sample.name)
        if (is.null(filtered.mtx.output.dir)) {
            write10xCounts(paste0("./", sample.name), 
                           af, overwrite = T)
            if (is.null(total.mtx)) {
                write10xCounts(paste0("./", sample.name), 
                               sf, overwrite = T)
                write10xCounts(paste0("./", sample.name), 
                               uf, overwrite = T)
            }
        }
        else {
            write10xCounts(paste0(filtered.mtx.output.dir), 
                           af, overwrite = T)
            if (is.null(total.mtx)) {
                write10xCounts(paste0(filtered.mtx.output.dir), 
                               sf, overwrite = T)
                write10xCounts(paste0(filtered.mtx.output.dir), 
                               uf, overwrite = T)
            }
        }
    } else {
        if (!is.null(total.mtx)) {
            seu <- suppressWarnings(CreateSeuratObject(counts = af, 
                                                       assay = "RNA", project = sample.name))
        }
        else {
            seu <- suppressWarnings(CreateSeuratObject(counts = af, 
                                                       assay = "RNA", project = sample.name))
            seu[["spliced_RNA"]] <- CreateAssayObject(sf)
            seu[["unspliced_RNA"]] <- CreateAssayObject(uf)
            DefaultAssay(seu) <- "RNA"
        }
        seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = paste(mt.pattern, 
                                                                         collapse = "|"))
        seu[["percent.cp"]] <- PercentageFeatureSet(seu, pattern = paste(cp.pattern, 
                                                                         collapse = "|"))
        if (HVG) {
            vfn <- HVGN
        }
        else {
            vfn <- nrow(seu)
        }
        suppressWarnings(seu <- SCTransform(seu, variable.features.n = vfn, 
                                            assay = "RNA", new.assay.name = "SCT", verbose = FALSE))
        DefaultAssay(seu) <- "SCT"
        if (!is.null(mt.pattern) && !is.null(cp.pattern) && !is.null(unwanted.genes)) {
            use.genes <- rownames(seu)[-c(grep(paste(mt.pattern, 
                                                     collapse = "|"), rownames(seu)), grep(paste(cp.pattern, 
                                                                                                 collapse = "|"), rownames(seu)), sort(match(unwanted.genes, 
                                                                                                                                             rownames(seu))))]
        }
        else if (!is.null(mt.pattern) && !is.null(cp.pattern) && 
                 is.null(unwanted.genes)) {
            use.genes <- rownames(seu)[-c(grep(paste(mt.pattern, 
                                                     collapse = "|"), rownames(seu)), grep(paste(cp.pattern, 
                                                                                                 collapse = "|"), rownames(seu)))]
        }
        else if (!is.null(mt.pattern) && is.null(cp.pattern) && 
                 !is.null(unwanted.genes)) {
            use.genes <- rownames(seu)[-c(grep(paste(mt.pattern, 
                                                     collapse = "|"), rownames(seu)), sort(match(unwanted.genes, 
                                                                                                 rownames(seu))))]
        }
        else if (!is.null(mt.pattern) && is.null(cp.pattern) && 
                 is.null(unwanted.genes)) {
            use.genes <- rownames(seu)[-c(grep(paste(mt.pattern, 
                                                     collapse = "|"), rownames(seu)))]
        }
        else if (is.null(mt.pattern) && is.null(cp.pattern) && 
                 !is.null(unwanted.genes)) {
            use.genes <- rownames(seu)[-c(sort(match(unwanted.genes, 
                                                     rownames(seu))))]
        }
        else if (is.null(mt.pattern) && !is.null(cp.pattern) && 
                 is.null(unwanted.genes)) {
            use.genes <- rownames(seu)[-c(grep(paste(cp.pattern, 
                                                     collapse = "|"), rownames(seu)))]
        }
        else if (is.null(mt.pattern) && !is.null(cp.pattern) && 
                 !is.null(unwanted.genes)) {
            use.genes <- rownames(seu)[-c(grep(paste(mt.pattern, 
                                                     collapse = "|"), rownames(seu)), grep(paste(cp.pattern, 
                                                                                                 collapse = "|"), rownames(seu)))]
        }
        else {
            use.genes <- rownames(seu)
        }
        seu <- RunPCA(seu, verbose = FALSE, approx = FALSE, npcs = 10, 
                      features = use.genes)
        nExp_poi <- round((doublet.rate/100) * ncol(seu))
        suppressMessages(suppressWarnings(seu <- doubletFinder_v3(seu, 
                                                                  PCs = 1:10, pN = 0.25, pK = 0.15, nExp = nExp_poi, 
                                                                  reuse.pANN = FALSE, sct = TRUE)))
        if (remove.doublet) {
            seu <- subset(seu, cells = colnames(seu)[which(seu@meta.data[, 
                                                                         grep("DF.classifications", colnames(seu@meta.data))] == 
                                                               "Singlet")])
            suppressWarnings(seu <- SCTransform(seu, variable.features.n = vfn, 
                                                assay = "RNA", new.assay.name = "SCT", verbose = FALSE))
            suppressWarnings(seu <- SCTransform(seu, variable.features.n = vfn, 
                                                assay = "spliced_RNA", new.assay.name = "spliced_SCT", 
                                                verbose = FALSE))
            suppressWarnings(seu <- SCTransform(seu, variable.features.n = vfn, 
                                                assay = "unspliced_RNA", new.assay.name = "unspliced_SCT", 
                                                verbose = FALSE))
        }
        DefaultAssay(seu) <- "SCT"
        seu@misc$percent_mt_raw <- list(pmt)
        seu@misc$log_nCount_raw <- list(lnc_gg)
        seu@misc$log_nFeature_raw <- list(lng_gg)
        seu@misc$parameters <- list(parameters)
        if (do.annotation) {
            if (species.name == "Arabidopsis thaliana") {
                seu <- cor.anno.at(seu = seu, dir_to_bulk = dir_to_bulk, 
                                   unwanted.genes = unwanted.genes, clustering_alg = clustering_alg, 
                                   res = res, mt.pattern = mt.pattern, cp.pattern = cp.pattern)
                feature_plot_at(seu = seu, res = res, doublet.rate = doublet.rate, 
                                dir_to_color_scheme = dir_to_color_scheme, 
                                save = paste0("./", sample.name, "/feature_plot.png"))
                print_HTML(parameters = parameters, cell_stats = cell_stats, 
                           seq_stats = seq_stats, sample_stats = sample.stats, 
                           dir = paste0("./", sample.name), sample.name = sample.name)
            }
            else if (species.name == "Oryza sativa") {
                seu <- cor.anno.os(seu = seu, dir_to_bulk = dir_to_bulk, 
                                   unwanted.genes = unwanted.genes, clustering_alg = clustering_alg, 
                                   res = res, mt.pattern = mt.pattern, cp.pattern = cp.pattern)
                feature_plot_os(seu = seu, res = res, doublet.rate = doublet.rate, 
                                dir_to_color_scheme = dir_to_color_scheme, 
                                save = paste0("./", sample.name, "/feature_plot.png"))
                print_HTML(parameters = parameters, cell_stats = cell_stats, 
                           seq_stats = seq_stats, sample_stats = sample.stats, 
                           dir = paste0("./", sample.name), sample.name = sample.name)
            }
        }
        else {
            seu <- RunPCA(seu, verbose = FALSE, approx = FALSE, 
                          npcs = 50, features = use.genes)
            suppressMessages(suppressWarnings(seu <- RunUMAP(seu, 
                                                             reduction = "pca", dims = 1:50, umap.method = "umap-learn", 
                                                             metric = "correlation")))
            suppressMessages(suppressWarnings(seu <- FindNeighbors(seu, 
                                                                   reduction = "pca", dims = 1:50)))
            suppressMessages(suppressWarnings(seu <- FindClusters(seu, 
                                                                  resolution = res, algorithm = clustering_alg)))
            feature_plot(seu = seu, res = res, doublet.rate = doublet.rate, 
                         save = paste0("./", sample.name, "/feature_plot.png"))
            print_HTML(parameters = parameters, cell_stats = cell_stats, 
                       seq_stats = seq_stats, sample_stats = sample.stats, 
                       dir = paste0("./", sample.name), sample.name = sample.name)
        }
        if (is.null(filtered.mtx.output.dir)) {
            write10xCounts(paste0(filtered.mtx.output.dir), 
                           af, overwrite = T)
            if (is.null(total.mtx)) {
                write10xCounts(paste0(filtered.mtx.output.dir), 
                               sf, overwrite = T)
                write10xCounts(paste0(filtered.mtx.output.dir), 
                               uf, overwrite = T)
            }
        }
        else {
            write10xCounts(paste0(filtered.mtx.output.dir), 
                           af, overwrite = T)
            if (is.null(total.mtx)) {
                write10xCounts(paste0(filtered.mtx.output.dir), 
                               sf, overwrite = T)
                write10xCounts(paste0(filtered.mtx.output.dir), 
                               uf, overwrite = T)
            }
        }
        saveRDS(seu, file = paste0(paste0("./", sample.name, 
                                          "/"), sample.name, "_COPILOT.rds"))
    }
}