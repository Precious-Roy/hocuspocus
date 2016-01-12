#' Unbiased estimate of the number of cell or gene clusters using the gap
#' statistic.
#' 
#' Takes ExpressionSet object and calculates the optimal number of kmeans, pam,
#' or hierarchical clusters for the samples or genes using the gap statistic.
#' 
#' @param cellData ExpressionSet object created with readCells (and preferably 
#'   transformed with prepCells).  It is also helpful to first run 
#'   reduceGenes_var and reduceGenes_pca.
#' @param gene_clust Boolean specifying whether the gap statistic should be
#'   calculated for the samples or genes.  TRUE calculates for the cells, FALSE
#'   for the genes.
#' @param fun Character string specifying whether the gap statistic should be
#'   calculated for kmeans, pam, or hierarchical clustering.  Possible values
#'   are kmeans, pam, or hclust. clustering methods to perform.  All three can
#'   be specified, or a subset of the three.
#' @param max_clust Integer specifying the maximum possible number of clusters
#'   in the dataset.  Set higher than the expected value. matrix for
#'   'hierarchical.'  Equivalent to the 'method' parameter within the dist
#'   function.
#' @param boot Integer specifying the number of bootstrap iterations to perform
#'   when calculating the gap statistic. 'hierarchical.'  Equivalent to the
#'   'method' parameter within the hclust function.
#' @param plot Boolean specifying whether a plot of the gap values vs the number
#'   of clusters should be produced.
#' @param save Boolean specifying whether the plot should be saved.
#' @param print Boolean specifying whether the optimal number of clusters should
#'   be printed in the terminal window.
#' @return The optimal number of clusters calculated from the gap statistic with
#'   the given parameters.  A new column is added to pData indicating the
#'   optimal number of cell or gene clusters for the chosen clustering method.
#' @export

gapStats <- function(cellData, gene_clust = FALSE, fun = "kmeans", max_clust = 25, boot = 100, plot = TRUE, save = FALSE, print = TRUE) {
    
    hierarchical <- function(x, k) list(cluster = cutree(hclust(dist(x), method = "ward"), k = k))
    
    if (.Platform$OS.type == "windows") {
        quartz <- function() windows()
    }
    
    if (cellData@logData$prepCells[1] == "No") {
      warning("It would be wise to run prepCells prior to gapStats.", call. = FALSE)
    }
    
    exprs_matrix <- exprs(cellData)
    
    if (gene_clust == TRUE) {
        exprs_matrix <- t(exprs_matrix)
    }
    
    if (fun == "kmeans") {
        cg <- clusGap(t(exprs_matrix), FUNcluster = kmeans, K.max = max_clust, B = boot)
    }
    
    if (fun == "pam") {
        cg <- clusGap(t(exprs_matrix), FUNcluster = pam, K.max = max_clust, B = boot)
    }
    
    if (fun == "hclust") {
        cg <- clusGap(t(exprs_matrix), FUNcluster = hierarchical, K.max = max_clust, B = boot)
    }
    
    if (!(fun == "hclust" | fun == "pam" | fun == "kmeans")) {
        stop("fun must be one of kmeans, pam, or hiearchical.", call. = FALSE)
    }
    
    gap <- cg$Tab[, "gap"]
    
    gapSE <- cg$Tab[, "SE.sim"]
    
    k_opt <- maxSE(gap, gapSE)
    
    
    if (plot == TRUE) {
        quartz()
        plot(cg)
        if (save == TRUE && gene_clust == FALSE) {
            quartz(type = "pdf", file = "Gap Stats Plot - Cells.pdf")
            plot(cg)
            dev.off
        }
        if (save == TRUE && gene_clust == TRUE) {
            quartz(type = "pdf", file = "Gap Stats Plot - Genes.pdf")
            plot(cg)
            dev.off
        }
    }
    
    if (gene_clust == FALSE) {
        out <- paste("kOptC_", fun, sep = "")
        #pData(cellData)[, out] <- ""
        #pData(cellData)[1, out] <- k_opt
        if (print == TRUE) {
            print(paste("Optimal number of cell clusters =", k_opt))
        }
    }
    
    if (gene_clust == TRUE) {
        out <- paste("kOptG_", fun, sep = "")
        #pData(cellData)[, out] <- ""
        #pData(cellData)[1, out] <- k_opt
        if (print == TRUE) {
            print(paste("Optimal number of gene clusters =", k_opt))
        }
    }
    
    cellData
    
} 
