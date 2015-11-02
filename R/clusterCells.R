#' Perform clustering analysis on the samples in the expression matrix.
#' 
#' Takes ExpressionSet object and clusters the samples using hierarchical,
#' k-means, or pam clustering.
#' 
#' @param cellData ExpressionSet object created with readCells (and preferably 
#'   transformed with prepCells).  It is also helpful to first run 
#'   reduceGenes_var and reduceGenes_pca.
#' @param k The number of desired clusters into which the samples are to be
#'   grouped.  The function calcGap can be used to obtain an unbiased estimate
#'   of the number of clusters.
#' @param methods Vector of character strings specifying which of the three
#'   clustering methods to perform.  All three can be specified, or a subset of
#'   the three.
#' @param hier_dist Character string specifying how to compute the distance
#'   matrix for 'hierarchical.'  Equivalent to the 'method' parameter within the
#'   dist function.
#' @param hier_clust Character string specifying the clustering method for
#'   'hierarchical.'  Equivalent to the 'method' parameter within the hclust
#'   function.
#' @return For each clustering method specified, a column is added to pData
#'   containing the cluster information for each sample (clusters are designated
#'   by roman numerals).

clusterCells <- function(cellData, k = 3, methods = c("hierarchical", "kmeans", "pam"), hier_dist = "euclidean", hier_clust = "ward") {
    
    if (!("prepCells" %in% colnames(pData(cellData)))) {
        warning("It would be wise to run prepCells prior to reduceGenes_pca", call. = FALSE)
    }
    
    x <- exprs(cellData)
    
    if ("hierarchical" %in% methods) {
        Hclust <- hclust(dist(t(x), method = hier_dist), method = hier_clust)
        tree <- cutree(Hclust, k)
        pData(cellData)$HC_Groups <- as.character(as.roman(tree))
    }
    
    if ("kmeans" %in% methods) {
        set.seed(25)
        km <- kmeans(t(x), centers = k, nstart = 25)$cluster
        pData(cellData)$KM_Groups <- as.character(as.roman(km))
    }
    
    if ("pam" %in% methods) {
        pm <- pam(t(x), k = k, cluster.only = TRUE)
        pData(cellData)$PAM_Groups <- as.character(as.roman(pm))
    }
    
    cellData
    
} 
