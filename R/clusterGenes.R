#' Perform clustering analysis on the genes in the expression matrix.
#' 
#' Takes ExpressionSet object and clusters the genes using hierarchical,
#' k-means, or pam clustering.
#' 
#' @param cellData ExpressionSet object created with readCells (and preferably 
#'   transformed with prepCells).  It is also helpful to first run 
#'   reduceGenes_var and reduceGenes_pca.
#' @param k The number of desired clusters into which the genes are to be
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
#' @export


clusterGenes <- function(cellData, k = 3, methods = c("hierarchical", "kmeans", "pam"), hier_dist = "euclidean", hier_clust = "ward") {
    
  if (cellData@logData$prepCells[1] == "No") {
    warning("It would be wise to run prepCells prior to clusterGenes.", call. = FALSE)
  }
    
    x <- exprs(cellData)
    
    if ("hierarchical" %in% methods) {
        Hclust <- hclust(dist(x, method = hier_dist), method = hier_clust)
        tree <- cutree(Hclust, k)
        fData(cellData)[row.names(x), "HC_Groups"] <- as.character(as.roman(tree))
    }
    
    if ("kmeans" %in% methods) {
        set.seed(25)
        km <- kmeans(x, centers = k, nstart = 25)$cluster
        fData(cellData)[row.names(x), "KM_Groups"] <- as.character(as.roman(km))
    }
    
    if ("pam" %in% methods) {
        pm <- pam(x, k = k, cluster.only = TRUE)
        fData(cellData)[row.names(x), "PAM_Groups"] <- as.character(as.roman(pm))
    }
    
    cellData
    
} 
