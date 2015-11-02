#' Prepare matrix of expression values for downstream analyses
#' 
#' Takes ExpressionSet object and (optionally) 1) normalizes the expression 
#' levels across samples using the median of the geometric means method found 
#' within the DESeq package, 2) transforms the values in the matrix of 
#' expression values to a log2 scale, 3) sets expression values below the set 
#' limit of detection (LOD) to zero prior to the log2 transformation, and 4) 
#' provides an option to eliminate potential batch effects specified in a 
#' phenoData column using the ComBat method within the sva package.
#' 
#' @param cellData ExpressionSet object created with readCells.
#' @param LOD Limit of detection (numeric).  All expression values below this 
#'   limit will be set to zero prior to the log2 transformation.  LOD is best 
#'   determined by spiking in a dilution series of RNA standards with known 
#'   concentration.
#' @param norm Boolean specifying whether to normalize the expression values of 
#'   each gene across all samples in the ExpressionSet object using the median 
#'   of the geometric means.  This helps with making expression values more 
#'   comparable across samples. See explanation of the 
#'   estimateSizeFactorsForMatrix function within the DESeq package 
#'   documentation for more explanation.
#' @param batchGroup Character string specifying the name of the column header 
#'   in pData where the batch group information is stored.  All samples within a
#'   particular batch group should have the same designation.  If no character
#'   string is specified, batch is not run.  E.g. if there are four samples and
#'   two batch groups, the information in the batch group column could look like
#'   c('Batch1','Batch1','Batch2','Batch2').
#' @return ExpressionSet object modified according to the optional parameters 
#'   specified above.  Note that it is highly recommended that you at least 
#'   log-transform the data prior to any further downstream analysis.
#' @export



prepCells <- function(cellData, LOD = 1, norm = TRUE, batchGroup) {
    
    if (("prepCells" %in% colnames(pData(cellData)))) {
        stop("You have already run prepCells on this data", call. = FALSE)
    }
    
    data <- exprs(cellData)
    
    # Run DESeq geometric mean normalization
    if (norm == TRUE) {
        loggeomeans <- rowMeans(log(data))
        factors <- apply(data, 2, function(cnts) exp(median((log(cnts) - loggeomeans)[is.finite(loggeomeans)])))
        data <- t(t(data)/factors)
    }
    
    # Spit warning if LOD<1
    
    if (LOD < 1) {
        warning("Setting LOD<1 is not advised", call. = FALSE)
    }
    
    # Look for anything less than LOD and set it to 1.
    
    data[data < LOD] <- 1
    
    
    # Convert data to log2 scale (anything that was below LOD becomes 0)
    
    log.data <- log2(data)
    
    exprs(cellData) <- log.data
    
    
    if (!missing(batchGroup)) {
        
        warning("Removed genes with zero expression across all cells in order to identify batch effects", call. = FALSE)
        
        pheno <- pData(cellData)
        
        modcombat <- model.matrix(~1, data = pheno)
        
        batch_names = pheno[, batchGroup]
        
        expr.genes <- which(rowSums(log.data) > 0)
        
        log.data <- log.data[expr.genes, ]
        
        combat_edata <- ComBat(dat = log.data, batch = batch_names, mod = modcombat, par.prior = TRUE)
        
        log.data <- combat_edata
        
        fData(cellData)$BatchGenes <- row.names(fData(cellData)) %in% row.names(log.data)
        
        exprs(cellData) <- log.data
        
    }
    
    pData(cellData)$prepCells <- c("")
    
    pData(cellData)$prepCells[1] <- "Yes"
    
    cellData
    
} 
