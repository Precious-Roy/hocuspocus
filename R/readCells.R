#' Load all experiment data into R
#' 
#' Takes separate .txt files containing assay, pheno, feature, and experiment 
#' data and combines them into a single ExpressionSet object (uses Bioconductor 
#' ExpressionSet class format).
#' 
#' @param assay Character string specifying the path and filename of a 
#'   tab-delimited .txt file containing a matrix of expression values (e.g. 
#'   counts, FPKMs, TPMs, etc.) with cells in columns and genes in rows.  Be 
#'   sure to include cell and gene names as header and row names!
#' @param pheno Character string specifying the path and filename tab-delimited 
#'   .txt file containing information about the cells in separate columns.  Be 
#'   sure that the row names in this file match the column names in assayData 
#'   (can be a subset)
#' @param feature Character string specifying the path and filename of a 
#'   tab-delimited .txt file containing information about the genes in separate 
#'   columns.  Be sure that the row names in this file match the row names in 
#'   assayData (can be a subset)
#' @param experiment Character string specifying the path and filename of a 
#'   tab-delimited .txt file containing additional information about all the 
#'   column headers in the phenoData file in separate columns.  Be sure that the
#'   row names in this file match the column headers in the phenoData file.
#' @return ExpressionSet object containg matrix of expression values, phenoData,
#'   assayData, and experimentData.
#' @export


readCells <- function(assay, pheno, feature, experiment) {
    
    if (missing(assay)) {
        stop("You must specify an assay file.", call. = FALSE)
    } else {
        aData <- read.table(assay, header = TRUE, sep = "\t", row.names = 1, as.is = TRUE, check.names = FALSE)
    }
    
    aData <- as.matrix(aData)
    params <- list(assayData = aData)
    
    if (!missing(pheno)) {
        pData <- read.table(pheno, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
        aData <- aData[, row.names(pData)]
        pData <- new("AnnotatedDataFrame", data = pData)
        params["phenoData"] <- pData
        params$assayData <- aData
    } else {
        pData <- data.frame(colnames(aData))
        pData$GroupID <- "Group1"
        pData$Group_Colors <- "grey"
        pData <- new("AnnotatedDataFrame", data = pData)
        params["phenoData"] <- pData
    }
    
    if (!missing(feature)) {
        fData <- read.table(feature, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
        aData <- aData[row.names(fData), ]
        fData <- new("AnnotatedDataFrame", data = fData)
        params["featureData"] <- fData
        params$assayData <- aData
    } else {
        fData <- data.frame(row.names(aData))
        colnames(fData) <- "GeneID"
        row.names(fData) <- fData$GeneID
        fData <- new("AnnotatedDataFrame", data = fData)
        params["featureData"] <- fData
    }
    
    if (!missing(experiment)) {
        eData <- read.table(experiment, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
        eData <- do.call(MIAME, as.list(eData))
        params["experimentData"] <- eData
    }
    
    cellData <- do.call(ExpressionSet, params)
    cellData
    
} 
