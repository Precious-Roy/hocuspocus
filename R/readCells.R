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


readCells <- function(assay, pheno, feature, experiment, color) {
    
    if (missing(assay)) {
        stop("You must specify an assay file.", call. = FALSE)
    } else {
        aData <- read.table(assay, header = TRUE, sep = "\t", row.names = 1, as.is = TRUE, check.names = FALSE)
    }
    
    aData <- as.matrix(aData)
    params <- list(Class="Hocus",assayData = aData)
    lData <- data.frame(outliers="No", prepCells="No", reduceGenes_var="No", reduceGenes_pca="No", stringsAsFactors=FALSE)
    lData <- new("AnnotatedDataFrame", data=lData)
    params$logData <- lData
    
    if (!missing(pheno)) {
        pData <- read.table(pheno, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
        aData <- aData[, row.names(pData)]
        pData_color<-pData
        pData <- new("AnnotatedDataFrame", data = pData)
        params["phenoData"] <- pData
        params$assayData <- aData
    } else {
        pData <- data.frame("GroupID"=1:length(colnames(aData)))
        row.names(pData)<-colnames(aData)
        pData$GroupID <- "Group1"
        pData_color<-pData
        pData <- new("AnnotatedDataFrame", data = pData)
        params$phenoData <- pData
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
        params$featureData <- fData
    }
    
    if (!missing(experiment)) { 
      eData <- read.table(experiment, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
      eData <- new("MIAME", other=as.list(eData))
      params$experimentData <- eData
    }
    
    if (!missing(color) && !missing(pheno)) {
      cData <- read.table(color, header = TRUE, sep = "\t", check.names = FALSE,stringsAsFactors=FALSE)
      if (!all(colnames(cData) %in% colnames(pData))) {
        stop("All column names in the colorData file must match the column names in the phenoData matrix.",call.=FALSE)
      }
      cData <- as.list(cData)
      cData <- lapply(cData, function(x) x[x != ""])
      if (!all(lengths(cData) == apply(data.frame(pData_color[,names(cData)]), 2, function(x) nlevels(as.factor(x))))){
        stop("The number of colors for each group in the colorData file must match the number of levels within the corresponding columns of the phenoData matrix.",call.=FALSE)
      }
      for (i in 1:length(cData)){
        names(cData[[i]]) <- unique(pData_color[,names(cData)[i]])
      }
      params$colorData <- cData
    }
    
    if (missing(pheno)){
      cData <- list(GroupID="grey")
      names(cData[[1]]) <- "Group1"
      params$colorData <- cData
    }
    
    
    #params$assayData <- assayDataNew("environment", exprs=aData)
    params$assayData <- assayDataNew(exprs=aData)
    cellData <- do.call(new,params)
    cellData
    
} 
