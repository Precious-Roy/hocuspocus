#' Reduce number of genes in expression matrix based on variance
#' 
#' Takes ExpressionSet object and (optionally) 1) eliminates genes that are 
#' expressed below the LOD in a specified number of samples, 2) eliminates genes
#' whose coefficient of variation (CV) across cells is below a set limit, and 3)
#' eliminates genes using the CV z-score/binning method within the Seurat 
#' package.
#' 
#' @param cellData ExpressionSet object created with readCells (and preferably 
#'   transformed with prepCells)
#' @param exprGenes Boolean specifying whether to remove genes that are 
#'   expressed above the expression threshold (exprThresh) in fewer than the 
#'   number of samples specified with cellThresh.
#' @param exprThresh Numeric specifying the the threshold for considering a gene
#'   to be expressed.  Note that if prepCells has been run, expression values 
#'   will be on a log2 scale, so exprThresh = 1 corresponds to an expression 
#'   level of 2 in untransformed space.
#' @param cellThresh Integer indicating the number of cells in which a gene has 
#'   to be expressed above the expression threshold (exprThresh) in order to be 
#'   kept.
#' @param varThresh Boolean specifying whether to remove genes whose CV across 
#'   samples is less than the threshold specified with cv.  If exprGenes and 
#'   varThresh are both TRUE, exprGenes is performed prior to varThresh.
#' @param cv Numeric specifying the CV threshold for varThresh.  Values ranging 
#'   from 0.1 to 0.5 are typically effective.
#' @param seuratThresh Boolean specifying whether the Seurat package's method of
#'   elminating genes by variance is to be used.  See the Seurat package 
#'   documentation for more details.  If exprThresh, varThresh, and seuratThresh
#'   are all TRUE, they are performed in that order.
#' @param z_cutoff Numeric specifying the z-score cutoff for the Seurat method.
#' @param ... Pass more options to the Seurat function mean.var.plot()
#' @return ExpressionSet object with genes removed from the expression matrix 
#'   according to the optional parameters specified above.  Note that the 
#'   original list of genes will still be present within fData.  Genes that pass
#'   filter will be stored in fData as TRUE, genes that do not pass filter will 
#'   be stored as FALSE.

reduceGenes_var <- function(cellData, exprGenes = TRUE, exprThresh = 1, cellThresh = 2, varThresh = TRUE, cv = 0.5, seuratThresh = FALSE, z_cutoff = 1.5, ...) {
    
    if (!("prepCells" %in% colnames(pData(cellData)))) {
        warning("It would be wise to run prepCells prior to reduceGenes_var", call. = FALSE)
    }
    
    log.data <- exprs(cellData)
    
    if (exprGenes == TRUE) {
        
        # Find genes expressed (>LOD, i.e. a two-fold change) in three or more cells and limit the gene list down to those genes
        
        expr.genes <- which(rowSums(log.data > exprThresh) > (cellThresh))
        
        log.data <- log.data[expr.genes, ]
        
        fData(cellData)$ExpressedGenes <- row.names(fData(cellData)) %in% row.names(log.data)
        
        exprs(cellData) <- log.data
        
    }
    
    if (varThresh == TRUE) {
        
        # Find genes with a c.v. between cells greater than varThresh and limit the gene list down to those genes
        
        co.var <- function(x) {
            sd(x)/mean(x)
        }
        
        log.data <- log.data[which(apply(log.data, 1, co.var) >= cv), ]
        
        fData(cellData)$CV_Genes <- row.names(fData(cellData)) %in% row.names(log.data)
        
        exprs(cellData) <- log.data
        
    }
    
    if (seuratThresh == TRUE) {
        
        nbt <- new("seurat", raw.data = data.frame(log.data))
        
        nbt <- setup(nbt, project = "NBT", min.cells = 3, min.genes = 1, is.expr = 1e-04)
        
        quartz()
        
        nbt <- mean.var.plot(nbt, y.cutoff = z_cutoff, x.low.cutoff = 1, x.high.cutoff = 14, fxn.x = expMean, fxn.y = logVarDivMean, ...)
        
        log.data <- log.data[nbt@var.genes, ]
        
        fData(cellData)$SeuratGenes <- row.names(fData(cellData)) %in% nbt@var.genes
        
        exprs(cellData) <- log.data
        
    }
    
    cellData
    
} 
