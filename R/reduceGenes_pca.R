#' Reduce number of genes in expression matrix using principal component 
#' analysis (PCA)
#' 
#' Takes ExpressionSet object, performs PCA on the transposed expression matrix,
#' then selects the specified number of most influential genes on the specified 
#' PCs.  Two methods are available for selecting the genes: 1) the genes with 
#' the highest absolute loading value across all the specified PCs are chosen, 
#' or 2) the genes whose loadings show the most significant positive and 
#' negative correlations with the specified PCs are chosen using the dimdesc 
#' function within the FactoMineR package.
#' 
#' @param cellData ExpressionSet object created with readCells (and preferably 
#'   transformed with prepCells).  It is also helpful to first run 
#'   reduceGenes_var.
#' @param corr Boolean when set to TRUE reduces genes using the dimedesc 
#'   function within the FactoMineR package.  With this method, an equal number 
#'   of genes from each of the specified PCs is chosen.  Also, within an 
#'   individual PC, an equal number of genes with positive and negative loadings
#'   is chosen.  E.g. if 200 genes on PCs 1 and 2 are specified, the top 50 
#'   genes with the most significant positive correlation and the top 50 genes 
#'   with the most signficant negative correlation on PC 1 and PC 2 will be 
#'   chosen.  If set to FALSE, the genes with the maximum absolute loading value
#'   across all the specified PCs are chosen.  Be warned that with this method, 
#'   genes from an individual PC can dominate.
#' @param PCs Vector of integers specifying which PCs to select genes from. 
#'   Running pcaMatrix on the full list of genes can be helpful in determining 
#'   which PCs to use.
#' @param genes Integer specifying the desired number of genes to select from 
#'   the pca analysis.  The number of genes in the expression matrix will be 
#'   reduced to this number.
#' @param center Boolean specifying whether to the center the data prior to PCA.
#'   This is generally recommended.
#' @param scale Boolean specifying whether the data should be scaled prior to 
#'   PCA.  This is generally not recommended unless samples have different units
#'   (e.g. some samples are counts and some are TPMs).
#' @param print Boolean specifying whether the results from the PCA analysis 
#'   should be displayed in the terminal window.
#' @param saveTable Boolean specifying whether the results from the PCA analysis
#'   should be saved in a .txt output file.
#' @return ExpressionSet object with genes removed from the expression matrix 
#'   according to the optional parameters specified above.  Note that the 
#'   original list of genes will still be present within fData.  Genes that pass
#'   filter will be stored in fData as TRUE, genes that do not pass filter will 
#'   be stored as FALSE.
#' @export


reduceGenes_pca <- function(cellData, corr = TRUE, PCs = c(1, 2, 3), genes = 300, center = TRUE, scale = FALSE, print = FALSE, saveTable = FALSE) {
    
  if (cellData@logData$prepCells[1] == "No") {
        warning("It would be wise to run prepCells prior to reduceGenes_pca.", call. = FALSE)
    }
    
    # In case reduceGenes_var hasn't been run, remove all genes with undetectable expression across all cells, or else PCA will fail
    
    log.data <- exprs(cellData)
    
    expr.genes <- which(rowSums(log.data > 0) > 0)
  
    log.data <- log.data[expr.genes, ]
    
    exprs(cellData) <- log.data
    
    
    if (corr == FALSE) {
        
        # Run PCA using prcomp and get specified number of top positive and negative eigenvalue genes from specified PCs
        
        PCA.allgenes <- prcomp(t(exprs(cellData)), center = center, scale. = scale)
        
        user_PCs <- data.frame(PCA.allgenes$rotation[, PCs])
        
        user_PCs_abs <- abs(user_PCs)
        max_vals <- apply(user_PCs_abs, 1, max)
        max_vals <- data.frame(max_vals)
        max_vals_index <- which(apply(user_PCs_abs, 1, max) == user_PCs_abs, arr.ind = TRUE)
        max_vals_index <- data.frame(max_vals_index)
        
        unstacked <- unstack(max_vals_index)
        raw_PC_vals <- c()
        
        for (i in 1:length(PCs)) {
            raw_vals <- data.frame(user_PCs[, i])
            raw_vals_subset <- data.frame(raw_vals[unstacked[[i]], ])
            raw_PC_vals <- c(raw_PC_vals, raw_vals_subset)
        }
        
        stacked_raw_PC_vals <- data.frame(stack(raw_PC_vals))
        
        max_vals_rearr <- data.frame(max_vals[row.names(max_vals_index), ])
        row.names(max_vals_rearr) <- row.names(max_vals_index)
        
        max_vals_rearr_index <- data.frame(row.names(max_vals_rearr), max_vals_rearr, stacked_raw_PC_vals[, 1], max_vals_index[, 2])
        
        o <- order(max_vals_rearr_index[, 2], decreasing = TRUE)
        ordered <- max_vals_rearr_index[o, ]
        colnames(ordered) <- c("Gene", "Abs value", "Raw value", "PC Index")
        
        shaved_genes_table <- ordered[1:genes, ]
        shaved_genes_list <- data.frame(ordered[1:genes, 1], stringsAsFactors = FALSE)
        
        # Write table containing top gene names on specified PCs, which genes are found on multiple PCs, the final list of genes
        
        if (print == TRUE) {
            print(shaved_genes_table[, 2:4])
        }
        
        if (saveTable == TRUE) {
            out1 <- write.table(shaved_genes_table, paste("Shaved Genes on PCs", paste(PCs, collapse = ","), "- Value Method", genes, "Genes.txt"), row.names = FALSE, 
                col.names = TRUE, sep = "\t", na = "")
        }
        
        
        # Reduce expression table by final list of genes
        
        shaved_expr_table <- exprs(cellData)[t(shaved_genes_list), ]
        
        PCA_genes_colname <- paste(genes, "PCA genes on PCs", paste(PCs, collapse = ","), "- Values")
        
        fData(cellData)[, PCA_genes_colname] <- row.names(fData(cellData)) %in% shaved_genes_list
        
        exprs(cellData) <- shaved_expr_table
        
    }
    
    # --------------------------------- Corr PCA-------------------------------------------
    
    
    if (corr == TRUE) {
        # Run PCA using FactoMineR, determine correlation of genes with each PC, and get specified number of top correlated genes from specified PCs
        
        PCA.allgenes <- FactoMineR::PCA(t(exprs(cellData)), scale.unit = scale, ncp = 100, graph = F)
        
        # Find genes most highly correlated with each specified PC
        
        dimension.PCA.allgenes <- FactoMineR::dimdesc(PCA.allgenes, axes = PCs, proba = 0.01)
        
        dim_PC_genes <- list()
        dim_PC_genes_nums <- list()
        
        for (i in 1:length(PCs)) {
            dim1 <- as.data.frame(dimension.PCA.allgenes[i])
            dim2 <- row.names(as.data.frame(dimension.PCA.allgenes[i]))
            dim_PC_genes[[i]] <- dim2
            dim1 <- data.frame(dim1, stringsAsFactors = FALSE)
            dim2 <- data.frame(dim2, stringsAsFactors = FALSE)
            row.names(dim1) <- row.names(dim2)
            dim_PC_genes_nums[[i]] <- data.frame(dim2, dim1)
        }
        
        g <- genes/(length(PCs) * 2)
        
        genes.corr.dims <- list()
        genes.nums.corr.dims <- list()
        PCnames <- c()
        
        for (i in 1:length(PCs)) {
            genes.corr.dims[[i]] <- c(head(dim_PC_genes[[i]], n = g), tail(dim_PC_genes[[i]], n = g))
            genes.nums.corr.dims[[i]] <- rbind(head(dim_PC_genes_nums[[i]], n = g), tail(dim_PC_genes_nums[[i]], n = g))
            PCnames <- c(PCnames, paste("Top", as.integer(g) * 2, "genes PC", PCs[i]))
            PCnames <- c(PCnames, paste("Top", as.integer(g) * 2, "corr-values PC", PCs[i]))
            PCnames <- c(PCnames, paste("Top", as.integer(g) * 2, "p-values PC", PCs[i]))
        }
        
        genes.corr.dims <- do.call(cbind, genes.corr.dims)
        genes.nums.corr.dims <- do.call(cbind, genes.nums.corr.dims)
        
        redundant.genes.corr.dims <- matrix(genes.corr.dims, ncol = 1)
        redundant.genes.corr.dims <- data.frame(table(redundant.genes.corr.dims))
        redundant.genes.corr.dims <- redundant.genes.corr.dims[redundant.genes.corr.dims$Freq > 1, ]
        
        unique.genes.corr.dims <- matrix(genes.corr.dims, ncol = 1)
        unique.genes.corr.dims <- unique(unique.genes.corr.dims)
        
        genes.nums.corr.dims <- data.frame(genes.nums.corr.dims)
        redundant.genes.corr.dims <- data.frame(redundant.genes.corr.dims[, 1])
        unique.genes.corr.dims <- data.frame(unique.genes.corr.dims)
        
        if (length(t(unique.genes.corr.dims)) > length(genes.nums.corr.dims[, 1])) {
            genes.nums.corr.dims[(length(genes.nums.corr.dims[, 1]) + 1):length(t(unique.genes.corr.dims)), ] <- NA
            redundant.genes.corr.dims[(length(t(redundant.genes.corr.dims)) + 1):length(t(unique.genes.corr.dims)), ] <- NA
        }
        
        PC.genes.summary <- cbind(genes.nums.corr.dims, redundant.genes.corr.dims, unique.genes.corr.dims)
        PC.genes.print <- genes.nums.corr.dims
        colnames(PC.genes.print) <- PCnames
        PCnames <- c(PCnames, "Duplicate Genes", "Final List")
        colnames(PC.genes.summary) <- PCnames
        
        # Write table containing top gene names on specified PCs, which genes are found on multiple PCs, the final list of genes
        
        if (print == TRUE) {
            print(PC.genes.print[1:(genes/length(PCs)), ])
        }
        
        if (saveTable == TRUE) {
            out1 <- write.table(PC.genes.summary, paste("Shaved Genes on PCs", paste(PCs, collapse = ","), "- Corr Method", genes, "Genes.txt"), row.names = FALSE, 
                col.names = TRUE, sep = "\t", na = "")
        }
        
        # Reduce expression table by final list of genes
        
        shaved_expr_table <- exprs(cellData)[t(unique.genes.corr.dims), ]
        
        PCA_genes_colname <- paste(genes, "PCA genes on PCs", paste(PCs, collapse = ","), "- Corr")
        
        fData(cellData)[, PCA_genes_colname] <- row.names(fData(cellData)) %in% t(unique.genes.corr.dims)
        
        exprs(cellData) <- shaved_expr_table
        
    }
    
    cellData@logData$reduceGenes_pca[1] <- "Yes"
    
    cellData
    
} 
