#' Generate a heatmap from the matrix of expression values.
#' 
#' Takes ExpressionSet object and creates a highly customizable heatmap, which
#' can display bars that contain information about samples and/or gene clusters.
#' Genes are plotted in rows and samples in columns.
#' 
#' @param cellData ExpressionSet object created with readCells (and preferably 
#'   transformed with prepCells).  It is also helpful to first run 
#'   reduceGenes_var and reduceGenes_pca to keep the number of displayed genes
#'   reasonable.
#' @param clusterCellsBy Character string that can take on the values
#'   'hierarchical' or 'groups' to indicate whether the samples should be
#'   clustered with hiearchical clustering or by specified groupings.
#' @param clusterGenesBy Character string that can take on the values
#'   'hierarchical' or 'groups' to indicate whether the genes should be
#'   clustered with hiearchical clustering or by specified groupings.
#' @param cellGroups Character string indicating the title of the column
#'   containing the sample group designations in pData.  Length should equal the
#'   number of samples in pData.  An example input would be 'KM_Groups', i.e.
#'   the title of one of the columns generated in pData by the clusterCells
#'   function.  Groups are organized from left-to-right according to the
#'   alphabetical order of the factor levels in cellGroups.
#' @param cellOrder Vector of integers indicating the order of the sample groups
#'   to be plotted from left-to-right on the heatmap.  The position in the
#'   vector indicates the order on the heatmap (i.e. first, second, third,
#'   etc.), and the integer specifies the original alphabetical position of the
#'   group.  E.g., for three groups plotted from left-to-right by default
#'   alphabetical order as, 'A', 'B', 'C', the vector c(3,1,2) would plot them
#'   as 'C', 'A', 'B'.
#' @param geneGroups Character string indicating the title of the column
#'   containing the gene group designations in fData.  Length should equal the
#'   total number of genes in fData.  All genes present in the expression table
#'   should have a group value in the fData column.  An example input would be
#'   'KM_Groups', i.e. the title of one of the columns generated in fData by the
#'   clusterGenes function.  Groups are organized from bottom-to-top according
#'   to the alphabetical order of the factor levels in geneGroups
#' @param geneOrder Vector of integers indicating the order of the gene groups
#'   to be plotted from bottom-to-top on the heatmap.  The position in the
#'   vector indicates the order on the heatmap (i.e. first, second, third,
#'   etc.), and the integer specifies the original alphabetical position of the
#'   group.  E.g., for three groups plotted from bottom-to-top by default
#'   alphabetical order as, 'A', 'B', 'C', the vector c(3,1,2) would plot them
#'   as 'C', 'A', 'B'.
#' @param center Numeric specifying a value for the center of the lookup table
#'   for expression values.  Default is to use the center of minimum and maximum
#'   expression values.
#' @param bars Vector of character strings specifying the titles of the columns
#'   in pData that contain the names of groups or values for each sample for
#'   which a color bar should be generated.  If the column contains characters
#'   or factors, each level will be plotted with a different color.  If the
#'   column is numeric, a color gradient spanning two colors will be created.
#' @param colors Vector of character strings specifying the titles of the
#'   columns in pData that contain the colors to be used for plotting bars.  The
#'   order of the title names should match the order of the title names in bars
#'   to which the colors correspond.  If the corresponding bar column is a
#'   character or factor, the number of colors specified should equal the number
#'   of factor levels in the bar column.  The order of colors corresponds to the
#'   order of factor levels in the bar column.  If the corresponding bar column
#'   is numeric, two colors should be specified.
#' @param logNumeric Boolean specifying whether the numeric values in bars
#'   should be transformed to log2 space.
#' @param mapColors Vector of character strings specifying the colors to be used
#'   for the lookup table of expression values.  Low to high values are
#'   specified from left-to-right in the vector.
#' @param save Boolean specifying whether to save the resultant heatmap as a
#'   .tiff file.
#' @return Highly customizeable heatmap plot in a new window.
#' @export

heatMap <- function(cellData, clusterCellsBy = "hierarchical", clusterGenesBy = "hierarchical", cellGroups, cellOrder, geneGroups, 
    geneOrder, center = NA, bars, colors, logNumeric = TRUE, mapColors = c("midnightblue", "dodgerblue3", "mistyrose", "red2", 
        "red4"), save = FALSE) {
    
    if (.Platform$OS.type == "windows") {
        quartz <- function() windows()
    }
    
    samples <- pData(cellData)
    
    exprs_table <- exprs(cellData)
    
    if (!("prepCells" %in% colnames(pData(cellData)))) {
        warning("It would be wise to run prepCells prior to heatMap.", call. = FALSE)
    }
    
    if (!missing(bars)) {
        
        warning("Make certain that bars are in corresponding order to colors", call. = FALSE)
        
        if (!missing(bars) && !all((bars %in% colnames(pData(cellData))))) {
            stop("One or more of the column names specified for bars is not found in phenoData. If absent, add at least one column to phenoData containg values (groups or numeric) for each cell.", 
                call. = FALSE)
        }
        
        if (!missing(colors) && !all((colors %in% colnames(pData(cellData))))) {
            stop("One or more of the column names specified for colors is not found in phenoData. If absent, add at least one column of colors to phenoData that is equal in length to the corresponding number of levels in each groups column.", 
                call. = FALSE)
        }
        
        groupData <- data.frame(samples[, bars])
        colnames(groupData) <- bars
        
        if (missing(colors)) {
            annCol = groupData
            annColors = NULL
        }
        
        if (!missing(colors)) {
            groupColors <- data.frame(samples[, colors])
            
            groupColors_list <- list()
            for (i in 1:length(colors)) {
                groupColors_temp <- data.frame(groupColors[, i])
                groupColors_temp <- as.character(groupColors_temp[!apply(is.na(groupColors_temp) | groupColors_temp == "", 1, 
                  all), ])
                
                if (is.factor(groupData[, i]) | is.character(groupData[, i])) {
                  names(groupColors_temp) <- as.vector(unique(groupData[, i]))
                  groupColors_temp <- groupColors_temp[order(names(groupColors_temp))]
                }
                
                if (is.numeric(groupData[, i]) & logNumeric == TRUE) {
                  groupData[, i][groupData[, i] < 1] <- 1
                  groupData[, i] <- log2(groupData[, i])
                }
                
                groupColors_list[[i]] <- groupColors_temp
            }
            
            names(groupColors) <- bars
            
            annCol <- groupData
            annColors <- groupColors_list
            
        }
    }
    
    if (missing(bars)) {
        
        annCol <- NULL
        annColors <- NULL
        
    }
    
    
    if (clusterCellsBy == "groups") {
        
        if (missing(cellGroups)) {
            stop("No cellGroups are specified!")
        }
        
        if (any(pData(cellData)[, cellGroups] == "" | is.null(pData(cellData)[, cellGroups]))) {
            stop("Not all cells have a group designation within the cellGroups column", call. = FALSE)
        }
        
        
        cell_groups <- data.frame(as.character(pData(cellData)[, cellGroups]))
        o <- order(cell_groups[, 1])
        
        if (!missing(cellOrder)) {
            cell_groups2 <- cell_groups
            comp1 <- as.vector(unique(cell_groups[, 1]))
            comp1 <- comp1[order(comp1)]
            for (i in 1:length(cellOrder)) {
                cell_groups2[cell_groups == comp1[cellOrder[i]]] <- comp1[i]
            }
            o <- order(cell_groups2[, 1])
        }
        
    }
    
    if (clusterGenesBy == "groups") {
        
        if (missing(geneGroups)) {
            stop("No geneGroups are specified!")
        }
        
        if (any(fData(cellData)[row.names(exprs_table), geneGroups] == "" | is.null(fData(cellData)[row.names(exprs_table), geneGroups]))) {
            stop("Not all genes have a group designation within one of the geneGroups columns", call. = FALSE)
        }
        
        genes <- fData(cellData)[row.names(exprs_table), geneGroups]
        gene_groups <- data.frame(factor(genes))
        colnames(gene_groups) <- "Gene Groups"
        row.names(gene_groups) <- row.names(exprs_table)
        og <- order(genes)
        
        if (!missing(geneOrder)) {
            genes2 <- genes
            comp2 <- as.vector(unique(genes))
            comp2 <- comp2[order(comp2)]
            for (i in 1:length(geneOrder)) {
                genes2[genes == comp2[geneOrder[i]]] <- comp2[i]
            }
            og <- order(genes2)
        }
    }
    
    
    shaved_expr_table <- exprs_table
    
    
    if (clusterCellsBy == "hierarchical" && clusterGenesBy == "hierarchical") {
        
        quartz()
        NMF::aheatmap(shaved_expr_table, Rowv = c("euclidean", "ward"), Colv = c("euclidean", "ward"), revC = F, annCol = annCol, 
            annColors = annColors, color = mapColors, breaks = center, width = 8, height = 10, fontsize = 10)
        
        
        if (save == TRUE) {
            a <- NMF::aheatmap(shaved_expr_table, Rowv = c("euclidean", "ward"), Colv = c("euclidean", "ward"), revC = F, annCol = annCol, 
                annColors = annColors, color = mapColors, breaks = center, width = 8, height = 10, fontsize = 10, file = paste("Heatmap.tiff"))
        }
    }
    
    if (clusterCellsBy == "groups" && clusterGenesBy == "hierarchical") {
        
        quartz()
        NMF::aheatmap(shaved_expr_table, Rowv = c("euclidean", "ward"), Colv = o, revC = F, annCol = annCol, annColors = annColors, 
            color = mapColors, breaks = center, width = 8, height = 10, fontsize = 10)
        
        
        if (save == TRUE) {
            a <- NMF::aheatmap(shaved_expr_table, Rowv = c("euclidean", "ward"), Colv = o, revC = F, annCol = annCol, annColors = annColors, 
                color = mapColors, breaks = center, width = 8, height = 10, fontsize = 10, file = paste("Heatmap.tiff"))
        }
    }
    
    
    
    if (clusterCellsBy == "groups" && clusterGenesBy == "groups") {
        
        quartz()
        NMF::aheatmap(shaved_expr_table, Rowv = og, Colv = o, revC = F, annRow = gene_groups, annCol = annCol, annColors = annColors, 
            color = mapColors, breaks = center, width = 8, height = 10, fontsize = 10)
        
        
        if (save == TRUE) {
            a <- NMF::aheatmap(shaved_expr_table, Rowv = og, Colv = o, revC = F, annRow = genes, annCol = annCol, annColors = annColors, 
                color = mapColors, breaks = center, width = 8, height = 10, fontsize = 10, file = paste("Heatmap.tiff"))
        }
    }
    
    if (clusterCellsBy == "hierarchical" && clusterGenesBy == "groups") {
        
        quartz()
        NMF::aheatmap(shaved_expr_table, Rowv = og, Colv = c("euclidean", "ward"), revC = F, annRow = gene_groups, annCol = annCol, 
            annColors = annColors, color = mapColors, breaks = center, width = 8, height = 10, fontsize = 10)
        
        
        if (save == TRUE) {
            a <- NMF::aheatmap(shaved_expr_table, Rowv = og, Colv = c("euclidean", "ward"), revC = F, annRow = gene_groups, annCol = annCol, 
                annColors = annColors, color = mapColors, breaks = center, width = 8, height = 10, fontsize = 10, file = paste("Heatmap.tiff"))
        }
    }
    
} 
