#' Generate a series of plots comparing the first 10 PCs of a PCA analysis.
#' 
#' Takes ExpressionSet object, performs PCA on the transposed expression matrix,
#' then plots the projections of the sample scores on the first 10 PCs.
#' 
#' @param cellData ExpressionSet object created with readCells (and preferably 
#'   transformed with prepCells).  It is also helpful to first run 
#'   reduceGenes_var.
#' @param scree Boolean specifying whether to generate a scree plot showing the 
#'   amount of variance contributed by each PC.
#' @param center Boolean specifying whether to the center the data prior to PCA.
#'   This is generally recommended.
#' @param scale Boolean specifying whether the data should be scaled prior to 
#'   PCA.  This is generally not recommended unless samples have different units
#'   (e.g. some samples are counts and some are TPMs).
#' @param groups Character string specifying the title of the column in pData 
#'   that contains the names of the groups to which each sample belongs.  The 
#'   dots representing each sample in the plots will be colored by group. The 
#'   column length should be the same length as the number of samples.
#' @param values Character string specifying the title of the column in pData 
#'   that contains a vector of numeric values to be plotted as a color gradient 
#'   on the dots in the plots.  The column length should be the same length as 
#'   the number of samples.  Gene and values cannot be specified simultaneously.
#'   If groups and values are both specified, gene will be used for coloring. 
#'   Set bubble to TRUE to display group information as well.
#' @param gene Character string specifying a gene whose expression values will 
#'   be plotted a color gradient on the dots in the plots.  One gene name should
#'   be specified, and the gene must be present within the expression table. 
#'   Gene and values cannot be specified simultaneously.  If groups and gene are
#'   both specified, gene will be used for coloring.  Set bubble to TRUE to 
#'   display group information as well.
#' @param colors Character string specifying the title of the column in pData
#'   that contains the colors to be used for plotting.  If groups is specified,
#'   length of colors should be equal to the number of levels in groups.  The
#'   order of colors corresponds to the order of factor levels in groups.  If
#'   values or gene is specified, length of colors should be 2, or else a
#'   default black-to-yellow color gradient will be used.
#' @param logNumeric Boolean specifying whether the numbers in values should be 
#'   transformed to log2 space.
#' @param refPC Character string specifying the PC to be used as the reference 
#'   PC that is plotted on the x-axis of every plot.
#' @param alpha Numeric specifying the transparency (from 0 to 1) level of the 
#'   dot colors.
#' @param dotsize Numeric specifying the size of the dots in the plots.
#' @param bubble Boolean specifying whether dots of different sizes should be
#'   plotted to display information about another variable on the plot.
#' @param bubbleSizes Vector of numerics specifying the dot sizes for each group
#'   within groups.  Length of bubbleSizes should be the same length as groups. 
#'   If bubbleSizes is not specified, default sizes will be generated according 
#'   to the levels in the groups column.  Useful for displaying group 
#'   information along with values or gene.  Additionally, can be used to 
#'   indicate additional information about groups.  E.g. if the levels in groups
#'   are 'A1', 'B1', 'C2', and 'D2', bubbleSizes could be c(3,3,6,6) to indicate
#'   the '1' and '2' components of the groups with different dot sizes.
#' @param print Boolean specifying whether the genes with the most positive and 
#'   negative loadings on each of the 10 PCs should be printed in the terminal 
#'   window.
#' @param printNum Integer specifying the number of genes from each PC to print.
#' @param save Boolean specifying whether to save the resultant plots as a .tiff
#'   file.
#' @return Matrix of PCA plots for the first 10 PCs.
#' @import ggplot2
#' @import gridExtra

pcaMatrix <- function(cellData, scree = FALSE, center = TRUE, scale = FALSE, groups, values, gene, colors, logNumeric = TRUE, refPC = "PC1", alpha = 0.7, dotsize = 2, 
    bubble = FALSE, bubbleSizes, print = FALSE, printNum = 50, save = FALSE) {
    
    if (.Platform$OS.type == "windows") {
        quartz <- function() windows()
    }
    
    if (!("prepCells" %in% colnames(pData(cellData)))) {
        warning("It would be wise to run prepCells prior to pcaMatrix", call. = FALSE)
    }
    
    if (!missing(groups) && !all((groups %in% colnames(pData(cellData))))) {
        stop("The column name specified for groups is not found in phenoData. If absent, add a column to phenoData containing the group for each cell, even if it only consists of 1 group.", 
            call. = FALSE)
    }
    
    if (!missing(values) && !all((values %in% colnames(pData(cellData))))) {
        stop("The column name specified for values is not found in phenoData. If absent, add a column of colors to phenoData that is equal in length to the number of levels in groups.", 
            call. = FALSE)
    }
    
    if (!missing(colors) && !all((colors %in% colnames(pData(cellData))))) {
        stop("The column name specified for colors is not found in phenoData. If absent, add a column of colors to phenoData that is equal in length to the number of levels in groups.", 
            call. = FALSE)
    }
    
    PCA.allgenes <- prcomp(t(exprs(cellData)), center = center, scale. = scale, save = FALSE)
    
    if (missing(colors)) {
        cell_colors <- "blue"
    }
    
    if (!missing(colors)) {
        cell_colors <- data.frame(pData(cellData)[, colors], stringsAsFactors = FALSE)
        cell_colors <- as.character(cell_colors[!apply(is.na(cell_colors) | cell_colors == "", 1, all), ])
    }
    
    if (scree == TRUE) {
        
        quartz()
        screeplot(PCA.allgenes, type = "lines", main = "Scree Plot")
        
        if (save == TRUE) {
            quartz(type = "pdf", file = "Scree Plot.pdf")
            screeplot(PCA.allgenes, type = "lines")
            dev.off
        }
    }
    
    
    comp2 <- data.frame(PCA.allgenes$x[, 1:10])
    comp.cell.type2 <- comp2
    
    if (missing(groups)) {
        comp.cell.type2["Cell_Type"] <- rep("All", nrow(pData(cellData)))
    }
    
    if (!missing(groups)) {
        groupz <- data.frame(pData(cellData)[, groups], stringsAsFactors = FALSE)
        comp.cell.type2["Cell_Type"] <- factor(groupz[, colnames(groupz)], levels = unique(groupz[, colnames(groupz)]), ordered = FALSE)
    }
    
    if (bubble == TRUE) {
        
        if (missing(bubbleSizes)) {
            sizes <- 1:length(levels(comp.cell.type2$Cell_Type))
        }
        
        if (!missing(bubbleSizes)) {
            sizes <- bubbleSizes
        }
        
    }
    
    if (bubble == FALSE) {
        sizes <- vector(length = length(cell_colors))
        sizes[1:length(cell_colors)] <- dotsize
    }
    
    if (missing(values) && missing(gene)) {
        
        if (bubble == TRUE) {
            warning("Bubble sizes currently reflect the number of levels in groups and are not adding additional information to the plot.", call. = FALSE)
        }
        
        
        g1 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC1", fill = "factor(Cell_Type,levels=unique(Cell_Type),ordered=FALSE)", size = "factor(Cell_Type,levels=unique(Cell_Type),ordered=FALSE)"), 
            shape = 21, color = "black", alpha = alpha) + scale_fill_manual(name = "Cell Type", values = cell_colors) + scale_size_manual(name = "Cell Type", values = sizes) + 
            guides(fill = guide_legend("Cell Type"), color = guide_legend("Cell Type"), size = guide_legend("Cell Type")) + theme(legend.key.size = grid::unit(0.35, 
            "cm"), legend.text = element_text(size = 8))
        g2 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC2", fill = "factor(Cell_Type,levels=unique(Cell_Type),ordered=FALSE)", size = "factor(Cell_Type,levels=unique(Cell_Type),ordered=FALSE)"), 
            shape = 21, color = "black", alpha = alpha) + scale_fill_manual(name = "Cell Type", values = cell_colors) + scale_size_manual(name = "Cell Type", values = sizes) + 
            guides(fill = guide_legend("Cell Type"), color = guide_legend("Cell Type"), size = guide_legend("Cell Type")) + theme(legend.key.size = grid::unit(0.35, 
            "cm"), legend.text = element_text(size = 8))
        g3 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC3", fill = "factor(Cell_Type,levels=unique(Cell_Type),ordered=FALSE)", size = "factor(Cell_Type,levels=unique(Cell_Type),ordered=FALSE)"), 
            shape = 21, color = "black", alpha = alpha) + scale_fill_manual(name = "Cell Type", values = cell_colors) + scale_size_manual(name = "Cell Type", values = sizes) + 
            guides(fill = guide_legend("Cell Type"), color = guide_legend("Cell Type"), size = guide_legend("Cell Type")) + theme(legend.key.size = grid::unit(0.35, 
            "cm"), legend.text = element_text(size = 8))
        g4 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC4", fill = "factor(Cell_Type,levels=unique(Cell_Type),ordered=FALSE)", size = "factor(Cell_Type,levels=unique(Cell_Type),ordered=FALSE)"), 
            shape = 21, color = "black", alpha = alpha) + scale_fill_manual(name = "Cell Type", values = cell_colors) + scale_size_manual(name = "Cell Type", values = sizes) + 
            guides(fill = guide_legend("Cell Type"), color = guide_legend("Cell Type"), size = guide_legend("Cell Type")) + theme(legend.key.size = grid::unit(0.35, 
            "cm"), legend.text = element_text(size = 8))
        g5 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC5", fill = "factor(Cell_Type,levels=unique(Cell_Type),ordered=FALSE)", size = "factor(Cell_Type,levels=unique(Cell_Type),ordered=FALSE)"), 
            shape = 21, color = "black", alpha = alpha) + scale_fill_manual(name = "Cell Type", values = cell_colors) + scale_size_manual(name = "Cell Type", values = sizes) + 
            guides(fill = guide_legend("Cell Type"), color = guide_legend("Cell Type"), size = guide_legend("Cell Type")) + theme(legend.key.size = grid::unit(0.35, 
            "cm"), legend.text = element_text(size = 8))
        g6 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC6", fill = "factor(Cell_Type,levels=unique(Cell_Type),ordered=FALSE)", size = "factor(Cell_Type,levels=unique(Cell_Type),ordered=FALSE)"), 
            shape = 21, color = "black", alpha = alpha) + scale_fill_manual(name = "Cell Type", values = cell_colors) + scale_size_manual(name = "Cell Type", values = sizes) + 
            guides(fill = guide_legend("Cell Type"), color = guide_legend("Cell Type"), size = guide_legend("Cell Type")) + theme(legend.key.size = grid::unit(0.35, 
            "cm"), legend.text = element_text(size = 8))
        g7 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC7", fill = "factor(Cell_Type,levels=unique(Cell_Type),ordered=FALSE)", size = "factor(Cell_Type,levels=unique(Cell_Type),ordered=FALSE)"), 
            shape = 21, color = "black", alpha = alpha) + scale_fill_manual(name = "Cell Type", values = cell_colors) + scale_size_manual(name = "Cell Type", values = sizes) + 
            guides(fill = guide_legend("Cell Type"), color = guide_legend("Cell Type"), size = guide_legend("Cell Type")) + theme(legend.key.size = grid::unit(0.35, 
            "cm"), legend.text = element_text(size = 8))
        g8 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC8", fill = "factor(Cell_Type,levels=unique(Cell_Type),ordered=FALSE)", size = "factor(Cell_Type,levels=unique(Cell_Type),ordered=FALSE)"), 
            shape = 21, color = "black", alpha = alpha) + scale_fill_manual(name = "Cell Type", values = cell_colors) + scale_size_manual(name = "Cell Type", values = sizes) + 
            guides(fill = guide_legend("Cell Type"), color = guide_legend("Cell Type"), size = guide_legend("Cell Type")) + theme(legend.key.size = grid::unit(0.35, 
            "cm"), legend.text = element_text(size = 8))
        g9 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC9", fill = "factor(Cell_Type,levels=unique(Cell_Type),ordered=FALSE)", size = "factor(Cell_Type,levels=unique(Cell_Type),ordered=FALSE)"), 
            shape = 21, color = "black", alpha = alpha) + scale_fill_manual(name = "Cell Type", values = cell_colors) + scale_size_manual(name = "Cell Type", values = sizes) + 
            guides(fill = guide_legend("Cell Type"), color = guide_legend("Cell Type"), size = guide_legend("Cell Type")) + theme(legend.key.size = grid::unit(0.35, 
            "cm"), legend.text = element_text(size = 8))
        
        if (save == TRUE) {
            PCA_grob <- arrangeGrob(g1, g2, g3, g4, g5, g6, g7, g8, g9, ncol = 3)
            ggsave(PCA_grob, file = paste("PCA Matrix.tiff"))
        }
        
        quartz()
        PCA_grid <- grid.arrange(g1, g2, g3, g4, g5, g6, g7, g8, g9, ncol = 3)
        
    }
    
    if (!missing(values)) {
        
        if (!missing(gene)) {
            stop("Cannot specify both values and gene.", call. = FALSE)
        }
        
        if (length(cell_colors) != 2) {
            warning("Specify a vector of two colors to control color, using default colors instead", call. = FALSE)
            cell_colors[1:2] <- c("black", "yellow")
        }
        
        valuez <- data.frame(pData(cellData)[, values], stringsAsFactors = FALSE)
        if (logNumeric == TRUE) {
            valuez <- log2(valuez)
        }
        comp.cell.type2["Values"] <- valuez
        
        if (bubble == TRUE) {
            g1 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC1", fill = "Values", size = "factor(Cell_Type,levels=unique(Cell_Type),ordered=FALSE)"), 
                shape = 21, color = "black", alpha = alpha) + scale_size_manual(name = "Groups", values = sizes) + guides(color = guide_legend("Values"), size = guide_legend("Groups")) + 
                scale_fill_gradient(low = cell_colors[1], high = cell_colors[2]) + theme(legend.key.size = grid::unit(0.35, "cm"), legend.text = element_text(size = 8))
            g2 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC2", fill = "Values", size = "factor(Cell_Type,levels=unique(Cell_Type),ordered=FALSE)"), 
                shape = 21, color = "black", alpha = alpha) + scale_size_manual(name = "Groups", values = sizes) + guides(color = guide_legend("Values"), size = guide_legend("Groups")) + 
                scale_fill_gradient(low = cell_colors[1], high = cell_colors[2]) + theme(legend.key.size = grid::unit(0.35, "cm"), legend.text = element_text(size = 8))
            g3 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC3", fill = "Values", size = "factor(Cell_Type,levels=unique(Cell_Type),ordered=FALSE)"), 
                shape = 21, color = "black", alpha = alpha) + scale_size_manual(name = "Groups", values = sizes) + guides(color = guide_legend("Values"), size = guide_legend("Groups")) + 
                scale_fill_gradient(low = cell_colors[1], high = cell_colors[2]) + theme(legend.key.size = grid::unit(0.35, "cm"), legend.text = element_text(size = 8))
            g4 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC4", fill = "Values", size = "factor(Cell_Type,levels=unique(Cell_Type),ordered=FALSE)"), 
                shape = 21, color = "black", alpha = alpha) + scale_size_manual(name = "Groups", values = sizes) + guides(color = guide_legend("Values"), size = guide_legend("Groups")) + 
                scale_fill_gradient(low = cell_colors[1], high = cell_colors[2]) + theme(legend.key.size = grid::unit(0.35, "cm"), legend.text = element_text(size = 8))
            g5 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC5", fill = "Values", size = "factor(Cell_Type,levels=unique(Cell_Type),ordered=FALSE)"), 
                shape = 21, color = "black", alpha = alpha) + scale_size_manual(name = "Groups", values = sizes) + guides(color = guide_legend("Values"), size = guide_legend("Groups")) + 
                scale_fill_gradient(low = cell_colors[1], high = cell_colors[2]) + theme(legend.key.size = grid::unit(0.35, "cm"), legend.text = element_text(size = 8))
            g6 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC6", fill = "Values", size = "factor(Cell_Type,levels=unique(Cell_Type),ordered=FALSE)"), 
                shape = 21, color = "black", alpha = alpha) + scale_size_manual(name = "Groups", values = sizes) + guides(color = guide_legend("Values"), size = guide_legend("Groups")) + 
                scale_fill_gradient(low = cell_colors[1], high = cell_colors[2]) + theme(legend.key.size = grid::unit(0.35, "cm"), legend.text = element_text(size = 8))
            g7 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC7", fill = "Values", size = "factor(Cell_Type,levels=unique(Cell_Type),ordered=FALSE)"), 
                shape = 21, color = "black", alpha = alpha) + scale_size_manual(name = "Groups", values = sizes) + guides(color = guide_legend("Values"), size = guide_legend("Groups")) + 
                scale_fill_gradient(low = cell_colors[1], high = cell_colors[2]) + theme(legend.key.size = grid::unit(0.35, "cm"), legend.text = element_text(size = 8))
            g8 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC8", fill = "Values", size = "factor(Cell_Type,levels=unique(Cell_Type),ordered=FALSE)"), 
                shape = 21, color = "black", alpha = alpha) + scale_size_manual(name = "Groups", values = sizes) + guides(color = guide_legend("Values"), size = guide_legend("Groups")) + 
                scale_fill_gradient(low = cell_colors[1], high = cell_colors[2]) + theme(legend.key.size = grid::unit(0.35, "cm"), legend.text = element_text(size = 8))
            g9 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC9", fill = "Values", size = "factor(Cell_Type,levels=unique(Cell_Type),ordered=FALSE)"), 
                shape = 21, color = "black", alpha = alpha) + scale_size_manual(name = "Groups", values = sizes) + guides(color = guide_legend("Values"), size = guide_legend("Groups")) + 
                scale_fill_gradient(low = cell_colors[1], high = cell_colors[2]) + theme(legend.key.size = grid::unit(0.35, "cm"), legend.text = element_text(size = 8))
            
            if (save == TRUE) {
                PCA_grob <- arrangeGrob(g1, g2, g3, g4, g5, g6, g7, g8, g9, ncol = 3)
                ggsave(PCA_grob, file = paste("PCA Matrix.tiff"))
            }
            
            quartz()
            PCA_grid <- grid.arrange(g1, g2, g3, g4, g5, g6, g7, g8, g9, ncol = 3)
            
        }
        
        if (bubble == FALSE) {
            g1 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC1", fill = "Values"), size = dotsize, shape = 21, color = "black", alpha = alpha) + 
                guides(color = guide_legend(title = "Values")) + scale_fill_gradient(low = cell_colors[1], high = cell_colors[2]) + theme(legend.key.size = grid::unit(0.35, 
                "cm"), legend.text = element_text(size = 8))
            g2 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC2", fill = "Values"), size = dotsize, shape = 21, color = "black", alpha = alpha) + 
                guides(color = guide_legend(title = "Values")) + scale_fill_gradient(low = cell_colors[1], high = cell_colors[2]) + theme(legend.key.size = grid::unit(0.35, 
                "cm"), legend.text = element_text(size = 8))
            g3 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC3", fill = "Values"), size = dotsize, shape = 21, color = "black", alpha = alpha) + 
                guides(color = guide_legend(title = "Values")) + scale_fill_gradient(low = cell_colors[1], high = cell_colors[2]) + theme(legend.key.size = grid::unit(0.35, 
                "cm"), legend.text = element_text(size = 8))
            g4 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC4", fill = "Values"), size = dotsize, shape = 21, color = "black", alpha = alpha) + 
                guides(color = guide_legend(title = "Values")) + scale_fill_gradient(low = cell_colors[1], high = cell_colors[2]) + theme(legend.key.size = grid::unit(0.35, 
                "cm"), legend.text = element_text(size = 8))
            g5 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC5", fill = "Values"), size = dotsize, shape = 21, color = "black", alpha = alpha) + 
                guides(color = guide_legend(title = "Values")) + scale_fill_gradient(low = cell_colors[1], high = cell_colors[2]) + theme(legend.key.size = grid::unit(0.35, 
                "cm"), legend.text = element_text(size = 8))
            g6 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC6", fill = "Values"), size = dotsize, shape = 21, color = "black", alpha = alpha) + 
                guides(color = guide_legend(title = "Values")) + scale_fill_gradient(low = cell_colors[1], high = cell_colors[2]) + theme(legend.key.size = grid::unit(0.35, 
                "cm"), legend.text = element_text(size = 8))
            g7 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC7", fill = "Values"), size = dotsize, shape = 21, color = "black", alpha = alpha) + 
                guides(color = guide_legend(title = "Values")) + scale_fill_gradient(low = cell_colors[1], high = cell_colors[2]) + theme(legend.key.size = grid::unit(0.35, 
                "cm"), legend.text = element_text(size = 8))
            g8 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC8", fill = "Values"), size = dotsize, shape = 21, color = "black", alpha = alpha) + 
                guides(color = guide_legend(title = "Values")) + scale_fill_gradient(low = cell_colors[1], high = cell_colors[2]) + theme(legend.key.size = grid::unit(0.35, 
                "cm"), legend.text = element_text(size = 8))
            g9 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC9", fill = "Values"), size = dotsize, shape = 21, color = "black", alpha = alpha) + 
                guides(color = guide_legend(title = "Values")) + scale_fill_gradient(low = cell_colors[1], high = cell_colors[2]) + theme(legend.key.size = grid::unit(0.35, 
                "cm"), legend.text = element_text(size = 8))
            
            if (save == TRUE) {
                PCA_grob <- arrangeGrob(g1, g2, g3, g4, g5, g6, g7, g8, g9, ncol = 3)
                ggsave(PCA_grob, file = paste("PCA Matrix.tiff"))
            }
            
            quartz()
            PCA_grid <- grid.arrange(g1, g2, g3, g4, g5, g6, g7, g8, g9, ncol = 3)
            
        }
        
    }
    
    if (!missing(gene)) {
        
        if (!missing(values)) {
            stop("Cannot specify both values and gene.", call. = FALSE)
        }
        
        if (!gene %in% row.names(exprs(cellData))) {
            stop("Specified gene is not present in expression table", call. = FALSE)
        }
        
        if (length(cell_colors) != 2) {
            warning("Specify a vector of two colors to control color, using default colors instead", call. = FALSE)
            cell_colors[1:2] <- c("black", "yellow")
        }
        
        genez <- data.frame(exprs(cellData)[gene, ], stringsAsFactors = FALSE)[, 1]
        comp.cell.type2[gene] <- genez
        
        if (bubble == TRUE) {
            
            g1 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC1", fill = gene, size = "factor(Cell_Type,levels=unique(Cell_Type),ordered=FALSE)"), 
                shape = 21, color = "black", alpha = alpha) + scale_size_manual(name = "Groups", values = sizes) + guides(color = guide_legend(title = gene), size = guide_legend("Groups")) + 
                scale_fill_gradient(low = cell_colors[1], high = cell_colors[2]) + theme(legend.key.size = grid::unit(0.35, "cm"), legend.text = element_text(size = 8))
            g2 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC2", fill = gene, size = "factor(Cell_Type,levels=unique(Cell_Type),ordered=FALSE)"), 
                shape = 21, color = "black", alpha = alpha) + scale_size_manual(name = "Groups", values = sizes) + guides(color = guide_legend(title = gene), size = guide_legend("Groups")) + 
                scale_fill_gradient(low = cell_colors[1], high = cell_colors[2]) + theme(legend.key.size = grid::unit(0.35, "cm"), legend.text = element_text(size = 8))
            g3 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC3", fill = gene, size = "factor(Cell_Type,levels=unique(Cell_Type),ordered=FALSE)"), 
                shape = 21, color = "black", alpha = alpha) + scale_size_manual(name = "Groups", values = sizes) + guides(color = guide_legend(title = gene), size = guide_legend("Groups")) + 
                scale_fill_gradient(low = cell_colors[1], high = cell_colors[2]) + theme(legend.key.size = grid::unit(0.35, "cm"), legend.text = element_text(size = 8))
            g4 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC4", fill = gene, size = "factor(Cell_Type,levels=unique(Cell_Type),ordered=FALSE)"), 
                shape = 21, color = "black", alpha = alpha) + scale_size_manual(name = "Groups", values = sizes) + guides(color = guide_legend(title = gene), size = guide_legend("Groups")) + 
                scale_fill_gradient(low = cell_colors[1], high = cell_colors[2]) + theme(legend.key.size = grid::unit(0.35, "cm"), legend.text = element_text(size = 8))
            g5 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC5", fill = gene, size = "factor(Cell_Type,levels=unique(Cell_Type),ordered=FALSE)"), 
                shape = 21, color = "black", alpha = alpha) + scale_size_manual(name = "Groups", values = sizes) + guides(color = guide_legend(title = gene), size = guide_legend("Groups")) + 
                scale_fill_gradient(low = cell_colors[1], high = cell_colors[2]) + theme(legend.key.size = grid::unit(0.35, "cm"), legend.text = element_text(size = 8))
            g6 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC6", fill = gene, size = "factor(Cell_Type,levels=unique(Cell_Type),ordered=FALSE)"), 
                shape = 21, color = "black", alpha = alpha) + scale_size_manual(name = "Groups", values = sizes) + guides(color = guide_legend(title = gene), size = guide_legend("Groups")) + 
                scale_fill_gradient(low = cell_colors[1], high = cell_colors[2]) + theme(legend.key.size = grid::unit(0.35, "cm"), legend.text = element_text(size = 8))
            g7 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC7", fill = gene, size = "factor(Cell_Type,levels=unique(Cell_Type),ordered=FALSE)"), 
                shape = 21, color = "black", alpha = alpha) + scale_size_manual(name = "Groups", values = sizes) + guides(color = guide_legend(title = gene), size = guide_legend("Groups")) + 
                scale_fill_gradient(low = cell_colors[1], high = cell_colors[2]) + theme(legend.key.size = grid::unit(0.35, "cm"), legend.text = element_text(size = 8))
            g8 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC8", fill = gene, size = "factor(Cell_Type,levels=unique(Cell_Type),ordered=FALSE)"), 
                shape = 21, color = "black", alpha = alpha) + scale_size_manual(name = "Groups", values = sizes) + guides(color = guide_legend(title = gene), size = guide_legend("Groups")) + 
                scale_fill_gradient(low = cell_colors[1], high = cell_colors[2]) + theme(legend.key.size = grid::unit(0.35, "cm"), legend.text = element_text(size = 8))
            g9 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC9", fill = gene, size = "factor(Cell_Type,levels=unique(Cell_Type),ordered=FALSE)"), 
                shape = 21, color = "black", alpha = alpha) + scale_size_manual(name = "Groups", values = sizes) + guides(color = guide_legend(title = gene), size = guide_legend("Groups")) + 
                scale_fill_gradient(low = cell_colors[1], high = cell_colors[2]) + theme(legend.key.size = grid::unit(0.35, "cm"), legend.text = element_text(size = 8))
            
            if (save == TRUE) {
                PCA_grob <- arrangeGrob(g1, g2, g3, g4, g5, g6, g7, g8, g9, ncol = 3)
                ggsave(PCA_grob, file = paste("PCA Matrix.tiff"))
            }
            
            quartz()
            PCA_grid <- grid.arrange(g1, g2, g3, g4, g5, g6, g7, g8, g9, ncol = 3)
            
        }
        
        if (bubble == FALSE) {
            g1 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC1", fill = gene), size = dotsize, shape = 21, color = "black", alpha = alpha) + 
                guides(color = guide_legend(title = gene)) + scale_fill_gradient(low = cell_colors[1], high = cell_colors[2]) + theme(legend.key.size = grid::unit(0.35, 
                "cm"), legend.text = element_text(size = 8))
            g2 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC2", fill = gene), size = dotsize, shape = 21, color = "black", alpha = alpha) + 
                guides(color = guide_legend(title = gene)) + scale_fill_gradient(low = cell_colors[1], high = cell_colors[2]) + theme(legend.key.size = grid::unit(0.35, 
                "cm"), legend.text = element_text(size = 8))
            g3 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC3", fill = gene), size = dotsize, shape = 21, color = "black", alpha = alpha) + 
                guides(color = guide_legend(title = gene)) + scale_fill_gradient(low = cell_colors[1], high = cell_colors[2]) + theme(legend.key.size = grid::unit(0.35, 
                "cm"), legend.text = element_text(size = 8))
            g4 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC4", fill = gene), size = dotsize, shape = 21, color = "black", alpha = alpha) + 
                guides(color = guide_legend(title = gene)) + scale_fill_gradient(low = cell_colors[1], high = cell_colors[2]) + theme(legend.key.size = grid::unit(0.35, 
                "cm"), legend.text = element_text(size = 8))
            g5 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC5", fill = gene), size = dotsize, shape = 21, color = "black", alpha = alpha) + 
                guides(color = guide_legend(title = gene)) + scale_fill_gradient(low = cell_colors[1], high = cell_colors[2]) + theme(legend.key.size = grid::unit(0.35, 
                "cm"), legend.text = element_text(size = 8))
            g6 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC6", fill = gene), size = dotsize, shape = 21, color = "black", alpha = alpha) + 
                guides(color = guide_legend(title = gene)) + scale_fill_gradient(low = cell_colors[1], high = cell_colors[2]) + theme(legend.key.size = grid::unit(0.35, 
                "cm"), legend.text = element_text(size = 8))
            g7 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC7", fill = gene), size = dotsize, shape = 21, color = "black", alpha = alpha) + 
                guides(color = guide_legend(title = gene)) + scale_fill_gradient(low = cell_colors[1], high = cell_colors[2]) + theme(legend.key.size = grid::unit(0.35, 
                "cm"), legend.text = element_text(size = 8))
            g8 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC8", fill = gene), size = dotsize, shape = 21, color = "black", alpha = alpha) + 
                guides(color = guide_legend(title = gene)) + scale_fill_gradient(low = cell_colors[1], high = cell_colors[2]) + theme(legend.key.size = grid::unit(0.35, 
                "cm"), legend.text = element_text(size = 8))
            g9 <- ggplot(comp.cell.type2) + geom_point(aes_string(x = refPC, y = "PC9", fill = gene), size = dotsize, shape = 21, color = "black", alpha = alpha) + 
                guides(color = guide_legend(title = gene)) + scale_fill_gradient(low = cell_colors[1], high = cell_colors[2]) + theme(legend.key.size = grid::unit(0.35, 
                "cm"), legend.text = element_text(size = 8))
            
            if (save == TRUE) {
                PCA_grob <- arrangeGrob(g1, g2, g3, g4, g5, g6, g7, g8, g9, ncol = 3)
                ggsave(PCA_grob, file = paste("PCA Matrix.tiff"))
            }
            
            quartz()
            PCA_grid <- grid.arrange(g1, g2, g3, g4, g5, g6, g7, g8, g9, ncol = 3)
            
        }
        
        
    }
    
    
    
    if (print == TRUE) {
        if (printNum > nrow(exprs(cellData))) {
            stop("There are fewer genes than printNum in the expression matrix.")
        }
        n <- printNum/2
        comp3 <- data.frame(PCA.allgenes$rotation[, 1:10])
        comp3_order <- apply(comp3, 2, order, decreasing = TRUE)
        terminal_output <- data.frame()
        for (i in 1:10) {
            terminal_output[1:nrow(comp3), i] <- row.names(comp3)[comp3_order[, i]]
        }
        colnames(terminal_output) <- colnames(comp3)
        terminal_output_final <- terminal_output[c(1:n, (nrow(terminal_output) - n):nrow(terminal_output)), ]
        print(terminal_output_final)
    }
    
} 
