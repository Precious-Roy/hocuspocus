#' Generate violin plots with a high density fit (i.e. blot plots).
#' 
#' Takes ExpressionSet object and generates blot plots for the specified groups 
#' and genes.
#' 
#' @param cellData ExpressionSet object created with readCells (and preferably 
#'   transformed with prepCells).
#' @param genes Vector of character strings specifying the names of genes to 
#'   plot.  All gene names must be genes that are present in the expression 
#'   matrix.
#' @param groups Character string specifying the title of the column in pData 
#'   that contains the names of the groups to which each sample belongs.  The 
#'   blots representing each group in the plots for each gene will be colored 
#'   accordingly. The column length should be the same length as the number of 
#'   samples.
#' @param colors  Character string specifying the title of the column in pData 
#'   that contains the colors to be used for plotting.  The number of colors 
#'   should be equal to the number of levels in groups.  The order of colors 
#'   corresponds to the order of factor levels in groups.  If nothing is 
#'   specified, random colors from the rainbow function are assigned.
#' @param cols Integer specifying the number of columns in which to arrange the 
#'   blot plots.  The number of cols should be less than or equal to the number 
#'   of genes.
#' @param singleGroup Character string specifying the name of a single group 
#'   within groups to plot.  Individual blot plots are generated for each gene 
#'   on a single set of axes.
#' @param singleColor Character string specifying the color to use for the 
#'   singleGroup blot plot.
#' @param order Boolean specifying whether the groups should be arranged in 
#'   alphabetical factor order prior to plotting.  If FALSE, the current order 
#'   of levels is used.
#' @param save Boolean specifying whether to save the resultant plot as a .tiff 
#'   file.
#' @return Blot plots for the specified groups and genes in a new window.
#' @export


blotPlot <- function(cellData, genes, groups = "GroupID", colors, cols = 3, singleGroup, singleColor = "red", order = TRUE, save = FALSE) {
    
    if (.Platform$OS.type == "windows") {
        quartz <- function() windows()
    }
    
    if (!("prepCells" %in% colnames(pData(cellData)))) {
        warning("It would be wise to run prepCells prior to blotPlot.", call. = FALSE)
    }
    
    if (missing(singleGroup) && missing(groups)) {
        stop("Either singleGroups or groups must be specified to generate a plot.", call. = FALSE)
    }
    
    if (!missing(singleGroup) && length(singleGroup) > 1) {
        stop("singleGroup should be length 1.", call. = FALSE)
    }
    
    log.data <- data.frame(exprs(cellData), stringsAsFactors = FALSE)
    
    samples <- pData(cellData)
    
    # Provide vector of genes you want to blot plot
    
    genes <- t(genes)
    
    if (!all(genes %in% row.names(log.data))) {
        stop("At least of the specified genes is not present in expression table", call. = FALSE)
    }
    
    
    # Get colors for groups
    
    if (missing(colors)) {
        group_colors <- palette(rainbow(length(levels(samples[, groups]))))
    }
    
    if (!missing(colors)) {
        group_colors <- data.frame(samples[, colors])
        group_colors <- as.character(group_colors[!apply(is.na(group_colors) | group_colors == "", 1, all), ])
        group_colors <- t(group_colors)
    }
    
    
    gene_refined_table <- log.data[genes, ]
    names(gene_refined_table) <- as.vector(samples[, groups])
    
    
    stacked_table <- stack(gene_refined_table)
    stacked_genes <- rep(row.names(gene_refined_table), length(row.names(samples)))
    
    if (!missing(singleGroup) && missing(groups)) {
        stop("A column of groups must be specified from a singleGroup can be chosen.", call. = FALSE)
    }
    if (!missing(singleGroup) && !missing(groups)) {
        plotme <- transform(stacked_table, group = gsub("\\..*$", "", stacked_table$ind), genez = stacked_genes)
        plotme$genez <- factor(plotme$genez)
        
        g1 <- ggplot(subset(plotme, (group == singleGroup)), aes(x = genez, y = values)) + geom_violin(adjust = 0.5, trim = F, 
            scale = "width", fill = singleColor, alpha = 0.3, color = "darkgray", width = 0.5) + geom_point(position = position_jitter(w = 0.1, 
            h = 0.1), size = 1.3, color = "grey26") + scale_fill_manual(values = group_colors) + theme(legend.position = "none", 
            axis.text.x = element_text(angle = 90, color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12), 
            axis.title.x = element_blank()) + labs(y = "Expression") + theme(axis.title.y = element_text(size = 17)) + theme(plot.background = element_blank(), 
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(), 
            strip.background = element_blank()) + theme(axis.line = element_line(color = "black"), axis.ticks.x = element_blank())
        if (save == TRUE) {
            ggsave(g1, filename = "DE Ute Diff Cells.tiff", height = 20, width = 8)
        }
    }
    
    if (!missing(groups) && missing(singleGroup)) {
        
        plotme <- transform(stacked_table, group = gsub("\\..*$", "", stacked_table$ind), genez = stacked_genes)
        plotme$group <- factor(plotme$group, levels = unique(plotme$group), ordered = FALSE)
        plotme$genez <- factor(plotme$genez, levels = unique(plotme$genez), ordered = FALSE)
        
        if (order == TRUE) {
          g <- plotme$group
          g <- factor(g[order(plotme$ind)],levels(g)[order(levels(g))])
          group_colors <- group_colors[order(levels(plotme$group))]
          plotme <- plotme[order(plotme$ind),c("values","ind","genez")]
          plotme$group <- g  
        }
        
        g1 <- ggplot(plotme, aes(x = group, y = values)) + geom_violin(adjust = 0.5, trim = F, scale = "width", aes(fill = group), 
            color = "black", alpha = 0.8, size = 0.5) + geom_point(position = position_jitter(w = 0.1, h = 0.1), size = 1.2, 
            color = "black") + scale_fill_manual(values = group_colors) + facet_wrap(~genez, nrow = round(length(genes)/cols)) + 
            theme(legend.position = "none", axis.text.x = element_text(angle = 50, hjust = 1, size = 15, color = "black"), axis.title.x = element_blank(), 
                axis.text.y = element_text(size = 12, color = "black"), strip.text.x = element_text(size = 13, face = "italic")) + 
            labs(y = "Expression") + theme(panel.background = element_rect(colour = "black", fill = "white"), strip.background = element_rect(color = "black", 
            size = 0.5), axis.title.y = element_text(size = 20, vjust = 1.5))
        if (save == TRUE) {
            ggsave(g1, filename = "DE Ute Diff Cells.tiff", height = 20, width = 8)
        }
        
    }
    
    quartz()
    g1
} 
