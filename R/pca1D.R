#'
#'
#'
#'
#' @export


pca1D <- function(cellData, center = TRUE, scale = FALSE, PC = 1, ICA=FALSE, groups, values, gene, colors, logNumeric = TRUE, alpha = 0.7, dotsize = 5, 
                    bubble = FALSE, bubbleSizes, save = FALSE) {
    
  if (.Platform$OS.type == "windows") {
    quartz <- function() windows()
  }
  
  if (cellData@logData$prepCells[1] == "No") {
    warning("It would be wise to run prepCells prior to pca2D.", call. = FALSE)
  }
  
  if (!missing(groups) && !all((groups %in% colnames(pData(cellData))))) {
    stop("The column name specified for groups is not found in phenoData. If absent, add a column to phenoData containing the group for each cell, even if it only consists of 1 group.", 
         call. = FALSE)
  }
  
  if (!missing(values) && !all((values %in% colnames(pData(cellData))))) {
    stop("The column name specified for values is not found in phenoData.", 
         call. = FALSE)
  }
  
  if (!missing(groups) && !all((names(cData(cellData)) %in% colnames(pData(cellData))))) {
    stop("Some or all of the column names in colorData do not have matching column names in phenoData.", 
         call. = FALSE)
  }

    
  samples <- pData(cellData)
  exprs_table <- exprs(cellData)
  
  # Obtain cell colors
  
  if (ICA==FALSE){
    pca1 <- prcomp(t(exprs_table), center = center, scale. = scale)
    comp1 <- data.frame(pca1$x[, PC])
    comp1["y"] <- seq(0, 0, length.out = length(samples[, groups]))
    comp.cell.type1 <- comp1
    
  }
  
  if (ICA==TRUE){
    ica1 <- fastICA::fastICA(t(exprs_table), PC)
    comp1 <- data.frame(ica1$S[, PC])
    comp1["y"] <- seq(0, 0, length.out = length(samples[, groups]))
    comp.cell.type1 <- comp1
  }
  
  if (missing(groups)) {
    cell_colors <- "blue"
    comp.cell.type1["Cell_Type"] <- rep("All", nrow(pData(cellData)))
  }
  
  if (!missing(groups)) {
    cell_colors <- data.frame(cellData@colorData[[groups]], stringsAsFactors = FALSE)
    cell_colors <- as.character(cell_colors[!apply(is.na(cell_colors) | cell_colors == "", 1, all), ])
    groupz <- data.frame(pData(cellData)[, groups], stringsAsFactors = FALSE)
    comp.cell.type1["Cell_Type"] <- factor(groupz[, colnames(groupz)], levels = unique(groupz[, colnames(groupz)]), ordered = FALSE)
  }
  
  
  if (bubble == TRUE) {
    
    if (missing(bubbleSizes)) {
      sizes <- 1.4^((4:(4 + length(levels(comp.cell.type1$Cell_Type)))))
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
  
      samp_PCA <- ggplot(comp.cell.type1) + geom_point(aes_string(x = colnames(comp1)[1], y = colnames(comp1)[2], fill = "factor(Cell_Type,levels=unique(Cell_Type),ordered=FALSE)", 
      size = "factor(Cell_Type,levels=unique(Cell_Type),ordered=FALSE)"), shape = 21, color = "black", alpha = alpha, position = position_jitter(height = 2)) + scale_fill_manual(name = "Groups", 
      values = cell_colors) + scale_size_manual(name = "Groups", values = sizes) + guides(fill = guide_legend("Groups"), color = guide_legend("Groups"), 
      size = guide_legend("Groups")) + theme_bw() + theme(legend.key = element_rect(colour = F), legend.background = element_rect(fill = NA)) + theme(legend.key = element_rect(colour = F), legend.background = element_rect(fill = NA)) + theme(axis.text.x = element_text(size = 25), 
      axis.text.y = element_blank(), axis.title.x = element_text(size = 35, vjust = -1), axis.title.y = element_blank(), legend.text = element_text(size = 20), 
      legend.title = element_text(size = 30, face = "bold"), plot.margin = grid::unit(c(2.1, 1.1, 1.1, 1.1), "cm"), legend.key.size = grid::unit(0.9, "cm"), axis.ticks.y = element_blank()) + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylim(-10, 10) + xlab(paste("PC", PC)) + scale_color_gradient2()
      
      
      quartz()
      print(samp_PCA)
      if (save == TRUE) {
        ggsave(samp_PCA, file = paste("PCA Plot.tiff"), width = 7, height = 6)
      }
    
  }
  
  if (!missing(values)) {
    
    if (!missing(gene)) {
      stop("Cannot specify both values and gene.", call. = FALSE)
    }
    
    if (!missing(colors)){
      cell_colors <- colors
    }
    
    if (length(cell_colors) != 2) {
      warning("Specify a vector of two colors to control color, using default colors instead.", call. = FALSE)
      cell_colors[1:2] <- c("black", "yellow")
    }
    
    valuez <- data.frame(pData(cellData)[, values], stringsAsFactors = FALSE)
    if (logNumeric == TRUE) {
      valuez <- log2(valuez)
    }
    comp.cell.type1["Values"] <- valuez
    
    if (bubble == TRUE) {
      
      samp_PCA <- ggplot(comp.cell.type1) + geom_point(aes_string(x = colnames(comp1)[1], y = colnames(comp1)[2], fill = "Values", size = "factor(Cell_Type,levels=unique(Cell_Type),ordered=FALSE)"), 
        shape = 21, color = "black", alpha = alpha, position = position_jitter(height = 2)) + scale_size_manual(name = "Groups", values = sizes) + guides(color = guide_legend("Values"), size = guide_legend("Groups")) + 
        theme_bw() + theme(legend.key = element_rect(colour = F), legend.background = element_rect(fill = NA)) + theme(legend.key = element_rect(colour = F), legend.background = element_rect(fill = NA)) + theme(legend.key = element_rect(colour = F), legend.background = element_rect(fill = NA)) + theme(axis.text.x = element_text(size = 25), 
        axis.text.y = element_blank(), axis.title.x = element_text(size = 35, vjust = -1), axis.title.y = element_blank(), legend.text = element_text(size = 20), 
        legend.title = element_text(size = 30, face = "bold"), plot.margin = grid::unit(c(2.1, 1.1, 1.1, 1.1), "cm"), legend.key.size = grid::unit(0.9, "cm"), axis.ticks.y = element_blank()) + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylim(-10, 10) + xlab(paste("PC", PC)) + scale_fill_gradient(low = cell_colors[1], mid = muted(cell_colors[2]), high = cell_colors[2], mid=max(valuez)/2) 
                                                                                                                                                                                                                                                 
      
      quartz()
      print(samp_PCA)
      if (save == TRUE) {
        ggsave(samp_PCA, file = paste("PCA Plot.tiff"), width = 7, height = 6)
      }
    }
    
    if (bubble == FALSE) {
      
      samp_PCA <- ggplot(comp.cell.type1) + geom_point(aes_string(x = colnames(comp1)[1], y = colnames(comp1)[2], fill = "Values"), size = dotsize, shape = 21, 
      color = "black", alpha = alpha, position = position_jitter(height = 2)) + guides(color = guide_legend(title = "Values")) + theme_bw() + theme(legend.key = element_rect(colour = F), 
      legend.background = element_rect(fill = NA)) + theme(legend.key = element_rect(colour = F), legend.background = element_rect(fill = NA)) + theme(legend.key = element_rect(colour = F), legend.background = element_rect(fill = NA)) + theme(axis.text.x = element_text(size = 25), 
      axis.text.y = element_blank(), axis.title.x = element_text(size = 35, vjust = -1), axis.title.y = element_blank(), legend.text = element_text(size = 20), 
      legend.title = element_text(size = 30, face = "bold"), plot.margin = grid::unit(c(2.1, 1.1, 1.1, 1.1), "cm"), legend.key.size = grid::unit(0.9, "cm"), axis.ticks.y = element_blank()) + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylim(-10, 10) + xlab(paste("PC", PC)) + scale_fill_gradient(low = cell_colors[1], mid = muted(cell_colors[2]), high = cell_colors[2], mid=max(valuez)/2)
      
      quartz()
      print(samp_PCA)
      if (save == TRUE) {
        ggsave(samp_PCA, file = paste("PCA Plot.tiff"), width = 7, height = 6)
      }
    }
  }
  
  if (!missing(gene)) {
    
    if (!missing(values)) {
      stop("Cannot specify both values and gene.", call. = FALSE)
    }
    
    
    if (!missing(gene) && !gene %in% row.names(exprs_table)) {
      stop("Specified gene is not present in expression table", call. = FALSE)
    }
    
    if (!missing(colors)){
      cell_colors <- colors
    }
    
    if (length(cell_colors) != 2) {
      warning("Specify a vector of two colors to control color, using default colors instead", call. = FALSE)
      cell_colors[1:2] <- c("black", "yellow")
    }
    
    genez <- data.frame(exprs_table[gene, ], stringsAsFactors = FALSE)[, 1]
    comp.cell.type1[gene] <- genez
    
    if (bubble == TRUE) {
      
      samp_PCA <- ggplot(comp.cell.type1) + geom_point(aes_string(x = colnames(comp1)[1], y = colnames(comp1)[2], fill = gene, size = "factor(Cell_Type,levels=unique(Cell_Type),ordered=FALSE)"), 
      shape = 21, color = "black", alpha = alpha, position = position_jitter(height = 2)) + scale_size_manual(name = "Groups", values = sizes) + guides(color = guide_legend(title = gene), 
      size = guide_legend("Groups")) + theme_bw() + theme(legend.key = element_rect(colour = F), legend.background = element_rect(fill = NA)) + theme(legend.key = element_rect(colour = F), legend.background = element_rect(fill = NA)) + theme(axis.text.x = element_text(size = 25), 
      axis.text.y = element_blank(), axis.title.x = element_text(size = 35, vjust = -1), axis.title.y = element_blank(), legend.text = element_text(size = 20), 
      legend.title = element_text(size = 30, face = "bold"), plot.margin = grid::unit(c(2.1, 1.1, 1.1, 1.1), "cm"), legend.key.size = grid::unit(0.9, "cm"), axis.ticks.y = element_blank()) + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylim(-10, 10) + xlab(paste("PC", PC)) + scale_fill_gradient(low = cell_colors[1], mid = muted(cell_colors[2]), high = cell_colors[2], mid=max(exprs_table[gene, ])/2)
      
      quartz()
      print(samp_PCA)
      if (save == TRUE) {
        ggsave(samp_PCA, file = paste("PCA Plot.tiff"), width = 7, height = 6)
      }
    }
    
    if (bubble == FALSE) {
      
      samp_PCA <- ggplot(comp.cell.type1) + geom_point(aes_string(x = colnames(comp1)[1], y = colnames(comp1)[2], fill = gene), size = dotsize, shape = 21, 
      color = "black", alpha = alpha, position = position_jitter(height = 2)) + guides(color = guide_legend(title = gene)) + theme_bw() + theme(legend.key = element_rect(colour = F), legend.background = element_rect(fill = NA)) + theme(legend.key = element_rect(colour = F), legend.background = element_rect(fill = NA)) + theme(axis.text.x = element_text(size = 25), 
      axis.text.y = element_blank(), axis.title.x = element_text(size = 35, vjust = -1), axis.title.y = element_blank(), legend.text = element_text(size = 20), 
      legend.title = element_text(size = 30, face = "bold"), plot.margin = grid::unit(c(2.1, 1.1, 1.1, 1.1), "cm"), legend.key.size = grid::unit(0.9, "cm"), axis.ticks.y = element_blank()) + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylim(-10, 10) + xlab(paste("PC", PC)) + scale_fill_gradient(low = cell_colors[1], mid = muted(cell_colors[2]), high = cell_colors[2], mid=max(exprs_table[gene, ])/2) 
      
      quartz()
      print(samp_PCA)
      if (save == TRUE) {
        ggsave(samp_PCA, file = paste("PCA Plot.tiff"), width = 7, height = 6)
      }
    }
  }
  samp_PCA
} 
