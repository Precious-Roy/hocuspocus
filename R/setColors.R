#'
#'
#'
#'
#' @export


setColors <- function(cellData, groups, colors, auto=FALSE, palette="rainbow"){
  
  if (missing(groups) && setAll==FALSE){
    stop("No group headers are specified and setAll is FALSE, no colors will be set.", call. = FALSE)
  }
  
  if (!all(groups %in% colnames(pData(cellData)))) {
    stop("All names in groups must match the column names in the phenoData matrix.",call.=FALSE)
  }
  
  if (!missing(colors) && auto==FALSE){
    if (!all(lengths(colors) == apply(data.frame(pData(cellData)[,groups]), 2, function(x) nlevels(as.factor(x))))){
      stop("The number of colors for each group in the colorData file must match the number of levels within the corresponding columns of the phenoData matrix.",call.=FALSE)
    }
    for (i in 1:length(colors)){
      cellData@colorData[[groups[i]]] <- colors[[i]]
      names(cellData@colorData[[groups[i]]]) <- unique(cellData@phenoData@data[,groups[i]])
    }
    
  }
  
  if (auto==TRUE){
    
    if (!missing(colors)){
      warning("Manual colors are specified and auto is TRUE, defaulting to auto colors.", call. = FALSE)
    }
    
    if (palette != "rainbow" && palette != "heat" && palette != "terrain" && palette != "topo" && palette != "cm" && palette != "greys"){
      stop("palette must be one of rainbow, heat, terrain, topo, cm", call.=FALSE)
    }
    
    if (palette == "rainbow"){
      for (i in 1:length(groups)){
        cellData@colorData[[groups[i]]] <- rainbow(length(unique(cellData@phenoData@data[,groups[i]])))[order(order(unique(as.vector(cellData@phenoData@data[,groups[i]]))))]
        names(cellData@colorData[[groups[i]]]) <- unique(cellData@phenoData@data[,groups[i]])
      }
    }
    
    if (palette == "heat"){
      for (i in 1:length(groups)){
        cellData@colorData[[groups[i]]] <- heat.colors(length(unique(cellData@phenoData@data[,groups[i]])))[order(order(unique(as.vector(cellData@phenoData@data[,groups[i]]))))]
        names(cellData@colorData[[groups[i]]]) <- unique(cellData@phenoData@data[,groups[i]])
      }
    }
    
    if (palette == "terrain"){
      for (i in 1:length(groups)){
        cellData@colorData[[groups[i]]] <- terrain.colors(length(unique(cellData@phenoData@data[,groups[i]])))[order(order(unique(as.vector(cellData@phenoData@data[,groups[i]]))))]
        names(cellData@colorData[[groups[i]]]) <- unique(cellData@phenoData@data[,groups[i]])
      }
    }
    
    if (palette == "topo"){
      for (i in 1:length(groups)){
        cellData@colorData[[groups[i]]] <- topo.colors(length(unique(cellData@phenoData@data[,groups[i]])))[order(order(unique(as.vector(cellData@phenoData@data[,groups[i]]))))]
        names(cellData@colorData[[groups[i]]]) <- unique(cellData@phenoData@data[,groups[i]])
      }
    }
    
    if (palette == "cm"){
      for (i in 1:length(groups)){
        cellData@colorData[[groups[i]]] <- cm.colors(length(unique(cellData@phenoData@data[,groups[i]])))[order(order(unique(as.vector(cellData@phenoData@data[,groups[i]]))))]
        names(cellData@colorData[[groups[i]]]) <- unique(cellData@phenoData@data[,groups[i]])
      }
    }
    
    if (palette == "greys"){
      for (i in 1:length(groups)){
        cellData@colorData[[groups[i]]] <- grey.colors(length(unique(cellData@phenoData@data[,groups[i]])),start=0.9,end=0)[order(order(unique(as.vector(cellData@phenoData@data[,groups[i]]))))]
        names(cellData@colorData[[groups[i]]]) <- unique(cellData@phenoData@data[,groups[i]])
      }
    }
    
  }
  
  cellData
  
}
