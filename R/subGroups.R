#' Select a subset of groups in phenoData
#'
#' Select a subset of specified groups in the phenoData slot of a Hocus object while retaining associated color info in colorData.
#'
#' @export


subGroups <- function(cellData,groups,names){
  
  cellData<-cellData[,row.names(pData(cellData))[pData(cellData)[,groups]==names]]
  
  for (i in 1:length(cellData@colorData)){
    cellData@colorData[[i]] <- cellData@colorData[[i]][names(cellData@colorData[[i]])==unique(pData(cellData)[,names(cellData@colorData[i])])]
  }
  
  cellData
  
}
