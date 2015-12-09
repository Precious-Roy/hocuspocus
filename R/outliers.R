outliers <- function(){
  
  loggeomeans <- rowMeans(log(data))
  factors <- apply(data, 2, function(cnts) exp(median((log(cnts) - loggeomeans)[is.finite(loggeomeans)])))
  data <- t(t(data)/factors)
  
  outliers <- data.frame(pData(ute)$Outliers,stringsAsFactors=FALSE)
  outliers <- as.character(outliers[!apply(is.na(outliers) | outliers == "", 1, all), ])
  
  t<-!(row.names(pData(cellData)) %in% outliers)
  
  cellData <- cellData[,t]
  
}


