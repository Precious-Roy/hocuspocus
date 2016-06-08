#' The Hocus class
#'
#' An extension of Bioconductor's ExpressionSet class, which stores all single cell expression data.
#'
#' @export


setClass( "Hocus", 
          contains = "ExpressionSet",
          slots = c(colorData = "list",
                    logData = "AnnotatedDataFrame"),
          prototype = prototype( new( "VersionedBiobase",
                                      versions = c( Biobase:::classVersion("ExpressionSet"), Hocus = "1.0.0" ) ))
)
