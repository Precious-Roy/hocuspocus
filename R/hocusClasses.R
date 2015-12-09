#'
#'
#'
#'
#' @export


setClass( "Hocus", 
          contains = "ExpressionSet",
          slots = c(colorData = "list",
                    logData = "AnnotatedDataFrame"),
          prototype = prototype( new( "VersionedBiobase",
                                      versions = c( classVersion("ExpressionSet"), Hocus = "1.0.0" ) ))
)
