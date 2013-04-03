print.shrink <- 
function(x, ...)
{
  if (!class(x) %in% "shrink") { stop("object is not of class shrink") }

  cat(paste("Shrinkage Factors (type=", x$type, ", method=", x$method, "):\n", sep=""))
  if (x$type %in% "global") { print(as.vector(x$shrinkage)) } else { print(x$shrinkage) }
 
  if (!is.null(x$postfit)) { 
    cat("\nPostfit:\n")
    print(x$postfit, ...)
  } else {
    cat("\nShrunken Regression Coefficients:\n")
    print(x$shrunken, ...)
  }    
}   
