coef.shrink <- 
  function(object, ...)
{
  if(!inherits(object, "shrink")) stop("'object' is not of class shrink")
  
  if(object$type == "all") {
    result <- vector(mode = "list", length = 2)
    names(result) <- c("global", "parameterwise")
    
    cat("GLOBAL SHRINKAGE\n")
    result[[1]] <- print(object$global$ShrunkenRegCoef)
    
    cat("\n\nPARAMETERWISE SHRINKAGE\n")
    result[[2]] <- print(object$parameterwise$ShrunkenRegCoef)
    
    if(!is.null(object$join)) {
      cat("\n\nPARAMETERWISE SHRINKAGE WITH JOIN OPTION\n")
      result$joint <- print(object$joint$ShrunkenRegCoef)
      if(!is.null(object$join)) cat("\njoint shrinkage was requested for:", sapply(object$join, function(x) paste(x, collapse="+")))
    }
  } else
    result <- print(object$ShrunkenRegCoef)
  
  invisible(result)
}   
