coef.shrink <- 
function(object, ...)
{
  if (!class(object) %in% "shrink") { stop("object is not of class shrink") }
  
  print(object$shrunken)
}   
