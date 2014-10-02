isrcs <- 
  function(fit) 
{ 
  # was the rms:rcs argument used in the fit?
  # [[1]] rcs is included in formula (right-hand side), but not as a variable name
  # [[2]] variables with rcs
  
  list(any(all.names(fit$call$formula[3]) == "rcs") && !any(all.vars(fit$call$formula[3]) == "rcs"),  
       regexpr(pattern = "rcs(", fixed =TRUE, text = names(fit$coefficients)) == TRUE)
}
