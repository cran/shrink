isrcs <- 
function(x) 
{ 
  # is the rms:rcs argument used in the respective model? (T/F)
  regexpr(pattern = "rcs", text = x) == TRUE 
}
