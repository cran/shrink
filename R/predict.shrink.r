predict.shrink <-
  function(object, newdata = NULL, type = c("link", "response", "lp", "risk", "expected", "terms"), 
           terms = NULL, na.action = na.pass, collapse, safe = FALSE, ...)   
           # se.fit = FALSE, dispersion = NULL, 
{                                       
  # class glm = c("link", "response", "terms")
  # class mfp = c("link", "response", "lp", "risk", "expected", "terms"),
  # class coxph = c("lp", "risk", "expected", "terms")
  # class lm = c("response", "terms")
    
  if(!inherits(object, "shrink")) stop("'object' is not of class shrink")
  if(inherits(object$fit, "mfp")) { fit <- object$fit$fit } else { fit <- object$fit }         

  fit$coefficients <- object$ShrunkenRegCoef

  if(inherits(object$fit, "mfp")) {  
    # code from mfp:predict.mfp Version 1.4.9 (Nov. 2012)
    if(object$fit$family$family == "Cox") {
        if(is.null(object$fit$terms)) terms = names(object$fit$assign)
    	if(!missing(newdata)) 
    		if(!missing(collapse)) 
    			getFromNamespace("predict.coxph", "survival")(fit, newdata = newdata, type = type, 
                                                        se.fit = FALSE, terms = terms, 
                                                        collapse = collapse, safe = safe, ...) 
    		else
    			getFromNamespace("predict.coxph", "survival")(fit, newdata = newdata, type = type, 
                                                        se.fit = FALSE, terms = terms, 
                                                        safe = safe, ...) 
    	else
    		if(!missing(collapse)) 
    			getFromNamespace("predict.coxph", "survival")(fit, type = type, se.fit = FALSE, 
                                                        terms = terms, collapse = collapse, 
                                                        safe = safe, ...) 
    		else
    			getFromNamespace("predict.coxph", "survival")(fit, type = type, se.fit = FALSE, 
                                                        terms = terms, safe = safe, ...) 
    } else {
    	if(!missing(newdata)) 
    		predict.glm(fit, newdata = newdata, type = type, se.fit = FALSE, dispersion = NULL, 
                    terms = terms, na.action = na.action, ...)
    	else 
    		predict.glm(fit, type = type, se.fit = FALSE, dispersion = NULL, terms = terms, 
                    na.action = na.action, ...)
    }
  } else {
    if(inherits(object$fit, c("glm", "lm"))) {    
    	if(!missing(newdata)) { 
        predict.glm(fit, newdata = newdata, type = type, se.fit = FALSE, dispersion = NULL, 
                    terms = terms, na.action = na.action, ...)
      } else {
        predict.glm(fit, type = type, se.fit = FALSE, dispersion = NULL, terms = terms, 
                    na.action = na.action, ...)
      }
    } else {
#    if (inherits(object$fit, "coxph")) {           
      if(!is.null(terms)) { terms <- names(fit$assign) }
      if(is.null(newdata)) { newdata <- data.frame(fit$x) }      
      predict(fit, newdata = newdata, type = type, se.fit = FALSE, na.action = na.action, 
              terms = terms, collapse,  ...)
    }
  }
}
