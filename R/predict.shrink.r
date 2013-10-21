predict.shrink <-
function(object, newdata = NULL, type = c("link", "response", "lp", "risk", "expected", "terms"), terms = NULL, 
         na.action = na.pass, collapse, safe = FALSE, ...)   # se.fit = FALSE, dispersion = NULL, 
{                                       
  # class glm = c("link", "response", "terms")
  # class mfp = c("link", "response", "lp", "risk", "expected", "terms"),
  # class coxph = c("lp", "risk", "expected", "terms")
  # class lm = c("response", "terms")
    
  if (class(object) != "shrink") { stop("object is not of class shrink") }

  if (sum(class(object$fit) %in% "mfp") == 1) { fit <- object$fit$fit } else
  if (sum(class(object$fit) %in% "mfp") == 0) { fit <- object$fit }         

  fit$coefficients <- object$shrunken

  # code from mfp:predict.mfp Version 1.4.9 (Nov. 2012)
  if (sum(class(object$fit) %in% "mfp") == 1) { 
    if (object$fit$family$family == "Cox") {
        if (is.null(object$fit$terms)) terms = names(object$fit$assign)
    	if (!missing(newdata)) 
    		if (!missing(collapse)) 
    			getFromNamespace("predict.coxph", "survival")(fit, newdata = newdata, type = type, se.fit = FALSE, terms = terms, collapse = collapse, safe = safe, ...) 
    		else
    			getFromNamespace("predict.coxph", "survival")(fit, newdata = newdata, type = type, se.fit = FALSE, terms = terms, safe = safe, ...) 
    	else
    		if (!missing(collapse)) 
    			getFromNamespace("predict.coxph", "survival")(fit, type = type, se.fit = FALSE, terms = terms, collapse = collapse, safe = safe, ...) 
    		else
    			getFromNamespace("predict.coxph", "survival")(fit, type = type, se.fit = FALSE, terms = terms, safe = safe, ...) 
    } else {
    	if (!missing(newdata)) 
    		predict.glm(fit, newdata = newdata, type = type, se.fit = FALSE, dispersion = NULL, terms = terms, na.action = na.action, ...)
    	else 
    		predict.glm(fit, type = type, se.fit = FALSE, dispersion = NULL, terms = terms, na.action = na.action, ...)
    }
  } else
  if (sum(class(object$fit) %in% "mfp") == 0) {
    if (sum(class(object$fit) %in% c("glm", "lm")) >= 1) {
    	if (!missing(newdata)) { 
        predict.glm(fit, newdata = newdata, type = type, se.fit = FALSE, dispersion = NULL, terms = terms, na.action = na.action, ...)
      } else {
        predict.glm(fit, type = type, se.fit = FALSE, dispersion = NULL, terms = terms, na.action = na.action, ...)
      }
    } else
    if (sum(class(object$fit) %in% "coxph") == 1) {
      if (!is.null(terms)) { terms <- names(fit$assign) }
      if (is.null(newdata)) { newdata <- data.frame(fit$x) }      
      predict(fit, newdata = newdata, type = type, se.fit = FALSE, na.action = na.action, terms = terms, collapse,  ...)
    }
  }
}
