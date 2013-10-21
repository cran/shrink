# R function 'shrink': obtain global, parameterwise, and joint post-estimation 
# shrinkage factors for regression coefficients from fit objects of class lm, 
# class glm with family=c("gaussian", "binomial"), class coxph, or class mfp 
# with family=c(cox, gaussian, binomial)
# Copyright (C) 2013 Daniela Dunkler, Georg Heinze
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License version 2
# as published by the Free Software Foundation
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

shrink <-
function(fit, type = "parameterwise", method = "jackknife", join = NULL)
{
  if (method == "d") { method <- "dfbeta" }
  if ((method == "j") | (method == "jack")) { method <- "jackknife" }  
  if ((method != "dfbeta") & (method != "jackknife")) { stop("method must be either jackknife or dfbeta") }
  
  if (type == "g") { type <- "global" }
  if ((type == "p") | (type == "pw")) { type <- "parameterwise" }  
  if ((type != "global") & (type != "parameterwise")) { stop("type must be either global or parameterwise") }

  if (is.null(fit$x) == TRUE | !is.matrix(fit$x) | ((class(fit)[1] %in% "coxph")) & !is.matrix(fit$x)) { 
    stop("recalculate the fit with x=TRUE") 
  }
  if (is.null(fit$y) == TRUE) { stop("recalculate the fit with y=TRUE") }
  

  varnames <- names(fit$coefficients)                                           # for rcs
  rcs.test <- sapply(varnames, isrcs)                                         
  if (sum(rcs.test) >= 1) {
    varnames.org <- varnames
    index <- ((1:length(varnames))[as.vector(rcs.test)])
    pre <- strsplit(varnames[index[1]], ")")[[1]][2]
    I <- 1
    for (i in 1:length(index))  {
      if (i != 1) { if (pre != strsplit(strsplit(varnames[index[i]], ")")[[1]][2], "'")[[1]][1]) { 
        I <- 1 
        pre <- strsplit(varnames[index[i]], ")")[[1]][2]
      } }
      varnames[index[i]] <- paste("rcs", I, ".", strsplit(varnames[index[i]], ")")[[1]][2], sep="")
      varnames[index[i]] <- strsplit(varnames[index[i]], "'")[[1]][1]
      I <- I + 1
    }
    names(fit$coefficients) <- dimnames(fit$x)[[2]] <- varnames 
    if (!is.null(join)) { 
      if (sum(all.vars(fit$call$formula[[3]]) %in% join) >= 1) {
        index2 <- join %in% all.vars(fit$call$formula[[3]])                     # rcs(var) is the only element in join level
        for (i in (1:length(index2))[index2]) { join[[i]] <- varnames[grep(join[[i]], varnames)] }
      }
    }
  }


  if (!is.null(join)) {
    if (!length(unlist(join)) == length(unique(unlist(join)))) { stop("each variable is allowed only once in join") }
    if (length(unlist(join)) == 1) { stop("need at least 2 variables in join") }
    if (!sum(dimnames(fit$x)[[2]] %in% unlist(join)) == length(unlist(join))) { stop("variables in fit$x and variables listed in join do not match") }
    if (type == "global") { stop("argument join is only applicable to type=parameterwise") }
  }
  
  
  n <- nrow(fit$x)
  callT <- match.call(expand.dots = TRUE) 


  # class: coxph, mfp (with family = cox)
  if (sum(class(fit) %in% "coxph") == 1) {
    intercept <- FALSE
        
    if (method == "dfbeta") {
      dfb <- residuals(fit, "dfbeta")
      bi <- rep(coefficients(fit), each = n) - dfb
    } else
    if (method == "jackknife") {
      bi <- matrix(0, n, ncol(fit$x))
      for(i in 1:n) bi[i,] <- coxph(fit$y[-i,] ~ fit$x[-i,])$coefficients
    }

    if (type == "parameterwise") {
      if (is.null(join))  { shrinkage <- fit$x * bi } else
      if (!is.null(join)) {
        shrinkage <- matrix(NA, nrow = n, ncol = ncol(fit$x) - length(unlist(join)) + length(join))
        Joined <- rep(0, ncol(fit$x))
        index <- rep(0, times = ncol(fit$x))
        names(index) <- dimnames(fit$x)[[2]]
        for (i in 1:length(join)) {
          joined <- dimnames(fit$x)[[2]] %in% join[[i]]
          index[joined] <- i
          Joined <- Joined + joined
          shrinkage[,i] <- sapply(1:n, function(X) fit$x[X, joined] %*% bi[X, joined])
        }
        if (sum(!dimnames(fit$x)[[2]]%in%unlist(join)) >=1 ) {                  # shrinkage factors for variables not listed in join
          shrinkage[, -c(1:length(join))] <- fit$x[, !Joined] * bi[, !Joined] 
          nam <- dimnames(fit$x)[[2]]
          dimnames(shrinkage)[[2]] <- c(paste("join", unlist(lapply(join, function(X) X[1])), sep="."), nam[!Joined])
        } else
        if (sum(!dimnames(fit$x)[[2]] %in% unlist(join))==0) { 
          dimnames(shrinkage) <- list(1:nrow(shrinkage), paste("join", unlist(lapply(join, function(X) X[1])), sep=".")) 
        }
      }
      if (is.vector(shrinkage)) { shrinkage <- matrix(shrinkage, nrow = n) } 
      if (is.matrix(shrinkage) & ncol(shrinkage)==1 & is.null(dimnames(shrinkage))) { 
        shrinkage <- matrix(shrinkage, nrow = n, dimnames = list(1:n, dimnames(fit$x)[[2]])) 
      }
      f <- coxph(formula(paste("fit$y~", paste(dimnames(shrinkage)[[2]], collapse="+"), sep="")), 
           data = data.frame(shrinkage))
    } else
    if (type == "global") {
      if (is.matrix(bi))  { shrinkage <- sapply(1:n, function(X) fit$x[X,] %*% bi[X,]) }
      if (!is.matrix(bi)) { shrinkage <- sapply(1:n, function(X) fit$x[X,] %*% bi[X]) }
      f <- coxph(fit$y ~ shrinkage)
    }
    res <- f$coefficients
    cov <- vcov(f)    

    if (type == "parameterwise" & (!is.null(join))) {                           # sort res
      res2 <- c(names(fit$coefficients)[!Joined], paste("join", rep(unlist(lapply(join, function(x) { x[1] } )), 
                unlist(lapply(join, function(x) { length(x) } ))), sep="."))
      res2 <- res[match(res2, names(res))]
      if ((sum(Joined == 0)) == 0) { names(res2) <- unlist(join) } else {
      names(res2)[-c(0:sum(Joined == 0))] <- unlist(join) }
      res <- res2[names(fit$coefficients)]
    }

    shr <- res * fit$coefficients

    if (type == "parameterwise" & is.null(join)) {                              # postfit 
      betamat <- matrix(rep(coefficients(fit), each = n), n, ncol(fit$x))
      modx <- bi/betamat * fit$x
      colnames(modx) <- names(res)

      modx <- cbind(fit$y, modx)
      f2 <- formula(paste(paste("Surv(", dimnames(modx)[[2]][1], ", ", dimnames(modx)[[2]][2], ")~", sep=""), 
                    paste(names(res), collapse="+"), sep=""))
      fit2 <- coxph(f2, data = data.frame(modx), x = TRUE)
      fit2$call$formula <- f2
    } 
  } else

  # class: lm, glm (with family = gaussian, binomial), mfp (with family = gaussian, binomial)
  if (sum(class(fit) %in% c("glm", "lm")) >= 1) {
    if (sum(class(fit) %in% "glm") == 1) { if (!fit$family[[1]] %in% c("gaussian", "binomial")) { 
      stop(paste("family", fit$family[[1]], "is not supported", sep=" ")) } 
    }
    
    if (sum(varnames %in% c("(Intercept)", "Intercept")) == 1) { intercept <- TRUE } else intercept <- FALSE
    varnames <- varnames[!(varnames %in% c("(Intercept)", "Intercept"))]

    if (sum(class(fit) %in% "mfp") == 1) {                                      # mfp inconsistency
      if (!identical(as.vector(fit$coefficients),as.vector(fit$fit$coefficients))) { 
        cat("Note, that coef(fit) != coef(fit$fit) in the mfp object. shrink uses coef(fit$fit).\n\n")
        fit.save <- fit
        varnames2 <- names(fit$fit$coefficients)
        varnames2 <- varnames2[!(varnames2 %in% c("(Intercept)", "Intercept"))]
  
        names(fit$fit$coefficients) <-  names(fit$coefficients)
        
        fit$coefficients <- fit$fit$coefficients
      }
    } else { varnames2 <- varnames }
    
    if (method == "dfbeta") {                                                   # intercept is not considered
      if (intercept) { dfb <- dfbeta(fit)[, -1] } else                          
      if (!intercept) { dfb <- dfbeta(fit) }                                    
      bi <- rep(fit$coefficients[varnames], each = n) - dfb
    } else
    if (method == "jackknife") {
      bi <- matrix(NA, n, length(varnames))
      if (intercept) {
        if (sum(class(fit) %in% "glm") == 1)  { 
          for(i in 1:n) bi[i,] <- glm(fit$y[-i] ~ fit$x[-i, varnames], family = fit$family[[1]])$coefficients[-1] 
        }
        if (!sum(class(fit) %in% "glm") == 1) { 
          for(i in 1:n) bi[i,] <- lm(fit$y[-i] ~ fit$x[-i, varnames])$coefficients[-1] 
        }
      } else
      if (!intercept) {
        if (sum(class(fit) %in% "glm") == 1)  { 
          for(i in 1:n) bi[i,] <- glm(fit$y[-i] ~ -1 + fit$x[-i, varnames], family = fit$family[[1]])$coefficients 
        }
        if (!sum(class(fit) %in% "glm") == 1) { 
          for(i in 1:n) bi[i,] <- lm(fit$y[-i] ~ -1 + fit$x[-i, varnames])$coefficients 
        }
      }
    }

    if (type == "parameterwise") {
      if (is.null(join))  { shrinkage <- fit$x[, varnames] * bi } else
      if (!is.null(join)) {
        shrinkage <- matrix(NA, nrow = n, ncol = length(varnames) - length(unlist(join)) + length(join))
        Joined <- rep(0, length(varnames))
        index <- rep(0, times = length(varnames))
        names(index) <- varnames
        for (i in 1:length(join)) {
          joined <- varnames %in% join[[i]]
          index[joined] <- i
          Joined <- Joined + joined
          if (intercept)  { shrinkage[,i] <- sapply(1:n, function(X) fit$x[X, c(FALSE, joined)] %*% bi[X, joined]) } else
          if (!intercept) { shrinkage[,i] <- sapply(1:n, function(X) fit$x[X, joined] %*% bi[X, joined]) } 
        }
        if (sum(!varnames%in%unlist(join)) >= 1) {
          if (intercept)  { shrinkage[, -c(1:length(join))] <- fit$x[, c(FALSE, !Joined)] * bi[, !Joined] } else
          if (!intercept) { shrinkage[, -c(1:length(join))] <- fit$x[, !Joined] * bi[, !Joined] }
          dimnames(shrinkage)[[2]] <- c(paste("join", unlist(lapply(join, function(X) X[1])), sep="."), varnames[!Joined])
        } else
        if (sum(!varnames %in% unlist(join))==0) { 
          dimnames(shrinkage) <- list(1:nrow(shrinkage), paste("join", unlist(lapply(join, function(X) X[1])), sep=".")) 
        }
      }
      if (is.vector(shrinkage)) { shrinkage <- matrix(shrinkage, nrow = n) } 
      if (is.matrix(shrinkage) & ncol(shrinkage)==1 & is.null(dimnames(shrinkage))) { 
        shrinkage <- matrix(shrinkage, nrow = n, dimnames = list(1:n, varnames)) 
      }
  
      if (intercept) {      
        if (sum(class(fit) %in% "glm") == 1)  { 
          f <- glm(formula(paste("fit$y~", paste(dimnames(shrinkage)[[2]], collapse = "+"), sep = "")), 
                   family = fit$family[[1]], data=data.frame(shrinkage))
        }    
        if (!sum(class(fit) %in% "glm") == 1) { 
          f <- glm(formula(paste("fit$y~", paste(dimnames(shrinkage)[[2]], collapse="+"), sep="")),  
                   family = "gaussian", data = data.frame(shrinkage))
        }    
      } else
      if (!intercept) {      
        if (sum(class(fit) %in% "glm") == 1)  { 
          f <- glm(formula(paste("fit$y~ -1 + ", paste(dimnames(shrinkage)[[2]], collapse="+"), sep = "")), 
                   family = fit$family[[1]], data = data.frame(shrinkage))
        }    
        if (!sum(class(fit) %in% "glm") == 1) { 
          f <- glm(formula(paste("fit$y~ -1 + ", paste(dimnames(shrinkage)[[2]], collapse="+"), sep="")), 
                   family = "gaussian", data = data.frame(shrinkage))
        }    
      } 
    } else
    if (type == "global") {
      if (is.matrix(bi))  { shrinkage <- sapply(1:n, function(X) fit$x[X, varnames] %*% bi[X,]) }
      if (!is.matrix(bi)) { shrinkage <- sapply(1:n, function(X) fit$x[X, varnames] %*% bi[X]) }

      if (intercept) { 
        if (sum(class(fit) %in% "glm") == 1)  { f <- glm(fit$y ~ shrinkage, family = fit$family[[1]]) }
        if (!sum(class(fit) %in% "glm") == 1) { f <- glm(fit$y ~ shrinkage, family = "gaussian") }
      } else 
      if (!intercept) { 
        if (sum(class(fit) %in% "glm") == 1)  { f <- glm(fit$y ~ -1 + shrinkage, family = fit$family[[1]]) }
        if (!sum(class(fit) %in% "glm") == 1) { f <- glm(fit$y ~ -1 + shrinkage, family = "gaussian") }
      }      
    }
    
    if (intercept) { res <- f$coefficients[-1] } else if (!intercept) { res <- f$coefficients }  # intercept is not shown in the results
    cov <- vcov(f)

    if (type == "parameterwise" & (!is.null(join))) {                           # sort res
      res2 <- c(varnames[!Joined], paste("join", rep(unlist(lapply(join, function(x) { x[1] } )), 
                unlist(lapply(join, function(x) { length(x) } ))), sep="."))
      res2 <- res[match(res2, names(res))]
      if ((sum(Joined == 0)) == 0) { names(res2) <- unlist(join) } else {
      names(res2)[-c(0:sum(Joined == 0))] <- unlist(join) }
      res <- res2[varnames]
    } 
    
    if (type != "parameterwise" | (type == "parameterwise" & !is.null(join)) | 
       (type == "parameterwise" & method=="jackknife") |
       (method == "dfbeta" & type == "parameterwise" & is.null(join) & intercept)) {
      shr <- rep(NA, length = length(fit$coefficients))
      names(shr) <- names(fit$coefficients)
      if (intercept) {
        shr[-1] <- res * fit$coefficients[-1]
        lp <- apply(t(t(fit$x[,-1]) * shr[-1]), 1, sum)                         # adapt intercept
        if (sum(class(fit) %in% "glm") == 1)  { 
          shr[1] <- glm(fit$y ~ offset(lp), family = fit$family[[1]])$coefficients[1] 
        }
        if (!sum(class(fit) %in% "glm") == 1) { 
          shr[1] <- glm(fit$y ~ offset(lp), family = "gaussian")$coefficients[1] 
        }
      } else
      if (!intercept) { shr <- res*fit$coefficients }
    }
    
    if (type == "parameterwise" & is.null(join) & !intercept) {                 # postfit
      betamat <- matrix(rep(fit$coefficients[varnames], each = n), n, length(varnames))
      modx <- bi / betamat * fit$x[,varnames]
      modx <- cbind(fit$y, modx)
      dimnames(modx)[[2]] <- c(all.vars(fit$call$formula)[1], names(res))
      if (intercept) {
        f2 <- formula(paste(paste(dimnames(modx)[[2]][1], "~", sep=""), paste(names(res), 
              collapse = "+"), sep = "")) 
      } else
      if (!intercept) {
        f2 <- formula(paste(paste(dimnames(modx)[[2]][1], "~ -1 +", sep=""), paste(names(res), 
              collapse = "+"), sep = "")) 
      }
      if (sum(class(fit) %in% "glm") == 1) { 
        fit2 <- glm(f2, family=fit$family[[1]], data=data.frame(modx)) 
        fit2$call$family <- fit$family[[1]]
      }
      if (!sum(class(fit) %in% "glm") == 1) {                                   # for class lm
        fit2 <- glm(f2, family = "gaussian", data = data.frame(modx)) 
      }
      fit2$call$formula <- f2
      shr <- fit2$coefficients
    }
  } else { stop(paste("class", class(fit), "is not supported", sep = " ")) }
  
  
  if ((sum(class(fit) %in% "mfp") == 1) & (sum(!(class(fit) %in% "coxph")) == 0)) {
    if (!identical(varnames, varnames2)) { fit <- fit.save }
  }
  
  if (type == "parameterwise" & is.null(join) & !intercept) { 
    res <- list(shrinkage = res, vcov.shrinkage = cov, shrunken = shr, 
           postfit = fit2, fit = fit, type = type, method = method, call = callT) 
  }
  
  if (!is.list(res)) { 
    res <- list(shrinkage = res, vcov.shrinkage = cov, 
           shrunken = shr, fit = fit, type = type, method = method, call = callT) 
  }
  res <- structure(res, class = "shrink")
  return(res)
}
