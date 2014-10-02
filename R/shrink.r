# R function 'shrink': obtain global, parameterwise, and joint post-estimation 
# shrinkage factors for regression coefficients from fit objects of class lm, 
# class glm with family=c("gaussian", "binomial"), class coxph, or class mfp 
# with family=c(cox, gaussian, binomial)
#
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
  function(fit, 
           type = c("parameterwise", "global", "all"), 
           method = c("jackknife", "dfbeta"), 
           join = NULL,
           notes = TRUE, 
           postfit = TRUE)
{
  if(!any(is.coxph <- inherits(fit, "coxph"), is.glm <- inherits(fit, "glm"), 
          is.lm <- inherits(fit, "lm")))     
      stop("'fit' is not of class 'coxph', 'glm' or 'lm'")
  
#  if(is.lm && length(method) == 2) method <- "dfbeta"
  type <- match.arg(type)
  method <- match.arg(method)
  
  if(!is.matrix(fit$x)) stop("recalculate 'fit' with 'x = TRUE'")
  if(is.null(fit$y)) stop("recalculate 'fit' with 'y = TRUE'")
  
  if(!is.null(join)) {
    if(!length(unlist(join)) == length(unique(unlist(join)))) stop("each variable is allowed only once in 'join'")
    if(type == "global") stop("argument 'join' is only applicable to 'type = parameterwise'")
  }
  
  if(is.coxph) { 
    if(type == "all") 
      results  <- list("global"        = shrink.coxph(fit, type = "global", method, join = NULL, notes, postfit),
                       "parameterwise" = shrink.coxph(fit, type = "parameterwise", method, join = NULL, notes, postfit),
                       "joint"         = if(!is.null(join)) shrink.coxph(fit, type = "parameterwise", method, join, notes, postfit))
    else results <- shrink.coxph(fit, type, method, join, notes, postfit)
  } else
  if(is.glm) {   
    if(type == "all") 
      results  <- list("global"        = shrink.glm(fit, type = "global", method, join = NULL, notes, postfit), 
                       "parameterwise" = shrink.glm(fit, type = "parameterwise", method, join = NULL, notes, postfit),
                       "joint"         = if(!is.null(join)) shrink.glm(fit, type = "parameterwise", method, join, notes, postfit))
    else results <- shrink.glm(fit, type, method, join, notes, postfit)
  } else {
    if(type == "all")
      results  <- list("global"        = shrink.lm(fit, type = "global", method, join = NULL, notes, postfit),
                       "parameterwise" = shrink.lm(fit, type = "parameterwise", method, join = NULL, notes, postfit),
                       "joint"         = if(!is.null(join)) shrink.lm(fit, type = "parameterwise", method, join, notes, postfit))
    else results <- shrink.lm(fit, type, method, join, notes, postfit)
  }

  if(type == "all" && is.null(join)) results[[3]] <- NULL
   
  results <- structure(c(results, fit = list(fit), type = type, method = method, 
                         join = list(join), call = match.call(expand.dots = TRUE)), 
                       class = "shrink")
  return(results)
}
