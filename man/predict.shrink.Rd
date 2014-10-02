\name{predict.shrink}

\alias{predict.shrink}

\title{ Predict Method for Objects of Class \code{shrink} }

\description{ 
  Obtains predictions from shrunken regression coefficients from an object of 
  class \code{shrink}.
  This class of objects is returned by the \code{shrink} function. Objects of this 
  class have methods for the functions \code{coef}, \code{predict}, \code{print}, 
  \code{summary}, and \code{vcov}.
}

\usage{ 
\method{predict}{shrink}(object, newdata = NULL, type = c("link", "response", "lp", "risk",
   "expected", "terms"), terms = NULL, na.action = na.pass, collapse, 
   safe = FALSE, ...) 
}

\arguments{
  \item{object}{ an object of class \code{shrink}. }
  \item{newdata}{ a data frame for which predictions are obtained, otherwise 
                  predictions are based on the data stored in \code{object}. }
  \item{type}{ the type of prediction required. }
  \item{terms}{ with \code{type = "terms"} by default all terms are returned. A 
                character vector specifies which terms are to be returned. }
  \item{na.action}{	function determining what should be done with missing values 
                    in \code{newdata}. The default is to include all observations. }
  \item{collapse}{ if \code{family = coxph} or \code{Cox}, an optional vector of 
                   subject identifiers. If specified, the output will contain 
                   one entry per subject rather than one entry per observation. }
  \item{safe}{ option from \code{predict.mfp}. }
  \item{...}{ additional arguments to be passed to methods. }
}

%%\details{  }

\value{ A vector or matrix of predictions. 
%%\item{fit}{ predictions. 
}

%%\references{ XXX }

\author{ Daniela Dunkler, Georg Heinze }

%%\note{ xxx }

\seealso{ \code{\link{shrink}}, \code{\link{coef.shrink}}, \code{\link{print.shrink}}, 
          \code{\link{summary.shrink}}, \code{\link{vcov.shrink}} }

%%\examples{}

\keyword{ utilities }
