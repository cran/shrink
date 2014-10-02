\name{coef.shrink}

\alias{coef.shrink}

\title{ Returns Shrunken Regression Coefficients from Objects of Class \code{shrink} }

\description{ 
  This class of objects is returned by the \code{shrink} function. Objects of this 
  class have methods for the functions \code{coef}, \code{predict}, \code{print}, 
  \code{summary}, and \code{vcov}.
}

\usage{                 
\method{coef}{shrink}(object, ...) 
}

\arguments{
  \item{object}{ object of class \code{shrink}. }
  \item{...}{ further arguments. }  
}

%%\details{ }

\value{ A vector with shrunken regression coefficients. }

%%\references{  }

\author{ Daniela Dunkler, Georg Heinze }

%%\note{  }

\seealso{ \code{\link{shrink}}, \code{\link{print.shrink}}, \code{\link{predict.shrink}}, 
          \code{\link{summary.shrink}}, \code{\link{vcov.shrink}} }

%%\examples{ }

\keyword{ utilities }
