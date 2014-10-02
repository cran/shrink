\name{vcov.shrink}

\alias{vcov.shrink}

\title{ Calculate Variance-Covariance Matrix of Shrinkage Factors for Objects of Class \code{shrink} }

\description{ 
  This class of objects is returned by the \code{shrink} function. Objects of this 
  class have methods for the functions \code{coef}, \code{predict}, \code{print}, 
  \code{summary}, and \code{vcov}.
}

\usage{
\method{vcov}{shrink}(object, digits = 6, ...) 
}

\arguments{
  \item{object}{ object of class \code{shrink}. }
  \item{digits}{ integer, used for number formatting with \code{\link[base]{signif}}(). }
  \item{...}{ further arguments. }  
}

%%\details{ }

\value{ A matrix of the estimated covariances between the obtained shrinkage factors. }

%%\references{  }

\author{ Daniela Dunkler, Georg Heinze }

%%\note{  }

\seealso{ \code{\link{shrink}},\code{\link{coef.shrink}}, \code{\link{predict.shrink}}, 
          \code{\link{print.shrink}}, \code{\link{summary.shrink}} }

%%\examples{ }

\keyword{ utilities }
