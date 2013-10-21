\name{predict.shrink}

\alias{predict.shrink}

\title{ Predict Method for Objects of Class \code{shrink} }

\description{ 
  Obtains predictions from shrunken regression coefficients from a 
  fit object of class \code{shrink}.
  This class of objects is returned by the \code{shrink} function. Objects of this 
  class have methods for the functions \code{print}, \code{predict}, and \code{coeff}. 
}

\usage{ 
\method{predict}{shrink}(object, newdata = NULL, type = c("link", "response", "lp", "risk",
"expected", "terms"), terms = NULL, na.action = na.pass, collapse, 
safe = FALSE, ...) 
}

\arguments{
  \item{object}{ an object of class \code{shrink}. }
  \item{newdata}{ a data frame for which predictions are obtained, otherwise 
                  predictions are based on the data used in the fit object. }
  \item{type}{ the type of prediction required. }
  \item{terms}{ with \code{type = "terms"} by default all terms are returned. A 
                character vector specifies which terms are to be returned. }
  \item{na.action}{	function determining what should be done with missing values 
                    in \code{newdata}. The default is to include all observations. }
  \item{collapse}{ if \code{family = coxph}, an optional vector of subject identifiers. 
                   If specified, the output will contain one entry per subject rather 
                   than one entry per observation. }
  \item{safe}{ option from \code{predict.mfp}. }
  \item{...}{ additional arguments to be passed to methods. }
}

%%\details{  }

\value{ A vector or matrix of predictions. 
  \item{fit}{ predictions. }
}

%%\references{ XXX }

\author{ Daniela Dunkler, Georg Heinze }

%%\note{ xxx }

\seealso{ \code{\link{shrink}} }

\examples{
library("mfp")
data("GBSG")
dat <- GBSG
set.seed(79)
dat <- dat[sample(seq_len(nrow(dat)), size = 400, replace = FALSE),]

fit1 <- mfp(Surv(rfst, cens) ~ fp(age, df = 4, select = 0.05) + 
            fp(prm, df = 4, select = 0.05), family = cox, data = dat)

dfbeta.global <- shrink(fit1, type = "global",  method = "dfbeta")
dfbeta.pw     <- shrink(fit1, type = "parameterwise", method = "dfbeta")
dfbeta.join   <- shrink(fit1, type = "parameterwise", method = "dfbeta", 
                        join=list(c("age.1", "age.2")))

# unshrunken
plot(20:80, predict(fit1, newdata = data.frame(age = 20:80, prm = 0), type = "lp") -
     predict(fit1, newdata = data.frame(age = 50, prm = 0), type = "lp"), xlab = "Age", 
     ylab = "Log hazard relative to 50 years", type = "l", lwd = 2)

# globally shrunken
lines(20:80, predict(dfbeta.global,newdata = data.frame(age = 20:80, prm = 0), type = "lp") - 
      predict(dfbeta.global, newdata = data.frame(age = 50, prm = 0), type = "lp"), 
      lty = 4, col = "blue", lwd = 2)

# jointly shrunken
lines(20:80, predict(dfbeta.join, newdata = data.frame(age = 20:80, prm = 0), type = "lp") -
      predict(dfbeta.join, newdata = data.frame(age = 50, prm = 0), type = "lp"), 
      lty = 3, col = "red", lwd = 2)

# parameterwise shrunken
lines(20:80, predict(dfbeta.pw, newdata = data.frame(age = 20:80, prm = 0), type = "lp") -
      predict(dfbeta.pw, newdata =data.frame(age = 50, prm = 0), type = "lp"), 
      lty = 2, col = "green", lwd = 2)

legend("topright", lty = c(1, 4, 3, 2), legend = c("No", "Global", "Joint", "Parameterwise"), 
       title = "SHRINKAGE", inset = 0.01, bty = "n", col = c("black", "blue", "red", "green"), 
       lwd = 2)  
}

\keyword{ utilities }
