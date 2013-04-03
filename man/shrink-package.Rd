\name{shrink-package}

\alias{shrink-package}

\docType{package}

\title{ Global, Parameterwise, and Joint Shrinkage of Regression Coefficients }

\description{ In prognostic modeling shrinkage of regression coefficients is needed to correct the 
 optimism in such a model caused by variable selection. While global shrinkage modifies all regression 
 coefficients by the same factor, parameterwise shrinkage factors differ between regression coefficients,
 their values depending on each variable's prognostic relevance. With highly correlated or semantically related 
 variables, such as several columns of a design matrix describing a nonlinear effect, a compromise between 
 global, and parameterwise shrinkage, termed 'joint shrinkage', is a useful extension. 
 A computational shortcut to resampling-based shrinkage factor estimation based on DFBETA residuals is 
 applied. 
 Global, parameterwise, and joint shrinkage for models fitted by \code{lm}, \code{glm}, \code{coxph}, and 
 \code{mfp} is available. }

\details{
\tabular{ll}{
Package: \tab shrink\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-04-03\cr
License: \tab GPL\cr
}
Functions contained in the shrink-package are
\tabular{ll}{
\code{ shrink } \tab a function to compute global, parameterwise, and joint shrinkage factors of \cr
\code{ } \tab regression coefficients of fit objects of class \code{lm}, \code{glm}, \code{coxph}, or \code{mfp}. \cr
\code{ predict.shrink } \tab obtains shrunken regression coefficients from a fit object of class \code{shrink}. \cr
\code{ print.shrink } \tab prints objects of class \code{shrink}. \cr
\code{ coef.shrink } \tab returns shrunken regression coefficients from objects of class \code{shrink}. \cr
}
}

\author{ Daniela Dunkler, Georg Heinze \cr{Maintainer: \email{daniela.dunkler@meduniwien.ac.at}} }

\references{ Sauerbrei W (1999) The use of resampling methods to simplify regression models in medial statistics. \emph{Applied Statistics} 48(3): 313-329. \cr
Verweij P, van Houwelingen J (1993) Cross-validation in survival analysis. \emph{Statistics in Medicine} 12(24): 2305-2314. }

\keyword{ package }
\keyword{ models }
\keyword{ regression }
\keyword{ survival }
\keyword{ nonlinear }

\seealso{ \code{\link{shrink}}, \code{\link{print.shrink}}, \code{\link{predict.shrink}}, \code{\link{coef.shrink}} }

\examples{
# with glm, family=binomial
set.seed(888)
intercept <- 1
beta <- c(0.5, 1.2)
n <- 1000 
x1 <- rnorm(n,1,1)
x2 <- rbinom(n, 1, 0.3)
linpred <- intercept + x1*beta[1] + x2*beta[2]
prob <- exp(linpred)/(1 + exp(linpred))
runis <- runif(n,0,1)
ytest <- ifelse(runis < prob,1,0)
simdat <- data.frame(cbind(y=ifelse(runis < prob, 1, 0), x1, x2))

fit2 <- glm(y ~ x1 + x2, family = binomial, data = simdat, x = TRUE)
summary(fit2)

dfbeta.pw <- shrink(fit2, type = "parameterwise", method = "dfbeta")
print(dfbeta.pw)
coef(dfbeta.pw)

#shrink(fit2, type = "parameterwise", method = "jackknife")

shrink(fit2, type = "parameterwise", method = "dfbeta", join = list(c("x1", "x2")))

#shrink(fit2, type = "parameterwise", method = "jackknife", join = list(c("x1", "x2")))

shrink(fit2, type = "global", method = "dfbeta")

#shrink(fit2, type = "global", method = "jackknife")

# For more examples see shrink
}
