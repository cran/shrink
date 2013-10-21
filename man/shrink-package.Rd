\name{shrink-package}

\alias{shrink-package}

\docType{package}

\title{ Global, Parameterwise, and Joint Post-Estimation Shrinkage }

\description{ In statistical modeling post-estimation shrinkage can be used to 
 correct for the overestimation of regression coefficients caused by variable 
 selection. While global shrinkage modifies all regression coefficients by the 
 same factor, parameterwise shrinkage factors differ between regression coefficients.
 With highly correlated or semantically related variables, such as several columns
 of a design matrix describing a nonlinear effect, parameterwise shrinkage factors
 are not interpretable and a compromise between global and parameterwise shrinkage, 
 termed 'joint shrinkage', is a useful extension. 
 A computational shortcut to resampling-based shrinkage factor estimation based 
 on DFBETA residuals is applied. 
 Global, parameterwise, and joint shrinkage for models fitted by \code{lm}, 
 \code{glm}, \code{coxph}, and \code{mfp} is available. 
}

\details{
\tabular{ll}{
Package: \tab shrink\cr
Type: \tab Package\cr
Version: \tab 1.1.\cr
Date: \tab 2013-10-21\cr
License: \tab GPL-2\cr
}
Functions contained in the \code{shrink}-package are
\tabular{ll}{
\code{ shrink } \tab a function to compute global, parameterwise, and joint post-estimation \cr 
\code{ } \tab shrinkage factors of fit objects of class \code{lm}, \code{glm}, \code{coxph}, or \code{mfp}. \cr
\code{ } \tab   \cr
\code{ coef.shrink } \tab returns shrunken regression coefficients from objects of class \code{shrink}. \cr
\code{ predict.shrink } \tab obtains predictions from shrunken regression coefficients from a fit object \cr
\code{ } \tab  of class \code{shrink}. \cr
\code{ print.shrink } \tab prints objects of class \code{shrink}. \cr
}

Data set contained in the \code{shrink}-package
\tabular{ll}{
\code{ deepvein } \tab deep vein thrombosis study \cr 
}
}

\author{ Daniela Dunkler, Georg Heinze \cr{Maintainer: \email{daniela.dunkler@meduniwien.ac.at}} }

\references{ Sauerbrei W (1999) The use of resampling methods to simplify regression models in medial statistics. \emph{Applied Statistics} \bold{48}(3): 313-329. \cr
Verweij P, van Houwelingen J (1993) Cross-validation in survival analysis. \emph{Statistics in Medicine} \bold{12}(24): 2305-2314. }

\keyword{ package }
\keyword{ models }
\keyword{ regression }
\keyword{ survival }
\keyword{ nonlinear }

\seealso{ \code{\link{shrink}}, \code{\link{coef.shrink}}, \code{\link{predict.shrink}}, \code{\link{print.shrink}} }

\examples{
# with glm, family = binomial
set.seed(888)
intercept <- 1
beta <- c(0.5, 1.2)
n <- 200 
x1 <- rnorm(n, mean = 1, sd = 1)
x2 <- rbinom(n, size = 1, prob = 0.3)
linpred <- intercept + x1*beta[1] + x2*beta[2]
prob <- exp(linpred)/(1 + exp(linpred))
runis <- runif(n,0,1)
ytest <- ifelse(runis < prob,1,0)
simdat <- data.frame(cbind(y = ifelse(runis < prob, 1, 0), x1, x2))

fit <- glm(y ~ x1 + x2, family = binomial, data = simdat, x = TRUE)
summary(fit)

global <- shrink(fit, type = "global", method = "dfbeta")
print(global)
coef(global)

shrink(fit, type = "parameterwise", method = "dfbeta")

shrink(fit, type = "parameterwise", method = "dfbeta", join = list(c("x1", "x2")))

#shrink(fit, type = "global", method = "jackknife")
#shrink(fit, type = "parameterwise", method = "jackknife")
#shrink(fit, type = "parameterwise", method = "jackknife", 
#       join = list(c("x1", "x2")))

# For more examples see shrink
}
