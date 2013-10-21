\name{shrink}

\alias{shrink}

\title{ Global, Parameterwise, and Joint Shrinkage of Regression Coefficients }

\description{ 
  Obtain global, parameterwise, and joint post-estimation shrinkage factors for 
  regression coefficients from fit objects of class \code{lm}, class \code{glm} with 
  \code{family = c("gaussian", "binomial")}, class \code{coxph}, or class \code{mfp} 
  with \code{family = c(cox, gaussian, binomial)}. 
}

\usage{ shrink(fit, type = "parameterwise", method = "jackknife", join = NULL) }

\arguments{
  \item{fit}{ a fit object of class \code{lm}, \code{glm}, \code{coxph} or \code{mfp}. 
              The fit object must have been called with \code{x = TRUE} (and 
              \code{y = TRUE} in case of \code{lm}). }
  \item{type}{ of shrinkage, either \code{"parameterwise"} (or \code{"p"}; default) 
               or \code{"global"} (or \code{"g"}) shrinkage. }
  \item{method}{ of shrinkage estimation, either \code{"jackknife"} (or \code{"j"}; 
                 default, based on leave-one-out resampling) or \code{"dfbeta"} 
                 (or \code{"d"}; excellent approximation based on DFBETA residuals). }
  \item{join}{ compute optional joint shrinkage factors for sets of specified columns 
               of the design matrix, if \code{type =} \code{"parameterwise"}. See details. }
}

\details{ While global shrinkage modifies all regression coefficients by the same 
factor, parameterwise shrinkage factors differ between regression coefficients. 
With highly correlated or semantically related variables, such as several columns 
of a design matrix describing a nonlinear effect, parameterwise shrinkage factors 
are not interpretable. Joint shrinkage of a set of columns of the design matrix 
will give one common shrinkage factor for this set.

Joint shrinkage factors may be useful when analysing highly correlated and/or 
semantically related columns of the design matrix, such as dummy variables 
corresponding to a categorical explanatory variable with more than two levels, 
two variables and their pairwise interaction term, or several transformations of 
an explantory variable enabling estimation of nonlinear effects. The analyst can 
define such 'joint' shrinkage factors by specifing the \code{join} option if 
\code{type = "parameterwise"}. \code{join} expects a list with at least one character 
vector including the names of the columns of the design matrix for which a joint 
shrinkage factor is requested; e. g. \code{join = list(c("dummy1"}, \code{"dummy2"}, 
\code{"dummy3"}), \code{c("main1"}, \code{"main2"}, \code{"interaction"}), 
\code{c("varX.fp1"}, \code{"varX.fp2"})). 

\code{shrink} also works for models incorporating restricted cubic splines 
computed with the \code{rcs} function from the \code{rms} library. A joint shrinkage
factor of variable \code{varX} transformed with \code{rcs} can be obtained by 
\code{join =} \code{list(c("varX"))}.

For fit objects of class \code{coxph} or \code{glm} with \code{family = "binomial"} 
the computational effort of estimating shrinkage factors may be greatly reduced 
by using \code{method =} \code{"dfbeta"} instead. However, for (very) small data sets 
\code{method = "jackknife"} may be of advantage, as the use of DFBETA residuals may 
underestimate the influence of some highly influential observations.

A shrunken intercept is estimated as follows: For all columns of the design matrix 
except for the intercept the shrinkage factors are multiplied with the respective 
regression coefficients and a linear predictor is computed. Then the shrunken 
intercept is estimated by modeling \code{fit$y ~} \code{offset(linear predictor)}. 

For regression models without an intercept, i.e. fit objects of class \code{coxph}, 
with \code{type = "parameterwise"}, the shrunken regression coefficients can be 
directly estimated. This postfit is retained in the \code{$postfit} slot of the 
\code{shrink} object.
}

\value{ 
\code{shrink} returns an object with the following components:
  \item{shrinkage}{ a vector of shrinkage factors of regression coefficients. }
  \item{vcov.shrinkage}{ a covariance matrix of shrinkage factors. }
  \item{shrunken}{ a vector with the shrunken regression coefficients. }
  \item{postfit}{ an optional postfit model with shrunken regression coefficients and associated \cr
                  standard errors if \code{type = "parameterwise"} and \code{join = NULL}. }
  \item{fit}{ the original (unshrunken) fit object. }
  \item{type}{ the requested shrinkage \code{type}. }
  \item{method}{ the requested shrinkage \code{method}. }  
  \item{call}{ the function call. }
}

\note{ For fit objects of class \code{mfp} with \code{family = binomial} or 
\code{gaussian} the regression coefficients of \code{fit} (obtained by \code{coef(fit)}) 
and \code{fit$fit} (\code{coef(fit$fit)}) may not always be identical, because of 
\code{mfp}'s pretransformation applied to the explanatory variables in the model. 
The \code{shrink} function uses the regression coefficients from \code{fit$fit} 
which correspond to the pretransformed explanatory variables. }
%For models of \code{family = binomial} the explanatory variable must be numeric in order to compute DFBETA residuals. 

\references{ Sauerbrei W (1999) The use of resampling methods to simplify regression 
models in medial statistics. \emph{Applied Statistics} \bold{48}(3): 313-329. \cr
Verweij P, van Houwelingen J (1993) Cross-validation in survival analysis. 
\emph{Statistics in Medicine} \bold{12}(24): 2305-2314. }

\author{ Daniela Dunkler, Georg Heinze }

\seealso{ \code{\link{coef.shrink}}, \code{\link{predict.shrink}}, \code{\link{print.shrink}} }

\examples{
# Example with mfp (family = cox)
library("mfp")
data("GBSG")
fit1 <- mfp(Surv(rfst, cens) ~ fp(age, df = 4, select = 0.05) + 
            fp(prm, df = 4, select = 0.05), family = cox, data = GBSG)

shrink(fit1, type = "global", method = "dfbeta")

dfbeta.pw <- shrink(fit1, type = "parameterwise", method = "dfbeta")
dfbeta.pw
cov2cor(dfbeta.pw$vcov.shrinkage)  
sqrt(diag(dfbeta.pw$vcov.shrinkage))

shrink(fit1, type = "parameterwise", method = "dfbeta", 
       join = list(c("age.1", "age.2"))) 

#shrink(fit1, type = "global", method = "jackknife")
#shrink(fit1, type = "parameterwise", method = "jackknife")
#shrink(fit1, type = "parameterwise", method = "jackknife", 
#       join = list(c("age.1", "age.2"))) 


# Example with rcs
library("rms")
fit2 <- coxph(Surv(rfst, cens) ~ rcs(age) + rcs(prm), data = GBSG, x = TRUE)

shrink(fit2, type = "global", method = "dfbeta")
shrink(fit2, type = "parameterwise", method = "dfbeta")
shrink(fit2, type = "parameterwise", method = "dfbeta", 
       join = list(c("age"), c("prm"))) 
  

# Examples with glm & mfp (family=binomial)
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
simdat <- data.frame(cbind(y = ifelse(runis < prob, 1, 0), x1, x2))

fit2 <- glm(y ~ x1 + x2, family = binomial, data = simdat, x = TRUE)
summary(fit2)

shrink(fit2, type = "global", method = "dfbeta")
shrink(fit2, type = "parameterwise", method = "dfbeta")
shrink(fit2, type = "parameterwise", method = "dfbeta", join = list(c("x1", "x2")))


utils::data("Pima.te", package="MASS")      
utils::data("Pima.tr", package="MASS")      
Pima <- rbind(Pima.te, Pima.tr)
Pima$type2 <- as.numeric(Pima$type)-1
fit3 <- mfp(type2 ~ npreg + glu + bmi + ped + fp(age, select = 0.05), 
            family = binomial, data = Pima) 
fit3

shrink(fit3, type = "global", method = "dfbeta")
shrink(fit3, type = "parameterwise", method="dfbeta")


# Examples with glm & mfp (family = gaussian) and lm 
utils::data("anorexia", package = "MASS")
contrasts(anorexia$Treat) <- contr.treatment(3, base = 2)
fit4 <- glm(Postwt ~ Prewt + Treat, family = gaussian, data = anorexia, x = TRUE)
fit4

shrink(fit4, type = "global", method = "dfbeta")
shrink(fit4, type = "parameterwise", method = "dfbeta")
shrink(fit4, type = "parameterwise", method = "dfbeta", 
       join = list(c("Treat1", "Treat3")))


fit5 <- lm(Postwt ~ Prewt + Treat, data =anorexia, x = TRUE, y = TRUE)
fit5

shrink(fit5, type = "global", method = "dfbeta")
shrink(fit5, type = "parameterwise", method = "dfbeta")
shrink(fit5, type = "parameterwise", method = "dfbeta", 
       join=list(c("Treat1", "Treat3")))

       
utils::data("GAGurine", package="MASS")
fit6 <- mfp(Age ~ fp(GAG, select = 0.05), family = gaussian, data = GAGurine)
fit6

shrink(fit6, type = "global", method = "dfbeta")
shrink(fit6, type = "parameterwise", method = "dfbeta")
shrink(fit6, type = "parameterwise", method = "dfbeta", 
       join = list(c("GAG.1", "GAG.2")))
}

\keyword{ models }
\keyword{ regression }
\keyword{ survival }
\keyword{ nonlinear }
