\name{shrink}

\alias{shrink}

\title{ Global, Parameterwise and Joint Shrinkage of Regression Coefficients }

\description{ 
  Obtain global, parameterwise and joint post-estimation shrinkage factors for 
  regression coefficients from fit objects of class \code{lm}, \code{glm}, 
  \code{coxph}, or (for R<=3.1.0) \code{mfp}. 
}

\usage{ shrink(fit, type = c("parameterwise", "global", "all"), method = c("jackknife", 
   "dfbeta"), join = NULL, notes = TRUE, postfit = TRUE) }

\arguments{
  \item{fit}{ a fit object of class \code{lm}, \code{glm}, \code{coxph}, or \code{mfp}. 
              The fit object must have been called with \code{x = TRUE} (and 
              \code{y = TRUE} in case of class \code{lm}). }
  \item{type}{ of shrinkage, either \code{"parameterwise"} (default), \code{"global"} shrinkage,
               or \code{"all"}. }
  \item{method}{ of shrinkage estimation, either \code{"jackknife"} (based on leave-one-out 
                 resampling, default) or \code{"dfbeta"} (excellent approximation based on 
                 DFBETA residuals). }
  \item{join}{ compute optional joint shrinkage factors for sets of specified columns 
               of the design matrix, if \code{type =} \code{"parameterwise"}. See details. }
  \item{notes}{ print notes. Default is TRUE. }
  \item{postfit}{ obtain fit with shrunken regression coefficients. This option is only available 
                  for models without an intercept. Default is TRUE. }
}

\details{ While global shrinkage modifies all regression coefficients by the same 
factor, parameterwise shrinkage factors differ between regression coefficients. 
With variables which are either highly correlated or associated with regard to 
contents, such as several columns of a design matrix describing a nonlinear 
effect, parameterwise shrinkage factors are not interpretable. Joint shrinkage 
of a set of such associated design variables will give one common shrinkage 
factor for this set.

Joint shrinkage factors may be useful when analysing highly correlated and/or 
such associated columns of the design matrix, e.g. dummy variables corresponding 
to a categorical explanatory variable with more than two levels, two variables 
and their pairwise interaction term, or several transformations of an explantory 
variable enabling estimation of nonlinear effects. The analyst can define 
'joint' shrinkage factors by specifing the \code{join} option if 
\code{type = "parameterwise"}. \code{join} expects a list with at least one character 
vector including the names of the columns of the design matrix for which a joint 
shrinkage factor is requested. For example the following specification of 
\code{join = list(c("dummy1"}, \code{"dummy2"}, \code{"dummy3"}), \code{c("main1"}, 
\code{"main2"}, \code{"interaction"}), \code{c("varX.fp1"}, \code{"varX.fp2"})) requests 
the joint shrinkage factors for a) \code{"dummy1"}, \code{"dummy2"} and \code{"dummy3"}, 
b) \code{"main1"}, \code{"main2"} and \code{"interaction"} and c) \code{"varX.fp1"} and 
\code{"varX.fp2"}. 

\code{shrink} also works for models incorporating restricted cubic splines 
computed with the \code{\link[rms]{rcs}} function from the \code{rms} library. A joint 
shrinkage factor of explanatory variable \code{varX} transformed with \code{\link[rms]{rcs}} 
can be obtained by \code{join =} \code{list(c("rcs(varX)"))} or by stating the names of the 
\code{\link[rms]{rcs}}-transformed variables as given in the respective fit object. (These 
two notations should not be mixed within one call to \code{shrink}.)

For linear regression models (\code{lm} or \code{glm} with \code{family = "gaussian"}) 
shrinkage factors obtained by Jackknife and the DFBETA approximation will be identical. For 
all other types of regression, the computational effort of estimating shrinkage factors may 
be greatly reduced by using \code{method =} \code{"dfbeta"} instead. However, for (very) 
small data sets \code{method = "jackknife"} may be of advantage, as the use of DFBETA 
residuals may underestimate the influence of some highly influential observations.

A shrunken intercept is estimated as follows: For all columns of the design matrix 
except for the intercept the shrinkage factors are multiplied with the respective 
regression coefficients and a linear predictor is computed. Then the shrunken 
intercept is estimated by modeling \code{fit$y ~} \code{offset(linear predictor)}. 

For regression models without an intercept, i.e., fit objects of class \code{coxph}, 
the shrunken regression coefficients can be directly estimated. This postfit is retained 
in the \code{$postfit} component of the \code{shrink} object.
}

\value{ 
\code{shrink} returns an object with the following components:
  \item{ShrinkageFactors}{ a vector of shrinkage factors of regression coefficients. }
  \item{ShrinkageFactorsVCOV}{ the covariance matrix of the shrinkage factors. }
  \item{ShrunkenRegCoef}{ a vector with the shrunken regression coefficients. }
  \item{postfit}{ an optional postfit model with shrunken regression coefficients and associated \cr
                  standard errors for models without an intercept. }
  \item{fit}{ the original (unshrunken) fit object. }
  \item{type}{ the requested shrinkage \code{type}. }
  \item{method}{ the requested shrinkage \code{method}. }  
  \item{join}{ the requested joint shrinkage factors. }
  \item{call}{ the function call. }
  
If \code{type = "all"} then the object returned by \code{shrink} additionally contains
  \item{global}{ a list with the following elements: \code{ShrinkageFactors}, 
                 \code{ShrinkageFactorsVCOV} and \code{ShrunkenRegCoef}. }
  \item{parameterwise}{ a list with the following elements: \code{ShrinkageFactors}, 
                 \code{ShrinkageFactorsVCOV} and \code{ShrunkenRegCoef}. }
  \item{joint}{ an optional list with the following elements: \code{ShrinkageFactors}, 
                 \code{ShrinkageFactorsVCOV} and \code{ShrunkenRegCoef}. }
}

\note{ For fit objects of class \code{mfp} with \code{family != cox} regression coefficients 
of \code{fit} (obtained by \code{coef(fit)}) and \code{fit$fit} may not always be identical, 
because of \code{mfp}'s pretransformation applied to the explanatory variables in the model. 
The \code{shrink} function uses a) the names as given in \code{names(coef(fit))} and b) the 
regression coefficients as given in \code{summary(fit)} which correspond to the pretransformed 
explanatory variables. }
%For models of \code{family = binomial} the explanatory variable must be numeric in order to compute DFBETA residuals. 

\references{ Sauerbrei W (1999) The use of resampling methods to simplify regression 
models in medial statistics. \emph{Applied Statistics} \bold{48}(3): 313-329. \cr
Verweij P, van Houwelingen J (1993) Cross-validation in survival analysis. 
\emph{Statistics in Medicine} \bold{12}(24): 2305-2314. }

\author{ Daniela Dunkler, Georg Heinze }

\seealso{ \code{\link{coef.shrink}}, \code{\link{predict.shrink}}, \code{\link{print.shrink}}, 
\code{\link{summary.shrink}}, \code{\link{vcov.shrink}} }

\examples{
## Example with Cox proportional hazards regression & rcs
library("survival")
library("rms")
data("GBSG")

fit1 <- coxph(Surv(rfst, cens) ~ rcs(age) + log(prm + 1), data = GBSG, x = TRUE)

shrink(fit1, type = "global", method = "dfbeta")
shrink(fit1, type = "parameterwise", method = "dfbeta")
shrink(fit1, type = "parameterwise", method = "dfbeta", 
       join = list(c("rcs(age)"))) 
shrink(fit1, type = "parameterwise", method = "dfbeta", 
       join = list(c("rcs(age)"), c("log(prm + 1)"))) 

dfbeta.pw <- shrink(fit1, type = "parameterwise", method = "dfbeta")
dfbeta.pw
dfbeta.pw$postfit

# correlations between shrinkage factors and standard errors of shrinkage factors
cov2cor(dfbeta.pw$ShrinkageFactorsVCOV)  
sqrt(diag(dfbeta.pw$ShrinkageFactorsVCOV))


## Examples with glm (family = binomial)
set.seed(888)
intercept <- 1
beta <- c(0.5, 1.2)
n <- 1000 
x1 <- rnorm(n, mean = 1, sd = 1)
x2 <- rbinom(n, size = 1, prob = 0.3)
linpred <- intercept + x1 * beta[1] + x2 * beta[2]
prob <- exp(linpred) / (1 + exp(linpred))
runis <- runif(n, 0, 1)
ytest <- ifelse(test = runis < prob, yes = 1, no = 0)
simdat <- data.frame(cbind(y = ifelse(runis < prob, 1, 0), x1, x2))

fit2 <- glm(y ~ x1 + x2, family = binomial, data = simdat, x = TRUE)
summary(fit2)

shrink(fit2, type = "global", method = "dfbeta")
shrink(fit2, type = "parameterwise", method = "dfbeta")
shrink(fit2, type = "parameterwise", method = "dfbeta", join = list(c("x1", "x2")))


utils::data("Pima.te", package="MASS")      
utils::data("Pima.tr", package="MASS")      
Pima <- rbind(Pima.te, Pima.tr)
fit3 <- glm(type ~ npreg + glu + bmi + ped + rcs(age), 
            family = binomial, data = Pima, x = TRUE) 
summary(fit3)                                 

shrink(fit3, type = "global", method = "dfbeta")
shrink(fit3, type = "parameterwise", method = "dfbeta")
shrink(fit3, type = "parameterwise", method = "dfbeta", join = list(c("rcs(age)")))


## Examples with glm (family = gaussian) and lm 
utils::data("anorexia", package = "MASS")
contrasts(anorexia$Treat) <- contr.treatment(n = 3, base = 2)
fit4 <- glm(Postwt ~ Prewt + Treat, family = gaussian, data = anorexia, x = TRUE)
fit4

shrink(fit4, type = "global", method = "dfbeta")
# which is identical to the more time-consuming jackknife approach:
# shrink(fit4, type = "global", method = "jackknife")

shrink(fit4, type = "parameterwise", method = "dfbeta")
shrink(fit4, type = "parameterwise", method = "dfbeta", 
       join = list(c("Treat1", "Treat3")))


fit5 <- lm(Postwt ~ Prewt + Treat, data = anorexia, x = TRUE, y = TRUE)
fit5

shrink(fit5, type = "global", method = "dfbeta")
shrink(fit5, type = "parameterwise", method = "dfbeta")
shrink(fit5, type = "parameterwise", method = "dfbeta", 
       join=list(c("Treat1", "Treat3")))


utils::data("GAGurine", package = "MASS")
fit6 <- glm(Age ~ rcs(GAG), family = gaussian, data = GAGurine, x = TRUE)
summary(fit6)   

shrink(fit6, type = "global", method = "dfbeta")
shrink(fit6, type = "parameterwise", method = "dfbeta")
shrink(fit6, type = "parameterwise", method = "dfbeta", join = list("rcs(GAG)"))

## For R up to 3.1.0: Example with mfp (family = cox)
#library("mfp")
#fit7 <- mfp(Surv(rfst, cens) ~ fp(age, df = 4, select = 0.05) + 
#              fp(prm, df = 4, select = 0.05), family = cox, data = GBSG)

#shrink(fit7, type = "global", method = "dfbeta")
#shrink(fit7, type = "parameterwise", method = "dfbeta", 
#       join = list(c("age.1", "age.2"))) 

#shrink(fit7, type = "global", method = "jackknife")
#shrink(fit7, type = "parameterwise", method = "jackknife")
#shrink(fit7, type = "parameterwise", method = "jackknife", 
#       join = list(c("age.1", "age.2"))) 

# obtain global, parameterwise and joint shrinkage with one call to 'shrink'
#shrink(fit7, type = "all", method = "dfbeta", 
#       join = list(c("age.1", "age.2")))
}

\keyword{ models }
\keyword{ regression }
\keyword{ survival }
\keyword{ nonlinear }
