\name{deepvein}

\alias{deepvein}

\title{ Deep Vein Thrombosis Study }

\usage{data(deepvein)}

\description{
  A data frame containing time to recurrence of thrombosis and several potential
  prognostic factors measured at baseline for 929 individuals with deep vein
  thrombosis or unprovoked pulmonary embolism. 147 events of recurrence were 
  observed during a median follow-up time of 37.8 months.
    
  The data are a modified and partly simulated version of the data set used by 
  Eichinger et al. (2010) and are available under a GPL-2 license.
}

\format{
  The data frame contains observations of 929 individuals and the following 
  variables:
  \describe{
    \item{pnr}{ patient number. }
    \item{time}{ time to recurrence of thrombosis or end of study in months. }
    \item{status}{ = 1 recurrence of thrombosis. }
    \item{sex}{ gender. }
    \item{fiimut}{ factor II G20210A mutation. }
    \item{fvleid}{ factor V Leiden mutation. }
    \item{log2ddim}{ log2-transformed D-dimer. }
    \item{bmi}{ body mass index. }
    \item{durther}{ duration of anticoagulation therapy. }
    \item{age}{ age in years. }
    \item{loc}{ location of first thrombosis: pulmonary embolism (PE), distal, 
                or proximal deep vein \cr thrombosis. }
  }
}

%%\source{  e.g. \url{http://lib.stat.cmu.edu/datasets/bodyfat}}

\references{
  Eichinger S, Heinze G, Jandeck LM, Kyrle P (2010). Risk assessment of recurrence 
  in patients with unprovoked deep vein thrombosis or pulmonary embolism - The 
  Vienna prediction model. \emph{Circulation}, \bold{121}(14): 1630-1636.
}

\examples{
data("deepvein")
summary(deepvein)
 
library("survival")
deepvein$loc <- relevel(deepvein$loc, ref = "PE")
fitfull <- coxph(Surv(time, status) ~ sex + fiimut + fvleid + log2ddim + 
                 bmi + durther + age + loc, data = deepvein, x = TRUE)
summary(fitfull)

# DFBETA approximation
shrink(fitfull, type = "global", method = "dfbeta")
shrink(fitfull, type = "parameterwise", method = "dfbeta")
shrink(fitfull, type = "parameterwise", method = "dfbeta", 
       join = list(c("locdistal", "locproximal")))

# jackknife estimate (time consuming)
#shrink(fitfull, type = "global", method = "jackknife")
#shrink(fitfull, type = "parameterwise", method = "jackknife")
#shrink(fitfull, type = "parameterwise", method = "jackknife", 
#       join = list(c("locdistal", "locproximal")))

# backward elimination based on AIC
fitselect <- step(fitfull, direction = "backward")
summary(fitselect)

shrink(fitselect, type = "global", method = "dfbeta")
shrink(fitselect, type = "parameterwise", method = "dfbeta")
shrink(fitselect, type = "parameterwise", method = "dfbeta", 
       join = list(c("locdistal", "locproximal")))

#shrink(fitselect, type = "global", method = "jackknife")
#shrink(fitselect, type = "parameterwise", method = "jackknife")
#shrink(fitselect, type = "parameterwise", method = "jackknife", 
#       join = list(c("locdistal", "locproximal")))
}

\keyword{datasets}
