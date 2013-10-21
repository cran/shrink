### tests
library(shrink)
library(rms)
library(mfp)
data(GBSG)

# cox (univariate, multivariate)
fit <- coxph(Surv(rfst, cens) ~ age, data = GBSG, x=TRUE)
a<-shrink(fit, method = "dfbeta", type = "parameterwise")
shrink(fit, method = "dfbeta", type = "parameterwise", join=list(c("age")))
c<-shrink(fit, method = "jackknife", type = "parameterwise")
d<-shrink(fit, method = "dfbeta", type = "global")
e<-shrink(fit, method = "jackknife", type = "global")

a
c
d
e

names(a)
a$shrinkage
a$vcov.shrinkage
a$shrunken
a$postfit
a$fit
a$type
a$method
a$call

names(c)
c$shrinkage
c$vcov.shrinkage
c$shrunken
c$postfit
c$fit
c$type
c$method
c$call

names(d)
d$shrinkage
d$vcov.shrinkage
d$shrunken
d$postfit
d$fit
d$type
d$method
d$call

names(e)
e$shrinkage
e$vcov.shrinkage
e$shrunken
e$postfit
e$fit
e$type
e$method
e$call

fit <- coxph(Surv(rfst, cens) ~ age + prm, data = GBSG, x=TRUE)
shrink(fit, method = "dfbeta", type = "parameterwise")
shrink(fit, method = "dfbeta", type = "parameterwise", join=list(c("age", "prm")))
shrink(fit, method = "jackknife", type = "parameterwise")
shrink(fit, method = "dfbeta", type = "global")
shrink(fit, method = "jackknife", type = "global")

# cox (univariate, multivariate, rcs)
fit <- coxph(Surv(rfst, cens) ~ rcs(age), data = GBSG, x=TRUE)
shrink(fit, method = "dfbeta", type = "parameterwise")
shrink(fit, method = "dfbeta", type = "parameterwise", join=list(c("rcs1.age", "rcs2.age", "rcs3.age", "rcs4.age")))
shrink(fit, method = "dfbeta", type = "parameterwise", join=list(c("age")))
shrink(fit, method = "jackknife", type = "parameterwise")
shrink(fit, method = "dfbeta", type = "global")
shrink(fit, method = "jackknife", type = "global")

fit <- coxph(Surv(rfst, cens) ~ rcs(age) + rcs(prm), data = GBSG, x=TRUE)
shrink(fit, method = "dfbeta", type = "parameterwise")
shrink(fit, method = "dfbeta", type = "parameterwise", join=list(c("age"), c("prm")))
shrink(fit, method = "dfbeta", type = "parameterwise", join=list(c("age", "prm")))
shrink(fit, method = "jackknife", type = "parameterwise")
shrink(fit, method = "dfbeta", type = "global")
shrink(fit, method = "jackknife", type = "global")

fit <- coxph(Surv(rfst, cens) ~ age + rcs(prm), data = GBSG, x=TRUE)
shrink(fit, method = "dfbeta", type = "parameterwise")
shrink(fit, method = "dfbeta", type = "parameterwise", join=list(c("prm")))
shrink(fit, method = "jackknife", type = "parameterwise")
shrink(fit, method = "dfbeta", type = "global")
shrink(fit, method = "jackknife", type = "global")

# cox and mfp
library(mfp)
data(GBSG)
fit <- mfp(Surv(rfst, cens) ~ fp(age, df = 4, select = 0.05) +  fp(prm, df = 4, select = 0.05), family = cox, data = GBSG)    # + fp(esm, df = 4, select = 1)
shrink(fit, method = "dfbeta", type = "parameterwise")
shrink(fit, method = "dfbeta", type = "global")
shrink(fit, method = "dfbeta", type = "parameterwise", join=list(c("age.1", "age.2", "prm.1")))

# glm, binomial
utils::data(Pima.te, package="MASS")      
utils::data(Pima.tr, package="MASS")      
Pima <- rbind(Pima.te, Pima.tr)
nrow(Pima)
fit <- glm(type ~ npreg + glu + rcs(bmi) + ped + age, family = binomial, data = Pima, x = TRUE) 
shrink(fit, method="dfbeta", type="global")
shrink(fit, method="dfbeta", type="parameterwise")
shrink(fit, method="dfbeta", type="parameterwise", join=list(c("npreg", "glu"), c("bmi")))        
shrink(fit, method="dfbeta", type="parameterwise", join=list(c("npreg", "glu", "bmi")))           # does not work
shrink(fit, method="dfbeta", type="parameterwise", join=list(c("npreg", "glu", "rcs1.bmi" ,   "rcs2.bmi",    "rcs3.bmi"  ,  "rcs4.bmi")))           # use instead
shrink(fit, method="dfbeta", type="parameterwise", join=list(c("npreg", "glu", "ped","age")))


################################################################################
# further examples
fit1a <- coxph(Surv(rfst, cens) ~ age + prm, data = GBSG, x=TRUE)
shrink(fit1a, method = "dfbeta", type = "parameterwise")
shrink(fit1a, method = "dfbeta", type = "parameterwise", join=list(c("age", "prm")))
shrink(fit1a, method = "jackknife", type = "parameterwise")
shrink(fit1a, method = "dfbeta", type = "global")
shrink(fit1a, method = "jackknife", type = "global")

library(mfp)
data(GBSG)
fit1 <- coxph(Surv(rfst, cens) ~ age + prm+htreat+tumsize+menostat, data = GBSG, x=TRUE)
shrink(fit1, method = "dfbeta", type = "parameterwise")
shrink(fit1, method = "jackknife", type = "parameterwise")
shrink(fit1, method = "dfbeta", type = "parameterwise", join = list(c("age", "prm"), c("htreat1", "tumsize", "menostat2")))
shrink(fit1, method = "j", type = "parameterwise", join = list(c("age", "prm"), c("htreat1", "tumsize", "menostat2")))
shrink(fit1, method = "dfbeta", type = "global")
shrink(fit1, method = "jackknife", type = "global")

fit1 <- coxph(Surv(rfst, cens) ~ age + prm+htreat+tumsize+menostat, data = GBSG)
shrink(fit1, method = "dfbeta", type = "parameterwise")

utils::data(Pima.te, package="MASS")      
utils::data(Pima.tr, package="MASS")      
Pima <- rbind(Pima.te, Pima.tr)
nrow(Pima)
fit3 <- glm(type ~ npreg + glu + bmi + ped + age, family = binomial, data = Pima, x = TRUE) 
shrink(fit3, "parameterwise", "dfbeta")
shrink(fit3, "parameterwise", "dfbeta", join=list(c("npreg", "glu")))
shrink(fit3, "parameterwise", "dfbeta", join=list(c("npreg", "glu", "bmi")))
shrink(fit3, "parameterwise", "dfbeta", join=list(c("npreg", "glu", "bmi","ped","age")))
shrink(fit3, "parameterwise", "dfbeta", join=list(c("npreg", "glu"),c("bmi","ped","age")))
shrink(fit3, "parameterwise", "dfbeta", join=list(c("npreg", "glu"),c("bmi"), c("ped","age")))

fit3 <- glm(type ~ npreg + glu + bmi + ped + age, family = binomial, data = Pima) 
shrink(fit3, "parameterwise", "dfbeta")

# cox and mfp
library(mfp)
data(GBSG)
fit1 <- mfp(Surv(rfst, cens) ~ fp(age, df = 4, select = 0.05) +
            fp(prm, df = 4, select = 0.05), family = cox, data = GBSG)    # + fp(esm, df = 4, select = 1)

dfbeta.pw <- shrink(fit1, method = "dfbeta", type = "parameterwise")
dfbeta.pw

names(dfbeta.pw)
dfbeta.pw$shrinkage
dfbeta.pw$vcov.shrinkage
dfbeta.pw$shrunken
dfbeta.pw$postfit
dfbeta.pw$fit
dfbeta.pw$method
dfbeta.pw$type
dfbeta.pw$call

cov2cor(dfbeta.pw$vcov.shrinkage)

jack.pw <- shrink(fit1, method = "jackknife", type = "parameterwise")
jack.pw
jack.pw$vcov.shrinkage

dfbeta.join <- shrink(fit1, method = "dfbeta", type = "parameterwise",
                    join = list(c("age.1", "age.2")))
dfbeta.join
dfbeta.join$vcov.shrinkage

jack.join <- shrink(fit1, method = "jackknife", type = "parameterwise",
                  join = list(c("age.1", "age.2")))
jack.join
jack.join$vcov.shrinkage

dfbeta.global <- shrink(fit1, method = "dfbeta", type = "global")
dfbeta.global
dfbeta.global$vcov.shrinkage

jack.global <- shrink(fit1, method = "jackknife", type = "global")
jack.global
jack.global$vcov.shrinkage


plot(20:80, predict(fit1, newdata=data.frame(age=20:80, prm=32.5), "lp") -
     predict(fit1, newdata=data.frame(age=50, prm=32.5), "lp"), xlab="Age",
     ylab="Log hazard relative to 50 years", type="l", lwd=2)
lines(20:80, predict(dfbeta.global,newdata=data.frame(age=20:80, prm=32.5),"lp") -
      predict(dfbeta.global, newdata=data.frame(age=50, prm=32.5), "lp"), lty=4, col="blue", lwd=2)
lines(20:80, predict(dfbeta.join, newdata=data.frame(age=20:80, prm=32.5),"lp") -
      predict(dfbeta.join, newdata=data.frame(age=50, prm=32.5), "lp"), lty=3, col="red", lwd=2)
lines(20:80, predict(dfbeta.pw, newdata=data.frame(age=20:80, prm=32.5), type="lp") -
      predict(dfbeta.pw, newdata=data.frame(age=50, prm=32.5), type="lp"), lty=2, col="green", lwd=2)
legend("topright", lty=c(1,4,3,2), legend=c("No", "Global", "Joint", "Parameterwise"),
       title="SHRINKAGE", inset=0.01, bty="n", col=c("black", "blue", "red", "green"), lwd=2)

#-------------------------------------------------------------------------------
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

fit2 <- glm(y ~ -1 + x1 + x2, family = binomial, data = simdat, x = TRUE)
summary(fit2)

#fit2 <- glm(y ~ -1 + x1, family = binomial, data = simdat, x = TRUE)

#fit=fit2; method="d"; type="p"; join=NULL; postfit=TRUE
#join = list(c("x1", "x2"))

dfbeta.pw <- shrink(fit2, method = "dfbeta", type = "parameterwise")  
dfbeta.pw

names(dfbeta.pw)
dfbeta.pw$shrinkage
dfbeta.pw$vcov.shrinkage
dfbeta.pw$shrunken
dfbeta.pw$postfit
dfbeta.pw$fit
dfbeta.pw$method
dfbeta.pw$type
dfbeta.pw$call

dfbeta.join <- shrink(fit2, method = "dfbeta", type = "parameterwise", join = list(c("x1", "x2")))
dfbeta.join
dfbeta.join$vcov.shrinkage

jack.pw <- shrink(fit2, method = "jackknife", type = "parameterwise")
jack.pw
jack.pw$vcov.shrinkage

dfbeta.global <- shrink(fit2, method = "dfbeta", type = "global")
dfbeta.global
dfbeta.global$vcov.shrinkage

jack.global <- shrink(fit2, method = "jackknife", type = "global")
jack.global
jack.global$vcov.shrinkage

plot(-2:4, predict(fit2, newdata=data.frame(x1=-2:4, x2=0), "response") -
     predict(fit2, newdata=data.frame(x1=0, x2=0), "response"), xlab="Age",
     ylab="Log hazard relative to 50 years", type="l", lwd=2)
lines(-2:4, predict(dfbeta.global,newdata=data.frame(x1=-2:4, x2=0), "response") -
      predict(dfbeta.global, newdata=data.frame(x1=0, x2=0), "response"), lty=4, col="blue", lwd=2)
lines(-2:4, predict(dfbeta.join, newdata=data.frame(x1=-2:4, x2=0),"response") -
      predict(dfbeta.join, newdata=data.frame(x1=0, x2=0), "response"), lty=3, col="red", lwd=2)
lines(-2:4, predict(dfbeta.pw, newdata=data.frame(x1=-2:4, x2=0), "response") -
      predict(dfbeta.pw, newdata=data.frame(x1=0, x2=0), type="response"), lty=2, col="green", lwd=2)
legend("topright", lty=c(1,4,3,2), legend=c("No", "Global", "Joint", "Parameterwise"),
       title="SHRINKAGE", inset=0.01, bty="n", col=c("black", "blue", "red", "green"), lwd=2)


#-------------------------------------------------------------------------------
utils::data(anorexia, package="MASS")
fit5 <- lm(Postwt ~ Prewt + Treat, data =anorexia, x = TRUE, y = TRUE)
fit5

dfbeta.pw <- shrink(fit5, method = "dfbeta", type = "parameterwise")
dfbeta.pw
dfbeta.pw$vcov.shrinkage

dfbeta.join <- shrink(fit5, method = "dfbeta", type = "parameterwise",
                    join=list(c("TreatCont", "TreatFT")))
dfbeta.join
dfbeta.join$vcov.shrinkage

jack.pw <- shrink(fit5, method = "jackknife", type = "parameterwise")
jack.pw
jack.pw$vcov.shrinkage

dfbeta.global <- shrink(fit5, method = "dfbeta", type = "global")
dfbeta.global
dfbeta.global$vcov.shrinkage

jack.global <- shrink(fit5, method = "jackknife", type = "global")
jack.global
jack.global$vcov.shrinkage

plot(70:95, predict(fit5, newdata=data.frame(Prewt=70:95, Treat="Cont"), type="response") -
     predict(fit5, newdata=data.frame(Prewt=82, Treat="Cont"), type="response"), type="l", lwd=2)
lines(70:95, predict(dfbeta.global,newdata=data.frame(Prewt=70:95, Treat="Cont"), type="response") -
      predict(dfbeta.global, newdata=data.frame(Prewt=82, Treat="Cont"), type="response"), lty=4, col="blue", lwd=2)
lines(70:95, predict(dfbeta.join, newdata=data.frame(Prewt=70:95, Treat="Cont"), type="response") -
      predict(dfbeta.join, newdata=data.frame(Prewt=82, Treat="Cont"), type="response"), lty=3, col="red", lwd=2)
lines(70:95, predict(dfbeta.pw, newdata=data.frame(Prewt=70:95, Treat="Cont"), type="response") -
      predict(dfbeta.pw, newdata=data.frame(Prewt=82, Treat="Cont"), type="response"), lty=2, col="green", lwd=2)

#-------------------------------------------------------------------------------

utils::data(Pima.te, package="MASS")      
utils::data(Pima.tr, package="MASS")      
Pima <- rbind(Pima.te, Pima.tr)
Pima$type2 <- as.numeric(Pima$type)-1
fit3 <- mfp(type2 ~ npreg + glu + bmi + ped + fp(age, select = 0.05), family = binomial, data = Pima) 
fit3

dfbeta.pw <- shrink(fit3, type = "parameterwise", method="dfbeta")
dfbeta.global <- shrink(fit3, type = "global", method = "dfbeta")

predict(fit3,  type="link")[1:10]
predict(dfbeta.global, type="link")[1:10]
predict(dfbeta.pw, type="link")[1:10]


utils::data(GAGurine, package="MASS")
fit6 <- mfp(Age ~ fp(GAG, select = 0.05), family = gaussian, data = GAGurine)
#fit6 <- glm(Age ~ GAG, family = gaussian, data = GAGurine, x=TRUE)

dfbeta.pw <- shrink(fit6, type = "parameterwise", method = "dfbeta")
dfbeta.global <- shrink(fit6, type = "global", method = "dfbeta")

predict(fit6,  type="link")[1:10] 
predict(dfbeta.global, type="link")[1:10]
predict(dfbeta.pw, type="link")[1:10]

# unshrunken
plot(3:55, predict(fit6, newdata=data.frame(GAG=3:55), "link") -
     predict(fit6, newdata=data.frame(GAG=10), "link"), xlab="GAG",
     ylab="Log hazard relative to 50 years", type="l", lwd=2)

# globally shrunken
lines(3:55, predict(dfbeta.global, newdata=data.frame(GAG=3:55), type="link") - 
      predict(dfbeta.global, newdata=data.frame(GAG=10), type="link"), lty=4, col="blue", lwd=2)

# parameterwise shrunken
lines(3:55, predict(dfbeta.pw, newdata=data.frame(GAG=3:55), type="link") -
      predict(dfbeta.pw, newdata=data.frame(GAG=10), type="link"), lty=2, col="green", lwd=2)

legend("topright", lty=c(1,4,2), legend=c("No", "Global", "Parameterwise"), 
       title="SHRINKAGE", inset=0.01, bty="n", col=c("black", "blue", "green"), lwd=2)


utils::data(anorexia, package = "MASS")
nrow(anorexia)
contrasts(anorexia$Treat) <- contr.treatment(3, base = 2)
anorexia2 <- data.frame(cbind(Postwt=anorexia[, "Postwt"], model.matrix(~ Prewt + Treat, data = anorexia)[, -1]))
fit4 <- glm(Postwt ~ -1 +Prewt + Treat1 + Treat3, family = gaussian, data = anorexia2, x = TRUE)
fit4

shrink(fit4, type = "parameterwise", method = "dfbeta")
shrink(fit4, type = "parameterwise", method = "jackknife")
shrink(fit4, type = "parameterwise", method = "dfbeta", join = list(c("Treat1", "Treat3")))
shrink(fit4, type = "parameterwise", method = "jackknife", join = list(c("Treat1", "Treat3")))
shrink(fit4, type = "global", method = "dfbeta")
shrink(fit4, type = "global", method = "jackknife")


fit4 <- glm(Postwt ~ -1 +Prewt + Treat, family = gaussian, data = anorexia, x = TRUE)
fit4
shrink(fit4, type = "parameterwise", method = "dfbeta")
shrink(fit4, type = "parameterwise", method = "jackknife")
shrink(fit4, type = "parameterwise", method = "dfbeta", join = list(c("TreatCont", "TreatFT")))
shrink(fit4, type = "parameterwise", method = "jackknife", join = list(c("TreatCont", "TreatFT")))
shrink(fit4, type = "global", method = "dfbeta")
shrink(fit4, type = "global", method = "jackknife")
