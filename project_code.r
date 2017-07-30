rm(list=ls())
ls()
a = read.csv("3000cases_NE_Adlt_AllDay_2.csv")
nrow(a)
attach(a)
data.frame(a$HHSIZE)
proj = data.frame(CNTTDTR
                       , DRVR
                       , HHSIZE
                       , HHVEHCNT
                       , ONEVEH
                       , TWOVEH
                       , THREEVEH
                       , DISTWRK
                       , WRKR
                       , WRKR_SALES
                       , WRKR_ADMIN
                       , WRKR_MFGCONST
                       , WRKR_PROTECH
                       , WRKTIMEFLEX
                       , HIGHINC
                       , MEDINC
                       , YOUNG
                       , OLD
                       , URB
                       , MALE
                       , Weekend, SMTRVLTIME)
proj = na.omit(proj)
nrow(proj)
head(proj)

pairs(~SMTRVLTIME+CNTTDTR+DISTWRK+ ONEVEH+TWOVEH+THREEVEH+URB+MEDINC+HIGHINC+
        YOUNG+OLD+URB+ MALE+ Weekend + WRKR_SALES
      + WRKR_ADMIN
      + WRKR_MFGCONST
      + WRKR_PROTECH
      + WRKTIMEFLEX,
      data=proj, main="Simple Scatterplot Matrix")
corrvariable= data.frame(log(SMTRVLTIME), CNTTDTR, DISTWRK, ONEVEH,TWOVEH,
                         THREEVEH,URB,MEDINC,HIGHINC, YOUNG
                         , OLD
                         , URB
                         , MALE
                         , Weekend)
k = cor(corrvariable, method="pearson")
as.matrix(k)
attach(proj)

null = lm(SMTRVLTIME ~ 1, data = proj)
full = lm(SMTRVLTIME ~ CNTTDTR
          + DRVR
          + HHSIZE
          + HHVEHCNT
          + ONEVEH
          + TWOVEH
          + THREEVEH
          + DISTWRK
          + WRKR
          + WRKR_SALES
          + WRKR_ADMIN
          + WRKR_MFGCONST
          + WRKR_PROTECH
          + WRKTIMEFLEX
          + HIGHINC
          + MEDINC
          + YOUNG
          + OLD
          + URB
          + MALE
          + Weekend)

library(leaps)
x  = nrow(proj)

step(null, scope = list(upper=full, lower=null), 
     data=data, direction="both", nbest = 2)
step(null, scope = list(upper=full, lower=null), 
     data=flu, direction="both", k = log(x))

testmod= lm((SMTRVLTIME) ~ (CNTTDTR) + DISTWRK + HIGHINC +MEDINC + URB + Weekend + MALE  + (WRKR) + THREEVEH)
summary(testmod)
extractAIC(testmod) 
extractAIC(testmod, k=log(2999))
library(lmtest)
bptest(testmod, data = proj)
abline(0,1)

qqnorm(resid(testmod))
plot(fitted(testmod), rstudent(testmod), main="Externally Studentized Residuals Plot Vs Fitted Value", 
            xlab= "Fitted Y", ylab="Ext Studentized Residuals" )
abline(0,0)
abline(2,0, col=4)
abline(-2,0, col=4)

#Leverage
p=length(coef(testmod))
n=nrow(proj)
hii = influence(testmod)$hat
par(mfrow=c(1,1))
plot(hii, main="Hat Diagonals")
abline(h=2*p/n, col="red")

#BC Transforation
library(car)
p1 <- powerTransform((SMTRVLTIME) ~  CNTTDTR + URB+ DISTWRK  
                     +MEDINC + MALE  + WRKR 
                     +THREEVEH, proj)
summary(p1)
p1$roundlam
m1 <- lm(bcPower(SMTRVLTIME, p1$roundlam) ~  CNTTDTR + URB+ DISTWRK  
         +MEDINC + MALE  + WRKR 
         +THREEVEH, data =proj)
summary(m1)
qqnorm(resid(m1), main="Normal Q-Q Plot of Residuals")
abline(0,1)
plot(fitted(m1), rstudent(m1), 
     main="Residual Vs Fitted Values Plot", 
     xlab="Fitted Values", ylab="Residuals")
abline(0,0, col=1)
abline(2,0, col=4)
abline(-2,0, col=4)
library(lmtest)
dwtest(m1,data=proj)
#INFLUNETIAL OUTLIERS
hii = influence(m1)$hat
p=length(coef(m1))
2*p/n
highleverage = hii[hii > 2*p/n]
highinfluentials

CooksD<-  cooks.distance(m1)
ri    <- rstandard(m1)
dfits  <- dffits(m1)
cov.rat<- covratio(m1)
p=length(coef(m1))
n=nrow(proj)

par(mfrow=c(1,1))
plot(hii, main="Hat Diagonals")
abline(h=2*p/n, col="red")

plot(CooksD, main="Cook's Distances")
abline(h=.8,col="red")

plot(abs(dfits), main="DFFITS")
abline(h=2*sqrt( p/n ), col="red" )
plot(abs(cov.rat-1), main="|COVRATIO - 1|")
abline(h=3*p/n, col="red")

highleverage = hii[hii > 2*p/n]
highleverage
highCooksD = CooksD[CooksD > 0.8]
highdffits = dffits[dffits > 2*sqrt( p/n )]
highcovratio = covratio[covratio > abs(cov.rat-1)]
sum.mat <- cbind( hii,CooksD, dfits, cov.rat)

proj = proj[-1145,]
nrow(proj)
pr= na.omit(pr)
sort(fitted(m1), decreasing=TRUE)
which.max( fitted(m1))
order(fitted(m1),decreasing=T)[1:10]
# unweighted regression
unwt.mod <- lm(log(SMTRVLTIME) ~ (CNTTDTR) + URB+ DISTWRK 
                 +MEDINC + MALE  + WRKR 
               +THREEVEH,data=proj)
res <- unwt.mod$residuals
yhat<- unwt.mod$fitted.values
# estimate sd
sd.mod <- lm(abs(res)~ (CNTTDTR) + URB+ DISTWRK 
             +MEDINC + MALE  + WRKR 
             +THREEVEH,data=proj)
sest <- sd.mod$fitted.values
wt <- 1/sest^2
# weighted regression
wt.mod <- lm(log(SMTRVLTIME) ~ (CNTTDTR) + URB+ DISTWRK  
             +MEDINC + MALE  + WRKR 
             +THREEVEH,data=proj)
resw <- wt.mod$residuals
ywhat <- wt.mod$fitted.values
par(mfrow=c(1,2))
plot(yhat, res, main="Unweighted Model", xlab="predicted", ylab="residuals")
abline(h=0, col="red")
plot(ywhat, sqrt(wt)*resw, main="Weighted Model", xlab="predicted", ylab="residuals")
abline(h=0, col="red")
qqnorm(resid(wt.mod), main="Normal Q-Q Plot of Residuals")
abline(0,1)
summary(wt.mod)
summary(m1)
sort(sqrt(wt)*resw,decreasing=F)[1:6]
sort(ywhat,decreasing=T)[1:10]
bptest(wt.mod, data=proj, studentize = FALSE)
library(alr3)
pureErrorAnova(wt.mod)
dwtest(wt.mod, data=proj)

#Hold-out
rm(list=ls())
ls()
hold = read.csv("1803cases_NE_Adlt_AllDay.csv")
attach(hold)
hold = na.omit(hold)
nrow(hold)
hold$predict = NULL
hold$predict = exp(3.479570 + 0.151877 *CNTTDTR - 0.097403* URB
                   + 0.014379 *DISTWRK  -0.084656 *MEDINC
                   +  0.057195  * MALE  -0.092133   *WRKR
                   +0.093728  * THREEVEH)
