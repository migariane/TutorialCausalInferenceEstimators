
#############################################################
# Luke Keele 1 & Miguel Angel Luque Fernandez 2
# 1 Georgetown University 
# 2 London School of Hygiene and Tropical Medicine
# Estimating the Average Treatment Effect (ATE)
# Using: 
# 1. Standardization (G-Formula), IP Weighting (G-Formula)
# 2. MSM
# 3. Double Robust Methods AIPTW 
# 4. TMLE
# 5. Matching 
################################################################

library(sandwich)
library(tmle)
library(MatchIt)
library(geepack)

rm(list=ls())

#data
L <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1)
A <- c(0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
Y <- c(1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1)
data <- as.data.frame(cbind(L,A,Y))

# Imbalance
# Level of X When Treated
mean(data$A[data$L==1], na.rm=TRUE)

# Level of X When Untreated
mean(data$A[data$L==0], na.rm=TRUE)

# Unadjusted Estimate
mean(data$Y[data$A==1], na.rm=TRUE) - mean(data$Y[data$A==0],na.rm=TRUE)

# Unadjusted Estimate Regression 
reg <- lm(Y ~ A, data=data); reg

# Naive Estimate Regression Adjustment 
naive <- lm(Y ~ A + L, data=data); naive

# G-Formula Non-Parametric = Standardization (Problem Dimensionallity)
pr.l <- prop.table(table(data$L))
pr.l

# Table of Classified Means
tab.out <- aggregate(Y ~ A + L, data, mean)
tab.out

# Non-parametric ATE 
((mean(data$Y[data$A==1 & data$L==1]) - mean(data$Y[data$A==0 & data$L==1]))*pr.l[2]) + 
(mean(data$Y[data$A==1 & data$L==0]) - mean(data$Y[data$A==0 & data$L==0]))*pr.l[1]

# G-Formula with Regression Adjustment (RA)
mod1  <- glm(Y ~  L, family="binomial", data=data[data$A==1,])
mod0  <- glm(Y ~  L, family="binomial", data=data[data$A==0,])
L <- data$L
PORA  <- cbind(Y1 = predict(mod1, newdata=data.frame(A = 1, L), type="response"),
               Y0 = predict(mod0, newdata=data.frame(A = 0, L), type="response"))
RAPO <- as.data.frame(PORA)
Y.1 <- RAPO$Y1
Y.0 <- RAPO$Y0
# G-computation ATE (95%CI Bootstraping?)
RApsi <- mean((Y.1) - (Y.0), na.rm=TRUE); RApsi

# Bootstrapping the G-comp esitmator
# G-computation wrapper function
library(boot)
attach(data)
g.comp=function(data,indices)
{
  dat=data[indices,]
  glm1  <- glm(Y ~ L, family="binomial", dat=dat[dat$A==1,])
  glm2  <- glm(Y ~ L, family="binomial", dat=dat[dat$A==0,])
  Y.1 = predict(glm1, newdata=data.frame(A = 1, L), type="response")
  Y.0 = predict(glm2, newdata=data.frame(A = 0, L), type="response")
  mean((Y.1) - mean(Y.0))
}
# Can get original estimate, by plugging in indices 1:n
g.comp(data,indices=1:nrow(data))
# Draw 200 bootstrap sample estimates of RD
boot.out=boot(data,g.comp,500)
# compute confidence intervals using percentile method
boot.ci(boot.out,type="perc",conf=0.95)
boot.ci(boot.out,type="norm",conf=0.95)

# IPW
# Estimate Propensity Score
# Joint Probability Distribution of X and D
table(data$L, data$A)
prop.table(table(data$L, data$A))

# P-scores by Hand
.30/pr.l[2]
.20/pr.l[1]
.05/pr.l[2]
.45/pr.l[1]

# Weights
1/(.30/pr.l[2])
1/(.20/pr.l[1])

1/(.05/pr.l[2])
1/(.45/pr.l[1])

# Via Logit
p.a <- glm(A ~ as.factor(L), data=data, family=binomial)
p.score <- ifelse(data$A == 0, 1 - predict(p.a, type = "response"),
                  predict(p.a, type = "response"))

# Table of Propensity Scores
table(p.score)

# Generate IP Weights                     
data$w <- 1/p.score
table(data$w)
summary(data$w)
sd(data$w)

# Estimate ATE
mean(data$w*as.numeric(data$A==1)*data$Y) - mean(data$w*as.numeric(data$A==0)*data$Y)

# Estimate ATE via modelling in MSM
(glm(Y ~ A + L, data=data, weights=data$w))

# Boostraping 
library(boot)
iptw.w = function(data,indices)
{
    dat=data[indices,]
    mean(dat$w*as.numeric(dat$A==1)*dat$Y) - mean(dat$w*as.numeric(dat$A==0)*dat$Y)
}
# Can get original estimate, by plugging in indices 1:n
iptw.w(data,indices=1:nrow(data))
# Draw 200 bootstrap sample estimates 
boot.out=boot(data,iptw.w,100)
# compute confidence intervals using percentile method
boot.ci(boot.out,type="perc",conf=0.95)
boot.ci(boot.out,type="norm",conf=0.95)

# Estimate ATE SLIPTW
#devtools::install_github("benkeser/drtmle")
library(drtmle)
library(SuperLearner)
attach(data)
L <- as.data.frame(L)
fit <- islptw(Y = Y, A = A, W = L, a_0 = c(0,1),
              SL_g = c("SL.glm","SL.glm.interaction","SL.step.interaction"),
              SL_Qr = "SL.glm","SL.glm.interaction","SL.step.interaction")
fit <- confint(fit); fit 
fit <- islptw(Y = Y, A = A, W = L, a_0 = c(0,1),
              SL_g = c("SL.glm"),
              SL_Qr = "SL.glm")
ci_ate <- confint(fit, contrast=c(-1,1)); ci_ate

# MSM
msm <- lm(Y  ~ A, data = data, weights = w)

# Now with Inference
SE <-sqrt(diag(vcovHC(msm, type="HC0"))) # robust standard errors
beta <- coef(msm)
lcl <- beta-1.96*SE 
ucl <- beta+1.96*SE
cbind(beta, lcl, ucl)[2,]

# Stablized Weights  Fix this.
mean(data$w*as.numeric(data$A==1)*data$Y)/mean( data$w*as.numeric(data$A==1)) - mean(data$w*as.numeric(data$A==0)*data$Y)/mean(data$w*as.numeric(data$A==0))

denom.fit <- glm(A ~ as.factor(L), 
                 family = binomial(), data = data)
denom.p <- predict(denom.fit, type = "response")

# estimation of numerator of ip weights
numer.fit <- glm(A ~ 1, family = binomial(), data = data)
summary(numer.fit)
numer.p <- predict(numer.fit, type = "response")

data$sw <- ifelse(data$A == 0, ((1-numer.p)/(1-denom.p)),
                  (numer.p/denom.p))

msm <- lm(Y  ~ A, data = data, weights = sw)

# Now with Inference
SE <-sqrt(diag(vcovHC(msm, type="HC0"))) # robust standard errors
beta <- coef(msm)
lcl <- beta-1.96*SE 
ucl <- beta+1.96*SE
cbind(beta, lcl, ucl)[2,]

# Double Robust Methods
# AIPTW 
# Potential Outcomes
mod  <- glm(Y ~ A + L, family="binomial", data=data)
PO   <- cbind(Yhat = predict(mod),
             Y1 = predict(mod, newdata=data.frame(A = 1, L), type="response"),
             Y0 = predict(mod, newdata=data.frame(A = 0, L), type="response"))
RA <- as.data.frame(PO)
Yhat <- RA$Yhat
Y.1a <- RA$Y1
Y.0a <- RA$Y0
# estimation of weights
g <- glm(A ~ L, family = binomial(), data = data)
gw <- predict(g, type = "response")
gws <- ifelse(data$A == 0, (-(1 - data$A)/(1 - gw)),(data$A/gw))
sum(gws)
AIPTW <- mean(gws*(data$Y - plogis(RA$Yhat)) + ((Y.1a) - (Y.0a)));AIPTW
# Inference Influence Curve
IC <- (gws*(data$Y - plogis(RA$Yhat)) + ((Y.1a) - (Y.0a)))-AIPTW
n  <- nrow(data)
varHat.IC <- var(IC)/n; varHat.IC
lci <- AIPTW-1.96*sqrt(varHat.IC)
uci <- AIPTW+1.96*sqrt(varHat.IC)
cat(AIPTW,lci,uci)

# Inference using Bootstrap
AIPTW.b = function(data,indices)
  {
  dat=data[indices,]
  mod  <- glm(Y ~ A + L, family="binomial", data=dat)
  Yhat = predict(mod)
  Y1 = predict(mod, newdata=data.frame(A = 1, L))
  Y0 = predict(mod, newdata=data.frame(A = 0, L))
  g <- glm(A ~ L, family="binomial",  data = data)
  gw <- predict(g,type="response")
  gws <- ifelse(A == 0, (-(1 - A)/(1 - gw)),(A/gw))
  mean(gws*(Y - plogis(Yhat)) + (plogis(Y1) - plogis(Y0)))
  }
# Can get original estimate, by plugging in indices 1:n
AIPTW.b(data,indices=1:nrow(data))
# Draw 200 bootstrap sample estimates 
boot.out=boot(data,AIPTW.b,100)
# compute confidence intervals using percentile method
boot.ci(boot.out,type="perc",conf=0.95)
boot.ci(boot.out,type="norm",conf=0.95)

# TMLE
set.seed(7777)
library(tmle)
w <- subset(data, select=c(L))
fittmle <- tmle(data$Y, data$A, W=w, family="binomial", Q.SL.library = c("SL.glm","SL.glm.interaction","SL.step.interaction","SL.gam","SL.randomForest"),
             g.SL.library = c("SL.glm","SL.glm.interaction","SL.step.interaction","SL.gam","SL.randomForest"))
fittmle
# Matching Methods
# Matching (Two methods: nearest neighbors with logit and optimal matching with mahalanobis)
match.1 <- matchit(Y ~ A + L, 
                   data = data, method = "nearest", distance = "logit", reestimate = TRUE)

match.2 <- matchit(Y ~ A + L,
                   data = data, method = "optimal", distance = "mahalanobis", reestimate = TRUE)

##### Compare the results from each of the matches
summary(match.1)
summary(match.2)

##### Create suitable data structures
m1 <- match.data(match.1)
m2 <- match.data(match.2)
mmat1 <- as.numeric(match.1$match.matrix)
mmat2 <- as.numeric(match.2$match.matrix)

m1data <- data.frame(rbind(
  cbind(1:(dim(m1)[1]/2), m1$Y[m1$A == 1], m1$A[m1$A == 1]),
  cbind(1:(dim(m1)[1]/2), data$Y[mmat1], data$A[mmat1]))
)
names(m1data) <- c("id", "Y", "A")

m2data <- data.frame(rbind(
  cbind(1:(dim(m2)[1]/2), m2$Y[m2$A == 1], m2$A[m2$A == 1]),
  cbind(1:(dim(m2)[1]/2), data$Y[mmat2], data$A[mmat2]))
)
names(m2data) <- c("id", "Y", "A")

##### Logistic regression based on matched pairs
zz.match1 <- geese(Y ~ A, family = binomial(link = "log"), id = id, data = m1data)
summary(zz.match1)

zz.match2 <- geese(Y ~ A, family = binomial(link = "log"), id = id, data = m2data)
summary(zz.match2)

## Add more covariates and use Non-parametric approach from LUKEE
