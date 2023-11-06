## Prerequisites ##

library(MASS)
library(survival)
library(flexsurv)
library(pec)
library(rms)
library(randomForestSRC)
library(alabama)

## Simulation setting ##

#path.test <- "C:/Users/haoli/Desktop/Dissertation/04_project1/03_simulations/02_data_generation/RC/test/test/"
#path.c <- "C:/Users/haoli/Desktop/Dissertation/04_project1/03_simulations/02_data_generation/RC/test/c/"
#path.cc <- "C:/Users/haoli/Desktop/Dissertation/04_project1/03_simulations/02_data_generation/RC/test/cc/"
#path.sc <- "C:/Users/haoli/Desktop/Dissertation/04_project1/03_simulations/02_data_generation/RC/test/sc/"
#path.results <- "C:/Users/haoli/Desktop/Dissertation/04_project1/03_simulations/03_results/01_RC_test/"

path.test <- "/nas/longleaf/home/haolin/dissertation1/round3/data/test/"
path.c <- "/nas/longleaf/home/haolin/dissertation1/round3/data/c/"
#path.r <- "/nas/longleaf/home/haolin/dissertation1/round3/data/r/"
#path.cc <- "/nas/longleaf/home/haolin/dissertation1/round3/data/cc/"
#path.sc <- "/nas/longleaf/home/haolin/dissertation1/round3/data/sc/"
path.results <- "/nas/longleaf/home/haolin/dissertation1/round3/results/"

p = 10
nsim = 500
batch = 1
seq = ((batch-1)*nsim+1):(batch*nsim)
tau = 0.0162
fold = 10
L = 4
gamma.weibull = 0.7

model = 'RC'
spec = 'linear'
cen.mod = 'ind'

ipcw.wt.cum = matrix(0, nrow = nsim, ncol = L)

exp.pa.cum = rep(0, nsim)
cox.pa.cum = rep(0, nsim)
gam.pa.cum = rep(0, nsim)
srf.pa.cum = rep(0, nsim)
ipcwsl.pa.cum = rep(0, nsim)


for (m in seq){
  cat(m)
  try({
  # 1. Prepare data for k-fold CV #
  
  dat = read.csv(paste0(path.c, 'dat', m, '.csv'))
  n = nrow(dat)
  #dat$subid = c(1:n)
  dat$fold = 0
  for (i in 1:n){
    dat$fold[i] = sample(c(1:fold), 1)
  }
  
  # 2. Model fitting of candidate learners and prediction for k-fold CV #
  
  # 2.1 Cox model #
  
  cox.survival.pred = data.frame(subid = rep(0, n), t1 = rep(0, n), t2 = rep(0, n), t3 = rep(0, n), t4 = rep(0, n),
                           t5 = rep(0, n), t6 = rep(0, n), t7 = rep(0, n), t8 = rep(0, n), t9 = rep(0, n))
  for (j in 1:fold){
    dat.train = dat[dat$fold!=j, ]
    dat.test = dat[dat$fold==j, ]
    cox.fit = cph(Surv(time, status) ~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10 , data = dat.train, surv=T)
    for (i in 1:nrow(dat.test)){
      subid = as.numeric(dat.test$subid[i])
      cox.survival.pred$subid[subid] = subid
      cox.survival.pred$t1[subid] = predictSurvProb(cox.fit, newdata=dat.test[i, 1:p], time = 1*tau/10)
      cox.survival.pred$t2[subid] = predictSurvProb(cox.fit, newdata=dat.test[i, 1:p], time = 2*tau/10)
      cox.survival.pred$t3[subid] = predictSurvProb(cox.fit, newdata=dat.test[i, 1:p], time = 3*tau/10)
      cox.survival.pred$t4[subid] = predictSurvProb(cox.fit, newdata=dat.test[i, 1:p], time = 4*tau/10)
      cox.survival.pred$t5[subid] = predictSurvProb(cox.fit, newdata=dat.test[i, 1:p], time = 5*tau/10)
      cox.survival.pred$t6[subid] = predictSurvProb(cox.fit, newdata=dat.test[i, 1:p], time = 6*tau/10)
      cox.survival.pred$t7[subid] = predictSurvProb(cox.fit, newdata=dat.test[i, 1:p], time = 7*tau/10)
      cox.survival.pred$t8[subid] = predictSurvProb(cox.fit, newdata=dat.test[i, 1:p], time = 8*tau/10)
      cox.survival.pred$t9[subid] = predictSurvProb(cox.fit, newdata=dat.test[i, 1:p], time = 9*tau/10)
    }
  }
  
  # 2.2 exponential model #
  
  exp.survival.pred = data.frame(subid = rep(0, n), t1 = rep(0, n), t2 = rep(0, n), t3 = rep(0, n), t4 = rep(0, n),
                                 t5 = rep(0, n), t6 = rep(0, n), t7 = rep(0, n), t8 = rep(0, n), t9 = rep(0, n))
  for (j in 1:fold){
    dat.train = dat[dat$fold!=j, ]
    dat.test = dat[dat$fold==j, ]
    exp.fit = flexsurvreg(Surv(time, status) ~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10, data = dat.train, dist="exp")
    for (i in 1:nrow(dat.test)){
      subid = as.numeric(dat.test$subid[i])
      exp.survival.pred$subid[subid] = subid
      exp.survival.pred$t1[subid] = predict(exp.fit, newdata=dat.test[i, 1:p], type = "survival", times = 1*tau/10)$.pred
      exp.survival.pred$t2[subid] = predict(exp.fit, newdata=dat.test[i, 1:p], type = "survival", times = 2*tau/10)$.pred
      exp.survival.pred$t3[subid] = predict(exp.fit, newdata=dat.test[i, 1:p], type = "survival", times = 3*tau/10)$.pred
      exp.survival.pred$t4[subid] = predict(exp.fit, newdata=dat.test[i, 1:p], type = "survival", times = 4*tau/10)$.pred
      exp.survival.pred$t5[subid] = predict(exp.fit, newdata=dat.test[i, 1:p], type = "survival", times = 5*tau/10)$.pred
      exp.survival.pred$t6[subid] = predict(exp.fit, newdata=dat.test[i, 1:p], type = "survival", times = 6*tau/10)$.pred
      exp.survival.pred$t7[subid] = predict(exp.fit, newdata=dat.test[i, 1:p], type = "survival", times = 7*tau/10)$.pred
      exp.survival.pred$t8[subid] = predict(exp.fit, newdata=dat.test[i, 1:p], type = "survival", times = 8*tau/10)$.pred
      exp.survival.pred$t9[subid] = predict(exp.fit, newdata=dat.test[i, 1:p], type = "survival", times = 9*tau/10)$.pred
    }
  }
  
  # 2.3 gamma AFT model #
  
  gam.survival.pred = data.frame(subid = rep(0, n), t1 = rep(0, n), t2 = rep(0, n), t3 = rep(0, n), t4 = rep(0, n),
                                 t5 = rep(0, n), t6 = rep(0, n), t7 = rep(0, n), t8 = rep(0, n), t9 = rep(0, n))
  for (j in 1:fold){
    dat.train = dat[dat$fold!=j, ]
    dat.test = dat[dat$fold==j, ]
    gam.fit = flexsurvreg(Surv(time, status) ~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10, data = dat.train, dist="gamma")
    for (i in 1:nrow(dat.test)){
      subid = as.numeric(dat.test$subid[i])
      gam.survival.pred$subid[subid] = subid
      gam.survival.pred$t1[subid] = predict(gam.fit, newdata=dat.test[i, 1:p], type = "survival", times = 1*tau/10)$.pred
      gam.survival.pred$t2[subid] = predict(gam.fit, newdata=dat.test[i, 1:p], type = "survival", times = 2*tau/10)$.pred
      gam.survival.pred$t3[subid] = predict(gam.fit, newdata=dat.test[i, 1:p], type = "survival", times = 3*tau/10)$.pred
      gam.survival.pred$t4[subid] = predict(gam.fit, newdata=dat.test[i, 1:p], type = "survival", times = 4*tau/10)$.pred
      gam.survival.pred$t5[subid] = predict(gam.fit, newdata=dat.test[i, 1:p], type = "survival", times = 5*tau/10)$.pred
      gam.survival.pred$t6[subid] = predict(gam.fit, newdata=dat.test[i, 1:p], type = "survival", times = 6*tau/10)$.pred
      gam.survival.pred$t7[subid] = predict(gam.fit, newdata=dat.test[i, 1:p], type = "survival", times = 7*tau/10)$.pred
      gam.survival.pred$t8[subid] = predict(gam.fit, newdata=dat.test[i, 1:p], type = "survival", times = 8*tau/10)$.pred
      gam.survival.pred$t9[subid] = predict(gam.fit, newdata=dat.test[i, 1:p], type = "survival", times = 9*tau/10)$.pred
    }
  }
  
  # 2.4 Survival random forest #
  
  srf.survival.pred = data.frame(subid = rep(0, n), t1 = rep(0, n), t2 = rep(0, n), t3 = rep(0, n), t4 = rep(0, n),
                            t5 = rep(0, n), t6 = rep(0, n), t7 = rep(0, n), t8 = rep(0, n), t9 = rep(0, n))
  for (j in 1:fold){
    dat.train = dat[dat$fold!=j, ]
    dat.test = dat[dat$fold==j, ]
    srf.fit = rfsrc(Surv(time, status) ~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10, data = dat.train)
    for (i in 1:nrow(dat.test)){
      srf.pred = data.frame(predict(srf.fit, newdata = dat.test[i,1:p])$survival[1,], predict(srf.fit, newdata = dat.test[i,1:p])$time.interest)
      subid = as.numeric(dat.test$subid[i])
      srf.pred$time1 = as.numeric(srf.pred[,2]-1*tau/10 <0)
      srf.pred$time2 = as.numeric(srf.pred[,2]-2*tau/10 <0)
      srf.pred$time3 = as.numeric(srf.pred[,2]-3*tau/10 <0)
      srf.pred$time4 = as.numeric(srf.pred[,2]-4*tau/10 <0)
      srf.pred$time5 = as.numeric(srf.pred[,2]-5*tau/10 <0)
      srf.pred$time6 = as.numeric(srf.pred[,2]-6*tau/10 <0)
      srf.pred$time7 = as.numeric(srf.pred[,2]-7*tau/10 <0)
      srf.pred$time8 = as.numeric(srf.pred[,2]-8*tau/10 <0)
      srf.pred$time9 = as.numeric(srf.pred[,2]-9*tau/10 <0)
      srf.survival.pred$subid[subid] = subid
      srf.survival.pred$t1[subid] = min(srf.pred[srf.pred$time1>0, 1])
      srf.survival.pred$t2[subid] = min(srf.pred[srf.pred$time2>0, 1])
      srf.survival.pred$t3[subid] = min(srf.pred[srf.pred$time3>0, 1])
      srf.survival.pred$t4[subid] = min(srf.pred[srf.pred$time4>0, 1])
      srf.survival.pred$t5[subid] = min(srf.pred[srf.pred$time5>0, 1])
      srf.survival.pred$t6[subid] = min(srf.pred[srf.pred$time6>0, 1])
      srf.survival.pred$t7[subid] = min(srf.pred[srf.pred$time7>0, 1])
      srf.survival.pred$t8[subid] = min(srf.pred[srf.pred$time8>0, 1])
      srf.survival.pred$t9[subid] = min(srf.pred[srf.pred$time9>0, 1])
    }
  }
  
  # 3. IPCW Super Learner #
  
  # 3.1 Indicators #
  
  ind = dat
  ind$t1 = as.numeric((ind$time > 1*tau/10))
  ind$t2 = as.numeric((ind$time > 2*tau/10))
  ind$t3 = as.numeric((ind$time > 3*tau/10))
  ind$t4 = as.numeric((ind$time > 4*tau/10))
  ind$t5 = as.numeric((ind$time > 5*tau/10))
  ind$t6 = as.numeric((ind$time > 6*tau/10))
  ind$t7 = as.numeric((ind$time > 7*tau/10))
  ind$t8 = as.numeric((ind$time > 8*tau/10))
  ind$t9 = as.numeric((ind$time > 9*tau/10))
  ind = subset(ind, select = c(subid, t1, t2, t3, t4, t5, t6, t7, t8, t9))
  
  # 3.2 Fraction in the front #
  
  cen.fit <- survfit(Surv(time, (1-status)) ~ 1 , data=dat)
  #cent.pred = summary(cen.fit, times = 0.01)$surv
  
  denom = dat
  denom$time1 = (as.numeric(denom$time <= 1*tau/10)*denom$time + as.numeric(denom$time > 1*tau/10)*1*tau/10) 
  denom$time2 = (as.numeric(denom$time <= 2*tau/10)*denom$time + as.numeric(denom$time > 2*tau/10)*2*tau/10) 
  denom$time3 = (as.numeric(denom$time <= 3*tau/10)*denom$time + as.numeric(denom$time > 3*tau/10)*3*tau/10) 
  denom$time4 = (as.numeric(denom$time <= 4*tau/10)*denom$time + as.numeric(denom$time > 4*tau/10)*4*tau/10) 
  denom$time5 = (as.numeric(denom$time <= 5*tau/10)*denom$time + as.numeric(denom$time > 5*tau/10)*5*tau/10) 
  denom$time6 = (as.numeric(denom$time <= 6*tau/10)*denom$time + as.numeric(denom$time > 6*tau/10)*6*tau/10) 
  denom$time7 = (as.numeric(denom$time <= 7*tau/10)*denom$time + as.numeric(denom$time > 7*tau/10)*7*tau/10) 
  denom$time8 = (as.numeric(denom$time <= 8*tau/10)*denom$time + as.numeric(denom$time > 8*tau/10)*8*tau/10) 
  denom$time9 = (as.numeric(denom$time <= 9*tau/10)*denom$time + as.numeric(denom$time > 9*tau/10)*9*tau/10) 
  
  denom$cen.pred1 = 0
  denom$cen.pred2 = 0
  denom$cen.pred3 = 0
  denom$cen.pred4 = 0
  denom$cen.pred5 = 0
  denom$cen.pred6 = 0
  denom$cen.pred7 = 0
  denom$cen.pred8 = 0
  denom$cen.pred9 = 0
  for (i in 1:n){
    denom$cen.pred1[i] = summary(cen.fit, times = denom$time1[i])$surv
    denom$cen.pred2[i] = summary(cen.fit, times = denom$time2[i])$surv
    denom$cen.pred3[i] = summary(cen.fit, times = denom$time3[i])$surv
    denom$cen.pred4[i] = summary(cen.fit, times = denom$time4[i])$surv
    denom$cen.pred5[i] = summary(cen.fit, times = denom$time5[i])$surv
    denom$cen.pred6[i] = summary(cen.fit, times = denom$time6[i])$surv
    denom$cen.pred7[i] = summary(cen.fit, times = denom$time7[i])$surv
    denom$cen.pred8[i] = summary(cen.fit, times = denom$time8[i])$surv
    denom$cen.pred9[i] = summary(cen.fit, times = denom$time9[i])$surv
  }
  
  denom$t1 = (1-as.numeric((denom$time <= 1*tau/10)&(denom$status == 0)))/denom$cen.pred1
  denom$t2 = (1-as.numeric((denom$time <= 2*tau/10)&(denom$status == 0)))/denom$cen.pred2
  denom$t3 = (1-as.numeric((denom$time <= 3*tau/10)&(denom$status == 0)))/denom$cen.pred3
  denom$t4 = (1-as.numeric((denom$time <= 4*tau/10)&(denom$status == 0)))/denom$cen.pred4
  denom$t5 = (1-as.numeric((denom$time <= 5*tau/10)&(denom$status == 0)))/denom$cen.pred5
  denom$t6 = (1-as.numeric((denom$time <= 6*tau/10)&(denom$status == 0)))/denom$cen.pred6
  denom$t7 = (1-as.numeric((denom$time <= 7*tau/10)&(denom$status == 0)))/denom$cen.pred7
  denom$t8 = (1-as.numeric((denom$time <= 8*tau/10)&(denom$status == 0)))/denom$cen.pred8
  denom$t9 = (1-as.numeric((denom$time <= 9*tau/10)&(denom$status == 0)))/denom$cen.pred9
  
  frac = subset(denom, select = c(subid, t1, t2, t3, t4, t5, t6, t7, t8, t9))
  
  # 3.3 Compute IPCW SL weights #
  
  fn <- function(x){
    fn = sum((as.matrix(ind[,-1])-(x[1]*as.matrix(exp.survival.pred[,-1]) + x[2]*as.matrix(gam.survival.pred[,-1]) + x[3]*as.matrix(cox.survival.pred[,-1])+
                                     x[4]*as.matrix(srf.survival.pred[,-1]) ))^2 * (as.matrix(frac[,-1]))) / (n)
    fn
  } 
  
  heq <- function(x) {
    h <- rep(NA, 1)
    h[1] <- x[1] + x[2] + x[3] + x[4] - 1
    h
  }
  
  hin <- function(x) {
    h <- rep(NA, 1)
    h[1] <- x[1]
    h[2] <- x[2]
    h[3] <- x[3]
    h[4] <- x[4]
    h
  }
  
  set.seed(12)
  p0 <- runif(4)
  ans <- constrOptim.nl(par=p0, fn=fn, heq=heq, hin=hin)
  ipcw.wt = ans$par
  
  ipcw.wt.cum[m,] = ipcw.wt
  
  
  # 5. Prediction Accuracy for candidate learners and super learner #
  
  # 5.0 Calculate the true survival probability #
  dat.pa = read.csv(paste0(path.test, 'dat', m, '.csv'))
  
  beta = matrix(c(0.1, 0.1, 0.3, -0.3, -0.2, -0.3, 0.01, -0.01, 0.25, -0.25), nrow = p)
  true.lin.pred <- as.matrix(dat.pa[,1:p]) %*% beta + 0.3*as.matrix(dat.pa[,1])*as.matrix(dat.pa[,3]) - 0.3*as.matrix(dat.pa[,3])^2 + 0.2*as.matrix(dat.pa[,1])*as.matrix(dat.pa[,6])*as.matrix(dat.pa[,9]) - 0.4*as.matrix(dat.pa[,2])*as.matrix(dat.pa[,4])*as.matrix(dat.pa[,5])+1.6
  
  survival.true <- function(t){
    survival =  exp(-exp(true.lin.pred)*t^gamma.weibull)
    return(survival)
  }
  true.survival = data.frame(subid = c(1:nrow(dat.pa)), t1 = rep(0, nrow(dat.pa)), t2 = rep(0, nrow(dat.pa)), 
                        t3 = rep(0, nrow(dat.pa)), t4 = rep(0, nrow(dat.pa)), t5 = rep(0, nrow(dat.pa)), 
                        t6 = rep(0, nrow(dat.pa)), t7 = rep(0, nrow(dat.pa)), t8 = rep(0, nrow(dat.pa)), 
                        t9 = rep(0, nrow(dat.pa)), t10 = rep(0, nrow(dat.pa)), t11 = rep(0, nrow(dat.pa)),
                        t12 = rep(0, nrow(dat.pa)), t13 = rep(0, nrow(dat.pa)), t14 = rep(0, nrow(dat.pa)))
  true.survival$t1 = survival.true(1*tau/15)
  true.survival$t2 = survival.true(2*tau/15)
  true.survival$t3 = survival.true(3*tau/15)
  true.survival$t4 = survival.true(4*tau/15)
  true.survival$t5 = survival.true(5*tau/15)
  true.survival$t6 = survival.true(6*tau/15)
  true.survival$t7 = survival.true(7*tau/15)
  true.survival$t8 = survival.true(8*tau/15)
  true.survival$t9 = survival.true(9*tau/15)
  true.survival$t10 = survival.true(10*tau/15)
  true.survival$t11 = survival.true(11*tau/15)
  true.survival$t12 = survival.true(12*tau/15)
  true.survival$t13 = survival.true(13*tau/15)
  true.survival$t14 = survival.true(14*tau/15)
  
  # 5.1 Exponential model #
  
  exp.survival.pa = data.frame(subid = rep(0, nrow(dat.pa)), t1 = rep(0, nrow(dat.pa)), t2 = rep(0, nrow(dat.pa)), 
                         t3 = rep(0, nrow(dat.pa)), t4 = rep(0, nrow(dat.pa)), t5 = rep(0, nrow(dat.pa)), 
                         t6 = rep(0, nrow(dat.pa)), t7 = rep(0, nrow(dat.pa)), t8 = rep(0, nrow(dat.pa)), 
                         t9 = rep(0, nrow(dat.pa)), t10 = rep(0, nrow(dat.pa)), t11 = rep(0, nrow(dat.pa)),
                         t12 = rep(0, nrow(dat.pa)), t13 = rep(0, nrow(dat.pa)), t14 = rep(0, nrow(dat.pa)))
  exp.fit.pa = flexsurvreg(Surv(time, status) ~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10, data = dat, dist="exp")
  for (i in 1:nrow(dat.pa)){
    subid = as.numeric(dat.pa$subid[i])
    exp.survival.pa$subid[subid] = subid
    exp.survival.pa$t1[subid] = predict(exp.fit.pa, newdata=dat.pa[i, 1:p], type = "survival", times = 1*tau/15)$.pred
    exp.survival.pa$t2[subid] = predict(exp.fit.pa, newdata=dat.pa[i, 1:p], type = "survival", times = 2*tau/15)$.pred
    exp.survival.pa$t3[subid] = predict(exp.fit.pa, newdata=dat.pa[i, 1:p], type = "survival", times = 3*tau/15)$.pred
    exp.survival.pa$t4[subid] = predict(exp.fit.pa, newdata=dat.pa[i, 1:p], type = "survival", times = 4*tau/15)$.pred
    exp.survival.pa$t5[subid] = predict(exp.fit.pa, newdata=dat.pa[i, 1:p], type = "survival", times = 5*tau/15)$.pred
    exp.survival.pa$t6[subid] = predict(exp.fit.pa, newdata=dat.pa[i, 1:p], type = "survival", times = 6*tau/15)$.pred
    exp.survival.pa$t7[subid] = predict(exp.fit.pa, newdata=dat.pa[i, 1:p], type = "survival", times = 7*tau/15)$.pred
    exp.survival.pa$t8[subid] = predict(exp.fit.pa, newdata=dat.pa[i, 1:p], type = "survival", times = 8*tau/15)$.pred
    exp.survival.pa$t9[subid] = predict(exp.fit.pa, newdata=dat.pa[i, 1:p], type = "survival", times = 9*tau/15)$.pred
    exp.survival.pa$t10[subid] = predict(exp.fit.pa, newdata=dat.pa[i, 1:p], type = "survival", times = 10*tau/15)$.pred
    exp.survival.pa$t11[subid] = predict(exp.fit.pa, newdata=dat.pa[i, 1:p], type = "survival", times = 11*tau/15)$.pred
    exp.survival.pa$t12[subid] = predict(exp.fit.pa, newdata=dat.pa[i, 1:p], type = "survival", times = 12*tau/15)$.pred
    exp.survival.pa$t13[subid] = predict(exp.fit.pa, newdata=dat.pa[i, 1:p], type = "survival", times = 13*tau/15)$.pred
    exp.survival.pa$t14[subid] = predict(exp.fit.pa, newdata=dat.pa[i, 1:p], type = "survival", times = 14*tau/15)$.pred
  }
  
  exp.pa <- sum((as.matrix(exp.survival.pa[,-1]) - as.matrix(true.survival[,-1]))^2)/nrow(dat.pa)/14
  exp.pa.cum[m] = exp.pa
  
  # 5.2 Gamma model #
  
  gam.survival.pa = data.frame(subid = rep(0, nrow(dat.pa)), t1 = rep(0, nrow(dat.pa)), t2 = rep(0, nrow(dat.pa)), 
                               t3 = rep(0, nrow(dat.pa)), t4 = rep(0, nrow(dat.pa)), t5 = rep(0, nrow(dat.pa)), 
                               t6 = rep(0, nrow(dat.pa)), t7 = rep(0, nrow(dat.pa)), t8 = rep(0, nrow(dat.pa)), 
                               t9 = rep(0, nrow(dat.pa)), t10 = rep(0, nrow(dat.pa)), t11 = rep(0, nrow(dat.pa)),
                               t12 = rep(0, nrow(dat.pa)), t13 = rep(0, nrow(dat.pa)), t14 = rep(0, nrow(dat.pa)))
  gam.fit.pa = flexsurvreg(Surv(time, status) ~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10, data = dat, dist="gamma")
  for (i in 1:nrow(dat.pa)){
    subid = as.numeric(dat.pa$subid[i])
    gam.survival.pa$subid[subid] = subid
    gam.survival.pa$t1[subid] = predict(gam.fit.pa, newdata=dat.pa[i, 1:p], type = "survival", times = 1*tau/15)$.pred
    gam.survival.pa$t2[subid] = predict(gam.fit.pa, newdata=dat.pa[i, 1:p], type = "survival", times = 2*tau/15)$.pred
    gam.survival.pa$t3[subid] = predict(gam.fit.pa, newdata=dat.pa[i, 1:p], type = "survival", times = 3*tau/15)$.pred
    gam.survival.pa$t4[subid] = predict(gam.fit.pa, newdata=dat.pa[i, 1:p], type = "survival", times = 4*tau/15)$.pred
    gam.survival.pa$t5[subid] = predict(gam.fit.pa, newdata=dat.pa[i, 1:p], type = "survival", times = 5*tau/15)$.pred
    gam.survival.pa$t6[subid] = predict(gam.fit.pa, newdata=dat.pa[i, 1:p], type = "survival", times = 6*tau/15)$.pred
    gam.survival.pa$t7[subid] = predict(gam.fit.pa, newdata=dat.pa[i, 1:p], type = "survival", times = 7*tau/15)$.pred
    gam.survival.pa$t8[subid] = predict(gam.fit.pa, newdata=dat.pa[i, 1:p], type = "survival", times = 8*tau/15)$.pred
    gam.survival.pa$t9[subid] = predict(gam.fit.pa, newdata=dat.pa[i, 1:p], type = "survival", times = 9*tau/15)$.pred
    gam.survival.pa$t10[subid] = predict(gam.fit.pa, newdata=dat.pa[i, 1:p], type = "survival", times = 10*tau/15)$.pred
    gam.survival.pa$t11[subid] = predict(gam.fit.pa, newdata=dat.pa[i, 1:p], type = "survival", times = 11*tau/15)$.pred
    gam.survival.pa$t12[subid] = predict(gam.fit.pa, newdata=dat.pa[i, 1:p], type = "survival", times = 12*tau/15)$.pred
    gam.survival.pa$t13[subid] = predict(gam.fit.pa, newdata=dat.pa[i, 1:p], type = "survival", times = 13*tau/15)$.pred
    gam.survival.pa$t14[subid] = predict(gam.fit.pa, newdata=dat.pa[i, 1:p], type = "survival", times = 14*tau/15)$.pred
  }
  
  gam.pa <- sum((as.matrix(gam.survival.pa[,-1]) - as.matrix(true.survival[,-1]))^2)/nrow(dat.pa)/14
  gam.pa.cum[m] = gam.pa
  
  # 5.3 Cox model #
  
  cox.survival.pa = data.frame(subid = rep(0, nrow(dat.pa)), t1 = rep(0, nrow(dat.pa)), t2 = rep(0, nrow(dat.pa)), 
                               t3 = rep(0, nrow(dat.pa)), t4 = rep(0, nrow(dat.pa)), t5 = rep(0, nrow(dat.pa)), 
                               t6 = rep(0, nrow(dat.pa)), t7 = rep(0, nrow(dat.pa)), t8 = rep(0, nrow(dat.pa)), 
                               t9 = rep(0, nrow(dat.pa)), t10 = rep(0, nrow(dat.pa)), t11 = rep(0, nrow(dat.pa)),
                               t12 = rep(0, nrow(dat.pa)), t13 = rep(0, nrow(dat.pa)), t14 = rep(0, nrow(dat.pa)))
  cox.fit.pa = cph(Surv(time, status) ~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10 , data = dat, surv=T)
  for (i in 1:nrow(dat.pa)){
    subid = as.numeric(dat.pa$subid[i])
    cox.survival.pa$subid[subid] = subid
    cox.survival.pa$t1[subid] = predictSurvProb(cox.fit.pa, newdata=dat.pa[i, 1:p], time = 1*tau/15)
    cox.survival.pa$t2[subid] = predictSurvProb(cox.fit.pa, newdata=dat.pa[i, 1:p], time = 2*tau/15)
    cox.survival.pa$t3[subid] = predictSurvProb(cox.fit.pa, newdata=dat.pa[i, 1:p], time = 3*tau/15)
    cox.survival.pa$t4[subid] = predictSurvProb(cox.fit.pa, newdata=dat.pa[i, 1:p], time = 4*tau/15)
    cox.survival.pa$t5[subid] = predictSurvProb(cox.fit.pa, newdata=dat.pa[i, 1:p], time = 5*tau/15)
    cox.survival.pa$t6[subid] = predictSurvProb(cox.fit.pa, newdata=dat.pa[i, 1:p], time = 6*tau/15)
    cox.survival.pa$t7[subid] = predictSurvProb(cox.fit.pa, newdata=dat.pa[i, 1:p], time = 7*tau/15)
    cox.survival.pa$t8[subid] = predictSurvProb(cox.fit.pa, newdata=dat.pa[i, 1:p], time = 8*tau/15)
    cox.survival.pa$t9[subid] = predictSurvProb(cox.fit.pa, newdata=dat.pa[i, 1:p], time = 9*tau/15)
    cox.survival.pa$t10[subid] = predictSurvProb(cox.fit.pa, newdata=dat.pa[i, 1:p], time = 10*tau/15)
    cox.survival.pa$t11[subid] = predictSurvProb(cox.fit.pa, newdata=dat.pa[i, 1:p], time = 11*tau/15)
    cox.survival.pa$t12[subid] = predictSurvProb(cox.fit.pa, newdata=dat.pa[i, 1:p], time = 12*tau/15)
    cox.survival.pa$t13[subid] = predictSurvProb(cox.fit.pa, newdata=dat.pa[i, 1:p], time = 13*tau/15)
    cox.survival.pa$t14[subid] = predictSurvProb(cox.fit.pa, newdata=dat.pa[i, 1:p], time = 14*tau/15)
  }
  
  cox.pa <- sum((as.matrix(cox.survival.pa[,-1]) - as.matrix(true.survival[,-1]))^2)/nrow(dat.pa)/14
  cox.pa.cum[m] = cox.pa
  
  # 5.5 random survival forest #
  
  srf.survival.pa = data.frame(subid = rep(0, nrow(dat.pa)), t1 = rep(0, nrow(dat.pa)), t2 = rep(0, nrow(dat.pa)), 
                          t3 = rep(0, nrow(dat.pa)), t4 = rep(0, nrow(dat.pa)), t5 = rep(0, nrow(dat.pa)), 
                          t6 = rep(0, nrow(dat.pa)), t7 = rep(0, nrow(dat.pa)), t8 = rep(0, nrow(dat.pa)), 
                          t9 = rep(0, nrow(dat.pa)), t10 = rep(0, nrow(dat.pa)), t11 = rep(0, nrow(dat.pa)),
                          t12 = rep(0, nrow(dat.pa)), t13 = rep(0, nrow(dat.pa)), t14 = rep(0, nrow(dat.pa)))
  srf.fit <-  rfsrc(Surv(time, status) ~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10, data = dat)
    for (i in 1:nrow(dat.pa)){
      srf.pred = data.frame(predict(srf.fit, newdata = dat.pa[i,1:p])$survival[1,], predict(srf.fit, newdata = dat.pa[i,1:p])$time.interest)
      subid = as.numeric(dat.pa$subid[i])
      srf.pred$time1 = as.numeric(srf.pred[,2]-1*tau/15 <0)
      srf.pred$time2 = as.numeric(srf.pred[,2]-2*tau/15 <0)
      srf.pred$time3 = as.numeric(srf.pred[,2]-3*tau/15 <0)
      srf.pred$time4 = as.numeric(srf.pred[,2]-4*tau/15 <0)
      srf.pred$time5 = as.numeric(srf.pred[,2]-5*tau/15 <0)
      srf.pred$time6 = as.numeric(srf.pred[,2]-6*tau/15 <0)
      srf.pred$time7 = as.numeric(srf.pred[,2]-7*tau/15 <0)
      srf.pred$time8 = as.numeric(srf.pred[,2]-8*tau/15 <0)
      srf.pred$time9 = as.numeric(srf.pred[,2]-9*tau/15 <0)
      srf.pred$time10 = as.numeric(srf.pred[,2]-10*tau/15 <0)
      srf.pred$time11 = as.numeric(srf.pred[,2]-11*tau/15 <0)
      srf.pred$time12 = as.numeric(srf.pred[,2]-12*tau/15 <0)
      srf.pred$time13 = as.numeric(srf.pred[,2]-13*tau/15 <0)
      srf.pred$time14 = as.numeric(srf.pred[,2]-14*tau/15 <0)
      srf.survival.pa$subid[subid] = subid
      srf.survival.pa$t1[subid] = min(srf.pred[srf.pred$time1>0, 1])
      srf.survival.pa$t2[subid] = min(srf.pred[srf.pred$time2>0, 1])
      srf.survival.pa$t3[subid] = min(srf.pred[srf.pred$time3>0, 1])
      srf.survival.pa$t4[subid] = min(srf.pred[srf.pred$time4>0, 1])
      srf.survival.pa$t5[subid] = min(srf.pred[srf.pred$time5>0, 1])
      srf.survival.pa$t6[subid] = min(srf.pred[srf.pred$time6>0, 1])
      srf.survival.pa$t7[subid] = min(srf.pred[srf.pred$time7>0, 1])
      srf.survival.pa$t8[subid] = min(srf.pred[srf.pred$time8>0, 1])
      srf.survival.pa$t9[subid] = min(srf.pred[srf.pred$time9>0, 1])
      srf.survival.pa$t10[subid] = min(srf.pred[srf.pred$time10>0, 1])
      srf.survival.pa$t11[subid] = min(srf.pred[srf.pred$time11>0, 1])
      srf.survival.pa$t12[subid] = min(srf.pred[srf.pred$time12>0, 1])
      srf.survival.pa$t13[subid] = min(srf.pred[srf.pred$time13>0, 1])
      srf.survival.pa$t14[subid] = min(srf.pred[srf.pred$time14>0, 1])
  }
  
  srf.pa <- sum((as.matrix(srf.survival.pa[,-1]) - as.matrix(true.survival[,-1]))^2)/nrow(dat.pa)/14
  srf.pa.cum[m] = srf.pa
  
  # 5.6 IPCW super learner #
  
  ipcwsl.survival.pa = ipcw.wt[1]*exp.survival.pa[,-1] + ipcw.wt[2]*gam.survival.pa[,-1] + ipcw.wt[3]*cox.survival.pa[,-1] + ipcw.wt[4]*srf.survival.pa[,-1]
  ipcwsl.pa <- sum((as.matrix(ipcwsl.survival.pa) - as.matrix(true.survival[,-1]))^2)/nrow(dat.pa)/14
  ipcwsl.pa.cum[m] = ipcwsl.pa
  })
}


write.csv(ipcw.wt.cum, file=paste0(path.results, 'c_large_ipcw_wt_', model, '_', spec, '_', cen.mod, '_', batch, '.csv'), row.names = F)
write.csv(ipcwsl.pa.cum, file=paste0(path.results, 'c_large_ipcwsl_pa_', model, '_', spec, '_', cen.mod, '_', batch, '.csv'), row.names = F)

write.csv(exp.pa.cum, file=paste0(path.results, 'c_large_exp_pa_', model, '_', spec, '_', cen.mod, '_', batch, '.csv'), row.names = F)
write.csv(gam.pa.cum, file=paste0(path.results, 'c_large_gam_pa_', model, '_', spec, '_', cen.mod, '_', batch, '.csv'), row.names = F)
write.csv(cox.pa.cum, file=paste0(path.results, 'c_large_cox_pa_', model, '_', spec, '_', cen.mod, '_', batch, '.csv'), row.names = F)
write.csv(srf.pa.cum, file=paste0(path.results, 'c_large_srf_pa_', model, '_', spec, '_', cen.mod, '_', batch, '.csv'), row.names = F)

mean(ipcwsl.pa.cum[1:nsim])
#mean(exp.pa.cum[1:nsim])
#mean(gam.pa.cum[1:nsim])
#mean(cox.pa.cum[1:nsim])
#mean(srf.pa.cum[1:nsim])
#mean(exp2.pa.cum[1:nsim])
#mean(gam2.pa.cum[1:nsim])
#mean(cox2.pa.cum[1:nsim])
#mean(srf2.pa.cum[1:nsim])


median(ipcwsl.pa.cum[1:nsim])
#median(exp.pa.cum[1:nsim])
#median(gam.pa.cum[1:nsim])
#median(cox.pa.cum[1:nsim])
#median(srf.pa.cum[1:nsim])
#median(exp2.pa.cum[1:nsim])
#median(gam2.pa.cum[1:nsim])
#median(cox2.pa.cum[1:nsim])
#median(srf2.pa.cum[1:nsim])


         








