## Prerequisites ##

library(MASS)
library(survival)
library(flexsurv)
library(pec)
library(rms)
library(randomForestSRC)
library(alabama)
library(CoxBoost)

## Simulation setting ##

#path.test <- "C:/Users/Haolin Li/Desktop/Dissertation/04_project1/07_simulation_new/data/test/"
#path.c <- "C:/Users/Haolin Li/Desktop/Dissertation/04_project1/07_simulation_new/data/c/"
#path.cc <- "C:/Users/Haolin Li/Desktop/Dissertation/04_project1/07_simulation_new/data/cc/"
#path.sc <- "C:/Users/Haolin Li/Desktop/Dissertation/04_project1/07_simulation_new/data/sc/"
#path.r <- "C:/Users/Haolin Li/Desktop/Dissertation/04_project1/07_simulation_new/data/r/"

path.test <- "/nas/longleaf/home/haolin/dissertation1/round6/data/test/"
#path.c <- "/nas/longleaf/home/haolin/dissertation1/round6/data/c/"
path.cc <- "/nas/longleaf/home/haolin/dissertation1/round6/data/cc/"
#path.sc <- "/nas/longleaf/home/haolin/dissertation1/round6/data/sc/"
path.results <- "/nas/longleaf/home/haolin/dissertation1/round6/results/"

p = 400
p.full = 140
Paste <- function(x, y = "") paste(paste0(x, y), collapse = " + ")
formula.cc = as.formula(paste0('Surv(time, status) ~ ', Paste("X", 1:400)))
formula.full = as.formula(paste0('Surv(time, status) ~ ', Paste("X", 1:140)))
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

cb.pa.cum = rep(0, nsim)
srf.pa.cum = rep(0, nsim)
cb2.pa.cum = rep(0, nsim)
srf2.pa.cum = rep(0, nsim)
ipcwsl.pa.cum = rep(0, nsim)


for (m in seq){
  cat(m)
  try({
  # 1. Prepare data for k-fold CV #
  
  dat = read.csv(paste0(path.cc, 'dat', m, '.csv'))
  n = nrow(dat)
  dat$subid = c(1:n)
  dat$fold = 0
  for (i in 1:n){
    dat$fold[i] = sample(c(1:fold), 1)
  }
  ncc = sum(is.na(dat$X400)==F)
  
  # 2. Model fitting of candidate learners and prediction for k-fold CV #
  
  # 2.1 CoxBoost model #
  
  cb.survival.pred = data.frame(subid = rep(NA, n), t1 = rep(NA, n), t2 = rep(NA, n), t3 = rep(NA, n), t4 = rep(NA, n),
                                t5 = rep(NA, n), t6 = rep(NA, n), t7 = rep(NA, n), t8 = rep(NA, n), t9 = rep(NA, n))
  for (j in 1:fold){
    dat.train = dat[dat$fold!=j, ]
    dat.train = dat.train[is.na(dat.train$X400)==F,]
    dat.test = dat[dat$fold==j, ]
    dat.test = dat.test[is.na(dat.test$X400)==F,]
    cb.fit = CoxBoost(time=dat.train$time, status=dat.train$status, x=as.matrix(dat.train[,1:p]), weights = dat.train$weight) 
    for (i in 1:nrow(dat.test)){
      subid = as.numeric(dat.test$subid[i])
      cb.survival.pred$subid[subid] = subid
      cb.survival.pred$t1[subid] = 1-predict(cb.fit, newdata = as.matrix(dat.test[i, 1:p]),  time = 1*tau/10, type = "CIF")
      cb.survival.pred$t2[subid] = 1-predict(cb.fit, newdata = as.matrix(dat.test[i, 1:p]),  time = 2*tau/10, type = "CIF")
      cb.survival.pred$t3[subid] = 1-predict(cb.fit, newdata = as.matrix(dat.test[i, 1:p]),  time = 3*tau/10, type = "CIF")
      cb.survival.pred$t4[subid] = 1-predict(cb.fit, newdata = as.matrix(dat.test[i, 1:p]),  time = 4*tau/10, type = "CIF")
      cb.survival.pred$t5[subid] = 1-predict(cb.fit, newdata = as.matrix(dat.test[i, 1:p]),  time = 5*tau/10, type = "CIF")
      cb.survival.pred$t6[subid] = 1-predict(cb.fit, newdata = as.matrix(dat.test[i, 1:p]),  time = 6*tau/10, type = "CIF")
      cb.survival.pred$t7[subid] = 1-predict(cb.fit, newdata = as.matrix(dat.test[i, 1:p]),  time = 7*tau/10, type = "CIF")
      cb.survival.pred$t8[subid] = 1-predict(cb.fit, newdata = as.matrix(dat.test[i, 1:p]),  time = 8*tau/10, type = "CIF")
      cb.survival.pred$t9[subid] = 1-predict(cb.fit, newdata = as.matrix(dat.test[i, 1:p]),  time = 9*tau/10, type = "CIF")
    }
  }
  
  # 2.4 Survival random forest #
  
  srf.survival.pred = data.frame(subid = rep(NA, n), t1 = rep(NA, n), t2 = rep(NA, n), t3 = rep(NA, n), t4 = rep(NA, n),
                            t5 = rep(NA, n), t6 = rep(NA, n), t7 = rep(NA, n), t8 = rep(NA, n), t9 = rep(NA, n))
  for (j in 1:fold){
    dat.train = dat[dat$fold!=j, ]
    dat.train = dat.train[is.na(dat.train$X400)==F,]
    dat.test = dat[dat$fold==j, ]
    dat.test = dat.test[is.na(dat.test$X400)==F,]
    srf.fit = rfsrc(formula.cc, case.wt = dat.train$weight, data = dat.train)
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
  
  # 2.5 CoxBoost model 2 #
  
  cb2.survival.pred = data.frame(subid = rep(NA, n), t1 = rep(NA, n), t2 = rep(NA, n), t3 = rep(NA, n), t4 = rep(NA, n),
                                t5 = rep(NA, n), t6 = rep(NA, n), t7 = rep(NA, n), t8 = rep(NA, n), t9 = rep(NA, n))
  for (j in 1:fold){
    dat.train = dat[dat$fold!=j, ]
    dat.test = dat[dat$fold==j, ]
    dat.test = dat.test[is.na(dat.test$X400)==F,]
    cb.fit = CoxBoost(time=dat.train$time, status=dat.train$status, x=as.matrix(dat.train[,1:p.full])) 
    for (i in 1:nrow(dat.test)){
      subid = as.numeric(dat.test$subid[i])
      cb2.survival.pred$subid[subid] = subid
      cb2.survival.pred$t1[subid] = 1-predict(cb.fit, newdata = as.matrix(dat.test[i, 1:p.full]),  time = 1*tau/10, type = "CIF")
      cb2.survival.pred$t2[subid] = 1-predict(cb.fit, newdata = as.matrix(dat.test[i, 1:p.full]),  time = 2*tau/10, type = "CIF")
      cb2.survival.pred$t3[subid] = 1-predict(cb.fit, newdata = as.matrix(dat.test[i, 1:p.full]),  time = 3*tau/10, type = "CIF")
      cb2.survival.pred$t4[subid] = 1-predict(cb.fit, newdata = as.matrix(dat.test[i, 1:p.full]),  time = 4*tau/10, type = "CIF")
      cb2.survival.pred$t5[subid] = 1-predict(cb.fit, newdata = as.matrix(dat.test[i, 1:p.full]),  time = 5*tau/10, type = "CIF")
      cb2.survival.pred$t6[subid] = 1-predict(cb.fit, newdata = as.matrix(dat.test[i, 1:p.full]),  time = 6*tau/10, type = "CIF")
      cb2.survival.pred$t7[subid] = 1-predict(cb.fit, newdata = as.matrix(dat.test[i, 1:p.full]),  time = 7*tau/10, type = "CIF")
      cb2.survival.pred$t8[subid] = 1-predict(cb.fit, newdata = as.matrix(dat.test[i, 1:p.full]),  time = 8*tau/10, type = "CIF")
      cb2.survival.pred$t9[subid] = 1-predict(cb.fit, newdata = as.matrix(dat.test[i, 1:p.full]),  time = 9*tau/10, type = "CIF")
    }
  }
  
  # 2.8 Survival random forest 2 #
  
  srf2.survival.pred = data.frame(subid = rep(NA, n), t1 = rep(NA, n), t2 = rep(NA, n), t3 = rep(NA, n), t4 = rep(NA, n),
                                 t5 = rep(NA, n), t6 = rep(NA, n), t7 = rep(NA, n), t8 = rep(NA, n), t9 = rep(NA, n))
  for (j in 1:fold){
    dat.train = dat[dat$fold!=j, ]
    dat.test = dat[dat$fold==j, ]
    dat.test = dat.test[is.na(dat.test$X400)==F,]
    srf.fit = rfsrc(formula.full, data = dat.train)
    for (i in 1:nrow(dat.test)){
      srf.pred = data.frame(predict(srf.fit, newdata = dat.test[i,1:p.full])$survival[1,], predict(srf.fit, newdata = dat.test[i,1:p.full])$time.interest)
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
      srf2.survival.pred$subid[subid] = subid
      srf2.survival.pred$t1[subid] = min(srf.pred[srf.pred$time1>0, 1])
      srf2.survival.pred$t2[subid] = min(srf.pred[srf.pred$time2>0, 1])
      srf2.survival.pred$t3[subid] = min(srf.pred[srf.pred$time3>0, 1])
      srf2.survival.pred$t4[subid] = min(srf.pred[srf.pred$time4>0, 1])
      srf2.survival.pred$t5[subid] = min(srf.pred[srf.pred$time5>0, 1])
      srf2.survival.pred$t6[subid] = min(srf.pred[srf.pred$time6>0, 1])
      srf2.survival.pred$t7[subid] = min(srf.pred[srf.pred$time7>0, 1])
      srf2.survival.pred$t8[subid] = min(srf.pred[srf.pred$time8>0, 1])
      srf2.survival.pred$t9[subid] = min(srf.pred[srf.pred$time9>0, 1])
    }
  }
  
  # 3. IPCW Super Learner #
  
  # 3.0 weights #
  wt = dat
  wt$t1 = wt$weight
  wt$t2 = wt$weight
  wt$t3 = wt$weight
  wt$t4 = wt$weight
  wt$t5 = wt$weight
  wt$t6 = wt$weight
  wt$t7 = wt$weight
  wt$t8 = wt$weight
  wt$t9 = wt$weight
  wt = subset(wt, select = c(subid, t1, t2, t3, t4, t5, t6, t7, t8, t9))
  
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
  
  idx = na.omit(cb.survival.pred$subid)
  
  cb.survival.pred = cb.survival.pred[idx,]
  srf.survival.pred = srf.survival.pred[idx,]
  cb2.survival.pred = cb2.survival.pred[idx,]
  srf2.survival.pred = srf2.survival.pred[idx,]
  
  ind = ind[idx,]
  frac = frac[idx,]
  wt = wt[idx, ]
  
  fn <- function(x){
    fn = sum((as.matrix(ind[,-1])-(x[1]*as.matrix(cb.survival.pred[,-1]) + x[2]*as.matrix(srf.survival.pred[,-1]) + 
                                     x[3]*as.matrix(cb2.survival.pred[,-1]) + x[4]*as.matrix(srf2.survival.pred[,-1]) ))^2 
             * (as.matrix(frac[,-1]))) / (n)
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
  
  beta = matrix(c(rep(0,p.full-2), c(0.1, 0.1, 0.3, -0.3, -0.2, -0.3, 0.01, -0.01, 0.25, -0.25), rep(0,p-p.full-8)), nrow=p)
  true.lin.pred <- as.matrix(dat.pa[,1:p]) %*% beta + 0.3*as.matrix(dat.pa[,139])*as.matrix(dat.pa[,141]) - 0.3*as.matrix(dat.pa[,141])^2 + 0.2*as.matrix(dat.pa[,139])*as.matrix(dat.pa[,144])*as.matrix(dat.pa[,147]) - 0.4*as.matrix(dat.pa[,140])*as.matrix(dat.pa[,142])*as.matrix(dat.pa[,143])+1.6#+ 0.3*X[,1]*X[,3] + 0.2*X[,3]^2 + 0.2*X[,1]*X[,6]*X[,9]
  
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
  
  # 5.0.1 Prepare CC data #
  
  dat.cc = dat[is.na(dat$X400)==F,]
  
  # 5.3 Cox model #
  
  cb.survival.pa = data.frame(subid = rep(NA, nrow(dat.pa)), t1 = rep(NA, nrow(dat.pa)), t2 = rep(NA, nrow(dat.pa)), 
                              t3 = rep(NA, nrow(dat.pa)), t4 = rep(NA, nrow(dat.pa)), t5 = rep(NA, nrow(dat.pa)), 
                              t6 = rep(NA, nrow(dat.pa)), t7 = rep(NA, nrow(dat.pa)), t8 = rep(NA, nrow(dat.pa)), 
                              t9 = rep(NA, nrow(dat.pa)), t10 = rep(NA, nrow(dat.pa)), t11 = rep(NA, nrow(dat.pa)),
                              t12 = rep(NA, nrow(dat.pa)), t13 = rep(NA, nrow(dat.pa)), t14 = rep(NA, nrow(dat.pa)))
  cb.fit.pa = CoxBoost(time=dat.cc$time, status=dat.cc$status, x=as.matrix(dat.cc[,1:p]), weights = dat.cc$weight) 
  for (i in 1:nrow(dat.pa)){
    subid = as.numeric(dat.pa$subid[i])
    cb.survival.pa$subid[subid] = subid
    cb.survival.pa$t1[subid] = 1-predict(cb.fit.pa, newdata = as.matrix(dat.pa[i, 1:p]),  time = 1*tau/15, type = "CIF")
    cb.survival.pa$t2[subid] = 1-predict(cb.fit.pa, newdata = as.matrix(dat.pa[i, 1:p]),  time = 2*tau/15, type = "CIF")
    cb.survival.pa$t3[subid] = 1-predict(cb.fit.pa, newdata = as.matrix(dat.pa[i, 1:p]),  time = 3*tau/15, type = "CIF")
    cb.survival.pa$t4[subid] = 1-predict(cb.fit.pa, newdata = as.matrix(dat.pa[i, 1:p]),  time = 4*tau/15, type = "CIF")
    cb.survival.pa$t5[subid] = 1-predict(cb.fit.pa, newdata = as.matrix(dat.pa[i, 1:p]),  time = 5*tau/15, type = "CIF")
    cb.survival.pa$t6[subid] = 1-predict(cb.fit.pa, newdata = as.matrix(dat.pa[i, 1:p]),  time = 6*tau/15, type = "CIF")
    cb.survival.pa$t7[subid] = 1-predict(cb.fit.pa, newdata = as.matrix(dat.pa[i, 1:p]),  time = 7*tau/15, type = "CIF")
    cb.survival.pa$t8[subid] = 1-predict(cb.fit.pa, newdata = as.matrix(dat.pa[i, 1:p]),  time = 8*tau/15, type = "CIF")
    cb.survival.pa$t9[subid] = 1-predict(cb.fit.pa, newdata = as.matrix(dat.pa[i, 1:p]),  time = 9*tau/15, type = "CIF")
    cb.survival.pa$t10[subid] = 1-predict(cb.fit.pa, newdata = as.matrix(dat.pa[i, 1:p]),  time = 10*tau/15, type = "CIF")
    cb.survival.pa$t11[subid] = 1-predict(cb.fit.pa, newdata = as.matrix(dat.pa[i, 1:p]),  time = 11*tau/15, type = "CIF")
    cb.survival.pa$t12[subid] = 1-predict(cb.fit.pa, newdata = as.matrix(dat.pa[i, 1:p]),  time = 12*tau/15, type = "CIF")
    cb.survival.pa$t13[subid] = 1-predict(cb.fit.pa, newdata = as.matrix(dat.pa[i, 1:p]),  time = 13*tau/15, type = "CIF")
    cb.survival.pa$t14[subid] = 1-predict(cb.fit.pa, newdata = as.matrix(dat.pa[i, 1:p]),  time = 14*tau/15, type = "CIF")
  }
  
  cb.pa <- sum((as.matrix(cb.survival.pa[,-1]) - as.matrix(true.survival[,-1]))^2)/nrow(dat.pa)/14
  cb.pa.cum[m] = cb.pa
  
  # 5.5 random survival forest #
  
  srf.survival.pa = data.frame(subid = rep(NA, nrow(dat.pa)), t1 = rep(NA, nrow(dat.pa)), t2 = rep(NA, nrow(dat.pa)), 
                          t3 = rep(NA, nrow(dat.pa)), t4 = rep(NA, nrow(dat.pa)), t5 = rep(NA, nrow(dat.pa)), 
                          t6 = rep(NA, nrow(dat.pa)), t7 = rep(NA, nrow(dat.pa)), t8 = rep(NA, nrow(dat.pa)), 
                          t9 = rep(NA, nrow(dat.pa)), t10 = rep(NA, nrow(dat.pa)), t11 = rep(NA, nrow(dat.pa)),
                          t12 = rep(NA, nrow(dat.pa)), t13 = rep(NA, nrow(dat.pa)), t14 = rep(NA, nrow(dat.pa)))
  srf.fit <-  rfsrc(formula.cc, case.wt = dat.cc$weight, data = dat.cc)
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
  
  # 5.7 Cox model 2 #
  
  cb2.survival.pa = data.frame(subid = rep(NA, nrow(dat.pa)), t1 = rep(NA, nrow(dat.pa)), t2 = rep(NA, nrow(dat.pa)), 
                              t3 = rep(NA, nrow(dat.pa)), t4 = rep(NA, nrow(dat.pa)), t5 = rep(NA, nrow(dat.pa)), 
                              t6 = rep(NA, nrow(dat.pa)), t7 = rep(NA, nrow(dat.pa)), t8 = rep(NA, nrow(dat.pa)), 
                              t9 = rep(NA, nrow(dat.pa)), t10 = rep(NA, nrow(dat.pa)), t11 = rep(NA, nrow(dat.pa)),
                              t12 = rep(NA, nrow(dat.pa)), t13 = rep(NA, nrow(dat.pa)), t14 = rep(NA, nrow(dat.pa)))
  cb.fit.pa = CoxBoost(time=dat$time, status=dat$status, x=as.matrix(dat[,1:p.full])) 
  for (i in 1:nrow(dat.pa)){
    subid = as.numeric(dat.pa$subid[i])
    cb2.survival.pa$subid[subid] = subid
    cb2.survival.pa$t1[subid] = 1-predict(cb.fit.pa, newdata = as.matrix(dat.pa[i, 1:p.full]),  time = 1*tau/15, type = "CIF")
    cb2.survival.pa$t2[subid] = 1-predict(cb.fit.pa, newdata = as.matrix(dat.pa[i, 1:p.full]),  time = 2*tau/15, type = "CIF")
    cb2.survival.pa$t3[subid] = 1-predict(cb.fit.pa, newdata = as.matrix(dat.pa[i, 1:p.full]),  time = 3*tau/15, type = "CIF")
    cb2.survival.pa$t4[subid] = 1-predict(cb.fit.pa, newdata = as.matrix(dat.pa[i, 1:p.full]),  time = 4*tau/15, type = "CIF")
    cb2.survival.pa$t5[subid] = 1-predict(cb.fit.pa, newdata = as.matrix(dat.pa[i, 1:p.full]),  time = 5*tau/15, type = "CIF")
    cb2.survival.pa$t6[subid] = 1-predict(cb.fit.pa, newdata = as.matrix(dat.pa[i, 1:p.full]),  time = 6*tau/15, type = "CIF")
    cb2.survival.pa$t7[subid] = 1-predict(cb.fit.pa, newdata = as.matrix(dat.pa[i, 1:p.full]),  time = 7*tau/15, type = "CIF")
    cb2.survival.pa$t8[subid] = 1-predict(cb.fit.pa, newdata = as.matrix(dat.pa[i, 1:p.full]),  time = 8*tau/15, type = "CIF")
    cb2.survival.pa$t9[subid] = 1-predict(cb.fit.pa, newdata = as.matrix(dat.pa[i, 1:p.full]),  time = 9*tau/15, type = "CIF")
    cb2.survival.pa$t10[subid] = 1-predict(cb.fit.pa, newdata = as.matrix(dat.pa[i, 1:p.full]),  time = 10*tau/15, type = "CIF")
    cb2.survival.pa$t11[subid] = 1-predict(cb.fit.pa, newdata = as.matrix(dat.pa[i, 1:p.full]),  time = 11*tau/15, type = "CIF")
    cb2.survival.pa$t12[subid] = 1-predict(cb.fit.pa, newdata = as.matrix(dat.pa[i, 1:p.full]),  time = 12*tau/15, type = "CIF")
    cb2.survival.pa$t13[subid] = 1-predict(cb.fit.pa, newdata = as.matrix(dat.pa[i, 1:p.full]),  time = 13*tau/15, type = "CIF")
    cb2.survival.pa$t14[subid] = 1-predict(cb.fit.pa, newdata = as.matrix(dat.pa[i, 1:p.full]),  time = 14*tau/15, type = "CIF")
  }
  
  cb2.pa <- sum((as.matrix(cb.survival.pa[,-1]) - as.matrix(true.survival[,-1]))^2)/nrow(dat.pa)/14
  cb2.pa.cum[m] = cb.pa
  
  # 5.8 random survival forest 2 #
  
  srf2.survival.pa = data.frame(subid = rep(NA, nrow(dat.pa)), t1 = rep(NA, nrow(dat.pa)), t2 = rep(NA, nrow(dat.pa)), 
                               t3 = rep(NA, nrow(dat.pa)), t4 = rep(NA, nrow(dat.pa)), t5 = rep(NA, nrow(dat.pa)), 
                               t6 = rep(NA, nrow(dat.pa)), t7 = rep(NA, nrow(dat.pa)), t8 = rep(NA, nrow(dat.pa)), 
                               t9 = rep(NA, nrow(dat.pa)), t10 = rep(NA, nrow(dat.pa)), t11 = rep(NA, nrow(dat.pa)),
                               t12 = rep(NA, nrow(dat.pa)), t13 = rep(NA, nrow(dat.pa)), t14 = rep(NA, nrow(dat.pa)))
  srf.fit <-  rfsrc(formula.full, data = dat)
  for (i in 1:nrow(dat.pa)){
    srf.pred = data.frame(predict(srf.fit, newdata = dat.pa[i,1:p.full])$survival[1,], predict(srf.fit, newdata = dat.pa[i,1:p.full])$time.interest)
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
    srf2.survival.pa$subid[subid] = subid
    srf2.survival.pa$t1[subid] = min(srf.pred[srf.pred$time1>0, 1])
    srf2.survival.pa$t2[subid] = min(srf.pred[srf.pred$time2>0, 1])
    srf2.survival.pa$t3[subid] = min(srf.pred[srf.pred$time3>0, 1])
    srf2.survival.pa$t4[subid] = min(srf.pred[srf.pred$time4>0, 1])
    srf2.survival.pa$t5[subid] = min(srf.pred[srf.pred$time5>0, 1])
    srf2.survival.pa$t6[subid] = min(srf.pred[srf.pred$time6>0, 1])
    srf2.survival.pa$t7[subid] = min(srf.pred[srf.pred$time7>0, 1])
    srf2.survival.pa$t8[subid] = min(srf.pred[srf.pred$time8>0, 1])
    srf2.survival.pa$t9[subid] = min(srf.pred[srf.pred$time9>0, 1])
    srf2.survival.pa$t10[subid] = min(srf.pred[srf.pred$time10>0, 1])
    srf2.survival.pa$t11[subid] = min(srf.pred[srf.pred$time11>0, 1])
    srf2.survival.pa$t12[subid] = min(srf.pred[srf.pred$time12>0, 1])
    srf2.survival.pa$t13[subid] = min(srf.pred[srf.pred$time13>0, 1])
    srf2.survival.pa$t14[subid] = min(srf.pred[srf.pred$time14>0, 1])
  }
  
  srf2.pa <- sum((as.matrix(srf2.survival.pa[,-1]) - as.matrix(true.survival[,-1]))^2)/nrow(dat.pa)/14
  srf2.pa.cum[m] = srf2.pa
  
  # 5.6 IPCW super learner #
  
  ipcwsl.survival.pa = ipcw.wt[1]*cb.survival.pa[,-1] + ipcw.wt[2]*srf.survival.pa[,-1] +
    ipcw.wt[3]*cb2.survival.pa[,-1] + ipcw.wt[4]*srf2.survival.pa[,-1]
  ipcwsl.pa <- sum((as.matrix(ipcwsl.survival.pa) - as.matrix(true.survival[,-1]))^2)/nrow(dat.pa)/14
  ipcwsl.pa.cum[m] = ipcwsl.pa
  })
}


write.csv(ipcw.wt.cum, file=paste0(path.results, 'cc_uw_w_ipcw_wt_', model, '_', spec, '_', cen.mod, '_', batch, '.csv'), row.names = F)
write.csv(ipcwsl.pa.cum, file=paste0(path.results, 'cc_uw_w_ipcwsl_pa_', model, '_', spec, '_', cen.mod, '_', batch, '.csv'), row.names = F)

#write.csv(cb.pa.cum, file=paste0(path.results, 'cc_w_cb_pa_', model, '_', spec, '_', cen.mod, '_', batch, '.csv'), row.names = F)
#write.csv(srf.pa.cum, file=paste0(path.results, 'cc_w_srf_pa_', model, '_', spec, '_', cen.mod, '_', batch, '.csv'), row.names = F)
#write.csv(cb2.pa.cum, file=paste0(path.results, 'cc_w_cb2_pa_', model, '_', spec, '_', cen.mod, '_', batch, '.csv'), row.names = F)
#write.csv(srf2.pa.cum, file=paste0(path.results, 'cc_w_srf2_pa_', model, '_', spec, '_', cen.mod, '_', batch, '.csv'), row.names = F)

mean(ipcwsl.pa.cum[1:nsim])
mean(cb.pa.cum[1:nsim])
mean(srf.pa.cum[1:nsim])
mean(cb2.pa.cum[1:nsim])
mean(srf2.pa.cum[1:nsim])

median(ipcwsl.pa.cum[1:nsim])
median(cb.pa.cum[1:nsim])
median(srf.pa.cum[1:nsim])
median(cb2.pa.cum[1:nsim])
median(srf2.pa.cum[1:nsim])

         








