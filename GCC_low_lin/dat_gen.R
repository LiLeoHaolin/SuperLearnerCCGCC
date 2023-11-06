## Prerequisites ##

library(MASS)

## Simulation setting ##

#path.test <- "C:/Users/haoli/Desktop/Dissertation/04_project1/03_simulations/02_data_generation/RC/test/test/"
#path.c <- "C:/Users/haoli/Desktop/Dissertation/04_project1/03_simulations/02_data_generation/RC/test/c/"
#path.cc <- "C:/Users/haoli/Desktop/Dissertation/04_project1/03_simulations/02_data_generation/RC/test/cc/"
#path.sc <- "C:/Users/haoli/Desktop/Dissertation/04_project1/03_simulations/02_data_generation/RC/test/sc/"
#path.r <- "C:/Users/haoli/Desktop/Dissertation/04_project1/03_simulations/02_data_generation/RC/test/r/"

path.test <- "/nas/longleaf/home/haolin/dissertation1/round1/data/test/"
path.c <- "/nas/longleaf/home/haolin/dissertation1/round1/data/c/"
path.cc <- "/nas/longleaf/home/haolin/dissertation1/round1/data/cc/"
path.sc <- "/nas/longleaf/home/haolin/dissertation1/round1/data/sc/"
path.r <- "/nas/longleaf/home/haolin/dissertation1/round1/data/r/"

gamma.weibull = 0.7
p = 10
p.full = 6
n = 1500
q = 0.1
q.comp = 0.154
q.case = 0.3
nsim = 500

cenc = rep(0, nsim)

## Data generation (training) ##

for (j in 1:nsim){
  cat(j)
  
  ar1_cor <- function(n, rho) {
    exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                      (1:n - 1))
    rho^exponent
  }
  
  
  X = as.matrix(mvrnorm(n = 2*n, matrix(0, nrow = p, ncol = 1), ar1_cor(p, 0.6)))
  
  beta = matrix(c(0.1, 0.1, 0.3, -0.3, -0.2, -0.3, 0.01, -0.01, 0.25, -0.25), nrow = p)
  
  true.lin.pred <- X %*% beta +1.5 #+ 0.3*X[,1]*X[,3] + 0.2*X[,3]^2 + 0.2*X[,1]*X[,6]*X[,9]
  v <- runif(n=2*n)
  Tlat <- (- log(v) / (1* exp(true.lin.pred)))^(1 / gamma.weibull)
 
  C = runif(2*n, min = 0, max = 0.03) # censoring times

  end <- rep(0.0162, 2*n)
  time <- pmin(Tlat, C, end) # follow-up times
  
  status <- as.numeric(Tlat == time) #event indicators
  
  
  dat.pre = data.frame(X, status, time)
  
  dat = dat.pre[1:n,]
  test = dat.pre[(n+1):(2*n),]

  dat$subid = 1:n
  test$subid = 1:n
  
  dat$sc = sample(c(1:0), n, replace=T, prob = c(q,1-q))
  dat$r = sample(c(1:0), n, replace=T, prob = c(q.comp,1-q.comp))
 
  cenc[j] = table(dat$status)[1]/n
  
  c <- subset(dat, select = -c(sc, r))
  write.csv(c,file=paste0(path.c, 'dat', j, '.csv'), row.names = F)
  
  sc <- subset(dat, select = -c(r))
  for (m in 1:n){
    if (sc$sc[m]==1){
      sc[m,c((p.full+1):p)] = sc[m,c((p.full+1):p)]
    }else{
      sc[m,c((p.full+1):p)] = NA
    }
  }
  sc <- subset(sc, select = -c(sc))
  write.csv(sc,file=paste0(path.sc, 'dat', j, '.csv'), row.names = F)
  
  r <- subset(dat, select = -c(sc))
  for (m in 1:n){
    if (r$r[m]==1){
      r[m,c((p.full+1):p)] = r[m,c((p.full+1):p)]
    }else{
      r[m,c((p.full+1):p)] = NA
    }
  }
  r <- subset(r, select = -c(r))
  write.csv(r,file=paste0(path.r, 'dat', j, '.csv'), row.names = F)
  
  ntilde = sum(dat$sc)
  
  dat$sc.case = sample(c(1:0), n, replace=T, prob = c(q.case,1-q.case))
  
  cc <- subset(dat, select = -c(r))
  for (m in 1:n){
    if ((dat$sc[m]==1)|((dat$status[m]==1)&(dat$sc.case[m]==1))){
      cc[m,c((p.full+1):p)] = cc[m,c((p.full+1):p)]
    }else{
      cc[m,c((p.full+1):p)] = NA
    }
  }
  cc$weight = 1*(cc$status==1)*(cc$sc==1)+n/ntilde*(cc$status==0)*(cc$sc==1)+sum(dat$status*(1-dat$sc))/sum(dat$status*(1-dat$sc)*(dat$sc.case))*(cc$status==1)*(cc$sc==0)*(cc$sc.case==1)
  cc <- subset(cc, select = -c(sc, sc.case))
  #cc$weight = 1*(cc$status==1)+n/ntilde*(cc$status==0)
  write.csv(cc,file=paste0(path.cc, 'dat', j, '.csv'), row.names = F)
  
  write.csv(test,file=paste0(path.test, 'dat', j, '.csv'), row.names = F)
  
}

mean(cenc)

