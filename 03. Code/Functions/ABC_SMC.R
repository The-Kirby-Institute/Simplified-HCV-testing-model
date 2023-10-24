# pre ABC SMC 
calc_distance <- function(model_dt, calib_dt){
  sre1 <- sqrt(sum(
    (c(as.numeric(unlist((model_dt[["N"]])))) - 
                       c(as.numeric(unlist(calib_dt[["N"]][ ,3]))))^2))
    

  
  sre2 <- sqrt(sum((c(unlist(model_dt[["frac"]])) - 
                      calib_dt[["frac"]]$realPop)^2))

  return(c(sre1, sre2))
}


# Perturbation kernel 
rK <- function(mean, sigma){   
  return(rtmvnorm(1,mean=mean, sigma=sigma, lower=lm.low, upper=lm.upp)) 
}

#  Identity function: H(x)= 1 if x=T
H <- function(x) as.numeric(x>0)

#  Test if prior is non zero
prior.non.zero<-function(par){
  prod(sapply(1:14, function(a) H(par[a]-lm.low[a])* H(lm.upp[a]-par[a])), na.rm = TRUE)
}

Norm.Eucl.dist<-function(p1,p2){
  sqrt(sum(((p1-p2)/(lm.upp-lm.low))^2)) }

#  Covariance based on M neighbours
getSigmaNeighbours<-function(M, theta, Theta){
  dist<- sapply(1:N, function(a) Norm.Eucl.dist(as.numeric(theta), as.numeric(Theta[a,])))
  temp<-data.frame(no=seq(1,N), dist)
  temp<-temp[order(temp$dist),]
  sigma<-cov(Theta[temp$no[1:(M+1)],])
  return(sigma)
}