

#### round2: to round values to specified digits
round2 = function(x, n) {
  posneg = sign(x)
  z = abs(x)*10^n
  z = z + 0.5
  z = trunc(z)
  z = z/10^n
  z*posneg
}


#### rmse: to calculate root mean squared errors as defined in the paper
rmse <- function(theta, est)
  sqrt(mean((theta - est)^2))


#### pval.2side: function for 2-sided p-value calculation 
####   following the generated statistics by the proposed resampling procedure 
pval.2side <- function(stat, stat.bt)
{
  len.ub <- length(which(stat.bt >= stat))
  len.lb <- length(which(stat.bt <= stat))
  
  round2(2*min(len.ub, len.lb)/length(stat.bt), 4)
}

#### PLS1: auxiliary function to fit a general PLS regression
PLS1 <- function(X, y, k) # PLS1, y is a univariate outcome
{
  R <- matrix(NA, nrow = ncol(X), ncol = k)
  V <- matrix(0, nrow = ncol(X), ncol = k)
  P <- matrix(NA, nrow = ncol(X), ncol = k)
  T <- matrix(NA, nrow = nrow(X), ncol = k)
  
  S <- cov(X,y)
  
  for (a in 1:k)
  {
    r <- S
    r <- r/sqrt(c(crossprod(r))) # r is normalized vector for linear transformation of X
    
    t <- X%*%r # X is tranformed without change the unit
    
    p <- t(X)%*%X%*%r/c(t(r)%*%t(X)%*%X%*%r)
    
    v <- p
    if (a>1) v <- v - V%*%(t(V)%*%p) #make v orthognal to all previous loadings
    v <- v/sqrt(c(crossprod(v))) #get normalized
    
    S <- S - v%*%(t(v)%*%S)
    
    T[,a] <- t
    R[,a] <- r
    V[,a] <- v
    P[,a] <- p
  }
  
  A <- solve(t(T)%*%T)%*%(t(T)%*%y)
  
  B <- R%*%solve(t(T)%*%T)%*%(t(T)%*%y)
  
  list(T=T, R=R, P=P, A=A, B=B)
}
