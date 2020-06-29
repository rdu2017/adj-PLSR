

#### adjPLS1: function for model fitting
####   k: the number of the components; 
####   avl: the indicies of X.all that have the assoicated outcome y

adjPLS1 <- function(X.all, y, k, avl=c(1:100))  
{
  X <- X.all[avl,]
  
  R <- matrix(NA, nrow = ncol(X), ncol = k)
  V <- matrix(0, nrow = ncol(X), ncol = k)
  P <- matrix(NA, nrow = ncol(X), ncol = k)
  T.all <- matrix(NA, nrow = nrow(X.all), ncol = k)
  
  S <- cov(X,y)
  
  for (a in 1:k)
  {
    r <- S
    r <- r/sqrt(c(crossprod(r))) # r is normalized vector for linear transformation of X
    
    t <- X.all%*%r # X is tranformed without change the unit
    
    p <- t(X.all)%*%X.all%*%r/c(t(r)%*%t(X.all)%*%X.all%*%r) # here use X.all for estimation of x-loadings
    
    v <- p
    if (a>1) v <- v - V%*%(t(V)%*%p) #make v orthognal to all previous loadings
    v <- v/sqrt(c(crossprod(v))) #get normalized
    
    S <- S - v%*%(t(v)%*%S)
    
    T.all[,a] <- t
    R[,a] <- r
    V[,a] <- v
    P[,a] <- p
  }
  
  T <- T.all[avl,]
  A <- solve(t(T)%*%T)%*%(t(T)%*%y)
  
  B <- R%*%solve(t(T)%*%T)%*%(t(T)%*%y)
  
  list(T=T, R=R, P=P, A=A, B=B, T.all=T.all)
}


#### stat.adjPLS1.wiBT: function for calculating the testing statistics
####   p: the number of covariates
####   n.bt: the number of Bootstrapping proceures for the variance estimate

stat.adjPLS1.wiBT <- function(p, k, n.bt, y, X.all, avl, seed=NULL)
{
  if (!is.null(seed)) set.seed(seed)
  
  out.adjPLS1 <- adjPLS1(X.all, y, k, avl=avl)
  
  B.adjPLS1.bt <- matrix(NA, ncol=n.bt, nrow=p)
  for (j in 1:n.bt)
  {
    avl.bt <- sample(avl, length(avl), replace=T) #resampling the part that has y
    oth.bt <- sample((1:nrow(X.all))[-avl], length((1:nrow(X.all))[-avl]), replace=T) #resampling the part that does not have y
    X.all.bt <- X.all[c(avl.bt, oth.bt),] #combine the two parts from having y to not having y
    
    y.bt <- y[match(avl.bt, avl)] #y.bt is ordered according to the resamplied indecies
    
    out.adjPLS1.bt <- adjPLS1(X.all.bt, y.bt, k, avl=1:length(avl.bt)) #
    B.adjPLS1.bt[,j] <- out.adjPLS1.bt$B
  }
  
  stat.adjPLS1 <- c(out.adjPLS1$B)/sqrt(diag(cov(t(B.adjPLS1.bt))))
  
  list(B=out.adjPLS1$B, stat.adjPLS1=stat.adjPLS1)
}

