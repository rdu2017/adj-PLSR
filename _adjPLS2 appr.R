

#### adjPLS2: function for model fitting
####   k: the number of the components; 
####   avl: the indicies of X.all that have the assoicated outcome y

adjPLS2.c <- function(X.all, Y, k, avl=c(1:100)) #avl: available part, the row indices of X.all 
{
  X <- X.all[avl,]
  
  R <- matrix(NA, nrow = ncol(X), ncol = k)
  V <- matrix(0, nrow = ncol(X), ncol = k)
  P <- matrix(NA, nrow = ncol(X), ncol = k)
  T.all <- matrix(NA, nrow = nrow(X.all), ncol = k)
  
  Q <- matrix(NA, nrow=ncol(Y), ncol=k) # Y linear combination;
  U <- matrix(NA, nrow=nrow(Y), ncol=k) # Y scores
  
  S <- cov(X,Y)
  
  for (a in 1:k)
  {
    q <- eigen(crossprod(S), symmetric = TRUE)$vectors[,1] #sysmmetric being assumed or not
    
    r <- S %*% q  #t(q)%*%t(S)%*%r
    r <- r/sqrt(c(crossprod(r))) # r is normalized for linear transformation of X
    
    t <- X.all%*%r # X is tranformed without change the unit
    
    p <- t(X.all)%*%X.all%*%r/c(t(r)%*%t(X.all)%*%X.all%*%r) # here use larger dataset for estimation of x-loadings
    
    u <- Y%*%q #by using the eigen func, the q has already been normalized
    
    v <- p
    if (a>1) v <- v - V%*%(t(V)%*%p) #make v orthognal to all previous loadings
    v <- v/sqrt(c(crossprod(v))) #get normalized
    
    S <- S - v%*%(t(v)%*%S)
    
    T.all[,a] <- t
    R[,a] <- r
    V[,a] <- v
    P[,a] <- p
    Q[,a] <- q
    U[,a] <- u
  }
  
  T <- T.all[avl,]
  A <- solve(t(T)%*%T)%*%(t(T)%*%Y)
  
  B <- R%*%solve(t(T)%*%T)%*%(t(T)%*%Y)
  
  L <- R%*%solve(t(T)%*%T)%*%t(T)
  
  list(T=T, R=R, P=P, Q=Q, U=U, A=A, B=B, L=L, T.all=T.all)
}


#### stat.adjPLS2.wiBT: function for calculating the testing statistics
####   p: the number of covariates
####   n.bt: the number of Bootstrapping proceures for the variance estimate

stat.adjPLS2.wiBT <- function(p, k, n.bt, Y, X.all, avl, seed=NULL)
{
  if (!is.null(seed)) set.seed(seed)
  
  out.adjPLS2 <- adjPLS2.c(X.all, Y, k, avl=avl)
  
  B.adjPLS2.bt <- list()
  for (j in 1:n.bt)
  {
    avl.bt <- sample(avl, length(avl), replace=T) #resampling the part that has Y
    oth.bt <- sample((1:nrow(X.all))[-avl], length((1:nrow(X.all))[-avl]), replace=T) #resampling the part that does not have y
    X.all.bt <- X.all[c(avl.bt, oth.bt),] #combine the two parts from having y to not having y
    
    Y.bt <- Y[match(avl.bt, avl),] #y.bt is ordered according to the resampled indicies
    
    out.adjPLS2.bt <- adjPLS2.c(X.all.bt, Y.bt, k, avl=1:length(avl.bt)) #
    B.adjPLS2.bt[[length(B.adjPLS2.bt)+1]] <- out.adjPLS2.bt$B
  }
  
  B.var <- matrix(NA, nrow=p, ncol=ncol(Y))
  for (i in 1:p)
    for (j in 1:ncol(Y))
    {
      B.var[i,j] <- var(sapply(B.adjPLS2.bt, function(x) x[i,j]))
    }
  
  stat.adjPLS2 <- out.adjPLS2$B/sqrt(B.var)
  
  stat.adjPLS2
}

