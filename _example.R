
rm(list = ls()) #clear out all the existing objects in the R environment

setwd('C:/RDU/Research/adj PLS/R codes/to post/')

library(pls)
source('_simu func.R')
source('_adjPLS1 appr.R')
source('_oth funcs.R')

n.obs <- 500
n.avl <- 100
avl <- 1:100

n.sim <- 5

n.bt <- 200
n.null <- 200

p <- 5 #dimension of x
k <- 2 #number of components

B.PLS1 <- NULL
B.adjPLS1 <- NULL

p.plsr <- matrix(NA, ncol=n.sim, nrow=p)
p.adjPLS1 <- matrix(NA, ncol=n.sim, nrow=p)

for (i in 1:n.sim)
{
  #i <- 1
  print(i)
  
  # q=1, p=5, k=2; avl=1:100 
  dat <- mkdata(n=n.obs,
                t.mu=c(1,1), t.sigma=c(2,2),
                P=matrix(c(5,2,0,0,0,0,0,0,1,1), ncol=2, byrow=F),
                g.sigma=rep(1,5),
                f.sigma=1,
                tA=matrix(c(1,3), ncol=1),
                seed=i)
  
  ## the par of the dataset that has complete observations: both covariates and outcome
  X <- dat$X[avl,]
  y <- dat$y[avl,]
  
  X.all <- dat$X #X.all has observations without outcome values
  
  ## model fitting and test statistic calculation
  compt.adjPLS1 <- stat.adjPLS1.wiBT(p=p, k=k, n.bt=n.bt, y=y, X.all=X.all, avl=avl, seed=i)
  
  ## estimation
  if (is.null(B.PLS1)) B.PLS1 <- PLS1(X=X, y=y, k=k)$B else 
    B.PLS1 <- cbind(B.PLS1, PLS1(X=X, y=y, k=k)$B) #estimate from a general PLSR
  
  if (is.null(B.adjPLS1)) B.adjPLS1 <- compt.adjPLS1$B else
    B.adjPLS1 <- cbind(B.adjPLS1, compt.adjPLS1$B) #estimate from adj PLSR
  
  ## observed statistics
  stat.adjPLS1 <- compt.adjPLS1$stat.adjPLS1 #test statistic from adj PLSR
  
  ## bootstrapping generated statistics
  stat.adjPLS1.bt <- matrix(NA, ncol=n.null, nrow=p)
  for (j in 1:n.bt)
  {
    # j <- 1
    print(c(i,j))
    
    seed.ij <- as.numeric(paste0(c(i,j), collapse = ''))
    
    # adjPLS1
    set.seed(seed.ij)
    avl.bt <- sample(1:nrow(X.all), length(y), replace=F)
    
    stat.adjPLS1.bt[,j] <- 
      stat.adjPLS1.wiBT(p=p, k=k, n.bt=n.bt, y=y, X.all=X.all, avl=avl.bt, seed=seed.ij)$stat.adjPLS1
  }
  
  ##### R pkg plsr func for testing with a general PLSR
  fit.plsr <- plsr(y ~ X, method='simpls', ncomp=k, validation = "LOO", jackknife = TRUE)
  tst.plsr <- jack.test(fit.plsr, ncomp=k)

  for (m in 1:p)
  {
    p.plsr[m,i] <- tst.plsr$pvalues[m]
    p.adjPLS1[m,i] <- pval.2side(stat.adjPLS1[m], stat.adjPLS1.bt[m,])
  }
}

