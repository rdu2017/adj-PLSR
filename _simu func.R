
library('MASS')

mkdata <- function(n=500, #number of the total observations to simulate
                   t.mu=c(1,1), t.sigma=c(1, 1), #t: k X 1 vector, T: n X k
                   P=matrix(c(5,2,0,0,0,0,0,0,1,1), ncol=2, byrow=F), #P: p X k, containing vector of p as x-loadings
                   g.mu=c(0,0,0,0,0), g.sigma=c(1, 1, 1, 1,1), #g: p X 1, residual vector
                   tA=matrix(c(1,3), ncol=1), #tA: k X q, coefficient 
                   f.mu=0, f.sigma=1, #f: q X 1
                   seed = NULL)
{
  if (!is.null(seed)) set.seed(seed) #set a ranodm seed for reproducible simulations if needed
  
  P <- apply(P, 2, function(x) x/sqrt(c(crossprod(x)))) #normalize the loading vector
  
  if (length(t.sigma)==1) T <- mvrnorm(n=n, mu=t.mu, Sigma=diag(t.sigma, nrow=1, ncol=1)) else
    T <- mvrnorm(n=n, mu=t.mu, Sigma=diag(t.sigma)) #simulate T
  
  if (length(g.sigma)==1) X <- T%*%t(P) + mvrnorm(n=n, mu=g.mu, Sigma=diag(g.sigma, nrow=1, ncol=1)) else
    X <- T%*%t(P) + mvrnorm(n=n, mu=g.mu, Sigma=diag(g.sigma)) #simulate X
  
  if (length(f.sigma)==1) y <- T%*%tA + mvrnorm(n=n, mu=f.mu, Sigma=diag(f.sigma, nrow=1, ncol=1)) else
    y <- T%*%tA + mvrnorm(n=n, mu=f.mu, Sigma=diag(f.sigma)) #simulate Y
  
  B <- t(t(P)/rowSums(t(P)))%*%(tA) #the true values in B
  
  list(T=T, X=X, y=y, B=B) #output a list 
}


