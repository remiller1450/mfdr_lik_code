library(ncvreg)
library(Matrix)

fir.linear <- function(fit, X, true.sig){
  p <- dim(fit$beta)[1]-1
  n.false <- p - 4          #### For this sim set there are 4 true variables
  S <- predict(fit, type="nvars")
  tau <- sqrt(fit$loss/(fit$n - S + 1))
  l <- fit$lam*fit$alpha
  EF <- pmin(2*n.false*pnorm(-l*sqrt(fit$n)/tau),S)
  FIR <- EF/S
  FIR[S==0] <- 0
  df <- cbind(EF=EF, S=S, FIR=FIR)
  structure(df, class="fir")
}
### Function used in calculating set S
nonzero <- function(x){
  return(sum(x != 0))
}


#findlamb <- function(tfirmat, fir){
#  clamb  <- min(which(tfirmat > fir))
#  return(clamb)
#}
findlamb <- function(tfirmat, fir){
  clamb  <- max(which(tfirmat < fir))
  return(clamb)
}


genData <- function(n, J, K=1, beta, family=c("gaussian","binomial","survival"), J0=ceiling(J/2), K0=K, SNR=1, sig = c("homogeneous","heterogeneous"), sig.g = c("homogeneous","heterogeneous"), rho = 0, rho.g = rho, corr=c("exchangeable", "autoregressive")) {
  family <- match.arg(family)
  sig <- match.arg(sig)
  sig.g <- match.arg(sig.g)
  corr <- match.arg(corr)
  
  ## Gen X, S
  if (corr=="exchangeable") {
    X <- genX(n=n, J=J, K=K, rho=rho, rho.g=rho.g)
  } else {
    RHO <- matrix(rho^(0:(J-1)), J, J, byrow=TRUE)
    S <- bandSparse(J, k=0:(J-1), diagonals=RHO, symmetric=TRUE)
    R <- chol(S)
    X <- as.matrix(matrix(rnorm(n*J), n, J) %*% R)
  }
  
  j <- rep(1:J,rep(K,J))
  
  ## Gen beta
  if (missing(beta) || length(beta)==1) {
    k <- rep(1:K,J)
    b <- (j <= J0) * (k <= K0)
    s <- c(1,-1)[1+j%%2] * c(1,-1)[1+k%%2]
    if (missing(beta)) {
      S <- matrix(rho, nrow=J*K, ncol=J*K)
      for (i in 1:J) S[(i-1)*K+1:K,(i-1)*K+1:K] <- rho.g
      diag(S) <- rep(1,J*K)
      if (sig=="heterogeneous") b <- b*j
      if (sig.g=="heterogeneous") b <- b*k
      b <- b*s
      beta <- b*sqrt(SNR)/sqrt(crossprod(b,S)%*%b)
    } else beta <- b*s*beta
  }
  
  ## Gen y
  y <- genY(X%*%beta, family=family, sigma=1)
  return(list(X=X,y=y,beta=beta,family=family,group=j))
}

## rho  : correlation across all explanatory variables
## rho.g: correlation within group (must be at least rho)
genX <- function(n, J, K=1, rho=0, rho.g=rho, corr=corr) {
  a <- sqrt(rho/(1-rho.g))
  b <- sqrt((rho.g-rho)/(1-rho.g))
  Z <- rnorm(n)
  ZZ <- t(matrix(rep(rnorm(n*J), rep(K,n*J)), ncol=n))
  ZZZ <- matrix(rnorm(n*J*K),nrow=n)
  return(matrix(as.numeric(a*Z + b*ZZ + ZZZ),nrow=n)/sqrt(1+a^2+b^2))
}

genY <- function(eta,family=c("gaussian","binomial","survival"),sigma=1,lam.0=0.001) {
  family=match.arg(family)
  n <- length(eta)
  if (family=="gaussian") y <- rnorm(n,mean=eta,sd=sigma)
  else if (family=="binomial")
  {
    pi. <- exp(eta)/(1+exp(eta))
    pi.[eta > log(.9999/.0001)] <- 1
    pi.[eta < log(.0001/.9999)] <- 0
    y <- rbinom(n,1,pi.)
  } else if (family == "survival")
  {
    haz <- lam.0*exp(eta)
    y <- rexp(n, haz)
  }
  return(y)
}



#### Parameters for sim
t <- 4
p <- 40
nn <- c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
true.sig <- 1
nreps <- 1000


#### Outer loop over N
result.lasso <- result.lasso2 <- result <- result2 <-  NULL
for (ii in 1:length(nn)){
  n <- nn[ii]
  b <- 10/sqrt(n)
  beta <- c(rep(c(b,-b),t/2), rep(0,p-t))
  
  ### Generate initial lambda seq
  lams <- NULL
  for (qq in 1:5){
    X <- matrix(rnorm(n*p),nrow = n, ncol = p)
    y <- X %*% beta + rnorm(n, sd = true.sig)
    fit.temp <- ncvreg(X,y, family = "gaussian", penalty = "lasso", nlambda = 3)
    lams <- c(lams, fit.temp$lambda)
  }
  lmax <- max(lams)
  lseq <- seq(1.1*lmax, .01*lmax, length.out = 300)
  
  
  ### Inner loop within a given N
  res.estimate.lasso <- res.estimate.MCP <- res.true.lasso <- res.true.MCP <-  res.estimate.MCP2 <- res.estimate.lasso2 <-NULL
  res.S.MCP <- res.S.lasso <- NULL
  
  for (i in 1:nreps){
    ### Generate Correlated data
    X.A <- genData(n, t, beta = c(rep(b, t/2), rep(-b, t/2)), J0=4, rho=0, family = "gaussian")
    X.C <- genData(n, p-t, rho=0.8, beta=0, corr="auto")
    X <- cbind(X.A$X, X.C$X)
    y <- X.A$y
    
    
    
    ### fit lasso/MCP models and find estimated FIR using both function tailored to this data and ncvreg function
    fit.lasso <- ncvreg(X,y, family = "gaussian", penalty = "lasso", lambda = lseq)
    fir.lasso <- fir.linear(fit.lasso, X, true.sig = true.sig)  
    fir.lasso2 <- mfdr(fit.lasso)
    
    fit.MCP <- ncvreg(X,y, family = "gaussian", penalty = "MCP", lambda = lseq)
    fir.MCP <- fir.linear(fit.MCP, X, true.sig = true.sig)
    fir.MCP2 <- mfdr(fit.MCP)
    
    ### find true false selection rates
    fi.lasso <- apply(fit.lasso$beta[(t+2):(p+1),], 2, nonzero)  #### first beta is intercept, next t are true vars, so t+2 is first false
    
    fi.MCP <- apply(fit.MCP$beta[(t+2):(p+1),], 2, nonzero)
    
    
    
    ### Store results for current iteration
    res.estimate.lasso <- cbind(res.estimate.lasso, fir.lasso[,1])
    res.estimate.MCP <- cbind(res.estimate.MCP, fir.MCP[,1])
    
    res.estimate.lasso2 <- cbind(res.estimate.lasso2, fir.lasso2[,1])
    res.estimate.MCP2 <- cbind(res.estimate.MCP2, fir.MCP2[,1])
    
    res.true.lasso <- cbind(res.true.lasso, fi.lasso)
    res.true.MCP <- cbind(res.true.MCP, fi.MCP)
    
    res.S.MCP <- cbind(res.S.MCP, fir.MCP[,2])
    res.S.lasso <- cbind(res.S.lasso, fir.lasso[,2])
  }
  
  ### sum across iterations 
  mres.est.lasso <- apply(res.estimate.lasso, 1, sum, na.rm = TRUE)
  mres.est.MCP <- apply(res.estimate.MCP, 1, sum, na.rm = TRUE)
  mres.est.lasso2 <- apply(res.estimate.lasso2,1, sum, na.rm = TRUE)
  mres.est.MCP2 <- apply(res.estimate.MCP2, 1, sum, na.rm = TRUE)
  mres.true.lasso <- apply(res.true.lasso, 1, sum, na.rm = TRUE)
  mres.true.MCP <- apply(res.true.MCP, 1, sum, na.rm = TRUE)
  mres.S.lasso <- apply(res.S.lasso, 1, sum, na.rm = TRUE)
  mres.S.MCP <- apply(res.S.MCP, 1, sum, na.rm = TRUE)
  
  ### Calculate mfdr as sum/sum at each lambda
  mfdr.est <- mres.est.MCP/mres.S.MCP
  mfdr.est2 <- mres.est.MCP2/mres.S.MCP
  mfdr.est.lasso <- mres.est.lasso/mres.S.lasso
  mfdr.est.lasso2 <- mres.est.lasso2/mres.S.lasso
  mfdr.true <- mres.true.MCP/mres.S.MCP
  mfdr.true.lasso <- mres.true.lasso/mres.S.lasso
  
  ### Store results after averaging
  result <- c(result, mfdr.est[findlamb(mfdr.true, fir = .1)])
  result.lasso <- c(result.lasso,   mfdr.est.lasso[findlamb(mfdr.true.lasso, fir = .1)])
  result2 <- c(result2, mfdr.true[findlamb(mfdr.est, fir = .1)])
  result.lasso2 <- c(result.lasso2,   mfdr.true.lasso[findlamb(mfdr.est.lasso, fir = .1)])
  
  save(result, result.lasso, n, file = "sim1_linear_cor.RData")
}



