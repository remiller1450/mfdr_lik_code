library(ncvreg)
library(survival)
library(Matrix)

fir.logistic <- function(fit, X){
  X <- X
  PW <-  predict(fit, X, type = 'response')
  S <- predict(fit, type="nvars")
  l <- fit$lambda*fit$alpha
  pz <- which(fit$penalty.factor != 0)
  p <- 40
  t <- 4
  EF <- WW <- NULL
  for (i in 1:length(l)){
    ### Calculate W matrix
    W <- diag(as.vector(PW[,i]*(1 - PW[,i])))
    
    ### Calculate x'Wx
    for (j in 1:ncol(X)){
      WW[j] <- t(fit$X[,j]) %*% W %*% fit$X[,j]
    }
    
    ### Calculate expected false inclusions
    EF[i] <- 2*sum(pnorm(-(nrow(X)*l[i])/sqrt(WW))[(t+1):p])
  }
  EF <- pmin(EF, S)
  FIR <- EF/S
  FIR[S==0] <- 0
  df <- data.frame(EF=EF, S=S, FIR=FIR)
  structure(df, class=c("fir", "data.frame"))
}



### Function used in calculating set S
nonzero <- function(x){
  return(sum(x != 0))
}

findlamb <- function(tfirmat, fir){
  clamb  <- min(which(tfirmat > fir))
  return(clamb)
}



genData <- function(n, J, K=1, beta, family=c("gaussian","binomial"), J0=ceiling(J/2), K0=K, SNR=1, sig = c("homogeneous","heterogeneous"), sig.g = c("homogeneous","heterogeneous"), rho = 0, rho.g = rho, corr=c("exchangeable", "autoregressive")) {
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

genY <- function(eta,family=c("gaussian","binomial"),sigma=1) {
  family=match.arg(family)
  n <- length(eta)
  if (family=="gaussian") y <- rnorm(n,mean=eta,sd=sigma)
  else if (family=="binomial")
  {
    pi. <- exp(eta)/(1+exp(eta))
    pi.[eta > log(.9999/.0001)] <- 1
    pi.[eta < log(.0001/.9999)] <- 0
    y <- rbinom(n,1,pi.)
  }
  return(y)
}



nn <- c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)

result <- result2 <- f.est <- NULL
result.lasso <- result2.lasso <- f.est.lasso <- NULL
for (i in 1:length(nn)){
  t <- 4
  p <- 40
  n <- nn[i]
  nsim <- 1000
  res.est <- matrix(NA, nrow = 200, ncol = nsim)
  res.true <- matrix(NA, nrow = 200, ncol = nsim)
  res.est.lasso <- matrix(NA, nrow = 200, ncol = nsim)
  res.true.lasso <- matrix(NA, nrow = 200, ncol = nsim)
  res.S <- matrix(NA, nrow = 200, ncol = nsim)
  res.S.lasso <- matrix(NA, nrow = 200, ncol = nsim)
  
  ### Gen lambda seq
  lmax <- lmax.lasso <- NULL
  for(q in 1:5){
    ### Generate Autocorrelated data
    X.A <- genData(n, t, beta = c(rep(14/sqrt(n), t/2), rep(-14/sqrt(n), t/2)), J0=4, rho=0, family = "binomial")
    X.C <- genData(n, p-t, rho=0.8, beta=0, corr="auto")
    yy <- X.A$y
    X <- cbind(X.A$X, X.C$X)
    ### Fit Model
    fit.lam <- ncvreg(X, yy, family = "binomial", returnX = TRUE, penalty = 'lasso', nlambda = 3)
    lmax <- c(lmax, max(fit.lam$lambda))
  }
  lseq <- seq(1.2*max(lmax),0.1*max(lmax), length.out = 200)
  lseq.lasso <- seq(1.3*max(lmax), 0.1*max(lmax), length.out = 200)
  
  for (k in 1:nsim){
    ### Generate Data
    X.A <- genData(n, t, beta = c(rep(14/sqrt(n), t/2), rep(-14/sqrt(n), t/2)), J0=4, rho=0, family = "binomial")
    X.C <- genData(n, p-t, rho=0.8, beta=0, corr="auto")
    yy <- X.A$y
    X <- cbind(X.A$X, X.C$X)
    
    ### Fit Model
    fit <- ncvreg(X, yy, returnX = TRUE, family = "binomial", penalty = 'MCP', lambda = lseq)
    fit.lasso <-  ncvreg(X, yy, returnX = TRUE, family = "binomial", penalty = 'lasso', lambda = lseq)
    
    
    
    ### Number of selections
    S <- apply(fit$beta, 2, nonzero) - 1     ### -1 for intercept 
    S.lasso <- apply(fit.lasso$beta, 2, nonzero) - 1  ### -1 for intercept 
    
    ### Calcutate true FIR
    tfd1 = apply(fit$beta[(t+2):(p+1),], 2, nonzero)  ### True number of false inclusions
    ti1 = apply(fit$beta[2:(t+1),], 2, nonzero)       ### Number of correct inclusions 
    
    
    tfd1.lasso = apply(fit.lasso$beta[(t+2):(p+1),], 2, nonzero)  ### True number of false inclusions
    ti1.lasso = apply(fit.lasso$beta[2:(t+1),], 2, nonzero)       ### Number of correct inclusions 
    
    
    
    ### Estimate FIR
    FI <- fir.logistic(fit, X)
    FI.lasso <- fir.logistic(fit.lasso, X)
    
    
    ### Make everything length 200 (model saturation can lead to it being shorter), everything after saturation is NA and excluded from future calculations
    if(length(FI$EF) < 200){
      FI.MCP <- c(FI$EF, rep(NA,(200-length(FI$EF))))
      S.MCP <- c(FI$S, rep(NA,(200-length(FI$S))))
    } else {
      FI.MCP <- FI$EF
      S.MCP <- FI$S
    }
    if(length(FI.lasso$EF) < 200){
      FI.LASSO <- c(FI.lasso$EF, rep(NA,(200-length(FI.lasso$EF))))
      S.LASSO <- c(FI.lasso$S, rep(NA,(200-length(FI.lasso$S))))
    } else {
      FI.LASSO <- FI.lasso$EF
      S.LASSO <- FI.lasso$S
    }
    if(length(tfd1) < 200){
      FI.MCP.TRUE <- c(tfd1, rep(NA,(200-length(tfd1))))
    } else {
      FI.MCP.TRUE <- tfd1
    }
    
    if(length(tfd1.lasso) < 200){
      FI.LASSO.TRUE <- c(tfd1.lasso, rep(NA,(200-length(tfd1.lasso))))
    } else {
      FI.LASSO.TRUE <- tfd1.lasso
    }    
    
    ### Store current iteration
    res.est[,k] <- FI.MCP
    res.true[,k] <- FI.MCP.TRUE
    res.est.lasso[,k] <- FI.LASSO
    res.true.lasso[,k] <- FI.LASSO.TRUE
    res.S[,k] <- S.MCP
    res.S.lasso[,k] <- S.LASSO
  }
  
  ### Sum across the iterations
  mres.est <- apply(res.est, 1, sum, na.rm = TRUE)
  mres.true <- apply(res.true, 1, sum, na.rm = TRUE)
  mres.est.lasso <- apply(res.est.lasso, 1, sum, na.rm = TRUE)
  mres.true.lasso <- apply(res.true.lasso, 1, sum, na.rm = TRUE)
  mres.S <- apply(res.S, 1, sum, na.rm = TRUE)
  mres.S.lasso <- apply(res.S.lasso, 1, sum, na.rm = TRUE)
  
  ### Calculate mfdr as sum/sum at each lambda
  mfdr.est <- mres.est/mres.S
  mfdr.est.lasso <- mres.est.lasso/mres.S.lasso
  mfdr.true <- mres.true/mres.S
  mfdr.true.lasso <- mres.true.lasso/mres.S.lasso
  
  result <- c(result, mfdr.est[findlamb(mfdr.true, fir = .1)])
  result.lasso <- c(result.lasso,   mfdr.est.lasso[findlamb(mfdr.true.lasso, fir = .1)])
  result2 <- c(result2, mfdr.true[findlamb(mfdr.est, fir = .1)])
  result2.lasso <- c(result2.lasso,   mfdr.true.lasso[findlamb(mfdr.est.lasso, fir = .1)])
  
  save(result, result.lasso, nn, file = "sim1_logistic_ind.RData")
}

