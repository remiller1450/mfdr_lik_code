#######################################################################################
#
#                                                                                     
#   Filename    :	mfdr_lik_reproduce.R										  
#   Input data files  :    ---                                                        
#   Output data files :    fig1.pdf, fig2.pdf, fig3.pdf, fig4.pdf, fig5.pdf,
#			   table1.pdf, table2.pdf 
#
#   Required R packages :  Rccp, Matrix, ncvreg, survival, ggplot2, covTest, 
#			   SelectiveInference, gridExtra, GEOQuery
#
#
########################################################################################


#install.packages("ncvreg")
#install.packages('covTest')
#install.packages('selectiveInference')
#install.packages('gridExtra')
#install.packages('ggplot2')
#install.packages('Rccp')
library(ggplot2)
library(ncvreg)
library(Matrix)
library(survival)
library(covTest)
library(selectiveInference)
library(Rccp)
library(gridExtra)


## Set seed
set.seed(123456)

########################################################################################
#
#    Setup custom functions to be used throughout simulations
#
########################################################################################

## C++ for calculating mfdr for cox models
Rcpp::cppFunction(
  'NumericVector mfdr_cox(int n,
  int L,
  int p,
  int t,
  double alpha,
  NumericVector lambda,
  NumericVector Eta,
  NumericVector d,
  NumericMatrix X,
  NumericVector S
) {
  NumericVector w(n);
  NumericVector rsk(n);
  NumericVector haz(n);
  NumericVector EF(L);
  double tau = 0;
  int nn = 0;
  NumericVector stt(1);
  NumericVector tmp(1);
  
  
  for (int l=0; l<L; l++) {
  for (int i=0; i<n; i++) haz[i] = exp(Eta[n*l+i]);
  rsk[n-1] = haz[n-1];
  for (int i=n-2; i>=0; i--) rsk[i] = rsk[i+1] + haz[i];
  for (int j=0; j<n; j++) {
  w[j] = 0;
  for (int i=0; i <= j; i++) {
  w[j] += d[i]*haz[j]/rsk[i]*(1-haz[j]/rsk[i]);
  }
  }
  tmp[0] = 0;
  for (int j=t; j<p; j++) {
  tau = 0;
  nn = n*j;
  for (int m=0;m<n;m++) tau += w[m] * pow(X[nn+m], 2);
  stt[0] = -sqrt(n)*lambda[l]*alpha/sqrt(tau/n);
  tmp += 2*pnorm(stt, 0.0, 1.0, 1, 0);
  EF[l] = tmp[0];
  }
  }
  
  return(EF/S); }'
)

## C++ function for calculating mfdr for logistic regression models
Rcpp::cppFunction(
  'NumericVector mfdr_log(int n,
  int L,
  int p,
  int t,
  double alpha,
  NumericVector lambda,
  NumericVector Eta,
  NumericMatrix X,
  NumericVector S
) {
  NumericVector w(n);
  NumericVector EF(L);
  double tau = 0;
  double pi = 0;
  int nn = 0;
  NumericVector stt(1);
  NumericVector tmp(1);
  
  
 for (int l=0; l<L; l++) {
  for (int i=0; i<n; i++) {
    pi = exp(Eta[n*l+i])/(1 + exp(Eta[n*l+i]));
    w[i] = pi*(1-pi);
  }
  tmp[0] = 0;
  for (int j=t; j<p; j++) {
    tau = 0;
    nn = n*j;
      for (int m=0;m<n;m++) tau += w[m] * pow(X[nn+m], 2);
        stt[0] = -sqrt(n)*lambda[l]*alpha/sqrt(tau/n);
        tmp += 2*pnorm(stt, 0.0, 1.0, 1, 0);
        EF[l] = tmp[0];
    }
  }
  
  return(EF/S); }'
)

### Short function to find lambda of interest (used later)
findlam <- function(x, thresh = .1){
  max(which(x <= thresh))
}

## Functions to generate data with correlated predictors
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

genY <- function(eta,family=c("gaussian","binomial","survival"),sigma=1,lam.0=1) {
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


########################################################################################
#
#    Simulation #1 (Produces Figure 1)
#
########################################################################################

### Initial setup
nn <- seq(200, 1600, by = 200)
p <- 100
t <- 10
nreps <- 1000

icox.res.est <- ilog.res.est <- ilin.res.est <- icox.res.true <- ilog.res.true <- ilin.res.true<- matrix(NA,100,length(nn))
ccox.res.est <- clog.res.est <- clin.res.est <- ccox.res.true <- clog.res.true <- clin.res.true<- matrix(NA,100,length(nn))
icox.res.est.lasso <- ilog.res.est.lasso <- ilin.res.est.lasso <- icox.res.true.lasso <- ilog.res.true.lasso <- ilin.res.true.lasso <- matrix(NA,100,length(nn))
ccox.res.est.lasso <- clog.res.est.lasso <- clin.res.est.lasso <- ccox.res.true.lasso <- clog.res.true.lasso <- clin.res.true.lasso <- matrix(NA,100,length(nn))

### Loop over sample sizes
for(a in 1:length(nn)){
  n <- nn[a]
  beta.cox <- c(rep(4/sqrt(n),t), rep(0, p-t))
  beta.log <- c(rep(10/sqrt(n),t), rep(0, p-t))
  beta.lin <- c(rep(6/sqrt(n),t), rep(0, p-t))
  
  ## Generate lambda sequence to use for this n
  ilog.max.lam <- icox.max.lam <- ilin.max.lam <- clog.max.lam <- ccox.max.lam <- clin.max.lam <- 0
  for(ll in 1:3){
    X <- matrix(rnorm(n*p), n, p)
    y.cox <- survival::Surv(rexp(n, rate = exp(X %*% beta.cox)), rep(1, n))
    y.log <- rbinom(n, size = 1, prob = exp(X %*% beta.log)/(1 + exp(X %*% beta.log)))
    y.lin <- rnorm(n, X%*%beta.lin)
    
    ## Correlated data
    X.cox <- genData(n, t, beta = beta.cox[1:t], J0=t, rho=0, family = "survival")
    X.log <- genData(n, t, beta = beta.log[1:t], J0=t, rho=0, family = "binomial")
    X.lin <- genData(n, t, beta = beta.lin[1:t], J0=t, rho=0, family = "gaussian")
    X.C <- genData(n, p-t, rho=0.8, beta=0, corr="auto")
    yc.cox <- Surv(X.cox$y, rep(1, n))
    yc.log <- X.log$y
    yc.lin <- X.lin$y
    XC.cox <- as.matrix(cbind(X.cox$X, X.C$X))
    XC.log <- as.matrix(cbind(X.log$X, X.C$X))
    XC.lin <- as.matrix(cbind(X.lin$X, X.C$X))
    
    ## Fit models
    icox.fit.temp <- ncvsurv(X, y.cox)
    ccox.fit.temp <- ncvsurv(XC.cox, yc.cox)
    ilog.fit.temp <- ncvreg(X, y.log, family = "binomial")
    clog.fit.temp <- ncvreg(XC.log, y.log, family = "binomial")
    ilin.fit.temp <- ncvreg(X, y.lin)
    clin.fit.temp <- ncvreg(XC.lin, y.lin)
    
    ## Extract max of each default lambda seq
    icox.max.lam <- max(icox.max.lam, max(icox.fit.temp$lambda))
    ilog.max.lam <- max(ilog.max.lam, max(ilog.fit.temp$lambda))
    ilin.max.lam <- max(ilin.max.lam, max(ilin.fit.temp$lambda))
    
    ccox.max.lam <- max(ccox.max.lam, max(ccox.fit.temp$lambda))
    clog.max.lam <- max(clog.max.lam, max(clog.fit.temp$lambda))
    clin.max.lam <- max(clin.max.lam, max(clin.fit.temp$lambda))
  }
  
  ### Create safe lambda sequences
  icox.lseq <- icox.lseq.lasso <- seq(1.1*icox.max.lam, .05*icox.max.lam, length.out = 100)
  ilog.lseq <- ilog.lseq.lasso <- seq(1.1*ilog.max.lam, .05*ilog.max.lam, length.out = 100)
  ilin.lseq <- ilin.lseq.lasso <- seq(1.1*ilin.max.lam, .05*ilin.max.lam, length.out = 100)
  
  ccox.lseq <- ccox.lseq.lasso <- seq(1.1*ccox.max.lam, .05*ccox.max.lam, length.out = 100)
  clog.lseq <- clog.lseq.lasso <- seq(1.1*clog.max.lam, .05*clog.max.lam, length.out = 100)
  clin.lseq <- clin.lseq.lasso <- seq(1.1*clin.max.lam, .05*clin.max.lam, length.out = 100)
  
  
  ### Simulate within each sample size
  icox.true <- ilog.true <- ilin.true <- icox.est <- ilog.est <- ilin.est <- matrix(NA, 100, nreps)
  ccox.true <- clog.true <- clin.true <- ccox.est <- clog.est <- clin.est <- matrix(NA, 100, nreps)
  icox.true.lasso <- ilog.true.lasso <- ilin.true.lasso <- icox.est.lasso <- ilog.est.lasso <- ilin.est.lasso <- matrix(NA, 100, nreps)
  ccox.true.lasso <- clog.true.lasso <- clin.true.lasso <- ccox.est.lasso <- clog.est.lasso <- clin.est.lasso <- matrix(NA, 100, nreps)
  for(i in 1:nreps){
    
    ## Uncorrelated data
    X <- matrix(rnorm(n*p), n, p)
    y.cox <- survival::Surv(rexp(n, rate = exp(X %*% beta.cox)), rep(1, n))
    y.log <- rbinom(n, size = 1, prob = exp(X %*% beta.log)/(1 + exp(X %*% beta.log)))
    y.lin <- rnorm(n, X%*%beta.lin)
    
    ## Correlated data
    X.cox <- genData(n, t, beta = beta.cox[1:t], J0=t, rho=0, family = "survival")
    X.log <- genData(n, t, beta = beta.log[1:t], J0=t, rho=0, family = "binomial")
    X.lin <- genData(n, t, beta = beta.lin[1:t], J0=t, rho=0, family = "gaussian")
    X.C <- genData(n, p-t, rho=0.8, beta=0, corr="auto")
    yc.cox <- Surv(X.cox$y, rep(1, n))
    yc.log <- X.log$y
    yc.lin <- X.lin$y
    XC.cox <- as.matrix(cbind(X.cox$X, X.C$X))
    XC.log <- as.matrix(cbind(X.log$X, X.C$X))
    XC.lin <- as.matrix(cbind(X.lin$X, X.C$X))
    
    ### MCP MODELS ###
    
    ## Fit cox model on ind data
    icox.fit <- ncvsurv(X, y.cox, returnX = TRUE, lambda = icox.lseq)
    mfcox.temp <- mfdr_cox(n = n, L = length(icox.fit$lambda), p = p, t = t, alpha = 1, d = icox.fit$fail,
                           lambda = icox.fit$lambda, X = icox.fit$X, Eta = icox.fit$Eta, S = predict(icox.fit, type = 'nvars'))
    mfcox.temp[is.infinite(mfcox.temp)] = 0
    mfcox.temp[mfcox.temp > 1] = 1
    if(length(mfcox.temp) == 100){
      icox.est[,i] <- mfcox.temp
    } else{icox.est[,i] <- c(mfcox.temp, rep(1, 100-length(mfcox.temp)))} ## If seq stopped early 
    tfcox.temp <- colSums(icox.fit$beta[(t+1):p,] != 0 )/colSums(icox.fit$beta != 0 )
    if(length(tfcox.temp) == 100){
      icox.true[,i] <- tfcox.temp
    } else {icox.true[,i] <- c(tfcox.temp, rep(1, 100-length(tfcox.temp)))}
    #icox.est[,i] <- mfdr(icox.fit)$mFDR
    #icox.true[,i] <- colSums(icox.fit$beta[(t+1):p,] != 0 )/colSums(icox.fit$beta != 0 )
    icox.true[is.na(icox.true)] <- 0
    
    ## Fit cox model on cor data
    ccox.fit <- ncvsurv(XC.cox, yc.cox, returnX = TRUE, lambda = ccox.lseq)
    cmfcox.temp <- mfdr_cox(n = n, L = length(ccox.fit$lambda), p = p, t = t, alpha = 1, d = ccox.fit$fail,
                           lambda = ccox.fit$lambda, X = ccox.fit$X, Eta = ccox.fit$Eta, S = predict(ccox.fit, type = 'nvars'))
    cmfcox.temp[is.infinite(cmfcox.temp)] = 0
    cmfcox.temp[cmfcox.temp > 1] = 1
    if(length(mfcox.temp) == 100){
      ccox.est[,i] <- cmfcox.temp
    } else{ccox.est[,i] <- c(cmfcox.temp, rep(1, 100-length(cmfcox.temp)))} ## If seq stopped early 
    ctfcox.temp <- colSums(ccox.fit$beta[(t+1):p,] != 0 )/colSums(ccox.fit$beta != 0 )
    if(length(ctfcox.temp) == 100){
      ccox.true[,i] <- ctfcox.temp
    } else {ccox.true[,i] <- c(ctfcox.temp, rep(1, 100-length(ctfcox.temp)))}
    ccox.true[is.na(ccox.true)] <- 0
    
    ## Fit logistic model on ind data
    ilog.fit <- ncvreg(X, y.log, family = "binomial", returnX = TRUE, lambda = ilog.lseq)
    mflog.temp  <- mfdr_log(n = n, L = length(ilog.fit$lambda), p = p, t = t, alpha = 1,
                                  lambda = ilog.fit$lambda, X = ilog.fit$X, Eta = ilog.fit$Eta,
                                  S = predict(ilog.fit, type = 'nvars'))
    mflog.temp[is.infinite(mflog.temp)] = 0
    mflog.temp[mflog.temp > 1] = 1
    if(length(mflog.temp) == 100){
      ilog.est[,i] <- mflog.temp
    } else{ilog.est[,i] <- c(mflog.temp, rep(1, 100-length(mflog.temp)))} ## If seq stopped early 
    tflog.temp <- colSums(ilog.fit$beta[(t+2):(p+1),] != 0 )/colSums(ilog.fit$beta[-1,] != 0 )
    if(length(tflog.temp) == 100){
      ilog.true[,i] <- tflog.temp
    } else {ilog.true[,i] <- c(tflog.temp, rep(1, 100-length(tflog.temp)))}
    ilog.true[is.na(ilog.true)] <- 0
    
    ## Fit logistic model on cor data
    clog.fit <- ncvreg(XC.log, yc.log, family = "binomial", returnX = TRUE, lambda = clog.lseq)
    cmflog.temp  <- mfdr_log(n = n, L = length(clog.fit$lambda), p = p, t = t, alpha = 1,
                                   lambda = clog.fit$lambda, X = clog.fit$X, Eta = clog.fit$Eta,
                                   S = predict(clog.fit, type = 'nvars'))
    cmflog.temp[is.infinite(cmflog.temp)] = 0
    cmflog.temp[cmflog.temp > 1] = 1
    if(length(cmflog.temp) == 100){
      clog.est[,i] <- cmflog.temp
    } else{clog.est[,i] <- c(cmflog.temp, rep(1, 100-length(cmflog.temp)))} ## If seq stopped early 
    ctflog.temp <- colSums(clog.fit$beta[(t+2):(p+1),] != 0 )/colSums(clog.fit$beta[-1,] != 0 )
    if(length(ctflog.temp) == 100){
      clog.true[,i] <- ctflog.temp
    } else {clog.true[,i] <- c(ctflog.temp, rep(1, 100-length(ctflog.temp)))}
    clog.true[is.na(clog.true)] <- 0
    
    ## Fit lm on ind data
    ilin.fit <- ncvreg(X, y.lin, returnX = TRUE, lambda = ilin.lseq)
    ilin.est[,i] <- mfdr(ilin.fit)$mFDR
    ilin.true[,i] <- colSums(ilin.fit$beta[(t+2):(p+1),] != 0 )/colSums(ilin.fit$beta[-1,] != 0 )
    ilin.true[is.na(ilin.true)] <- 0
    
    ## Fit lm on cor data
    clin.fit <- ncvreg(XC.lin, yc.lin, returnX = TRUE, lambda = clin.lseq)
    clin.est[,i] <- mfdr(clin.fit)$mFDR
    clin.true[,i] <- colSums(clin.fit$beta[(t+2):(p+1),] != 0 )/colSums(clin.fit$beta[-1,] != 0 )
    clin.true[is.na(clin.true)] <- 0
    
    ### LASSO MODELS ###
    
    ## Fit cox model on ind data
    icox.fit.lasso <- ncvsurv(X, y.cox, returnX = TRUE, lambda = icox.lseq.lasso, penalty = "lasso")
    mfcox.temp.lasso  <- mfdr_cox(n = n, L = length(icox.fit.lasso$lambda), p = p, t = t, alpha = 1, d = icox.fit.lasso$fail,
                           lambda = icox.fit.lasso$lambda, X = icox.fit.lasso$X, Eta = icox.fit.lasso$Eta,
                           S = predict(icox.fit.lasso, type = 'nvars'))
    mfcox.temp.lasso[is.infinite(mfcox.temp.lasso)] = 0
    mfcox.temp.lasso[mfcox.temp.lasso > 1] = 1
    if(length(mfcox.temp.lasso) == 100){
      icox.est.lasso[,i] <- mfcox.temp.lasso
    } else{icox.est.lasso[,i] <- c(mfcox.temp.lasso, rep(1, 100-length(mfcox.temp.lasso)))} ## If seq stopped early 
    tfcox.temp.lasso <- colSums(icox.fit.lasso$beta[(t+1):p,] != 0 )/colSums(icox.fit.lasso$beta != 0 )
    if(length(tfcox.temp.lasso) == 100){
      icox.true.lasso[,i] <- tfcox.temp.lasso
    } else {icox.true.lasso[,i] <- c(tfcox.temp.lasso, rep(1, 100-length(tfcox.temp.lasso)))}
    #icox.est[,i] <- mfdr(icox.fit)$mFDR
    #icox.true[,i] <- colSums(icox.fit$beta[(t+1):p,] != 0 )/colSums(icox.fit$beta != 0 )
    icox.true.lasso[is.na(icox.true.lasso)] <- 0
    
    ## Fit cox model on cor data
    ccox.fit.lasso <- ncvsurv(XC.cox, yc.cox, returnX = TRUE, lambda = ccox.lseq.lasso, penalty = "lasso")
    cmfcox.temp.lasso  <- mfdr_cox(n = n, L = length(ccox.fit.lasso$lambda), p = p, t = t, alpha = 1, d = ccox.fit.lasso$fail,
                                  lambda = ccox.fit.lasso$lambda, X = ccox.fit.lasso$X, Eta = ccox.fit.lasso$Eta,
                                  S = predict(ccox.fit.lasso, type = 'nvars'))
    cmfcox.temp.lasso[is.infinite(cmfcox.temp.lasso)] = 0
    cmfcox.temp.lasso[cmfcox.temp.lasso > 1] = 1
    if(length(mfcox.temp.lasso) == 100){
      ccox.est.lasso[,i] <- cmfcox.temp.lasso
    } else{ccox.est.lasso[,i] <- c(cmfcox.temp.lasso, rep(1, 100-length(cmfcox.temp.lasso)))} ## If seq stopped early 
    ctfcox.temp.lasso <- colSums(ccox.fit.lasso$beta[(t+1):p,] != 0 )/colSums(ccox.fit.lasso$beta != 0 )
    if(length(ctfcox.temp.lasso) == 100){
      ccox.true.lasso[,i] <- ctfcox.temp.lasso
    } else {ccox.true.lasso[,i] <- c(ctfcox.temp.lasso, rep(1, 100-length(ctfcox.temp.lasso)))}
    ccox.true.lasso[is.na(ccox.true.lasso)] <- 0
    
    ## Fit logistic model on ind data
    ilog.fit.lasso <- ncvreg(X, y.log, family = "binomial", returnX = TRUE, lambda = ilog.lseq.lasso, penalty = "lasso")
    mflog.temp.lasso  <- mfdr_log(n = n, L = length(ilog.fit.lasso$lambda), p = p, t = t, alpha = 1,
                                   lambda = ilog.fit.lasso$lambda, X = ilog.fit.lasso$X, Eta = ilog.fit.lasso$Eta,
                                   S = predict(ilog.fit.lasso, type = 'nvars'))
    mflog.temp.lasso[is.infinite(mflog.temp.lasso)] = 0
    mflog.temp.lasso[mflog.temp.lasso > 1] = 1
    if(length(mflog.temp.lasso) == 100){
      ilog.est.lasso[,i] <- mflog.temp.lasso
    } else{ilog.est.lasso[,i] <- c(mflog.temp.lasso, rep(1, 100-length(mflog.temp.lasso)))} ## If seq stopped early 
    tflog.temp.lasso <- colSums(ilog.fit.lasso$beta[(t+2):(p+1),] != 0 )/colSums(ilog.fit.lasso$beta[-1,] != 0 )
    if(length(tflog.temp.lasso) == 100){
      ilog.true.lasso[,i] <- tflog.temp.lasso
    } else {ilog.true.lasso[,i] <- c(tflog.temp.lasso, rep(1, 100-length(tflog.temp.lasso)))}
    ilog.true.lasso[is.na(ilog.true.lasso)] <- 0
    
    ## Fit logistic model on cor data
    clog.fit.lasso <- ncvreg(XC.log, yc.log, family = "binomial", returnX = TRUE, lambda = clog.lseq.lasso, penalty = "lasso")
    cmflog.temp.lasso  <- mfdr_log(n = n, L = length(clog.fit.lasso$lambda), p = p, t = t, alpha = 1,
                                   lambda = clog.fit.lasso$lambda, X = clog.fit.lasso$X, Eta = clog.fit.lasso$Eta,
                                   S = predict(clog.fit.lasso, type = 'nvars'))
    cmflog.temp.lasso[is.infinite(cmflog.temp.lasso)] = 0
    cmflog.temp.lasso[cmflog.temp.lasso > 1] = 1
    if(length(cmflog.temp.lasso) == 100){
      clog.est.lasso[,i] <- cmflog.temp.lasso
    } else{clog.est.lasso[,i] <- c(cmflog.temp.lasso, rep(1, 100-length(cmflog.temp.lasso)))} ## If seq stopped early 
    ctflog.temp.lasso <- colSums(clog.fit.lasso$beta[(t+2):(p+1),] != 0 )/colSums(clog.fit.lasso$beta[-1,] != 0 )
    if(length(ctflog.temp.lasso) == 100){
      clog.true.lasso[,i] <- ctflog.temp.lasso
    } else {clog.true.lasso[,i] <- c(ctflog.temp.lasso, rep(1, 100-length(ctflog.temp.lasso)))}
    clog.true.lasso[is.na(clog.true.lasso)] <- 0
    
    ## Fit lm on ind data
    ilin.fit.lasso <- ncvreg(X, y.lin, returnX = TRUE, lambda = ilin.lseq.lasso, penalty = "lasso")
    ilin.est.lasso[,i] <- mfdr(ilin.fit.lasso)$mFDR
    ilin.true.lasso[,i] <- colSums(ilin.fit.lasso$beta[(t+2):(p+1),] != 0 )/colSums(ilin.fit.lasso$beta[-1,] != 0 )
    ilin.true.lasso[is.na(ilin.true.lasso)] <- 0
    
    ## Fit lm on cor data
    clin.fit.lasso <- ncvreg(XC.lin, yc.lin, returnX = TRUE, lambda = clin.lseq.lasso, penalty = "lasso")
    clin.est.lasso[,i] <- mfdr(clin.fit.lasso)$mFDR
    clin.true.lasso[,i] <- colSums(clin.fit.lasso$beta[(t+2):(p+1),] != 0 )/colSums(clin.fit.lasso$beta[-1,] != 0 )
    clin.true.lasso[is.na(clin.true.lasso)] <- 0
    
  }
  
  ## Aggregate mean mfdr for each lambda value and store for that sample size
  
  ## MCP MODELS ##
  icox.res.est[,a] <- apply(icox.est, 1, mean)
  icox.res.true[,a] <- apply(icox.true, 1, mean)
  
  ccox.res.est[,a] <- apply(ccox.est, 1, mean)
  ccox.res.true[,a] <- apply(ccox.true, 1, mean)
  
  ilog.res.est[,a] <- apply(ilog.est, 1, mean)
  ilog.res.true[,a] <- apply(ilog.true, 1, mean)
  
  clog.res.est[,a] <- apply(clog.est, 1, mean)
  clog.res.true[,a] <- apply(clog.true, 1, mean)
  
  ilin.res.est[,a] <- apply(ilin.est, 1, mean)
  ilin.res.true[,a] <- apply(ilin.true, 1, mean)
  
  clin.res.est[,a] <- apply(clin.est, 1, mean)
  clin.res.true[,a] <- apply(clin.true, 1, mean)
  
  ## LASSO MODELS ##
  icox.res.est.lasso[,a] <- apply(icox.est.lasso, 1, mean)
  icox.res.true.lasso[,a] <- apply(icox.true.lasso, 1, mean)
  
  ccox.res.est.lasso[,a] <- apply(ccox.est.lasso, 1, mean)
  ccox.res.true.lasso[,a] <- apply(ccox.true.lasso, 1, mean)
  
  ilog.res.est.lasso[,a] <- apply(ilog.est.lasso, 1, mean)
  ilog.res.true.lasso[,a] <- apply(ilog.true.lasso, 1, mean)
  
  clog.res.est.lasso[,a] <- apply(clog.est.lasso, 1, mean)
  clog.res.true.lasso[,a] <- apply(clog.true.lasso, 1, mean)
  
  ilin.res.est.lasso[,a] <- apply(ilin.est.lasso, 1, mean)
  ilin.res.true.lasso[,a] <- apply(ilin.true.lasso, 1, mean)
  
  clin.res.est.lasso[,a] <- apply(clin.est.lasso, 1, mean)
  clin.res.true.lasso[,a] <- apply(clin.true.lasso, 1, mean)
  
}

## Process the results into estimated/true mfdr ratio targeting "thr%" mfdr for each value of n
thr <- .2

## MCP MODELS ##
icox.idx <- apply(icox.res.true, 2, findlam, thresh = thr)
icox.estbylam <- ((p-t)/p)*diag(icox.res.est[icox.idx,])

ccox.idx <- apply(ccox.res.true, 2, findlam, thresh = thr)
ccox.estbylam <- ((p-t)/p)*diag(ccox.res.est[ccox.idx,])

ilog.idx <- apply(ilog.res.true, 2, findlam, thresh = thr)
ilog.estbylam <- ((p-t)/p)*diag(ilog.res.est[ilog.idx,])

clog.idx <- apply(clog.res.true, 2, findlam, thresh = thr)
clog.estbylam <- ((p-t)/p)*diag(clog.res.est[clog.idx,])

ilin.idx <- apply(ilin.res.true, 2, findlam, thresh = thr)
ilin.estbylam <- ((p-t)/p)*diag(ilin.res.est[ilin.idx,])

clin.idx <- apply(clin.res.true, 2, findlam, thresh = thr)
clin.estbylam <- ((p-t)/p)*diag(clin.res.est[clin.idx,])

## LASSO MODELS ##
icox.idx.lasso <- apply(icox.res.true.lasso, 2, findlam, thresh = thr)
icox.estbylam.lasso <- ((p-t)/p)*diag(icox.res.est.lasso[icox.idx.lasso,])

ccox.idx.lasso <- apply(ccox.res.true.lasso, 2, findlam, thresh = thr)
ccox.estbylam.lasso <- ((p-t)/p)*diag(ccox.res.est.lasso[ccox.idx.lasso,])

ilog.idx.lasso <- apply(ilog.res.true.lasso, 2, findlam, thresh = thr)
ilog.estbylam.lasso <- ((p-t)/p)*diag(ilog.res.est.lasso[ilog.idx.lasso,])

clog.idx.lasso <- apply(clog.res.true.lasso, 2, findlam, thresh = thr)
clog.estbylam.lasso <- ((p-t)/p)*diag(clog.res.est[clog.idx.lasso,])

ilin.idx.lasso <- apply(ilin.res.true.lasso, 2, findlam, thresh = thr)
ilin.estbylam.lasso <- ((p-t)/p)*diag(ilin.res.est.lasso[ilin.idx.lasso,])

clin.idx.lasso <- apply(clin.res.true.lasso, 2, findlam, thresh = thr)
clin.estbylam.lasso <- ((p-t)/p)*diag(clin.res.est.lasso[clin.idx.lasso,])



## Format all different models and cor settings to plot in ggplot
lii1 <- data.frame(val = ilin.estbylam.lasso, reg = "Linear", ID = "Independent - lasso", n  = nn, Penalty = "lasso", Setting = "Independent")
lii2 <- data.frame(val = ilin.estbylam, reg = "Linear", ID = "Independent - MCP", n  = nn, Penalty = "MCP", Setting = "Independent")

lic1 <- data.frame(val = clin.estbylam.lasso, reg = "Linear", ID = "Correlated - lasso", n = nn, Penalty = "lasso", Setting = "Correlated")
lic2 <- data.frame(val = clin.estbylam, reg = "Linear", ID = "Correlated - MCP", n  = nn, Penalty = "MCP", Setting = "Correlated")

ci1 <- data.frame(val = icox.estbylam.lasso, reg = "Cox", ID = "Independent - lasso", n  = nn, Penalty = "lasso", Setting = "Independent")
ci2 <- data.frame(val = icox.estbylam, reg = "Cox", ID = "Independent - MCP", n  = nn, Penalty = "MCP", Setting = "Independent")

cc1 <- data.frame(val = ccox.estbylam.lasso, reg = "Cox", ID = "Correlated - lasso", n  = nn, Penalty = "lasso", Setting = "Correlated")
cc2 <- data.frame(val = ccox.estbylam, reg = "Cox", ID = "Correlated - MCP", n  = nn, Penalty = "MCP", Setting = "Correlated")

li1 <- data.frame(val = ilog.estbylam.lasso, reg = "Logistic", ID = "Independent - lasso", n  = nn, Penalty = "lasso", Setting = "Independent")
li2 <- data.frame(val = ilog.estbylam, reg = "Logistic", ID = "Independent - MCP", n  = nn, Penalty = "MCP", Setting = "Independent")

lc1 <- data.frame(val = clog.estbylam.lasso, reg = "Logistic", ID = "Correlated - lasso", n  = nn, Penalty = "lasso", Setting = "Correlated")
lc2 <- data.frame(val = clog.estbylam, reg = "Logistic", ID = "Correlated - MCP", n  = nn, Penalty = "MCP", Setting = "Correlated")

plot.mat <- rbind(lc1,lc2,li1,li2,cc1,cc2,ci1,ci2,lic1,lic2,lii1,lii2)

###################################
############## FIGURE 1 ###########
###################################
pp <- ggplot(plot.mat, aes(x=n, y=val, colour = Setting, group=ID, linetype = Penalty))  +
  geom_line(size = 1.5) + geom_hline(yintercept= thr, linetype="dotted", size = 1.2) + 
  scale_x_continuous("n") + scale_y_continuous(name = "Estimated mFDR", limits = c(thr-thr/4,thr*3.3))

pdf("fig1.pdf",width=7.5,height=3)
(pp + facet_grid(. ~ reg))
dev.off()

########################################################################################
#
#    Simulation #2 (Produces Figures 2-4 and Table 1)
#
########################################################################################

### Setup
n <- 400
p <- 1000
t <- 10
b <- seq(0.05, 1.25, by = .05)
nreps <- 100
epv <- 20  #### This var is how many features should be selected in stage 1 of SS, not actual EPV for stage 2

### Objects for storage, rows are replication for a given beta (cols are diff betas)
pv <- rep(NA, p)
pt2.coxA <- pt2.coxB <- pt2.coxC <- matrix(NA, nrow = nreps, ncol = length(b)) 
mfdrsel.coxA <- mfdrsel.coxB <- mfdrsel.coxC <- matrix(NA, nrow = nreps, ncol = length(b)) 
cvsel.coxA <- cvsel.coxB <- cvsel.coxC <- matrix(NA, nrow = nreps, ncol = length(b)) 
uni.coxA <- uni.coxB <- uni.coxC <- matrix(NA, nrow = nreps, ncol = length(b)) 

pt2.logA <- pt2.logB <- pt2.logC <- matrix(NA, nrow = nreps, ncol = length(b)) 
mfdrsel.logA <- mfdrsel.logB <- mfdrsel.logC <- matrix(NA, nrow = nreps, ncol = length(b)) 
cvsel.logA <- cvsel.logB <- cvsel.logC <- matrix(NA, nrow = nreps, ncol = length(b)) 
covt.logA <- covt.logB <- covt.logC <-  matrix(NA, nrow = nreps, ncol = length(b)) 
uni.logA <- uni.logB <- uni.logC <- matrix(NA, nrow = nreps, ncol = length(b)) 


##### Outer loop over beta
for (k in 1:length(b)) {
  bb <- numeric(p)
  bb[(0:(t-1))*10+1] <- c(rep(b[k],t/2),rep(-b[k],t/2))  #### 10 true variables (type A), spaced by 10 spots
  
  id.A <- which(bb != 0)
  id.B <- 1:100
  id.B <- id.B[-id.A]
  id.C <- 101:p
  
  ### Inner loop within a beta
  for (i in 1:nreps){
    ### Generate data
    ## Logistic
    D1 <- genData(n, J=10, J0=10, K=10, K0=1, rho=0, rho.g=0.5, beta=bb[1:100], family = "binomial")  #### 9 correlated vars for each true (type B)
    D2 <- genData(n, p - 100, rho=0.8, beta=0, corr="auto")
    
    y.log <- D1$y
    X.log <- cbind(D1$X, D2$X)
    
    ## Cox
    D1 <- genData(n, J=10, J0=10, K=10, K0=1, rho=0, rho.g=0.5, beta=bb[1:100], family = "survival")  #### 9 correlated vars for each true (type B)
    D2 <- genData(n, p - 100, rho=0.8, beta=0, corr="auto")
    
    y.cox <- Surv(D1$y, rep(1,n))
    X.cox <- cbind(D1$X, D2$X)
    
    ###############################################
    ############# Univariate Testing ##############
    ###############################################
    
    upv.cox <- upv.log <- rep(NA, p)
    for (j in 1:p){
      fit.lm.cox <- coxph(y.cox ~ X.cox[,j])
      upv.cox[j] <- summary(fit.lm.cox)$coefficients[1,5]
      fit.lm.log <- glm(y.log ~ X.log[,j], family = 'binomial')
      upv.log[j] <- summary(fit.lm.log)$coefficients[2,4]
    }
    aupv.cox <- p.adjust(upv.cox, method = "fdr")
    aupv.log <- p.adjust(upv.log, method = "fdr")
    
    uni.coxA[i,k] <- sum(aupv.cox[id.A] < .1)
    uni.coxB[i,k] <- sum(aupv.cox[id.B] < .1)
    uni.coxC[i,k] <- sum(aupv.cox[id.C] < .1)
    
    uni.logA[i,k] <- sum(aupv.log[id.A] < .1)
    uni.logB[i,k] <- sum(aupv.log[id.B] < .1)
    uni.logC[i,k] <- sum(aupv.log[id.C] < .1)
    
    
    #############################################
    ############## Sample spliting ##############
    #############################################
    
    #### Split the data
    idx <- sample(n, n/2)
    
    X1.log <- X.log[idx,]
    X2.log <- X.log[-idx,]
    X1.cox <- X.cox[idx,]
    X2.cox <- X.cox[-idx,]
    y1.log <- y.log[idx]
    y2.log <- y.log[-idx]
    y1.cox <- y.cox[idx]
    y2.cox <- y.cox[-idx]
    
    
    #### Part 1
    cv.cox <- ncvsurv(X1.cox, y1.cox, penalty = 'lasso')
    cv.log <- ncvreg(X1.log, y1.log, family = "binomial", penalty = 'lasso')
    
    ### Part 2 (cox)
    S <- predict(cv.cox, type = "nvars")
    sel.cox <- which(cv.cox$beta[,max(which(S <= epv))] != 0)
    if (length(sel.cox) > 0 ){
      XX <- X2.cox[,sel.cox]
      ols.fit <- coxph(y2.cox ~ as.matrix(XX))
      pv.selected <- 2*pnorm(-abs(ols.fit$coefficients/sqrt(diag(vcov(ols.fit)))))
      pv.selected <- p.adjust(pv.selected, method = "fdr")
      pv[sel.cox] <- pv.selected
      pv[-sel.cox] <- 1
      rpv.cox <- pv
    } else {
      rpv.cox <- rep(1,p)
    } 
    pt2.coxA[i,k] <- sum(rpv.cox[id.A] < .1)
    pt2.coxB[i,k] <- sum(rpv.cox[id.B] < .1)
    pt2.coxC[i,k] <- sum(rpv.cox[id.C] < .1)
    
    ### Part 2 (log)
    S <- predict(cv.log, type = "nvars")
    sel.log <- which(cv.log$beta[-1,max(which(S <= epv))] != 0)
    if (length(sel.log) > 0 ){
      XX <- X2.log[,sel.log]
      ols.fit <- glm(y2.log ~ as.matrix(XX), family = "binomial")
      pv.selected <- 2*pnorm(-abs(ols.fit$coefficients[-1]/sqrt(diag(vcov(ols.fit)))[-1]))
      pv.selected <- p.adjust(pv.selected, method = "fdr")
      pv[sel.log] <- pv.selected
      pv[-sel.log] <- 1
      rpv.log <- pv
    } else {
      rpv.log <- rep(1,p)
    } 
    pt2.logA[i,k] <- sum(rpv.log[id.A] < .1)
    pt2.logB[i,k] <- sum(rpv.log[id.B] < .1)
    pt2.logC[i,k] <- sum(rpv.log[id.C] < .1)
    
    ###################################################
    ##################  MFDR ##########################
    ###################################################
    
    fit.cox <- ncvsurv(X.cox, y.cox, penalty = 'lasso', returnX = TRUE)
    fit.log <- ncvreg(X.log, y.log, family = "binomial", penalty = 'lasso', returnX = TRUE)
    
    step.cox <- max(which(mfdr(fit.cox)[,3] < 0.1))
    step.log <- max(which(mfdr(fit.log)[,3] < 0.1))
    
    mfdrsel.coxA[i,k] <- sum(coef(fit.cox)[id.A, step.cox] != 0)
    mfdrsel.coxB[i,k] <- sum(coef(fit.cox)[id.B, step.cox] != 0)
    mfdrsel.coxC[i,k] <- sum(coef(fit.cox)[id.C, step.cox] != 0)
    
    mfdrsel.logA[i,k] <- sum(coef(fit.log)[id.A+1, step.log] != 0)  ### Need +1 to skip intercept
    mfdrsel.logB[i,k] <- sum(coef(fit.log)[id.B+1, step.log] != 0)
    mfdrsel.logC[i,k] <- sum(coef(fit.log)[id.C+1, step.log] != 0)
    
    
    ###################################################
    ############# Cross Validation ####################
    ###################################################
    cv.cox <- cv.ncvsurv(X.cox, y.cox, penalty = 'lasso')
    cv.log <- cv.ncvreg(X.log, y.log, family = "binomial", penalty = 'lasso')
    
    cvsel.coxA[i,k] <- sum(coef(cv.cox, s = cv.cox$lambda.min)[id.A] != 0)
    cvsel.coxB[i,k] <- sum(coef(cv.cox, s = cv.cox$lambda.min)[id.B] != 0)
    cvsel.coxC[i,k] <- sum(coef(cv.cox, s = cv.cox$lambda.min)[id.C] != 0)
    
    cvsel.logA[i,k] <- sum(coef(cv.log, s = cv.log$lambda.min)[id.A+1] != 0)
    cvsel.logB[i,k] <- sum(coef(cv.log, s = cv.log$lambda.min)[id.B+1] != 0)
    cvsel.logC[i,k] <- sum(coef(cv.log, s = cv.log$lambda.min)[id.C+1] != 0)
    
    
    ###################################################
    ############# Selective Inference #################
    ###################################################
    ### Logistic
    logfit.lars <- lars.glm(as.matrix(X.log), y.log, family = "binomial")
    ctest.log <- tryCatch(covTest(logfit.lars, as.matrix(X.log), y.log), error=function(e) NULL)
    
    if (is.null(ctest.log)){
      covt.logA[i,k] <- NA
      covt.logB[i,k] <- NA
      covt.logC[i,k] <- NA
    } else {
      cov.cut <- forwardStop(ctest.log$results[!is.na(ctest.log$results[,3]),3], .1)
      if (cov.cut < 1){
        covt.logA[i,k] <- 0
        covt.logB[i,k] <- 0
        covt.logC[i,k] <- 0
      } else {
        covt.logA[i,k] <- sum(ctest.log$results[1:cov.cut,1] %in% id.A) 
        covt.logB[i,k] <- sum(ctest.log$results[1:cov.cut,1] %in% id.B)
        covt.logC[i,k] <- sum(ctest.log$results[1:cov.cut,1] %in% id.C) 
      }
    }
  }
}



plot.power <- rbind(data.frame(beta = b,val = colMeans(pt2.coxA),Method = "SampleSplitting", Type = "A", reg = "Cox"),
                    data.frame(beta = b,val = colMeans(mfdrsel.coxA),Method = "mFDR", Type = "A", reg = "Cox"),
                    data.frame(beta = b,val = colMeans(pt2.coxB),Method = "SampleSplitting", Type = "B", reg = "Cox"),
                    data.frame(beta = b,val = colMeans(mfdrsel.coxB),Method = "mFDR", Type = "B", reg = "Cox"),
                    data.frame(beta = b,val = colMeans(pt2.coxC),Method = "SampleSplitting", Type = "C", reg = "Cox"),
                    data.frame(beta = b,val = colMeans(mfdrsel.coxC),Method = "mFDR", Type = "C", reg = "Cox"),
                    data.frame(beta = b,val = colMeans(pt2.logA),Method = "SampleSplitting", Type = "A", reg = "Logistic"),
                    data.frame(beta = b,val = colMeans(mfdrsel.logA),Method = "mFDR", Type = "A", reg = "Logistic"),
                    data.frame(beta = b,val = colMeans(pt2.logB),Method = "SampleSplitting", Type = "B", reg = "Logistic"),
                    data.frame(beta = b,val = colMeans(mfdrsel.logB),Method = "mFDR", Type = "B", reg = "Logistic"),
                    data.frame(beta = b,val = colMeans(pt2.logC),Method = "SampleSplitting", Type = "C", reg = "Logistic"),
                    data.frame(beta = b,val = colMeans(mfdrsel.logC),Method = "mFDR", Type = "C", reg = "Logistic"),
                    data.frame(beta = b,val = colMeans(covt.logA, na.rm = TRUE),Method = "CovTest", Type = "A", reg = "Logistic"),
                    data.frame(beta = b,val = colMeans(covt.logB, na.rm = TRUE),Method = "CovTest", Type = "B", reg = "Logistic"),
                    data.frame(beta = b,val = colMeans(covt.logC, na.rm = TRUE),Method = "CovTest", Type = "C", reg = "Logistic"))


plot.cv <- rbind(data.frame(beta = b,val = colMeans(cvsel.coxA),Method = "CV", Type = "A", reg = "Cox"),
                 data.frame(beta = b,val = colMeans(cvsel.coxB),Method = "CV", Type = "B", reg = "Cox"),
                 data.frame(beta = b,val = colMeans(cvsel.coxC),Method = "CV", Type = "C", reg = "Cox"),
                 data.frame(beta = b,val = colMeans(cvsel.logA),Method = "CV", Type = "A", reg = "Logistic"),
                 data.frame(beta = b,val = colMeans(cvsel.logB),Method = "CV", Type = "B", reg = "Logistic"),
                 data.frame(beta = b,val = colMeans(cvsel.logC),Method = "CV", Type = "C", reg = "Logistic"))


plot.uni <- rbind(data.frame(beta = b,val = colMeans(uni.coxA),Method = "univariate", Type = "A", reg = "Cox"),
                  data.frame(beta = b,val = colMeans(mfdrsel.coxA),Method = "mFDR", Type = "A", reg = "Cox"),
                  data.frame(beta = b,val = colMeans(uni.coxB),Method = "univariate", Type = "B", reg = "Cox"),
                  data.frame(beta = b,val = colMeans(mfdrsel.coxB),Method = "mFDR", Type = "B", reg = "Cox"),
                  data.frame(beta = b,val = colMeans(uni.coxC),Method = "univariate", Type = "C", reg = "Cox"),
                  data.frame(beta = b,val = colMeans(mfdrsel.coxC),Method = "mFDR", Type = "C", reg = "Cox"),
                  data.frame(beta = b,val = colMeans(uni.logA),Method = "univariate", Type = "A", reg = "Logistic"),
                  data.frame(beta = b,val = colMeans(mfdrsel.logA),Method = "mFDR", Type = "A", reg = "Logistic"),
                  data.frame(beta = b,val = colMeans(uni.logB),Method = "univariate", Type = "B", reg = "Logistic"),
                  data.frame(beta = b,val = colMeans(mfdrsel.logB),Method = "mFDR", Type = "B", reg = "Logistic"),
                  data.frame(beta = b,val = colMeans(uni.logC),Method = "univariate", Type = "C", reg = "Logistic"),
                  data.frame(beta = b,val = colMeans(mfdrsel.logC),Method = "mFDR", Type = "C", reg = "Logistic"))


###################################
############## FIGURE 2 ###########
###################################
pdf("fig2.pdf",width=6,height=3)
plot.power <- plot.power[complete.cases(plot.power) & plot.power$Type == "A",]
plot.power <- plot.power[(plot.power$reg == "Cox" & plot.power$beta < .8) | (plot.power$reg == "Logistic") ,]
p1 <- ggplot(plot.power, aes(x=beta, y=val, colour = Method, group=Method))  +
  geom_line(size = 1.5) + scale_x_continuous("Signal Strength (b)") + scale_y_continuous('Causal Feature Selections')
p1 + facet_grid(. ~ reg, scales = "free")
dev.off()

###################################
############## FIGURE 3 ###########
###################################
pdf("fig3.pdf",width=6,height=3)
plot.cv <- plot.cv[complete.cases(plot.cv),]
levels(plot.cv$Type) <- c('Causal "A"','Correlated "B"','Noise "C"')
p2 <- ggplot(plot.cv, aes(x=beta, y=val, colour = Type, group=Type))  +  scale_fill_discrete(name="Type") + 
  geom_line(size = 1.5) + scale_x_continuous("Signal Strength (b)", breaks = seq(0,1.2,by=.2)) + 
  scale_y_continuous("Feature Selections", breaks = seq(0,70,by=10)) +
  labs(col = "Feature Type")
p2 + facet_grid( ~ reg, scales = "free")
dev.off()

###################################
############## FIGURE 4 ###########
###################################
pdf("fig4.pdf",width=6.5,height=4.5)
plot.uni <- plot.uni[complete.cases(plot.uni),]
plot.uni <- plot.uni[(plot.uni$reg == "Cox" & plot.uni$beta < .8) | (plot.uni$reg == "Logistic") ,]
levels(plot.uni$Type) <- c('Causal "A"','Correlated "B"','Noise "C"')
levels(plot.uni$Method) <- c("Univariate","mFDR","SampleSplitting","CovTest")
p1 <- ggplot(plot.uni, aes(x=beta, y=val, colour = Method, group=Method))  +
  geom_line(size = 1.5) + scale_x_continuous("Signal Strength (b)") + scale_y_continuous("Feature Selections")
p1 + facet_grid(Type ~ reg, scales = "free")
dev.off()

###################################
############## TABLE 1  ###########
###################################

uni.inrange <- plot.uni[which(plot.uni$b %in% seq(.25, .75, by = .1)),]
uni.inrangeA <- uni.inrange[uni.inrange$Method == "Univariate" & uni.inrange$Type == 'Causal "A"', ]
uni.inrangeC <-  uni.inrange[uni.inrange$Method == "Univariate" & uni.inrange$Type == 'Noise "C"', ]
log.uni.AC <- (uni.inrangeA$val/uni.inrangeC$val)[uni.inrangeA$reg == "Logistic"]
cox.uni.AC <- (uni.inrangeA$val/uni.inrangeC$val)[uni.inrangeA$reg == "Cox"]

mfdr.inrangeA <- uni.inrange[uni.inrange$Method == "mFDR" & uni.inrange$Type == 'Causal "A"', ]
mfdr.inrangeC <-  uni.inrange[uni.inrange$Method == "mFDR" & uni.inrange$Type == 'Noise "C"', ]
log.mfdr.AC <- (mfdr.inrangeA$val/mfdr.inrangeC$val)[mfdr.inrangeA$reg == "Logistic"]
cox.mfdr.AC <- (mfdr.inrangeA$val/mfdr.inrangeC$val)[mfdr.inrangeA$reg == "Cox"]

tab.AC <- data.frame(b = unique(uni.inrange$beta),
                     Univariate_Cox = cox.uni.AC, lasso_mFDR_Cox = cox.mfdr.AC,
                     Univariate_Logistic = log.uni.AC, lasso_mFDR_Logistic = log.mfdr.AC)

pdf("tab1.pdf",width=8,height=4)
grid.table(round(tab.AC, 2))
dev.off()


###############################################################################
###########                                                   #################
###########                  SHEDDEN CASE STUDY               #################
###########                                                   #################
###############################################################################

### Packages
#source("https://bioconductor.org/biocLite.R")
#biocLite("GEOquery")
library(GEOquery)

## Load the Shedden case study data
shedden <- load("shedden.RData")

## Pre-processing on Shedden data
colnames(X) <- fData[,1]
ZZ <- model.matrix(~ factor(Sex) + factor(Race) + factor(AdjChemo) + factor(SmHist) + factor(Margin) + factor(Grade), Z)
D <- cbind(ZZ[,-1], Z$Age ,X) ## Combined environment and gene model matrix
pf <- c(rep(0, ncol(ZZ)),rep(1, ncol(X))) ## Penalty factor setup to allow environmental variables into the model unpenalized

## Fit CV models
set.seed(1234)
cv.shedden.lasso <- cv.ncvsurv(D, S, penalty.factor = pf, penalty = 'lasso', returnX = TRUE)
cv.shedden.mcp <-  cv.ncvsurv(D, S, penalty.factor = pf, returnX = TRUE)


## Apply the mFDR method to each model
mfdr.shedden.lasso <- mfdr(cv.shedden.lasso$fit)
mfdr.shedden.mcp <- mfdr(cv.shedden.mcp$fit)

## Find lambda that limits mFDR < 10%
mfdr.shedden.lasso.lam10 <- cv.shedden.lasso$lambda[max(which(mfdr.shedden.lasso$mFDR < .1))]
mfdr.shedden.mcp.lam10 <- cv.shedden.mcp$lambda[max(which(mfdr.shedden.mcp$mFDR < .1))]


## Find CVE 1se and CV 1se lambda value
lasso.cve.1se <- cv.shedden.lasso$cve[which(cv.shedden.lasso$lambda == cv.shedden.lasso$lambda.min)] + cv.shedden.lasso$cvse[which(cv.shedden.lasso$lambda == cv.shedden.lasso$lambda.min)]
lasso.lam.1se <- cv.shedden.lasso$lambda[which.min(cv.shedden.lasso$cve > lasso.cve.1se)]

mcp.cve.1se <- cv.shedden.mcp$cve[which(cv.shedden.mcp$lambda == cv.shedden.mcp$lambda.min)] + cv.shedden.mcp$cvse[which(cv.shedden.mcp$lambda == cv.shedden.mcp$lambda.min)]
mcp.lam.1se <- cv.shedden.mcp$lambda[which.min(cv.shedden.mcp$cve > mcp.cve.1se)]

## mFDR summaries at each lambda in table
shedden.CV.S <- mfdr.shedden.lasso[which(cv.shedden.lasso$lambda == cv.shedden.lasso$lambda.min),]
shedden.CV1se.S <- mfdr.shedden.lasso[which(cv.shedden.lasso$lambda == lasso.lam.1se),]
shedden.mfdr.S <- mfdr.shedden.lasso[which(cv.shedden.lasso$lambda == mfdr.shedden.lasso.lam10),]

shedden.CV.S.mcp <- mfdr.shedden.mcp[which(cv.shedden.mcp$lambda == cv.shedden.mcp$lambda.min),]
shedden.CV1se.S.mcp <- mfdr.shedden.mcp[which(cv.shedden.mcp$lambda == mcp.lam.1se),]
shedden.mfdr.S.mcp <- mfdr.shedden.mcp[which(cv.shedden.mcp$lambda == mfdr.shedden.mcp.lam10),]

tab.temp <- rbind(shedden.CV.S, shedden.CV1se.S, shedden.mfdr.S, shedden.CV.S.mcp, shedden.CV1se.S.mcp,shedden.mfdr.S.mcp)

## CVE and MCE at each lambda in table
shedden.CV.CVE <- cv.shedden.lasso$cve[which(cv.shedden.lasso$lambda == cv.shedden.lasso$lambda.min)]
shedden.CV1se.CVE <- cv.shedden.lasso$cve[which(cv.shedden.lasso$lambda == lasso.lam.1se)]
shedden.mfdr.CVE <- cv.shedden.lasso$cve[which(cv.shedden.lasso$lambda == mfdr.shedden.lasso.lam10)]

shedden.CV.CVE.mcp <- cv.shedden.mcp$cve[which(cv.shedden.mcp$lambda == cv.shedden.mcp$lambda.min)]
shedden.CV1se.CVE.mcp <- cv.shedden.mcp$cve[which(cv.shedden.mcp$lambda == mcp.lam.1se)]
shedden.mfdr.CVE.mcp <- cv.shedden.mcp$cve[which(cv.shedden.mcp$lambda == mfdr.shedden.mcp.lam10)]

shedden.CVEs <- c(shedden.CV.CVE, shedden.CV1se.CVE, shedden.mfdr.CVE, shedden.CV.CVE.mcp, shedden.CV1se.CVE.mcp, shedden.mfdr.CVE.mcp)

tab.shedden <- data.frame(penalty = c(rep("lasso",3), rep("MCP",3)), 
                          desc = c("CV", "CV(1se)", "mFDR", "CV", "CV(1se)", "mFDR"),
                          lambda = c(cv.shedden.lasso$lambda.min, lasso.lam.1se, mfdr.shedden.lasso.lam10, 
                                     cv.shedden.mcp$lambda.min, mcp.lam.1se, mfdr.shedden.mcp.lam10),
                          EF = tab.temp$EF, S = tab.temp$S, mFDR = tab.temp$mFDR,
                          CVE = shedden.CVEs)


###############################################################################
###########                                                   #################
###########                   SPIRA CASE STUDY                #################
###########                                                   #################
###############################################################################

## Load the Spira case study data
spira <- load("spira-geo.RData")

## Pre-processing on spira data
ys <- pData(geo[[1]])[,1]
y <- as.numeric(regexpr('NOT', ys) == -1)
X2 <- t(exprs(geo[[1]]))  
X <- apply(X2, 2, as.numeric)  ## Make gene expressions numeric

## Fit CV models
set.seed(12345)
cv.spira.lasso <- cv.ncvreg(X, y, family = 'binomial', penalty = 'lasso', returnX = TRUE)
cv.spira.mcp <- cv.ncvreg(X, y, family = 'binomial', returnX = TRUE)

## Apply the mFDR method to each model
mfdr.spira.lasso <- mfdr(cv.spira.lasso$fit)
mfdr.spira.mcp <- mfdr(cv.spira.mcp$fit)

## Find lambda that limits mFDR < 10%
mfdr.spira.lasso.lam10 <- cv.spira.lasso$lambda[max(which(mfdr.spira.lasso$mFDR < .1))]
mfdr.spira.mcp.lam10 <- cv.spira.mcp$lambda[max(which(mfdr.spira.mcp$mFDR < .1))]


## Find CVE 1se and CV 1se lambda value
lasso.cve.1se <- cv.spira.lasso$cve[which(cv.spira.lasso$lambda == cv.spira.lasso$lambda.min)] + cv.spira.lasso$cvse[which(cv.spira.lasso$lambda == cv.spira.lasso$lambda.min)]
lasso.lam.1se <- cv.spira.lasso$lambda[which.min(cv.spira.lasso$cve > lasso.cve.1se)]

mcp.cve.1se <- cv.spira.mcp$cve[which(cv.spira.mcp$lambda == cv.spira.mcp$lambda.min)] + cv.spira.mcp$cvse[which(cv.spira.mcp$lambda == cv.spira.mcp$lambda.min)]
mcp.lam.1se <- cv.spira.mcp$lambda[which.min(cv.spira.mcp$cve > mcp.cve.1se)]

## mFDR summaries at each lambda in table
spira.CV.S <- mfdr.spira.lasso[which(cv.spira.lasso$lambda == cv.spira.lasso$lambda.min),]
spira.CV1se.S <- mfdr.spira.lasso[which(cv.spira.lasso$lambda == lasso.lam.1se),]
spira.mfdr.S <- mfdr.spira.lasso[which(cv.spira.lasso$lambda == mfdr.spira.lasso.lam10),]

spira.CV.S.mcp <- mfdr.spira.mcp[which(cv.spira.mcp$lambda == cv.spira.mcp$lambda.min),]
spira.CV1se.S.mcp <- mfdr.spira.mcp[which(cv.spira.mcp$lambda == mcp.lam.1se),]
spira.mfdr.S.mcp <- mfdr.spira.mcp[which(cv.spira.mcp$lambda == mfdr.spira.mcp.lam10),]

tab.temp <- rbind(spira.CV.S, spira.CV1se.S, spira.mfdr.S, spira.CV.S.mcp, spira.CV1se.S.mcp,spira.mfdr.S.mcp)

## CVE and MCE at each lambda in table
spira.CV.CVE <- cv.spira.lasso$cve[which(cv.spira.lasso$lambda == cv.spira.lasso$lambda.min)]
spira.CV1se.CVE <- cv.spira.lasso$cve[which(cv.spira.lasso$lambda == lasso.lam.1se)]
spira.mfdr.CVE <- cv.spira.lasso$cve[which(cv.spira.lasso$lambda == mfdr.spira.lasso.lam10)]

spira.CV.CVE.mcp <- cv.spira.mcp$cve[which(cv.spira.mcp$lambda == cv.spira.mcp$lambda.min)]
spira.CV1se.CVE.mcp <- cv.spira.mcp$cve[which(cv.spira.mcp$lambda == mcp.lam.1se)]
spira.mfdr.CVE.mcp <- cv.spira.mcp$cve[which(cv.spira.mcp$lambda == mfdr.spira.mcp.lam10)]

spira.CVEs <- c(spira.CV.CVE, spira.CV1se.CVE, spira.mfdr.CVE, spira.CV.CVE.mcp, spira.CV1se.CVE.mcp, spira.mfdr.CVE.mcp)

spira.CV.MCE <- cv.spira.lasso$pe[which(cv.spira.lasso$lambda == cv.spira.lasso$lambda.min)]
spira.CV1se.MCE <- cv.spira.lasso$pe[which(cv.spira.lasso$lambda == lasso.lam.1se)]
spira.mfdr.MCE <- cv.spira.lasso$pe[which(cv.spira.lasso$lambda == mfdr.spira.lasso.lam10)]

spira.CV.MCE.mcp <- cv.spira.mcp$pe[which(cv.spira.mcp$lambda == cv.spira.mcp$lambda.min)]
spira.CV1se.MCE.mcp <- cv.spira.mcp$pe[which(cv.spira.mcp$lambda == mcp.lam.1se)]
spira.mfdr.MCE.mcp <- cv.spira.mcp$pe[which(cv.spira.mcp$lambda == mfdr.spira.mcp.lam10)]

spira.MCEs <- c(spira.CV.MCE, spira.CV1se.MCE, spira.mfdr.MCE, spira.CV.MCE.mcp, spira.CV1se.MCE.mcp, spira.mfdr.MCE.mcp)


tab.spira <- data.frame(penalty = c(rep("lasso",3), rep("MCP",3)), 
                        desc = c("CV", "CV(1se)", "mFDR", "CV", "CV(1se)", "mFDR"),
                        lambda = c(cv.spira.lasso$lambda.min, lasso.lam.1se, mfdr.spira.lasso.lam10, 
                                   cv.spira.mcp$lambda.min, mcp.lam.1se, mfdr.spira.mcp.lam10),
                        EF = tab.temp$EF, S = tab.temp$S, mFDR = tab.temp$mFDR,
                        CVE = spira.CVEs, MCE = spira.MCEs)


###################################
############## FIGURE 5 ###########
###################################

pdf("fig5.pdf",width=6.5,height=3)
par(mfrow = c(1,2))
plot(mfdr.shedden.lasso, type = "EF", log.l = TRUE, xlim = c(-1.7, -2.45), ylim = c(0,70), lwd = 2)
plot(mfdr.shedden.lasso, type = "mFDR", log.l = TRUE, xlim = c(-1.7, -2.45), ylim = c(0,1), selected=FALSE, lwd = 2)
dev.off()

###################################
############## TABLE 2  ###########
###################################

pdf("tab2.pdf",width=15,height=4)
grid.table(cbind("Shedden", tab.shedden, "Spira", tab.spira))
dev.off()
