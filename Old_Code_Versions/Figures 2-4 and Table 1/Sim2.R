#install.packages("devtools")
#library(devtools)
#install_github("pbreheny/ncvreg")
library(ncvreg)
library(survival)
#install.packages('covTest')
library(covTest)
#install.packages('selectiveInference')
library(selectiveInference)


### Function used to generate data
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



### Fixed parms
n <- 400
p <- 1000
t <- 10
b <- seq(0.05, 1.25, by = .05)
nreps <- 300
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
    
    ### Cox  --- Not yet implemented by covTest
    #coxfit.lars <- lars.glm(as.matrix(X), y.cox[,1], status = y.cox[,2], family = "cox")
    #ctest.cox <- tryCatch(covTest(coxfit.lars, as.matrix(X), y.cox[,1], status = y.cox[,2]), error=function(e) NULL)
    
    #if (is.null(ctest.cox)){
    #  covt.log[i,k] <- NA
    #  covf.log[i,k] <- NA
    #} else {
    #  cov.cut.cox <- forwardStop(ctest.cox$results[!is.na(ctest.cox$results[,3]),3], .1)
    #  if (cov.cut.cox < 1){
    #    covt.log[i,k] <- 0
    #    covf.log[i,k] <- 0
    #  } else {
    #    covt.log[i,k] <- sum(ctest.cox$results[1:cov.cut.cox,1] %in% 1:t) 
    #    covf.log[i,k] <- sum(ctest.cox$results[1:cov.cut.cox,1] %in% (t+1):p) 
    #  }
    #}
    
    
  }
  
  #### Save before moving to next beta value
  save(b, t, uni.coxA, uni.coxB, uni.coxC,
       pt2.coxA, pt2.coxB, pt2.coxC,
       mfdrsel.coxA, mfdrsel.coxB, mfdrsel.coxC,
       cvsel.coxA, cvsel.coxB, cvsel.coxC,
       uni.logA, uni.logB, uni.logC,
       pt2.logA, pt2.logB, pt2.logC, 
       mfdrsel.logA, mfdrsel.logB, mfdrsel.logC,
       cvsel.logA, cvsel.logB, cvsel.logC,
       covt.logA, covt.logB, covt.logC, file = "sim2.RData")
}
