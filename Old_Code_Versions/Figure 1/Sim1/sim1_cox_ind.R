library(ncvreg)
library(survival)


fir.surv <- function(fit){
  X <- fit$X
  d <- fit$fail
  S <- predict(fit, type="nvars")
  l <- fit$lambda*fit$alpha
  pf <- which(fit$penalty.factor != 0)
  p <- 40
  t <- 4
  EF <- WW <- NULL
  for (i in 1:length(l)){
    ### Calculate W matrix
    R <- rev(cumsum(rev(exp(fit$Eta[,i]))))
    Rk <- 1/R
    Rk2 <- 1/(R^2)
    RR <- cumsum(Rk)
    RR2 <- cumsum(Rk2)
    W <- diag(exp(fit$Eta[,i])*RR - (exp(fit$Eta[,i])^2)*RR2)
    
    ### Calculate x'Wx
    for (j in 1:ncol(fit$X)){
      WW[j] <- t(X[,j]) %*% W %*% X[,j]
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

nn <- c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)

result <- result2 <- f.est <- lambind <- NULL
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
    X <- matrix(rnorm(n*p), nrow = n, ncol = p)
    B <- c(rep(5/sqrt(n), t/2), rep(-5/sqrt(n), t/2), rep(0, p - t))
    y <- rexp(n, rate = exp(X %*% B))  ### proportional hazards, everyone with a constant baseline hazard of 0.001
    nc <- rep(1, n)  ### Use for no censoring
    yy <- Surv(y, nc)
    ### Fit Cox Model
    fit.lam <- ncvsurv(X, yy, returnX = TRUE, penalty = 'MCP', nlambda = 3)
    lmax <- c(lmax, max(fit.lam$lambda))
  }
  lseq <- seq(1.2*max(lmax),0.05*max(lmax), length.out = 200)
  lseq.lasso <- seq(1.2*max(lmax), 0.05*max(lmax), length.out = 200)
  
  for (k in 1:nsim){
    X <- matrix(rnorm(n*p), nrow = n, ncol = p)
    B <- c(rep(5/sqrt(n), t/2), rep(-5/sqrt(n), t/2), rep(0, p - t))
    y <- rexp(n, rate = exp(X %*% B))  ### proportional hazards, everyone with a constant baseline hazard of 0.001
    nc <- rep(1, n)  ### Use for no censoring
    yy <- Surv(y, nc)
    
    ### Fit Cox Model
    fit <- ncvsurv(X, yy, returnX = TRUE, penalty = 'MCP', lambda = lseq)
    fit.lasso <- ncvsurv(X, yy, returnX = TRUE, penalty = 'lasso', lambda = lseq.lasso)
    
    
    ### Number of selections
    S <- apply(fit$beta, 2, nonzero)
    S.lasso <- apply(fit.lasso$beta, 2, nonzero)
    
    
    ### Calcutate true FIR
    tfd1 = apply(fit$beta[(t+1):p,], 2, nonzero)  ### True number of false inclusions
    ti1 = apply(fit$beta[1:t,], 2, nonzero)       ### Number of correct inclusions 

    
    tfd1.lasso = apply(fit.lasso$beta[(t+1):p,], 2, nonzero)  ### True number of false inclusions
    ti1.lasso = apply(fit.lasso$beta[1:t,], 2, nonzero)       ### Number of correct inclusions 

    
    ### Estimate FIR
    FI <- fir.surv(fit)
    FI.lasso <- fir.surv(fit.lasso)
    
    
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
  lambind <- c(lambind, findlamb(mres.est, fir = .1))
  
  save(result, result.lasso, nn, file = "sim1_cox_ind.RData")
}
