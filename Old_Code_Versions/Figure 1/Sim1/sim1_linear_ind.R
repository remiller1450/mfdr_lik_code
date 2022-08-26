library(ncvreg)
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

#findlamb <- function(tfirmat, fir){  ### Alternative
#  clamb  <- min(which(tfirmat > fir))
#  return(clamb)
#}
findlamb <- function(tfirmat, fir){
  clamb  <- max(which(tfirmat < fir))
  return(clamb)
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
### Generate IID data
X <- matrix(rnorm(n*p),nrow = n, ncol = p)
y <- X %*% beta + rnorm(n, sd = true.sig)



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

save(result, result.lasso, nn, file = "sim1_linear_ind.RData")
}

