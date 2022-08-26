########################################################################
################ SIMULATION 1 -- GENERATE FIGURE 1 #####################
########################################################################
## This file is used to generate Figure 1
## The six .R code files must be run first to obtain simulation results,
## alternatively you can use the appropriate .RData files.
########################################################################

#install.packages(ggplot2) #Required package for plot
library(ggplot2)

wd1 <- getwd()
wd2 <- paste0(wd1,"/Data1")
setwd(wd2)

## Load all of the sim results for linear, logistic, etc.
load(file = "sim1_linear_ind.RData")
lii1 <- data.frame(val = result.lasso, reg = "Linear", ID = "Independent - lasso", n  = nn, Penalty = "lasso", Setting = "Independent")
lii2 <- data.frame(val = result, reg = "Linear", ID = "Independent - MCP", n  = nn, Penalty = "MCP", Setting = "Independent")

load(file = "sim1_linear_cor.RData")
lic1 <- data.frame(val = result.lasso, reg = "Linear", ID = "Correlated - lasso", n = nn, Penalty = "lasso", Setting = "Correlated")
lic2 <- data.frame(val = result, reg = "Linear", ID = "Correlated - MCP", n  = nn, Penalty = "MCP", Setting = "Correlated")

load(file = "sim1_cox_ind.RData")
ci1 <- data.frame(val = result.lasso, reg = "Cox", ID = "Independent - lasso", n  = nn, Penalty = "lasso", Setting = "Independent")
ci2 <- data.frame(val = result, reg = "Cox", ID = "Independent - MCP", n  = nn, Penalty = "MCP", Setting = "Independent")

load(file = "sim1_cox_cor.RData")
cc1 <- data.frame(val = result.lasso, reg = "Cox", ID = "Correlated - lasso", n  = nn, Penalty = "lasso", Setting = "Correlated")
cc2 <- data.frame(val = result, reg = "Cox", ID = "Correlated - MCP", n  = nn, Penalty = "MCP", Setting = "Correlated")

load(file = "sim1_logistic_ind.RData")
li1 <- data.frame(val = result.lasso, reg = "Logistic", ID = "Independent - lasso", n  = nn, Penalty = "lasso", Setting = "Independent")
li2 <- data.frame(val = result, reg = "Logistic", ID = "Independent - MCP", n  = nn, Penalty = "MCP", Setting = "Independent")

load(file = "sim1_logistic_cor.RData")
lc1 <- data.frame(val = result.lasso, reg = "Logistic", ID = "Correlated - lasso", n  = nn, Penalty = "lasso", Setting = "Correlated")
lc2 <- data.frame(val = result, reg = "Logistic", ID = "Correlated - MCP", n  = nn, Penalty = "MCP", Setting = "Correlated")


### Create plot
plot.mat <- rbind(lc1,lc2,li1,li2,cc1,cc2,ci1,ci2,lic1,lic2,lii1,lii2)
p1 <- ggplot(plot.mat, aes(x=n, y=val, colour = Setting, group=ID, linetype = Penalty))  +
  geom_line(size = 1.5) + geom_hline(yintercept= 0.1, linetype="dotted", size = 1.2) + 
  scale_x_continuous("n") + scale_y_continuous(name = "Estimated mFDR", limits = c(0,.46))

setwd(wd1)
pdf("fig1.pdf",width=7.5,height=3)
(p1 + facet_grid(. ~ reg))
dev.off()