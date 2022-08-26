########################################################################
#################### GENERATE FIGURES 2,3,4 ############################
########################################################################
## This file is used to generate Figures 2,3,4
## The .R code files must be run first to obtain simulation results,
## alternatively you can load the provided .RData file
########################################################################
#install.packages("ggplot2")
library(ggplot2)

#### Note - make sure sim.RData is in the current working directly
load(file = "sim2.RData")

plot.power <- rbind(data.frame(beta = b,val = colMeans(pt2.coxA),Method = "SampleSplitting", Type = "A", reg = "Cox"),
                    #data.frame(beta = b,val = colMeans(uni.coxA),Method = "univariate", Type = "A", reg = "Cox"),
                    data.frame(beta = b,val = colMeans(mfdrsel.coxA),Method = "mFDR", Type = "A", reg = "Cox"),
                    data.frame(beta = b,val = colMeans(pt2.coxB),Method = "SampleSplitting", Type = "B", reg = "Cox"),
                    #data.frame(beta = b,val = colMeans(uni.coxB),Method = "univariate", Type = "B", reg = "Cox"),
                    data.frame(beta = b,val = colMeans(mfdrsel.coxB),Method = "mFDR", Type = "B", reg = "Cox"),
                    data.frame(beta = b,val = colMeans(pt2.coxC),Method = "SampleSplitting", Type = "C", reg = "Cox"),
                    #data.frame(beta = b,val = colMeans(uni.coxC),Method = "univariate", Type = "C", reg = "Cox"),
                    data.frame(beta = b,val = colMeans(mfdrsel.coxC),Method = "mFDR", Type = "C", reg = "Cox"),
                    data.frame(beta = b,val = colMeans(pt2.logA),Method = "SampleSplitting", Type = "A", reg = "Logistic"),
                    #data.frame(beta = b,val = colMeans(uni.logA),Method = "univariate", Type = "A", reg = "Logistic"),
                    data.frame(beta = b,val = colMeans(mfdrsel.logA),Method = "mFDR", Type = "A", reg = "Logistic"),
                    data.frame(beta = b,val = colMeans(pt2.logB),Method = "SampleSplitting", Type = "B", reg = "Logistic"),
                    #data.frame(beta = b,val = colMeans(uni.logB),Method = "univariate", Type = "B", reg = "Logistic"),
                    data.frame(beta = b,val = colMeans(mfdrsel.logB),Method = "mFDR", Type = "B", reg = "Logistic"),
                    data.frame(beta = b,val = colMeans(pt2.logC),Method = "SampleSplitting", Type = "C", reg = "Logistic"),
                    #data.frame(beta = b,val = colMeans(uni.logC),Method = "univariate", Type = "C", reg = "Logistic"),
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
# Power comparison of mFDR, SS, CovTest
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
# Failures of CV
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
###### Comparison with univariate testing
pdf("fig4.pdf",width=6.5,height=4.5)
plot.uni <- plot.uni[complete.cases(plot.uni),]
plot.uni <- plot.uni[(plot.uni$reg == "Cox" & plot.uni$beta < .8) | (plot.uni$reg == "Logistic") ,]
levels(plot.uni$Type) <- c('Causal "A"','Correlated "B"','Noise "C"')
levels(plot.uni$Method) <- c("Univariate","mFDR","SampleSplitting","CovTest")
p1 <- ggplot(plot.uni, aes(x=beta, y=val, colour = Method, group=Method))  +
  geom_line(size = 1.5) + scale_x_continuous("Signal Strength (b)") + scale_y_continuous("Feature Selections")
p1 + facet_grid(Type ~ reg, scales = "free")
dev.off()
