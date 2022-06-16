library(metafor)
library(tidyverse)
dat <- read_csv("data.csv")

## IV Method (Arcsine Square Root Transformation)
ma <- rma(xi = xi, ni = ni, measure = "PAS", method = "REML", data = dat)
summary(ma)
forest(ma, slab = paste(dat$author, dat$year), 
       atransf = transf.iarcsin, refline = NA,
       ilab = cbind(xi, ni), ilab.xpos = c(-0.25, -0.1),
       xlim = c(-1.2, 1.4), cex = 0.75, mlab = "",
       xlab = "Proportion",
       header = "Authors and Year")
op <- par(cex = 0.75, font = 2)
text(c(-0.25, -0.1), 18, c("e", "n"))
text(-1.2, -1, pos = 4, cex = 1, bquote(paste("RE Model (Q = ",
                                              .(formatC(ma$QE, digits = 2, format = "f")), ", df = ", .(ma$k - ma$p),
                                              ", p = ", .(formatC(ma$QEp, digits = 2, format = "f")), "; ", I^2, " = ",
                                              .(formatC(ma$I2, digits = 1, format = "f")), "%)")))
title("Meta-analysis Using Inverse Variance Method with Arcsine Square Root Transformation")
funnel(ma, atransf = transf.iarcsin, xlab = "Proportion")
title("Funnel Plot")

## Diagnostics and Sensitivity Analysis
plot(influence(ma))
sensdat <- dat[-16, ]
ma.sens <- rma(xi = xi, ni = ni, measure = "PAS", method = "REML", data = sensdat)
summary(ma.sens)

## GLMM with Logit Link
ma.glmm <- rma.glmm(xi = xi, ni = ni, measure = "PLO", data = dat)
summary(ma.glmm)
forest(ma.glmm, slab = paste(dat$author, dat$year), 
       atransf = transf.ilogit, refline = NA,
       ilab = cbind(xi, ni), ilab.xpos = c(-7.5, -6.5),
       xlim = c(-12.5, 2.5), cex = 0.75,
       xlab = "Proportion",
       header = "Authors and Year")
op <- par(cex = 0.75, font = 2)
text(c(-7.5, -6.5), 18, c("e", "n"))
title("Meta-analysis Using Generalized Linear Mixed Model with Logit Link")

## Subgroup Analysis: Time since Tested Positive
newdat <- na.omit(dat)
mlabfun <- function(text, res) {
  list(bquote(paste(.(text),
                    " (Q = ", .(formatC(res$QE, digits = 2, format = "f")),
                    ", df = ", .(res$k - res$p),
                    ", p ", .(metafor:::.pval(res$QEp, digits = 2, showeq = TRUE, sep = " ")), "; ",
                    I^2, " = ", .(formatC(res$I2, digits = 1, format = "f")), "%, ",
                    tau^2, " = ", .(formatC(res$tau2, digits = 2, format = "f")), ")")))}
ma.total <- rma(xi = xi, ni = ni, measure = "PAS", method = "REML", data = newdat)
forest(ma.total, slab = paste(newdat$author, newdat$year), 
       atransf = transf.iarcsin, refline = NA,
       ilab = cbind(xi, ni), ilab.xpos = c(-0.3, -0.1),
       ylim = c(-2.25, 24),
       xlim = c(-2.3, 1.4), cex = 0.75, 
       mlab = mlabfun("RE Model for All Studies", ma.total),
       order = time,
       rows = c(17:20, 3:4, 9:12),
       psize = 1,
       xlab = "Proportion",
       header = "Authors and Year")
op <- par(cex = 0.75, font = 2)
text(c(-0.3, -0.1), 23, c("e", "n"))
par(font = 4)
text(-2.3, c(21, 13, 5), pos = 4, c("< 2 weeks", "2 weeks - 2 months", "> 2 months"))
ma.sub1 <- rma(xi = xi, ni = ni, measure = "PAS", method = "REML", subset = (time == "< 2 weeks"), data = newdat)
ma.sub2 <- rma(xi = xi, ni = ni, measure = "PAS", method = "REML", subset = (time == "2 weeks - 2 months"), data = newdat)
ma.sub3 <- rma(xi = xi, ni = ni, measure = "PAS", method = "REML", subset = (time == "> 2 months"), data = newdat)
par(op)
addpoly(ma.sub1, row = 15.5, mlab = mlabfun("RE Model for Subgroup", ma.sub1))
addpoly(ma.sub2, row = 7.5, mlab = mlabfun("RE Model for Subgroup", ma.sub2))
addpoly(ma.sub3, row = 1.5, mlab = mlabfun("RE Model for Subgroup", ma.sub3))
ma.sub <- rma(xi = xi, ni = ni, mods = ~ time, measure = "PAS", method = "REML", data=newdat)
summary(ma.sub)
text(-2.3, -2.25, pos = 4, cex = 0.75, bquote(paste("Test for Subgroup Differences: ",
                                                    Q[M], " = ", .(formatC(ma.sub$QM, digits = 2, format = "f")), ", df = ",
                                                    .(ma.sub$p - 1),", p = ", .(formatC(ma.sub$QMp, digits = 2, format = "f")))))
title("Subgroup Analysis: Time since Tested Positive")

## Meta-regression: HDI
mr <- rma(xi = xi, ni = ni, mods = ~ hdi, measure = "PAS", method = "REML", data = dat)
summary(mr)
regplot(mr, atransf = transf.iarcsin, xlab = "HDI", ylab = "Proportion")
title("Bubble Plot: HDI")