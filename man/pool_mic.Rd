\name{pool_mic}
\alias{pool_mic}
\title{
Pool parameter estimates
}
\description{
Combine multiple parameter estimates from imputed data by mic into an overall estimate and variance, incorporating both the within and between imputation variability
}
\usage{
pool_mic(est.mat, se.mat, padj.method="holm")
}
\arguments{
  \item{est.mat}{a matrix of multiple parameter estimates}
  \item{se.mat}{a matrix of multiple parameter estimates' standard errors}
  \item{padj.method}{a method for multiple testing correction}
}
\value{
  A data frame consisting of Wald statistics, p-values, and adjusted p-values
}
\references{
  Sohn, M.B., Scheible, K., and Gill, S.R. \emph{A semi-parametric multiple imputation method for high-sparse, high-dimensional, compositional data} (Submitted)
}
\author{
  Michael B. Sohn

  Maintainer: Michael B. Sohn <michael_sohn@urmc.rochester.edu>
}
\examples{
\dontrun{
library(tidyverse); library(emdbook); library(LinDA)
# Simulate datasets for two groups using zero-inflated negative binomial models
n.smpl <- 50; n.taxa <- 100; prop.str.zeros <- 0.5; da.taxa.prop <- 0.2
n.low.abun <- round(n.taxa*0.6)
n.mid.abun <- round(n.taxa*0.3)
n.high.abun <- n.taxa - n.low.abun - n.mid.abun
sim.mean1 <- sim.mean2 <- c(sample(1:3, n.low.abun, rep=TRUE),
                            sample(4:10, n.mid.abun, rep=TRUE),
                            sample(11:100, n.high.abun, rep=TRUE))
sim.wt <- sample(1:20, size=n.smpl, rep=TRUE)
sim.size <- sample(seq(0.5, 2, 0.1), n.taxa, rep=TRUE)
M1 <- rzinbinom(n = n.smpl*n.taxa,
                mu = as.vector(outer(sim.wt, sim.mean1, "*")),
                size = rep(sim.size, each=n.smpl),
                zprob = prop.str.zeros) %>%
  matrix(., nrow=n.smpl, ncol=n.taxa)
da.indx <- sample(1:n.taxa, round(n.taxa*da.taxa.prop))
sim.mean2[da.indx] <- sim.mean1[da.indx]*sample(c(seq(0.1, 0.5, 0.05), seq(2, 10, 0.5)),
                                                round(n.taxa*da.taxa.prop), rep=T)
M2 <- rzinbinom(n = n.smpl*n.taxa,
                mu = as.vector(outer(sim.wt, sim.mean2, "*")),
                size = rep(sim.size, each=n.smpl),
                zprob = prop.str.zeros) %>%
  matrix(., nrow=n.smpl, ncol=n.taxa)
colnames(M1) <- colnames(M2) <- paste0("T", 1:n.taxa)

# Combine M1 and M2 to construct the final profile
M.f <- rbind(M1, M2)
rownames(M.f) <- paste0("S", 1:nrow(M.f))

# Assign group memberships
y <- c(rep(0, n1), rep(1, n2))

# Get differentially abundant taxa
is.true.da.taxa <- (sim.mean1 != sim.mean2)
true.da.taxa <- paste0("T", which(is.true.da.taxa))

# Construct a design matrix
X <- model.matrix(~y)

# Run mic
rslt <- mic(M.f, X=X, pdenom=0.30, nimp=200)

# Run DA analysis on imputed profiles using LinDA
da.rslt.mic.fc <- da.rslt.mic.fc.se <- list()
for(i in 1:length(rslt)){
  M.tmp <- rslt[[i]] %>% t()
  colnames(M.tmp) <- paste0("S", 1:ncol(M.tmp))
  Meta.dat <- data.frame(sid=colnames(M.tmp), arm=as.factor(y))
  linda.out <- linda(M.tmp, Meta.dat, formula='~arm')
  da.rslt.mic.fc[[i]] <- linda.out$output$arm1$log2FoldChange
  da.rslt.mic.fc.se[[i]] <- linda.out$output$arm1$lfcSE
}

# Pool the estimates
est.mat <- simplify2array(da.rslt.mic.fc)
se.mat <- simplify2array(da.rslt.mic.fc.se)
rownames(est.mat) <- colnames(M.f)
pool.rslt <- pool_mic(est.mat, se.mat)
}
}
