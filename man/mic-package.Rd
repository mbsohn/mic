\name{mic-package}
\alias{mic-package}
\docType{package}
\title{
A semi-parametric multiple imputation method for high-sparse, high-dimensional, compositional data
}
\description{
Generate multiple imputed compositional data
}
\author{
Michael B. Sohn

Maintainer: Michael B. Sohn <michael_sohn@urmc.rochester.edu>
}
\references{
Sohn, M.B., Scheible, K., and Gill, S.R. \emph{A semi-parametric multiple imputation method for high-sparse, high-dimensional, compositional data} (Submitted)
}
\examples{
\dontrun{
library(tidyverse); library(emdbook)
# Simulate a dataset using a zero-inflated negative binomial model
n.smpl <- 100; n.taxa <- 200
n.low.abun <- round(n.taxa*0.6)
n.mid.abun <- round(n.taxa*0.3)
n.high.abun <- n.taxa - n.low.abun - n.mid.abun
sim.mean <- c(sample(1:3, n.low.abun, rep=TRUE),
              sample(4:10, n.mid.abun, rep=TRUE),
              sample(11:100, n.high.abun, rep=TRUE))
sim.wt <- sample(1:20, size=n.smpl, rep=TRUE)
sim.size <- sample(seq(0.5, 2, 0.1), n.taxa, rep=TRUE)
M <- rzinbinom(n = n.smpl*n.taxa,
               mu = as.vector(outer(sim.wt, sim.mean, "*")),
               size = rep(sim.size, each=n.smpl),
               zprob = prop.str.zeros) %>%
     matrix(., nrow=n.smpl, ncol=n.taxa)

# Run mic
rslt <- mic(M, nimp=10)
}
}
