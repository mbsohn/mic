\name{mic}
\alias{mic}
\title{
A multiple imputation method for compositional microbiome data
}
\description{
Generate multiple imputed compositional data
}
\usage{
mic(M, G=NULL, wt=NULL, nimp=100, pdenom=0.5, imp.model="mean",
    sequential=FALSE, cutoff=0.1, pmin=2)
}
\arguments{
  \item{M}{a matrix of compositional data}
  \item{G}{a vector of group membership}
  \item{wt}{a weight for each sample}
  \item{nimp}{the number of imputations}
  \item{pdenom}{the proportion of taxa whose sum is used for denominator in ratio}
  \item{imp.model}{currently available imputation models: mean or linear}
  \item{sequential}{whether or not sequentially impute components}
  \item{cutoff}{the minimum proportion of zeros in a component not to be included in imputation}
  \item{pmin}{a factor determining a small value that replaces zeros}
}
\value{
  A list of imputed compositional matrices
}
\references{
  Sohn, M.B., Gill, S.R., and Scheible, K. \emph{A multiple imputation method for compositional microbiome data} (Submitted)
}
\author{
  Michael B. Sohn

  Maintainer: Michael B. Sohn <michael_sohn@urmc.rochester.edu>
}
\examples{
\dontrun{
library(tidyverse); library(emdbook)
# Simulate a dataset using a zero-inflated negative binomial model
n.smpl <- 100; n.taxa <- 200; prop.str.zeros <- 0.5
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
               zprob = prop.str.zeros) \%>\%
     matrix(., nrow=n.smpl, ncol=n.taxa)

# Run mic
rslt <- mic(M, nimp=10)
}
}
