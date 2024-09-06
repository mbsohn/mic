README
================
Michael B. Sohn
09/05/2024

## MIC: <ins>M</ins>ultiple <ins>I</ins>mputation Method for High-Sparse, High-Dimensional, <ins>C</ins>ompositional Data

The *mic* function generates a number of imputed compositional data
using a random selection and amalgamation procedure. For detailed
information about the arguments, please see the documentation for
*mic()*.

### Installation

install.packages(“devtools”)

devtools::install_github(“mbsohn/mic”)

### Example: Generate 200 imputed profiles

``` r
library(mic); library(tidyverse); library(emdbook)
# Simulate datasets for two groups using zero-inflated negative binomial models
set.seed(2024)
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
y <- c(rep(0, n.smpl), rep(1, n.smpl))
# Get differentially abundant taxa
is.true.da.taxa <- (sim.mean1 != sim.mean2)
true.da.taxa <- paste0("T", which(is.true.da.taxa))
# Construct a design matrix
X <- model.matrix(~y)
# Run mic
rslt <- mic(M.f, X=X, pdenom=0.30, nimp=200)
```

### Differential abundance analysis for two groups using LinDA

``` r
library(LinDA)
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
cat("DA taxa: ", true.da.taxa, "\n")
```

    ## DA taxa:  T10 T15 T26 T29 T38 T39 T43 T47 T48 T51 T53 T56 T64 T70 T77 T80 T83 T87 T92 T97

``` r
print(pool.rslt)
```

    ##                W            p         adjp
    ## T1    0.91415577 1.803175e-01 1.000000e+00
    ## T2   -0.77324213 2.196895e-01 1.000000e+00
    ## T3    0.84597425 1.987836e-01 1.000000e+00
    ## T4   -0.80630567 2.100333e-01 1.000000e+00
    ## T5   -0.16962303 4.326533e-01 1.000000e+00
    ## T6   -1.01732868 1.544985e-01 1.000000e+00
    ## T7   -0.17553073 4.303313e-01 1.000000e+00
    ## T8    2.65095793 4.013192e-03 3.250685e-01
    ## T9    1.09323292 1.371458e-01 1.000000e+00
    ## T10   2.63758404 4.174947e-03 3.339957e-01
    ## T11  -1.31471316 9.430316e-02 1.000000e+00
    ## T12  -0.40457137 3.428963e-01 1.000000e+00
    ## T13  -0.44541258 3.280108e-01 1.000000e+00
    ## T14  -0.07602804 4.696984e-01 1.000000e+00
    ## T15   2.87528873 2.018291e-03 1.654999e-01
    ## T16   1.58951323 5.597229e-02 1.000000e+00
    ## T17  -1.54039723 6.173178e-02 1.000000e+00
    ## T18   0.83898665 2.007384e-01 1.000000e+00
    ## T19   2.31185835 1.039275e-02 7.898486e-01
    ## T20   0.49497676 3.103083e-01 1.000000e+00
    ## T21  -0.37908753 3.523114e-01 1.000000e+00
    ## T22  -0.29180447 3.852181e-01 1.000000e+00
    ## T23  -0.43783989 3.307512e-01 1.000000e+00
    ## T24   0.43573465 3.315146e-01 1.000000e+00
    ## T25  -0.17474887 4.306385e-01 1.000000e+00
    ## T26   6.88035032 2.985278e-12 2.776308e-10
    ## T27   0.40689194 3.420437e-01 1.000000e+00
    ## T28   1.00458762 1.575477e-01 1.000000e+00
    ## T29  10.20086086 9.826139e-25 9.826139e-23
    ## T30  -0.71025544 2.387729e-01 1.000000e+00
    ## T31   1.44828278 7.376899e-02 1.000000e+00
    ## T32  -0.71280816 2.379822e-01 1.000000e+00
    ## T33   0.62576450 2.657347e-01 1.000000e+00
    ## T34   1.63195835 5.134413e-02 1.000000e+00
    ## T35  -1.90074489 2.866772e-02 1.000000e+00
    ## T36  -0.36330059 3.581902e-01 1.000000e+00
    ## T37  -0.60110900 2.738837e-01 1.000000e+00
    ## T38  -1.75382540 3.973022e-02 1.000000e+00
    ## T39   5.78245383 3.680939e-09 3.202417e-07
    ## T40   1.58731168 5.622104e-02 1.000000e+00
    ## T41  -1.55626229 5.982287e-02 1.000000e+00
    ## T42   0.05153071 4.794513e-01 1.000000e+00
    ## T43   9.05906359 6.578761e-20 6.381398e-18
    ## T44   2.45247250 7.093911e-03 5.462312e-01
    ## T45  -1.94269807 2.602632e-02 1.000000e+00
    ## T46   0.29589432 3.836554e-01 1.000000e+00
    ## T47   6.07344016 6.259933e-10 5.633940e-08
    ## T48   2.54028909 5.538044e-03 4.319674e-01
    ## T49  -0.26207452 3.966320e-01 1.000000e+00
    ## T50  -0.83487365 2.018944e-01 1.000000e+00
    ## T51   5.90942477 1.716522e-09 1.527705e-07
    ## T52   0.63747898 2.619064e-01 1.000000e+00
    ## T53   4.81474857 7.369267e-07 6.263877e-05
    ## T54   2.96770811 1.500146e-03 1.260122e-01
    ## T55  -0.92481426 1.775313e-01 1.000000e+00
    ## T56   9.11772397 3.835926e-20 3.759208e-18
    ## T57   1.77954683 3.757508e-02 1.000000e+00
    ## T58  -0.50880663 3.054439e-01 1.000000e+00
    ## T59  -0.99681787 1.594265e-01 1.000000e+00
    ## T60  -0.04199612 4.832509e-01 1.000000e+00
    ## T61   0.34950504 3.633551e-01 1.000000e+00
    ## T62   0.96995602 1.660342e-01 1.000000e+00
    ## T63   0.57254457 2.834765e-01 1.000000e+00
    ## T64   8.48846110 1.046957e-17 9.946087e-16
    ## T65   2.54708064 5.431416e-03 4.290819e-01
    ## T66  -1.10661782 1.342296e-01 1.000000e+00
    ## T67  -1.60208504 5.456840e-02 1.000000e+00
    ## T68   0.26336194 3.961358e-01 1.000000e+00
    ## T69   0.83855450 2.008597e-01 1.000000e+00
    ## T70  -6.67847203 1.207231e-11 1.110653e-09
    ## T71   1.47503431 7.010166e-02 1.000000e+00
    ## T72   0.29995501 3.821057e-01 1.000000e+00
    ## T73  -1.44799198 7.380965e-02 1.000000e+00
    ## T74  -0.25649610 3.987839e-01 1.000000e+00
    ## T75   0.73051372 2.325381e-01 1.000000e+00
    ## T76  -0.17686004 4.298092e-01 1.000000e+00
    ## T77   5.35597818 4.254739e-08 3.659075e-06
    ## T78   1.33316331 9.123911e-02 1.000000e+00
    ## T79   0.05817025 4.768065e-01 1.000000e+00
    ## T80  -7.54514096 2.258993e-14 2.123453e-12
    ## T81  -1.22627847 1.100470e-01 1.000000e+00
    ## T82  -0.20815807 4.175528e-01 1.000000e+00
    ## T83   5.88753582 1.959981e-09 1.724783e-07
    ## T84  -0.01630492 4.934956e-01 1.000000e+00
    ## T85  -1.65542535 4.891912e-02 1.000000e+00
    ## T86   0.29002648 3.858980e-01 1.000000e+00
    ## T87   8.61932942 3.367367e-18 3.232673e-16
    ## T88  -0.97053205 1.658907e-01 1.000000e+00
    ## T89   1.15702889 1.236303e-01 1.000000e+00
    ## T90  -1.44902498 7.366530e-02 1.000000e+00
    ## T91  -0.34035707 3.667938e-01 1.000000e+00
    ## T92   9.47687546 1.310055e-21 1.296955e-19
    ## T93   1.39263691 8.186481e-02 1.000000e+00
    ## T94   0.03344306 4.866606e-01 1.000000e+00
    ## T95  -1.97599217 2.407783e-02 1.000000e+00
    ## T96  -0.65302185 2.568711e-01 1.000000e+00
    ## T97   6.29985038 1.489665e-10 1.355596e-08
    ## T98   0.55725565 2.886764e-01 1.000000e+00
    ## T99   2.94445934 1.617597e-03 1.342606e-01
    ## T100 -2.29622082 1.083163e-02 8.123722e-01
