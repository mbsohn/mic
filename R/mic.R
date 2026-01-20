mic <- function(M, G=NULL, wt=NULL, nimp=100, pdenom=0.5, imp.model="mean",
                sequential=FALSE, cutoff=0.1, pmin=2){
  if(is.null(wt)) wt <- rep(1, nrow(M))
  if(is.null(G)) G <- rep(1, nrow(M))
  n.G <- length(unique(G))
  M <- as.data.frame(M)
  min.nonzero <- min(M[M>0])
  if(min.nonzero > 0.999){
    subzero <- 0.5
  } else{
    subzero <- min.nonzero/pmin
  }
  split_M <- split(M, G)
  split_wt <- split(wt, G)
  split_miM <- list()
  for(gr.indx in 1:n.G){
    M1ex.indx <- which(!apply(split_M[[gr.indx]], 2, function(x) sum(x>0)/length(x)) > cutoff)
    if(length(M1ex.indx)>0){
      M1ex <- (split_M[[gr.indx]])[,M1ex.indx]
      M1ex[M1ex == 0] <- subzero
      M1 <- (split_M[[gr.indx]])[,-M1ex.indx]
      M1 <- M1[,names(sort(apply(M1, 2, function(x) sum(x>0)), decreasing=TRUE))]
      wt1 <- split_wt[[gr.indx]]
    } else{
      M1ex <- NULL
      M1 <- split_M[[gr.indx]]
      M1 <- M1[,names(sort(apply(M1, 2, function(x) sum(x>0)), decreasing=TRUE))]
      wt1 <- split_wt[[gr.indx]]
    }
    if(imp.model == "mean"){
      if(sequential == TRUE){
        split_miM[[gr.indx]] <- mic.core.s.0(M1, M1ex, wt1, nimp, pdenom)
      } else{
        split_miM[[gr.indx]] <- mic.core.n.0(M1, M1ex, wt1, nimp, pdenom)
      }
    } else if(imp.model == "linear"){
      if(sequential == TRUE){
        split_miM[[gr.indx]] <- mic.core.s.1(M1, M1ex, wt1, nimp, pdenom)
      } else{
        split_miM[[gr.indx]] <- mic.core.n.1(M1, M1ex, wt1, nimp, pdenom)
      }
    } else{
      stop("Only the mean or simple linear model is currently supported!")
    }
  }
  names(split_miM) <- paste0("Group", 1:n.G)
  return(split_miM)
}

### MIC Core for Non-Sequential Imputation: Mean Model
mic.core.n.0 <- function(M1, M1ex, wt1, nimp, pdenom){
  n.cores <- detectCores(); mic.cl <- makeCluster(n.cores-1)
  registerDoParallel(cl=mic.cl)
  n.taxa <- ncol(M1)
  miM1 <- list()
  rbind.M1 <- foreach(i=1:nimp, .combine="rbind") %dopar% {
    tmp.M1 <- M1
    tx.indx <- 1:n.taxa
    for(j in 1:n.taxa){
      ran.den.indx <- sample(tx.indx[-j], round(n.taxa*pdenom))
      m1j <- M1[,j]
      mdj <- rowSums(M1[,ran.den.indx])
      rj <- m1j/mdj
      rj.zero.indx <- which(rj==0 | is.infinite(rj))
      if(length(rj.zero.indx)>0){
        yj <- log(rj[-rj.zero.indx])
        wj <- wt1[-rj.zero.indx]
        ln_rj <- weighted.mean(yj, w=wj)
        m1j[rj.zero.indx] <- exp(ln_rj) * mdj[rj.zero.indx]
        tmp.M1[,j] <- ifelse(m1j==0, NA, m1j)
      }
    }
    if(!is.null(M1ex)){
      tmp.M1 <- cbind(tmp.M1, M1ex)
    }
    tmp.M1 <- sweep(tmp.M1, 1, rowSums(tmp.M1), "/")
    tmp.M1 <- as.matrix(tmp.M1[,order(colnames(tmp.M1))])
  }
  stopCluster(cl=mic.cl)
  miM1 <- lapply(split(rbind.M1, rep(c(1:nimp), each=(nrow(rbind.M1)/nimp))), matrix, nrow(rbind.M1)/nimp)
  miM1 <- lapply(miM1, function(x) {colnames(x) <- colnames(rbind.M1); x})
  na.indx <- which(sapply(miM1, function(x) any(is.na(x))) == TRUE)
  l.na.indx <- length(na.indx)
  if(l.na.indx == nimp){
    stop("Too many zeros!!! Try to increase the value of parameter `pdenom'.")
  }
  if(l.na.indx > 0){
    miM1 <- miM1[-na.indx]
  }
  names(miM1) <- paste0("ImpProfile", 1:length(miM1))
  class(miM1) <- "mic"
  return(miM1)
}

### MIC Core for Non-Sequential Imputation: Simple Regression Model
mic.core.n.1 <- function(M1, M1ex, wt1, nimp, pdenom){
  n.cores <- detectCores(); mic.cl <- makeCluster(n.cores-1)
  registerDoParallel(cl=mic.cl)
  n.taxa <- ncol(M1)
  miM1 <- list()
  rbind.M1 <- foreach(i=1:nimp, .combine="rbind") %dopar% {
    tmp.M1 <- M1
    tx.indx <- 1:n.taxa
    for(j in 1:n.taxa){
      ran.den.indx <- sample(tx.indx[-j], round(n.taxa*pdenom))
      m1j <- M1[,j]
      mdj <- rowSums(M1[,ran.den.indx])
      mnj <- rowSums(M1[,-c(j, ran.den.indx)])
      rj <- m1j/mdj
      rxj <- mnj/mdj
      rj.zero.indx <- which(rj==0 | is.infinite(rj))
      if(length(rj.zero.indx)>0){
        yj <- log(rj[-rj.zero.indx])
        xj <- log(rxj[-rj.zero.indx])
        xj[is.infinite(xj)] <- NA
        wj <- wt1[-rj.zero.indx]
        lm_ln_rj <- lm(yj ~ xj, weights=wj, na.action=na.omit)
        ln_rxj <- data.frame(xj=log(rxj[rj.zero.indx]))
        ln_rj <- predict(lm_ln_rj, newdata=ln_rxj)
        ln_rj[is.infinite(ln_rj)] <- weighted.mean(yj, w=wj)
        m1j[rj.zero.indx] <- exp(ln_rj) * mdj[rj.zero.indx]
        tmp.M1[,j] <- ifelse(m1j==0, NA, m1j)
      }
    }
    if(!is.null(M1ex)){
      tmp.M1 <- cbind(tmp.M1, M1ex)
    }
    tmp.M1 <- sweep(tmp.M1, 1, rowSums(tmp.M1), "/")
    tmp.M1 <- as.matrix(tmp.M1[,order(colnames(tmp.M1))])
  }
  stopCluster(cl=mic.cl)
  miM1 <- lapply(split(rbind.M1, rep(c(1:nimp), each=(nrow(rbind.M1)/nimp))), matrix, nrow(rbind.M1)/nimp)
  miM1 <- lapply(miM1, function(x) {colnames(x) <- colnames(rbind.M1); x})
  na.indx <- which(sapply(miM1, function(x) any(is.na(x))) == TRUE)
  l.na.indx <- length(na.indx)
  if(l.na.indx == nimp){
    stop("Too many zeros!!! Try to increase the value of parameter `pdenom'.")
  }
  if(l.na.indx > 0){
    miM1 <- miM1[-na.indx]
  }
  names(miM1) <- paste0("ImpProfile", 1:length(miM1))
  class(miM1) <- "mic"
  return(miM1)
}

### MIC Core for Sequential Imputation: Mean Model
mic.core.s.0 <- function(M1, M1ex, wt1, nimp, pdenom){
  n.cores <- detectCores(); mic.cl <- makeCluster(n.cores-1)
  registerDoParallel(cl=mic.cl)
  n.taxa <- ncol(M1)
  miM1 <- list()
  rbind.M1 <- foreach(i=1:nimp, .combine="rbind") %dopar% {
    tmp.M1 <- M1
    tx.indx <- 1:n.taxa
    for(j in 1:n.taxa){
      ran.den.indx <- sample(tx.indx[-j], round(n.taxa*pdenom))
      m1j <- tmp.M1[,j]
      mdj <- rowSums(tmp.M1[,ran.den.indx])
      rj <- m1j/mdj
      rj.zero.indx <- which(rj==0 | is.infinite(rj))
      if(length(rj.zero.indx)>0){
        yj <- log(rj[-rj.zero.indx])
        wj <- wt1[-rj.zero.indx]
        ln_rj <- weighted.mean(yj, w=wj)
        m1j[rj.zero.indx] <- exp(ln_rj) * mdj[rj.zero.indx]
        tmp.M1[,j] <- ifelse(m1j==0, NA, m1j)
      }
    }
    if(!is.null(M1ex)){
      tmp.M1 <- cbind(tmp.M1, M1ex)
    }
    tmp.M1 <- sweep(tmp.M1, 1, rowSums(tmp.M1), "/")
    tmp.M1 <- as.matrix(tmp.M1[,order(colnames(tmp.M1))])
  }
  stopCluster(cl=mic.cl)
  miM1 <- lapply(split(rbind.M1, rep(c(1:nimp), each=(nrow(rbind.M1)/nimp))), matrix, nrow(rbind.M1)/nimp)
  miM1 <- lapply(miM1, function(x) {colnames(x) <- colnames(rbind.M1); x})
  na.indx <- which(sapply(miM1, function(x) any(is.na(x))) == TRUE)
  l.na.indx <- length(na.indx)
  if(l.na.indx == nimp){
    stop("Too many zeros!!! Try to increase the value of parameter `pdenom'.")
  }
  if(l.na.indx > 0){
    miM1 <- miM1[-na.indx]
  }
  names(miM1) <- paste0("ImpProfile", 1:length(miM1))
  class(miM1) <- "mic"
  return(miM1)
}

### MIC Core for Sequential Imputation: Simple Regression Model
mic.core.s.1 <- function(M1, M1ex, wt1, nimp, pdenom){
  n.cores <- detectCores(); mic.cl <- makeCluster(n.cores-1)
  registerDoParallel(cl=mic.cl)
  n.taxa <- ncol(M1)
  miM1 <- list()
  rbind.M1 <- foreach(i=1:nimp, .combine="rbind") %dopar% {
    tmp.M1 <- M1
    tx.indx <- 1:n.taxa
    for(j in 1:n.taxa){
      ran.den.indx <- sample(tx.indx[-j], round(n.taxa*pdenom))
      m1j <- tmp.M1[,j]
      mdj <- rowSums(tmp.M1[,ran.den.indx])
      mnj <- rowSums(tmp.M1[,-c(j, ran.den.indx)])
      rj <- m1j/mdj
      rxj <- mnj/mdj
      rj.zero.indx <- which(rj==0 | is.infinite(rj))
      if(length(rj.zero.indx)>0){
        yj <- log(rj[-rj.zero.indx])
        xj <- log(rxj[-rj.zero.indx])
        xj[is.infinite(xj)] <- NA
        wj <- wt1[-rj.zero.indx]
        lm_ln_rj <- lm(yj ~ xj, weights=wj, na.action=na.omit)
        ln_rxj <- data.frame(xj=log(rxj[rj.zero.indx]))
        ln_rj <- predict(lm_ln_rj, newdata=ln_rxj)
        ln_rj[is.infinite(ln_rj)] <- weighted.mean(yj, w=wj)
        m1j[rj.zero.indx] <- exp(ln_rj) * mdj[rj.zero.indx]
        tmp.M1[,j] <- ifelse(m1j==0, NA, m1j)
      }
    }
    if(!is.null(M1ex)){
      tmp.M1 <- cbind(tmp.M1, M1ex)
    }
    tmp.M1 <- sweep(tmp.M1, 1, rowSums(tmp.M1), "/")
    tmp.M1 <- as.matrix(tmp.M1[,order(colnames(tmp.M1))])
  }
  stopCluster(cl=mic.cl)
  miM1 <- lapply(split(rbind.M1, rep(c(1:nimp), each=(nrow(rbind.M1)/nimp))), matrix, nrow(rbind.M1)/nimp)
  miM1 <- lapply(miM1, function(x) {colnames(x) <- colnames(rbind.M1); x})
  na.indx <- which(sapply(miM1, function(x) any(is.na(x))) == TRUE)
  l.na.indx <- length(na.indx)
  if(l.na.indx == nimp){
    stop("Too many zeros!!! Try to increase the value of parameter `pdenom'.")
  }
  if(l.na.indx > 0){
    miM1 <- miM1[-na.indx]
  }
  names(miM1) <- paste0("ImpProfile", 1:length(miM1))
  class(miM1) <- "mic"
  return(miM1)
}

pool_mic <- function(est.mat, wiv.mat, padj.method="holm"){
  n.imp <- ncol(est.mat); Eest <- rowMeans(est.mat)
  VW <- apply(wiv.mat, 1, function(x) mean(x))
  VB <- apply(est.mat, 1, function(x) sum((x-mean(x))^2)/(n.imp-1))
  VT <- VW + (1+1/n.imp)*VB
  pool.Wald.value <- Eest/sqrt(VT)
  p.val <- pnorm(abs(pool.Wald.value), lower.tail=FALSE)
  p.adj <- p.adjust(p.val, method=padj.method)
  rslt <- data.frame(est=Eest, se=sqrt(VT), W=pool.Wald.value, p=p.val, adjp=p.adj)
  rownames(rslt) <- rownames(est.mat)
  return(rslt)
}

print.mic <- function(x, ...){
  cat(length(x), " Imputed Profiles")
}

