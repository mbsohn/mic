mic <- function(M, X=NULL, wt=NULL, nimp=100, pdenom=0.5){
  if(is.null(wt)){
    lib.size <- rowSums(M)
    wt <- lib.size/sum(lib.size)
  }
  if(!(all(X[,1] == 1))) stop("The first column of X must be 1s!!!")
  n.cores <- detectCores(); optimem.cl <- makeCluster(n.cores-1)
  registerDoParallel(cl=optimem.cl)
  n.taxa <- ncol(M)
  miM <- list()
  rbind.M <- foreach(i=1:nimp, .combine="rbind") %dopar% {
    tmp.M <- M
    tx.indx <- 1:n.taxa
    for(j in 1:n.taxa){
      m1j <- NULL
      ran.den.indx <- sample(tx.indx[-j], round(n.taxa*pdenom))
      xmj <- rowSums(tmp.M[,ran.den.indx])
      rj <- tmp.M[,j]/xmj
      rj.zero.indx <- which(rj==0 | is.infinite(rj))
      if(length(rj.zero.indx)>0){
        if(is.null(X)){
          m1j <- weighted.mean(log(rj[-rj.zero.indx]), w=wt[-rj.zero.indx])
        } else{
          Ws <- diag(wt[-rj.zero.indx])
          Xs <- X[-rj.zero.indx,,drop=FALSE]
          Ss.inv <- solve(crossprod(Xs, Ws%*%Xs))
          beta1 <- Ss.inv%*%crossprod(Xs, Ws)%*%log(rj[-rj.zero.indx])
          m1j <- X[rj.zero.indx,]%*%beta1
        }
        xij <- exp(m1j+log(xmj[rj.zero.indx]))
        tmp.M[rj.zero.indx,j] <- ifelse(xij==0, NA, xij)
      }
    }
    tmp.M <- proportions(tmp.M, margin=1)
  }
  stopCluster(cl=optimem.cl)
  miM <- lapply(split(rbind.M, rep(c(1:nimp), each=(nrow(rbind.M)/nimp))), matrix, nrow(rbind.M)/nimp)
  na.indx <- which(sapply(miM, function(x) any(is.na(x))) == TRUE)
  l.na.indx <- length(na.indx)
  if(l.na.indx == nimp){
    stop("Too many zeros!!! Try to increase the value of parameter `pdenom'.")
  }
  if(l.na.indx > 0){
    miM <- miM[-na.indx]
  }
  names(miM) <- paste0("ImpProfile", 1:length(miM))
  miM <- lapply(miM, function(x) {colnames(x) <- colnames(M); x})
  class(miM) <- "mic"
  return(miM)
}

pool_mic <- function(est.mat, se.mat, padj.method="holm"){
  n.imp <- ncol(est.mat); Eest <- rowMeans(est.mat)
  VW <- apply(se.mat, 1, function(x) mean(x^2))
  VB <- apply(est.mat, 1, function(x) sum((x-mean(x))^2)/(n.imp-1))
  VT <- VW + (1+1/n.imp)*VB
  pool.Wald.value <- Eest/sqrt(VT)
  p.val <- pnorm(abs(pool.Wald.value), lower.tail=FALSE)
  p.adj <- p.adjust(p.val, method=padj.method)
  rslt <- data.frame(W=pool.Wald.value, p=p.val, adjp=p.adj)
  rownames(rslt) <- rownames(est.mat)
  return(rslt)
}

print.mic <- function(x, ...){
  cat(length(x), "Imputed Profiles")
}
