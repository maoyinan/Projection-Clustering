
greedyXWVfixed <- function(y,X,W,V,U,vni,epm,M,IGa,IGb,sigdelp,sigbetap,nClusters) {
  
  bef <- proc.time() 
  
  # Fit 1 comp model #
  k <- 1
  n <- length(vni)
  nsum <- sum(vni)
  p <- dim(X)[2]
  s1 <- dim(W)[2]
  s2 <- dim(V)[2]
  d <- dim(U)[2]
  g <- dim(epm)[2]
  
  # Priors #
  sigbeta <- array(0,dim=c(k,p,p))
  for (j in 1:k) {sigbeta[j,,] <- sigbetap*diag(p)}
  alpa <- rep(IGa,k)
  lama <- rep(IGb,k)
  alpb <- rep(IGa,k)
  lamb <- rep(IGb,k)
  alp <- matrix(IGa,g,k)
  lam <- matrix(IGb,g,k)
  sigdel <- sigdelp*diag((k-1)*d)
  qp <- matrix(1,nrow=n,ncol=k)
  fitpre <- VLLMf(y,X,U,sigbeta,alpa,lama,alpb,lamb,alp,
                  lam,sigdel,vni,epm,qp,1.0e-5)
  # View(fitpre$qp)
  
  dsp <- NULL
  dif <- 10
  tot <- seq(1,n,1) 
  
  # save("image_greedy.RData")
  
  while (dif > 1 | k<nClusters) {
    
    # k <- length(fitpre$alpaq)
    postp <- rep(0,n)
    for(i in 1:n){
      for(j in 1:k){
        if (fitpre$qp[i,j] == max(fitpre$qp[i,])) (postp[i]<-j)}} 
    gp <- array(0,dim=c(k,n,M))
    for (j in 1:k) {
      for (m in 1:M) {
        current <- tot[postp==j]
        len <- length(current)
        sam <- sample(c(1,2),len,replace=TRUE,prob=c(0.5,0.5))
        gp[j,which(postp==j),m] <- sam }}
    
    
    # (K <- dim(gp)[1])
    fitore <- fitpre
    fits <- list(NULL)
    (k <- k+1)
    
    # Priors #
    sigbeta <- array(0,dim=c(k,p,p))
    for (j in 1:k) {sigbeta[j,,] <- sigbetap*diag(p)}
    alpa <- rep(IGa,k)
    lama <- rep(IGb,k)
    alpb <- rep(IGa,k)
    lamb <- rep(IGb,k)
    alp <- matrix(IGa,g,k)
    lam <- matrix(IGb,g,k)
    sigdel <- sigdelp*diag((k-1)*d)
    
    if (length(dsp)>=1) {sp <- seq(1,k-1,1)[-dsp]} else
      {sp <- seq(1,k-1,1)} 
    lbds <- cbind(sp,rep(0,length(sp)))
    num <- 0
    for (cp in sp) {
      num <- num+1
      fitold <- -10000000
      nredz <- 0
      for (m in 1:M) {
        fitcan <- VLLMfp(y,X,U,sigbeta,alpa,lama,alpb,lamb,alp,
                         lam,sigdel,vni,epm,1,fitpre,gp,cp,m) 
        cat(cp,m,fitcan$lb,"\n")
        
        # if (fitcan$lb > fitold) {
          fits[[cp]] <- fitcan
          lbds[num,2] <- fitcan$lb
          fitold <- fitcan$lb}# }
      
      if (any( colSums(fits[[cp]]$qp) <= 1)) {
        dsp <- c(dsp,cp)
        lbds[num,2] <-  -10000000
        }
      cat(cp,dsp,"\n")
    }
    lbdsord <- lbds[,2]
    fitpre <- fits[[lbds[which(lbdsord==max(lbdsord))[1],1]]]
    
    #Priors#
    sigbeta <- array(0,dim=c(k,p,p))
    for (j in 1:k) {sigbeta[j,,] <- sigbetap*diag(p)}
    alpa <- rep(IGa,k)
    lama <- rep(IGb,k)
    alpb <- rep(IGa,k)
    lamb <- rep(IGb,k)
    alp <- matrix(IGa,g,k)
    lam <- matrix(IGb,g,k)
    sigdel <- sigdelp*diag((k-1)*d)
    fitfinal <- VLLMf(y,X,U,sigbeta,alpa,lama,alpb,lamb,alp,
                      lam,sigdel,vni,epm,fitpre$qp,1.0e-5,fitpre)
    
    dsptemp <- dsp
    dsp <- NULL
    for (dd in dsptemp){
      dg <- which(fitpre$lamaq==fitore$lamaq[dd])
      if (max(abs(fitpre$qp[,dg]-fitfinal$qp[,dg])) <= 0.5){
        dsp <- c(dsp,dg)} }
    
    if ( any(colSums(fitfinal$qp)<=1) ){dsp <- c(dsp,which(colSums(fitfinal$qp)<=1))}
    
    fitpre <- fitfinal
    (lbpre <- fitore$lbadj)
    (lbnext <- fitfinal$lbadj)
    (dif <- lbnext-lbpre)
    cat(dim(fitpre$mubetaq)[2],fitpre$lbadj,dsp,"\n") 
  }
  
  aft <- proc.time()
  dur <- aft-bef
  list(dur=dur,fitpre=fitpre)}



