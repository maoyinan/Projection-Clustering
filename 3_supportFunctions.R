Q.func <- function(idxA, idxB, dat,Gamma, G, id){
  nSub= length(unique(id))

  GA <- G[idxA,idxA]
  GA1 <- solve(GA)
  if(length(idxB)>0){
    GB <- G[idxB,idxB]
    GAB <- G[idxA,idxB]
    ZB <- as.matrix(dat[,paste0('Z',idxB)])
  }
  ZA <- as.matrix(dat[,paste0('Z',idxA)])

  Q <- array(NA, dim=c(nSub, length(idxA), length(idxA)))
  for(i in seq(nSub)){
    idx <- which(id==i)
    Gammai <- Gamma*diag(1,nrow=length(idx),ncol=length(idx))
    ZiA <- ZA[idx,]
    if(length(idxB)>0){
      ZiB <- ZB[idx,]
      term1 <- ZiA+ZiB%*%t(GAB)%*%GA1
      term2 <- ZiB%*%(GB-t(GAB)%*%solve(GA)%*%GAB)%*%t(ZiB) + Gammai
      Q[i,,] <- t(term1) %*% solve(term2) %*% term1

    }
    else{
      Q[i,,] <- Gamma*diag(1,nrow=length(idxA),ncol=length(idxA))

    }

  }
  Q
}

# k means algorithm -------------------------------------------------------

dK.func0 <- function(nCluster, clusterVec, Q, bAMatrix,regQ){
  dK <- matrix(NA, nrow=nCluster, ncol=ncol(bAMatrix))
  for(j in seq(nCluster)){
    ind <- which(clusterVec==j)
    if(length(ind)>1){
      term1 <- apply(Q[ind,,], 2:3, sum)
      term2 <- apply(sapply(ind, function(i)Q[i,,]%*%bAMatrix[i,]),1,sum)
    } else {
      term1 <- Q[ind,,]
      term2 <- Q[ind,,]%*%bAMatrix[ind,]
    }
    term11 <-
      tryCatch(expr = solve(term1),
               error = function(c) {
                 # matlib::Ginv(term1, tol=regQ)
                 solve(term1+diag(regQ,nrow(term1),ncol(term1)))
               })

    dK[j,] <- term11%*%term2
  }
  dK
}
dK.func1 <- function(nCluster, clusterVec, ls_Q, ls_bAMatrix, regQ){
  # used for summing over posterior samples
  sum_Q <- apply(ls_Q, 1:3, sum)
  sum_Qb <- array(0, dim=dim(ls_Q)[1:2])
  for(k in seq(dim(ls_Q)[1])){
    for(i in seq(nDraw))
      sum_Qb[k,] <- sum_Qb[k,] + c(ls_Q[k,,,i] %*% ls_bAMatrix[k,,i])
  }

  dK <- matrix(NA, nrow=nCluster, ncol=dim(ls_Q)[2])
  for(j in seq(nCluster)){
    ind <- which(clusterVec==j)
    if(length(ind)>1){
      term1 <- apply(sum_Q[ind,,], 2:3, sum)
      term2 <- apply(sapply(ind, function(i)sum_Qb[i,]),1,sum)
    } else {
      term1 <- sum_Q[ind,,]
      term2 <- sum_Qb[ind,]
    }
    term11 <-
      tryCatch(expr = solve(term1),
               error = function(c) {
                 # matlib::Ginv(term1, tol=regQ)
                 solve(term1+diag(regQ,nrow(term1),ncol(term1)))
               })

    dK[j,] <- term11%*%term2
  }
  dK
}

KL.func <- function(diA, biA, Qi){
  1/2*t(diA-biA)%*%Qi%*%(diA-biA)
}

C.func0 <- function(nCluster, dK, Q, bAMatrix){
  clusterVec1 <- numeric(nrow(bAMatrix))
  termMatrix <- matrix(NA,nrow=nrow(bAMatrix),ncol=nCluster)
  for(i in seq_along(clusterVec1)){
    biA <- bAMatrix[i,]
    term <- numeric(nCluster)
    for(j in seq(nCluster)){
      term[j] <- KL.func(dK[j,],biA, Q[i,,])
      termMatrix[i,j] <- term[j]
    }
    clusterVec1[i] <- which(term==min(term))
  }
  clusterVec1
}
C.func1 <- function(nCluster, dK, ls_Q, ls_bAMatrix){
  nSub <- dim(ls_Q)[1]
  clusterVec1 <- numeric(nSub)
  termMatrix <- matrix(NA,nrow=nSub,ncol=nCluster)
  for(i in seq(nSub)){
    term <- numeric(nCluster)
    for(j in seq(nCluster)){
      for(k in seq(nDraw))
        term[j] <-  term[j] + KL.func(dK[j,],ls_bAMatrix[i,,k], ls_Q[i,,,k])

      # term[j] <- KL.func(dK[j,],biA, Q[i,,])
      termMatrix[i,j] <- term[j]
    }
    clusterVec1[i] <- which(term==min(term))
  }
  clusterVec1
}
# cluster.initialize <- function(dK0,clusterVec0,nCluster,Q,bAMatrix, seed){
#   set.seed(seed)
#   dK <- matrix(NA, nrow=nCluster, ncol=ncol(bAMatrix))
#   iCluster <- seq(nCluster-1)[which(table(clusterVec0)>1)[1]]
#   if(nCluster>2) dK[seq(nCluster-2),] <- dK0[setdiff(seq(nCluster-1), iCluster),]
#   # split the first cluster of previous iteration (who has more than 1 points)
#   idx <- which(clusterVec0==iCluster)
#   dK[seq(nCluster-1,nCluster),] <- bAMatrix[sort(sample(idx, 2, replace = F)),]
#   clusterVec0 <- C.func(nCluster,dK, Q, bAMatrix)
#   clusterVec0
# }
cluster.initialize <- function(clusterVec0,nCluster,Q,bAMatrix,regQ, seed){
  set.seed(seed)
  clusterVec <- clusterVec0
  temp <- which(table(clusterVec0)>1)
  if(length(temp)==1)  {
    iCluster <- temp
    # split the first cluster of previous iteration (who has more than 1 points)
    idx <- which(clusterVec0==iCluster)
    split <- cluster.initialize0(2,Q[idx,,], bAMatrix[idx,],seed)
    clusterVec[idx[split==2]] <- nCluster
  }
  else {
    clusterMat <- matrix(NA, nrow=length(clusterVec),ncol=length(temp))
    kl <- numeric(length(temp))
    for(j in seq_along(kl)){
      iCluster <- temp[j]

      # split the first cluster of previous iteration (who has more than 1 points)
      idx <- which(clusterVec0==iCluster)
      split <- cluster.initialize0(2,Q[idx,,], bAMatrix[idx,],seed)
      clusterVec[idx[split==2]] <- nCluster
      clusterMat[,j] <- clusterVec
      dK <- dK.func0(nCluster, clusterVec, Q, bAMatrix,regQ)
      kl[j] <- KLK.func(dK,clusterVec,bAMatrix,Q)
    }
    ii <- which(kl==min(kl))
    if(length(ii)==1) clusterVec <- clusterMat[,ii]
    else clusterVec <- clusterMat[,sample(ii,1)]

  }
  clusterVec
}
cluster.initialize0 <- function(nCluster,Q,bAMatrix, seed){
  set.seed(seed)
  dK <- bAMatrix[sort(sample(seq(nrow(bAMatrix)), nCluster, replace = F)),]
  clusterVec0 <- C.func0(nCluster,dK, Q, bAMatrix)
  clusterVec0
}

cluster.iter <- function(clusterVec0,nIter,nCluster,Q,bAMatrix,regQ, seed, do.plot){
  ls_Q <- Q
  ls_bAMatrix <- bAMatrix
  if(length(dim(Q))==4){
    Q <- ls_Q[,,,1]
    bAMatrix <- ls_bAMatrix[,,1]
    dK.func <- dK.func1
    C.func <- C.func1
  }
  else{
    dK.func <- dK.func0
    C.func <- C.func0
  }

  clusterMatrix <- matrix(NA, nrow= nrow(bAMatrix), ncol= nIter)
  dKArray <- array(NA, dim=c(nCluster, ncol(bAMatrix), nIter))
  if(is.null(clusterVec0))
    clusterVec <- cluster.initialize0(nCluster,Q,bAMatrix, seed)
  else clusterVec <- cluster.initialize(clusterVec0,nCluster, Q,bAMatrix,regQ, seed)

  dK <- dK.func(nCluster, clusterVec, ls_Q, ls_bAMatrix,regQ)
  clusterMatrix[,1] <- clusterVec
  dKArray[,,1] <- dK

  if(do.plot){
    plot(bAMatrix[,1],bAMatrix[,2], col=clusterVec, main=sprintf('Iter 1 with %d clusters',nCluster))
    text(dK[,1],dK[,2],seq(nCluster), col=seq(nCluster), cex=2)
  }

  for(i in seq(2,nIter)){
    #update
    clusterVec1 <- C.func(nCluster,dK, ls_Q, ls_bAMatrix)
    # if(sum(clusterVec1-clusterMatrix[,i-1])==0) break
    if(sum(clusterVec1-clusterMatrix[,i-1])==0 | length(unique(clusterVec1))<nCluster) break

    clusterVec <- clusterVec1
    dK <- dK.func(nCluster, clusterVec, ls_Q, ls_bAMatrix, regQ)
    clusterMatrix[,i] <- clusterVec
    dKArray[,,i] <- dK

    if(do.plot){
      plot(bAMatrix[,1],bAMatrix[,2], col=clusterVec, main=sprintf('Iter %d with %d clusters',i,nCluster))
      text(dK[,1],dK[,2],seq(nCluster), col=seq(nCluster), cex=2)
    }
  }
  clusterMatrix[,nIter] <- clusterVec
  dKArray[,,nIter] <- dK
  list(clusterMatrix,dKArray)
}

KLK.func <- function(dK,clusterVec,bAMatrix,Q){
  nSub <- nrow(bAMatrix)
  KLKvec <- numeric(nSub)
  for(i in seq(nSub)){
    KLKvec[i] <- KL.func(dK[clusterVec[i],], bAMatrix[i,], Q[i,,])
  }
  KLK <- sum(KLKvec)
  KLK
}

Kcluster.pick <- function(nIter, bAMatrix,Q,regQ, do.plot, seed){
  nSub <- nrow(bAMatrix)
  KL <- numeric(nSub)
  nCluster <- 1
  clusterVec <- rep(1,nSub)
  dK <- dK.func0(nCluster, clusterVec, Q, bAMatrix,regQ)
  KL[1] <- KLK.func(dK,clusterVec,bAMatrix,Q)

  clusterMatrix <- matrix(NA, nrow= nSub, ncol= nSub)
  clusterMatrix[,1] <- clusterVec

  for(nCluster in seq(2,nSub-1)){
    # for(nCluster in seq(2,41)){
    temp <- cluster.iter(clusterVec,nIter,nCluster,Q,bAMatrix,regQ, seed,do.plot)
    dK <- temp[[2]][,,nIter]
    clusterVec <- temp[[1]][,nIter]
    KL[nCluster] <- KLK.func(dK,clusterVec,bAMatrix,Q)
    clusterMatrix[,nCluster] <- clusterVec

  }
  list(KL=KL, clusterMatrix=clusterMatrix)
}

# cluster poeterior distribution ------------------------------------------
# require fit object from stan (df_of_draws)

cluster.posterior <- function(df_of_draws,nDraw, dat, nBasis,idxA, idxB, id, nCluster,nIter, regQ, seed){

  nSub= length(unique(id))
  idxDraw <- sample(nrow(df_of_draws),nDraw)

  foreach(idx = idxDraw) %dopar% {
    bAMatrix <- matrix(unlist(
      df_of_draws[idx,sprintf('beta[%d,%d]',rep(seq(nSub), each=nBasis),seq(nBasis))]),
      nrow=nSub,byrow = T)[,idxA]
    Gamma <- unlist(df_of_draws[idx,'sigmaepsilon'])^2
    G <-  matrix(unlist(
      df_of_draws[idx,sprintf('VarBeta[%d,%d]',rep(seq(nBasis), each=nBasis),seq(nBasis))]),
      nrow=nBasis,byrow = T)
    Q <- Q.func(idxA, idxB, dat,Gamma, G, id)

    temp <- cluster.iter(NULL,nIter,nCluster,Q,bAMatrix,regQ, seed,do.plot=F)
    dK <- temp[[2]][,,nIter]
    clusterVec <- temp[[1]][,nIter]

    list(KLK.func(dK,clusterVec,bAMatrix,Q),clusterVec, dK, Q, bAMatrix)
  } -> temp
  ls_Q <- array(unlist(lapply(temp,'[[',4)), dim=c(nSub,length(idxA),length(idxA), nDraw))
  ls_bAMatrix <- array(unlist(lapply(temp,'[[',5)), dim=c(nSub, length(idxA), nDraw))
  temp1 <- cluster.iter(NULL,nIter,nCluster,ls_Q,ls_bAMatrix,regQ, seed,do.plot=F)
  clusterVec <- temp1[[1]][,nIter]

  list( # must check if element sequence is right
    KL=unlist(lapply(temp,'[[',1)),
    clusterMatrix=matrix(unlist(lapply(temp,'[[',2)),ncol=nDraw,byrow=F),
    dKArray=array(unlist(lapply(temp,'[[',3)), dim=c(nCluster, length(idxA), nDraw)),
    clusterVec=clusterVec
  )
}
cluster.posterior0 <- function(df_of_draws,nDraw, dat, nBasis,idxA, idxB, id, nCluster,nIter, regQ, seed){

  nSub= length(unique(id))
  idxDraw <- sample(nrow(df_of_draws),nDraw)
  KL <- numeric(nDraw)
  clusterMatrix <- matrix(NA, nrow= nSub, ncol= nDraw)
  dKArray <- array(NA, dim=c(nCluster, length(idxA), nDraw))

  for(iter in seq(nDraw)){
    if(iter%%(nDraw/10)==0) cat(sprintf('%d out of %d\n',iter,nDraw))
    idx <- idxDraw[iter]
    bAMatrix <- matrix(unlist(
      df_of_draws[idx,sprintf('beta[%d,%d]',rep(seq(nSub), each=nBasis),seq(nBasis))]),
      nrow=nSub,byrow = T)[,idxA]
    Gamma <- unlist(df_of_draws[idx,'sigmaepsilon'])^2
    G <-  matrix(unlist(
      df_of_draws[idx,sprintf('VarBeta[%d,%d]',rep(seq(nBasis), each=nBasis),seq(nBasis))]),
      nrow=nBasis,byrow = T)
    Q <- Q.func(idxA, idxB, dat,Gamma, G, id)

    temp <- cluster.iter(NULL,nIter,nCluster,Q,bAMatrix,regQ, seed,do.plot=F)
    dK <- temp[[2]][,,nIter]
    clusterVec <- temp[[1]][,nIter]

    KL[iter] <- KLK.func(dK,clusterVec,bAMatrix,Q)
    clusterMatrix[,iter] <- clusterVec
    dKArray[,,iter] <- dK
  }

  list(KL=KL, clusterMatrix=clusterMatrix, dKArray=dKArray)
}


# compare different projection dimension ---------------------------

my.fit <- function(idxA, idxB, dat,G,bAMatrix,b0,  id){
  nSub= length(unique(id))
  GA <- G[idxA,idxA]
  GA1 <- solve(GA)
  if(length(idxB)>0){
    GAB <- G[idxA,idxB]
    ZB <- as.matrix(dat[,paste0('Z',idxB)])
    term <- t(GAB)%*%GA1
  }
  ZA <- as.matrix(dat[,paste0('Z',idxA)])

  y <- numeric(length(id))
  for(i in seq(length(y))){
    if(length(idxB)>0){
      y[i] <- b0+as.numeric((ZA[i,] + ZB[i,]%*%term) %*% bAMatrix[id[i],])
    }
    else
      y[i] <- b0+as.numeric((ZA[i,]) %*% bAMatrix[id[i],])

  }
  y
}

estimate_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}
