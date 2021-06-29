# posterior mean of parameters -----------

post.mean <- function(df_of_draws, dat, nBasis){
  nSub <- length(unique(dat$ID))
  b0 <- mean(df_of_draws[,'beta0'])
  bMatrix <- matrix(
    apply(df_of_draws[,sprintf('beta[%d,%d]',rep(seq(nSub), each=nBasis),seq(nBasis))],
          2, mean),
    nrow=nSub,byrow = T)
  G <-  matrix(
    apply(df_of_draws[,sprintf('VarBeta[%d,%d]',rep(seq(nBasis), each=nBasis),seq(nBasis))],
          2, mean),
    nrow=nBasis,byrow = T)
  Gamma <-
    mean(df_of_draws[,'sigmaepsilon'])^2

  list(b0=b0,bMatrix=bMatrix,G=G,Gamma=Gamma)
}

# fitted values based on posterior mean ----------------------

pc.fit <- function(ls_par, dat, ls_idxA,seed){
  mat_fitted <- data.frame(
    matrix(NA,nrow=nrow(dat),ncol=length(ls_idxA),dimnames=list(c(),sprintf('y_hatA%d',seq(length(ls_idxA))))))
  id <-  as.numeric(dat$ID)

  b0 <- ls_par$b0
  bMatrix <- ls_par$bMatrix
  G <- ls_par$G
  Gamma <- ls_par$Gamma


  foreach(t = seq(1,length(ls_idxA))) %dopar% {
    idxA <- ls_idxA[[t]]
    idxB <- setdiff(ls_idxA[[1]],idxA)
    set.seed(seed)

    cat(sprintf('Set %d: calculating fitted curve at MAP\n',t))

    bAMatrix <- bMatrix[,idxA]

    mat_fitted_t <- my.fit(idxA, idxB, dat, G,bAMatrix,b0, id)
    list(mat_fitted_t)
  }->temp
  list(
    mat_fitted=matrix(unlist(lapply(temp,'[[',1)),ncol=length(ls_idxA),byrow=F,dimnames=list(c(),sprintf('y_hatA%d',seq(length(ls_idxA)))))
  )
}

# cluster number by KL method --------------------------------

pc.KL <- function(ls_par, dat, ls_idxA, nIter, thKL, regQ, seed){
  id <-  as.numeric(dat$ID)
  nSub <- length(unique(id))

  bMatrix <- ls_par$bMatrix
  G <- ls_par$G
  Gamma <- ls_par$Gamma


  foreach(t = seq(1,length(ls_idxA))) %dopar% {
    idxA <- ls_idxA[[t]]
    idxB <- setdiff(ls_idxA[[1]],idxA)
    set.seed(seed)

    cat(sprintf('Set %d: optimizing cluster number by k means\n',t))

    bAMatrix <- bMatrix[,idxA]
    Q <- Q.func(idxA, idxB, dat,Gamma, G,id)
    Kcluster.obj <- Kcluster.pick(nIter, bAMatrix,Q, regQ, do.plot=F, seed=seed)
    threshold <- Kcluster.obj$KL[1]*thKL # find ncluster such that KL<threshold
    nClusters_t <- which(Kcluster.obj$KL<threshold)[1]
    KLs_t <- Kcluster.obj$KL
    cluster_t <- Kcluster.obj$clusterMatrix

    list(nClusters_t, KLs_t, cluster_t)
  }->temp

  list(
    nClusters=unlist(lapply(temp,'[[',1)),
    KLs=lapply(temp,'[[',2),
    cluster0=lapply(temp,'[[',3)
  )
}

# draw from posterior with chosen nCluster --------------------

pc.pair <- function(df_of_draws, dat, nBasis, ls_idxA, nIter, nDraw,nClusters, regQ, seed){
  id <-  as.numeric(dat$ID)
  nSub <- length(unique(id))
  if(length(nClusters)==1) nClusters <- rep(nClusters,length(ls_idxA))

  foreach(t = seq(1,length(ls_idxA))) %dopar% {
    idxA <- ls_idxA[[t]]
    idxB <- setdiff(ls_idxA[[1]],idxA)
    nCluster <- nClusters[t]

    set.seed(seed)

    cat(sprintf('Set %d: drawing samples from posterior\n',t))
    post.obj <-
      cluster.posterior(df_of_draws,nDraw, dat, nBasis,idxA, idxB, id, nCluster,nIter, regQ, seed)

    # # posterior of KL
    # hist(post.obj$KL,main=sprintf('Posterior of KL with %d clusters',nCluster))
    #
    # # posterior of centroid
    # plot(post.obj$dKArray[,1,],post.obj$dKArray[,2,],pch='.', main=sprintf('Posterior of centroids with %d clusters',nCluster))

    # posterior probability of 2 subjects in 1 cluster
    tb <- matrix(NA, nrow= nSub, ncol=nSub)
    for(i in seq(2,nSub)){
      for(j in seq(i-1)){
        tb[i,j] <- mean(post.obj$clusterMatrix[i,]==post.obj$clusterMatrix[j,])
      }
    }
    ls_prob_t <- tb
    list(ls_prob_t,
         post.obj$clusterMatrix, post.obj$clusterVec)
  }->temp

  list(
    ls_prob=lapply(temp,'[[',1),
    ls_clust=lapply(temp,'[[',2),
    ls_clust0=lapply(temp,'[[',3)
  )
}

# cluster number by bootstrap on the fitted curve ---------------------------------

kmeans.predict <- function(x, centers) {

  # compute squared euclidean distance from each sample to each cluster center
  tmp <- sapply(seq_len(nrow(x)),
                function(i) apply(centers, 1,
                                  function(v) sum((x[i, ]-v)^2)))
  if(!is.matrix(tmp)) tmp <- t(tmp)
  max.col(-t(tmp))  # find index of min distance
}

nCluster.boot <- function(mat_fitted, dat, nB, nCl, seed){
  set.seed(seed)
  nSet <- ncol(mat_fitted)
  sBtb <- matrix(NA,nrow=nCl-1, ncol=nSet)
  nClusters = cutoff <- numeric(nSet)

  for(j in seq(nSet)){
    xx <- dat %>% # formatting observations for distance
      select(ID,t) %>%
      bind_cols(data.frame(y=mat_fitted[,j])) %>%
      pivot_wider(names_from = t, values_from=y) %>%
      select(-ID) %>%
      as.matrix
    nSub <- nrow(xx)
    sB <- sapply(seq(2,nCl), function(cl){
      print(cl)
      mean(sapply(seq(nB),function(b){
        idx1 <- sample(nSub,nSub,replace=T)
        idx2 <- sample(nSub,nSub,replace=T)
        out1 <- kmeans(xx[idx1,], cl)
        out2 <- kmeans(xx[idx2,], cl)
        fitted1 <- fitted(out1, method='classes')
        fitted12 <- kmeans.predict(xx[idx2,], out1$centers)
        fitted2 <- fitted(out2, method='classes')
        fitted21 <- kmeans.predict(xx[idx1,], out2$centers)

        # empirical clustering distance
        idx_pair <- matrix(c(rep(1:nSub,nSub),rep(1:nSub,each=nSub)),nrow=2, ncol=nSub^2, byrow = T)
        d <- mean(apply(idx_pair,2,function(x){
          (fitted1[x[1]]==fitted12[x[2]] & fitted21[x[1]]!=fitted2[x[2]] )+
            ( fitted1[x[1]]!=fitted12[x[2]] & fitted21[x[1]]==fitted2[x[2]])
        }))
      }))
    })
    sBtb[,j] <- sB
    cutoff[j] <- max(sB)/2
    idx <- setdiff(which(sB<=cutoff[j]), 1)
    nClusters[j] <- idx+1
  }
  ls_sB <- list(sBtb=sBtb, cutoff=cutoff, nClusters=nClusters)
  ls_sB
}

