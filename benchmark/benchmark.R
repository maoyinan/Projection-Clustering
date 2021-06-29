library(tidyverse)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("BHC")
library(BHC)
library(TSclust)
library(kml) # Genolini 2010
library(mclust) # no time dependency
library(fossil) # rand index
library(kableExtra)
select <- dplyr::select
paths <- paste0('../',
                c('eg2_crop',
                  'eg3_yeast',
                  'eg4_sleep',
                  'eg5_active',
                  'eg5m_active_x1x2',
                  'eg5g_active_gap'
                ))

# mean of KL
pc.rand4 <- function(true_cluster, pcOutput, df_groups,nDraw, ls_idxA, rand.fcn){
  arr_cluster <- array(unlist(pcOutput$ls_clust0), 
                       dim=c(nrow(df_groups), length(ls_idxA)))
  apply(arr_cluster, 2, function(y)  {
    ret <- rand.fcn(true_cluster, y)
    ret
  })
}

# some classes distinguishing better than others
# adjusted rand index
# weighted index prefer one pair of grouping 


# benchmark clustering methods using rand index  ------------------------------------------------

benchmark.rand <- function(ls_benchmark, true_cluster, rand.fcn){
  ret <- rep(NA,length(ls_benchmark))
  for(i in seq_along(ls_benchmark)){
    if(!is.null(ls_benchmark[[i]])) ret[i] <- rand.fcn(true_cluster,ls_benchmark[[i]])
  }
  ret
}

sd <- 12345
nData <- 4

for(i in seq(nData)){
  cat(i,'\n')
  load(file.path(paths[i],'dat'))  
  dat1 <- dat %>% 
    select(ID, t, Record) %>% 
    pivot_wider(names_from = ID, values_from= Record) %>% 
    select(-t) %>% 
    as.data.frame
  dat %>% select(Group, ID) %>% unique %>% 
    select(Group) %>% unlist %>% as.factor %>% as.numeric -> true_cluster
  
  nCluster <- length(unique(true_cluster))
  
  set.seed(sd)
  IP.dis <- diss(dat1, "INT.PER")#"DWT"
  IP.hclus <- cutree(hclust(IP.dis), k = nCluster)
  
  diffs <- rep(1, ncol(dat1))
  logs <- rep(F, ncol(dat1))
  hc.hclus <- NA
  hc.hclus <- tryCatch({
    dpred <- diss(dat1, "PRED", h = nCluster, B = 1200, logarithms = logs, differences = diffs)
    cutree( hclust(dpred$dist), k = nCluster)
  },
  error=function(cond){
    message('No values returned')
    message(cond)
  })
  
  dat2 <- t(dat1)
  if(length(unique(c(dat2)))>10){
    hc2 <- bhc(dat2, rownames(dat2), 0, seq(ncol(dat2)), "time-course",
               numReps=1, noiseMode=0, numThreads=1, verbose=TRUE) # account for correlation in time series (must be continuous): squared exponential covariance
  }
  else hc2 <- bhc(dat2, rownames(dat2), 0, seq(ncol(dat2)), "multinomial", # used for discrete with 2+ categories data
                  numReps=1, noiseMode=0, numThreads=1, verbose=TRUE)
  bhc.hc <- cutree(as.hclust(hc2), k=nCluster)
  
  cld1 <- clusterLongData(dat2)
  kml(cld1,nCluster,10)
  kml.c <- as.numeric(getClusters(cld1,nCluster,1))
  
  BIC <- mclustBIC(dat2)
  mod1 <- Mclust(dat2, x = BIC, G = nCluster) # chosen cluster number always the minimum of G range
  mclust.c <- mod1$classification
  
  ls_benchmark <- list(IP.hclus=IP.hclus, hc.hclus=hc.hclus, bhc.hc=bhc.hc, kml.c=kml.c, mclust.c=mclust.c)
  save(ls_benchmark, file=sprintf('ls_benchmark%d',i))
}


tb_rand <- array(NA, dim = c(nData, 5,2), dimnames = list(
  paste0('Eg',seq(nData)), 
  c('HC_dist', 'HC_pred', 'BHC', 'KML', 'Mclust'),
  c('rand','adj_rand')))

for(i in seq(nData)){
  load(sprintf('ls_benchmark%d',i))
  load(file.path(paths[i],'dat'))  
  dat %>% select(Group, ID) %>% unique %>% 
    select(Group) %>% unlist %>% as.factor %>% as.numeric -> true_cluster
  for(j in 1:2){
    rand.fcn <- list(rand.index,adj.rand.index)[[j]]
    z1 <- benchmark.rand(ls_benchmark, true_cluster, rand.fcn)
    tb_rand[i, ,j] <- z1
  }
}
tb_rand_benchmark <- tb_rand

# projection clustering rand index --------------------------------------

k=4 
nData <- 6
{
  pc.rand <- get(paste0('pc.rand',k))
  
  tb_rand <- array(NA, c(nData, 4, 2), dimnames = list(
    paste0('Eg', seq(nData)+1), 
    paste0('PC', seq(4)),
    c('rand','adj_rand')))
  z0 <- numeric(4)
  
  for(i in seq(nData)){
    load(file.path(paths[i],sprintf('pcOutputEg%dFixed',i+1)))  
    load(file.path(paths[i],'dat'))  
    
    dat %>% select(Group, ID) %>% unique %>% 
      select(Group) %>% unlist %>% as.factor %>% as.numeric -> true_cluster
    
    nCluster <- length(unique(true_cluster))
    for(j in 1:2){
      rand.fcn <- list(rand.index,adj.rand.index)[[j]]
      
      z0 <- pc.rand(true_cluster, out_pc0, df_groups,nDraw, ls_idxA, rand.fcn)
      tb_rand[i, ,j] <- z0
    }
  }
  print(tb_rand)
}
tb_rand
save(tb_rand_benchmark,tb_rand,file='tb_rand')
# tb1 for comparing examples 1-3 using all methods
# tb2 for example 4-6 using projective clustering


# one more benchmark clustering method ---------------------
# variational computation

source('VariationalComputeClustering/VGA Alg 3 (full centering).R')
source('VariationalComputeClustering/greedyXWVfixed.R')


sd <- 12345
nData <- 6

tb3 <- matrix(NA,nrow=nData,ncol=2)
vaClusters <- list()

for(i in seq(nData)){
  load(file.path(paths[i],'dat'))  
  
  dat1 <- dat %>% 
    select(ID, t, Record) %>% 
    pivot_wider(names_from = ID, values_from= Record) %>% 
    select(-t) %>% 
    as.data.frame
  dat %>% select(Group, ID) %>% unique %>% 
    select(Group) %>% unlist %>% as.factor %>% as.numeric -> true_cluster
  nClusters <- length(unique(true_cluster))
  data <- t(dat1)
  
  # construct y #
  n <- dim(data)[1]
  np <- dim(data)[2]
  vni <- rep(np,n)
  vni <- apply(data,1,function(x)sum(!is.na(x)))
  nsum <- sum(vni) 
  y <- NULL
  for (k in 1:n) {
    vv <- data[k,1:np]
    y <- c(y,vv)}
  idx <- which(is.na(y))
  # Covariates #
  X <- NULL
  for (k in 1:n) {
    Xi <- diag(np)
    X <- rbind(X,Xi)}
  if(length(idx)) {
    y <- y[-idx]
    X <- X[-idx,]
  }
  W <- X
  V <- X
  time <- seq(1,n,1)
  U <- cbind(1,time,time^2,time^3)
  scle<-function(x) { -1+2*(x-min(x))/(max(x)-min(x)) }
  U[,2:4]<- apply(U[,2:4],2,scle)
  
  # specify data missing as 0, available as 1
  # epm <- matrix(1,nrow=n,ncol=np)
  epm <- matrix(0,nrow=n,ncol=np)
  for (j in 1:n) {
    for (t in 1:np) {
      epm[j,t] <- length(na.omit(data[j,t]))}}
  
  # Applying VGA using Algorithm 3 (full centering) #
  # source('VariationalComputeClustering/greedyXWVfixed.R')
  M <- 5
  sigdelp <- 1000
  sigbetap <- 10000
  IGa <- 2
  IGb <- round(2*(fitdistr(lm(y~X-1)$residuals,'t',m=0,df=4)$estimate)^2,2)
  set.seed(sd)
  # mixturemodel <- greedyXWV(y,X,W,V,U,vni,epm,M,IGa,IGb,sigdelp,sigbetap)
  # mixturemodel$dur
  
  mixturemodel <- greedyXWVfixed(y,X,W,V,U,vni,epm,M,IGa,IGb,sigdelp,sigbetap,nClusters)
  mixturemodel$dur
  fit <- mixturemodel$fitpre
  
  # plot(seq(1,np,1),data[1,],type='l',ylim=c(min(data),max(data)))
  # for(i in 1:n){points(seq(1,np,1),data[i,],type='l')}
  # xx <- fit$mubetaq
  # for(i in 1:ncol(xx)){points(seq(1,np,1),xx[,i],type='l',col='red')}
  # xx <- fit$murhoq
  # for(i in 1:nrow(xx)){points(seq(1,np,1),xx[i,],type='l',col='blue')}
  predClusters <- apply(fit$qp,1,function(x)which(x==max(x)))
  plot(seq(1,np,1),data[1,],type='l',ylim=c(min(data,na.rm = T),max(data,na.rm = T)))
  for(k in 1:n){points(seq(1,np,1),data[k,],type='l',col=c('red','blue','green','purple','orange')[predClusters[k]])}
  
  vaClusters[[i]] <- predClusters
  rd <- rand.index(true_cluster, predClusters)
  tb3[i,1] <- rd
  rd <- adj.rand.index(true_cluster, predClusters)
  tb3[i,2] <- rd
}

# output in latex table
load('tb_rand')
tb1 <- cbind(
  Example=paste0('Eg',rep(2:4,2)),
  rbind(tb_rand_benchmark[,,1],tb_rand_benchmark[,,2])[-c(4,8),],
  data.frame(VC=c(tb3[2:4,])),
  rbind(tb_rand[,,1],tb_rand[,,2])[c(1:3,7:9),]
)
tb2 <- cbind(
  Example=paste0('Eg5',rep(c('','M','G'),2)),
  data.frame(VC=c(tb3[4:6,])),
  rbind(tb_rand[,,1],tb_rand[,,2])[c(4:6,10:12),]
  )

kable(tb1, "latex",digits = 2,
             label = 'rand1',row.names = F,
             caption = 
               "Rand and adjusted Rand indices for different clustering methods for examples 2-4. The methods compared are described in the text.") %>% 
  kable_styling() %>% 
  pack_rows("Rand Index", 1, 3) %>% 
  pack_rows("Adjusted Rand Index", 4, 6) %>% 
  write( file='rand1.tex')

kable(tb2, "latex",digits = 2,
      label = 'rand2',row.names = F,
      caption = 
        "Rand and adjusted Rand indices of different clustering methods for Example 5.  
	The methods compared are described in the text.  Eg5 is the case of the original data, Eg5M introduces 
	additional fixed effects in the model (accelerometer data in two other directions) 
	and Eg5G treats 10% of the original observations as missing") %>% 
  kable_styling() %>% 
  pack_rows("Rand Index", 1, 3) %>% 
  pack_rows("Adjusted Rand Index", 4, 6) %>% 
  write( file='rand2.tex')


# # benchmark with composite index --------------------------------------
# 
# library(fpc)
# tb_rand <- matrix(NA, nrow=nData, ncol= 4, dimnames = list(
#   paste0('Eg', seq(nData)), 
#   paste0('PC', seq(4))))
# z0 <- numeric(4)
# 
# clustermethod=c("kmeansCBI","hclustCBI")
# # A clustering method can be used more than once, with different
# # parameters
# clustermethodpars <- list()
# clustermethodpars[[2]] <- list()
# clustermethodpars[[2]]$method <- "average"
# # Last element of clustermethodpars needs to have an entry!
# methodname <- c("kmeans","average")
# 
# # indices:
# # CH
# # ASW
# # Dunn
# # Pearson Γ Prediction strength Bootstab
# # CVNN
# indexname <-  c(
#   'ch',
#   'asw',
#   'dunn',
#   'pearsongamma',
#   'cvnnd', # smaller is better
#   'corrected.rand')
# 
# for(i in seq(nData)){
#   load(paste0('pcOutput',i))
#   load(paste0('dat',i))
#   load(paste0('ls_benchmark',i))
#   
#   dat <- dat %>% 
#     select(ID,Time,Record) %>% 
#     pivot_wider(id_cols='ID', names_from='Time', values_from='Record') %>% 
#     select(-ID) %>% 
#     as.matrix
#   dat1 <- dist(dat)
#   
#   df_groups %>% 
#     select(Group) %>% unlist %>% as.factor %>% as.numeric -> true_cluster
#   
#   nCluster <- length(unique(true_cluster))
#   
#   arr_cluster0 <- array(unlist(ls_benchmark), 
#         dim=c(nrow(df_groups), length(ls_benchmark)))
#   arr_cluster <- array(unlist(pcOutput$ls_clust0), 
#                        dim=c(nrow(df_groups), length(ls_idxA)))
#   arr_cluster <- cbind(arr_cluster0,arr_cluster)
#   
#   z0 <-  mapply(function(j){
#     clust <- arr_cluster[,j]
#     outcstat <- cqcluster.stats(dat1,clust,alt.clustering=true_cluster)
#     # names(outcstat)
#     unlist(outcstat[indexname])
#   }, 1:ncol(arr_cluster))
#   
#   
#   
#   clustermethod=c("kmeansCBI","hclustCBI")
#   # A clustering method can be used more than once, with different
#   # parameters
#   clustermethodpars <- list()
#   clustermethodpars[[2]] <- list()
#   clustermethodpars[[2]]$method <- "average"
#   # Last element of clustermethodpars needs to have an entry!
#   methodname <- c("kmeans","average")
#   cbs <-  clusterbenchstats(dat,G=nCluster,clustermethod=clustermethod,
#                             methodname=methodname,distmethod=rep(FALSE,2),
#                             clustermethodpars=clustermethodpars,nnruns=1,kmruns=1,fnruns=1,avenruns=1)
#   
#  # xx <- cbs$qstat
#  # xx$statistics
#  # weights=c(1,0,0,0,0,1,0,1,0,1,0,1,0,0,1,1)
#  # xx$statistics[weights==1]
#  #  print(cbs$qstat,aggregate=TRUE,weights=c(1,0,0,0,0,1,0,1,0,1,0,1,0,0,1,1))
#   tb_rand[i, ] <- z0
# }
# 
# # k means with fitted -----------------------------------------------------
# StatMatch::mahalanobis.dist
# rand.fcn <- rand.index
# source('../7_kmeansFunction.R')
# regQ= 1e-6
# k=4 # the function to integrate clustering results
# nData <- 5
# {
#   pc.rand <- get(paste0('pc.rand',k))
#   
#   tb_rand <- matrix(NA, nrow=nData, ncol= 4, dimnames = list(
#     paste0('Eg', seq(nData)), 
#     paste0('PC', seq(4))))
#   z0 <- numeric(4)
#   ls_dK <- NULL
#   par(mfrow=c(1,2))
#   for(i in seq(nData)){
#     load(paste0('pcOutput',i))
#     load(paste0('dat',i))
#     
#     mat_fitted <- pcOutput$mat_fitted
#     xx <- dat %>% # formatting observations for distance
#       select(ID,t) %>% 
#       bind_cols(data.frame(y=mat_fitted[,1])) %>% 
#       pivot_wider(names_from = t, values_from=y) %>% 
#       select(-ID) %>% 
#       as.matrix
#     # xx <- as.matrix(iris[,-5])
#     # vc <- diag(1, nrow=ncol(xx))
#     vc <- cov(xx) +diag(regQ,ncol(xx))
#     temp <- Kmeans.simu(xx, cov.func, dist.func, nIter=10,nK=30, seed)
#     dK <- temp$mse/ncol(xx)
#     plot(dK,type='l', main=paste('Example',i))
#     # jK <- diff(dK^(-ncol(xx)/2))
#     jK <- diff(dK^(-2))
#     plot(jK,type='l', main=paste('Example',i))
#     
#     nCluster <- which(jK==max(jK,na.rm=T))
#     ls_dK[[i]] <- list(dK, jK, nCluster)
#     # tb_rand[i,1] <- rand.fcn(true_cluster,temp$clusterMatrix[,nCluster])
#   }
# }
# 
# dK <- pcOutput$KLs[[1]][1:50]/nSub/nBasis
# plot(dK,type='l', main=paste('Example',i))
# jK <- diff(dK^(-5))
# plot(jK,type='l', main=paste('Example',i))
