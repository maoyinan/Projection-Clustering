library(tidyverse)
library(grid)

nClusters <- 4
nSub <- 30*nClusters # subjects
nT <- 40 # time points per subject
dat <- data.frame(
  ID=rep(seq(nSub),each=nT),
  Group=rep(seq(nClusters),each=nSub/nClusters*nT),
  Time=rep(seq(nT),nSub)
) %>%
  mutate(t= (Time-min(Time))/diff(range((Time))) )
dat$Z1 <- 1
for(i in c(1:9)){
  dat[,paste0('Z',i+1)] <- cos(pi*i*dat$t)
}
nBasis <- 10

# simulate
set.seed(12345)
epsilon <- rnorm(nSub*nT,0,1/10)
i1wi <- sample(c(1:3),nSub,replace = T)+1
i2wi <- sample(c(7:9),nSub,replace = T)+1

dat$Record <- NA
for(j in 1:nSub){# filter ID==j
  idx <- which(dat$ID==j)
  gr <- dat$Group[idx[1]] # which cluster determines beta
  beta1 <- ifelse(gr%in%c(1,2),1,.1)
  beta2 <- ifelse(gr%in%c(1,3),1,.1)
  z1 <- dat[idx,paste0('Z',i1wi[j])]
  z2 <- dat[idx,paste0('Z',i2wi[j])]
  yi <- beta1*z1+beta2*z2 + epsilon[idx]
  dat$Record[idx] <- yi
}
dat$Record <- as.numeric(scale(dat$Record))


ls_groups <- list(
  '1'= c('2','3'),
  '2'= c('1','4'),
  '3'= c('1','4'),
  '4'= c('3','2')
)

png('plotSeriesEg1.png',height=8,width=12,units='cm',res=300,pointsize=10)
dat %>%
  ggplot(aes(x=t,y=Record,group=ID,color=as.factor(Group))) +
  geom_line()+
  facet_wrap(~factor(Group))+
  theme_bw()+
  theme(legend.background = ,legend.position = "none")+
  scale_x_continuous('Time', breaks=seq(0,1,.2),labels=c('0',seq(.2,.8,.2),'1'))
dev.off()

ls_idxA <- list(
  seq(nBasis),
  1:4,
  5:7,
  8:10
)
nIter <- 10
seed <- 100
nDraw <- 1000
thKL <- 1/10 # threshold of KL to choose cluster number

df_groups <- dat %>%
  mutate(Sub=ID) %>%
  select(Sub,ID,Group) %>%
  unique

save(nBasis, nIter, seed, nDraw, thKL,  df_groups, ls_groups, dat,ls_idxA, file='dat')


# linear mixed model ------------------------------------
source('../1_main.R')

# manually chosen nClusters
nClusters <- 4 # move choose cluster code here
ls_prob <- pc.pair(df_of_draws, dat, nBasis, ls_idxA, nIter, nDraw,nClusters, regQ, seed)$ls_prob

save(ls_par, mat_fitted, nClusters, ls_prob,file='pcOutputEg1')

save(df_of_draws,file='D:/df_of_drawsEg1')


# plots -------------------------------------------------------------------
load('pcOutputEg1')
# Plot connection probability
subTitles <-
  c('All frequencies','Low frequencies','Medium frequencies','High frequencies')

allColors <- c('darkseagreen', 'coral3','steelblue','sandybrown')
cutPoints <- c(.5,.8)

# Rearrange subjects by classes so that subjects in same group are put together for easy visulisation later
for(i in seq_along(ls_prob)){
  ls_prob[[i]] <- ls_prob[[i]][df_groups$ID,df_groups$ID]
}
ls_labels <- rep(list(1:4),4)
ls_labels[[4]] <- c(1,3,2,4)
png('plotClusterProbEg1.png',height=15,width=20,units='cm',res=300,pointsize=10)
cluster.plot(ls_prob,df_groups,ls_groups,cutPoints, fileName, subTitles, allColors, seed=123,thin=1, ls_labels)
dev.off()

# tabulate links count
tb2 <- NULL
nSub <- nrow(df_groups)
for(t in seq_along(ls_prob)){
  tb <- ls_prob[[t]]
  tb1 <- data.frame(as.matrix(tb,dimnames=list(1:nSub,1:nSub))) %>%
    bind_cols(row=1:nSub) %>%
    pivot_longer(-row,names_to='col',values_to='prob') %>%
    mutate(col=as.numeric(gsub('X','',col))) %>%
    filter(!is.na(prob)) %>%
    filter(prob>.8) %>%
    mutate(row=df_groups$Group[row],col=df_groups$Group[col],
           Link=paste(row,col,sep='-'))
  tb1 <- tb1 %>% count(row,col) %>%
    pivot_wider(names_from='col',values_from='n')
  if(nrow(tb1)>0) tb2 <- tb2 %>%
    bind_rows(bind_cols(PC=rep(t,nrow(tb1)),tb1))
}
# tb2 <- tb2[,-1]
# colnames(tb2)[1] <- 'Cluster'
#
# options(knitr.kable.NA = '')
# kable(tb2, "latex",digits = 0,
#       label = 'links0',row.names = F,
#       caption =
#         "Frequency table of clustering for links higher than 0.8 probability, based on four projection clustering schemes") %>%
#   kable_styling() %>%
#   pack_rows("PC1", 1, 4) %>%
#   pack_rows("PC2", 5, 8) %>%
#   # pack_rows("PC3", 11, 15) %>%
#   pack_rows("PC4", 9, 12) %>%
#   write( file='links_eg1.tex')
#

