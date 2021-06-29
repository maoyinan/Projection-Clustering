library(tidyverse)
library(grid)

#' Downloaded from UCR archive https://www.cs.ucr.edu/%7Eeamonn/time_series_data_2018/
#' Random draw 30 series from 5 clusters

dat1 <- read.csv("Crop_TRAIN.tsv", sep = "\t", header = F)
groups <- table(dat1$V1)
set.seed(123)
groups <- sort(sample(seq_along(groups), 5, replace=FALSE))
dat1 <- dat1 %>%
  filter(V1%in%groups) %>%
  group_by(V1) %>%
  slice(1:30)

dat <- dat1 %>%
  ungroup() %>%
  mutate(ID=row_number()) %>%
  gather('Time','Record',V2:V47) %>%
  mutate(Time=as.numeric(gsub('V','',Time))-1) %>%
  arrange(V1,ID,Time) %>%
  mutate(Group=factor(V1),
         t= Time/max(Time),
         Record= as.numeric(scale(Record))) %>%
  select(-V1) %>%
  as.data.frame()

nBasis <- 31 # number of random effects
dat$Z1 <- 1
for(i in 1:(nBasis-1)){
  dat[,paste0('Z',i+1)] <- cos(pi*i*dat$t)
}

png('plotSeriesEg2.png',height=8,width=12,units='cm',res=300,pointsize=10)
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
  5:11,
  12:31
)
nIter <- 10
seed <- 100
nDraw <- 1000
thKL <- 1/10 # threshold of KL to choose cluster number

df_groups <- dat %>%
  mutate(Sub=ID) %>%
  select(Sub,ID,Group) %>%
  unique
ls_groups <- NULL
save(nBasis, nIter, seed, nDraw, thKL,  df_groups, ls_groups, dat,ls_idxA, file='dat')

# linear mixed model ------------------------------------
source('../1_main.R')

# nClusters chosen by KL method
out_KL <- pc.KL(ls_par, dat, ls_idxA, nIter, thKL, regQ, seed=1)
KLs <- out_KL$KLs

# nClusters chosen by bootstrap method
out_boot <- nCluster.boot(mat_fitted, dat, nB=100, nCl=30, seed=1)
nClusters <- out_boot$nClusters

out_pc <- pc.pair(df_of_draws, dat, nBasis, ls_idxA, nIter, nDraw,nClusters, regQ, seed)
ls_prob <- out_pc$ls_prob
ls_clust0 <- out_pc$ls_clust0
save(ls_par, mat_fitted, nClusters,out_KL, out_boot, out_pc,file='pcOutputEg2')

out_pc0 <- pc.pair(df_of_draws, dat, nBasis, ls_idxA, nIter, nDraw, length(unique(df_groups$Group)), regQ, seed)
save(mat_fitted, out_pc0,file='pcOutputEg2Fixed')


save(df_of_draws,file='D:/df_of_drawsEg2')

# plots -------------------------------------------------------------------

# Check fitted curve with projection onto different dimensions

xText <- 'Time'
yText <- 'Record'
plotRow <- length(unique(df_groups$Group))
colNames <- c('Record',colnames(mat_fitted)) # the columns to plot fitted curve
plotColors <- c("#030303", "#FA8072", "#87CEFA", "#DDA0DD", "#9ACD32")
legendLabels <-  c(
  'a'='Observation',
  'c'='All frequencies',
  'd1'='Low frequencies',
  'd2'='Medium frequencies',
  'd3'='High frequencies')

png('plotFitEg2.png',height=20,width=24,units='cm',res=300,pointsize=10)
fitted.plot(dat,mat_fitted,colNames, plotColors, legendLabels, plotRow,xText,yText, nPerGroup=4)
dev.off()


# Plot k means KL

subTitles <-
  c('All frequencies','Low frequencies','Medium frequencies','High frequencies')
xText <- 'K: number of clusters'
yText <- 'KL'
plotRow <- 2

png('plotKLEg2.png',height=12,width=20,units='cm',res=300,pointsize=10)
kmeans.plot(KLs,thKL,plotRow, subTitles,xText,yText)
dev.off()

# Plot bootstrap cluster number
xText <- 'K: number of clusters'
yText <- 'Instability'
plotRow <- 2
ls_sB <- out_boot
subTitles <- c('All frequencies','Low frequencies','Medium frequencies','High frequencies')
png('plotKmeansEg2.png',height=10,width=15,units='cm',res=150,pointsize=10)
kmeans_fitted.plot(ls_sB,plotRow, subTitles,xText,yText)
dev.off()


# Plot connection probability
subTitles <-
  c('All frequencies','Low frequencies','Medium frequencies','High frequencies')

allColors <- c('darkseagreen', 'coral3','steelblue','sandybrown')
cutPoints <- c(.5,.8)

# Rearrange subjects by classes so that subjects in same group are put together for easy visulisation later
for(i in seq_along(ls_prob)){
  ls_prob[[i]] <- ls_prob[[i]][df_groups$ID,df_groups$ID]
}
ls_labels <- rep(list(c('3',14,19,10,15)),4)
png('plotClusterProbEg2.png',height=15,width=20,units='cm',res=300,pointsize=10)
cluster.plot(ls_prob, df_groups, ls_groups, cutPoints, fileName, subTitles, allColors,seed=123,thin=1,ls_labels)
dev.off()
png('plotClusterProb10percEg2.png',height=15,width=20,units='cm',res=300,pointsize=10)
cluster.plot(ls_prob, df_groups, ls_groups, cutPoints, fileName, subTitles, allColors,seed=123,thin=.1,ls_labels)
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
tb2 <- tb2[,-1]
colnames(tb2)[1] <- 'Cluster'

options(knitr.kable.NA = '')
kable(tb2, "latex",digits = 0,
      label = 'links1',row.names = F,
      caption =
        "Frequency table of clustering for links higher than 0.8 probability, based on four projection clustering schemes") %>%
  kable_styling() %>%
  pack_rows("PC1", 1, 5) %>%
  pack_rows("PC2", 6, 10) %>%
  pack_rows("PC3", 11, 15) %>%
  pack_rows("PC4", 16, 20) %>%
  write( file='linksEg2.tex')

