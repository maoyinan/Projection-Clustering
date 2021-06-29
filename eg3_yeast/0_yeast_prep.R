# http://genome-www.stanford.edu/cellcycle/
# https://arxiv.org/pdf/1112.4675.pdf
library(tidyverse)
# library(R.matlab)
# library(readtext)
library(stringr)
library(grid)
filter <- dplyr::filter

#' ### Import dataset
#' 512 ticks per second; Signal: Fp2-F4	1 tick per sample; 32.76 adu/uV; 12-bit ADC, zero at 0; baseline is 0
path <- '/Users/maoyinan/Documents/DataCellCycle'
path <- 'D:/Data/DataCellCycle'
dat0 <- read.csv(file.path(path,'CellCycle98.csv'))
dat <- NULL
for(i in seq(18)){
  file <- sprintf('cycle/%d.txt',i)
  dat1 <- read.table(file.path(path,file), header=F, sep = '\t',quote="",  skip=2)
  colnames(dat1)[3] <- 'ORF'
  dat <- bind_rows(
    dat, bind_cols(
      dat0 %>%
        select(Group= Peak,ORF) %>%
        left_join(dat1 %>% select(ORF,Record=V17), by='ORF'),
      data.frame(t=i)))
}
dat <- dat %>%
  filter(ORF%in%(dat %>% count(ORF) %>% filter(n==18&ORF!='') %>% select(ORF) %>% unlist), !is.na(Record)) %>%
  select(Group,ORF) %>%
  unique %>%
  group_by(Group) %>%
  slice(1:30) %>%
  ungroup %>%
  select(ORF) %>%
  left_join(dat, by='ORF')  %>%
  mutate(Record=log(Record),
         t=(t-min(t))/max(t),
         Record=as.numeric(scale(Record)))

png('plotSeriesEg3.png',height=8,width=12,units='cm',res=300,pointsize=10)
dat %>%
  ggplot(aes(x=t,y=Record,group=ID,color=as.factor(Group))) +
  geom_line()+
  facet_wrap(~factor(Group))+
  theme_bw()+
  theme(legend.background = ,legend.position = "none")+
  scale_x_continuous('Time', breaks=seq(0,1,.2),labels=c('0',seq(.2,.8,.2),'1'))
dev.off()

groups <- c('G1','S','S/G2', 'G2/M','M/G1')
ls_groups <- list(
  'G1'= c('M/G1','S'),
  'S'= c('G1','S/G2'),
  'S/G2'= c('S','G2/M'),
  'G2/M'= c('S/G2','M/G1'),
  'M/G1'= c('G2/M','G1')
)
df_groups <- data.frame(Group=groups) %>%
  left_join(
    dat %>%
      select(Group,ORF) %>%
      unique %>%
      mutate(ID=row_number())
  ) %>%
  mutate(Sub=row_number())

df_groups %>%
  count(Group)

dat <- dat %>%
  right_join(df_groups %>%
               select(-Group), by='ORF')
dat %>% count(ID)

#' Prepare parameters for projection clustering
dat$Z1 <- 1
for(i in 1:18){
  dat[,paste0('Z',i+1)] <- cos(pi*i*dat$t)
}
nBasis <- 19 # number of random effects

#' Check observations compared with most frequent Fourier basis
dat %>%
  filter(ID%in%seq(30)) %>%
  ggplot(aes(x=t,y=Record,group=Sub)) +
  geom_line(aes(color=as.factor(Group)))+
  geom_line(aes_string('t', paste0('Z',nBasis)), color='black') +
  facet_wrap(~Sub,nrow=5) +
  theme_bw()

ls_idxA <- list(
  seq(nBasis),
  1:7,
  8:13,
  14:19
)
nIter <- 10
seed <- 100
nDraw <- 1000
thKL <- 1/10 # threshold of KL to choose cluster number

save(nBasis, nIter, seed, nDraw, thKL,  df_groups, ls_groups, dat,ls_idxA, file='dat')


# linear mixed model ------------------------------------
source('../1_main.R')
#
# # nClusters chosen by KL method
# out_KL <- pc.KL(ls_par, dat, ls_idxA, nIter, thKL, regQ, seed=1)
# nClusters <- out_KL$nClusters
# KLs <- out_KL$KLs

# nClusters chosen by bootstrap method
out_boot <- nCluster.boot(mat_fitted, dat, nB=100, nCl=30, seed=1)
nClusters <- out_boot$nClusters


out_pc <- pc.pair(df_of_draws, dat, nBasis, ls_idxA, nIter, nDraw,nClusters, regQ, seed)
ls_prob <- out_pc$ls_prob
ls_clust0 <- out_pc$ls_clust0
save(ls_par, mat_fitted, nClusters, out_boot, out_pc,file='pcOutputEg3')


out_pc0 <- pc.pair(df_of_draws, dat, nBasis, ls_idxA, nIter, nDraw, length(unique(df_groups$Group)), regQ, seed)
save(mat_fitted, out_pc0,file='pcOutputEg3Fixed')

save(df_of_draws,file='D:/df_of_drawsEg3')

# plots -------------------------------------------------------------------

# # Check fitted curve with projection onto different dimensions
#
# xText <- 'Time'
# yText <- 'Record'
# plotRow <- length(unique(df_groups$Group))
# colNames <- c('Record',colnames(mat_fitted)) # the columns to plot fitted curve
# plotColors <- c("#030303", "#FA8072", "#87CEFA", "#DDA0DD", "#9ACD32")
# legendLabels <-  c(
#   'a'='Observation',
#   'c'='All frequencies',
#   'd1'='Low frequencies',
#   'd2'='Medium frequencies',
#   'd3'='High frequencies')
#
# png('plotFitEg1.png',height=20,width=24,units='cm',res=300,pointsize=10)
# fitted.plot(dat,mat_fitted,colNames, plotColors, legendLabels, plotRow,xText,yText, nPerGroup=4)
# dev.off()
#
#
# # Plot k means KL
#
# subTitles <-
#   c('All frequencies','Low frequencies','Medium frequencies','High frequencies')
# xText <- 'K: number of clusters'
# yText <- 'KL'
# plotRow <- 2
#
# png('plotKLEg1.png',height=12,width=20,units='cm',res=300,pointsize=10)
# kmeans.plot(KLs,thKL,plotRow, subTitles,xText,yText)
# dev.off()

# Plot bootstrap cluster number
xText <- 'K: number of clusters'
yText <- 'Instability'
plotRow <- 2
ls_sB <- out_boot
subTitles <- c('All frequencies','Low frequencies','Medium frequencies','High frequencies')
png('plotKmeansEg3.png',height=10,width=15,units='cm',res=150,pointsize=10)
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

png('plotClusterProbEg3.png',height=15,width=20,units='cm',res=300,pointsize=10)
cluster.plot(ls_prob, df_groups, ls_groups, cutPoints, fileName, subTitles, allColors,seed=123,thin=1,ls_labels=NULL)
dev.off()
png('plotClusterProb10percEg3.png',height=15,width=20,units='cm',res=300,pointsize=10)
cluster.plot(ls_prob, df_groups, ls_groups, cutPoints, fileName, subTitles, allColors,seed=123,thin=.1,ls_labels=NULL)
dev.off()
