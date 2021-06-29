# https://archive.ics.uci.edu/ml/datasets/Activity+Recognition+from+Single+Chest-Mounted+Accelerometer

library(tidyverse)
filter <- dplyr::filter

#' ### Import dataset
path <- '/Users/maoyinan/Documents/Activity/'
path <- 'D:/Data/Activity/'

dat0 <- NULL
for(i in seq(15)){
  dat1 <- read.csv(file.path(path,sprintf('%d.csv',i)), header=F)
  dat1 <- dat1 %>%
    group_by(V1) %>%
    slice(1) %>%
    ungroup
  dat0 <- rbind(dat0,
                dat1)
}
dat <- dat0 %>%
  mutate(Subject=c(0,diff(V5)),
         Subject=ifelse(Subject!=0, 1, 0),
         Subject=cumsum(Subject)+1,
         Group= V5,
         t=V1,
         Record=V4)
dat %>% count(Subject)
dat <- dat %>%
  filter(
    Group>0,
    Subject%in% (
    dat %>% count(Subject) %>% filter(n>10) %>% select(Subject) %>% unlist)) %>%
  select(Subject, Group, t, Record) %>%
  group_by(Subject) %>%
  slice(1:50) %>%
  mutate(t= t-min(t)) %>%
  ungroup %>%
  mutate(Record=as.numeric(scale(Record)),
         t=t/(max(t)-min(t)))
dat %>% group_by(Subject) %>% summarise(max(t)) %>% View
dat <- dat %>% filter(
  Subject%in%(dat %>% group_by(Subject) %>% summarise(x=max(t)) %>% filter(x==0.1) %>% select(Subject) %>% unlist)) %>%
  group_by(Subject) %>%
  mutate(t= t-min(t)) %>%
  ungroup %>%
  mutate(t=t/(max(t)-min(t)))


#' Check data vs grouping
dat %>%
  ggplot(aes(x=t,y=Record,group=Subject,color=Subject)) +
  geom_line()+
  facet_wrap(~Group) +
  theme_bw()

groups <- sort(unique(dat$Group))
# --- 1: Working at Computer
# --- 2: Standing Up, Walking and Going updown stairs
# --- 3: Standing
# --- 4: Walking
# --- 5: Going UpDown Stairs
# --- 6: Walking and Talking with Someone
# --- 7: Talking while Standing
ls_groups <- list(
  '1'= c('2'),
  '2'= c('1','3'),
  '3'= c('2','7'),
  '4'= c('2','6'),
  '5'= c('4','6'),
  '6'= c('2','4'),
  '7'=c('3')
)
df_groups <- data.frame(Group=groups) %>%
  left_join(
    dat %>%
      select(Group,Subject) %>%
      unique %>%
      mutate(ID=row_number())
  ) %>%
  mutate(Sub=row_number())

df_groups %>%
  count(Group)

dat <- dat %>%
  left_join(df_groups %>% select(-Group), by='Subject')

#' Prepare parameters for projection clustering
dat$Z1 <- 1
for(i in 1:30){
  dat[,paste0('Z',i+1)] <- cos(pi*i*dat$t)
}
nBasis <- 31 # number of random effects

png('plotSeriesEg5.png',height=8,width=12,units='cm',res=300,pointsize=10)
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
  1:6,
  7:11,
  12:31
)
nIter <- 10
seed <- 100
nDraw <- 1000
thKL <- 1/10 # threshold of KL to choose cluster number

save(nBasis, nIter, seed, nDraw, thKL,  df_groups, ls_groups, dat,ls_idxA, file='dat')



# linear mixed model ------------------------------------
source('../1_main.R')

# nClusters chosen by bootstrap method
out_boot <- nCluster.boot(mat_fitted, dat, nB=100, nCl=30, seed=1)
nClusters <- out_boot$nClusters

out_pc <- pc.pair(df_of_draws, dat, nBasis, ls_idxA, nIter, nDraw,nClusters, regQ, seed)
ls_prob <- out_pc$ls_prob
ls_clust0 <- out_pc$ls_clust0
save(ls_par, mat_fitted, nClusters,out_boot, out_pc,file='pcOutputEg5')


out_pc0 <- pc.pair(df_of_draws, dat, nBasis, ls_idxA, nIter, nDraw, length(unique(df_groups$Group)), regQ, seed)
save(mat_fitted, out_pc0,file='pcOutputEg5Fixed')

save(df_of_draws,file='D:/df_of_drawsEg5')

# plots -------------------------------------------------------------------

# Plot bootstrap cluster number
xText <- 'K: number of clusters'
yText <- 'Instability'
plotRow <- 2
ls_sB <- out_boot
subTitles <- c('All frequencies','Low frequencies','Medium frequencies','High frequencies')
png('plotKmeansEg5.png',height=10,width=15,units='cm',res=150,pointsize=10)
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

png('plotClusterProbEg5.png',height=15,width=20,units='cm',res=300,pointsize=10)
cluster.plot(ls_prob, df_groups, ls_groups, cutPoints, fileName, subTitles, allColors,seed=123,thin=1,ls_labels=NULL)
dev.off()
png('plotClusterProb10percEg5.png',height=15,width=20,units='cm',res=300,pointsize=10)
cluster.plot(ls_prob, df_groups, ls_groups, cutPoints, fileName, subTitles, allColors,seed=123,thin=.1,ls_labels=NULL)
dev.off()


