# https://physionet.org/content/capslpdb/1.0.0/
library(tidyverse)
library(grid)
library(fda)

filter <- dplyr::filter
select <- dplyr::select

#' ### Import dataset
#' 512 ticks per second; Signal: Fp2-F4	1 tick per sample; 32.76 adu/uV; 12-bit ADC, zero at 0; baseline is 0
path <- 'D:/Data/sleep/capslpdb/'
file <- 'brux2.mat'
dat0 <- R.matlab::readMat(file.path(path,'brux2.mat'))$temp %>% as.data.frame
dim(dat0)
head(dat0)
colnames(dat0) <- c('t','Record')
# segments of 30 seconds
set.seed(1)
dat <- dat0 %>%
  mutate(segment=floor(t/30)) %>%
  filter(segment%in%sample(seq(max(segment)), 100, replace=F))

dat %>% count(segment)



# hypnogram (0=W=wake, S1-S4=sleep stages, 5=R=REM, 7=MT=body movements)
hyp <- R.matlab::readMat(file.path(path,paste0('hyp_',file)))$hyp
head(hyp)
colnames(hyp) <- c('Group','t')
max(hyp[,'t'])

# merge eeg record with sleep types
hyp <- hyp %>%
  as.data.frame %>%
  mutate(gp1= lag(Group, default=0), cut=gp1!=Group, stack=cumsum(cut))
hyp %>%
  group_by(stack) %>%
  summarise(t1=min(t),t2=max(t),rg= diff(range(t))) %>%
  # filter(rg>100) %>% # find segments of hypnogram with >100s window
  left_join(hyp %>% select(Group,stack) %>% unique, by='stack') -> tb_group
tb_group %>% count(Group)
tb_group %>%
  filter(rg>0) %>%
  # group_by(Group) %>%
  # slice(1:10) %>%
  # ungroup %>%
  mutate(idx=1:n(),
         Group=as.character(Group)) -> tb_group
tb_group %>%
  count(Group)

group.find <- function(vals){
  x <- mapply(function(a,b) between(vals, a, b),
              tb_group$t1, tb_group$t2)
  if(sum(x)==0) return(NA)
  else return(which(x))
}
dat1 <- dat %>%
  rowwise() %>%
  mutate(idx=group.find(t)) %>%
  ungroup
dat <- dat1 %>%
  filter(!is.na(idx), !is.na(Record)) %>%
  left_join(tb_group, by='idx') %>%
  mutate(Subject=segment)

dat <- dat %>%
  filter(!segment%in% (dat %>% count(segment) %>% filter(n==1) %>% select(segment) %>% unlist))

tb_group1 <- dat %>%
  group_by(Subject) %>%
  summarise(t1=min(t),t2=max(t)) %>%
  left_join(dat %>% select(Subject, Group) %>% unique, by='Subject')

fs <- 256 # sampling frequency
hz <- 40 # look at frequency spectrum for below hz
x_hz <- seq(hz)
w_hz <- diff(x_hz)[1]

# fast fourier transformation segments for frequency spectrum
mapply(function(low,high){
  cat('Calculating time interval', low,high,'\n')
  dat %>%
    filter(between(t, low, high)) %>%
    select(Record) %>%
    unlist->x
  y <- fft(x)
  n <- length(y)
  t_power = (0:(n-1))*(fs/n)     # frequency range
  power = abs(y)^2/n    # power of the DFT
  idx <- which(t_power<=hz)
  # plot(t_power[idx],power[idx], type='l')
  y1 <- sapply(x_hz, function(h){
    mean(power[between(t_power, h-w_hz/2, h+w_hz/2)], na.rm = T)
  })
  # plot(x_hz,y1, type='b')
  return(y1)
}, tb_group1$t1, tb_group1$t2) -> temp

lineColors <- colorRamp(c("red", "blue"))(seq(0, 1, len = 6))
plot(temp[,1], type='l', col=lineColors[as.numeric(tb_group1$Group[1])+1,])
for(i in 1:ncol(temp)) lines(temp[,i], col=lineColors[as.numeric(tb_group1$Group[i])+1,])

dat <- tb_group1 %>%
  bind_cols(as.data.frame(t(temp))) %>%
  select(Subject, Group, V1:V30) %>%
  pivot_longer(cols=V1:V30, names_to='Herz', values_to='Power') %>%
  mutate(Herz=as.numeric(gsub('V','',Herz)), Group=as.factor(Group)) %>%
  mutate(Record=as.numeric(scale(log(Power))),
         t=(Herz-min(Herz))/diff(range(Herz)))
dat <- na.omit(dat)
dat %>%
  ggplot(aes(x=t, y=Record, group=Subject, color=Group)) +
  geom_line() +
  facet_wrap(~Group)


df_groups <- dat %>%
  select(Subject, Group) %>%
  unique %>%
  mutate(ID=row_number()) %>%
  arrange(Group) %>%
  mutate(Sub=row_number())

dat %>% count(Subject)
dat <- dat %>%
  left_join(df_groups %>% select(-Group), by='Subject')


# # bspline basis
basis.create <- function(nBasis, t){
  obj <- create.bspline.basis(nbasis = nBasis)
  # obj <- create.fourier.basis(nBasis = 6)
  # plot(obj)
  fitted <- predict(obj, t)
  # plot(dat$t,fitted[,6],type='l')
  fitted
}
gr <- dat$t[1:30]
a <- basis.create(30, gr)
plot(gr, a[,1],type='l')
for(j in seq(2,6))
  lines(gr,a[,j])

nBasis <- 31
dat$Z1 <- 1
dat[,paste0('Z',2:nBasis)] <- basis.create(nBasis-1, dat$t)


ls_idxA <- list(
  seq(nBasis),
  1:6,
  7:16,
  17:31
)

for(i in seq_along(ls_idxA)){
  idx <- ls_idxA[[i]]
  plot(dat$t[1:30], dat$Z2[1:30],type='l')
  for(j in idx[-1])
    lines(dat$t[1:30], unlist(dat[1:30,paste0('Z',j)]))
}



png('plotSeriesEg4.png',height=8,width=12,units='cm',res=300,pointsize=10)
dat %>%
  ggplot(aes(x=t,y=Record,group=ID,color=as.factor(Group))) +
  geom_line()+
  facet_wrap(~factor(Group))+
  theme_bw()+
  theme(legend.background = ,legend.position = "none")+
  scale_x_continuous('Frequency', breaks=seq(0,1,.2),labels=c('0',seq(.2,.8,.2),'1'))
dev.off()

nIter <- 10
seed <- 100
nDraw <- 1000
thKL <- 1/10 # threshold of KL to choose cluster number

# close by groups  (0=W=wake, S1-S4=sleep stages, 5=R=REM, 7=MT=body movements)
ls_groups <- list(
  '0'= c('1'),
  '1'= c('0','2'),
  '2'= c('1','3'),
  '3'= c('2','4'),
  '4'= c('3','5'),
  '5'= c('4')
)
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
save(ls_par, mat_fitted, nClusters, out_boot, out_pc,file='pcOutputEg4')


out_pc0 <- pc.pair(df_of_draws, dat, nBasis, ls_idxA, nIter, nDraw, length(unique(df_groups$Group)), regQ, seed)
save(mat_fitted, out_pc0,file='pcOutputEg4Fixed')

save(df_of_draws,file='D:/df_of_drawsEg4')

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
png('plotKmeansEg4.png',height=10,width=15,units='cm',res=150,pointsize=10)
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

png('plotClusterProbEg4.png',height=15,width=20,units='cm',res=300,pointsize=10)
cluster.plot(ls_prob, df_groups, ls_groups, cutPoints, fileName, subTitles, allColors,seed=123,thin=1,ls_labels=NULL)
dev.off()
