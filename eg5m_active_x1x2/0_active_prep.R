# https://archive.ics.uci.edu/ml/datasets/Activity+recognition+with+healthy+older+people+using+a+batteryless+wearable+sensor

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
         Record=V4,
         x1=V2,
         x2=V3
         )
dat %>% count(Subject)
dat <- dat %>%
  filter(
    Group>0,
    Subject%in% (
    dat %>% count(Subject) %>% filter(n>10) %>% select(Subject) %>% unlist)) %>%
  select(Subject, Group, t, Record, x1, x2) %>%
  group_by(Subject) %>%
  slice(1:50) %>%
  mutate(t= t-min(t)) %>%
  ungroup %>%
  mutate(t=t/(max(t)-min(t)),
         Record=as.numeric(scale(Record)),
         x1= as.numeric(scale(x1)),
         x2= as.numeric(scale(x2)))
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
ls_groups <- NULL
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
dat$Z32 <- dat$x1
dat$Z33 <- dat$x2
nBasis <- 33

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
  c(1:6,32:33),
  c(7:11,32:33),
  12:33
)
nIter <- 10
seed <- 100
nDraw <- 1000
thKL <- 1/5 # threshold of KL to choose cluster number

save(nBasis, nIter, seed, nDraw, thKL,  df_groups, ls_groups, dat,ls_idxA, file='dat')


# linear mixed model ------------------------------------
source('../1_main.R')

out_pc0 <- pc.pair(df_of_draws, dat, nBasis, ls_idxA, nIter, nDraw, length(unique(df_groups$Group)), regQ, seed)
save(mat_fitted, out_pc0,file='pcOutputEg6Fixed')

save(df_of_draws,file='D:/df_of_drawsEg6')

