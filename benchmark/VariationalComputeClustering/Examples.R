
#-----------------------------#
# Example on time course data #
#-----------------------------#
 
orf <- read.table(file='ORF_DATA.txt',sep="\t",header=TRUE,fill=TRUE)
dim(orf)
genes <- as.character(orf[1:800,1])
genes <- sort(genes)
which(genes=="YEL077C")
genes <- genes[-185]
which(genes=="YKL096W-A")
genes[427] <- "YKL097W-A"

yeast <- read.table(file='CDC_DATA.txt',sep="\t",header=TRUE,fill=TRUE)
yeast <- yeast[, -(2:5)]
yeast <- yeast[, -(20:74)]
dim(yeast)
vec <- rep(0,799)
for (i in 1:799){
vec[i] <- which(yeast[,1]==genes[i])}
yeast <- yeast[vec,]
yeast <- na.omit(yeast)
dim(yeast)

# construct y #
n <- 612
vni <- rep(18,n)
nsum <- sum(vni)
data <- matrix(0,nrow=612,ncol=18)
for (j in 1:18) {
data[,j] <- yeast[,j+1]}
y <- rep(0,nsum)
for (i in 1:n){
y[(18*(i-1)+1):(18*i)] <- data[i,]}

# Covariates #
X <- matrix(0,nrow=nsum,ncol=2)
W <- matrix(1,nrow=nsum,ncol=1)
V <- matrix(0,nrow=nsum,ncol=18)
U <- matrix(1,nrow=n,ncol=1)

# construct X,W,V #
Xi <- matrix(0,nrow=18,ncol=2)
for (l in 0:17) {Xi[(l+1),] <- c(cos(7*l*2*pi/53),sin(7*l*2*pi/53))}
for (i in 1:n) { 
X[((i-1)*18+1):(i*18),] <- Xi
V[((i-1)*18+1):(i*18),] <- diag(18)  }
epm <- matrix(18,nrow=n,ncol=1)


# Run VGA using Algorithm 1(no centering) #
M <- 5
IGa <- 2
IGb <- round(2*(fitdistr(lm(y~X-1)$residuals,'t', m=0,df = 4)$estimate)^2,2)
sigdelp <- 1000
sigbetap <- 1000
set.seed(826)
mixturemodel <- greedy(y,X,W,V,U,vni,epm,M,IGa,IGb,sigdelp,sigbetap)


# Applying merge moves #
fit <- mixturemodel$fitpre
m1 <- 12
m2 <- 13
tol <- 1.0e-5
type <- 'all'
IGa <- 2
IGb <- round(2*(fitdistr(lm(y~X-1)$residuals,'t', m=0,df = 4)$estimate)^2,2)
sigdelp <- 1000
sigbetap <- 1000
fitmerge <- MLMMmerge(y,X,W,V,U,IGa,IGb,sigbetap,sigdelp,vni,epm,tol,fit,m1,m2,type)
           


#------------------------------------------#
# Example on completely synthetic data set #
#------------------------------------------#

synthetic <- read.table(file='syn_sine_5_mult1.txt',header=TRUE)
y <- matrix(t(as.matrix(synthetic[4:83])),nrow=32000,ncol=1,byrow=FALSE)
y <- as.vector(y)
n <- 400
vni <- rep(80,400)

X <- NULL
for (i in 1:n) {
Xi <- matrix(0,nrow=80,ncol=20)
count <- 0
for (t in 1:20) {
 for (r in 1:4) {
  count <- count+1
  Xi[count,t] <- 1 } }
X <- rbind(X,Xi)}
W <- X

V <- NULL
for (i in 1:n) {V <- rbind(V,diag(80))}
U <- matrix(1,n,1)
epm <- matrix(4,n,20)


# Run VGA using Algorithm 2 (partial centering) #
M <- 5
IGa <- 2
IGb <- round(2*(fitdistr(lm(y~X-1)$residuals,'t',m=0,df=4)$estimate)^2,2)
sigdelp <- 1000
sigbetap <- 1000
set.seed(3)
mixturemodel <- greedyXW(y,X,W,V,U,vni,epm,M,IGa,IGb,sigdelp,sigbetap)


# Applying merge moves #
fit <- mixturemodel$fitpre
m1 <- 5
m2 <- 6
tol <- 1.0e-5
type <- 'all'
IGa <- 2
IGb <- round(2*(fitdistr(lm(y~X-1)$residuals,'t',m=0,df=4)$estimate)^2,2)
sigdelp <- 1000
sigbetap <- 1000
fitmerge <- MLMMpmerge(y,X,W,V,U,IGa,IGb,sigbetap,sigdelp,vni,epm,tol,fit,m1,m2,type)
             


#-------------------------------------#
# Example on yeast galactose data set #
#-------------------------------------#

gbm <- read.table(file='gal205.txt',sep="\t",header=TRUE)
gbm <- gbm[,-c(1,2)]
gbmnew <- as.matrix(gbm[,2:81])

n <- 205
datam <- array(gbmnew[,1:80], dim=c(n,4,20))
data <- array(0, dim=c(n,20,4))
for (i in 1:n) {
data[i,,] <- t(datam[i,,])}

vni <- rep(0,n)
y <- NULL
for (i in 1:n) {
vv <- as.vector(na.omit(gbmnew[i,1:80]))
y <- c(y,vv)
vni[i] <- length(vv)}

X <- NULL
for (i in 1:n) {
Xi <- matrix(0,nrow=vni[i],ncol=20)
count <- 0
 for (t in 1:20) {
  for (r in 1:4) {
  if (is.finite(data[i,t,r])) {
   count <- count+1
   Xi[count,t] <- 1 } }}
X <- rbind(X,Xi)}
W <- X

V <- NULL
for (i in 1:n) {
Vi <- matrix(0,nrow=vni[i],ncol=80)
count <- 0
 for (l in 1:80) {
  if (is.finite(gbm[i,(l+1)])) {
   count <- count+1
   Vi[count,l] <- 1 } }
V <- rbind(V,Vi)}

U <- matrix(,nrow=n,ncol=4)
for (i in 1:n){
for (j in 1:4){
if (gbm[i,1]==j) {U[i,j] <- 1} else U[i,j] <- 0 }}

epm <- matrix(0,nrow=n,ncol=20)
for (i in 1:n) {
for (t in 1:20) {
epm[i,t] <- length(na.omit(data[i,t,]))}}


# Applying VGA using Algorithm 2 (partial centering) #
M <- 5
sigdelp <- 1000
sigbetap <- 1000
IGa <- 2
IGb <- round(2*(fitdistr(lm(y~X-1)$residuals,'t',m=0,df=4)$estimate)^2,2)
set.seed(2)
mixturemodel <- greedyXW(y,X,W,V,U,vni,epm,M,IGa,IGb,sigdelp,sigbetap)


# Applying merge moves #
fit <- mixturemodel$fitpre
m1 <- 4
m2 <- 7
tol <- 1.0e-5
type <- 'all'
IGa <- 2
IGb <- round(2*(fitdistr(lm(y~X-1)$residuals,'t',m=0,df=4)$estimate)^2,2)
sigdelp <- 1000
sigbetap <- 1000
fitmerge <- MLMMpmerge(y,X,W,V,U,IGa,IGb,sigbetap,sigdelp,vni,epm,tol,fit,m1,m2,type)
             

#-----------------------------------#
# Example on water temperature data #
#-----------------------------------#

obs <- read.csv(file='temperature_data.csv', header=TRUE)
str(obs)
names(obs)
on <- length(unique(obs$Date))
obs$oday <- seq(1,on,1)
dateref <- obs[,-c(2:12)]
obsdata <- obs[,-c(1,13)]
data <-matrix(nrow=dim(obsdata)[1],ncol=dim(obsdata)[2])
for (j in 1: dim(obsdata)[2]) {
data[,j] <- obsdata[,j]}
obsdepth <- c(0.5,2,4,6,8,10,12,14,16,18,'bot')
plot(seq(1,11,1),data[1,],type='l',ylim=c(min(data),max(data)))
for(i in 1:on){points(seq(1,11,1),data[i,],type='l')}

# construct y #
n <- dim(data)[1]
vni <- rep(11,n)
nsum <- sum(vni) 
y <- NULL
for (i in 1:n) {
vv <- data[i,1:11]
y <- c(y,vv)}

# Covariates #
X <- NULL
for (i in 1:n) {
Xi <- diag(11)
X <- rbind(X,Xi)}
W <- X
V <- X
time <- seq(1,n,1)
U <- cbind(1,time,time^2,time^3)
scle<-function(x) { -1+2*(x-min(x))/(max(x)-min(x)) }
U[,2:4]<- apply(U[,2:4],2,scle)
epm <- matrix(1,nrow=on,ncol=11)


# Applying VGA using Algorithm 3 (full centering) #
M <- 5
sigdelp <- 1000
sigbetap <- 10000
IGa <- 2
IGb <- round(2*(fitdistr(lm(y~X-1)$residuals,'t',m=0,df=4)$estimate)^2,2)
set.seed(2)
mixturemodel <- greedyXWV(y,X,W,V,U,vni,epm,M,IGa,IGb,sigdelp,sigbetap)
mixturemodel$dur

# Applying merge moves #
fit <- mixturemodel$fitpre

xx <- fit$mubetaq
for(i in 1:ncol(xx)){points(seq(1,11,1),xx[,i],type='l',col='red')}
xx <- fit$murhoq
for(i in 1:nrow(xx)){points(seq(1,11,1),xx[i,],type='l',col='blue')}
predClusters <- apply(fit$qp,1,function(x)which(x==max(x)))
plot(seq(1,11,1),data[1,],type='l',ylim=c(min(data),max(data)))
for(i in 1:on){points(seq(1,11,1),data[i,],type='l',col=c('red','blue','green','purple','orange')[predClusters[i]])}


m1 <- 1
m2 <- 2
tol <- 1.0e-5
type <- 'all'
sigdelp <- 1000
sigbetap <- 10000
IGa <- 2
IGb <- round(2*(fitdistr(lm(y~X-1)$residuals,'t',m=0,df=4)$estimate)^2,2)
fitmerge <- MLMMfmerge(y,X,W,V,U,IGa,IGb,sigbetap,sigdelp,vni,epm,tol,fit,m1,m2,type)
          


        


