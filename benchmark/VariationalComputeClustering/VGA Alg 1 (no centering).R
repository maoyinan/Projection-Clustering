library(nnet) # for multinom() function
library(mclust) # for adjustedRandIndex() function
library(R.utils) # for insert() function
library(MASS) # for fitdistr() function
library(micEcon) # for insertRow() function

#function to compute t(x)%*%A^(-1)%*%x
quadinv <- function(x,A) sum(x*solve(A,x))

#function to compute t(x)%*%A%*%x
quad <- function(x,A) sum(x*(A%*%x))

#function to compute trace of A^(-1)%*%B
tri <- function(A,B) sum(diag(solve(A,B))) 

#function to compute trace of A
tr <- function(A) sum(diag(A))

#--------------------------------------------------------------------------#
# Function to check that update of mudelq leads to increase in lower bound #
#--------------------------------------------------------------------------#

Lmudelq <- function(mudelq,qp,sigdel,U){
L <- -0.5*quadinv(mudelq,sigdel)
k <- dim(qp)[2]
d <- dim(U)[2]
pp <- matrix(0,n,k)
ppnew <- matrix(0,n,k)
for (j in 2:k) { pp[,j] <- (U%*%mudelq[((j-2)*d+1):((j-1)*d)]) }
for (j in 1:k) {ppnew[,j] <- 1/rowSums(exp(pp-pp[,j]))}
pp <- ppnew
A <- log(pp/qp)
A[is.na(A)] <- 0
A[is.infinite(A)] <- 0
L <- L +sum(qp*A)
list(L=L) }

#--------------------------------------------------------------------#
# Function for adjustment of lower bound due to normal approximation #
#--------------------------------------------------------------------#

vlbadj <- function(U,alpa,sigdel,sigdelq) {
d <- dim(U)[2]
k <- length(alpa)
lbadj <- 0.5*(k-1)*d*log(2*pi)+0.5*determinant(sigdelq)$modulus[1]-
         0.5*tri(sigdel,sigdelq)+0.5*d*(k-1)
list(lbadj=lbadj)
}

#------------------------------------------------#
# Functions to calculate variational lower bound #
#------------------------------------------------#

vlb <- function(y,X,W,V,U,sigbeta,alpa,lama,alpb,lamb,alp,lam,sigdel,
                mubetaq,sigbetaq,muaq,sigaq,alpaq,lamaq,mubq,sigbq,
                alpbq,lambq,alpq,lamq,mudelq,qp,vni,epm) {
n <- length(vni)
nsum <- sum(vni)
g <- dim(epm)[2]
p <- dim(X)[2]
d <- dim(U)[2]
s1 <- dim(W)[2]
s2 <- dim(V)[2]
k <- length(alpa)

t1 <- matrix(0,g,k)
for (l in 1:g){
for (j in 1:k){
t1[l,j] <- sum(epm[,l]*qp[,j]) }}

lb <- ( 0.5*(k*(p+s2)+n*s1-nsum*log(2*pi))
 +sum(-(0.5*s2+alpb)*log(lambq)+alpb*log(lamb)-lgamma(alpb)
      +lgamma(alpbq)+alpbq-lamb*alpbq/lambq)
 +sum(-(0.5*s1*colSums(qp)+alpa)*log(lamaq)+alpa*log(lama)
      -lgamma(alpa)+lgamma(alpaq)+alpaq-lama*alpaq/lamaq
      +digamma(alpaq)*(alpa-alpaq+0.5*s1*colSums(qp)))  
 +sum(alp*log(lam)-lgamma(alp)+lgamma(alpq)+alpq-lam*alpq/lamq
     +digamma(alpq)*(alp-alpq+0.5*t1)-(0.5*t1+alp)*log(lamq)) )

for (j in 1:k) {
lb <- ( lb+0.5*determinant(solve(sigbeta[j,,],sigbetaq[j,,]))$modulus[1] 
        -0.5*quadinv(mubetaq[,j],sigbeta[j,,])-0.5*tri(sigbeta[j,,],sigbetaq[j,,])
        +0.5*determinant(as.matrix(sigbq[j,,]))$modulus[1]
        -0.5*alpbq[j]/lambq[j]*(crossprod(mubq[,j])+ tr(sigbq[j,,])) ) }

for (i in 1:n) {
lb <- lb + 0.5*determinant(as.matrix(sigaq[i,,]))$modulus[1] }

lt <- matrix(nrow=n,ncol=k)
for (i in 1:n) {
 for (j in 1:k) {
  yi <- y[(sum(vni[1:i-1])+1):(sum(vni[1:i]))]
  Xi <- as.matrix(X[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
  Wi <- as.matrix(W[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
  Vi <- as.matrix(V[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
  t2 <- yi-Xi%*%mubetaq[,j]-Wi%*%muaq[i,]-Vi%*%mubq[,j]
  t3 <- rep(alpq[,j]/lamq[,j],epm[i,])
lt[i,j] <- ( -0.5*alpaq[j]/lamaq[j]*(crossprod(muaq[i,])+ tr(as.matrix(sigaq[i,,])))
            -0.5*(crossprod(t3*t2,t2) + sum(sigbetaq[j,,]*crossprod(t3*Xi,Xi))
   + sum(sigaq[i,,]*crossprod(t3*Wi,Wi))+ sum(sigbq[j,,]*crossprod(t3*Vi,Vi)) ) )}}         
lb <- lb + sum(qp*lt)

if (k >1) {
lb <- lb-0.5*determinant(sigdel)$modulus[1]-0.5*(k-1)*d*log(2*pi)-0.5*quadinv(mudelq,sigdel)

pp <- matrix(0,n,k)
ppnew <- matrix(0,n,k)
for (j in 2:k) { pp[,j] <- (U%*%mudelq[((j-2)*d+1):((j-1)*d)]) }
for (j in 1:k) {ppnew[,j] <- 1/rowSums(exp(pp-pp[,j]))}
pp <- ppnew
A <- log(pp/qp)
A[is.na(A)] <- 0
A[is.infinite(A)] <- 0
lb <- lb+sum(qp*A)}
list(lb=lb) }


vlbshort <- function(y,X,W,V,U,sigbeta,alpa,lama,alpb,lamb,alp,lam,sigdel,
                mubetaq,sigbetaq,muaq,sigaq,alpaq,lamaq,mubq,sigbq,
                alpbq,lambq,alpq,lamq,mudelq,qp,vni,epm) {
n <- length(vni)
nsum <- sum(vni)
g <- dim(epm)[2]
p <- dim(X)[2]
d <- dim(U)[2]
s1 <- dim(W)[2]
s2 <- dim(V)[2]
k <- length(alpa)

lb <- ( 0.5*(k*(p+s2)+n*s1-nsum*log(2*pi))
 +sum(-alpbq*log(lambq)+alpb*log(lamb)-lgamma(alpb)+lgamma(alpbq))
 +sum(-alpaq*log(lamaq)+alpa*log(lama)-lgamma(alpa)+lgamma(alpaq))  
 +sum(alp*log(lam)-lgamma(alp)+lgamma(alpq)-alpq*log(lamq)) )

for (j in 1:k) {
lb <- ( lb+0.5*determinant(solve(sigbeta[j,,],sigbetaq[j,,]))$modulus[1] 
        -0.5*quadinv(mubetaq[,j],sigbeta[j,,])-0.5*tri(sigbeta[j,,],sigbetaq[j,,])
        +0.5*determinant(as.matrix(sigbq[j,,]))$modulus[1] ) }

for (i in 1:n) {
lb <- lb + 0.5*determinant(as.matrix(sigaq[i,,]))$modulus[1] }

if (k >1) {
lb <- lb-0.5*determinant(sigdel)$modulus[1]-0.5*(k-1)*d*log(2*pi)-0.5*quadinv(mudelq,sigdel)

pp <- matrix(0,n,k)
ppnew <- matrix(0,n,k)
for (j in 2:k) { pp[,j] <- (U%*%mudelq[((j-2)*d+1):((j-1)*d)]) }
for (j in 1:k) {ppnew[,j] <- 1/rowSums(exp(pp-pp[,j]))}
pp <- ppnew
A <- log(pp/qp)
A[is.na(A)] <- 0
A[is.infinite(A)] <- 0
lb <- lb+sum(qp*A)}
list(lb=lb) }


vlbfix <- function(n,p,s2,s1,g,nsum,sigbeta,alpa,lama,alpb,lamb,alp,lam,mubetaq,
                 sigbetaq,alpaq,lamaq,mubq,sigbq,alpbq,lambq,alpq,lamq,qp,set) {
k <- length(alpa)
fixed <- seq(1,k,1)[-set]
alpbf <- alpb[fixed]
lambf <- lamb[fixed]
alpaf <- alpa[fixed]
lamaf <- lama[fixed]
alpf <- matrix(alp[,fixed],g,length(fixed))
lamf <- matrix(lam[,fixed],g,length(fixed))
alpbqf <- alpbq[fixed]
lambqf <- lambq[fixed]
alpaqf <- alpaq[fixed]
lamaqf <- lamaq[fixed]
alpqf <- matrix(alpq[,fixed],g,length(fixed))
lamqf <- matrix(lamq[,fixed],g,length(fixed))
qpf <- matrix(qp[,fixed],n,length(fixed))

lb <- 0.5*(k*(p+s2)+n*s1-nsum*log(2*pi))+
      sum(-alpbqf*log(lambqf)+alpbf*log(lambf)-lgamma(alpbf)
           +lgamma(alpbqf)+alpbqf-lambf*alpbqf/lambqf)+
      sum(-alpaqf*log(lamaqf)+alpaf*log(lamaf)-lgamma(alpaf)
           +lgamma(alpaqf)+alpaqf-lamaf*alpaqf/lamaqf)+
     sum(alpf*log(lamf)-lgamma(alpf)+lgamma(alpqf)+alpqf-lamf*alpqf/lamqf
           -alpqf*log(lamqf)) 

for (j in fixed) {
lb <- ( lb+0.5*determinant(solve(sigbeta[j,,],sigbetaq[j,,]))$modulus[1] 
        -0.5*quadinv(mubetaq[,j],sigbeta[j,,])-0.5*tri(sigbeta[j,,],sigbetaq[j,,])
        +0.5*determinant(as.matrix(sigbq[j,,]))$modulus[1]
        -0.5*alpbq[j]/lambq[j]*(crossprod(mubq[,j])+ tr(sigbq[j,,])) ) }

list(lb=lb) }


vlbadd <- function(y,X,V,U,sigbeta,alpa,lama,alpb,lamb,alp,lam,sigdel,
                 mubetaq,sigbetaq,muaq,sigaq,alpaq,lamaq,mubq,sigbq,
                 alpbq,lambq,alpq,lamq,mudelq,qp,vni,epm,set) {
n <- length(vni)
nsum <- sum(vni)
g <- dim(epm)[2]
p <- dim(X)[2]
d <- dim(U)[2]
s2 <- dim(V)[2]
k <- length(alpa)
alpbc <- alpb[set]
lambc <- lamb[set]
alpac <- alpa[set]
lamac <- lama[set]
alpc <- matrix(alp[,set],g,length(set))
lamc <- matrix(lam[,set],g,length(set))
alpbqc <- alpbq[set]
lambqc <- lambq[set]
alpaqc <- alpaq[set]
lamaqc <- lamaq[set]
alpqc <- matrix(alpq[,set],g,length(set))
lamqc <- matrix(lamq[,set],g,length(set))
qpc <- matrix(qp[,set],n,length(set))
fixed <- seq(1,k,1)[-set]

lb <- sum(-alpbqc*log(lambqc)+alpbc*log(lambc)-lgamma(alpbc)+lgamma(alpbqc))+
      sum(-alpaqc*log(lamaqc)+alpac*log(lamac)-lgamma(alpac)+lgamma(alpaqc))+
      sum(-alpqc*log(lamqc)+alpc*log(lamc)-lgamma(alpc)+lgamma(alpqc)) 

for (j in set) {
lb <- ( lb+0.5*determinant(solve(sigbeta[j,,],sigbetaq[j,,]))$modulus[1] 
        -0.5*quadinv(mubetaq[,j],sigbeta[j,,])-0.5*tri(sigbeta[j,,],sigbetaq[j,,])
        +0.5*determinant(as.matrix(sigbq[j,,]))$modulus[1] )  }

for (i in 1:n) {
lb <- lb + 0.5*determinant(as.matrix(sigaq[i,,]))$modulus[1] }

if (k >1) {
lb <- lb-0.5*determinant(sigdel)$modulus[1]-0.5*(k-1)*d*log(2*pi)-0.5*quadinv(mudelq,sigdel)

pp <- matrix(0,n,k)
ppnew <- matrix(0,n,k)
for (j in 2:k) { pp[,j] <- (U%*%mudelq[((j-2)*d+1):((j-1)*d)]) }
for (j in 1:k) {ppnew[,j] <- 1/rowSums(exp(pp-pp[,j]))}
pp <- ppnew
A <- log(pp/qp)
A[is.na(A)] <- 0
A[is.infinite(A)] <- 0
lb <- lb+sum(qp*A)}

lt <- matrix(0,n,k)
for (i in 1:n) {
 for (j in fixed) {
  yi <- y[(sum(vni[1:i-1])+1):(sum(vni[1:i]))]
  Xi <- as.matrix(X[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
  Wi <- as.matrix(W[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
  Vi <- as.matrix(V[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
  t2 <- yi-Xi%*%mubetaq[,j]-Wi%*%muaq[i,]-Vi%*%mubq[,j]
  t3 <- rep(alpq[,j]/lamq[,j],epm[i,])
lt[i,j] <- ( -0.5*alpaq[j]/lamaq[j]*(crossprod(muaq[i,])+ tr(as.matrix(sigaq[i,,])))
            -0.5*(crossprod(t3*t2,t2) + sum(sigbetaq[j,,]*crossprod(t3*Xi,Xi))
   + sum(sigaq[i,,]*crossprod(t3*Wi,Wi))+ sum(sigbq[j,,]*crossprod(t3*Vi,Vi)) ) )}} 
lb <- lb + sum(qp*lt)

list(lb=lb) }


 #--------------#
 # VA Algorithm #
 #--------------#

MLMM <- function(y,X,W,V,U,sigbeta,alpa,lama,alpb,lamb,alp,lam,sigdel,vni,epm,qp,tol,stfit=NULL){
n <- length(vni)
nsum <- sum(vni)
p <- dim(X)[2]
d <- dim(U)[2]
g <- dim(epm)[2]
s1 <- dim(W)[2]
s2 <- dim(V)[2]
k <- length(alpa)

# Deterministic Updates #
alpbq <- alpb + s2/2

#Initialization
if(is.null(stfit)) {
mubetaq <- matrix(0,p,k)
sigbetaq <- array(0,dim=c(k,p,p)) 
muaq <- matrix(0,n,s1)
sigaq <- array(0,dim=c(n,s1,s1)) 
mubq <- matrix(0,s2,k) 
sigbq <- array(0,dim=c(k,s2,s2)) 
lambq <- alpbq
alpaq <- rep(1,k)
lamaq <- rep(1,k)
alpq <- matrix(1,g,k)
lamq <- matrix(1,g,k)
if (k>1) {mudelq <- rep(0,(k-1)*d)} else {mudelq <- 0}
lbold <- -1000000
} else {
mubetaq <- stfit$mubetaq
sigbetaq <- stfit$sigbetaq
muaq <- stfit$muaq
sigaq <- stfit$sigaq
mubq <- stfit$mubq
sigbq <- stfit$sigbq
alpaq <- stfit$alpaq
lamaq <- stfit$lamaq
lambq <- stfit$lambq
alpq <- stfit$alpq
lamq <- stfit$lamq
mudelq <- stfit$mudelq
qp <- stfit$qp
lbold <- stfit$lb}

dif <- 1000
while (dif > tol) {

# update sigbetaq and mubetaq #
ymb <- matrix(y-crossprod(t(V),mubq),ncol=k)
for (j in 1:k) {
t4 <- matrix(0,p,p)
t5 <- matrix(0,p,1)
for (i in 1:n) { 
  ymbi <- ymb[(sum(vni[1:i-1])+1):(sum(vni[1:i])),j]
  Xi <- as.matrix(X[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
  Wi <- as.matrix(W[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
  t3 <- rep(qp[i,j]*alpq[,j]/lamq[,j],epm[i,])
t4 <- t4 + crossprod(t3*Xi,Xi)
t5 <- t5 + crossprod(t3*Xi,(ymbi-Wi%*%muaq[i,])) }
sigbetaq[j,,] <- solve(solve(sigbeta[j,,])+ t4)
mubetaq[,j] <- crossprod(sigbetaq[j,,],t5)  }

# update sigaq, muaq #
for (i in 1:n) {
 yi <- as.matrix(y[(sum(vni[1:i-1])+1):(sum(vni[1:i]))])
 Xi <- as.matrix(X[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
 Wi <- as.matrix(W[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
 Vi <- as.matrix(V[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
 t1 <- rep(0,vni[i])
 t2 <- matrix(0,s1,1)
for (j in 1:k) {
t3 <- rep(alpq[,j]/lamq[,j],epm[i,])
t1 <- t1 + qp[i,j]*t3
t2 <- t2 + qp[i,j]*crossprod(t3*Wi,(yi-Xi%*%mubetaq[,j]-Vi%*%mubq[,j])) }
sigaq[i,,] <- solve(sum(alpaq/lamaq*qp[i,])*diag(s1) + crossprod(t1*Wi,Wi))
muaq[i,] <- crossprod(sigaq[i,,],t2)  }

# update sigbq, mubq #
for (j in 1:k) {
t7 <- matrix(0,nrow=s2,ncol=s2)
t6 <- matrix(0,nrow=s2,ncol=1)
for (i in 1:n) {
 yi <- as.matrix(y[(sum(vni[1:i-1])+1):(sum(vni[1:i]))])
 Xi <- as.matrix(X[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
 Wi <- as.matrix(W[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
 Vi <- as.matrix(V[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
 t3 <- rep(alpq[,j]/lamq[,j],epm[i,])
t7 <- t7 + qp[i,j]*crossprod(t3*Vi,Vi)
t6 <- t6 + qp[i,j]*crossprod(t3*Vi,yi-Xi%*%mubetaq[,j]-Wi%*%muaq[i,])} 
sigbq[j,,] <- solve(alpbq[j]/lambq[j]*diag(s2)+t7)
mubq[,j] <- crossprod(sigbq[j,,],t6) }

#update mudelq
if (k>1) 
{L1 <- Lmudelq(mudelq,qp,sigdel,U)$L
mudelqnew <- as.vector(t(coef(multinom(qp~U-1,trace=FALSE,decay=1/sigdel[1,1])))) 
L2 <- Lmudelq(mudelqnew,qp,sigdel,U)$L 
if (L2 > L1) {mudelq <- mudelqnew} } 

#update qij for each i,j
if (k>1) {
pp <- matrix(0,n,k)
ppnew <- matrix(0,n,k)
for (j in 2:k) { pp[,j] <- (U%*%mudelq[((j-2)*d+1):((j-1)*d)]) }
for (j in 1:k) {ppnew[,j] <- 1/rowSums(exp(pp-pp[,j]))}
pp <- ppnew

C <- matrix(0,n,k)
for (i in 1:n) {
 for (j in 1:k) {
yi <- y[(sum(vni[1:i-1])+1):(sum(vni[1:i]))]
Xi <- as.matrix(X[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
Wi <- as.matrix(W[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
Vi <- as.matrix(V[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
t2 <- yi-Xi%*%mubetaq[,j]-Wi%*%muaq[i,]-Vi%*%mubq[,j]
t3 <- rep(alpq[,j]/lamq[,j],epm[i,])
C[i,j] <- (-0.5*s1*(log(lamaq[j])-digamma(alpaq[j]))
           -0.5*alpaq[j]/lamaq[j]*(crossprod(muaq[i,])+ tr(as.matrix(sigaq[i,,])))
           -0.5*sum((log(lamq[,j])- digamma(alpq[,j]))*epm[i,])
           -0.5*(crossprod(t3*t2,t2)+ sum(sigbetaq[j,,]*crossprod(t3*Xi,Xi))
+ sum(sigaq[i,,]*crossprod(t3*Wi,Wi))+sum(sigbq[j,,]*crossprod(t3*Vi,Vi)) ) ) 
}}
qp <- pp*exp(C)/matrix(rep(rowSums(pp*exp(C)),times=k),ncol=k) }

# update alpaq, lamaq #
alpaq <- alpa +0.5*s1*colSums(qp)
for (j in 1:k) {
t8 <- 0
for (i in 1:n) {
t8 <- t8 + qp[i,j]*(crossprod(muaq[i,])+ tr(as.matrix(sigaq[i,,]))) }
lamaq[j] <- lama[j]+0.5*t8  }

# update lambq #
for (j in 1:k) {
lambq[j] <- lamb[j]+0.5*(crossprod(mubq[,j])+ tr(sigbq[j,,])) }

# update alpq, lamq #
t1 <- matrix(0,g,k)
for (l in 1:g){
for (j in 1:k){
t1[l,j] <- sum(0.5*epm[,l]*qp[,j]) }}
alpq <- alp+t1

t0 <- matrix(0,nrow=g,ncol=k)
for (i in 1:n) {
yi <- y[(sum(vni[1:i-1])+1):(sum(vni[1:i]))]
Xi <- as.matrix(X[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
Vi <- as.matrix(V[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
Wi <- as.matrix(W[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
t1 <- (matrix(rep(yi-Wi%*%muaq[i,],k),ncol=k)-Xi%*%mubetaq-Vi%*%mubq)^2
t2 <- matrix(0,nrow=vni[i],ncol=k)
for (j in 1:k) {t2[,j] <- diag(Xi%*%sigbetaq[j,,]%*%t(Xi)+Vi%*%sigbq[j,,]%*%t(Vi))}
t1 <- t1+matrix(rep(diag(Wi%*%sigaq[i,,]%*%t(Wi)),k),ncol=k)+t2
ag <- aggregate(t1,list(rep(seq(1,g,1),epm[i,])),sum)
if (k==1) {
ag <- ag[,2]
if (any(epm[i,]==0)) {for (ff in which(epm[i,]==0)) {ag <- insert(ag,ff,0)}}
} else {
ag <- as.matrix(ag[,2:(k+1)])
if (any(epm[i,]==0)) {for (ff in which(epm[i,]==0)) {ag <- insertRow(ag,ff,rep(0,k))}}}
t0 <-t0 +0.5*matrix(rep(qp[i,],g),ncol=k,byrow=TRUE)*ag }
lamq <- lam+t0

lb <- vlbshort(y,X,W,V,U,sigbeta,alpa,lama,alpb,lamb,alp,lam,sigdel,mubetaq,sigbetaq,
       muaq,sigaq,alpaq,lamaq,mubq,sigbq,alpbq,lambq,alpq,lamq,mudelq,qp,vni,epm)$lb     

dif <- (lb-lbold)/abs(lb)
lbold <- lb
cat(lb,k,round(colSums(qp),1),dif,"\n") }

if (k>1) {
gcc <- multinom(qp ~ U-1, decay=1/sigdel[1,1],Hess=TRUE,trace=FALSE)
sigdelq <- solve(gcc$Hess)
lbadj <- lb + vlbadj(U,alpa,sigdel,sigdelq)$lbadj 
} else {
sigdelq=NULL
lbadj=lb}

cat(k,lb,lbadj,"\n")

list (mubetaq=mubetaq,sigbetaq=sigbetaq,muaq=muaq,sigaq=sigaq,mubq=mubq,
sigbq=sigbq,alpaq=alpaq,lamaq=lamaq,alpbq=alpbq,lambq=lambq,
alpq=alpq,lamq=lamq,qp=qp,mudelq=mudelq,lb=lb,lbadj=lbadj,sigdelq=sigdelq)
}


#-----------------------#
# VA Alg for greedy alg #
#-----------------------#

MLMMp <- function(y,X,W,V,U,sigbeta,alpa,lama,alpb,lamb,alp,lam,sigdel,
                 vni,epm,tol,fitpre,gp,cp,m){
n <- length(vni)
nsum <- sum(vni)
p <- dim(X)[2]
d <- dim(U)[2]
g <- dim(epm)[2]
s1 <- dim(W)[2]
s2 <- dim(V)[2]
k <- length(alpa)

# Deterministic Updates #
alpbq <- alpb + s2/2

#Initialization
muaq <- fitpre$muaq
sigaq <- fitpre$sigaq
mubetaq <- cbind(fitpre$mubetaq,fitpre$mubetaq[,cp])
mubq <- cbind(fitpre$mubq,fitpre$mubq[,cp])
lambq <- c(fitpre$lambq,fitpre$lambq[cp])
alpaq <- c(fitpre$alpaq,fitpre$alpaq[cp])
lamaq <- c(fitpre$lamaq,fitpre$lamaq[cp])
alpq <- cbind(fitpre$alpq,fitpre$alpq[,cp])
lamq <- cbind(fitpre$lamq,fitpre$lamq[,cp])
sigbetaq <- array(0,dim=c(k,p,p))
sigbq <- array(0,dim=c(k,s2,s2))
for (j in 1:(k-1)){
sigbetaq[j,,] <- fitpre$sigbetaq[j,,]
sigbq[j,,] <- fitpre$sigbq[j,,] }
sigbetaq[k,,] <- fitpre$sigbetaq[cp,,]
sigbq[k,,] <- fitpre$sigbq[cp,,]
mudelq <- rep(0,(k-1)*d)
qp <- cbind(fitpre$qp,rep(0,n))
for (i in 1:n) {
if (gp[cp,i,m]==2) {
qp[i,k] <- qp[i,cp] 
qp[i,cp] <- 0}}

dif <- 10
lbold <- -1000000
count <- 0

while (dif > tol) {
count <- count+1

# update sigbetaq and mubetaq #
for (j in c(cp,k)) {
t4 <- matrix(0,nrow=p,ncol=p)
t5 <- matrix(0,nrow=p,ncol=1)
for (i in 1:n) { 
  yi <- as.matrix(y[(sum(vni[1:i-1])+1):(sum(vni[1:i]))])
  Xi <- as.matrix(X[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
  Wi <- as.matrix(W[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
  Vi <- as.matrix(V[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
  t3 <- rep(alpq[,j]/lamq[,j],epm[i,])
t4 <- t4 + qp[i,j]*crossprod(t3*Xi,Xi)
t5 <- t5 + qp[i,j]*crossprod(t3*Xi,(yi-Wi%*%muaq[i,]-Vi%*%mubq[,j]))}
sigbetaq[j,,] <- solve(solve(sigbeta[j,,])+ t4)
mubetaq[,j] <- crossprod(sigbetaq[j,,],t5)  }

# update sigaq, muaq #
for (i in 1:n) {
 yi <- as.matrix(y[(sum(vni[1:i-1])+1):(sum(vni[1:i]))])
 Xi <- as.matrix(X[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
 Wi <- as.matrix(W[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
 Vi <- as.matrix(V[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
 t1 <- rep(0,vni[i])
 t2 <- matrix(0,nrow=s1,ncol=1)
for (j in 1:k) {
t3 <- rep(alpq[,j]/lamq[,j],epm[i,])
t1 <- t1 + qp[i,j]*t3
t2 <- t2 + qp[i,j]*crossprod(t3*Wi,(yi-Xi%*%mubetaq[,j]-Vi%*%mubq[,j]))}
sigaq[i,,] <- solve(sum(alpaq/lamaq*qp[i,])*diag(s1) + crossprod(t1*Wi,Wi))
muaq[i,] <- crossprod(sigaq[i,,],t2)  }

# update sigbq, mubq #
for (j in c(cp,k)) {
t7 <- matrix(0,nrow=s2,ncol=s2)
t6 <- matrix(0,nrow=s2,ncol=1)
for (i in 1:n) {
 yi <- as.matrix(y[(sum(vni[1:i-1])+1):(sum(vni[1:i]))])
 Xi <- as.matrix(X[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
 Wi <- as.matrix(W[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
 Vi <- as.matrix(V[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
 t3 <- rep(alpq[,j]/lamq[,j],epm[i,])
t7 <- t7 + qp[i,j]*crossprod(t3*Vi,Vi)
t6 <- t6 + qp[i,j]*crossprod(t3*Vi,yi-Xi%*%mubetaq[,j]-Wi%*%muaq[i,])} 
sigbq[j,,] <- solve(alpbq[j]/lambq[j]*diag(s2)+t7)
mubq[,j] <- crossprod(sigbq[j,,],t6) }

#update mudelq
if (k>1) 
{L1 <- Lmudelq(mudelq,qp,sigdel,U)$L
mudelqnew <- as.vector(t(coef(multinom(qp~U-1,trace=FALSE,decay=1/sigdel[1,1])))) 
L2 <- Lmudelq(mudelqnew,qp,sigdel,U)$L 
if (L2 > L1) {mudelq <- mudelqnew} } 

#update qij 
if (k>1) {
pp <- matrix(0,nrow=n,ncol=k)
ppnew <- matrix(0,nrow=n,ncol=k)
for (j in 2:k) { pp[,j] <- (U%*%mudelq[((j-2)*d+1):((j-1)*d)]) }
for (j in 1:k) {ppnew[,j] <- 1/rowSums(exp(pp-pp[,j]))}
pp <- ppnew

C <- matrix(0,n,k)
for (i in 1:n) {
 for (j in c(cp,k)) {
yi <- y[(sum(vni[1:i-1])+1):(sum(vni[1:i]))]
Xi <- as.matrix(X[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
Wi <- as.matrix(W[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
Vi <- as.matrix(V[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
t2 <- yi-Xi%*%mubetaq[,j]-Wi%*%muaq[i,]-Vi%*%mubq[,j]
t3 <- rep(alpq[,j]/lamq[,j],epm[i,])
C[i,j] <- (-0.5*s1*(log(lamaq[j])-digamma(alpaq[j]))
           -0.5*alpaq[j]/lamaq[j]*(crossprod(muaq[i,])+ tr(as.matrix(sigaq[i,,])))
           -0.5*sum((log(lamq[,j])- digamma(alpq[,j]))*epm[i,])
           -0.5*(crossprod(t3*t2,t2)+ sum(sigbetaq[j,,]*crossprod(t3*Xi,Xi))
                   + sum(sigaq[i,,]*crossprod(t3*Wi,Wi))
                   + sum(sigbq[j,,]*crossprod(t3*Vi,Vi)) ) ) }}
Q1 <- rowSums(pp[,c(cp,k)]*exp(C[,c(cp,k)]))
Q2 <- rowSums(as.matrix(qp[,-c(cp,k)],n,(k-2)))
Q3 <- rowSums(as.matrix(qp[,c(cp,k)],n,(k-2)))
qp[,c(cp,k)] <- pp[,c(cp,k)]*exp(C[,c(cp,k)])/matrix(rep(Q1+(Q1*Q2/Q3),times=2),ncol=2)
for (i in 1:n) {if (Q3[i]==0 | Q1[i]==0) {qp[i,c(cp,k)] <- c(0,0)}}}

# update alpaq, lamaq #
alpaq[c(cp,k)] <- alpa[c(cp,k)] +0.5*s1*colSums(qp)[c(cp,k)]
for (j in c(cp,k)) {
t8 <- 0
for (i in 1:n) {
t8 <- t8 + qp[i,j]*(crossprod(muaq[i,])+ tr(as.matrix(sigaq[i,,]))) }
lamaq[j] <- lama[j]+0.5*t8  }

# update lambq #
for (j in c(cp,k)) {
lambq[j] <- lamb[j]+0.5*(crossprod(mubq[,j])+ tr(sigbq[j,,])) }

# update alpq, lamq #
t1 <- matrix(0,nrow=g,ncol=k)
for (l in 1:g){
for (j in c(cp,k)){
t1[l,j] <- sum(0.5*epm[,l]*qp[,j]) }}
alpq[,c(cp,k)] <- alp[,c(cp,k)]+t1[,c(cp,k)]

t0 <- matrix(0,g,2)
for (i in 1:n) {
yi <- y[(sum(vni[1:i-1])+1):(sum(vni[1:i]))]
Xi <- as.matrix(X[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
Wi <- as.matrix(W[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
Vi <- as.matrix(V[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
t1 <- (matrix(rep(yi-Wi%*%muaq[i,],2),ncol=2)-Xi%*%mubetaq[,c(cp,k)]-Vi%*%mubq[,c(cp,k)])^2
t2 <- matrix(0,nrow=vni[i],ncol=2)
num <- 0
for (j in c(cp,k)) {
num <- num+1
t2[,num] <- diag(Vi%*%sigbq[j,,]%*%t(Vi)+Xi%*%sigbetaq[j,,]%*%t(Xi))}
t1 <- t1+matrix(rep(diag(Wi%*%sigaq[i,,]%*%t(Wi)),2),ncol=2)+t2
ag <- aggregate(t1,list(rep(seq(1,g,1),epm[i,])),sum)
ag <- cbind(ag[,2],ag[,3])
if (any(epm[i,]==0)) {for (ff in which(epm[i,]==0)) {ag <- insertRow(ag,ff,rep(0,2))}}
t0 <-t0 +0.5*matrix(rep(qp[i,c(cp,k)],g),ncol=2,byrow=TRUE)*ag }
lamq[,c(cp,k)] <- lam[,c(cp,k)]+t0

if (count==1){lbinit <- vlbfix(n,p,s2,s1,g,nsum,sigbeta,alpa,lama,alpb,lamb,alp,lam,mubetaq,
                 sigbetaq,alpaq,lamaq,mubq,sigbq,alpbq,lambq,alpq,lamq,qp,c(cp,k))$lb}
lbaft <- vlbadd(y,X,V,U,sigbeta,alpa,lama,alpb,lamb,alp,lam,sigdel,
                 mubetaq,sigbetaq,muaq,sigaq,alpaq,lamaq,mubq,sigbq,
                 alpbq,lambq,alpq,lamq,mudelq,qp,vni,epm,c(cp,k))$lb
lb <- lbinit+lbaft

dif <- lb-lbold
lbold <- lb
prop <- colSums(qp)
cat(lb,round(prop[c(cp,k)],1),dif,"\n") }

list (mubetaq=mubetaq,sigbetaq=sigbetaq,muaq=muaq,sigaq=sigaq,
mubq=mubq,sigbq=sigbq,alpaq=alpaq,lamaq=lamaq,alpbq=alpbq,lambq=lambq,alpq=alpq,
lamq=lamq,qp=qp,mudelq=mudelq,lb=lb)}


MLMMpa <- function(y,X,W,V,U,sigbeta,alpa,lama,alpb,lamb,alp,lam,sigdel,
                 vni,epm,qp,alpq,lamq,alpaq,lamaq,alpbq,lambq,
                 mubetaq,sigbetaq,mubq,sigbq,set,tol){
n <- length(vni)
nsum <- sum(vni)
p <- dim(X)[2]
d <- dim(U)[2]
g <- dim(epm)[2]
s1 <- dim(W)[2]
s2 <- dim(V)[2]
k <- length(alpa)

# Deterministic Updates #
alpbq <- alpb + s2/2

#Initialization
muaq <- matrix(0,nrow=n,ncol=s1)
sigaq <- array(0,dim=c(n,s1,s1))
mudelq <- rep(0,(k-1)*d)

dif <- 10
lbold <- -100000
count <- 0

while (dif > tol) {
count <- count+1

# update sigbetaq and mubetaq #
for (j in set) {
t4 <- matrix(0,nrow=p,ncol=p)
t5 <- matrix(0,nrow=p,ncol=1)
for (i in 1:n) { 
  yi <- as.matrix(y[(sum(vni[1:i-1])+1):(sum(vni[1:i]))])
  Xi <- as.matrix(X[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
  Wi <- as.matrix(W[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
  Vi <- as.matrix(V[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
  t3 <- rep(alpq[,j]/lamq[,j],epm[i,])
t4 <- t4 + qp[i,j]*crossprod(t3*Xi,Xi)
t5 <- t5 + qp[i,j]*crossprod(t3*Xi,(yi-Wi%*%muaq[i,]-Vi%*%mubq[,j]))}
sigbetaq[j,,] <- solve(solve(sigbeta[j,,])+ t4)
mubetaq[,j] <- crossprod(sigbetaq[j,,],t5)  }

# update sigaq, muaq #
for (i in 1:n) {
 yi <- as.matrix(y[(sum(vni[1:i-1])+1):(sum(vni[1:i]))])
 Xi <- as.matrix(X[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
 Wi <- as.matrix(W[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
 Vi <- as.matrix(V[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
 t1 <- rep(0,vni[i])
 t2 <- matrix(0,nrow=s1,ncol=1)
for (j in 1:k) {
t3 <- rep(alpq[,j]/lamq[,j],epm[i,])
t1 <- t1 + qp[i,j]*t3
t2 <- t2 + qp[i,j]*crossprod(t3*Wi,(yi-Xi%*%mubetaq[,j]-Vi%*%mubq[,j]))}
sigaq[i,,] <- solve(sum(alpaq/lamaq*qp[i,])*diag(s1) + crossprod(t1*Wi,Wi))
muaq[i,] <- crossprod(sigaq[i,,],t2)  }

# update sigbq, mubq #
for (j in set) {
t7 <- matrix(0,nrow=s2,ncol=s2)
t6 <- matrix(0,nrow=s2,ncol=1)
for (i in 1:n) {
 yi <- as.matrix(y[(sum(vni[1:i-1])+1):(sum(vni[1:i]))])
 Xi <- as.matrix(X[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
 Wi <- as.matrix(W[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
 Vi <- as.matrix(V[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
 t3 <- rep(alpq[,j]/lamq[,j],epm[i,])
t7 <- t7 + qp[i,j]*crossprod(t3*Vi,Vi)
t6 <- t6 + qp[i,j]*crossprod(t3*Vi,yi-Xi%*%mubetaq[,j]-Wi%*%muaq[i,])} 
sigbq[j,,] <- solve(alpbq[j]/lambq[j]*diag(s2)+t7)
mubq[,j] <- crossprod(sigbq[j,,],t6) }

#update mudelq
if (k>1) 
{L1 <- Lmudelq(mudelq,qp,sigdel,U)$L
mudelqnew <- as.vector(t(coef(multinom(qp~U-1,trace=FALSE,decay=1/sigdel[1,1])))) 
L2 <- Lmudelq(mudelqnew,qp,sigdel,U)$L 
if (L2 > L1) {mudelq <- mudelqnew} } 

#update qij 
if (k>1) {
pp <- matrix(0,nrow=n,ncol=k)
ppnew <- matrix(0,nrow=n,ncol=k)
for (j in 2:k) { pp[,j] <- (U%*%mudelq[((j-2)*d+1):((j-1)*d)]) }
for (j in 1:k) {ppnew[,j] <- 1/rowSums(exp(pp-pp[,j]))}
pp <- ppnew

C <- matrix(0,nrow=n,ncol=k)
for (i in 1:n) {
 for (j in set) {
yi <- y[(sum(vni[1:i-1])+1):(sum(vni[1:i]))]
Xi <- as.matrix(X[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
Wi <- as.matrix(W[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
Vi <- as.matrix(V[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
t2 <- yi-Xi%*%mubetaq[,j]-Wi%*%muaq[i,]-Vi%*%mubq[,j]
t3 <- rep(alpq[,j]/lamq[,j],epm[i,])
C[i,j] <- (-0.5*s1*(log(lamaq[j])-digamma(alpaq[j]))
           -0.5*alpaq[j]/lamaq[j]*(crossprod(muaq[i,])+ tr(as.matrix(sigaq[i,,])))
           -0.5*sum((log(lamq[,j])- digamma(alpq[,j]))*epm[i,])
           -0.5*(crossprod(t3*t2,t2)+ sum(sigbetaq[j,,]*crossprod(t3*Xi,Xi))
                   + sum(sigaq[i,,]*crossprod(t3*Wi,Wi))
                   + sum(sigbq[j,,]*crossprod(t3*Vi,Vi)) ) ) }}
Q1 <- rowSums(pp[,set]*exp(C[,set]))
Q2 <- rowSums(as.matrix(qp[,-set],nrow=n,ncol=(k-2)))
Q3 <- rowSums(as.matrix(qp[,set],nrow=n,ncol=(k-2)))
qp[,set] <- pp[,set]*exp(C[,set])/matrix(rep(Q1+(Q1*Q2/Q3),times=length(set)),ncol=length(set))
for (i in 1:n) {if (Q3[i]==0 | Q1[i]==0) {qp[i,set] <- c(0,0)}}}

# update alpaq, lamaq #
alpaq[set] <- alpa[set] +0.5*s1*colSums(qp)[set]
for (j in set) {
t8 <- 0
for (i in 1:n) {
t8 <- t8 + qp[i,j]*(crossprod(muaq[i,])+ tr(as.matrix(sigaq[i,,]))) }
lamaq[j] <- lama[j]+0.5*t8  }

# update lambq #
for (j in set) {
lambq[j] <- lamb[j]+0.5*(crossprod(mubq[,j])+ tr(sigbq[j,,])) }

# update alpq, lamq #
t1 <- matrix(0,nrow=g,ncol=k)
for (l in 1:g){
for (j in set){
t1[l,j] <- sum(0.5*epm[,l]*qp[,j]) }}
alpq[,set] <- alp[,set]+t1[,set]

sel <- length(set)
t0 <- matrix(0,nrow=g,ncol=sel)
for (i in 1:n) {
yi <- y[(sum(vni[1:i-1])+1):(sum(vni[1:i]))]
Xi <- as.matrix(X[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
Vi <- as.matrix(V[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
Wi <- as.matrix(W[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
t1 <- (matrix(rep(yi-Wi%*%muaq[i,],sel),ncol=sel)-Xi%*%mubetaq[,set]-Vi%*%mubq[,set])^2
t2 <- matrix(0,nrow=vni[i],ncol=sel)
num <- 0
for (j in set) {
num <- num+1
t2[,num] <- diag(Vi%*%sigbq[j,,]%*%t(Vi)+Xi%*%sigbetaq[j,,]%*%t(Xi))}
t1 <- t1+matrix(rep(diag(Wi%*%sigaq[i,,]%*%t(Wi)),sel),ncol=sel)+t2
ag <- aggregate(t1,list(rep(seq(1,g,1),epm[i,])),sum)
ag <- as.matrix(ag[,seq(2,by=1,length.out=sel)])
if (any(epm[i,]==0)) {for (ff in which(epm[i,]==0)) {ag <- insertRow(ag,ff,rep(0,sel))}}
t0 <-t0 +0.5*matrix(rep(qp[i,set],g),ncol=sel,byrow=TRUE)*ag }
lamq[,set] <- lam[,set]+t0

if (count==1) {lbinit <- vlbfix(n,p,s2,s1,g,nsum,sigbeta,alpa,lama,alpb,lamb,alp,lam,mubetaq,
                 sigbetaq,alpaq,lamaq,mubq,sigbq,alpbq,lambq,alpq,lamq,qp,set)$lb}

lbaft <- vlbadd(y,X,V,U,sigbeta,alpa,lama,alpb,lamb,alp,lam,sigdel,
                 mubetaq,sigbetaq,muaq,sigaq,alpaq,lamaq,mubq,sigbq,
                 alpbq,lambq,alpq,lamq,mudelq,qp,vni,epm,set)$lb
lb <- lbinit+lbaft
   
dif <- (lb-lbold)/abs(lb)
lbold <- lb

cat(lb,round(colSums(qp),1),dif,"\n") }

if (k>1) {
gcc <- multinom(qp ~ U-1, decay=1/sigdel[1,1],Hess=TRUE,trace=FALSE)
sigdelq <- solve(gcc$Hess)
lbadj <- lb + vlbadj(U,alpa,sigdel,sigdelq)$lbadj 
} else {
sigdelq=NULL
lbadj=lb}

cat(k,lb,lbadj,"\n")

list (mubetaq=mubetaq,sigbetaq=sigbetaq,muaq=muaq,sigaq=sigaq,mubq=mubq,
sigbq=sigbq,alpaq=alpaq,lamaq=lamaq,alpbq=alpbq,lambq=lambq,alpq=alpq,
lamq=lamq,qp=qp,mudelq=mudelq,lb=lb,sigdelq=sigdelq,lbadj=lbadj)
}


#--------#
# Greedy #
#--------#

greedy <- function(y,X,W,V,U,vni,epm,M,IGa,IGb,sigdelp,sigbetap) {

bef <- proc.time() 

# Fit 1 comp model #
k <-1
n <- length(vni)
nsum <- sum(vni)
p <- dim(X)[2]
s1 <- dim(W)[2]
s2 <- dim(V)[2]
d <- dim(U)[2]
g <- dim(epm)[2]

# Priors #
sigbeta <- array(0,dim=c(k,p,p))
for (j in 1:k) {sigbeta[j,,] <- sigbetap*diag(p)}
alpa <- rep(IGa,k)
lama <- rep(IGb,k)
alpb <- rep(IGa,k)
lamb <- rep(IGb,k)
alp <- matrix(IGa,g,k)
lam <- matrix(IGb,g,k)
sigdel <- sigdelp*diag((k-1)*d)
qp <- matrix(1,n,1)

fitpre <- MLMM(y,X,W,V,U,sigbeta,alpa,lama,alpb,lamb,alp,
                lam,sigdel,vni,epm,qp,1.0e-5)

dsp <- NULL       # dsp (do not split) is the set of comp we do not attempt to split
dif <- 10         # set initial value of dif
tot <- seq(1,n,1) # total no. of data points

while (dif > 1) {

k <- length(fitpre$alpaq) # count current no. of mixture components

# Randomly partition each components into 2 subcomponents #
postp <- rep(0,n)
for(i in 1:n){
 for(j in 1:k){
if (fitpre$qp[i,j] == max(fitpre$qp[i,])) (postp[i]<-j)}} 
gp <- array(0,dim=c(k,n,M))
for (j in 1:k) {
for (m in 1:M) {
current <- tot[postp==j]
len <- length(current)
sam <- sample(c(1,2),len,replace=TRUE,prob=c(0.5,0.5))
gp[j,which(postp==j),m] <- sam }}


(K <- dim(gp)[1])
fitore <- fitpre
fits <- list(NULL)
k <- K+1 # no. of mixture components increased by 1

# Priors # 
sigbeta <- array(0,dim=c(k,p,p))
for (j in 1:k) {sigbeta[j,,] <- sigbetap*diag(p)}
alpa <- rep(IGa,k)
lama <- rep(IGb,k)
alpb <- rep(IGa,k)
lamb <- rep(IGb,k)
alp <- matrix(IGa,nrow=g,ncol=k)
lam <- matrix(IGb,nrow=g,ncol=k)
sigdel <- sigdelp*diag((k-1)*d)

if (length(dsp)>=1) {sp <- seq(1,K,1)[-dsp]} else {sp <- seq(1,K,1)} 
lbds <- cbind(sp,rep(0,length(sp)))  # sp is the set of comps we attempt splitting
num <- 0

for (cp in sp) {
num <- num+1
fitold <- -10000000
nredz <- 0

for (m in 1:M) {
fitcan <- MLMMp(y,X,W,V,U,sigbeta,alpa,lama,alpb,lamb,alp,
                lam,sigdel,vni,epm,1,fitpre,gp,cp,m) 
cat(cp,m,fitcan$lb,"\n")

if (fitcan$lb > fitold) {
fits[[cp]] <- fitcan
lbds[num,2] <- fitcan$lb
fitold <- fitcan$lb} }

if (any( colSums(fits[[cp]]$qp)<=1)) {dsp <- c(dsp,cp)}
cat(cp,dsp,"\n")
}

lbdsord <- sort(lbds[,2],decreasing=TRUE)
ap <- 1
for (jj in 1:length(sp)) {

jm <- which(lbds[,2]==lbdsord[jj])
cp <- as.numeric(lbds[jm,1])
rp <- which(fitpre$lamaq==fitore$lamaq[cp])

if  (all(colSums(fits[[cp]]$qp)>=1)) {

# Initialization #
qp <- fitpre$qp[,-rp]
qp <- cbind(qp,fits[[cp]]$qp[,c(cp,(K+1))])

if (k==2) {
alpq <- fits[[cp]]$alpq
lamq <- fits[[cp]]$lamq
} else {
alpq <- matrix(fitpre$alpq[,-rp],ncol=length(fitpre$alpaq)-1)
alpq <- cbind(alpq,matrix(fits[[cp]]$alpq[,c(cp,(K+1))],ncol=2))
lamq <- matrix(fitpre$lamq[,-rp],ncol=length(fitpre$alpaq)-1)
lamq <- cbind(lamq,matrix(fits[[cp]]$lamq[,c(cp,(K+1))],ncol=2)) }

alpaq <- fitpre$alpaq[-rp]
alpaq <- c(alpaq,fits[[cp]]$alpaq[c(cp,(K+1))])
lamaq <- fitpre$lamaq[-rp]
lamaq <- c(lamaq,fits[[cp]]$lamaq[c(cp,(K+1))])
alpbq <- fitpre$alpbq[-rp]
alpbq <- c(alpbq,fits[[cp]]$alpbq[c(cp,(K+1))])
lambq <- fitpre$lambq[-rp]
lambq <- c(lambq,fits[[cp]]$lambq[c(cp,(K+1))])
mubetaq <- fitpre$mubetaq[,-rp]
mubetaq <- cbind(mubetaq,fits[[cp]]$mubetaq[,c(cp,(K+1))])
mubq <- fitpre$mubq[,-rp]
mubq <- cbind(mubq,fits[[cp]]$mubq[,c(cp,(K+1))])
muetaq <- fits[[cp]]$muetaq
sigetaq <- fits[[cp]]$sigetaq
k <- length(fitpre$alpaq)+1

if (k==2) {
sigbetaq <- fits[[cp]]$sigbetaq
sigbq <- fits[[cp]]$sigbq
} else {
sigbetaq <- array(0,dim=c(k,p,p))
sigbetaq[1:(k-2),,]<- fitpre$sigbetaq[-cp,,]
sigbetaq[(k-1):k,,] <- fits[[cp]]$sigbetaq[c(cp,(K+1)),,]
sigbq <- array(0,dim=c(k,s2,s2))
sigbq[1:(k-2),,]<- fitpre$sigbq[-cp,,]
sigbq[(k-1):k,,] <- fits[[cp]]$sigbq[c(cp,(K+1)),,] }

# Priors #
sigbeta <- array(0,dim=c(k,p,p))
for (j in 1:k) {sigbeta[j,,] <- sigbetap*diag(p)}
alpa <- rep(IGa,k)
lama <- rep(IGb,k)
alpb <- rep(IGa,k)
lamb <- rep(IGb,k)
alp <- matrix(IGa,g,k)
lam <- matrix(IGb,g,k)
sigdel <- sigdelp*diag((k-1)*d)

set <- seq((k-2*ap+1),k,1) # set of components to update
fitprenew <- MLMMpa(y,X,W,V,U,sigbeta,alpa,lama,alpb,lamb,alp,
                lam,sigdel,vni,epm,qp,alpq,lamq,alpaq,lamaq,
                alpbq,lambq,mubetaq,sigbetaq,mubq,sigbq,set,1.0e-5)
if (fitprenew$lbadj>fitpre$lbadj) {
fitpre <- fitprenew
ap <- ap + 1} else {break} }} 

# Priors #
k <- length(fitpre$alpaq)
sigbeta <- array(0,dim=c(k,p,p))
for (j in 1:k) {sigbeta[j,,] <- sigbetap*diag(p)}
alpa <- rep(IGa,k)
lama <- rep(IGb,k)
alpb <- rep(IGa,k)
lamb <- rep(IGb,k)
alp <- matrix(IGa,g,k)
lam <- matrix(IGb,g,k)
sigdel <- sigdelp*diag((k-1)*d)

fitfinal <- MLMM(y,X,W,V,U,sigbeta,alpa,lama,alpb,lamb,alp,
                lam,sigdel,vni,epm,qp,1.0e-5,fitpre)

dsptemp <- dsp
dsp <- NULL
for (dd in dsptemp){
dg <- which(fitpre$lamaq==fitore$lamaq[dd])
if (max(abs(fitpre$qp[,dg]-fitfinal$qp[,dg])) <= 0.5){
dsp <- c(dsp,dg)} }

if ( any(colSums(fitfinal$qp)<=1) ){dsp <- c(dsp,which(colSums(fitfinal$qp)<=1))}

fitpre <- fitfinal
(lbpre <- fitore$lbadj)
(lbnext <- fitfinal$lbadj)
(dif <- lbnext-lbpre)
cat(dim(fitpre$mubetaq)[2],fitpre$lbadj,dsp,dsptemp,"\n") }

aft <- proc.time()
dur <- aft-bef
list(dur=dur,fitpre=fitpre)}


MLMMmerge <- function(y,X,W,V,U,IGa,IGb,sigbetap,sigdelp,vni,epm,tol,fit,m1,m2,type){
n <- length(vni)
nsum <- sum(vni)
p <- dim(X)[2]
d <- dim(U)[2]
g <- dim(epm)[2]
s1 <- dim(W)[2]
s2 <- dim(V)[2]
k <- dim(fit$mubetaq)[2]-1

sigbeta <- array(0,dim=c(k,p,p))
for (j in 1:k) {sigbeta[j,,] <- sigbetap*diag(p)}
alpa <- rep(IGa,k)
lama <- rep(IGb,k)
alpb <- rep(IGa,k)
lamb <- rep(IGb,k)
alp <- matrix(IGa,g,k)
lam <- matrix(IGb,g,k)
sigdel <- sigdelp*diag((k-1)*d)

# Deterministic Updates #
alpbq <- alpb + s2/2

#Initialization
mubetaq <- fit$mubetaq[,-m2]
sigbetaq <- fit$sigbetaq[-m2,,]
muaq <- fit$muaq
sigaq <- fit$sigaq
mubq <- fit$mubq[,-m2]
sigbq <- fit$sigbq[-m2,,]
alpaq <- fit$alpaq[-m2]
lamaq <- fit$lamaq[-m2]
alpbq <- fit$alpbq[-m2]
lambq <- fit$lambq[-m2]
alpq <- matrix((fit$alpq[,-m2]),g,k)
lamq <- matrix((fit$lamq[,-m2]),g,k)
qp <- fit$qp
qp <- qp[,-m2]
if (m1<m2) {
if (type=='single') {set <- m1} else {set <- seq(1,k,1)}
qp[,m1] <- fit$qp[,m1]+fit$qp[,m2]
} else  {
if (type=='single') {set <- m1-1} else {set <- seq(1,k,1)}
qp[,(m1-1)] <- fit$qp[,m1]+fit$qp[,m2]}
if (k>1) {mudelq <- rep(0,(k-1)*d)} else {mudelq <- 0}
lbold <- -1000000

count <- 0
dif <- 10
while (dif > tol) {
count <- count+1

# update sigbetaq and mubetaq #
for (j in set) {
t4 <- matrix(0,nrow=p,ncol=p)
t5 <- matrix(0,nrow=p,ncol=1)
for (i in 1:n) { 
  yi <- as.matrix(y[(sum(vni[1:i-1])+1):(sum(vni[1:i]))])
  Xi <- as.matrix(X[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
  Wi <- as.matrix(W[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
  Vi <- as.matrix(V[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
  t3 <- rep(alpq[,j]/lamq[,j],epm[i,])
t4 <- t4 + qp[i,j]*crossprod(t3*Xi,Xi)
t5 <- t5 + qp[i,j]*crossprod(t3*Xi,(yi-Wi%*%muaq[i,]-Vi%*%mubq[,j]))}
sigbetaq[j,,] <- solve(solve(sigbeta[j,,])+ t4)
mubetaq[,j] <- crossprod(sigbetaq[j,,],t5)  }

# update sigaq, muaq #
for (i in 1:n) {
 yi <- as.matrix(y[(sum(vni[1:i-1])+1):(sum(vni[1:i]))])
 Xi <- as.matrix(X[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
 Wi <- as.matrix(W[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
 Vi <- as.matrix(V[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
 t1 <- rep(0,vni[i])
 t2 <- matrix(0,nrow=s1,ncol=1)
for (j in 1:k) {
t3 <- rep(alpq[,j]/lamq[,j],epm[i,])
t1 <- t1 + qp[i,j]*t3
t2 <- t2 + qp[i,j]*crossprod(t3*Wi,(yi-Xi%*%mubetaq[,j]-Vi%*%mubq[,j]))}
sigaq[i,,] <- solve(sum(alpaq/lamaq*qp[i,])*diag(s1) + crossprod(t1*Wi,Wi))
muaq[i,] <- crossprod(sigaq[i,,],t2)  }

# update sigbq, mubq #
for (j in set) {
t7 <- matrix(0,nrow=s2,ncol=s2)
t6 <- matrix(0,nrow=s2,ncol=1)
for (i in 1:n) {
 yi <- as.matrix(y[(sum(vni[1:i-1])+1):(sum(vni[1:i]))])
 Xi <- as.matrix(X[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
 Wi <- as.matrix(W[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
 Vi <- as.matrix(V[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
 t3 <- rep(alpq[,j]/lamq[,j],epm[i,])
t7 <- t7 + qp[i,j]*crossprod(t3*Vi,Vi)
t6 <- t6 + qp[i,j]*crossprod(t3*Vi,yi-Xi%*%mubetaq[,j]-Wi%*%muaq[i,])} 
sigbq[j,,] <- solve(alpbq[j]/lambq[j]*diag(s2)+t7)
mubq[,j] <- crossprod(sigbq[j,,],t6) }

#update mudelq
if (k>1) 
{L1 <- Lmudelq(mudelq,qp,sigdel,U)$L
mudelqnew <- as.vector(t(coef(multinom(qp~U-1,trace=FALSE,decay=1/sigdel[1,1])))) 
L2 <- Lmudelq(mudelqnew,qp,sigdel,U)$L 
if (L2 > L1) {mudelq <- mudelqnew} } 

#update qij 
if (k>1) {
pp <- matrix(0,nrow=n,ncol=k)
ppnew <- matrix(0,nrow=n,ncol=k)
for (j in 2:k) { pp[,j] <- (U%*%mudelq[((j-2)*d+1):((j-1)*d)]) }
for (j in 1:k) {ppnew[,j] <- 1/rowSums(exp(pp-pp[,j]))}
pp <- ppnew

C <- matrix(0,nrow=n,ncol=k)
for (i in 1:n) {
 for (j in set) {
yi <- y[(sum(vni[1:i-1])+1):(sum(vni[1:i]))]
Xi <- as.matrix(X[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
Wi <- as.matrix(W[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
Vi <- as.matrix(V[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
t2 <- yi-Xi%*%mubetaq[,j]-Wi%*%muaq[i,]-Vi%*%mubq[,j]
t3 <- rep(alpq[,j]/lamq[,j],epm[i,])
C[i,j] <- (-0.5*s1*(log(lamaq[j])-digamma(alpaq[j]))
           -0.5*alpaq[j]/lamaq[j]*(crossprod(muaq[i,])+ tr(as.matrix(sigaq[i,,])))
           -0.5*sum((log(lamq[,j])- digamma(alpq[,j]))*epm[i,])
           -0.5*(crossprod(t3*t2,t2)+ sum(sigbetaq[j,,]*crossprod(t3*Xi,Xi))
                   + sum(sigaq[i,,]*crossprod(t3*Wi,Wi))
                   + sum(sigbq[j,,]*crossprod(t3*Vi,Vi)) ) ) }}
Q1 <- rowSums(as.matrix(pp[,set]*exp(C[,set]),nrow=n))
Q2 <- rowSums(as.matrix(qp[,-set],nrow=n,ncol=(k-2)))
Q3 <- rowSums(as.matrix(qp[,set],nrow=n,ncol=(k-2)))
qp[,set] <- pp[,set]*exp(C[,set])/matrix(rep(Q1+(Q1*Q2/Q3),times=length(set)),ncol=length(set))
for (i in 1:n) {if (Q3[i]==0 | Q1[i]==0) {qp[i,set] <- c(0,0)}}}

# update alpaq, lamaq #
alpaq[set] <- alpa[set] +0.5*s1*colSums(qp)[set]
for (j in set) {
t8 <- 0
for (i in 1:n) {
t8 <- t8 + qp[i,j]*(crossprod(muaq[i,])+ tr(as.matrix(sigaq[i,,]))) }
lamaq[j] <- lama[j]+0.5*t8  }

# update lambq #
for (j in set) {
lambq[j] <- lamb[j]+0.5*(crossprod(mubq[,j])+ tr(sigbq[j,,])) }

# update alpq, lamq #
t1 <- matrix(0,nrow=g,ncol=k)
for (l in 1:g){
for (j in set){
t1[l,j] <- sum(0.5*epm[,l]*qp[,j]) }}
alpq[,set] <- alp[,set]+t1[,set]

sel <- length(set)
t0 <- matrix(0,nrow=g,ncol=sel)
for (i in 1:n) {
yi <- y[(sum(vni[1:i-1])+1):(sum(vni[1:i]))]
Xi <- as.matrix(X[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
Vi <- as.matrix(V[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
Wi <- as.matrix(W[(sum(vni[1:i-1])+1):(sum(vni[1:i])),])
t1 <- (matrix(rep(yi-Wi%*%muaq[i,],sel),ncol=sel)-Xi%*%mubetaq[,set]-Vi%*%mubq[,set])^2
t2 <- matrix(0,nrow=vni[i],ncol=sel)
num <- 0
for (j in set) {
num <- num+1
t2[,num] <- diag(Vi%*%sigbq[j,,]%*%t(Vi)+Xi%*%sigbetaq[j,,]%*%t(Xi))}
t1 <- t1+matrix(rep(diag(Wi%*%sigaq[i,,]%*%t(Wi)),sel),ncol=sel)+t2
ag <- aggregate(t1,list(rep(seq(1,g,1),epm[i,])),sum)
ag <- as.matrix(ag[,seq(2,by=1,length.out=sel)])
if (any(epm[i,]==0)) {for (ff in which(epm[i,]==0)) {ag <- insertRow(ag,ff,rep(0,sel))}}
t0 <-t0 +0.5*matrix(rep(qp[i,set],g),ncol=sel,byrow=TRUE)*ag }
lamq[,set] <- lam[,set]+t0

if (count==1) {lbinit <- vlbfix(n,p,s2,s1,g,nsum,sigbeta,alpa,lama,alpb,lamb,alp,lam,mubetaq,
                 sigbetaq,alpaq,lamaq,mubq,sigbq,alpbq,lambq,alpq,lamq,qp,set)$lb}

lbaft <- vlbadd(y,X,V,U,sigbeta,alpa,lama,alpb,lamb,alp,lam,sigdel,
                 mubetaq,sigbetaq,muaq,sigaq,alpaq,lamaq,mubq,sigbq,
                 alpbq,lambq,alpq,lamq,mudelq,qp,vni,epm,set)$lb
lb <- lbinit+lbaft
   
dif <- (lb-lbold)/abs(lb)
lbold <- lb

cat(lb,k,round(apply(qp,2,sum),1),dif,"\n") }

if (k>1) {
gcc <- multinom(qp ~ U-1, decay=1/sigdel[1,1],Hess=TRUE,trace=FALSE)
sigdelq <- solve(gcc$Hess)
lbadj <- lb + vlbadj(U,alpa,sigdel,sigdelq)$lbadj 
} else {
sigdelq=NULL
lbadj=lb}

cat(k,lb,lbadj,"\n")

list (mubetaq=mubetaq,sigbetaq=sigbetaq,muaq=muaq,sigaq=sigaq,mubq=mubq,
sigbq=sigbq,alpaq=alpaq,lamaq=lamaq,alpbq=alpbq,lambq=lambq,
alpq=alpq,lamq=lamq,qp=qp,mudelq=mudelq,lb=lb,lbadj=lbadj,sigdelq=sigdelq)
}




