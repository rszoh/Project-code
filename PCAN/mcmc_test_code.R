#!/bin/bash
### Test code
#####################
require("far")
library("truncnorm")
require("MCMCpack")
library("dlm")
library("mvtnorm")
library("matlab")
n=10
p=20
q=25
k1=5
k2=7
d<-min(k1,k2)
al0=1
bet0=1
nu0=1
S0=.1
mute=rep(10,k1)
mulb=rep(10,k2)
sigte=rep(.1,p)
siglb=rep(.1,q)
tau<-rep(.1,min(k1,k2))
sigte0 <- rep(10,k1)
siglb0 <- rep(10,k2)
pial<-rep(.5,k1)
pibl <-rep(.5,k2)
alp=1
betp=1

## simulate A and B (lower triangular matrices)
triangl <-function(k,d,pam1,pam2) ##pam1=al1,al1 (parameter of the beta), pam2=al,bet parameters of the gamma
{
  mat <- matrix(0,nrow=k,ncol=d)
  pivl <- rbeta(d,shape1=pam1[1],shape2=pam1[2])
  sig <- rgamma(n=d,shape=pam2[1],scale=pam2[2]) 
  for(i in 1:k)
  {
    for(j in 1:min(i,d)) 
    {
      bin <-rbinom(n=1,size=1,prob=pivl[j])  
      if(j < i)  
      {mat[i,j] <- bin*rnorm(n=1,sd=sqrt(1/sig[j])) }
      else
      {
        mat[i,i] <- rtruncnorm(n=1,a=0,mean=0,sd=sqrt(1/sig[j]))
      }   
    }
  }
  return(mat)
} 

A <- triangl(k1,min(k1,k2),c(1,1),c(1,1))
B <- triangl(k2,min(k1,k2),c(1,1),c(1,1))

## Simulate phi1 and phi2
sim_psimat<-function(k,p,psi_v) ## output a p by k matrix ( k < p)
{
  repeat
  {
    val=runif(k*p)
    prob<-cumsum(c(psi_v^2, 2*psi_v*(1-psi_v), (1-psi_v)^2))
    lw <- (val <= prob[1])
    mid<- (val <= prob[2] & val > prob[1])
    hig<- (val > prob[2])
    ot <- numeric(k*p)
    if(any(lw)) ot[lw]<- -sqrt(1/psi_v)
    if(any(mid)) ot[mid]<- 0
    if(any(hig)) ot[hig]<- sqrt(1/psi_v)
    u <- matrix(ot,nrow=p,ncol=k)
    test <- (prod(abs(La.svd(u)$d) > 1e-08) == 0) 
    if(!test) { break }
  }
  out=orthonormalization(u,basis=F,norm=T)
  return(out)    
}

phi1 <- t(sim_psimat(k=k1,p=p,psi_v=.5)) ## output a p by k matrix ( k < p)
phi2 <- t(sim_psimat(k=k2,p=q,psi_v=.5)) ## output a p by k matrix ( k < p)

## correlation function
corfc <-function(A,B,phi1,phi2,sig1,sig2)
{
cx1= cbind(crossprod(t(A)%*%phi1) +diag(sig1),crossprod(t(A)%*%phi1,t(B)%*%phi2))   
cx2=cbind(crossprod(t(B)%*%phi2,t(A)%*%phi1), crossprod(t(B)%*%phi2)+diag(sig2)) 
return(cov2cor(rbind(cx1,cx2)))
}  

dim(corfc(A,B,phi1,phi2,sigte,siglb)[1:p,(p+1):(p+q)])
######################################################################
##simuate Date sets
Z=matrix(rnorm(n=n*d),nrow=d)
tet=t(t(phi1)%*%A%*%Z + t(phi1)%*%matrix(mute,ncol=1)%*%ones(1,ncol(Z)) + matrix(rnorm(n=n*p,sd=sigte[1]),nrow=p))
lamb=t(t(phi2)%*%B%*%Z + t(phi2)%*%matrix(mulb,ncol=1)%*%ones(1,ncol(Z)) + matrix(rnorm(n=n*q,sd=siglb[1]),nrow=q))

X=apply(exp(tet),c(1,2),rpois,n=1)
Y=apply(exp(lamb),c(1,2),rpois,n=1)

##
tau_sim<-NULL
# sigte_sim <-NULL
# siglb<-NULL
# mute <-NULL
# mulb <-NULL
# A_sim <-NULL
# B_sim <-NULL
# Z_sim <-NULL
# tet_sim<-NULL
# lamb_sim<-NULL
out <-list("tau"=NULL,"sigte"=NULL,"siglb"=NULL,"mute"=NULL,"mulb"=NULL,"A"=NULL,"B"=NULL,"Z"=NULL,"tet"=NULL,"lamb"=NULL,"pial"=NULL,"pibl"=NULL)
idf=NULL
for(i in 1:100)
{
idf <- c(idf,i)
#tau <- tauj(al0,bet0,A,B) 
tau <- apply(matrix(1:ncol(A),ncol=1),1,tauj,al0=al0,bet0=bet0,A=A,B=B)
out$tau=rbind(out$tau,tau)
write.csv(out$tau,file="tau.csv",append=T)
#sigte <-sigtej(nu0,S0,A,Z,tet,mute)
#sigte <- apply(matrix(1:ncol(tet),ncol=1),1,sigtej,nu0=nu0,S0=S0,A=A,Z=Z,tet=tet,mute=mute)
sigte <- apply(matrix(1:ncol(tet),ncol=1),1,sigtej,nu0=nu0,S0=S0,A=A,Z=Z,tet=tet,mute=mute,phi1=phi1)
out$sigte=rbind(out$sigte,sigte)
write.csv(out$sigte,file="sigte.csv")
#siglb <- siglbj(nu0,S0,B,Z,lamb,mulb)
siglb <- apply(matrix(1:ncol(lamb),ncol=1),1,siglbj,nu0=nu0,S0=S0,B=B,Z=Z,lamb=lamb,mulb=mulb,phi2=phi2)
out$siglb=rbind(out$siglb,siglb)
write.csv(out$siglb,file="siglb.csv")
#mute <- mutej(sigte,sigte0,A,Z,tet)
mute <- mutefc(sigte,sigte0,A,Z,tet,phi1)
out$mute=rbind(out$mute,mute)
write.csv(out$mute,file="mutecsv")
#mulb <- mulbl(siglb,siglb0,A,Z,lamb) ## l=1, ... ,q
mulb <-mulbfc(siglb,siglb0,B,Z,lamb,phi2)
out$mulb=rbind(out$mulb,mulb)
write.csv(out$mulb,file="mulb.csv")
#Z<-Zi(sigte,siglb,tet,lamb,A,B,mute,mulb)
Z <- Zi(sigte,siglb,tet,lamb,A,B,mute,mulb,phi1,phi2)
out$Z <- rbind(out$Z,as.vector(Z))
write.csv(out$Z,file="Z.csv")
# for(l in 1:d)
# {
# for(j in l:p)  
# {A[j,l] <- ajl(sigte,mute,A,Z,tet,tau,j,l)}
# for(k in l:q)
# {B[k,l] <- bkl(siglb,mulb,B,Z,lamb,tau,k,l)}
# }
ida <- matrix(cbind(rep(1:nrow(A),each=ncol(A)), rep(c(1:ncol(A)),nrow(A))),ncol=2)
A <- matrix(apply(ida,1,ajl,sigte=sigte,mute=mute,A=A,Z=Z,tet=tet,tau=tau,phi1=phi1,pial=pial),ncol=ncol(A),byrow=T)
out$A=rbind(out$A,as.vector(A))
write.csv(out$A,file="A.csv")

idb <- matrix(cbind(rep(1:nrow(B),each=ncol(B)), rep(c(1:ncol(B)),nrow(B))),ncol=2)
B <- matrix(apply(idb,1,bkl,siglb=siglb,mulb=mulb,B=B,Z=Z,lamb=lamb,tau=tau,phi2=phi2,pibl=pibl),ncol=ncol(B),byrow=T)
out$B=rbind(out$B,as.vector(B))
write.csv(out$B,file="B.csv")
tet <- adap_rej(An=A,Xn=X,Zn=Z,sig_n=sigte,mu=mute,phim=phi1) 
out$tet=rbind(out$tet,as.vector(tet))
write.csv(out$tet,file="teta.csv")
lamb <- adap_rej(An=B,Xn=Y,Zn=Z,sig_n=siglb,mu=mulb,phim=phi2) 
out$lamb=rbind(out$lamb,as.vector(lamb))
write.csv(out$lamb,file="lamb.csv")
#cp <- cp+1
#if(cp > 100) {break}
pial <- apply(matrix(1:ncol(A),ncol=1),1,pifc,alp=alp,betp=betp,A=A) ## update
out$pial=rbind(out$pial,pial)
write.csv(out$pial,file="pi_al.csv")

pibl <- apply(matrix(1:ncol(B),ncol=1),1,pifc,alp=alp,betp=betp,A=B) ## update
out$pibl=rbind(out$pibl,pibl)
write.csv(out$pibl,file="pi_bl.csv")
}

cat("Done !!!! \n ")
#traceplot(mcmc(out$A[,1]))
#traceback()
##  mcmc 
