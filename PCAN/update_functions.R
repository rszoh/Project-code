#### Update parameters functions
### Define all the inputs values
### A:matrix of p * d, B: matrix of q*d, tau = vector of length d
### mute:vector of length p, mulb: vector of length q
### tet: matrix of n*p and lamb: matrix of n*q
### sigte: vector of length p, siglb: vector of length q
###
library(dlm)
tauj <-function(al0,bet0,A,B,j) ### j: the column of A
{
#J <-ncol(A)
#val <-numeric(J)
#for(j in 1:J)
#{
maj <- sum( (1:nrow(A)) > j)
mbj <- sum( (1:nrow(B)) > j)
aln <-.5*(maj+mbj) + al0
betn <- (1/bet0) + .5*(sum(A[,j]^2) + sum(B[,j]^2))
val <-rgamma(n=1,shape=aln,scale=1/betn)
#}
return(val)  
}
apply(matrix(1:ncol(A),ncol=1),1,tauj,al0=al0,bet0=bet0,A=A,B=B)


sigtej <-function(nu0,S0,A,Z,tet,mute,j) ## tet is a n*p matrix ,Z is a d*n matrix
{
p <-ncol(tet)
n <- nrow(tet)
#val <-numeric(p)
#for(j in 1:p)
#{
scl <- (nu0*S0 + sum((tet[,j] - t(Z)%*%matrix(A[j,],ncol=1) - mute[j])^2))
val <- scl / rchisq(n=1,df=n+nu0) 
#}
return(val)
}
apply(matrix(1:ncol(tet),ncol=1),1,sigtej,nu0=nu0,S0=S0,A=A,Z=Z,tet=tet,mute=mute)

siglbj <-function(nu0,S0,B,Z,lamb,mulb,j) ## lamb is a n*q matrix ,Z is a d*n matrix
{
  n <- nrow(lamb)
  q <- ncol(lamb)
  val <-numeric(q)
  #for(j in 1:q)
  #{
  scl <- (nu0*S0 + sum( (lamb[,j] - t(Z)%*%matrix(B[j,],ncol=1)- mulb[j])^2))
  val <- scl / rchisq(n=1,df=n+nu0)
  #}
  return(val)
}
apply(matrix(1:ncol(lamb),ncol=1),1,siglbj,nu0=nu0,S0=S0,B=B,Z=Z,lamb=lamb,mulb=mulb)

mutej <-function(sigte,sigte0,A,Z,tet,j) ## j=1, ...,p
{
n=nrow(tet)
p<-ncol(tet)
val <-numeric(p)
#for(j in 1:p)
#{
sig <- ((n/sigte[j]) + (1/sigte0))^(-1)
mn <- (1/sigte[j])*sig*sum(tet[,j] - t(Z)%*%matrix(A[j,],ncol=1))
val <-rnorm(n=1,mean=mn,sd=sqrt(sig))
#}
return(val)
}
apply(matrix(1:ncol(tet),ncol=1),1,mutej,sigte=sigte,sigte0=sigte0,A=A,Z=Z,tet=tet)

 mulbl <-function(siglb,siglb0,B,Z,lamb,l) ## l=1, ... ,q
{
  n=nrow(lamb)
  q <- ncol(lamb)
  val <-numeric(q)
  #for(l in 1:q)
  #{
  sig <- ((n/siglb[l]) + (1/siglb0))^(-1)
  mn <- (1/siglb[l])*sig*sum(lamb[,l] - t(Z)%*%matrix(B[l,],ncol=1))
  val<-rnorm(n=1,mean=mn,sd=sqrt(sig))
  #}
  return(val)
}
apply(matrix(1:ncol(lamb),ncol=1),1,mulbl,siglb=siglb,siglb0=siglb0,B=B,Z=Z,lamb=lamb)

ajl <-function(sigte,mute,A,Z,tet,tau,id) ##j=1,...,p and l=1,...,d; Z is d*n matrix; id (j,l) vector
{
j=id[1]
l=id[2]
sig <- ((1/sigte[j])*sum(Z[l,]^2) + tau[l])^{-1}   
mun <- (1/sigte[j])*sig*sum(tet[,j] - mute[j] - t(Z[-l,]%*%matrix(A[j,-l],ncol=1)))
if(j > l)
{val <-rnorm(n=1,mean=mun,sd=sqrt(sig))}
else if(j ==l)
{val <- rtruncnorm(n=1,a=0,mean=mun,sd=sqrt(sig))}
else 
{val <- 0}
return(val)
}
ajl(sigte,mute,A,Z,tet,tau,c(2,1))
ida <- matrix(cbind(rep(1:nrow(A),each=ncol(A)), rep(c(1:ncol(A)),nrow(A))),ncol=2)
matrix(apply(ida,1,ajl,sigte=sigte,mute=mute,A=A,Z=Z,tet=tet,tau=tau),ncol=ncol(A),byrow=T)

############################
bkl <-function(siglb,mulb,B,Z,lamb,tau,id) ##k=1,...,p and l=1,...,d; Z is d*n matrix
{
  k<-id[1]
  l<-id[2]
  sig <- ((1/siglb[k])*sum(Z[l,]^2) + tau[l])^{-1}   
  mun <-(1/siglb[k])*sig*sum(lamb[,k] - mulb[k] - t(Z[-l,]%*%matrix(B[k,-l],ncol=1)))
if(k > l)
  {val <-rnorm(n=1,mean=mun,sd=sqrt(sig))}
else if(k==l)
{val <- rtruncnorm(n=1,a=0,mean=mun,sd=sqrt(sig))}
  else { val <- 0 }
return(val)
}

### 
Zi <-function(sigte,siglb,tet,lamb,A,B,mute,mulb,i)
{
sig = solve(diag(rep(1,ncol(A))) + t(A)%*%diag(1/sigte)%*%A + t(B)%*%diag(1/siglb)%*%B)
mn<-NULL
#for(i in 1:nrow(tet))
mn =as.vector((t(matrix(tet[i,] - mute,ncol=1))%*%diag(1/sigte)%*%A + t(matrix(lamb[i,] - mulb,ncol=1))%*%diag(1/siglb)%*%B)%*%sig)
val <- rmvnorm(n=1,mean=mn,sigma=sig)
return(val)
#return(rmvnorm(n=1,mean=as.vector(mn),sigma=sig))  
}


####
### likelihood function 

likhd<-function(tet,Z,X,sig,A,mu) # Z(d*1) vector; tet(p*1) vector; X(p*1) vector, sig(p*1)vector; A(p*d),mu(p*1)
{
  #gam_m<-diag(gam,nrow=length(gam))
  sigm <- diag(sig)
  res=-(1/2)*t(tet - (A%*%matrix(Z,ncol=1) + mu + X*sig))%*%diag(1/sig)%*%(tet - (A%*%Z + mu + X*sig)) - sum(exp(tet)) ### like lihood minus the exp(-exp(tet)) term
  return(res)  
}
#### rejection ratio
########################################################################################################
## this function uses the arsm
#stop("Lets see what we got! \n " )
adap_rej<-function(An,Xn,Zn,sig_n,mu)
{
  ot=matrix(0,ncol=nrow(An),nrow=nrow(Xn))
  for(j in 1:nrow(Xn))   
  {
    ot[j,] <- arms(runif(nrow(An),-30,30),likhd,function(tet,...) (min(tet) >-30)*(max(tet)<30),n, Z=Zn[,j],A=An,sig=sig_n,X=Xn[j,],mu=mu)[1,] 
  }   
  return(ot)
}




