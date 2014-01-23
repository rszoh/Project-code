### Prior Checks code
## function to simulate Orthogonal Matrices
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
## simulate orthogonal matrices
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

### lower triangular mat
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

## simulate A
A <- triangl(k1,d,c(1,10),c(2,1))
B <- triangl(k2,d,c(1,10),c(2,1))

### orthoganl matrix
mat1 <- sim_psimat(k=k1,p,psi_v=.5) ## output a p by k matrix ( k < p)
mat2 <- sim_psimat(k=k2,q,psi_v=.5) ## output a p by k matrix ( k < p)
## Approximate distribution of correlation estimates
N=10000
out <-NULL

for(i in 1:N)
{
Vm<- riwish(d+2,diag(rep(1,d)))
A <- triangl(k1,d,c(1,20),c(1,.1)); B <- triangl(k2,d,c(1,20),c(1,.1))
var_te <- crossprod(t(mat1%*%A)) + diag(sigte)
var_lb <- crossprod(t(mat2%*%B)) + diag(siglb)
cov_lbte <- crossprod(t(mat1%*%A),t(mat2%*%B))
fn1 <- cov2cor(rbind(cbind(var_te,cov_lbte),cbind(t(cov_lbte),var_lb)))
out=rbind(out,fn1[1,(p+1):(p+q)])
}  
#head(out[,1:5])
par(mfrow=c(3,3))
truehist(out[,1])
truehist(out[,2])
truehist(out[,3])
truehist(out[,4])
truehist(out[,5])
truehist(out[,6])
truehist(out[,7])
truehist(out[,8])
truehist(out[,9])
####
x=seq(.01,.999,length.out=100)
y=x/(1/5 + x)
par(mfrow=c(1,1))
plot(x,y,type="l",lwd=4)











