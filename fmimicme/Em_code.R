#### EM algorithm functions
#### Dat is a data matrix with the following columns
#### This is expecting data as follows:
#### col1: i(individual id); col2: j(outcome Outc); col3: d(Group), col4: k(Time), col5: y(Response), col6: Zid(error free covariate)
#### col 7 wid (Heat): average of heat production over time
#### Zdat is data matrix that contains the col1: i(individual); col2: d(Group),; col3: Zid 
library(matlab)
data_summary <-function(dat,Zdat) ## this function outputs some data summary
{                            ## sample size, # of outcomes, # of groups, # of time points for each animal, Zid list of values
samps <- unique(dat[,1]) ## obtaines the animal ids
i_val <- length(samps)

outs <- unique(dat[,2]) ## obtains outcomes ids
j_val <- length(outs)

groups <- unique(dat[,3]) ## obtains groups ids
d_val <- length(groups)

### find the length of time points for each animal
Mijk <-array(NA,c(i_val,j_val,d_val))
Did <-list() ### Did all observed data for individual i
Zid <-list()
for(d in 1:d_val)
{
for(i in 1:i_val)
{
  for(j in 1:j_val)
  {
    idx <- (1:nrow(dat))[dat[,1]==samps[i] & dat[,2]== outs[j] & dat[,3]==groups[d]]  
    Mijk[i,j,d] <- length(dat[idx,"Time"])
    }
  }
}

for(i in 1:i_val)
{
Did[[i]] <- c(dat[dat[,"Id"]==i,"y"],mean(dat[idx,"Wid"])) 
idx1 <- (1:nrow(Zdat))[Zdat[,1]==samps[i]]
Zid[[i]] <- Zdat[idx1] 
}
ouf<-list("i"=i_val,"j"=j_val,"d"=d_val,"Mijd"=Mijd,"Did"=Did,"Zid"=Zid) ## returns a list 
return(ouf)
} 


### Computes covariance matrices bit by peaces
## Mijd is an array with the obs for each comb of i, j, d
## Rijd[[i]][[j]][[d]] is basis function values evalued at time tijdk
##lambj is list of matrices, list of length J
## sigj is a vector of length J
## kap1jd (kap1_{jd}) is a matrix with nrow=J and ncol=D
## bet2jd (bet2_{jd}) is a matrix with nrow=J and ncol=D
## sigdnu (sig_{dnu}) is a vector of length D
## id0=c(i0,j0,d0); 
## Cjd is a list of a  matrix of J rows and D columns 
## kap2jd is a matrix of J and D
## al0d: is a vector of length d
## al3d: is a list of length d
# id0=c(i,d) 
# kap2jd is a list of list kap2jd[[j]][[d]] is a a vector
# Cjd is list Cjd[[j]][[d]] is a vector
# al0d is a list al0d[[d]] is a vector
# Qijd is a list Qijd[[i]][[j]][[d]] is a matrix
# Zid is a list Zid[[i]][[d]]
##### variance between y_ijd 
##id0=c(i,j,d) for y_ijd

covvaryijd <-function(Mijd,Rijd,lambj,sigj,kap1jd,bet2jd,sigdnu,id0,id1) # computes the mijd*mijd covariance matrix
{ 
i0=id0[1]; j0=id0[2]; d0 <-id0[3];
i1=id1[1]; j1=id1[2]; d1 <-id1[3];
if(i0==i1 & d0==d1)
{  
out <- (kap1jd[j0,d0]*kap1jd[j1,d1]*(sigdnu[d1]^2) + beta2[j0,d0]*beta2[j1,d1])*ones(c(Mijd[i0,j0,d0],Mijd[i1,j1,d1])) + (j0==j1)*(sigj[j]^2*eye(Mijd[i,j,d])
+ Rijd[[i0]][[j0]][[d0]]%*%lambj[[j]]%*%t(Rijd[[i0]][[j0]][[d0]]))
}
else
{
out <- 0*ones(c(Mijd[i0,j0,d0],Mijd[i1,j1,d1]))
}
return(out)
}

##### covariance between y_ijd and w_id
##id0=c(i,j,d) for y_ijd
## id1=c(i,d) for w_id
covYWijd <-function(Mijd,kap1jd,sigdnu,id0,id1) ## covariance between y_ijd and w_id
{
i=id0[1]; j=id0[2]; d <-id0[3];  
i1<-id1[1];d1 <-id1[2]
if(i==i1 & d == d1)
{out<- kap1jd[j,d]^2*sigdnu[d]^2*ones(Mijd[id0])}
else
{ out <- 0*ones(Mijd[id0])}
return(out)
}

##### variance of wid
covwid <-function(sigu,sigdnu,id0) ## id0=c(i,j d)
{
  d=id0[3] 
return(sigu^2 + sigdnu[d])  
}

### covariance between Y-ijd and X_id
##id0=c(i,j,d) for y_ijd
## id1=c(i,d) for X_id
covYXijd <-function(Rijd,lambj,Mijd,kap1jd,sigdnu,id0,id1)
{ 
out <- prod(id0[-2] == id1)*kap1jd[id0[2:3]]*sigdnu[id0[3]]^2*ones(Mijd[id0]) + prod(id0[-2] != id1)*ones(Mijd[id0])*0  
return(out)
}
##id0=c(i,j,d) for y_ijd
## id1=c(i,j,d) for Mis_ijd
covYMisijd <- function(Mijd,kap1jd,sigdnu,bet2jd,Rijd,lambj,id0,id1)
{ 
  ot1 <- prod(id0[-2] == id1)*kap1jd[id0[2:3]]*sigdnu[id0[3]]^2*ones(c(Mijd[id0],1)) ### cov_Y_X
  ot2 <- prod(id0[-2]==id1)*bet2jd[id0[-1]]*ones(c(Mijd[id0],1)) ## cov_Y_eta
  ot3 <- prod(id0==id1)*Rijd[[id0[1]]][[id0[2]]][[id0[3]]]%*%lambj[[id0[2]]]  ## cov_Y_b
  out=cbind(ot1,cbind(ot2,ot3))
  return(out)
}
### Covariance between W_id and the vector of missing obs
# id0=c(i,j,d)
# id1=c(i,j,d)
covWMisijd <-function(sigdnu,lambj,id0,id1)
{
ot1 <- prod(id0[-2] == id1)*sigdnu[id0[3]]^2
ot2 <- 0
ot3 <- rep(0,nrow(lambj[[id0[3]]])) ## determine the appropriate length of b_ijd
out=c(ot1,ot2,ot3)
return(out)  
}

### covariance between mis Sig_11
# id0=c(i,j,d)
Sig11 <-function(sigdnu,lambj,id0)
{
dm <- nrow(lambj[[id0[2]]]) + 2
sig <-matrix(0,nrow=dm,ncol=dm)
sig[1,1] <- sigdnu[id0[3]]^2
sig[2,2] <- 1
sig[3:dm,3:dm] <- lambj[[id0[2]]]
return(sig)
}
### mu2: mean of the obseverd values
# id0=c(i,d) 
# kap2jd is a list of list kap2jd[[j]][[d]] is a a vector
# Cjd is list Cjd[[j]][[d]] is a vector
# al0d is a list al0d[[d]] is a vector
# Qijd is a list Qijd[[i]][[j]][[d]] is a matrix
# Zid is a list Zid[[i]][[d]]
###
mu2id <-function(Cjd,Qijd,kap2jd,Mijd,Zid,J,al0d,al3d,id0)
{
outy=NULL
vc <- as.vector(Zid[[id0[1]]][[id0[2]]]) 
for(j in 1:J)
{
ix=c(id0[1],j,id0[2])
vx <- rep(sum((vc)*kap2jd[[j]][[id0[2]]]),Mijd[ix]) + Qijd[[id0[1]]][[j]][[id0[2]]]%*%Cjd[[j]][[id0[2]]]
outy = c(outy,vx)    
}
vw = alod[id0[2]] + sum(vc*al3d[[id0[2]]])
return(c(vc,vw))
}
#################################################################################################
#################################################################################################




