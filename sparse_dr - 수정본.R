library(mclust)
library(tidyverse)
library(dr)
#library(bdsmatrix)
library(flextable)
library(glmnet)


#----------------------------------------------------------------------------------------------
# supporting functions
#----------------------------------------------------------------------------------------------

# center a vector
center_vec<-function(v)  v - mean(v)

# normalize a vector
norm<-function(v)  
{ 
  sumv2<-sum(v^2)
  if(sumv2 == 0) sumv2<-1
  v/sqrt(sumv2)
}

# Gram-Schmidt orthonormalization
orthnormal<-function(X)
{
  X<-as.matrix(X)
  n<-nrow(X)
  p<-ncol(X)
  
  W<-NULL
  if(p > 1) {
    W<-cbind(W, X[,1])
    for(k in 2:p) {
      gw<-rep(0, n)
      for(i in 1:(k-1)) {
        gki<-as.vector((t(W[,i]) %*% X[,k])/(t(W[,i]) %*% W[,i]))
        gw<-gw + gki * W[,i]
      }
      W<-cbind(W, X[,k] - gw)
    }
  } else {
    W<-cbind(W, X[,1])
  }
  
  W<-apply(W, 2, norm)
  W
}

# covariance matrix 
cov.x<-function(X)
{
  Xc<-apply(X, 2, center_vec)
  t(Xc) %*% Xc / nrow(Xc)
}

# square-root of a matrix
mat.sqrt<-function(A)
{
  ei<-eigen(A)
  d<-ei$values
  d<-(d+abs(d))/2
  d2<-sqrt(d)
  ans<-ei$vectors %*% diag(d2) %*% t(ei$vectors)
  return(ans)
}

# square-root-inverse of a matrix
mat.sqrt.inv<-function(A)
{
  ei<-eigen(A)
  d<-ei$values
  d<-(d+abs(d))/2
  d2<-1 / sqrt(d)
  d2[d == 0]<-0
  ans<-ei$vectors %*% diag(d2) %*% t(ei$vectors)
  return(ans)
}

# angle between two spaces
angles<-function(B1, B2)
{
  if(!is.matrix(B1)) B1<-as.matrix(B1)
  if(!is.matrix(B2)) B2<-as.matrix(B2)
  
  if(ncol(B1) >= ncol(B2)) {
    B<-B1; B.hat<-B2
  } else {
    B<-B2; B.hat<-B1
  }
  
  P1<-B %*% solve(t(B) %*% B) %*% t(B)
  if(ncol(B.hat) == 1) {
    nume<-as.vector(t(B.hat) %*% P1 %*% B.hat)
    deno<-as.vector(t(B.hat) %*% B.hat)
    ratio<-nume / deno
  } else {
    BtB<-t(B.hat) %*% B.hat
    ei<-eigen(BtB)
    BtB2<-ei$vectors %*% diag(1/sqrt(ei$values)) %*% t(ei$vectors)
    M<-BtB2 %*% t(B.hat) %*% P1 %*% B.hat %*% BtB2
    ratio<-abs(eigen(M)$values[nrow(M)])
  }
  ans<-acos(sqrt(ratio))/pi * 180
  if(ans > 90) ans<-180 - ans
  return(ans)
}

################################################################
#          Power of matrix
################################################################
matpower=function(a,alpha){
  a = (a + t(a))/2
  tmp = eigen(a)
  return(tmp$vectors%*%diag((tmp$values)^alpha)%*%
           t(tmp$vectors))}


################################################################
#           discretize
################################################################
discretize=function(y,h){
  n=length(y);m=round(n/h)
  y=y+.00001*mean(y)*rnorm(n)
  yord = y[order(y)]
  divpt=numeric();for(i in 1:(h-1)) divpt = c(divpt,yord[i*m+1])
  y1=rep(0,n);y1[y<divpt[1]]=1;y1[y>=divpt[h-1]]=h
  for(i in 2:(h-1)) y1[(y>=divpt[i-1])&(y<divpt[i])]=i
  return(y1)
}

################################################################
#          distance between subspaces
################################################################
dist=function(v1,v2){
  v1=as.matrix(v1);v2=as.matrix(v2)
  if(dim(v1)[1]>1){
    p1 <- v1%*%matpower(t(v1)%*%v1,-1)%*%t(v1)
    p2 <- v2%*%matpower(t(v2)%*%v2,-1)%*%t(v2)
  }
  if(dim(v1)[1]==1){
    p1=v1%*%t(v1)/c(t(v1)%*%v1)
    p2=v2%*%t(v2)/c(t(v2)%*%v2)}
  d <- sqrt(sum((p1-p2)*(p1-p2)))
  return(d)
}

##############################################################
#                   symmtrize a matrix
############################################################## 
symmetry = function(a){
  return((a + t(a))/2)}

#----------------------------------------------------------------------------------------------
# functions for sparse sdr wrt G inner product
#----------------------------------------------------------------------------------------------

comp.sdr<-function(X, y, method, d, nslices)
{
  out.m<-switch(method,
                pc     = comp.sdr.pc(X, y, d), 
                sir    = comp.sdr.sir(X, y, d, nslices),
                save   = comp.sdr.save(X, y, d, nslices),
                phdres = comp.sdr.phdres(X, y, d),
                dr     = comp.sdr.dr(X, y, d, nslices)
  )
  
  ans<-list(m=out.m$m, G=out.m$G, beta.sdr=out.m$beta.sdr)
  return(ans)
}



comp.sdr.pc<-function(X, y, d)
{
  # m matrix
  Xc<-apply(X, 2, center_vec)
  m<-mat.sqrt(t(Xc) %*% Xc)
  
  # G matrix 
  G<-diag(1, ncol(X))
  
  # pc
  v<-eigen(cov.x(X))$vectors[, 1:d]
  if(d == 1) v<-matrix(v, ncol=1)
  
  # return
  ans<-list(m=m, G=G, beta.sdr=v)
  return(ans)
}



comp.sdr.sir<-function(X, y, d, nslices)
{
  # parameters
  n<-nrow(X)
  p<-ncol(X)
  Sigma.x<-cov.x(X)
  Sigma.x2<-mat.sqrt(Sigma.x)
  Sigma.x.inv2<-mat.sqrt.inv(Sigma.x)
  
  # standardize X
  Z<-apply(X, 2, center_vec) %*% Sigma.x.inv2
  
  # slice y
  sy<-dr.slices(y, nslices)
  nslices<-sy$nslices
  
  # compute sdr kernel matrix
  M.sir.z<-matrix(0, nrow=p, ncol=p)
  for(s in 1:nslices) {
    Z.s<-Z[sy$slice.indicator == s, ]
    if(sy$slice.sizes[s] == 1) Z.s<-matrix(Z.s, nrow=1)
    Z.sm<-as.vector(apply(Z.s, 2, mean))
    M.sir.z<-M.sir.z + (sy$slice.sizes[s]/n) * Z.sm %*% t(Z.sm) #sirmat과 동일한 candidate matrix
  }
  M.sir<-Sigma.x2 %*% M.sir.z %*% Sigma.x2
  #m<-mat.sqrt(n*M.sir)
  m<-mat.sqrt(M.sir) #M^(1/2) -> mi를 구하기 위해(ssdr을 위해서)
  
  # compute sdr estimate w/o shrinkage(그냥 sdr을 위해서)
  v<-eigen(M.sir.z)$vectors[,1:d]
  if(d == 1) v<-matrix(v, ncol=1)
  beta.sdr<-Sigma.x.inv2 %*% v
  beta.sdr<-apply(beta.sdr, 2, norm)
  
  # return
  ans<-list(m=m, G=Sigma.x, beta.sdr=beta.sdr)
  return(ans)
}



comp.sdr.save<-function(X, y, d, nslices)
{
  # parameters
  n<-nrow(X)
  p<-ncol(X)
  Sigma.x<-cov.x(X)
  Sigma.x2<-mat.sqrt(Sigma.x)
  Sigma.x.inv2<-mat.sqrt.inv(Sigma.x)
  
  # standardize X
  Z<-apply(X, 2, center_vec) %*% Sigma.x.inv2
  
  # slice y
  sy<-dr.slices(y, nslices)
  nslices<-sy$nslices
  
  # compute sdr kernel matrix
  M.save.z<-matrix(0, nrow=p, ncol=p)
  for(s in 1:nslices) {
    Z.s<-Z[sy$slice.indicator == s, ]
    if(sy$slice.sizes[s] == 1) Z.s<-matrix(Z.s, nrow=1)
    iVz<-diag(1, p) - cov.x(Z.s) 
    M.save.z<-M.save.z + (sy$slice.sizes[s]/n) * iVz %*% iVz 
  }
  M.save<-Sigma.x2 %*% M.save.z %*% Sigma.x2
  m<-mat.sqrt(M.save)
  
  # compute sdr estimate w/o shrinkage
  v<-eigen(M.save.z)$vectors[,1:d]
  if(d == 1) v<-matrix(v, ncol=1)
  beta.sdr<-Sigma.x.inv2 %*% v
  beta.sdr<-apply(beta.sdr, 2, norm)
  
  # return
  ans<-list(m=m, G=Sigma.x, beta.sdr=beta.sdr)
  return(ans)
}



comp.sdr.phdres<-function(X, y, d)
{
  # parameters
  n<-nrow(X)
  p<-ncol(X)
  Sigma.x<-cov.x(X)
  Sigma.x2<-mat.sqrt(Sigma.x)
  Sigma.x.inv2<-mat.sqrt.inv(Sigma.x)
  
  # standardize X
  Z<-apply(X, 2, center_vec) %*% Sigma.x.inv2
  
  # residual
  e<-resid(lm(as.vector(y)~Z))
  
  # compute sdr kernel matrix
  M.phd.z<-matrix(0, nrow=p, ncol=p)
  for(s in 1:n) {
    M.phd.z<-M.phd.z + e[s] * Z[s,] %*% t(Z[s,])
  }
  M.phd.z<-M.phd.z / n
  M.phd<-Sigma.x2 %*% (M.phd.z %*% t(M.phd.z)) %*% Sigma.x2
  m<-mat.sqrt(M.phd)
  
  # compute sdr estimate w/o shrinkage
  v<-eigen(M.phd.z %*% t(M.phd.z))$vectors[,1:d]
  if(d == 1) v<-matrix(v, ncol=1)
  beta.sdr<-Sigma.x.inv2 %*% v
  beta.sdr<-apply(beta.sdr, 2, norm)
  
  # return
  ans<-list(m=m, G=Sigma.x, beta.sdr=beta.sdr)
  return(ans)
}


comp.sdr.dr <- function(X,y,d,nslices){
  p=ncol(X);n=nrow(X)
  Sigma.x=cov.x(X)
  Sigma.x2=mat.sqrt(Sigma.x)
  signrt=matpower(var(X),-1/2)
  xc=t(t(X)-apply(X,2,mean))
  xst=xc%*%signrt
  ydis=discretize(y,nslices)
  ylabel=unique(ydis)
  prob=numeric() 
  for(i in 1:nslices) prob=c(prob,length(ydis[ydis==ylabel[i]])/n)
  vxy = array(0,c(p,p,nslices));exy=numeric()
  for(i in 1:nslices) {
    vxy[,,i]=var(xst[ydis==ylabel[i],])
    exy=rbind(exy,apply(xst[ydis==ylabel[i],],2,mean))}
  mat1 = matrix(0,p,p);mat2 = matrix(0,p,p)
  for(i in 1:nslices){
    mat1 = mat1+prob[i]*(vxy[,,i]+exy[i,]%*%t(exy[i,]))%*%
      (vxy[,,i]+exy[i,]%*%t(exy[i,]))
    mat2 = mat2+prob[i]*exy[i,]%*%t(exy[i,])}
  out = 2*mat1+2*mat2%*%mat2+2*sum(diag(mat2))*mat2-2*diag(p)
  M.dr = Sigma.x2 %*% out %*% Sigma.x2 
  m = mat.sqrt(M.dr)
  ans <- list(m=m, G=Sigma.x, beta.sdr =signrt%*%eigen(out)$vectors[,1:d])
  return(ans)
}


# solve beta using lasso estimation
solve.beta<-function(x, y, G2, lambda1, lambda2)
{  
  # transform data to an L1 problem
  x.star<-rbind(x, sqrt(lambda2) * G2) 
  y.star<-c(y, rep(0, ncol(x)))
  
  fit2 <- glmnet(x.star, y.star, family="gaussian")
  beta.est <- coef(fit2, s=lambda1)[-1,]
  
  # return
  return(beta.est)
}



ssdr.lambda<-function(X, y, method=c("pc", "sir", "save", "phdres", "dr"), d=1, nslices=5, lambda1=NULL, lambda2=1e-6, max.iter=200, eps.conv=1e-3)
{
  # parameters
  n<-nrow(X)
  p<-ncol(X)
  method<-match.arg(method)
  
  # compute sdr components
  out.m<-comp.sdr(X, y, method, d, nslices)
  m<-out.m$m
  M<-t(m) %*% m   
  G<-out.m$G
  G2<-mat.sqrt(G)
  G2.inv<-mat.sqrt.inv(G2)
  
  # initial estimate of alpha and beta
  alpha<-out.m$beta.sdr      
  beta<-alpha      
  for(i in 1:d) {
    ym<-m %*% alpha[,i]
    beta[,i]<-solve.beta(m, ym, G2, lambda1[i], lambda2)
    #fit1 <- glmnet(m, ym, family="gaussian")
    #beta[,i] <- coef(fit1, s=lambda1[i])[-1,]
  }
  
  # iteration
  iter<-0
  beta.n<-apply(beta, 2, norm)
  diff.conv<-1
  while((iter < max.iter) & (diff.conv[iter+1] > eps.conv)){
    z<-svd(G2.inv %*% M %*% beta)
    alpha<-G2.inv %*% (z$u) %*% t(z$v)
    for(i in 1:d) {
      ym<-m %*% alpha[,i]
      beta[,i]<-solve.beta(m, ym, G, lambda1[i], lambda2)
      #fit1 <- glmnet(m, ym, family="gaussian")
      #beta[,i] <- coef(fit1, s=lambda1[i])[-1,]
    }
    
    beta.n.new<-apply(beta, 2, norm)
    diff.conv<-c(diff.conv, max(abs(beta.n.new - beta.n)))
    beta.n<-beta.n.new
    
    iter<-iter + 1
  }
  
  # compute objective value
  comp1<-m %*% solve(G) - m %*% beta %*% t(beta)
  rss1<-sum(diag(comp1 %*% G %*% t(comp1)))
  
  z<-svd(G2.inv %*% M %*% beta)
  alpha<-G2.inv %*% (z$u) %*% t(z$v)
  comp2<-m %*% solve(G) - m %*% beta %*% t(alpha)
  rss2<-sum(diag(comp2 %*% G %*% t(comp2)))
  
  p.e<-sum(as.vector(beta) != 0)
  
  # normalize beta
  beta<-apply(beta, 2, norm)
  
  # return
  ans<-list(beta=beta, beta0=out.m$beta.sdr, alpha=alpha, diff.conv=diff.conv, iter=iter, p.e=p.e, rss1=rss1, rss2=rss2)
  return(ans)
}


ssdr.wrap<-function(X, y, method=c("pc", "sir", "save", "phdres","dr"), d=1, nslices=5, s1.range=NULL, s2.range=NULL, max.iter=200, eps.conv=1e-3)
{
  # parameters 
  n<-nrow(X)
  
  # compute criteria for all s1 and s2
  crit.all<-list()
  for(j in 1:length(s2.range)) {
    s2<-s2.range[j]
    
    crit<-NULL
    for(k in 1:length(s1.range)) {
      s1<-s1.range[k]
      
      out1<-ssdr.lambda(X, y, method=method, d=d, nslices=nslices, lambda1=rep(s1, d), lambda2=s2)
      
      aic1<-n*out1$rss1 + 2 * out1$p.e
      bic1<-n*out1$rss1 + log(n) * out1$p.e
      out.crit<-c(aic1, bic1)
      crit<-cbind(crit, out.crit) 
    }
    
    crit.all[[j]]<-crit
  }
  
  # locate optimal s 
  s1.min<-NULL
  ct.min<-NULL
  for(j in 1:length(s2.range)) {
    s1.min.j<-ct.min.j<-NULL
    for(l in 1:length(out.crit)) {
      s1.min.j<-c(s1.min.j, s1.range[order(crit.all[[j]][l,])[1]])
      ct.min.j<-c(ct.min.j, min(crit.all[[j]][l,], na.rm=T))
    }
    s1.min<-rbind(s1.min, s1.min.j)
    ct.min<-rbind(ct.min, ct.min.j)
  }
  rownames(s1.min)<-as.character(s2.range); colnames(s1.min)<-c("aic1", "bic1")
  rownames(ct.min)<-as.character(s2.range); colnames(ct.min)<-c("aic1", "bic1")
  
  # beta estimate with given optimal s
  beta.est.all<-NULL
  s12.est.all<-NULL
  for(l in 1:length(out.crit)) {
    pos<-order(ct.min[, l])[1]
    s1<-s1.min[pos, l]
    s2<-s2.range[pos]
    out1<-ssdr.lambda(X, y, method=method, d=d, nslices=nslices, lambda1=rep(s1, d), lambda2=s2)
    beta.est.all<-cbind(beta.est.all, out1$beta)
    s12.est.all <-cbind(s12.est.all,  c(s1, s2))
  }
  rownames(s12.est.all)<-c("s1", "s2"); colnames(s12.est.all)<-c("aic1", "bic1")
  
  # return
  ans<-list(beta.est.all=beta.est.all, beta.est0=out1$beta0, s12.est.all=s12.est.all, crit.all=crit.all, s1.min=s1.min, ct.min=ct.min)
  return(ans)
}


###########################################################################
###########################################################################



set.seed(2022)
############################################
########## simulation
########## y = sign(b1'X)*log(|b2'X+5|)+0.2*e
########## (1)
############################################

n = 200
p = 20
sigma = 0.2
x = matrix(rnorm(n*p),n,p)
dim(x)

############################################
########## (1)
############################################

beta.true=matrix(0,p,2)
beta.true[1:4,1] = 1; beta.true[17:20,2]=1; beta.true
beta.true.norm <- apply(beta.true, 2, norm)
dim(beta.true)

y1 = sign(beta.true[,1]%*%t(x))*log(abs(beta.true[,2]%*%t(x)+5)) + sigma*rnorm(n)
dim(y1)

############### choose the order of beta ################
beta.order <- function(beta){
  new.beta <- beta
  if (abs(beta[1,1]) < abs(beta[1,2])){
    new.beta <- matrix(c(beta[,2], beta[,1]), ncol=2)
  }
  return(new.beta)
}


##  first trial
lambda1 = c(0.01,0.01)
ssdr.result1 <- ssdr.lambda(x, t(y1), method="dr", d=2, nslices=4, lambda1, lambda2=0.01, max.iter=200, eps.conv=1e-7)
beta = beta.order(ssdr.result1$beta)
beta
cor((x%*%beta.true)[,1], (x%*%ssdr.result1$beta)[,1]) # cor
cor((x%*%beta.true)[,2], (x%*%ssdr.result1$beta)[,2])
mean(((x%*%beta.true.norm)[,1] - (x%*%ssdr.result1$beta)[,1])^2) # mse
mean(((x%*%beta.true.norm)[,2] - (x%*%ssdr.result1$beta)[,2])^2)

## choose lambda : s12.est.all의 결과 보면 됨
ssdr.wrap.result1 <- ssdr.wrap(x, t(y1), method="dr", d=2, nslices=4, s1.range=range(0.001, 0.01, 0.05, 0.1), s2.range=c(0.01, 0.1, 0.5), max.iter=200, eps.conv=1e-7)
ssdr.wrap.result1 # 0.001, 0.01


lambda1 = c(0.001,0.001)
ssdr.result1 <- ssdr.lambda(x, t(y1), method="dr", d=2, nslices=4, lambda1, lambda2=0.01, max.iter=200, eps.conv=1e-7)
ssdr.result1$beta


# simulation 1
beta.true=matrix(0,p,2)
beta.true[1:4,1] = 1; beta.true[17:20,2]=1
beta.true.norm <- apply(beta.true, 2, norm)
p.sdr <- c(); p.dr <- c()
sdr.cor1 <- c(); sdr.cor2 <- c()
dr.cor1 <- c(); dr.cor2 <- c()
sdr.mse1 <- c(); sdr.mse2 <- c()
dr.mse1 <- c(); dr.mse2 <- c()
niter <- 200
set.seed(2022)
for (i in 1:niter){
  y = sign(x%*%beta.true[,1])*log(abs(x%*%beta.true[,2]+5)) + sigma*rnorm(n)
  ssdr.result1 <- ssdr.lambda(x, y, method="dr", d=2, nslices=4, lambda1=c(0.01, 0.01), lambda2=0.01, max.iter=200, eps.conv=1e-3)
  dr.beta <- beta.order(comp.sdr.dr(x, y, d=2, nslices=4)$beta.sdr)
  beta <- beta.order(ssdr.result1$beta)
  
  # number of p
  p.sdr <- c(p.sdr, 40-ssdr.result1$p.e)
  p.dr <- c(p.dr, sum(dr.beta==0))
  
  # corr with true beta
  sdr.cor1 <- c(sdr.cor1, abs(cor((x%*%beta.true)[,1], (x%*%beta)[,1])))
  dr.cor1 <- c(dr.cor1, abs(cor((x%*%beta.true)[,1], (x%*%dr.beta)[,1])))
  sdr.cor2 <- c(sdr.cor2, abs(cor((x%*%beta.true)[,2], (x%*%beta)[,2])))
  dr.cor2 <- c(dr.cor2, abs(cor((x%*%beta.true)[,2], (x%*%dr.beta)[,2])))
  
  # mse with true x
  sdr.mse1 <- c(sdr.mse1, mean(((x%*%beta.true.norm)[,1] - (x%*%beta)[,1])^2))
  dr.mse1 <- c(dr.mse1, mean(((x%*%beta.true.norm)[,1] - (x%*%dr.beta)[,1])^2))
  sdr.mse2 <- c(sdr.mse2, mean(((x%*%beta.true.norm)[,2] - (x%*%beta)[,2])^2))
  dr.mse2 <- c(dr.mse2, mean(((x%*%beta.true.norm)[,2] - (x%*%dr.beta)[,2])^2))
}

result.sim1 <- data.frame(method = c("dr", "sparse dr"), p = c(mean(p.dr), mean(p.sdr)), corr1 = c(mean(dr.cor1), mean(sdr.cor1)), corr2 = c(mean(dr.cor2), mean(sdr.cor2)), mse1 = c(mean(dr.mse1), mean(sdr.mse1)), mse2 = c(mean(dr.mse2), mean(sdr.mse2)))
flextable(result.sim1)


############################################
########## second model
############################################

beta.true=matrix(0,p,2)
beta.true[1:2,1] = 1; beta.true[3:4,1] = 0.1
beta.true[19:20,2]=1; beta.true[17:18,2]=0.1; beta.true
beta.true.norm <- apply(beta.true, 2, norm)
dim(beta.true)

y2 = sign(beta.true[,1]%*%t(x))*log(abs(beta.true[,2]%*%t(x)+5)) + sigma*rnorm(n)
dim(y2)

## first trial
lambda1 = c(0.01,0.01)
ssdr.result2 <- ssdr.lambda(x, t(y2), method="dr", d=2, nslices=5, lambda1, lambda2=0.1, max.iter=200, eps.conv=1e-7)
beta = beta.order(ssdr.result2$beta)
beta
cor((x%*%beta.true)[,1], (x%*%ssdr.result2$beta)[,1]) # cor
cor((x%*%beta.true)[,2], (x%*%ssdr.result2$beta)[,2])
mean(((x%*%beta.true.norm)[,1] - (x%*%ssdr.result2$beta)[,1])^2) # mse
mean(((x%*%beta.true.norm)[,2] - (x%*%ssdr.result2$beta)[,2])^2)

## choose lambda
ssdr.wrap.result2 <- ssdr.wrap(x, t(y2), method="dr", d=2, nslices=5, s1.range=range(0.001, 0.01, 0.05, 0.1, 0.5), s2.range=c(0.001, 0.01, 0.1, 0.5), max.iter=200, eps.conv=1e-7)
ssdr.wrap.result2 # 0.05, 0.01


# simulation 2
beta.true=matrix(0,p,2)
beta.true[1:2,1] = 1; beta.true[3:4,1] = 0.1
beta.true[19:20,2]=1; beta.true[17:18,2]=0.1
beta.true.norm <- apply(beta.true, 2, norm)
p.sdr <- c(); p.dr <- c()
sdr.cor1 <- c(); sdr.cor2 <- c()
dr.cor1 <- c(); dr.cor2 <- c()
sdr.mse1 <- c(); sdr.mse2 <- c()
dr.mse1 <- c(); dr.mse2 <- c()
niter <- 200

for (i in 1:niter){
  y = sign(x%*%beta.true[,1])*log(abs(x%*%beta.true[,2]+5)) + sigma*rnorm(n)
  ssdr.result1 <- ssdr.lambda(x, y, method="dr", d=2, nslices=4, lambda1=c(0.01, 0.01), lambda2=0.01, max.iter=200, eps.conv=1e-3)
  dr.beta <- beta.order(comp.sdr.dr(x, y, d=2, nslices=4)$beta.sdr)
  beta <- beta.order(ssdr.result1$beta)
  
  # number of p
  p.sdr <- c(p.sdr, 40-ssdr.result1$p.e)
  p.dr <- c(p.dr, sum(dr.beta==0))
  
  # corr with true beta
  sdr.cor1 <- c(sdr.cor1, abs(cor((x%*%beta.true)[,1], (x%*%beta)[,1])))
  dr.cor1 <- c(dr.cor1, abs(cor((x%*%beta.true)[,1], (x%*%dr.beta)[,1])))
  sdr.cor2 <- c(sdr.cor2, abs(cor((x%*%beta.true)[,2], (x%*%beta)[,2])))
  dr.cor2 <- c(dr.cor2, abs(cor((x%*%beta.true)[,2], (x%*%dr.beta)[,2])))
  
  # mse with true x
  sdr.mse1 <- c(sdr.mse1, mean(((x%*%beta.true.norm)[,1] - (x%*%beta)[,1])^2))
  dr.mse1 <- c(dr.mse1, mean(((x%*%beta.true.norm)[,1] - (x%*%dr.beta)[,1])^2))
  sdr.mse2 <- c(sdr.mse2, mean(((x%*%beta.true.norm)[,2] - (x%*%beta)[,2])^2))
  dr.mse2 <- c(dr.mse2, mean(((x%*%beta.true.norm)[,2] - (x%*%dr.beta)[,2])^2))
}

result.sim2 <- data.frame(method = c("dr", "sparse dr"), p = c(mean(p.dr), mean(p.sdr)), corr1 = c(mean(dr.cor1), mean(sdr.cor1)), corr2 = c(mean(dr.cor2), mean(sdr.cor2)), mse1 = c(mean(dr.mse1), mean(sdr.mse1)), mse2 = c(mean(dr.mse2), mean(sdr.mse2)))
flextable(result.sim2)



############################################
########## third model
############################################

beta.true=matrix(0,p,2)
beta.true[1:10,1] = 1
beta.true[11:20,2]= 1; beta.true
beta.true.norm <- apply(beta.true, 2, norm)
dim(beta.true)

y3 = sign(beta.true[,1]%*%t(x))*log(abs(beta.true[,2]%*%t(x)+5)) + sigma*rnorm(n)
dim(y3)



## first trial
lambda1 = c(0.01,0.01)
ssdr.result3 <- ssdr.lambda(x, t(y3), method="dr", d=2, nslices=5, lambda1, lambda2=0.1, max.iter=200, eps.conv=1e-7)
beta = beta.order(ssdr.result3$beta)
beta
cor((x%*%beta.true)[,1], (x%*%beta)[,1]) # cor
cor((x%*%beta.true)[,2], (x%*%beta)[,2])
mean(((x%*%beta.true.norm)[,1] - (x%*%beta)[,1])^2) # mse
mean(((x%*%beta.true.norm)[,2] - (x%*%beta)[,2])^2)

## choose lambda
ssdr.wrap.result1 <- ssdr.wrap(x, t(y3), method="dr", d=2, nslices=5, s1.range=range(0.001, 0.01, 0.1, 0.5), s2.range=c(0.01, 0.1, 0.5), max.iter=200, eps.conv=1e-7)
ssdr.wrap.result1 # 0.001, 0.01


# simulation 3
beta.true=matrix(0,p,2)
beta.true[1:10,1] = 1
beta.true[11:20,2]= 1
beta.true.norm <- apply(beta.true, 2, norm)
p.sdr <- c(); p.dr <- c()
sdr.cor1 <- c(); sdr.cor2 <- c()
dr.cor1 <- c(); dr.cor2 <- c()
sdr.mse1 <- c(); sdr.mse2 <- c()
dr.mse1 <- c(); dr.mse2 <- c()
niter <- 200

for (i in 1:niter){
  y = sign(x%*%beta.true[,1])*log(abs(x%*%beta.true[,2]+5)) + sigma*rnorm(n)
  ssdr.result1 <- ssdr.lambda(x, y, method="dr", d=2, nslices=4, lambda1=c(0.01, 0.01), lambda2=0.01, max.iter=200, eps.conv=1e-3)
  dr.beta <- beta.order(comp.sdr.dr(x, y, d=2, nslices=4)$beta.sdr)
  beta <- beta.order(ssdr.result1$beta)
  
  # number of p
  p.sdr <- c(p.sdr, 40-ssdr.result1$p.e)
  p.dr <- c(p.dr, sum(dr.beta==0))
  
  # corr with true beta
  sdr.cor1 <- c(sdr.cor1, abs(cor((x%*%beta.true)[,1], (x%*%beta)[,1])))
  dr.cor1 <- c(dr.cor1, abs(cor((x%*%beta.true)[,1], (x%*%dr.beta)[,1])))
  sdr.cor2 <- c(sdr.cor2, abs(cor((x%*%beta.true)[,2], (x%*%beta)[,2])))
  dr.cor2 <- c(dr.cor2, abs(cor((x%*%beta.true)[,2], (x%*%dr.beta)[,2])))
  
  # mse with true x
  sdr.mse1 <- c(sdr.mse1, mean(((x%*%beta.true.norm)[,1] - (x%*%beta)[,1])^2))
  dr.mse1 <- c(dr.mse1, mean(((x%*%beta.true.norm)[,1] - (x%*%dr.beta)[,1])^2))
  sdr.mse2 <- c(sdr.mse2, mean(((x%*%beta.true.norm)[,2] - (x%*%beta)[,2])^2))
  dr.mse2 <- c(dr.mse2, mean(((x%*%beta.true.norm)[,2] - (x%*%dr.beta)[,2])^2))
}

result.sim3 <- data.frame(method = c("dr", "sparse dr"), p = c(mean(p.dr), mean(p.sdr)), corr1 = c(mean(dr.cor1), mean(sdr.cor1)), corr2 = c(mean(dr.cor2), mean(sdr.cor2)), mse1 = c(mean(dr.mse1), mean(sdr.mse1)), mse2 = c(mean(dr.mse2), mean(sdr.mse2)))
flextable(result.sim3)


###################################################################################
#################### 논문 3.2 ####################################################

set.seed(2022)
############################################
########## simulation
########## y = cos(2*b1'X)-cos(b2'X)+0.5*e
########## (1)
############################################

n = 100
p = 10
sigma = 0.5
x = matrix(rnorm(n*p),n,p)
dim(x)


beta.true=matrix(0,p,2)
beta.true[1,1] = 1; beta.true[2,2]=1; beta.true
beta.true.norm <- apply(beta.true, 2, norm)
dim(beta.true)

y2 = cos(2*beta.true[,1]%*%t(x))-cos(beta.true[,2]%*%t(x)) + sigma*rnorm(n)
dim(y2)

############### choose the order of beta ################
beta.order <- function(beta){
  new.beta <- beta
  if (abs(beta[1,1]) < abs(beta[1,2])){
    new.beta <- matrix(c(beta[,2], beta[,1]), ncol=2)
  }
  return(new.beta)
}


##  first trial
lambda1 = c(0.01,0.01)
ssdr.result1 <- ssdr.lambda(x, t(y2), method="dr", d=2, nslices=4, lambda1, lambda2=0.01, max.iter=100, eps.conv=1e-7)
beta = beta.order(ssdr.result1$beta)
beta
cor((x%*%beta.true)[,1], (x%*%ssdr.result1$beta)[,1]) # cor
cor((x%*%beta.true)[,2], (x%*%ssdr.result1$beta)[,2])
mean(((x%*%beta.true.norm)[,1] - (x%*%ssdr.result1$beta)[,1])^2) # mse
mean(((x%*%beta.true.norm)[,2] - (x%*%ssdr.result1$beta)[,2])^2)

## choose lambda : s12.est.all의 결과 보면 됨
ssdr.wrap.result1 <- ssdr.wrap(x, t(y2), method="dr", d=2, nslices=4, s1.range=range(0.001, 0.01, 0.05, 0.1), s2.range=c(0.01, 0.1, 0.5), max.iter=100, eps.conv=1e-7)
ssdr.wrap.result1 # 0.001, 0.01


lambda1 = c(0.001,0.001)
ssdr.result1 <- ssdr.lambda(x, t(y2), method="dr", d=2, nslices=4, lambda1, lambda2=0.01, max.iter=100, eps.conv=1e-7)
ssdr.result1$beta


# simulation 1
beta.true=matrix(0,p,2)
beta.true[1,1] = 1; beta.true[2,2]=1
beta.true.norm <- apply(beta.true, 2, norm)
p.sdr <- c(); p.dr <- c()
sdr.cor1 <- c(); sdr.cor2 <- c()
dr.cor1 <- c(); dr.cor2 <- c()
sdr.mse1 <- c(); sdr.mse2 <- c()
dr.mse1 <- c(); dr.mse2 <- c()
niter <- 100
set.seed(2022)
for (i in 1:niter){
  y = cos(2*beta.true[,1]%*%t(x))-cos(beta.true[,2]%*%t(x)) + sigma*rnorm(n)
  ssdr.result1 <- ssdr.lambda(x, y, method="dr", d=2, nslices=4, lambda1=c(0.01, 0.01), lambda2=0.01, max.iter=100, eps.conv=1e-3)
  dr.beta <- beta.order(comp.sdr.dr(x, y, d=2, nslices=4)$beta.sdr)
  beta <- beta.order(ssdr.result1$beta)
  
  # number of p
  p.sdr <- c(p.sdr, 20-ssdr.result1$p.e)
  p.dr <- c(p.dr, sum(dr.beta==0))
  
  # corr with true beta
  sdr.cor1 <- c(sdr.cor1, abs(cor((x%*%beta.true)[,1], (x%*%beta)[,1])))
  dr.cor1 <- c(dr.cor1, abs(cor((x%*%beta.true)[,1], (x%*%dr.beta)[,1])))
  sdr.cor2 <- c(sdr.cor2, abs(cor((x%*%beta.true)[,2], (x%*%beta)[,2])))
  dr.cor2 <- c(dr.cor2, abs(cor((x%*%beta.true)[,2], (x%*%dr.beta)[,2])))
  
  # mse with true x
  sdr.mse1 <- c(sdr.mse1, mean(((x%*%beta.true.norm)[,1] - (x%*%beta)[,1])^2))
  dr.mse1 <- c(dr.mse1, mean(((x%*%beta.true.norm)[,1] - (x%*%dr.beta)[,1])^2))
  sdr.mse2 <- c(sdr.mse2, mean(((x%*%beta.true.norm)[,2] - (x%*%beta)[,2])^2))
  dr.mse2 <- c(dr.mse2, mean(((x%*%beta.true.norm)[,2] - (x%*%dr.beta)[,2])^2))
}

result.sim1 <- data.frame(method = c("dr", "sparse dr"), p = c(mean(p.dr), mean(p.sdr)), corr1 = c(mean(dr.cor1), mean(sdr.cor1)), corr2 = c(mean(dr.cor2), mean(sdr.cor2)), mse1 = c(mean(dr.mse1), mean(sdr.mse1)), mse2 = c(mean(dr.mse2), mean(sdr.mse2)))
flextable(result.sim1)



