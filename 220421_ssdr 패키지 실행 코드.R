library(devtools)
install_github('GijeongNoh/ssdr')
library(ssdr)

###################패키지로 실행시켜보기####################################
library(flextable)
library(glmnet)

n = 300
p = 40
sigma = 0.5
x = matrix(rnorm(n*p),n,p)
dim(x)

############################################
########## (1)
############################################

beta.true=matrix(0,p,2)
beta.true[1,1] = 1; beta.true[2,2]=1; beta.true
beta.true.norm <- apply(beta.true, 2, norm)
dim(beta.true)
y1 = sin(beta.true[,1]%*%t(x))+sin(beta.true[,2]%*%t(x)) + sigma*rnorm(n)
dim(y1)


##  first trial
lambda1 = c(0.01,0.01)
ssdr.result1 <- ssdr.lambda(x, t(y1), method="dr", d=2, nslices=4, lambda1, lambda2=0.01, max.iter=100, eps.conv=1e-7)
beta = beta.order(ssdr.result1$beta)
beta


ssdr.wrap.result1 <- ssdr.wrap(x, t(y1), method="dr", d=2, nslices=4, s1.range=range(0.001, 0.01, 0.05, 0.1), s2.range=c(0.01, 0.1, 0.5), max.iter=100, eps.conv=1e-7)
ssdr.wrap.result1 # 0.001, 0.01


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
  y = sin(x %*% beta.true[,1])+sin(x%*% beta.true[,2]) + sigma*rnorm(n)
  ssdr.result1 <- ssdr.lambda(x, y, method="dr", d=2, nslices=4, lambda1=c(0.01, 0.01), lambda2=0.01, max.iter=100, eps.conv=1e-3)
  dr.beta <- beta.order(comp.sdr.dr(x, y, d=2, nslices=4)$beta.sdr)
  beta <- beta.order(ssdr.result1$beta)
  
  # number of p
  p.sdr <- c(p.sdr, 80-ssdr.result1$p.e)
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

result.sim1 <- data.frame(method = c("dr", "sparse dr"), p = c(mean(p.dr), mean(p.sdr)), corr1 = c(mean(dr.cor1), mean(sdr.cor1)), corr2 = c(mean(dr.cor2), mean(sdr.cor2, na.rm=TRUE)), mse1 = c(mean(dr.mse1), mean(sdr.mse1)), mse2 = c(mean(dr.mse2), mean(sdr.mse2)))
flextable(result.sim1)
