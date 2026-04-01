library(nnet)
#library(devtools)
#devtools::install_github("Polarolouis/blockmodels", force = TRUE)
library(blockmodels)
library("sbm")


n=200
p = 3
K = 3
X = matrix(1,n,p+1)
for(j in 2:(p+1)){
  X[,j] = rnorm(n,0,1/2)
}

Beta = round(matrix(abs(rnorm((p+1)*K,0,3)),(p+1),K))
Beta[,1] = 0
Sign_Beta  = matrix(0,p+1,K)
for (i in 1:(p+1))
  {
  for (k in 1:K){
    Sign_Beta[i,k]=(-1)^(k+i)}
}
Beta  = Beta*Sign_Beta


#------------------------
  
TAU = exp(X%*% Beta)
Z = sapply(1:n,function(i){sample(1:3,size=1,replace=FALSE,TAU[i,])})



myTau_1 <- TAU*exp(matrix(rnorm(n*K,0,1),n,K))

nu = 1000; 
myTau_2 = t(sapply(1:n,function(i){
  u =  rdirichlet(1,nu*TAU[i,]/sum(TAU[i,])); 
  return(u)}
  )
)


myTau_1 <- myTau_1/myTau_1[,1]
myTau_2 <- myTau_2/myTau_2[,1]

head(TAU)
head(myTau_1)
head(myTau_2)

myBeta_lin <- solve(t(X)%*%X)%*%t(X)%*%log(TAU)
myBeta_1_lin <- solve(t(X)%*%X)%*%t(X)%*%log(myTau_1)
myBeta_2_lin <- solve(t(X)%*%X)%*%t(X)%*%log(myTau_2)




##### 

res_estim <- multinom(TAU ~ X[,2] + X[,3]+X[,4])
res_estim_myTau1 <- multinom(myTau_1 ~ X[,2] + X[,3] + X[,4])
res_estim_myTau2 <- multinom(myTau_2 ~ X[,2] + X[,3] + X[,4])

res_estim_Z <- multinom(Z ~ X[,2] + X[,3]+X[,4])


#Beta_estim_multinom <- cbind(rep(0,p+1),t(summary(res_estim)$coefficients)) 
#Beta_estim_multinom_Z <- cbind(rep(0,p+1),t(summary(res_estim_Z)$coefficients)) 
# plot((X%*%Beta_estim_multinom), (X%*%Beta_estim_multinom_Z), main='Compar TAU / Z')



Beta_estim_multinom_2<- cbind(rep(0,p+1),t(summary(res_estim_myTau2)$coefficients)) 
Beta_estim_multinom_1<- cbind(rep(0,p+1),t(summary(res_estim_myTau1)$coefficients)) 


#---------------------------------------- 
par(mfrow=c(2,2))


plot(log(TAU), (X%*%myBeta_lin), main='Exact softmax')
points(log(TAU), (X%*%Beta_estim_multinom),col='green')
abline(0,1)

plot(log(myTau_1), (X%*%myBeta_1_lin), main='Presque softmax',col='red')
points(log(myTau_1), (X%*%Beta_estim_multinom_1),col='green')
abline(0,1)

plot(log(myTau_2), (X%*%Beta_estim_multinom_2),col='green',main = "Dirichlet",xlim=range(X%*%Beta_estim_multinom_2))
points(log(myTau_2),(X%*%myBeta_2_lin),col='red')
abline(0,1)


par(mfrow=c(2,2))

plot(Beta, myBeta_lin, main='Exact softmax',type='p',col='red')
points(Beta,Beta_estim_multinom, main='Exact softmax',type='p',col='green')
abline(0,1)


plot(Beta, myBeta_1_lin, main="Presque softmax",type='p',col='red')
points(Beta,Beta_estim_multinom_1,type='p',col='green')
abline(0,1)

plot(Beta,myBeta_2_lin, main="Dirichlet",type='p',col='red')
points(Beta,Beta_estim_multinom_2,type='p',col='green')
abline(0,1)

print(c(sum((Beta-myBeta_2_lin)^2),sum((Beta-Beta_estim_multinom_2)^2)))
print(c(sum((Beta-myBeta_1_lin)^2),sum((Beta-Beta_estim_multinom_1)^2)))
print(c(sum((Beta-myBeta_lin)^2),sum((Beta-Beta_estim_multinom)^2)))





############" SBM 
myTau_1 <- exp(X%*% Beta + matrix(rnorm(n*K,0,1.4),n,K))
myTau_1 <- myTau_1/myTau_1[,1]
Z = sapply(1:n,function(i){sample(1:3,size=1,replace=FALSE,TAU[i,])})

alpha = round(matrix(rbeta(K*K,0.9,0.9),K,K),1)
alpha = (alpha+t(alpha))/2
ord <- order(apply(alpha,1,sum),decreasing = TRUE)
alpha = alpha[ord,ord]

Y <- matrix(rbinom(n*n,1,alpha[Z,Z]),n,n)

res_SBM <- estimateSimpleSBM(Y,model="bernoulli")
alpha_estim <- res_SBM$connectParam$mean
ord_2 <- order(apply(alpha_estim,1,sum),decreasing = TRUE)
tau_estim <- res_SBM$probMemberships[,ord_2]
Z_estim <-apply(res_SBM$indMemberships[,ord_2],1,which.max)


Beta_estim_SBM <- solve(t(X)%*%X)%*%t(X)%*%log(tau_estim/tau_estim[,1])
Beta_estim_SBM_2 <- cbind(rep(0,p+1),t(summary(multinom(tau_estim ~ X[,2] + X[,3]+X[,4]))$coefficients))
Beta_estim_SBM_3 <- cbind(rep(0,p+1),t(summary(multinom(Z_estim ~ X[,2] + X[,3]+X[,4]))$coefficients))






par(mfrow=c(1,1))
plot(Beta,Beta_estim_multinom)
abline(a=0,b=1)
points(Beta,Beta_estim_SBM,col='red')
points(Beta,Beta_estim_SBM_2,col='blue',pch=2)
points(Beta,Beta_estim_SBM_3,col='green',pch=3)

#plot(log(tau_estim_2),log(TAU))

R_SBM = sum((Beta-Beta_estim_SBM)^2)
R_SBM_2 = sum((Beta-Beta_estim_SBM_2)^2)
R_SBM_3 = sum((Beta-Beta_estim_SBM_3)^2)
print(c(R_SBM,R_SBM_2,R_SBM_3))


table(Z,Z_estim)


