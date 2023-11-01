# install.packages("Rlab") 
install.packages("coda")
library(coda)
#############################
#example 2主程序
#（利用theta=0.7 和 lambda=0.5（即u=0）产生500个数据进行拟合）
##############################
library(Rlab)

#以下生成数据    #generate the data of x,y

m <- 15000    #length of chain
n=500     # the number of x, y
x.ss=2/3

set.seed(214)
x<-runif(n, min = -1, max = 1)

F<-function(theta,u,x){(exp(u)/(1+exp(u)))*(pnorm(theta*x,0,1))+(1/(1+exp(u)))*(exp(theta*x*pi/sqrt(3))/(1+exp(theta*x*pi/sqrt(3))))}

theta.t=0.7
lambda.t=0.5
Q<-function(z){log(z/(1-z))}
u.t=Q(lambda.t)
p=rep(1:n)
p<-F(theta.t, u.t, x)

y=rep(1:n)

for(i in 1:n){
  y[i]=rbern(1,p[i])
}
#########################
#h的真实值
h0=F(theta.t, u.t, x.ss)
h0
###############################
#以下准备开始生成MCMC链
start.time <- Sys.time()
beta0<-c(2,6)    # initial parameter value : theta和u的初始值
prop.s<-c(1, 10)  # sd of proposal normal 方差
beta <- matrix(nrow=m, ncol=2)
acc.prob <-c(0,0)   # 通过个数
current.beta<-beta0 

for (t in 1:m){
  for (j in 1:2){
    prop.beta<-current.beta
    prop.beta[j]<-rnorm(1, current.beta[j], prop.s[j] ) 
    #都用正态生成, [1] theta  [2] u
    
    prop.eta<-F(prop.beta[1], prop.beta[2], x)
    cur.eta<-F(current.beta[1], current.beta[2], x)
    
    prop.loga=prod((prop.eta)^y, (1-prop.eta)^(1-y))
    current.loga=prod((cur.eta)^y, (1-cur.eta)^(1-y))
    # for(i in 1:n){prop.loga[i]=prod( (prop.eta)^y[i], (1-prop.eta)^(1-y[i]))}
    
    loga1=prop.loga/current.loga
    loga2=(dnorm(prop.beta[1],0,1)*(exp(prop.beta[2])/(1+exp(prop.beta[2]))^2))/(dnorm(current.beta[1],0,1)*(exp(current.beta[2])/(1+exp(current.beta[2]))^2))
    
    loga=min(1, loga1*loga2)
    
    u<-runif(1)
    if (u< loga) {
      current.beta<-prop.beta
      acc.prob[j] <- acc.prob[j]+1
    }
  }
  beta[t,]<-current.beta  #m行两列
}
print(acc.prob/m)
###compute ess and ese of theta


#以上形成完整m=55000的一条MCMC链

################

#减少一些点，burn=50000，形成剩下5000个点的链

#以上形成完整m=55000的一条MCMC链

################

#减少一些点，burn=50000，形成剩下5000个点的链

burn=10000

u.f=beta[,2][(burn+1):m]          #burn后u链
theta.f=beta[,1][(burn+1):m]       #burn后theta链
s.f=exp(u.f)/(1+exp(u.f))          #burn后lambda链





#x.ss=2/3            #取x值为2/3
v.f=F(theta.f, u.f, x.ss)     #v.f是长度为5000的概率链
summary(v.f)            #列出h概率链的性质
summary(s.f)            #列出lambda链的性质
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken






###########################
#The end!

###########################
#The end!
#relation<-lm(v.f~s.f)