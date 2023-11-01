library(coda)
#############################
#example 2主程序
#（利用theta=0.7 和 lambda=0.5（即u=0）产生500个数据进行拟合）
##############################
##get the starting time
library(Rlab)

#以下生成数据    #generate the data of x,y

m <- 100000    #length of chain
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
prop.s<-c(1,10)  # sd of proposal normal 方差
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

burn=95000

u.f=beta[,2][(burn+1):m]          #burn后u链
theta.f=beta[,1][(burn+1):m]       #burn后theta链
s.f=exp(u.f)/(1+exp(u.f))          #burn后lambda链





#x.ss=2/3            #取x值为2/3
v.f=F(theta.f, u.f, x.ss)     #v.f是长度为5000的概率链
summary(v.f)            #列出h概率链的性质
summary(s.f)            #列出lambda链的性质


###compute the all time needed

#theta_ess <- effectiveSize(as.mcmc(v.f))
#theta_ese <- sd(v.f)/sqrt(theta_ess)

###compute ess and ese of u
#u_ess <- effectiveSize(as.mcmc(s.f))
#u_ese <- sd(s.f)/sqrt(u_ess)  
#xl = min(s.f)
#xu = max(s.f)

#index = seq(1, 90000, by = 10)
#y_sample = v.f[index]
#x_sample = s.f[index]
#relation = lm(y_sample ~ x_sample)
#summary(relation)
######
#dat <- structure(list(theta_s = c(y_sample), lambda_s = c(x_sample)), Names = c("theta_s", "lambda_s"),
                 #class = "data.frame")

#lmObject <- lm(theta_s ~ lambda_s, data = dat)

#newdat <- data.frame(lambda_s = c(xu, xl))

#predict(lmObject, newdat, se.fit = TRUE, interval = "confidence", level = 0.90)








#############################################
#画关于h和h真实值h0的直方图
h=v.f
h0=F(theta.t, u.t, x.ss)
#hist(h,col="grey")
#par(new=T)
#abline(v=h0,col="blue")
#蓝色竖线是h0
#############################################
#画关于lambda和lambda真实值lambda.t的直方图
lambda=s.f
#hist(lambda,col="grey")
#par(new=T)
#abline(v=lambda.t,col="blue")
#蓝色竖线是lambda.t
#####################################
#####
#来个loess
#plx=predict(loess(v.f~s.f),se=T)
#################
#对生成的MCMC链按lambda的大小排序

ab1 <- theta.f
ab2<- s.f
new=data.frame(ab1,ab2)
line<-new[order(new$ab2),]

line.theta=t(line[,1])     #即为排序好后的theta链 
line.lambda=t(line[,2])   #排序好的lambda链

#########################

#把5000个数据分成34组求各组最值（每组147个数据）（分34组是因为组数等于数据数目的立方根）

row.t=147
col.t=34
rc=row.t*col.t

sf=matrix(line.lambda[1:rc], nrow=row.t)
ssf=t(sf)     #转置，每一行是一组lam值，共34行，每行147个。

###############

tf=matrix(line.theta[1:rc], nrow=row.t)  #theta分为34组
ttf=t(tf)    #每一行是一组theta值，共34行，每行147个。

vff=matrix(1,nrow=row.t, ncol=col.t)
vf=t(vff)

for (j in 1:col.t){
  vf[j,]=F(ttf[j,], log(ssf[j,]/(1-ssf[j,])), x.ss)
}   
#vf的每一行即为每一组概率值,共34行

Evff.sort=matrix(1,nrow=1,ncol=col.t)
Evf.sort=t(Evff.sort)
for (j in 1:col.t){
  Evf.sort[j,]=sum(vf[j,])/row.t
}     
# Evf.sort[j,]就是每一组的概率期望。
Evf.sort
sort(Evf.sort)
###############################
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken
#画图
####################################
#（h关于lambda 链的散点图）

figure12=function(){
  plot(s.f,v.f,xlim=c(0,1),ylim=c(0.60,0.80),xlab="lambda",ylab="h")
  par(new=T)
  plot(loess.smooth(s.f,v.f)$x,loess.smooth(s.f,v.f)$y,col="blue",xlim=c(0,1),ylim=c(0.60,0.80),pch=21,xlab="lambda",ylab="h",lty=1,type="l")
  par(new=T)
  abline(0,0)
  legend('bottomleft', c("loess","mcmc"), col=c("blue","black"),lty=c(1,0),
         pch=c(-1,1))
}
figure12()

loess.smooth(s.f,v.f)
seeres = loess.smooth(s.f,v.f)$y
