###unchanged dataset x,y

install.packages("Rlab")
library(Rlab)

set.seed(214)
xset <- runif(500, -1, 1)
yset <- c()
ng <-function(u,theori,x){(exp(u)/(1+exp(u)))*(pnorm(theori*x,0,1))+(1/(1+exp(u)))*(exp(theori*x*pi/sqrt(3))/(1+exp(theori*x*pi/sqrt(3))))}
for (i in 1:500){
  uo<-log(0.5/(1-0.5))
  para<- ng(uo,0.7,xset[i])
  y<-rbern(1,para)
  yset=append(yset,y)
}
head(xset)
head(yset)

##
#seedlist = c(100)#,200,500)
#sdlist = c(0.1067892)#2*0.1067892,3*0.1067892)#c(0.1067892,2,3,2*0.1067892,3*0.1067892)
#for (sed in seedlist){
  #print(i)
  #for (sd in sdlist){

start.time <- Sys.time()
##x is fixed value 2/3#####This x is different from fdfunc (likelihood)
M = 5000#10000
x = 2/3
#for checking convenience, we use q(theta) as N(1.15,0.8), because the original theta is sampling from N(0,1)
set.seed(214)
#theta=rnorm(M,0.7424202,0.1067892)
#qtheta=dnorm(theta,0.7424202,0.1067892)
theta=rnorm(M,0.7314861,0.1067703)
qtheta=dnorm(theta,0.7314861,0.1067703)

###compute Berger's weight directly
#w = target/qtheta
#wrest = max(w)/mean(w)
#############
#hfunc<-function(u){(exp(u)/(1+exp(u)))*(pnorm(theta*x,0,1))+(1/(1+exp(u)))*(exp(theta*x*pi/sqrt(3))/(1+exp(theta*x*pi/sqrt(3))))}##define the original hfunc here : 把原版h的function定义好
hfunc<-function(lambda){lambda*(pnorm(theta*x,0,1))+(1-lambda)*(exp(theta*x*pi/sqrt(3))/(1+exp(theta*x*pi/sqrt(3))))}##define the original hfunc here
####define the likelihood function ####定义单个似然函数 f(D|theta, lambda)######
#hsecf<-function(u,xi){(exp(u)/(1+exp(u)))*(pnorm(theta*xi,0,1))+(1/(1+exp(u)))*(exp(theta*xi*pi/sqrt(3))/(1+exp(theta*xi*pi/sqrt(3))))}###only use for likehood (convenience)
hsecf<-function(lambda,xi){lambda*(pnorm(theta*xi,0,1))+(1-lambda)*(exp(theta*xi*pi/sqrt(3))/(1+exp(theta*xi*pi/sqrt(3))))}###only use for likehood (convenience)
#generate the xset and yset of size 500 for computing likelihood function.
#generate the xset and yset of size 500 for computing likelihood function.
# define nature lambda_0 and theta_0 for sampling y
#lamori <- 0.5
#theori <- 0.7
###we have xset and yset til now, the length of xset and yset is 500



##joint all of parameter function will be likelihood function
product<-1
fdfunc<-function(lambda){
  for (i in 1:500){
    xi=xset[i]
    y=yset[i]
    #likefun<-((lambda)*(pnorm(theta*xi*1.81,0,1))+(1/(1+(lambda/(1-lambda))))*(exp(theta*xi*1.81)/(1+exp(theta*xi*1.81))))^y*((1- (lambda)*(pnorm(theta*xi*1.81,0,1))-(1/(1+(lambda/(1-lambda))))*(exp(theta*xi*1.81)/(1+exp(theta*xi*1.81))))^(1-y))
    product<-product * hsecf(lambda,xi)^y*(1-hsecf(lambda,xi))^(1-y)       #likfun(lambda,xi,y)
    #product <- product * (hfunc(u)^y)*((1-hfunc(u))^(1-y))
    #print(product)
  }
  return(product)
}


#####define the f(theta | lambda)
numerator<-function(lambda){
  nume <- hfunc(lambda)*dnorm(theta)*fdfunc(lambda)/qtheta##this is numerator: h(theta)*(f(theta | labmda, D) / q(theta)
  return(nume)
}

denominator<-function(lambda){
  deno <- dnorm(theta)*fdfunc(lambda)/qtheta##this is denominator:(f(theta | labmda, D) / q(theta), which is W(theta)
  return(deno)
}

resapp<-c()
weight<-c()####now try to compute the weight to each lambda
grid <-seq(0, 1, by = 0.1)#lambda
#grid <- log(lambda/(1-lambda))# u 
for(i in 1:11){
  res = sum(numerator(grid[i]))/sum(denominator(grid[i]))##this is the result of E(h(theta) | lambda, D)
  we<-na.omit(denominator(grid[i]))
  we_res<-max(we)/mean(we)
  weight<-append(weight,we_res)
  resapp<-append(resapp,res)
}
resee=na.omit(resapp)
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken###compute the all time needed



# Print the elapsed time



sdcomp<-function(lambda){
  ###compute ng and xi on variance formula
  ng<-(denominator(lambda)-mean(denominator(lambda)))/mean(denominator(lambda))
  Eng<-mean(ng)
  xi<-(numerator(lambda)-mean(numerator(lambda)))/mean(numerator(lambda))#[1,2,3,4,5,6]-3
  Exi<-mean(xi)
  ####compute the E(xi) and E(yi) on the formula
  Enu<-mean(numerator(lambda))
  Ede<-mean(denominator(lambda))
  result<-(Enu^2/Ede^2)*mean((ng-xi)^2)/M
  
  return(sqrt(result))
  
  
  
}



sdr<-c()
for(i in 1:10){
  res = sdcomp(grid[i])
  sdr<-append(sdr,res)
}
print(paste("Seed = ", sed,"SD = ",sd))
print(min(resee))
print(max(resee))
print(weight)
print("SDR")
print(sdr)

print("_____________________________")

plot(s.f,v.f,xlab="lambda",ylab="h")
lines(grid, resee,col = "red",lwd=2)
legend("topright", legend = c("MCMC","MC-ISU"), col = c("black", "red"), lwd = c(NA, 2), pch = c(1, NA), cex = 0.7)




#relation = lm(y_sample ~ x_sample)
#summary(relation)
######
dat <- structure(list(theta_s = c(resapp), lambda_s = c(grid)), Names = c("theta_s", "lambda_s"),
                 class = "data.frame")

lmObject <- lm(theta_s ~ lambda_s, data = dat)

newdat <- data.frame(lambda_s = c(0, 1))

predict(lmObject, newdat, se.fit = TRUE, interval = "confidence", level = 0.90)

fit_mcl = 0.6984672
fit_mch = 0.6994005
se_mcl = 0.0005367917
se_mch = 0.0005399235
fit_ISl = 0.6976417
fit_ISh = 0.6992977
se_ISl = 2.440159e-05
se_ISh = 2.440159e-05

Z_h = (fit_mch - fit_ISh)/sqrt(se_mch^2 + se_ISh^2)
Z_l = (fit_mcl - fit_ISl)/sqrt(se_mcl^2 + se_ISl^2)

CI_mcl = c(0.6975842,0.6993502)
CI_mch = c(0.6985123,0.7002887)
CI_ISl = c(0.697597,0.6976865)
CI_ISh = c(0.699253,0.6993425)

