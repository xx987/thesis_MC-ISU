library(mvtnorm)

#M=3500
M=3500
xv = 1.65
set.seed(214)
#set.seed(100)

mean_v <- c( -38.37779, 59.93598)#c(-36.64903, 55.67207)
cov_v <- matrix(c(5.128097^2,-0.9654197, -0.9654197, 3.062571^2), nrow = 2)#matrix(c(2.252578^2, -0.9707999, -0.9707999, 1.247183^2), nrow = 2)
samples <- rmvnorm(M, mean = mean_v, sigma = cov_v)
sigm <- rnorm(M, 4.859392,0.1110831)
alpha<-samples[,1]
beta<-samples[,2]
#sig<- 0.1142345 #mean(sampled_sigma)
hfunc <- function(lambda){(samples[,1] + samples[,2]*xv)*(lambda/(1+lambda))}
#h_prior  <- function(lambda){dmvnorm(samples, mean = c(0,0), sigma = lambda*sig*coe)}
prop_den <- dmvnorm(samples,mean=mean_v, sigma = cov_v)
prop_sig <- dnorm(sigm,4.859392,0.1110831)
start.time <- Sys.time()




#######This is the f(Beta, sigma^2 | D, lambda)\
post_beta <- function(alpha, beta, lambda, sigma){
  result<-dmvnorm(c(beta, alpha),(lambda/(lambda+1))*c(b_hat),(lambda/(lambda+1))*(sigma^2)*coe)
  return(result)
}

#post_beta<-function(alpha,beta,gs,sigma){ #this is prior p(beta, alpha | lambda)
#result<-dmvnorm(c(alpha, beta),c(0,0),gs*sigma*coe)
#return(result)

#}


post_sigma <- function(sigma,lambda){
  result<- dgamma(1/sigma^2, length(x)/2, (svar/2)+ (1/(2*(lambda+1))) * (c(b_hat)%*%cx %*% t(cx)%*%c(b_hat)))
  return(result)
}




poster <- function(lambda,alp,bet,sig) {
  
  #result<-dmvt(c(alp,bet), (lambda/(1+lambda))*c(b_hat), ((svar/2)+ (1/(2*(lambda+1))) * (c(b_hat)%*%cx %*% t(cx)%*%c(b_hat)))[1]*coe/7.5, 15, log=FALSE)
  result<-post_beta(alp, bet, lambda, sig)*post_sigma(sig,lambda)
  return(result)
}


########



nuo = 0
numerator<-function(lambda){
  
  for (i in 1:M){
    
    nuo<-nuo+hfunc(lambda)*poster(lambda,alpha[i],beta[i],sigm[i])/(prop_den*prop_sig)
    #nuo<-hfunc(lambda)*poster(lambda)/prop_den
    #nuo<-hfunc(lambda)*poster(lambda,alpha,beta,sigm)/(prop_den*prop_sig)
  }
  return(nuo)
}

deo=0
denominator<-function(lambda){
  for (i in 1:M){
    deo<-deo+poster(lambda,alpha[i],beta[i],sigm[i])/(prop_den*prop_sig)
    #deo<-poster(lambda)/prop_den
    #deo<-poster(lambda,alpha,beta,sigm)/(prop_den*prop_sig)
  }
  return(deo)
}

resapp<-c()
retsee<-c()
grid <-seq(38, 50, by = 1)#lambda
weight<-c()
for(i in 1:13){
  res = sum(numerator(grid[i]))/sum(denominator(grid[i]))##this is the result of E(h(theta) | lambda, D)
  #res = numerator(grid[i])/denominator(grid[i])
  we<-na.omit(denominator(grid[i]))
  we_res<-max(we)/mean(we)
  weight<-append(weight,we_res)
  resapp<-append(resapp,res)
}

resee=na.omit(resapp)
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken
for(i in 1:13){
  ret = (grid[i]/(1+grid[i]))*(b_hat[1]*xv+b_hat[2])
  retsee <-append(retsee,ret)
}

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
for(i in 1:13){
  res = sdcomp(grid[i])
  sdr<-append(sdr,res)
}


retsee
resee
max(retsee)-min(retsee)
max(resee)-min(resee)

lambda = seq(38, 50, length.out = length(sampled_alpha))
err = rnorm(length(sampled_alpha),0, 4.8)
h = sampled_beta*1.8 + sampled_alpha+err
#plot(lambda,h)



plot(lambda,h,ylab = "Estimated_y")
lines(grid, resee,col = "green",lwd=2)
lines(grid, retsee,col = "red",lwd=2)
legend("topright", legend = c("MCMC", "h","MC-ISU"), col = c("black", "red","green"), lwd = c(NA, 2,3), pch = c(1, NA,NA), cex = 0.7)


