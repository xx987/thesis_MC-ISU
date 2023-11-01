
M=1500
xv = 1.8 #1.5，1.57，1.65， 1.8
set.seed(214)
start.time<-Sys.time()


mean_v <- c(0,0)
cov_v <- matrix(15*0.56*c(5.473653,-9.035177, -9.035177, 14.980733), nrow = 2)#matrix(c(2.252578^2, -0.9707999, -0.9707999, 1.247183^2), nrow = 2)

samples <- rmvnorm(M, mean = mean_v, sigma = cov_v)

sigm<-1/rgamma(M,1,0.1)

alpha<-samples[,1]
beta<-samples[,2]
hfunc <- function(two,one){(two + one*xv)}
prop_den <- dmvnorm(samples,mean=mean_v, sigma = cov_v)
#prop_sig <- 1/sigm/log(7/3)#dnorm(sigm,4.859392,0.1110831)
prop_sig <-dgamma(sigm,1,0.1)


post_beta <- function(alpha, beta, lambda, sigma){
  result<-dmvnorm(c(beta, alpha),(lambda/(lambda+1))*c(b_hat),(lambda/(lambda+1))*(sigma^2)*coe)
  return(result)
}

#post_beta<-function(alpha,beta,gs,sigma){ #this is prior p(beta, alpha | lambda)
#result<-dmvnorm(c(alpha, beta),c(0,0),gs*sigma*coe)
#return(result)

#}

inv_bet<-c(-39.06196, 61.27219)
post_sigma <- function(sigma,lambda){
  result<- dgamma(1/sigma^2, length(x)/2, (svar/2)+ (1/(2*(lambda+1))) * (inv_bet%*%cx %*% t(cx)%*%inv_bet))
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
    
    nuo<-nuo+hfunc(samples[,1][i],samples[,2][i])*poster(lambda,alpha[i],beta[i],sigm[i])/(prop_den[i]*prop_sig[i])
    #nuo<-hfunc(lambda)*poster(lambda)/prop_den
    #nuo<-hfunc(lambda)*poster(lambda,alpha,beta,sigm)/(prop_den*prop_sig)
  }
  return(nuo)
}

deo=0

denominator<-function(lambda){
  for (i in 1:M){
    
    #wei_d<-append(wei_d,1)#poster(lambda,alpha[i],beta[i],sigm[i])/(prop_den[i]*prop_sig[i]))
    deo<-deo+poster(lambda,alpha[i],beta[i],sigm[i])/(prop_den[i]*prop_sig[i])
    #deo<-poster(lambda)/prop_den
    #deo<-poster(lambda,alpha,beta,sigm)/(prop_den*prop_sig)
  }
  return(deo)
}







llamb = as.integer(range(sampled_lambda)[1])
ulamb = as.integer(range(sampled_lambda)[2])




resapp<-c()
#theo_v<-c()
grid <-seq(llamb, ulamb, length.out =15)#by = 1#lambda

for(i in 1:15){
  res = numerator(grid[i])/denominator(grid[i])##this is the result of E(h(theta) | lambda, D)
  resapp<-append(resapp,res)
}

mc_isu=na.omit(resapp)
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken

wei_d <- c()
all_v <- c()
for (i in 1:15){
  wei_d <- append(wei_d,max(all_v)/mean(all_v))
  for(j in 1:M){
    all_v<-append(all_v,poster(grid[i],alpha[j],beta[j],sigm[j])/(prop_den[j]*prop_sig[j]))
  }}
#for(i in 1:15){
#ret = (grid[i]/(1+grid[i]))*(b_hat[1]*xv+b_hat[2])
#theo_v <-append(theo_v,ret)
#}

ng<-c()
xi<-c()
sdcomp<-function(lambda){
  ###compute ng and xi on variance formula
  #ng<-(denominator(lambda)-mean(denominator(lambda)))/mean(denominator(lambda))
  #ng<-(denominator(lambda)-(denominator(lambda))/M)#/(denominator(lambda)/M)
  #ng<-append(ng,(lambda-(denominator(lambda))/M)/(denominator(lambda)/M))
  #Eng<-mean(ng)
  Eng<-ng/M
  #xi<-(numerator(lambda)-mean(numerator(lambda)))/mean(numerator(lambda))#[1,2,3,4,5,6]-3
  #Exi<-mean(xi)
  #xi<-(numerator(lambda)-(numerator(lambda))/M)#/(numerator(lambda)/M)#[1,2,3,4,5,6]-3
  #xi<-append(xi,(lambda-(numerator(lambda))/M)/(numerator(lambda)/M))
  Exi<-xi/M
  ####compute the E(xi) and E(yi) on the formula
  #Enu<-mean(numerator(lambda))
  #Ede<-mean(denominator(lambda))
  Enu<-numerator(lambda)/M#numerator does not contain h for divided
  Ede<-denominator(lambda)/M
  #result<-(Enu^2/Ede^2)*mean((ng-xi)^2)/M
  #result<-(Enu^2/Ede^2)*((ng-xi)^2)/M/M
  #result<-(Enu^2/Ede^2)*((ng-xi)^2)/M
  #result<-(Enu^2/Ede^2)*((ng-xi))/M
  #range(var(ng-xi))
  result<-Enu^2/Ede^2
  return(result)
}

for(i in 1:15){
  ng<-append(ng,(grid[i]-(denominator(grid[i]))/M)/(denominator(grid[i])/M))
  xi<-append(xi,(grid[i]-(numerator(grid[i]))/M)/(numerator(grid[i])/M))
}


sdr<-c()
for(i in 1:15){
  res = sdcomp(grid[i])/var(ng-xi)/M
  sdr<-append(sdr,res)
}


