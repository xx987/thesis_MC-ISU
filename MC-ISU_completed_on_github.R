library(mvtnorm)

#
#This is dataset, and we should make x to be n*2 matrix which is (x,1) for matrix manipulation
x <- c(1.47, 1.50, 1.52, 1.55, 1.57, 1.60, 1.63, 1.65, 1.68, 1.70, 1.73, 1.75, 1.78, 1.80, 1.83)
onec <-  rep(1, times = 15)
y <- c(52.21, 53.12, 54.48, 55.84, 57.20, 58.57, 59.93, 61.29, 63.11, 64.47, 66.28, 68.10, 69.92, 72.19, 74.46)
cx <- rbind(x,onec)
coe <- solve(cx %*% t(cx))##this is (x^T*x)^-1
# Set the number of iterations and burn-in period for MCMC
num_iterations <- 20000#num_iterations <- 30000
#burn_in <- 0#5000#burn_in <- 1000
b_hat <- coe%*% cx %*% y
svar <-t((y-t(cx)%*%b_hat))%*% ((y-t(cx)%*%b_hat))

set.seed(214)
start.time <- Sys.time()
#quantile(1/rgamma(1000,1/2,7.5), probs = c(0.35,0.65))
# Initialize empty vectors for storing the sampled parameter values
sampled_alpha <- numeric(num_iterations)
sampled_beta <- numeric(num_iterations)
sampled_lambda <- numeric(num_iterations)
sampled_sigma <- numeric(num_iterations)

# Initialize the starting values for alpha, beta, and sigma
current_alpha <- -39#-35#-30
current_beta <- 50#55#55
current_lambda<-10
current_sigma <- 3 #rgamma(1, shape = 1/2, rate = 1/2)

prior_lam <- function(gs) {###this is the prior of g, which is p(g) or p(lambda)
  #result<-gs^(-3/2)*exp(-length(x)/(2*gs))
  result<- dgamma(1/gs, 1/2, length(x)/2)##inverse gamma distribution
  return(result)
}

post_beta <- function(alpha, beta, lambda, sigma){
  result<-dmvnorm(c(beta, alpha),(lambda/(lambda+1))*c(b_hat),(lambda/(lambda+1))*(sigma^2)*coe)
  #result<-dmvnorm(c(beta, alpha),(lambda/(lambda+1))*c_b,(lambda/(lambda+1))*(sigma^2)*coe)
  return(result)
}

#post_beta<-function(alpha,beta,gs,sigma){ #this is prior p(beta, alpha | lambda)
#result<-dmvnorm(c(alpha, beta),c(0,0),gs*sigma*coe)
#return(result)

#}

inv_bet<-c(-39.06196, 61.27219)
post_sigma <- function(sigma,lambda){
  #result<- dgamma(1/sigma^2, length(x)/2, (svar/2)+ (1/(2*(lambda+1))) * (c(b_hat)%*%cx %*% t(cx)%*%c(b_hat)))
  result<- dgamma(1/sigma^2, length(x)/2, (svar/2)+ (1/(2*(lambda+1))) * (inv_bet%*%cx %*% t(cx)%*%inv_bet))
  return(result)
}



# Perform MCMC sampling
for (i in 1:num_iterations) {
  # Propose new values for alpha, beta, and sigma
  proposed_alpha <- rnorm(1, mean = current_alpha, sd =3)
  proposed_beta <- rnorm(1, mean = current_beta, sd = 3)
  proposed_lambda <- runif(1,min = 20, max = 70)
  #proposed_sigma <- rgamma(1, shape = 1/2, rate = 0.5)#1/rgamma(1, shape = length(x)/2, rate = (svar/2)+ (1/(2*(lambda+1))) * (c(b_hat)%*%cx %*% t(cx)%*%c(b_hat)))
  #proposed_sigma <- runif(1,min = 0, max=5)
  proposed_sigma <- rgamma(1, 2, 1.5)#runif(1,min = 0, max=10000)
  
  
  log_likelihood_proposed <-prior_lam(proposed_lambda) * post_beta(proposed_alpha, proposed_beta,proposed_lambda,proposed_sigma)*post_sigma(proposed_sigma,proposed_lambda)
  
  
  
  log_likelihood_current <- prior_lam(current_lambda) * post_beta(current_alpha, current_beta,current_lambda,current_sigma)*post_sigma(current_sigma,current_lambda)
  flag = log_likelihood_proposed/log_likelihood_current
  
  # Calculate the acceptance ratio
  if (!is.nan(flag)){
    acceptance_ratio <- min(1,log_likelihood_proposed/log_likelihood_current)}
  
  # Accept or reject the proposed values
  if (runif(1) < acceptance_ratio) {
    current_alpha <- proposed_alpha
    current_beta <- proposed_beta
    current_lambda <- proposed_lambda
    current_sigma <- proposed_sigma
  }
  sampled_alpha[i] <- current_alpha
  sampled_beta[i] <- current_beta
  sampled_lambda[i] <- current_lambda
  sampled_sigma[i]<-current_sigma
  
  # Store the current values after burn-in
}


# Print the posterior means
sampled_alpha=sampled_alpha[10000:num_iterations]
sampled_beta=sampled_beta[10000:num_iterations]
sampled_lambda=sampled_lambda[10000:num_iterations]
sampled_sigma=sampled_sigma[10000:num_iterations]


mean_beta = mean(sampled_beta)
mean_alpha = mean(sampled_alpha)
mean_sig = mean(sampled_sigma)
sd_sig = sd(sampled_sigma)













#M=3500
M=1500
xv = 1.8 #1.5，1.57，1.65， 1.8
set.seed(214)
#set.seed(100)

mean_v <- c(mean_alpha,mean_beta)#c(-36.64903, 55.67207)

data <- data.frame(al= sampled_alpha,one = sampled_beta)
cov_v <- cov(data)


#cov_v <- matrix(c(4.351714^2,-0.9900172, -0.9900172, 2.710438^2), nrow = 2)#matrix(c(2.252578^2, -0.9707999, -0.9707999, 1.247183^2), nrow = 2)
samples <- rmvnorm(M, mean = mean_v, sigma = cov_v)
sigm <- rnorm(M, mean_sig,sd_sig)
alpha<-samples[,1]
beta<-samples[,2]
#sig<- 0.1142345 #mean(sampled_sigma)
hfunc <- function(two,one){(two + one*xv)}
#h_prior  <- function(lambda){dmvnorm(samples, mean = c(0,0), sigma = lambda*sig*coe)}
prop_den <- dmvnorm(samples,mean=mean_v, sigma = cov_v)
prop_sig <- dnorm(sigm,mean_sig,sd_sig)





#######This is the f(Beta, sigma^2 | D, lambda)\
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
  #we<-na.omit(denominator(grid[i]))
  #we<-wei_d
  #we_res<-max(we)/mean(we)
  #weight<-append(weight,we_res)
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

theo_v = c(68.22343,70.62655)
plot(x_sample,y_sample,ylab = "Estimated_y")#plot MCMC via Wei and Jiang's method
lines(grid, mc_isu,col = "green",lwd=2)
lines(c(20,70),theo_v,col = "red",lwd=2)
lines(loess.smooth(x_sample,y_sample),col = "orange",lwd=2)
legend("topright", legend = c("MCMC", "h","MC-ISU","loess"), col = c("black", "red","green","orange"), lwd = c(NA, 2,3,4), pch = c(1, NA,NA,NA), cex = 0.4)




