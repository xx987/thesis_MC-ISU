library(mvtnorm)
#
set.seed(214)
#This is dataset, and we should make x to be n*2 matrix which is (x,1) for matrix manipulation
x <- c(1.47, 1.50, 1.52, 1.55, 1.57, 1.60, 1.63, 1.65, 1.68, 1.70, 1.73, 1.75, 1.78, 1.80, 1.83)
cal <- c(runif(15, 2,5))# artifical k/calories taken per day
onec <-  rep(1, times = 15)
y <- c(52.21, 53.12, 54.48, 55.84, 57.20, 58.57, 59.93, 61.29, 63.11, 64.47, 66.28, 68.10, 69.92, 72.19, 74.46)
cx <- rbind(x,cal,onec)
coe <- solve(cx %*% t(cx))##this is (x^T*x)^-1
# Set the number of iterations and burn-in period for MCMC
num_iterations <- 20000#20000#num_iterations <- 30000
#burn_in <- 0#5000#burn_in <- 1000
b_hat <- coe%*% cx %*% y
svar <-t((y-t(cx)%*%b_hat))%*% ((y-t(cx)%*%b_hat))
inv_bet <- c(b_hat[3],b_hat[2],b_hat[1])





time_c <- c()
inf_exp <- c()
sup_exp <- c()
msdr <- c()
sup_sdr<-c()
inf_sdr <- c()
inf_the <-c()
sup_the <- c()

##########Start to do 100 iterations for testing############
for(l in 1:100){

start.time<-Sys.time()



#quantile(1/rgamma(1000,1/2,7.5), probs = c(0.35,0.65))
# Initialize empty vectors for storing the sampled parameter values
sampled_alpha <- numeric(num_iterations)
sampled_beta_two <- numeric(num_iterations)
sampled_beta_one <- numeric(num_iterations)
sampled_lambda <- numeric(num_iterations)
sampled_sigma <- numeric(num_iterations)

# Initialize the starting values for alpha, beta, and sigma
current_alpha <- -30#-35#-30
current_beta_two <-0
current_beta_one <- 60#55#55
current_lambda<-10
current_sigma <- 3 #rgamma(1, shape = 1/2, rate = 1/2)

prior_lam <- function(gs) {###this is the prior of g, which is p(g) or p(lambda)
  #result<-gs^(-3/2)*exp(-length(x)/(2*gs))
  result<- dgamma(1/gs, 1/2, length(x)/2)##inverse gamma distribution
  return(result)
}


post_beta <- function(alpha,beta_two, beta_one, lambda, sigma){
  result<-dmvnorm(c(beta_one, beta_two, alpha),(lambda/(lambda+1))*c(b_hat),(lambda/(lambda+1))*(sigma^2)*coe)
  #result<-dmvnorm(c(beta_one, beta_two, alpha),(lambda/(lambda+1))*c_b,(lambda/(lambda+1))*(sigma^2)*coe)
  #result<-dmvnorm(c(alpha,beta_two,beta_one),(lambda/(lambda+1))*c(b_hat),(lambda/(lambda+1))*(sigma^2)*coe)
  return(result)
}

#post_beta<-function(alpha,beta,gs,sigma){ #this is prior p(beta, alpha | lambda)
#result<-dmvnorm(c(alpha, beta),c(0,0),gs*sigma*coe)
#return(result)

#}

#c(t(b_hat))
b_hat
#inv_bet <- c(-3.767022e+01,-3.202679e-04,6.117137e+01)
post_sigma <- function(sigma,lambda){
  #result<- dgamma(1/sigma^2, length(x)/2, (svar/2)+ (1/(2*(lambda+1))) * (c(b_hat)%*%cx %*% t(cx)%*%c(b_hat)))
  result<- dgamma(1/sigma^2, length(x)/2, (svar/2)+ (1/(2*(lambda+1))) * (inv_bet%*%cx %*% t(cx)%*%inv_bet))
  return(result)
}


see_ratio<-c()
see_p<-c()
see_cp<-c()
# Perform MCMC sampling
for (i in 1:num_iterations) {
  # Propose new values for alpha, beta, and sigma
  proposed_beta_one <- rnorm(1, mean = current_beta_one, sd =5)
  proposed_beta_two <- rnorm(1, mean = current_beta_two, sd =5)
  proposed_alpha <- rnorm(1, mean = current_alpha, sd = 5)
  proposed_lambda <- runif(1,min = 20, max = 70)
  #proposed_sigma <- rgamma(1, shape = 1/2, rate = 0.5)#1/rgamma(1, shape = length(x)/2, rate = (svar/2)+ (1/(2*(lambda+1))) * (c(b_hat)%*%cx %*% t(cx)%*%c(b_hat)))
  #proposed_sigma <- runif(1,min = 0, max=5)
  proposed_sigma <- rgamma(1, 2, 1.5)#
  
  
  log_likelihood_proposed <-prior_lam(proposed_lambda) * post_beta(proposed_alpha,proposed_beta_two, proposed_beta_one,proposed_lambda,proposed_sigma)*post_sigma(proposed_sigma,proposed_lambda)
  #see_p<-append(see_p,c(proposed_alpha,proposed_beta_two, proposed_beta_one,proposed_lambda,proposed_sigma))
  
  
  log_likelihood_current <- prior_lam(current_lambda) * post_beta(current_alpha,current_beta_two, current_beta_one,current_lambda,current_sigma)*post_sigma(current_sigma,current_lambda)
  #see_cp<-append(see_cp,c(current_alpha,current_beta_two, current_beta_one,current_lambda,current_sigma))
  
  
  flag = log_likelihood_proposed/log_likelihood_current
  
  # Calculate the acceptance ratio
  if (!is.nan(flag)){
    acceptance_ratio <- min(1,log_likelihood_proposed/log_likelihood_current)}
  #see_ratio<-append(see_ratio,acceptance_ratio)
  # Accept or reject the proposed values
  if (runif(1) < acceptance_ratio) {
    current_beta_one <- proposed_beta_one
    current_beta_two <- proposed_beta_two
    
    current_alpha <- proposed_alpha
    current_lambda <- proposed_lambda
    current_sigma <- proposed_sigma
  }
  sampled_beta_one[i] <- current_beta_one
  sampled_beta_two[i]<- current_beta_two
  sampled_alpha[i] <- current_alpha
  sampled_lambda[i] <- current_lambda
  sampled_sigma[i]<-current_sigma
  
  # Store the current values after burn-in
}
#end.time <- Sys.time()
#time.taken <- round(end.time - start.time,2)
#time.taken



# Print the posterior means
sampled_beta_one=sampled_beta_one[10000:num_iterations]
sampled_beta_two=sampled_beta_two[10000:num_iterations]
sampled_alpha=sampled_alpha[10000:num_iterations]
sampled_lambda=sampled_lambda[10000:num_iterations]
sampled_sigma=sampled_sigma[10000:num_iterations]

mbone = mean(sampled_beta_one)
mbtwo = mean(sampled_beta_two)
malp = mean(sampled_alpha)
#mean(sampled_lambda)
msig = mean(sampled_sigma)
sdsig = sd(sampled_sigma)






#max(sampled_lambda)
#min(sampled_lambda)



















##################################XIS+++++++++++++ISU+++++++++++++++++++++++
M=1500
xv = 1.8 #1.5，1.57，1.65， 1.8
#set.seed(214)
#set.seed(100)

mean_v <- c(mbone,mbtwo,malp)#c(-36.64903, 55.67207)
data <- data.frame(one = sampled_beta_one,two = sampled_beta_two,al= sampled_alpha)
cov_v <- cov(data)#matrix(c(2.252578^2, -0.9707999, -0.9707999, 1.247183^2), nrow = 2)
samples <- rmvnorm(M, mean = mean_v, sigma = cov_v)
sigm <- rnorm(M,msig,sdsig)
alpha<-samples[,3]
beta_t<-samples[,2]
beta_o<-samples[,1]
#sig<- 0.1142345 #mean(sampled_sigma)
hfunc <- function(lambda){(samples[,3] + 2.509275*samples[,2] + samples[,1]*xv)*(lambda/(1+lambda))}
#h_prior  <- function(lambda){dmvnorm(samples, mean = c(0,0), sigma = lambda*sig*coe)}
prop_den <- dmvnorm(samples,mean=mean_v, sigma = cov_v)
prop_sig <- dnorm(sigm,msig,sdsig)
#start.time <- Sys.time()




#######This is the f(Beta, sigma^2 | D, lambda)\
post_beta <- function(alpha,beta_two, beta_one, lambda, sigma){
  result<-dmvnorm(c(beta_one, beta_two, alpha),(lambda/(lambda+1))*c(b_hat),(lambda/(lambda+1))*(sigma^2)*coe)
  #result<-dmvnorm(c(beta_one, beta_two, alpha),(lambda/(lambda+1))*c_b,(lambda/(lambda+1))*(sigma^2)*coe)
  #result<-dmvnorm(c(alpha,beta_two,beta_one),(lambda/(lambda+1))*c(b_hat),(lambda/(lambda+1))*(sigma^2)*coe)
  return(result)
}
#post_beta<-function(alpha,beta,gs,sigma){ #this is prior p(beta, alpha | lambda)
#result<-dmvnorm(c(alpha, beta),c(0,0),gs*sigma*coe)
#return(result)

#}

#inv_bet <- c(-26.46931,-0.004570573,63.67152)
post_sigma <- function(sigma,lambda){
  result<- dgamma(1/sigma^2, length(x)/2, (svar/2)+ (1/(2*(lambda+1))) * (inv_bet%*%cx %*% t(cx)%*%inv_bet))
  return(result)
}




poster <- function(lambda, alp, beta_t,beta_o,sig) {
  
  #result<-dmvt(c(alp,bet), (lambda/(1+lambda))*c(b_hat), ((svar/2)+ (1/(2*(lambda+1))) * (c(b_hat)%*%cx %*% t(cx)%*%c(b_hat)))[1]*coe/7.5, 15, log=FALSE)
  result<-post_beta(alp, beta_t,beta_o, lambda, sig)*post_sigma(sig,lambda)
  return(result)
}


########



nuo = 0
numerator<-function(lambda){
  
  for (i in 1:M){
    
    nuo<-nuo+hfunc(lambda)*poster(lambda,alpha[i],beta_t[i], beta_o[i], sigm[i])/(prop_den*prop_sig)
    #nuo<-hfunc(lambda)*poster(lambda)/prop_den
    #nuo<-hfunc(lambda)*poster(lambda,alpha,beta,sigm)/(prop_den*prop_sig)
  }
  return(nuo)
}

deo=0
denominator<-function(lambda){
  for (i in 1:M){
    deo<-deo+poster(lambda,alpha[i],beta_t[i], beta_o[i], sigm[i])/(prop_den*prop_sig)
    #deo<-poster(lambda)/prop_den
    #deo<-poster(lambda,alpha,beta,sigm)/(prop_den*prop_sig)
  }
  return(deo)
}


#inf and sup of sampled_lambda for computing posterior expectation
llamb = as.integer(range(sampled_lambda)[1])
ulamb = as.integer(range(sampled_lambda)[2])

resapp<-c()
theo_v<-c()
grid <-seq(llamb, ulamb, length.out =15)#by = 1#lambda
weight<-c()
for(i in 1:15){
  res = sum(numerator(grid[i]))/sum(denominator(grid[i]))##this is the result of E(h(theta) | lambda, D)
  #res = numerator(grid[i])/denominator(grid[i])
  #we<-na.omit(denominator(grid[i]))
  #we_res<-max(we)/mean(we)
  #weight<-append(weight,we_res)
  resapp<-append(resapp,res)
}

mc_isu=na.omit(resapp)
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken
for(i in 1:15){
  ret = (grid[i]/(1+grid[i]))*(b_hat[1]*xv+b_hat[2]*2.509275+b_hat[3])#b_hat[2]*2.509275+)#*10+b_hat[3])
  theo_v <-append(theo_v,ret)
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
for(i in 1:15){
  res = sdcomp(grid[i])
  sdr<-append(sdr,res)
}

time_c <- append(time_c, time.taken)

inf_exp <- append(inf_exp, range(mc_isu)[1])
sup_exp <- append(sup_exp, range(mc_isu)[2])
msdr <- c(msdr,mean(sdr))

inf_sdr <- c(inf_sdr,range(sdr)[1])
sup_sdr<-c(sup_sdr,range(sdr)[2])

inf_the <-c(inf_the, range(theo_v)[1])
sup_the <- c(sup_the, range(theo_v)[2])
}



time_c 
inf_exp 
sup_exp 
msdr 
sup_sdr
inf_sdr
inf_the 
sup_the 





err = rnorm(length(sampled_alpha),0, 1.2)#sqrt(mean(sampled_sigma)))
h = sampled_beta_one*xv +sampled_beta_two*2.509275 +sampled_alpha+err
lambda = seq(20, 70, length.out =length(h))
range(loess.smooth(lambda,h)$y)
sd(loess.smooth(lambda,h)$y)


theo_v
mc_isu
