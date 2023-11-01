library(mvtnorm)
#This is dataset, and we should make x to be n*2 matrix which is (x,1) for matrix manipulation
x <- c(1.47, 1.50, 1.52, 1.55, 1.57, 1.60, 1.63, 1.65, 1.68, 1.70, 1.73, 1.75, 1.78, 1.80, 1.83)
onec <-  rep(1, times = 15)
y <- c(52.21, 53.12, 54.48, 55.84, 57.20, 58.57, 59.93, 61.29, 63.11, 64.47, 66.28, 68.10, 69.92, 72.19, 74.46)
cx <- rbind(x,onec)
coe <- solve(cx %*% t(cx))##this is (x^T*x)^-1
# Set the number of iterations and burn-in period for MCMC
num_iterations <- 15000#num_iterations <- 30000
burn_in <- 5000#burn_in <- 1000
b_hat <- coe%*% cx %*% y
svar <-t((y-t(cx)%*%b_hat))%*% ((y-t(cx)%*%b_hat))

set.seed(214)

# Initialize empty vectors for storing the sampled parameter values
sampled_alpha <- numeric(num_iterations - burn_in)
sampled_beta <- numeric(num_iterations - burn_in)
sampled_lambda <- numeric(num_iterations - burn_in)
sampled_sigma <- numeric(num_iterations - burn_in)

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



# Perform MCMC sampling
for (i in 1:num_iterations) {
  # Propose new values for alpha, beta, and sigma
  proposed_alpha <- rnorm(1, mean = current_alpha, sd =3)
  proposed_beta <- rnorm(1, mean = current_beta, sd = 3)
  proposed_lambda <- runif(1,min = 10, max = 50)
  #proposed_sigma <- rgamma(1, shape = 1/2, rate = 0.5)#1/rgamma(1, shape = length(x)/2, rate = (svar/2)+ (1/(2*(lambda+1))) * (c(b_hat)%*%cx %*% t(cx)%*%c(b_hat)))
  proposed_sigma <- runif(1,min = 0, max=5)
  
  
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
sampled_alpha=sampled_alpha[5000:num_iterations]
sampled_beta=sampled_beta[5000:num_iterations]
sampled_lambda=sampled_lambda[5000:num_iterations]
sampled_sigma=sampled_sigma[5000:num_iterations]

mean(sampled_alpha)
mean(sampled_beta)
mean(sampled_lambda)
mean(sampled_sigma)

sd(sampled_alpha)
sd(sampled_beta)
cor(sampled_alpha,sampled_beta)

max(sampled_lambda)
min(sampled_lambda)