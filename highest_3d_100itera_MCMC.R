
library(mvtnorm)
library(coda)#



set.seed(214)
#This is dataset, and we should make x to be n*2 matrix which is (x,1) for matrix manipulation
x <- c(1.47, 1.50, 1.52, 1.55, 1.57, 1.60, 1.63, 1.65, 1.68, 1.70, 1.73, 1.75, 1.78, 1.80, 1.83)
cal <- c(runif(15, 2,5))# artifical k/calories taken per day
onec <-  rep(1, times = 15)
y <- c(52.21, 53.12, 54.48, 55.84, 57.20, 58.57, 59.93, 61.29, 63.11, 64.47, 66.28, 68.10, 69.92, 72.19, 74.46)
cx <- rbind(x,cal,onec)
coe <- solve(cx %*% t(cx))##this is (x^T*x)^-1
# Set the number of iterations and burn-in period for MCMC
num_iterations <- 55000#20000#num_iterations <- 30000
#burn_in <- 0#5000#burn_in <- 1000
b_hat <- coe%*% cx %*% y
svar <-t((y-t(cx)%*%b_hat))%*% ((y-t(cx)%*%b_hat))
inv_bet <- c(b_hat[3],b_hat[2],b_hat[1])




time_mc <- c()
inf_mc_exp <- c()
sup_mc_exp <- c()
msdr_mc_exp <- c()
sup_sdr_mc <-c()
inf_sdr_mc <- c()


#inf_the <-c()
#sup_the <- c()
######

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
current_alpha <--30#-35#-30
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
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
#time.taken



# Print the posterior means
sampled_beta_one=sampled_beta_one[10000:num_iterations]
sampled_beta_two=sampled_beta_two[10000:num_iterations]
sampled_alpha=sampled_alpha[10000:num_iterations]
sampled_lambda=sampled_lambda[10000:num_iterations]
sampled_sigma=sampled_sigma[10000:num_iterations]

sampled_lambda/(1+sampled_lambda)



mean_sig = mean(sampled_sigma)
sup_lam = range(sampled_lambda)[1]
inf_lam = range(sampled_lambda)[2]

xv=1.8
err = rnorm(length(sampled_alpha),0, sqrt(mean_sig))#sqrt(mean(sampled_sigma)))
h_mc = (sampled_beta_one*xv +sampled_beta_two*2.509275 +sampled_alpha)#+err)#*sampled_lambda/(1+sampled_lambda)
lambda = seq(sup_lam, inf_lam, length.out =length(h))





theta_ess <- effectiveSize(as.mcmc(h_mc))


index = as.integer(seq(1, 45000, length.out = as.integer(theta_ess)))
y_sample = h_mc[index]
x_sample = lambda[index]


lmObject <- loess(y_sample ~ x_sample)

#newdat <- data.frame(lambda_s = c(xu, xl))

se.est <- predict(lmObject, se=TRUE)$se
####


###set a simple SD computing when lambda is consistent. 
#####count how many times non-converage of 100 iterations. 
###MC-ISU still could give answers, but MCMC may not. 



time_mc <- append(time_mc,time.taken)
inf_mc_exp <- append(inf_mc_exp,range(loess.smooth(x_sample,y_sample )$y)[1])
sup_mc_exp <- append(sup_mc_exp,range(loess.smooth(x_sample,y_sample )$y)[2])
msdr_mc_exp <- append(msdr_mc_exp,mean(se.est))
inf_sdr_mc <-append(inf_sdr_mc,range(se.est)[1])
sup_sdr_mc <- append(sup_sdr_mc,range(se.est)[2])
}


time_mc 
inf_mc_exp
sup_mc_exp
msdr_mc_exp 
inf_sdr_mc
sup_sdr_mc 




filename <- "rinf_sdr.txt"
write.table(rinf_sdr, file = filename)



data1 <- rsup_sdr
data2 <- sup_sdr_mc

# 将数据放入一个列表
data <- list(Data1 = data1, Data2 = data2)

# 创建箱线图
boxplot(data, 
        main = "Sup sdr, black line is sd(sup_mc_exp) Boxplots",  # 图的标题
        names = c("MC-ISU", "MCMC"),  # 每个箱线图的标签
        col = c("lightblue", "orange"),  # 箱线图的颜色
        border = "brown",  # 箱线的边框颜色
        horizontal = FALSE,  # 是否水平绘制
        notch = FALSE,  # 是否显示凹陷
        ylim = c(0, 1.2)  # Y 轴的范围
)
abline(h = sd(sup_mc_exp[is.finite(sup_mc_exp)]), col = 'black',lty = 1)
#70.62566
#68.22343
