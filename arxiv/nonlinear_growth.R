
dat <- list(
  "N" = 27,  
  "x" =
    c(1, 1.5, 1.5, 1.5, 2.5, 4, 5, 5, 7, 8, 8.5, 9, 9.5, 9.5, 10, 
      12, 12, 13, 13, 14.5, 15.5, 15.5, 16.5, 17, 22.5, 29, 31.5),
  "Y" =
    c(1.8, 1.85, 1.87, 1.77, 2.02, 2.27, 2.15, 2.26, 2.47, 2.19, 
      2.26, 2.4, 2.39, 2.41, 2.5, 2.32, 2.32, 2.43, 2.47, 2.56, 2.65, 
      2.47, 2.64, 2.56, 2.7, 2.72, 2.57))


stanmodel <- "
data {
int<lower=0> N; 
real x[N]; 
real Y[N]; 
} 
parameters {
real alpha; 
real beta;  
real<lower=.5,upper= 1> lambda; // orginal gamma in the JAGS example  
real<lower=0> tau; 
} 
transformed parameters {
real sigma; 
sigma <- 1 / sqrt(tau); 
} 
model {
real m[N];
for (i in 1:N) 
m[i] <- alpha - beta * pow(lambda, x[i]);

Y ~ normal(m, sigma); 

alpha ~ normal(0.0, 1000); 
beta ~ normal(0.0, 1000); 
lambda ~ uniform(.5, 1); 
tau ~ gamma(.0001, .0001); 
}
generated quantities{
real Y_mean[N]; 
real Y_pred[N]; 
for(i in 1:N){
# Posterior parameter distribution of the mean
Y_mean[i] <- alpha - beta * pow(lambda, x[i]);
# Posterior predictive distribution
Y_pred[i] <- normal_rng(Y_mean[i], sigma);   
}
}
"




library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
fit <- stan(model_code = stanmodel, 
            model_name = "GrowthCurve", 
            data = dat)

Y_mean <- extract(fit, "Y_mean")
Y_mean_cred <- apply(Y_mean$Y_mean, 2, quantile, c(0.05, 0.95))
Y_mean_mean <- apply(Y_mean$Y_mean, 2, mean)

Y_pred <- extract(fit, "Y_pred")
Y_pred_cred <- apply(Y_pred$Y_pred, 2, quantile, c(0.05, 0.95))
Y_pred_mean <- apply(Y_pred$Y_pred, 2, mean)

plot(dat$Y ~ dat$x, xlab="x", ylab="Y", 
     ylim=c(1.6, 2.8), main="Non-linear Growth Curve")
lines(dat$x, Y_mean_mean)
points(dat$x, Y_pred_mean, pch=19)
lines(dat$x, Y_mean_cred[1,], col=4)
lines(dat$x, Y_mean_cred[2,], col=4)
lines(dat$x, Y_pred_cred[1,], col=2)
lines(dat$x, Y_pred_cred[2,], col=2)
legend(x="bottomright", bty="n", lwd=2, lty=c(NA, NA, 1, 1,1),
       legend=c("observation", "prediction", "mean prediction",
                "90% mean cred. interval", "90% pred. cred. interval"),
       col=c(1,1,1,4,2),  pch=c(1, 19, NA, NA, NA))
fit