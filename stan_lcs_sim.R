library(lavaan)
sim.mod <- "

#latent variables
lX1 =~ 1*X1
lX2 =~ 1*X2
lX3 =~ 1*X3
lX4 =~ 1*X4
lX5 =~ 1*X5
lX6 =~ 1*X6



#autoregressions
lX2 ~ 1*lX1; lX3 ~ 1*lX2; lX4 ~ 1*lX3; lX5 ~ 1*lX4; lX6 ~ 1*lX5

#change - delta; d
dX1 =~ 1*lX2; dX2 =~ 1*lX3; dX3 =~ 1*lX4; dX4 =~ 1*lX5; dX5 =~ 1*lX6

#intercept and slope
inV =~ 1*lX1;
slope =~ 1*dX1 + 1*dX2 + 1*dX3 + 1*dX4 + 1*dX5


#manifest means @0
X1 ~ 0*1; X2 ~0*1; X3 ~ 0*1; X4 ~ 0*1; X5 ~ 0*1; X6 ~ 0*1

#slope and intercept means
slope ~ 1*1;
inV ~ 1*1;

#Latent variances and covariance
slope ~~ 0.5*slope;
inV ~~ 0.5*inV;
slope ~~ 0.1*inV;

#means and vars @0
lX1 ~ 0*1; lX2 ~0*1; lX3 ~ 0*1; lX4 ~ 0*1;lX5 ~ 0*1; lX6 ~ 0*1
dX1 ~ 0*1; dX2 ~0*1; dX3 ~ 0*1; dX4 ~0*1; dX5 ~ 0*1

lX1 ~~ 0*lX1; lX2 ~~ 0*lX2; lX3 ~~ 0*lX3; lX4 ~~ 0*lX4;lX5 ~~ 0*lX5; lX6 ~~ 0*lX6
dX1 ~~ 0*dX1; dX2 ~~ 0*dX2; dX3 ~~ 0*dX3; dX4 ~~ 0*dX4; dX5 ~~ 0*dX5

#auVo-proportions
dX1 ~ 0*lX1; dX2 ~ 0*lX2; dX3 ~ .2*lX3;dX4 ~ .2*lX4;dX5 ~ .2*lX5;

#residuals equal
X1 ~~ 1*X1; X2 ~~ 1*X2; X3 ~~ 1*X3; X4 ~~ 1*X4;X5 ~~ 1*X5;X6 ~~ 1*X6;
"

set.seed(1)
dat <- simulateData(sim.mod,model.type="lavaan")


lds_1 <- "
#latent variables
lX1 =~ 1*X1
lX2 =~ 1*X2
lX3 =~ 1*X3
lX4 =~ 1*X4
lX5 =~ 1*X5
lX6 =~ 1*X6



#autoregressions
lX2 ~ 1*lX1; lX3 ~ 1*lX2; lX4 ~ 1*lX3; lX5 ~ 1*lX4; lX6 ~ 1*lX5

#change - delta; d
dX1 =~ 1*lX2; dX2 =~ 1*lX3; dX3 =~ 1*lX4; dX4 =~ 1*lX5; dX5 =~ 1*lX6

#intercept and slope
inV =~ 1*lX1;
slope =~ 1*dX1 + 1*dX2 + 1*dX3 + 1*dX4 + 1*dX5


#manifest means @0
X1 ~ 0*1; X2 ~0*1; X3 ~ 0*1; X4 ~ 0*1; X5 ~ 0*1; X6 ~ 0*1

#slope and intercept means
slope ~ 1;
inV ~ 1;

#Latent variances and covariance
slope ~~ slope;
inV ~~ inV;
slope ~~ inV;

#means and vars @0
lX1 ~ 0*1; lX2 ~0*1; lX3 ~ 0*1; lX4 ~ 0*1;lX5 ~ 0*1; lX6 ~ 0*1
dX1 ~ 0*1; dX2 ~0*1; dX3 ~ 0*1; dX4 ~0*1; dX5 ~ 0*1

lX1 ~~ 0*lX1; lX2 ~~ 0*lX2; lX3 ~~ 0*lX3; lX4 ~~ 0*lX4;lX5 ~~ 0*lX5; lX6 ~~ 0*lX6
dX1 ~~ 0*dX1; dX2 ~~ 0*dX2; dX3 ~~ 0*dX3; dX4 ~~ 0*dX4; dX5 ~~ 0*dX5

#auVo-proportions
dX1 ~ beta*lX1; dX2 ~ beta*lX2; dX3 ~ beta*lX3;dX4 ~ beta*lX4;dX5 ~ beta*lX5;

#residuals equal
X1 ~~ resid*X1; X2 ~~ resid*X2; X3 ~~ resid*X3; X4 ~~ resid*X4;X5 ~~ resid*X5;X6 ~~ resid*X6;
"

fit.lds <- lavaan(lds_1, data=dat)
summary(fit.lds)

inspect(fit.lds,"gradient")
modindices(fit.lds)
# stan

library(blavaan);library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


X <- dat
data = list()
#data$alpha <- c(1,1,2,2)
X = as.matrix(X)
N <- nrow(X)

dat2 <- list(
  N = N,
  X = X,
  t = 6)



stan.dual <-"
data{
int N; // sample size
int t; 
vector[t] X[N]; // data matrix of order [N,P]
}

parameters{
//vector[2] FS[N]; // factor scores, matrix of order [N,D]
matrix[N,2] FS;
real sigma; // variance for each variable
vector[2] M;
real pi;

cholesky_factor_corr[2] L_Omega;
vector<lower=0>[2] L_sigma;
}

transformed parameters{
matrix[N,t] mu;
matrix[N,t-1] d;


mu[,1] = FS[,1];

for (tt in 2:t){
d[,tt-1] = pi*mu[,tt-1] + FS[,2];
mu[,tt] = d[,tt-1]+mu[,tt-1];
}


}

model{
matrix[2,2] L_Sigma;
M[1] ~ normal(20,5);
M[2] ~ normal(-2,5);

pi ~ normal(0.4,.1);
sigma ~ gamma(2,2);

L_Sigma = diag_pre_multiply(L_sigma, L_Omega);
L_Omega ~ lkj_corr_cholesky(2);
L_sigma ~ cauchy(0, 2.5);


for (i in 1:N){
FS[i,] ~ multi_normal_cholesky(M, L_Sigma);

for (tt in 1:t){
X[i,tt] ~ normal(mu[i,tt], pow(sigma,0.5));
}
}
}
generated quantities {
cov_matrix[2] phi;
vector[N] log_lik;
phi = diag_pre_multiply(L_sigma, L_Omega) * diag_pre_multiply(L_sigma, L_Omega)';
for (n in 1:N)
log_lik[n] = normal_lpdf(X[n,]| mu[n,], sigma);
}"

#init.function <- function(){
 # list(pi = rnorm(1,0.4,.1))
#}

system.time(dual.out <- stan(model_code=stan.dual,iter=200,
                             #init=init.function,
                             control = list(max_treedepth=20,adapt_delta=0.8),
                             data = dat2,chains=3,cores=3,
                             pars=c("sigma","M","phi","pi","log_lik")))

pairs(dual.out,pars = c("M","pi","phi","sigma"))
dual.out

