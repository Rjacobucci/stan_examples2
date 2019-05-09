library(blavaan);library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

wisc <- read.table("C:/Users/RJacobucci/Documents/GitHub/EDM_Labs/2015/wisc4vpe.dat")
wisc <- read.table("C:/Users/jacobucc/Documents/GitHub/EDM_Labs/2015/wisc4vpe.dat")
wisc <- read.table("/Users/rjacobuc/Documents/GitHub/EDM_Labs/2015/wisc4vpe.dat")
names(wisc)<- c("V1","V2","V4","V6","P1","P2","P4", "P6", "Moeducat")


summary(wisc)
id = wisc[,1] < 5
wisc2 = wisc[!id,]


library(lavaan)

lds_linear <- "

#latent variables
lV1 =~ 1*V1
lV2 =~ 1*V2
lV4 =~ 1*V4
lV6 =~ 1*V6



#autoregressions
lV2 ~ 1*lV1; lV4 ~ 1*lV2; lV6 ~ 1*lV4

#change - delta; d
dV1 =~ 1*lV2; dV2 =~ 1*lV4; dV3 =~ 1*lV6

#intercept and slope
inV =~ 1*lV1;

# match lgm
slope =~ 1*dV1 + 2*dV2 + 2*dV3


#manifest means @0
V1 ~ 0*1; V2 ~0*1; V4 ~ 0*1; V6 ~ 0*1

#slope and intercept means
slope ~ 1;
inV ~ 1;

#Latent variances and covariance
slope ~~ slope;
inV ~~ inV;
slope ~~ inV;

#means and vars @0
lV1 ~ 0*1; lV2 ~0*1; lV4 ~ 0*1; lV6 ~ 0*1
dV1 ~ 0*1; dV2 ~0*1; dV3 ~ 0*1

lV1 ~~ 0*lV1; lV2 ~~ 0*lV2; lV4 ~~ 0*lV4; lV6 ~~ 0*lV6
dV1 ~~ 0*dV1; dV2 ~~ 0*dV2; dV3 ~~ 0*dV3

#residuals equal
V1 ~~ resid*V1; V2 ~~ resid*V2; V4 ~~ resid*V4; V6 ~~ resid*V6;
"

fit.linear <- lavaan(lds_linear, data=wisc2)
summary(fit.linear,standardized=TRUE,fit=TRUE)


lds_1 <- "

#latent variables
lV1 =~ 1*V1
lV2 =~ 1*V2
lV4 =~ 1*V4
lV6 =~ 1*V6



#autoregressions
lV2 ~ 1*lV1; lV4 ~ 1*lV2; lV6 ~ 1*lV4

#change - delta; d
dV1 =~ 1*lV2; dV2 =~ 1*lV4; dV3 =~ 1*lV6

#intercept and slope
inV =~ 1*lV1;

# match lgm
slope =~ 1*dV1 + 2*dV2 + 2*dV3


#manifest means @0
V1 ~ 0*1; V2 ~0*1; V4 ~ 0*1; V6 ~ 0*1

#slope and intercept means
slope ~ 1;
inV ~ 1;

#Latent variances and covariance
slope ~~ slope;
inV ~~ inV;
slope ~~ inV;

#means and vars @0
lV1 ~ 0*1; lV2 ~0*1; lV4 ~ 0*1; lV6 ~ 0*1
dV1 ~ 0*1; dV2 ~0*1; dV3 ~ 0*1

lV1 ~~ 0*lV1; lV2 ~~ 0*lV2; lV4 ~~ 0*lV4; lV6 ~~ 0*lV6
dV1 ~~ 0*dV1; dV2 ~~ 0*dV2; dV3 ~~ 0*dV3

#auVo-proportions
dV1 ~ beta*lV1 + start(.4)*lV1;
dV2 ~ beta*lV2 + start(.4)*lV2; 
dV3 ~ beta*lV4 + start(.4)*lV4;

#residuals equal
V1 ~~ resid*V1; V2 ~~ resid*V2; V4 ~~ resid*V4; V6 ~~ resid*V6;
"

fit.lds <- lavaan(lds_1, data=wisc2,optim.method="L-BFGS-B")
summary(fit.lds,standardized=TRUE,fit=TRUE)
fitmeasures(fit.lds)
lavInspect(fit.lds,"gradient") # didn't converge


# differences to differences

lds_2 <- "

#latent variables
lV1 =~ 1*V1
lV2 =~ 1*V2
lV4 =~ 1*V4
lV6 =~ 1*V6



#autoregressions
lV2 ~ 1*lV1; lV4 ~ 1*lV2; lV6 ~ 1*lV4

#change - delta; d
dV1 =~ 1*lV2; dV2 =~ 1*lV4; dV3 =~ 1*lV6

#intercept and slope
inV =~ 1*lV1;

# match lgm
slope =~ 1*dV1 + 2*dV2 + 2*dV3


#manifest means @0
V1 ~ 0*1; V2 ~0*1; V4 ~ 0*1; V6 ~ 0*1

#slope and intercept means
slope ~ 1;
inV ~ 1;

#Latent variances and covariance
slope ~~ slope;
inV ~~ inV;
slope ~~ inV;

#means and vars @0
lV1 ~ 0*1; lV2 ~0*1; lV4 ~ 0*1; lV6 ~ 0*1
dV1 ~ 0*1; dV2 ~0*1; dV3 ~ 0*1

lV1 ~~ 0*lV1; lV2 ~~ 0*lV2; lV4 ~~ 0*lV4; lV6 ~~ 0*lV6
dV1 ~~ 0*dV1; dV2 ~~ 0*dV2; dV3 ~~ 0*dV3

#auVo-proportions
dV2 ~ beta*dV1; dV3 ~ beta*dV2; 

#residuals equal
V1 ~~ resid*V1; V2 ~~ resid*V2; V4 ~~ resid*V4; V6 ~~ resid*V6;
"

fit.lds2 <- lavaan(lds_2, data=wisc2,optim.method="BFGS")
summary(fit.lds2,standardized=TRUE,fit=TRUE)

lavInspect(fit.lds2,"gradient")


library(rstan)

X <- wisc2[,c("V1","V2","V4","V6")]
data = list()
#data$alpha <- c(1,1,2,2)
X = as.matrix(X)
N <- nrow(X)

dat2 <- list(
  N = N,
  X = X,
  t = 4,
  alpha=c(1,1,2,2))

# proportional change


stan.prop <-"
data{
int N; // sample size
int t; 
vector[t] alpha;
vector[t] X[N]; // data matrix of order [N,P]
}

parameters{

real sigma; // variance for each variable
real pi;
vector[N] intt;

real intt_mean;
real intt_var;
}

transformed parameters{
matrix[N,t] mu;
matrix[N,t-1] d;


mu[,1] = intt;

for (tt in 2:t){
d[,tt-1] = pi*mu[,tt-1];
mu[,tt] = d[,tt-1]+mu[,tt-1];
}


}

model{

pi ~ normal(0.4,5);
sigma ~ gamma(2,2);
intt_var ~ gamma(2,2);
intt_mean ~ normal(20,5);



for (i in 1:N){
intt[i] ~ normal(intt_mean,intt_var);

for (tt in 1:t){
X[i,tt] ~ normal(mu[i,tt], pow(sigma,0.5));
}
}
}
generated quantities {
vector[N] log_lik;
for (n in 1:N)
log_lik[n] = normal_lpdf(X[n,]| mu[n,], sigma);
}"



system.time(prop.out <- stan(model_code=stan.prop,iter=200,
                             #init=init.function,
                             #control = list(max_treedepth=20,adapt_delta=0.99),
                             data = dat2,chains=3,cores=3,
                             pars=c("sigma","intt_mean","intt_var","pi","log_lik")))
pairs(lcs.out5)
pairs(prop.out,pars = c("pi","intt_mean","intt_var","sigma"))
prop.out
traceplot(lcs.out5,"M")



lik.prop <- extract_log_lik(prop.out)
loo.prop = loo(lik.prop)
loo.prop

# proportional with vector of pi


stan.prop2 <-"
data{
int N; // sample size
int t; 
vector[t] alpha;
vector[t] X[N]; // data matrix of order [N,P]
}

parameters{

real sigma; // variance for each variable
vector[t-1] pi;
vector[N] intt;

real intt_mean;
real intt_var;
}

transformed parameters{
matrix[N,t] mu;
matrix[N,t-1] d;


mu[,1] = intt;

for (tt in 2:t){
d[,tt-1] = pi[tt-1]*mu[,tt-1];
mu[,tt] = d[,tt-1]+mu[,tt-1];
}


}

model{

pi ~ normal(0.4,.01);
sigma ~ gamma(2,2);
intt_var ~ gamma(2,2);
intt_mean ~ normal(20,5);



for (i in 1:N){
intt[i] ~ normal(intt_mean,intt_var);

for (tt in 1:t){
X[i,tt] ~ normal(mu[i,tt], pow(sigma,0.5));
}
}
}
generated quantities {
vector[N] log_lik;
for (n in 1:N)
log_lik[n] = normal_lpdf(X[n,]| mu[n,], sigma);
}"



system.time(prop.out2 <- stan(model_code=stan.prop2,iter=200,
                             #init=init.function,
                             #control = list(max_treedepth=20,adapt_delta=0.99),
                             data = dat2,chains=3,cores=3,
                             pars=c("sigma","intt_mean","intt_var","pi","log_lik")))

pairs(prop.out2,pars = c("pi","intt_mean","intt_var","sigma"))
prop.out2


lik.prop2 <- extract_log_lik(prop.out2)
loo.prop = loo(lik.prop2)
loo.prop2

# linear change


stan.linear <-"
data{
int N; // sample size
int t; 
vector[t] alpha;
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
d[,tt-1] = alpha[tt]*FS[,2];
mu[,tt] = d[,tt-1]+mu[,tt-1];
}

}

model{
matrix[2,2] L_Sigma;
M[1] ~ normal(20,5);
M[2] ~ normal(-2,5);

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

system.time(linear.out <- stan(model_code=stan.linear,iter=200,
                             #init=init.function,
                             #control = list(max_treedepth=20),
                             data = dat2,chains=3,cores=3,
                             pars=c("sigma","M","phi","log_lik")))

linear.out

library(loo)
lik.linear <- extract_log_lik(linear.out)
loo.linear = loo(lik.linear)
loo.linear
# dual change


stan.dual <-"
data{
int N; // sample size
int t; 
vector[t] alpha;
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
d[,tt-1] = pi*mu[,tt-1] + alpha[tt]*FS[,2];
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

init.function <- function(){
  list(pi = rnorm(1,0.4,.1))
}

system.time(dual.out <- stan(model_code=stan.dual,iter=200,
                             #init=init.function,
                             control = list(max_treedepth=20,adapt_delta=0.99),
                             data = dat2,chains=3,cores=3,
                             pars=c("sigma","M","phi","pi","log_lik")))
pairs(lcs.out5)
pairs(dual.out,pars = c("M","pi","phi","sigma"))
lcs.out5
traceplot(lcs.out5,"M")



lik.dual <- extract_log_lik(dual.out)
loo.dual = loo(lik.dual)
loo.dual


# dual with vector of pi


stan.dual2 <-"
data{
int N; // sample size
int t; 
vector[t] alpha;
vector[t] X[N]; // data matrix of order [N,P]
}

parameters{
//vector[2] FS[N]; // factor scores, matrix of order [N,D]
matrix[N,2] FS;
real sigma; // variance for each variable
vector[2] M;
vector[t-1] pi;

cholesky_factor_corr[2] L_Omega;
vector<lower=0>[2] L_sigma;
}

transformed parameters{
matrix[N,t] mu;
matrix[N,t-1] d;


mu[,1] = FS[,1];

for (tt in 2:t){
d[,tt-1] = pi[tt-1]*mu[,tt-1] + alpha[tt]*FS[,2];
mu[,tt] = d[,tt-1]+mu[,tt-1];
}


}

model{
matrix[2,2] L_Sigma;
M[1] ~ normal(20,5);
M[2] ~ normal(-2,5);

pi ~ normal(0.4,.01);
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

init.function <- function(){
  list(pi = rnorm(1,0.4,.1))
}

system.time(dual.out2 <- stan(model_code=stan.dual2,iter=200,
                             #init=init.function,
                             #control = list(max_treedepth=20,adapt_delta=0.99),
                             data = dat2,chains=3,cores=3,
                             pars=c("sigma","M","phi","pi","log_lik")))

pairs(dual.out2,pars = c("M","pi","phi","sigma"))
dual.out2
traceplot(dual.out2,"pi")


lik.dual2 <- extract_log_lik(dual.out2)
loo.dual2 = loo(lik.dual2)
loo.dual2


# estimate mean of pi


stan.dual3 <-"
data{
int N; // sample size
int t; 
vector[t] alpha;
vector[t] X[N]; // data matrix of order [N,P]
}

parameters{
//vector[2] FS[N]; // factor scores, matrix of order [N,D]
matrix[N,2] FS;
real sigma; // variance for each variable
vector[2] M;
vector[t-1] pi;
real pi_mean;

cholesky_factor_corr[2] L_Omega;
vector<lower=0>[2] L_sigma;
}

transformed parameters{
matrix[N,t] mu;
matrix[N,t-1] d;


mu[,1] = FS[,1];

for (tt in 2:t){
d[,tt-1] = pi[tt-1]*mu[,tt-1] + alpha[tt]*FS[,2];
mu[,tt] = d[,tt-1]+mu[,tt-1];
}


}

model{
matrix[2,2] L_Sigma;
M[1] ~ normal(20,5);
M[2] ~ normal(-2,5);

pi ~ normal(pi_mean,.01);
pi_mean ~ normal(0.3,.1);
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

init.function <- function(){
  list(pi = rnorm(1,0.4,.1))
}

system.time(dual.out3 <- stan(model_code=stan.dual3,iter=200,
                              #init=init.function,
                              #control = list(max_treedepth=20,adapt_delta=0.99),
                              data = dat2,chains=3,cores=3,
                              pars=c("sigma","M","phi","pi","log_lik")))

pairs(dual.out3,pars = c("M","pi","phi","sigma"))
dual.out3
traceplot(dual.out3,"pi")


lik.dual3 <- extract_log_lik(dual.out3)
loo.dual3 = loo(lik.dual3)
loo.dual3


# estimate variance


stan.dual4 <-"
data{
int N; // sample size
int t; 
vector[t] alpha;
vector[t] X[N]; // data matrix of order [N,P]
}

parameters{
//vector[2] FS[N]; // factor scores, matrix of order [N,D]
matrix[N,2] FS;
real sigma; // variance for each variable
vector[2] M;
vector[t-1] pi;
real pi_mean;
real pi_var;

cholesky_factor_corr[2] L_Omega;
vector<lower=0>[2] L_sigma;
}

transformed parameters{
matrix[N,t] mu;
matrix[N,t-1] d;


mu[,1] = FS[,1];

for (tt in 2:t){
d[,tt-1] = pi[tt-1]*mu[,tt-1] + alpha[tt]*FS[,2];
mu[,tt] = d[,tt-1]+mu[,tt-1];
}


}

model{
matrix[2,2] L_Sigma;
M[1] ~ normal(20,5);
M[2] ~ normal(-2,5);

pi ~ normal(pi_mean,pi_var);
pi_mean ~ normal(0,1);
sigma ~ gamma(2,2);
pi_var ~ gamma(.1,10);

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

init.function <- function(){
  list(pi = rnorm(3,0.4,.1))
}

system.time(dual.out4 <- stan(model_code=stan.dual4,iter=200,
                              init=init.function,
                              control = list(max_treedepth=20,adapt_delta=0.99),
                              data = dat2,chains=3,cores=3,
                              pars=c("sigma","M","phi","pi","pi_mean","pi_var","log_lik")))

pairs(dual.out4,pars = c("M","pi","phi","sigma"))
dual.out4

traceplot(dual.out4,"pi_var")

lik.dual4 <- extract_log_lik(dual.out4)
loo.dual4 = loo(lik.dual4)
loo.dual4

# try changes to changes


stan.changes <-"
data{
int N; // sample size
int t; 
vector[t] alpha;
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
d[,tt-1] = alpha[tt]*FS[,2];
mu[,tt] = d[,tt-1]+mu[,tt-1];
}

for (tt in 3:t){
d[,tt-1] = pi*d[,tt-2] + alpha[tt]*FS[,2];
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

init.function <- function(){
  list(pi = rnorm(1,0.4,.1))
}

system.time(changes.out <- stan(model_code=stan.changes,iter=200,
                             #init=init.function,
                             #control = list(max_treedepth=20,adapt_delta=0.99),
                             data = dat2,chains=3,cores=3,
                             pars=c("sigma","M","phi","pi","log_lik")))
pairs(changes.out)
pairs(changes.out,pars = c("M","pi","phi","sigma"))
changes.out
traceplot(changes.out,"M")



lik.changes <- extract_log_lik(changes.out)
loo.changes = loo(lik.changes)
loo.changes


# exponential

stan.exponential <-"
data{
int N; // sample size
int t; 
vector[t] X[N]; // data matrix of order [N,P]
}

parameters{
matrix[N,2] FS;
real sigma; // variance for each variable
vector[2] M;
real alpha;

cholesky_factor_corr[2] L_Omega;
vector<lower=0>[2] L_sigma;
}

transformed parameters{
matrix[N,t] mu;

for(i in 1:N){
for (tt in 1:t){
mu[i,tt] = FS[i,1] + FS[i,2]*(exp(alpha*(tt-1)));
}
}

}

model{
matrix[2,2] L_Sigma;
M[1] ~ normal(20,100);
M[2] ~ normal(0,100);

alpha ~ normal(0,2);
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

init.function <- function(){
  list(M = rnorm(2,c(15,-2),1))
}

exponential.out <- stan(model_code=stan.exponential,iter=200,
                        #init=init.function,
                        control = list(max_treedepth=20,adapt_delta=0.8),
                        data = dat2,chains=1,cores=1,
                        pars=c("sigma","M","phi","alpha","log_lik"))
exponential.out
pairs(exponential.out,pars=c("sigma","M","phi","alpha"))
traceplot(exponential.out,"alpha")
traceplot(exponential.out,"M")

