library(blavaan)

wisc <- read.table("C:/Users/RJacobucci/Documents/GitHub/EDM_Labs/2015/wisc4vpe.dat")
wisc <- read.table("C:/Users/jacobucc/Documents/GitHub/EDM_Labs/2015/wisc4vpe.dat")
wisc <- read.table("/Users/rjacobuc/Documents/GitHub/EDM_Labs/2015/wisc4vpe.dat")
names(wisc)<- c("V1","V2","V4","V6","P1","P2","P4", "P6", "Moeducat")


library(lavaan)



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
dV1 ~ beta*lV1; dV2 ~ beta*lV2; dV3 ~ beta*lV4;

#residuals equal
V1 ~~ resid*V1; V2 ~~ resid*V2; V4 ~~ resid*V4; V6 ~~ resid*V6;
"

fit.lds <- lavaan(lds_1, data=wisc)
summary(fit.lds,standardized=TRUE,fit=TRUE)





library(rstan)

X <- wisc[,c("V1","V2","V4","V6")]
data = list()
#data$alpha <- c(1,1,2,2)
X = as.matrix(X)
N <- nrow(X)

dat2 <- list(
  N = N,
  X = X,
  t = 4,
  alpha=c(1,1,2,2))


lcs.stan <-"
data{
int N; // sample size
int t; 
vector[t] alpha;
vector[t] X[N]; // data matrix of order [N,P]
}

parameters{
//vector[2] FS[N]; // factor scores, matrix of order [N,D]
matrix[N,2] FS;
cholesky_factor_corr[2] Rho; // correlation matrix between factors
real sigma; // variance for each variable
vector[2] M;
real pi;
}

transformed parameters{
matrix[N,t] mu;
matrix[N,t-1] d;
matrix[2,2] Sd_d;
Sd_d = multiply_lower_tri_self_transpose(Rho);


  mu[,1] = FS[,1];

  for (tt in 2:t){
    d[,tt-1] = pi*mu[,tt-1] + alpha[tt]*FS[,2];
    mu[,tt] = d[,tt-1]+mu[,tt-1];
  }

}

model{

M ~ normal(0,10);
Sd_d ~ lkj_corr_cholesky(2.0);
pi ~ normal(0,.5);
sigma ~ gamma(2,2);


X[,1] ~ normal(mu[,1],pow(sigma,0.5));

for (i in 1:N){
FS[i] ~ multi_normal(M, Sd_d);



for (tt in 2:t){
X[i,tt] ~ normal(mu[i,tt], pow(sigma,0.5));
}
}

}
"

system.time(lcs.out<-stan(model_code=lcs.stan,iter=200,
              data = dat2,chains=1,
              pars=c("sigma","M","Sd_d","pi")))

lcs.out



lcs.stan2 <-"
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

cov_matrix[2] Sd_d;
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

M ~ normal(0,10);
Sd_d ~ wishart(2,diag_matrix(rep_vector(1.0, 2)));
pi ~ normal(0,.5);
sigma ~ gamma(2,2);


X[,1] ~ normal(mu[,1],pow(sigma,0.5));

for (i in 1:N){
FS[i] ~ multi_normal(M, Sd_d);

for (tt in 2:t){
X[i,tt] ~ normal(mu[i,tt], pow(sigma,0.5));
}
}

}
"

lcs.out2 <- stan(model_code=lcs.stan2,iter=200,
              data = dat2,chains=1,
              pars=c("sigma","M","Sd_d","pi"))
# 209




lcs.stan3 <-"
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
M ~ normal(0,10);

pi ~ normal(0,.5);
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
//matrix[2,2] phi;
//phi = diag_matrix(L_sigma) * L_Omega * diag_matrix(L_sigma);
corr_matrix[2] Omega;
cov_matrix[2] phi;
Omega = L_Omega * L_Omega';
phi = quad_form_diag(Omega, L_sigma);
}"

system.time(lcs.out3 <- stan(model_code=lcs.stan3,iter=200,
                 data = dat2,chains=1,
                 pars=c("sigma","M","phi","pi")))

lcs.out3




lcs.stan4 <-"
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
M ~ normal(0,10);

pi ~ normal(0,.5);
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
phi = diag_pre_multiply(L_sigma, L_Omega) * diag_pre_multiply(L_sigma, L_Omega)';
}"

system.time(lcs.out4 <- stan(model_code=lcs.stan4,iter=200,
                 data = dat2,chains=1,
                 pars=c("sigma","M","phi","pi")))

lcs.out4


lcs.stan5 <-"
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
M ~ normal(0,10);

pi ~ normal(0,.5);
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
phi = diag_pre_multiply(L_sigma, L_Omega) * diag_pre_multiply(L_sigma, L_Omega)';
}"

system.time(lcs.out5 <- stan(model_code=lcs.stan5,iter=200,
                             control = list(max_treedepth=20),
                             data = dat2,chains=3,
                             pars=c("sigma","M","phi","pi")))

lcs.out5


