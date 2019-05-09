dual.pop <- '
#true scores
lx1 =~ 1*x1
lx2 =~ 1*x2
lx3 =~ 1*x3
lx4 =~ 1*x4

#autoregressions
lx2 ~ 1*lx1
lx3 ~ 1*lx2
lx4 ~ 1*lx3

#change scores
dx2 =~ 1*lx2
dx3 =~ 1*lx3
dx4 =~ 1*lx4

#proportional change
dx2 ~ 0.2*lx1
dx3 ~ 0.2*lx2
dx4 ~ 0.2*lx3

#intercept
i =~ 1*lx1
i ~ 10*1
i ~~ 1*i

#slope
s =~ 1*dx2 + 1*dx3 + 1*dx4
s ~ 1*1
s ~~ 1*s

#cov
i ~~ 0*s

#means
x1 ~ 0*1
x2 ~ 0*1
x3 ~ 0*1
x4 ~ 0*1

lx1 ~ 0*1
lx2 ~ 0*1
lx3 ~ 0*1
lx4 ~ 0*1

dx2 ~ 0*1
dx3 ~ 0*1
dx4 ~ 0*1

#variances
lx1 ~~ 0*lx1
lx2 ~~ 0*lx2
lx3 ~~ 0*lx3
lx4 ~~ 0*lx4

dx2 ~~ 0*dx2
dx3 ~~ 0*dx3
dx4 ~~ 0*dx4

x1 ~~ 1*x1
x2 ~~ 1*x2
x3 ~~ 1*x3
x4 ~~ 1*x4

#covariances
dx2 ~~ 0*dx3 + 0*dx4 
dx3 ~~ 0*dx4 

i ~~ 0*dx2 + 0*dx3 + 0*dx4 
'

data = lavaan::simulateData(dual.pop, sample.nobs=500)


dual.mod <- '
#true scores
lx1 =~ 1*x1
lx2 =~ 1*x2
lx3 =~ 1*x3
lx4 =~ 1*x4

#autoregressions
lx2 ~ 1*lx1
lx3 ~ 1*lx2
lx4 ~ 1*lx3

#change scores
dx2 =~ 1*lx2
dx3 =~ 1*lx3
dx4 =~ 1*lx4

#proportional change
dx2 ~ beta*lx1
dx3 ~ beta*lx2
dx4 ~ beta*lx3

#intercept
i =~ 1*lx1
i ~ 1
i ~~ i

#slope
s =~ 1*dx2 + 1*dx3 + 1*dx4
s ~ 1
s ~~ s

#cov
i ~~ s

#means
x1 ~ 0*1
x2 ~ 0*1
x3 ~ 0*1
x4 ~ 0*1

lx1 ~ 0*1
lx2 ~ 0*1
lx3 ~ 0*1
lx4 ~ 0*1

dx2 ~ 0*1
dx3 ~ 0*1
dx4 ~ 0*1

#variances
lx1 ~~ 0*lx1
lx2 ~~ 0*lx2
lx3 ~~ 0*lx3
lx4 ~~ 0*lx4

dx2 ~~ 0*dx2
dx3 ~~ 0*dx3
dx4 ~~ 0*dx4

x1 ~~ e*x1
x2 ~~ e*x2
x3 ~~ e*x3
x4 ~~ e*x4

#covariances
dx2 ~~ 0*dx3 + 0*dx4 
dx3 ~~ 0*dx4 

i ~~ 0*dx2 + 0*dx3 + 0*dx4 
'

fit = lavaan::lavaan(dual.mod, data=data)

summary(fit, fit.measures=T)

model = "
#Factor structure
p =~ L01*x1 + L02*x2 + L03*x3 + L04*x4 
r =~ L11*x1 + L12*x2 + L13*x3 + L14*x4
s =~ 1*x1 + 1*x2 + 1*x3 + 1*x4
#s =~ L21*x1 + L22*x2 + L23*x3 + L24*x4


#Means
x1 + x2 + x3 + x4 ~ 0*1
p ~ beta_p*1
r ~ 0*1
s ~ beta_s*1

#Variances
x1 ~~ sigma2_u*x1
x2 ~~ sigma2_u*x2
x3 ~~ sigma2_u*x3
x4 ~~ sigma2_u*x4

p ~~ p
r ~~ 0*r
s ~~ s

#Covariances
 p ~~ s

#Phantom variables
#pbeta_p =~ 0*x1
pbeta_r =~ 0*x1
#pbeta_s =~ 0*x1

#pbeta_p ~~ 0*pbeta_p
pbeta_r ~~ 0*pbeta_r
#pbeta_p ~~ 0*pbeta_r

#pbeta_p ~ beta_p*1
pbeta_r ~ beta_r*1
#pbeta_s ~ beta_s*1

#Constraints
L01 == exp(beta_r*0)
L02 == exp(beta_r*1)
L03 == exp(beta_r*2)
L04 == exp(beta_r*3)

L11 == 0*beta_p*exp(beta_r*0)
L12 == 1*beta_p*exp(beta_r*1)
L13 == 2*beta_p*exp(beta_r*2)
L14 == 3*beta_p*exp(beta_r*3)
"

model1 = "
#Factor structure
p =~ L01*x1 + L02*x2 + L03*x3 + L04*x4 
s =~ 1*x1 + 1*x2 + 1*x3 + 1*x4

#Means
x1 + x2 + x3 + x4 ~ 0*1
p ~ beta_p*1
s ~ beta_s*1

#Variances
x1 ~~ sigma2_u*x1
x2 ~~ sigma2_u*x2
x3 ~~ sigma2_u*x3
x4 ~~ sigma2_u*x4

p ~~ p
s ~~ s

#Covariances
p ~~ s

#Phantom variables
pbeta_r =~ 0*x1
pbeta_r ~~ 0*pbeta_r
pbeta_r ~ beta_r*1

#Constraints
L01 == exp(beta_r*0)
L02 == exp(beta_r*1)
L03 == exp(beta_r*2)
L04 == exp(beta_r*3)
"

fitexp = lavaan::lavaan(model, data=data)
summary(fitexp,fit=T)
summary(fit)
fitmeasures(fitexp)["rmsea"]
fit
fitexp

fitexp1 = lavaan::lavaan(model1, data=data)
summary(fitexp1,fit=T)
fitmeasures(fitexp1)["rmsea"]

# stan exponential
library(rstan)
dat2 <- list(
  N = nrow(data),
  X = data,
  t = 4)

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


# other parameterization



stan.exponential2 <-"
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
mu[i,tt] = FS[i,1] + FS[i,2]*(1+exp(-alpha*(tt-1)));
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

exponential.out2 <- stan(model_code=stan.exponential2,iter=200,
                        init=init.function,
                        control = list(max_treedepth=20,adapt_delta=0.8),
                        data = dat2,chains=1,cores=1,
                        pars=c("sigma","M","phi","alpha","log_lik"))
exponential.out2
pairs(exponential.out2,pars=c("sigma","M","phi","alpha"))
traceplot(exponential.out2,"alpha")
traceplot(exponential.out2,"M")

