library(rstan); library(lavaan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
HS <- HolzingerSwineford1939[complete.cases(HolzingerSwineford1939),]


mod <- "
f1 =~ NA*x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9
f1 ~ sex + grade + ageyr
f1~~1*f1
"
out <- sem(mod, HS,meanstructure=T)
summary(out,rsq=T)

X <- HS[,7:15]
cov <- HS[,c(2,3,6)]


X = as.matrix(X)
cov= as.matrix(cov)
N <- nrow(X)

dat <- list(
  N = N,
  X = X,
  cov = cov)


mod.stan <-"
data{
int N; // sample size
matrix[N,9] X; // data matrix of order [N,P]
matrix[N,3] cov; // data matrix of order [N,P]
}

parameters{
vector[N] FS; // factor scores, matrix of order [N,D]
vector<lower=0>[9] sigma;
cholesky_factor_cov[3] rho_cov; // correlation matrix between factors
vector[8] lam;
vector[9] alpha;
vector[3] beta;
real<lower=0> psi;
}

transformed parameters{

vector[9] mu[N];
vector[N] mu2;
matrix[3,3] sd_cov;
sd_cov = multiply_lower_tri_self_transpose(rho_cov);


for(i in 1:N){
mu[i,1] = alpha[1] + 1*FS[i];
mu[i,2] = alpha[2] + lam[1]*FS[i];
mu[i,3] = alpha[3] + lam[2]*FS[i];
mu[i,4] = alpha[4] + lam[3]*FS[i];
mu[i,5] = alpha[5] + lam[4]*FS[i];
mu[i,6] = alpha[6] + lam[5]*FS[i];
mu[i,7] = alpha[7] + lam[6]*FS[i];
mu[i,8] = alpha[8] + lam[7]*FS[i];
mu[i,9] = alpha[9] + lam[8]*FS[i];

mu2[i] =  beta[1]*cov[i,1] + beta[2]*cov[i,2] + beta[3]*cov[i,3];

}
}

model{
sigma ~ gamma(2,2);
//beta ~ normal(0,1);
alpha ~ normal(0,1);
psi ~ gamma(2,2);
rho_cov ~ lkj_corr_cholesky(2.0);

for(i in 1:N){  
  for(j in 1:9){
    X[i,j] ~ normal(mu[i,j],pow(sigma[j],0.5));
  }
FS[i] ~ normal(mu2[i],psi);
cov[i,] ~ multi_normal(rep_vector(0,3), sd_cov);
}
}
generated quantities{
// can use mean() and sd()
real rsq2;
vector[N] z_mu;
vector[N] z_FS;
z_mu = (mu2 - mean(mu2))/sd(mu2);
z_FS = (FS - mean(FS))/sd(FS);
rsq2 = pow((z_mu' * z_FS)/N,2);
}
"
library(rstan)

fa.model=stan(model_code=mod.stan,
              data = dat,chains=3,iter=800,
              pars=c("rsq2", "sigma","lam","beta","psi","alpha"))

print(fa.model)

plot(fa.model,plotfun="trace",pars="rsq2")

tidy_stan(fa.model)
