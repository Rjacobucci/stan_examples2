library(rstan);library(lavaan);library(loo)

#install_github('nathanvan/rstanmulticore')
library(rstanmulticore)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

HS =data.frame(scale(HolzingerSwineford1939[,7:15]))
D <-3
P =9
N <-301



HS.model <- '
f1 =~ NA*x1 + x2 + x3 
f2 =~ NA*x4 + x5 + x6
f3 =~ NA*x7 + x8 +  x9
f1~~1*f1
f2~~1*f2
f3~~1*f3
'


fit.HS1 = cfa(HS.model,HS,orthogonal=T,meanstructure=TRUE) # do.fit=F
summary(fit.HS1,fit.measures=T)


fa.data <-list(P=P,N=N,X=as.matrix(HS),D=D)

mod.stan <-"
data{
  int N; // sample size
  int P; // number of variables
  int D; // number of dimensions
  vector[P] X[N]; // data matrix of order [N,P]
}

parameters{
  vector[P] b; // intercepts
  vector<lower=0>[P] lam; // factor loadings
  vector[D] FS[N]; // factor scores, matrix of order [N,D]
  corr_matrix[D] Rho; // correlation matrix between factors
  vector<lower=0,upper=100>[P] var_p; // variance for each variable
}

transformed parameters{
  vector[D] M;
  vector<lower=0, upper=1000>[D] Sd_d; // sd of factors
  vector[P] mu[N];
  matrix[D,D] Ld;
  
  for (m in 1:D) {
    M[m] =0;
    Sd_d[m] =1;}
  
  Ld =diag_matrix(Sd_d) * cholesky_decompose(Rho);
  
  for(i in 1:N){
    mu[i,1] =b[1] + lam[1]*FS[i,1];
    mu[i,2] =b[2] + lam[2]*FS[i,1];
    mu[i,3] =b[3] + lam[3]*FS[i,1];
    mu[i,4] =b[4] + lam[4]*FS[i,2];
    mu[i,5] =b[5] + lam[5]*FS[i,2];
    mu[i,6] =b[6] + lam[6]*FS[i,2];
    mu[i,7] =b[7] + lam[7]*FS[i,3];
    mu[i,8] =b[8] + lam[8]*FS[i,3];
    mu[i,9] =b[9] + lam[9]*FS[i,3];
  }
}

model{
  
  b ~ normal(0, 100);
  lam ~ normal(0.4,var_p);
  
  for(i in 1:N){    
    X[i] ~ normal(mu[i],var_p);
    FS[i] ~ multi_normal_cholesky(M, Ld);
  }
  
}
generated quantities {
  vector[N] log_lik;
  for (n in 1:N){
    log_lik[n] = normal_lpdf(X[n,]| mu[n], var_p); 
  }
}

"

fa.model=stan(model_code=mod.stan,
                data = fa.data,chains=1,
              pars=c("lam","b","var_p","log_lik"))

mod <- stan_model(model_code=mod.stan)

fa.vb <- vb(mod,fa.data,pars=c("lam","b","var_p"),tol_rel_obj=.001)
fa.vb

fa.model

lik <- extract_log_lik(fa.model)
loo(lik)#;waic(lik)
