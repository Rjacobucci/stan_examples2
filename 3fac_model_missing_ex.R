library(rstan);library(lavaan)


# https://gist.github.com/ashiklom/d790474834ef1c67fde6466728ea71f4



dat =as.matrix(data.frame(scale(HolzingerSwineford1939[,7:15])))
nmiss <- 100
miss <- sample.int(prod(dim(dat)), size = nmiss)
dat[miss] <- NA


HS.model <- '
f1 =~ NA*x1 + x2 + x3 
f2 =~ NA*x4 + x5 + x6
f3 =~ NA*x7 + x8 +  x9
f1~~1*f1
f2~~1*f2
f3~~1*f3
'
fit.HS1 = cfa(HS.model,dat,orthogonal=T,meanstructure=TRUE,missing="fiml") # do.fit=F
summary(fit.HS1,fit.measures=T)





dat_complete <- dat[!is.na(dat)]

ind_pres <- which(!is.na(dat), arr.ind = TRUE)
ind_miss <- which(is.na(dat), arr.ind = TRUE)


mod.data <- list(Nrow = nrow(dat),
                 Ncol = ncol(dat),
                 Ncomp = length(dat_complete),
                 Nmiss = sum(is.na(dat)),
                 dat_complete = dat_complete,
                 ind_pres = ind_pres,
                 ind_miss = ind_miss,
                 D =3)








mod.stan <-"
data{
   int<lower=0> Nrow;
    int<lower=0> Ncol;
    int<lower=0> Ncomp; // Number of non-missing values
    int<lower=0> Nmiss; // Number of missing values
    real dat_complete[Ncomp];   // Vector of non-missing values
    int ind_pres[Ncomp, 2];     // Matrix (row, col) of non-missing value indices
    int ind_miss[Nmiss, 2];     // Matrix (row, col) of missing value indices
    int D; // number of dimensions
}

parameters{
  vector[Ncol] b; // intercepts
  vector<lower=0>[Ncol] lam; // factor loadings
  vector[D] FS[Nrow]; // factor scores, matrix of order [N,D]
  corr_matrix[D] Rho; // correlation matrix between factors
  vector<lower=0,upper=100>[Ncol] var_p; // variance for each variable


  // Vector containing stochastic nodes (for filling missing values
    real Xmiss[Nmiss]; 
}

transformed parameters{

  vector[Ncol] X[Nrow];   // The data with interpolated missing values

 
  vector<lower=0, upper=1000>[D] Sd_d; // sd of factors
  vector[D] M;
  vector[Ncol] mu[Nrow];
  matrix[D,D] Ld;
  
  for (m in 1:D) {
    M[m] =0;
    Sd_d[m] =1;}
  
  Ld =diag_matrix(Sd_d) * cholesky_decompose(Rho);
  
  for(i in 1:Nrow){
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

 for(n in 1:Ncomp) {
        X[ind_pres[n,1]][ind_pres[n,2]] = dat_complete[n];
}

for(n in 1:Nmiss){
X[ind_miss[n,1]][ind_miss[n,2]] = Xmiss[n];
}
}

model{
  
  b ~ normal(0, 100);
  lam ~ normal(0.4,var_p);
  
  for(i in 1:Nrow){    
    X[i] ~ normal(mu[i],var_p);
    FS[i] ~ multi_normal_cholesky(M, Ld);
  }
  
}
generated quantities {
  vector[Nrow] log_lik;
for (n in 1:Nrow){
log_lik[n] = normal_lpdf(X[n,]| mu[n], var_p); 
}
}
"

fa.model=stan(model_code=mod.stan,iter=5000,warmup=2000,
                data = mod.data,chains=1,
              pars=c("lam","b","var_p","Sd_d","M"))
print(fa.model)

lik <- extract_log_lik(fa.model)
loo(lik)#;waic(lik)
