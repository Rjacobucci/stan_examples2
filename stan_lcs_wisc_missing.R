library(blavaan)

wisc <- read.table("C:/Users/RJacobucci/Documents/GitHub/EDM_Labs/2015/wisc4vpe.dat")
wisc <- read.table("C:/Users/jacobucc/Documents/GitHub/EDM_Labs/2015/wisc4vpe.dat")
names(wisc)<- c("V1","V2","V4","V6","P1","P2","P4", "P6", "Moeducat")




dat <- as.matrix(wisc[,c("V1","V2","V4","V6")])
nmiss <- 50
miss <- sample.int(prod(dim(dat)), size = nmiss)
dat[miss] <- NA




dat_complete <- dat[!is.na(dat)]

ind_pres <- which(!is.na(dat), arr.ind = TRUE)
ind_miss <- which(is.na(dat), arr.ind = TRUE)


mod.data <- list(N = nrow(dat),
                 t = ncol(dat),
                 Ncomp = length(dat_complete),
                 Nmiss = sum(is.na(dat)),
                 dat_complete = dat_complete,
                 ind_pres = ind_pres,
                 ind_miss = ind_miss,
                 alpha=c(1,1,2,2))





lcs.stan <-"
data{
int N; // sample size
int t; 
vector[t] alpha;
int<lower=0> Ncomp; // Number of non-missing values
int<lower=0> Nmiss; // Number of missing values
real dat_complete[Ncomp];   // Vector of non-missing values
int ind_pres[Ncomp, 2];     // Matrix (row, col) of non-missing value indices
int ind_miss[Nmiss, 2];     // Matrix (row, col) of missing value indices
}

parameters{
vector[2] FS[N]; // factor scores, matrix of order [N,D]

cholesky_factor_cov[2] Rho; // correlation matrix between factors
real sigma; // variance for each variable
vector[2] M;
real pi;
// Vector containing stochastic nodes (for filling missing values
    real Xmiss[Nmiss]; 
}

transformed parameters{
vector[t] X[N];   // The data with interpolated missing values
vector[t] mu[N];
vector[t-1] d[N];
matrix[2,2] Sd_d;
Sd_d = multiply_lower_tri_self_transpose(Rho);


for (i in 1:N){
  mu[i,1] = FS[i,1];

  for (tt in 2:t){
    d[i,tt-1] = pi*mu[i,tt-1] + alpha[tt]*FS[i,2];
    mu[i,tt] = d[i,tt-1]+mu[i,tt-1];
  }
}

for(n in 1:Ncomp) {
        X[ind_pres[n,1]][ind_pres[n,2]] = dat_complete[n];
}

for(n in 1:Nmiss){
X[ind_miss[n,1]][ind_miss[n,2]] = Xmiss[n];
}

}

model{
//M ~ normal(0,10);
//Rho ~ lkj_corr_cholesky(2.0);
//pi ~ normal(0,.5);
//sigma ~ gamma(2,2);


for (i in 1:N){
  FS[i] ~ multi_normal(M, Sd_d);
  X[i,1] ~ normal(mu[i,1],pow(sigma,0.5));
  

  for (tt in 2:t){
    X[i,tt] ~ normal(mu[i,tt], pow(sigma,0.5));
  }
}
}
generated quantities {
  vector[N] log_lik;
for (n in 1:N){
log_lik[n] = normal_lpdf(X[n,]| mu[n], sigma); 
}
}
"

lcs.out=stan(model_code=lcs.stan,iter=5000,warmup=2000,
              data = mod.data,chains=1,
              pars=c("sigma","M","Sd_d","pi"))

lcs.out
