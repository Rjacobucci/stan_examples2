library(blavaan)

wisc <- read.table("C:/Users/RJacobucci/Documents/GitHub/EDM_Labs/2015/wisc4vpe.dat")
wisc <- read.table("C:/Users/jacobucc/Documents/GitHub/EDM_Labs/2015/wisc4vpe.dat")
names(wisc)<- c("V1","V2","V4","V6","P1","P2","P4", "P6", "Moeducat")




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
vector[2] FS[N]; // factor scores, matrix of order [N,D]

cholesky_factor_cov[2] Rho; // correlation matrix between factors
real sigma; // variance for each variable
vector[2] M;
real pi;
}

transformed parameters{
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
"

lcs.out=stan(model_code=lcs.stan,iter=5000,warmup=2000,
              data = dat2,chains=1,
              pars=c("sigma","M","Sd_d","pi"))

lcs.out
