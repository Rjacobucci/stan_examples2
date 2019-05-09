library(rstan)



wisc <- read.table("C:/Users/RJacobucci/Documents/GitHub/EDM_Labs/2015/wisc4vpe.dat")
wisc <- read.table("C:/Users/jacobucc/Documents/GitHub/EDM_Labs/2015/wisc4vpe.dat")
wisc <- read.table("C:/Users/Ross/Documents/GitHub/EDM_Labs/2015/wisc4vpe.dat")
names(wisc)<- c("V1","V2","V4","V6","P1","P2","P4", "P6", "Moeducat")







library(lavaan)
library(semPlot)

model <- ' i0_1 =~ 1*V1 + 1*V2 + 1*V4 + 1*V6
s0_1 =~ 1*V1 + 2*V2 + 4*V4 + 6*V6 
#residuals equal
#V1 ~~ resid*V1; V2 ~~ resid*V2; V4 ~~ resid*V4; V6 ~~ resid*V6;
'
fit <- growth(model, data=wisc)
summary(fit)


model <- ' i0_1 =~ 1*V1 + 1*V2 + 1*V4 + 1*V6
s0_1 =~ 1*V1 + 2*V2 + 4*V4 + 6*V6 
#residuals equal
V1 ~~ resid*V1; V2 ~~ resid*V2; V4 ~~ resid*V4; V6 ~~ resid*V6;
i0_1 ~0
s0_1 ~ 0
i0_1~~1*i0_1;
'
fit <- growth(model, data=wisc)
summary(fit)





X <- wisc[,c("V1","V2","V4","V6")]
data = list()
#data$alpha <- c(1,1,2,2)
X = as.matrix(X)
N <- nrow(X)

dat2 <- list(
  N = N,
  X = X,
  t = 4,
  tt = rep(0,4))


mod.stan <-"
data{
int N; // sample size
int t; 
vector[4] tt;
vector[t] X[N]; // data matrix of order [N,P]
}

parameters{
vector[2] FS[N]; // factor scores, matrix of order [N,D]
cholesky_factor_cov[2] Rho; // correlation matrix between factors
cholesky_factor_cov[4] Rho2; // correlation matrix between factors
vector[2] M;
//vector[4] sig;
//real alpha;


//vector[4] lv_sigma[N]; 
}

transformed parameters{

vector[t] mu[N];
matrix[2,2] Sd_d;
matrix[4,4] Sd_d2;
Sd_d = multiply_lower_tri_self_transpose(Rho);
Sd_d2 = multiply_lower_tri_self_transpose(Rho2);



for(i in 1:N){
mu[i,1] = FS[i,1] + 1*FS[i,2] + sigma[i,1];
mu[i,2] = FS[i,1] + 2*FS[i,2] + sigma[i,2]; 
mu[i,3] = FS[i,1] + 4*FS[i,2] + sigma[i,3];
mu[i,4] = FS[i,1] + 6*FS[i,2] + sigma[i,4];

//sigma[i,1] = lv_sigma[i,1];
//sigma[i,2] = lv_sigma[i,2];
//sigma[i,3] = lv_sigma[i,3];
//sigma[i,4] = lv_sigma[i,4];

//for(j in 1:4){
//sigma[i,1] = alpha*init_sigma[i,1];
//sigma[i,2] = alpha*sigma[i,1];
//sigma[i,3] = alpha*sigma[i,2];
//}
}

}

model{
M[1] ~ normal(20,10);
M[2] ~ normal(5,10);
Rho ~ lkj_corr_cholesky(2.0);
Rho2 ~ lkj_corr_cholesky(5.0);
//alpha ~ normal(0,1);


for(i in 1:N){    
for(j in 1:4){
X[i,j] ~ normal(mu[i,j],.001); // fix variance to .0000001
}

sigma[i] ~ multi_normal(tt,Sd_d2);

FS[i] ~ multi_normal(M, Sd_d);
}
}
"

fa.model=stan(model_code=mod.stan,iter=10000,warmup=2000,
              data = dat2,chains=1,
              pars=c("M","Sd_d","Sd_d2"))

fa.model
