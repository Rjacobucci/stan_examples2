library(rstan)



wisc <- read.table("C:/Users/RJacobucci/Documents/GitHub/EDM_Labs/2015/wisc4vpe.dat")
wisc <- read.table("C:/Users/jacobucc/Documents/GitHub/EDM_Labs/2015/wisc4vpe.dat")
names(wisc)<- c("V1","V2","V4","V6","P1","P2","P4", "P6", "Moeducat")







library(lavaan)
library(semPlot)

model <- ' i0_1 =~ 1*V1 + 1*V2 + 1*V4 + 1*V6
s0_1 =~ 1*V1 + 2*V2 + 4*V4 + 6*V6 
#residuals equal
V1 ~~ resid*V1; V2 ~~ resid*V2; V4 ~~ resid*V4; V6 ~~ resid*V6;'
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
  t = 4)


mod.stan <-"
data{
int N; // sample size
int t; 
vector[t] X[N]; // data matrix of order [N,P]
}

parameters{
real sigma; // variance for each variable
vector[2] beta;
  vector<lower=0>[2] sigma_u;
  cholesky_factor_corr[2] L_u;
  matrix[2,N] z_u;
}

transformed parameters{
vector[t] mu[N];

  matrix[2,N] FS;
  FS = diag_pre_multiply(sigma_u,L_u) * z_u; //subj random effects


  for(i in 1:N){
    mu[i,1] = beta[1] + FS[1,i] + 1*(FS[2,i]+beta[2]);
    mu[i,2] = beta[1] + FS[1,i] + 2*(FS[2,i]+beta[2]);
    mu[i,3] = beta[1] + FS[1,i] + 4*(FS[2,i]+beta[2]);
    mu[i,4] = beta[1] + FS[1,i] + 6*(FS[2,i]+beta[2]);
  }
}

model{
  L_u ~ lkj_corr_cholesky(2.0);
  to_vector(z_u) ~ normal(0,1);
  
  for(i in 1:N){    
    X[i] ~ normal(mu[i],sigma);
  }
}
"

fa.model=stan(model_code=mod.stan,
                data = dat2,chains=1,
              pars=c("sigma","L_u","sigma_u","beta"))

fa.model
