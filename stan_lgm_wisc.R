library(rstan)



wisc <- read.table("C:/Users/RJacobucci/Documents/GitHub/EDM_Labs/2015/wisc4vpe.dat")
wisc <- read.table("C:/Users/jacobucc/Documents/GitHub/EDM_Labs/2015/wisc4vpe.dat")
wisc <- read.table("C:/Users/Ross/Documents/GitHub/EDM_Labs/2015/wisc4vpe.dat")
wisc <- read.table("/Users/rjacobuc/Documents/GitHub/EDM_Labs/2015/wisc4vpe.dat")
names(wisc)<- c("V1","V2","V4","V6","P1","P2","P4", "P6", "Moeducat")



change1 = wisc[,2] - wisc[,1]
ind1 = change1 > 8.67

change2 = wisc[,3] - wisc[,2]

ind2 = change2>10.2

table(ind1,ind2)


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
vector[2] FS[N]; // factor scores, matrix of order [N,D]
cholesky_factor_cov[2] Rho; // correlation matrix between factors
real sigma; // variance for each variable
vector[2] M;
}

transformed parameters{

vector[t] mu[N];
matrix[2,2] Sd_d;
Sd_d = multiply_lower_tri_self_transpose(Rho);

  for(i in 1:N){
    mu[i,1] = FS[i,1] + 1*FS[i,2];
    mu[i,2] = FS[i,1] + 2*FS[i,2];
    mu[i,3] = FS[i,1] + 4*FS[i,2];
    mu[i,4] = FS[i,1] + 6*FS[i,2];
  }
}

model{
  M ~ normal(0,2);
  Rho ~ lkj_corr_cholesky(2.0);
  
  for(i in 1:N){    
    X[i] ~ normal(mu[i],pow(sigma,0.5));
    FS[i] ~ multi_normal(M, Sd_d);
  }
}
"

fa.model=stan(model_code=mod.stan,
                data = dat2,chains=1,
                pars=c("sigma","M","Sd_d"))

fa.model
