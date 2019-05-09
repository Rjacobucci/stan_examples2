library(blavaan)

wisc <- read.table("C:/Users/RJacobucci/Documents/GitHub/EDM_Labs/2015/wisc4vpe.dat")
wisc <- read.table("C:/Users/jacobucc/Documents/GitHub/EDM_Labs/2015/wisc4vpe.dat")
names(wisc)<- c("V1","V2","V4","V6","P1","P2","P4", "P6", "Moeducat")


X <- wisc[,c("V1","V2","V4","V6")]
data = list()
#data$alpha <- c(1,1,2,2)
X = as.matrix(X)
N <- nrow(X)

dat <- list(
  N = N,
  X = X,
  t = 4,
  R = diag(2))








stanmodel <- "
data {
  int N; 
  int t; 
  matrix[2,2] R;
  vector[t] X[N]; 
} 
parameters {
  vector[2] beta;  
  //real pi;
  real<lower=0> sigma;
  cov_matrix[2] phi;
} 
transformed parameters {
  vector[2] s[N];
  vector[t] y[N];
  vector[t-1] d[N];

for (i in 1:N){
  y[i,1] = s[i,1];

  for (tt in 2:t){
    d[i,tt-1] = s[i,2];
    y[i,tt] = d[i,tt-1]+y[i,tt-1];
  }
}
  
} 

model {

for (i in 1:N){
  s[i,1:2] ~ multi_normal(beta, phi);
  X[i,1] ~ normal(y[i,1],sigma);
  

  for (tt in 2:t){
    X[i,tt] ~ normal(y[i,tt], sigma);
  }
}

sigma ~ gamma(2,2);
//pi~normal(0,2);
beta[1]~normal(0,100);
beta[2]~normal(0,100);
phi~ inv_wishart(2,R);

}
"

init.list <- list()
beta <- c(20,4)
sigma <- 13
phi <- matrix(c(20,3,3,1),2,2)
init.list[[1]]  <- list(beta=beta,sigma=sigma,phi=phi)


library(rstan)

fit <- stan(model_code = stanmodel, model_name = "LCS", chains=1,
          data = dat,init=init.list,
        pars=c("beta","sigma"))

print(fit)








###### -----------------------------------------

###### -----------------------------------------

# latent growth curve

dat2 <- list(
  N = N,
  X = X,
  t = 4,
  R = diag(2))

stanmodel2 <- "
data {
int N; 
int t; 
matrix[2,2] R;
vector[t] X[N]; 
} 
parameters {
vector[2] beta;  
//real pi;
real<lower=0> sigma;
cov_matrix[2] phi;
} 
transformed parameters {
vector[2] s[N];
vector[t] mu[N];

for (i in 1:N){

for (tt in 1:t){
mu[i,tt] = s[i,1] + tt*s[i,2];
}
}

} 

model {

for (i in 1:N){
s[i,1:2] ~ multi_normal(beta, phi);

for (tt in 1:t){
X[i,tt] ~ normal(mu[i,tt], sigma);
}
}

sigma ~ gamma(3,3);
//pi~normal(0,2);
beta[1]~normal(0,100);
beta[2]~normal(0,100);
phi~ inv_wishart(2,R);

}
"



library(rstan)

fit2 <- stan(model_code = stanmodel2, model_name = "LGCM", chains=1,
            data = dat2,
            pars=c("beta","sigma"))

print(fit2)
