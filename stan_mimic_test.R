library(rstan); library(lavaan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
HS <- HolzingerSwineford1939[complete.cases(HolzingerSwineford1939),]

id <- sample(1:nrow(HS),150)
HS.train <- HS[id,]
HS.test <- HS[-id,]



mod <- "
f1 =~ NA*x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9
f1 ~ sex + grade + ageyr
f1~~1*f1
"
out <- sem(mod, HS.train,meanstructure=T)
summary(out,rsq=T)

X <- HS.train[,7:15]
cov <- HS.train[,c(2,3,6)]
X.test <- HS.test[,7:15]
cov.test <- HS.test[,c(2,3,6)]


X = as.matrix(X)
cov= as.matrix(cov)
X.test = as.matrix(X.test)
cov.test= as.matrix(cov.test)
N <- nrow(X)

dat <- list(
  N = N,
  X = X,
  cov = cov,
  X_test = X.test,
  cov_test = cov.test)


setwd("C:/Users/rjacobuc/Documents/Github/stan_examples")


# Diagonal residual covariance matrix
mod.stan <- stan_model(stanc_ret = stanc(file = "stan_mimic_mod_test.stan"))


fa.model=sampling(mod.stan,iter=100,
              data = dat,chains=1,
              pars=c("rsq1","rsq2", "sigma","lam","beta","psi","psi2","alpha","alpha_test"))

print(fa.model)

plot(fa.model,plotfun="trace",pars="rsq2")

tidy_stan(fa.model)
