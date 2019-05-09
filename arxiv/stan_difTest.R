library(mirt);library(rstan)

#simulate data where group 2 has a smaller slopes and more extreme intercepts
set.seed(12345)
a1 <- a2 <- matrix(abs(rnorm(15,1,.3)), ncol=1)
d1 <- d2 <- matrix(rnorm(15,0,.7),ncol=1)
a2[1:2, ] <- a1[1:2, ]/3
d1[c(1,3), ] <- d2[c(1,3), ]/4
head(data.frame(a.group1 = a1, a.group2 = a2, d.group1 = d1, d.group2 = d2))
itemtype <- rep('2PL', nrow(a1))
N <- 1000

dataset1 <- simdata(a1, d1, N, itemtype)
dataset2 <- simdata(a2, d2, N, itemtype, mu = .1, sigma = matrix(1.5))
dat <- rbind(dataset1, dataset2)
group <- c(rep('D1', N), rep('D2', N))


itemnames <- colnames(dat)
model_anchor <- multipleGroup(dat, model = 1, group = group,
                              invariance = c(itemnames[11:15], 'free_means', 'free_var'))
anchor <- DIF(model_anchor, c('a1', 'd'), items2test = 1:10)
anchor

dat <- as.data.frame(dat)
dat$id <- 1:nrow(dat)
dat$group <- ifelse(group=="D1",0,1)
df<-reshape2::melt(dat,id=c("id","group"))
colnames(df)<-c("Person","Group","Item","Response")

df<-plyr::arrange(df, Person) # sort by person
head(df)


y <- df$Response
jj<-df$Person
x <- df$Group
ii<-as.numeric(df$Item)
I <- 15 # Items
J <- nrow(dat) #Person
N=I*J
irt.data<-list(J=J,I=I,N=N,jj=jj,ii=ii,y=y,x=x)




# try stan
# items 1-3 have DIF, c(1,2) for a, c(1,3) for d

#https://github.com/danielcfurr/example-models/blob/master/education/hierarchical_2pl/hierarchical_2pl.stan

mod.stan <-"
data {
  int<lower=1> I;               // # items
int<lower=1> J;               // # persons
int<lower=1> N;               // # observations
int<lower=1, upper=I> ii[N];  // item for n
int x[N];                     // covariate
int<lower=1, upper=J> jj[N];  // person for n
int<lower=0, upper=1> y[N];   // correctness for n
}
parameters {
vector[J] theta;              // abilities
vector[2] xi[I];              // alpha/beta pair vectors
vector[2] mu;                 // vector for alpha/beta means
vector<lower=0>[2] tau;       // vector for alpha/beta residual sds
cholesky_factor_corr[2] L_Omega;
}
transformed parameters {
vector[I] alpha;
vector[I] beta;
for (i in 1:I) {
alpha[i] = exp(xi[i,1]);
beta[i] = xi[i,2];
}
}
model {
matrix[2,2] L_Sigma;
L_Sigma = diag_pre_multiply(tau, L_Omega);
L_Omega ~ lkj_corr_cholesky(4);
mu[1] ~ normal(0,1);
tau[1] ~ exponential(.1);
mu[2] ~ normal(0,5);
tau[2] ~ exponential(.1);
//epsilon ~ normal(0,1);

for (j in 1:J){
  theta[j] ~ normal(gamma * x[j],1);
}


for (i in 1:I){
xi[i] ~ multi_normal_cholesky(mu, L_Sigma);
y ~ bernoulli_logit(alpha[i] .* (theta[jj] - beta[i]));
}


generated quantities {
corr_matrix[2] Omega;
Omega <- multiply_lower_tri_self_transpose(L_Omega);
}
"


sim_fit <- stan(model_code=mod.stan, pars=c("Omega","mu","alpha","beta","gamma"),
                data = irt.data, chains = 1, iter = 500)
sim_fit

