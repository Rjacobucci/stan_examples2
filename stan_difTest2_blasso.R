
library(mirt)
library(rstan)

#simulate data where group 2 more extreme intercepts
set.seed(whichRep)
I <- 15
lsd_a <- .5
sd_b <- 1
lmu_a <- 0
mu_b <- 0
uformDiff <- 0 ### Uniform Diff Effect
a1 <- a2 <- matrix(rlnorm(I,lmu_a,lsd_a), ncol=1)
b1 <- b2 <- matrix(rnorm(I,mu_b,sd_b),ncol=1)
d1 <- -a1*b1
d2 <- -a2*b2
d1[1:2] <- d1[1:2] - uformDiff
## Leave out the non-uniform diff for now.
## a2[1:2, ] <- a1[1:2, ]/3 

=======
library(mirt);library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#simulate data where group 2 has a smaller slopes and more extreme intercepts
set.seed(12345)
a1 <- a2 <- matrix(abs(rnorm(15,1,.3)), ncol=1)
d1 <- d2 <- matrix(rnorm(15,0,.7),ncol=1)
a2[1:2, ] <- a1[1:2, ]/3
d1[c(1,3), ] <- d2[c(1,3), ]/4
>>>>>>> origin/master
head(data.frame(a.group1 = a1, a.group2 = a2, d.group1 = d1, d.group2 = d2))
itemtype <- rep('2PL', nrow(a1))
N <- 1000

<<<<<<< HEAD
thetaMu <- .1
thetaSd <- sqrt(1.5)
dataset1 <- simdata(a1, d1,N, itemtype = itemtype)
dataset2 <- simdata(a2, d2, N, itemtype, mu = thetaMu, sigma = matrix(thetaSd))
=======
dataset1 <- simdata(a1, d1, N, itemtype)
dataset2 <- simdata(a2, d2, N, itemtype, mu = .1, sigma = matrix(1.5))
>>>>>>> origin/master
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
x <- ifelse(group=="D1",0,1)
ii<-as.numeric(df$Item)
I <- 15 # Items
J <- nrow(dat) #Person
N=I*J
Ik <- ifelse((df$Item == "Item_1" | df$Item == "Item_2" | df$Item == "Item_3"), 1, 0)
irt.data<-list(J=J,I=I,N=N,jj=jj,x=x,ii=ii,y=y,Ik=Ik)




# try stan
# items 1-3 have DIF, c(1,2) for a, c(1,3) for d

#https://github.com/danielcfurr/example-models/blob/master/education/hierarchical_2pl/hierarchical_2pl.stan

mod.stan <-"
data {

int<lower=1> I;               // # questions
int<lower=1> J;               // # persons
int<lower=1> N;               // # observations
int<lower=1, upper=I> ii[N];  // question for n
int<lower=1, upper=J> jj[N];  // person for n
int<lower=0, upper=1> y[N];   // correctness for n
real x[J];                    // covariate for person j
int<lower=0, upper=1> Ik[N];   // Indicator for item k

}
parameters {
real <lower=0>lambda;
real <lower=0>tau;
real <lower=0>ptau;
vector<lower=0>[I] alpha;     // discrimination for item i
vector[I] beta;               // difficulty for item i
real gamma;                        // regression coefficient of x
vector[J] epsilon;                 // error term in the regression model

vector[I] delta;                   // DIF parameter for item k
}
transformed parameters{
real<lower=0> pen;
real<lower=0> ppsi;

ppsi = pow(lambda,-1);
pen = ppsi*ptau;
}
model {
vector[N] eta;   
vector[J] theta;              // ability for person j
alpha ~ lognormal(0.5,1);


lambda ~ gamma(1,.05);
tau ~ gamma(1,.05);
ptau ~ gamma(1,tau/2);
delta ~ normal(0,pen);

beta ~ normal(0,1);
epsilon ~ normal(0,1);
for (j in 1:J)
theta[j] = (gamma * x[j]) + epsilon[j];
for (n in 1:N)
eta[n] = (alpha[ii[n]]) *
(theta[jj[n]] - (beta[ii[n]] + delta[ii[n]] * x[jj[n]]));
y ~ bernoulli_logit(eta);
}
"

sim_fit <- stan(model_code=mod.stan, pars=c("alpha","beta","gamma","delta"),
                data = irt.data, chains = 1, iter = 1000)
sim_fit





# try out shit


library(brms)

dat <- data.frame(response=y,id=jj,item=ii)
get_prior(response ~ item + (item|id),data=dat)

brm.out <- brm(response ~ item + (item|id),data=dat,chains=1,family=bernoulli())
brm.out
=======
vector[I] delta1;                   // DIF parameter for item k
vector<lower=0>[I] lambda;
vector<lower=0>[I] tau;
vector<lower=0>[I] ptau;
}
transformed parameters{
vector<lower=0>[I] pen;
vector<lower=0>[I] ppsi;
for(p in 1:(I)){
ppsi[p] = pow(lambda[p],-1);
pen[p] = ppsi[p]*ptau[p];
}
}
model{
vector[N] eta;   
vector[J] theta;              // ability for person j
alpha ~ lognormal(0.5,1);
lambda ~ gamma(1,.05);
tau ~ gamma(1,.05);
ptau ~ gamma(1,tau/2);
for(p in 1:I){
  delta1[p] ~ normal(0,pen[p]);
}
beta ~ normal(0,10);
epsilon ~ normal(0,1);
for (j in 1:J){
theta[j] = (gamma * x[j]) + epsilon[j];
}

for (n in 1:N){
eta[n] = (alpha[ii[n]]) * (theta[jj[n]] - (beta[ii[n]] + delta1[ii[n]] * x[jj[n]]));
}

y ~ bernoulli_logit(eta);
}
"


sim_fit <- stan(model_code=mod.stan, pars=c("alpha","beta","gamma","delta1"),
                data = irt.data, chains = 2, iter = 500,control=list(adapt_delta=0.95))
sim_fit

saveRDS(sim_fit,"balasso_dif_ex.rds")

>>>>>>> origin/master
