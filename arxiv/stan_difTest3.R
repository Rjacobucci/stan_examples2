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
vector<lower=0>[I] alpha;     // discrimination for item i
vector[I] beta;               // difficulty for item i
real gamma;                        // regression coefficient of x
vector[J] epsilon;                 // error term in the regression model
real delta;                   // DIF parameter for item k
}
model {
vector[N] eta;   
vector[J] theta;              // ability for person j
alpha ~ lognormal(0.5,1);
beta ~ normal(0,10);
epsilon ~ normal(0,1);
for (j in 1:J)
theta[j] = (gamma * x[j]) + epsilon[j];
for (n in 1:N)
eta[n] = alpha[ii[n]] * (theta[jj[n]] - (beta[ii[n]] + delta * Ik[n] * x[jj[n]]));
y ~ bernoulli_logit(eta);
}
"


sim_fit <- stan(model_code=mod.stan, pars=c("alpha","beta","gamma","delta"),
                data = irt.data, chains = 1, iter = 500)
sim_fit

