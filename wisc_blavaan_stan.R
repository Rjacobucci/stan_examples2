library(blavaan);library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

wisc <- read.table("C:/Users/RJacobucci/Documents/GitHub/EDM_Labs/2015/wisc4vpe.dat")
wisc <- read.table("C:/Users/jacobucc/Documents/GitHub/EDM_Labs/2015/wisc4vpe.dat")
wisc <- read.table("/Users/rjacobuc/Documents/GitHub/EDM_Labs/2015/wisc4vpe.dat")
names(wisc)<- c("V1","V2","V4","V6","P1","P2","P4", "P6", "Moeducat")



library(blavaan)



lds_1 <- "

#latent variables
lV1 =~ 1*V1
lV2 =~ 1*V2
lV4 =~ 1*V4
lV6 =~ 1*V6



#autoregressions
lV2 ~ 1*lV1; lV4 ~ 1*lV2; lV6 ~ 1*lV4

#change - delta; d
dV1 =~ 1*lV2; dV2 =~ 1*lV4; dV3 =~ 1*lV6

#intercept and slope
inV =~ 1*lV1;

# match lgm
slope =~ 1*dV1 + 2*dV2 + 2*dV3


#manifest means @0
V1 ~ 0*1; V2 ~0*1; V4 ~ 0*1; V6 ~ 0*1

#slope and intercept means
slope ~ 1;
inV ~ 1;

#Latent variances and covariance
slope ~~ slope;
inV ~~ inV;
slope ~~ inV;

#means and vars @0
lV1 ~ 0*1; lV2 ~0*1; lV4 ~ 0*1; lV6 ~ 0*1
dV1 ~ 0*1; dV2 ~0*1; dV3 ~ 0*1

lV1 ~~ 0*lV1; lV2 ~~ 0*lV2; lV4 ~~ 0*lV4; lV6 ~~ 0*lV6
dV1 ~~ 0*dV1; dV2 ~~ 0*dV2; dV3 ~~ 0*dV3

#auVo-proportions
dV1 ~ beta*lV1
dV2 ~ beta*lV2
dV3 ~ beta*lV4

#residuals equal
V1 ~~ resid*V1; V2 ~~ resid*V2; V4 ~~ resid*V4; V6 ~~ resid*V6;
"

bl.out = blavaan(lds_1,data=wisc,burnin=50,adapt=25,sample=200,target="stan")
summary(bl.out)
