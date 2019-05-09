library(rstan);library(blavaan)

wisc <- read.table("C:/Users/RJacobucci/Documents/GitHub/EDM_Labs/2015/wisc4vpe.dat")
wisc <- read.table("/Users/rjacobuc/Documents/GitHub/EDM_Labs/2015/wisc4vpe.dat")
names(wisc)<- c("V1","V2","V4","V6","P1","P2","P4", "P6", "Moeducat")


bivariate_mod <- "

###
### Verbal Univariate
###

lV1 =~ 1*V1
lV2 =~ 1*V2
lV4 =~ 1*V4
lV6 =~ 1*V6

lV2 ~ 1*lV1; lV4 ~ 1*lV2; lV6 ~ 1*lV4

dV1 =~ 1*lV2; dV2 =~ 1*lV4; dV3 =~ 1*lV6

inv1 =~ 1*lV1;
slope1 =~ 1*dV1 + 2*dV2 + 2*dV3

V1 ~ 0*1; V2 ~0*1; V4 ~ 0*1; V6 ~ 0*1

slope1 ~ 1;
inv1 ~ 1;

slope1 ~~ slope1;
inv1 ~~ inv1;
slope1 ~~ inv1;

lV1 ~ 0*1; lV2 ~0*1; lV4 ~ 0*1; lV6 ~ 0*1
dV1 ~ 0*1; dV2 ~0*1; dV3 ~ 0*1

lV1 ~~ 0*lV1; lV2 ~~ 0*lV2; lV4 ~~ 0*lV4; lV6 ~~ 0*lV6
dV1 ~~ 0*dV1; dV2 ~~ 0*dV2; dV3 ~~ 0*dV3

dV1 ~ beta1*lV1; dV2 ~ beta1*lV2; dV3 ~ beta1*lV4;

V1 ~~ resid1*V1; V2 ~~ resid1*V2; V4 ~~ resid1*V4; V6 ~~ resid1*V6;

###
### Performance Univariate Model
###

lP1 =~ 1*P1
lP2 =~ 1*P2
lP4 =~ 1*P4
lP6 =~ 1*P6

lP2 ~ 1*lP1; lP4 ~ 1*lP2; lP6 ~ 1*lP4

dP1 =~ 1*lP2; dP2 =~ 1*lP4; dP3 =~ 1*lP6

inv2 =~ 1*lP1;
slope2 =~ 1*dP1 + 2*dP2 + 2*dP3

P1 ~ 0*1; P2 ~0*1; P4 ~ 0*1; P6 ~ 0*1

slope2 ~ 1;
inv2 ~ 1;

slope2 ~~ slope2;
inv2 ~~ inv2;
slope2 ~~ inv2;

lP1 ~ 0*1; lP2 ~0*1; lP4 ~ 0*1; lP6 ~ 0*1
dP1 ~ 0*1; dP2 ~0*1; dP3 ~ 0*1

lP1 ~~ 0*lP1; lP2 ~~ 0*lP2; lP4 ~~ 0*lP4; lP6 ~~ 0*lP6
dP1 ~~ 0*dP1; dP2 ~~ 0*dP2; dP3 ~~ 0*dP3

dP1 ~ beta2*lP1; dP2 ~ beta2*lP2; dP3 ~ beta2*lP4;

P1 ~~ resid2*P1; P2 ~~ resid2*P2; P4 ~~ resid2*P4; P6 ~~ resid2*P6;


####
#### Bivariate Model
####

inv1~~slope2
inv1~~inv2
slope1~~slope2
slope1~~inv2

V1 ~~ cov1*P1;V2 ~~ cov1*P2;V4 ~~ cov1*P4;V6 ~~ cov1*P6

dP1 ~ coup1*lV1;dP2 ~ coup1*lV2;dP3 ~ coup1*lV4
dV1 ~ coup2*lP1;dV2 ~ coup2*lP2;dV3 ~ coup2*lP4
"

fit1 <- lavaan(bivariate_mod, data=wisc,control=list(iter.max=100000))
summary(fit1)

