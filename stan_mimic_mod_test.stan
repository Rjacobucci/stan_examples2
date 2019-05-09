data{
int N; // sample size
matrix[N,9] X; // data matrix of order [N,P]
matrix[N,3] cov; // data matrix of order [N,P]
matrix[N,9] X_test; // data matrix of order [N,P]
matrix[N,3] cov_test; // data matrix of order [N,P]
}

parameters{
vector[N] FS; // factor scores, matrix of order [N,D]
vector[N] FS_test; // factor scores, matrix of order [N,D]
vector<lower=0>[9] sigma;
vector[8] lam;
vector[9] alpha;
vector[9] alpha_test;
vector[3] beta;
real<lower=0> psi;
real<lower=0> psi2;
}

transformed parameters{
vector[8] lam_test;
vector[N] mu2_test;
vector<lower=0>[9] sigma_test;
vector[9] mu[N];
vector[9] mu_test[N];
vector[N] mu2;
vector[3] beta_test;

for(j in 1:8){
lam_test[j] = lam[j];
}
for(j in 1:9){
sigma_test[j] = sigma[j];
}
for(j in 1:3){
beta_test[j] = beta[j];
}


for(i in 1:N){
mu[i,1] = alpha[1] + 1*FS[i];
mu[i,2] = alpha[2] + lam[1]*FS[i];
mu[i,3] = alpha[3] + lam[2]*FS[i];
mu[i,4] = alpha[4] + lam[3]*FS[i];
mu[i,5] = alpha[5] + lam[4]*FS[i];
mu[i,6] = alpha[6] + lam[5]*FS[i];
mu[i,7] = alpha[7] + lam[6]*FS[i];
mu[i,8] = alpha[8] + lam[7]*FS[i];
mu[i,9] = alpha[9] + lam[8]*FS[i];

mu2[i] =  beta[1]*cov[i,1] + beta[2]*cov[i,2] + beta[3]*cov[i,3];
mu2_test[i] =  beta_test[1]*cov_test[i,1] + beta_test[2]*cov_test[i,2] + beta_test[3]*cov_test[i,3];
}
for(i in 1:N){
mu_test[i,1] = alpha_test[1] + 1*FS_test[i];
mu_test[i,2] = alpha_test[2] + lam_test[1]*FS_test[i];
mu_test[i,3] = alpha_test[3] + lam_test[2]*FS_test[i];
mu_test[i,4] = alpha_test[4] + lam_test[3]*FS_test[i];
mu_test[i,5] = alpha_test[5] + lam_test[4]*FS_test[i];
mu_test[i,6] = alpha_test[6] + lam_test[5]*FS_test[i];
mu_test[i,7] = alpha_test[7] + lam_test[6]*FS_test[i];
mu_test[i,8] = alpha_test[8] + lam_test[7]*FS_test[i];
mu_test[i,9] = alpha_test[9] + lam_test[8]*FS_test[i];

}
}

model{
sigma ~ gamma(2,2);
alpha ~ normal(3,10);
alpha_test ~ normal(3,10);
psi ~ gamma(2,2);
psi2 ~ gamma(2,2);

for(i in 1:N){  
  for(j in 1:9){
    X[i,j] ~ normal(mu[i,j],pow(sigma[j],0.5));
    X_test[i,j] ~ normal(mu_test[i,j],pow(sigma_test[j],0.5));
  }
FS[i] ~ normal(mu2[i],psi);
FS_test[i] ~ normal(mu2_test[i],psi2);
}
}
generated quantities{
// can use mean() and sd()
real rsq1;
vector[N] z_mu1;
vector[N] z_FS1;
real rsq2;
vector[N] z_mu2;
vector[N] z_FS2;
z_mu1 = (mu2 - mean(mu2))/sd(mu2);
z_FS1 = (FS - mean(FS))/sd(FS);
rsq1 = pow((z_mu1' * z_FS1)/N,2);
z_mu2 = (mu2_test - mean(mu2_test))/sd(mu2_test);
z_FS2 = (FS_test - mean(FS_test))/sd(FS_test);
rsq2 = pow((z_mu2' * z_FS2)/N,2);
}