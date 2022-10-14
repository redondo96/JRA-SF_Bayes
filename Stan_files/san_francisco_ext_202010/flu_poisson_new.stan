functions {
  vector SEIR(real time, vector y, real[] params) {
    vector[8] dydt;
    real gamma_1;
    real DR;
    real gamma_2;
    real RR2;
    real RR12;
    real RR11;
    real ER;
    real HR6;
    real HR5;
    real HR4;
    real HR3;
    real HR2;
    real HR1;
    real N;
    real delta;
    real MR;
    real PR;
    real AR;
    real lambda;
    real IR;
    real BR;
    gamma_1 = params[4]*params[5];
    DR = y[4]*params[5];
    gamma_2 = 1/((1/gamma_1)-(1/params[5]));
    RR2 = y[5]*gamma_2;
    RR12 = y[4]*gamma_1;
    RR11 = y[3]*gamma_1;
    ER = y[2]*0.5263157895;
    HR6 = 4.56621e-05*y[6];
    HR5 = 4.56621e-05*y[5];
    HR4 = 4.56621e-05*y[4];
    HR3 = 4.56621e-05*y[3];
    HR2 = 4.56621e-05*y[2];
    HR1 = 4.56621e-05*y[1];
    N = y[1]+y[2]+y[4]+y[3]+y[5]+y[6];
    delta = (0.067/(1-0.067))*(4.56621e-05+gamma_2);
    MR = y[5]*delta;
    PR = params[3]*ER;
    AR = (1-params[3])*ER;
    lambda = ((y[4]+y[5]+params[2]*y[3])*params[1])/N;
    IR = lambda*y[1];
    BR = 4.56621e-05*N;
    dydt[1] = BR-IR-HR1;
    dydt[2] = IR-ER-HR2;
    dydt[3] = AR-RR11-HR3;
    dydt[4] = PR-DR-RR12-HR4;
    dydt[5] = DR-RR2-MR-HR5;
    dydt[6] = RR11+RR12+RR2-HR6;
    dydt[7] = MR;
    dydt[8] = DR;
    return dydt;
  }
}
data {
  int<lower = 1> n_obs;
  int<lower = 1> n_params;
  int<lower = 1> n_difeq;
  int y[n_obs];
  real t0;
  real ts[n_obs];
}
parameters {
  real<lower = 0>            beta;
  real<lower = 0, upper = 1> q;
  real<lower = 0, upper = 1> rho;
  real<lower = 0, upper = 1> omega;
  real<lower = 0, upper = 1> alpha;
  real<lower = 0, upper = 55000> I0;
  real<lower = 0, upper = 55000> E0;
}
transformed parameters{
  vector[n_difeq] o[n_obs]; // Output from the ODE solver
  real x[n_obs];
  vector[n_difeq] x0;
  real params[n_params];
  x0[1] = 550000 - E0 - I0;
  x0[2] = E0;
  x0[3] = 0;
  x0[4] = I0;
  x0[5] = 0;
  x0[6] = 0;
  x0[7] = 0;
  x0[8] = I0;
  params[1] = beta;
  params[2] = q;
  params[3] = rho;
  params[4] = omega;
  params[5] = alpha;
  o = ode_rk45(SEIR, x0, t0, ts, params);
  x[1] =  o[1, 8]  - x0[8];
  for (i in 1:n_obs-1) {
    x[i + 1] = o[i + 1, 8] - o[i, 8] + 1e-5;
  }
}
model {
  omega   ~ uniform(0.01, 0.99);
  beta    ~ lognormal(0.5, 1);
  q       ~ beta(0.75, 30);
  rho     ~ beta(10, 40);
  alpha   ~ lognormal(-0.7, 0.2);
  I0      ~ lognormal(1, 1);
  E0      ~ lognormal(1, 1);
  y       ~ poisson(x);
}
generated quantities {
  real log_lik;
  log_lik = poisson_lpmf(y | x);
}
