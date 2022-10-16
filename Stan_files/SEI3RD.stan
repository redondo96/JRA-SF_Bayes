functions {
  vector X_model(real time, vector y, array[] real params) {
    vector[8] dydt;
    real E_to_I;
    real I_to_J;
    real E_to_A;
    real N;
    real births;
    real S_deaths;
    real E_deaths;
    real I_deaths;
    real J_deaths;
    real A_deaths;
    real R_deaths;
    real gamma1;
    real gamma2;
    real C_in;
    real S_to_E;
    real J_to_R;
    real A_to_R;
    real I_to_R;
    real var_delta;
    real J_to_D;
    E_to_I = params[4]*0.5263157895*y[2];
    I_to_J = params[5]*y[4];
    E_to_A = (1-params[4])*0.5263157895*y[2];
    N = y[3]+y[2]+y[4]+y[5]+y[6]+y[1];
    births = 4.56621e-05*N;
    S_deaths = 4.56621e-05*y[1];
    E_deaths = 4.56621e-05*y[2];
    I_deaths = 4.56621e-05*y[4];
    J_deaths = 4.56621e-05*y[5];
    A_deaths = 4.56621e-05*y[3];
    R_deaths = 4.56621e-05*y[6];
    gamma1 = params[1]*params[5];
    gamma2 = 1/((1/gamma1)-(1/params[5]));
    C_in = I_to_J;
    S_to_E = params[2]*y[1]*(y[4]+y[5]+params[3]*y[3])/N;
    J_to_R = gamma2*y[5];
    A_to_R = gamma1*y[3];
    I_to_R = gamma1*y[4];
    var_delta = (0.067/(1-0.067))*(4.56621e-05+gamma2);
    J_to_D = var_delta*y[5];
    dydt[1] = births-S_to_E-S_deaths;
    dydt[2] = S_to_E-E_to_I-E_to_A-E_deaths;
    dydt[3] = E_to_A-A_to_R-A_deaths;
    dydt[4] = E_to_I-I_to_J-I_deaths-I_to_R;
    dydt[5] = I_to_J-J_to_D-J_to_R-J_deaths;
    dydt[6] = J_to_R+A_to_R+I_to_R-R_deaths;
    dydt[7] = J_to_D;
    dydt[8] = C_in;
    return dydt;
  }
}
data {
  int<lower = 1> n_obs;
  int<lower = 1> n_params;
  int<lower = 1> n_difeq;
  array[n_obs] int y;
  real t0;
  array[n_obs] real ts;
}
parameters {
  real<lower = 0, upper = 1> par_omega;
  real<lower = 0> par_beta;
  real<lower = 0, upper = 1> par_q;
  real<lower = 0, upper = 1> par_rho;
  real<lower = 0, upper = 1> par_alpha;
  real<lower = 0> E0;
  real<lower = 0> I0;
  real<lower = 0> inv_phi;
}
transformed parameters{
  array[n_obs] vector[n_difeq] x; // Output from the ODE solver
  array[n_params] real params;
  vector[n_difeq] x0; // init values
  array[n_obs] real delta_x_1;
  real phi;
  phi = 1 / inv_phi;
  x0[1] = 550000 - E0 - I0; // S
  x0[2] = E0; // E
  x0[3] = 0; // A
  x0[4] = I0; // I
  x0[5] = 0; // J
  x0[6] = 0; // R
  x0[7] = 0; // D
  x0[8] = 0; // C
  params[1] = par_omega;
  params[2] = par_beta;
  params[3] = par_q;
  params[4] = par_rho;
  params[5] = par_alpha;
  x = ode_rk45(X_model, x0, t0, ts, params);
  delta_x_1[1] =  x[1, 8] - x0[8] + 1e-5;
  for (i in 1:n_obs-1) {
    delta_x_1[i + 1] = x[i + 1, 8] - x[i, 8] + 1e-5;
  }
}
model {
  par_omega ~ beta(2, 2);
  par_beta ~ lognormal(0.5, 1);
  par_q ~ beta(2, 2);
  par_rho ~ beta(2, 2);
  par_alpha ~ beta(2, 2);
  E0 ~ lognormal(1, 1);
  I0 ~ lognormal(1, 1);
  inv_phi ~ exponential(5);
  y ~ neg_binomial_2(delta_x_1, phi);
}
generated quantities {
  real log_lik;
  log_lik = neg_binomial_2_lpmf(y | delta_x_1, phi);
}
