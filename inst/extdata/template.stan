//
data {
  int<lower=0> R;         // Number of areas
  int<lower=0> T;         // Number of time steps

  int N[R, T];            // Number of cases

  real Nsum[R, T];

  int distrib;

  x_distmat_x             // matrix[R, R] MIJ;
  x_popmat_x              // matrix[R, R] popmat;
  x_adjmat_x              // matrix[R, R] adjmat;
  x_con_mat_x
}

parameters {
  real<lower=0> gamma;           //Local transmission
  x_overdispersion_x             //real<lower=0> beta;
  x_spat_int_par_x               //real<lower=0> alpha_spat;          //spatial interaction
  x_adj_par_x                    //real<lower=0> alpha_adj;      //adjacency term

  x_con_par_x
  real<lower=0> epsilon;         //error
  x_power_for_dist_x              //real<lower=0> k;            //distance exponent
}

transformed parameters {

  matrix[R, T] foi;             //force of incection matrix

  x_gravity_law_x


  // calculate force of infection matrix
  foi = (diag_matrix(rep_vector(gamma, R)) x_gravity_interaction_x x_adjacency_interaction_x x_con_mat_interaction_x) *
        to_matrix(Nsum) +
        epsilon;
}




model {
  x_priors_grav_x   //gamma ~ gamma(1, 1); alpha_spat ~ gamma(1, 1); k ~ normal(1, 1);
  x_priors_adj_x    //alpha_adj  ~ gamma(1, 1);
  x_priors_con_x
  //tau1 ~ normal(1, 1);
  //tau2 ~ normal(1, 1);


  x_priors_over_x    //  beta ~ gamma(1, 1);

  for (t in 2:T) {
    for (i in 1:R) {
      if (distrib == 1){
        N[i, t] ~ neg_binomial_2(foi[i, t], 1 / beta);
      }
      if (distrib == 0)
      {
        N[i, t] ~ poisson(foi[i, t]);
      }
    }
  }
}



generated quantities {

  vector[R * (T-1)] log_lik;
  if (distrib == 1){
      for (i in 1:R) {
        for (t in 2:T) {
          log_lik[(t-2) * R + i] = neg_binomial_2_lpmf(N[i, t] | foi[i, t], 1 / beta);
            }
        }
    }

  if (distrib == 0){
    for (i in 1:R) {
      for (t in 2:T) {
        log_lik[(t-2) * R + i] = poisson_lpmf(N[i, t] | foi[i, t]);
          }
      }
    }
}

