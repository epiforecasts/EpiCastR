#' Construct stan model
#' @param model_path path to stan template to build model
#' @param interactions vector with numbers to indicate the types of interaction to include in the model: 1. Gravity model. 2. Gravity model with population density. 3. Power law (no population info). 4. Adjacency model.
#' @export
construct_stan_model <- function(model_path, interactions=c(1,4))

{


  includeparam = list(
  x_distmat_x = c("x_distmat_x", "matrix[R, R] MIJ; "),
  x_popmat_x = c("x_popmat_x", "matrix[R, R] popmat;"),
  x_adjmat_x = c("x_adjmat_x","matrix[R, R] adjmat;") ,
  x_spat_int_par_x = c("x_spat_int_par_x", "real<lower=0> alpha_spat; "),
  x_adj_par_x = c("x_adj_par_x", "real<lower=0> alpha_adj; "),
  x_con_par_x = c("x_con_par_x", "real<lower=0> alpha_con; "),
  x_power_for_dist_x = c("x_power_for_dist_x", "real<lower=0> k; "),
  x_con_mat_x = c("x_con_mat_x", "matrix[R, R] con_mat;"),



  x_gravity_law_x = c("x_gravity_law_x",

    "real s[2];                    //to offer 1 or grav
     matrix[R, R] space;           //spatial interaction matrix

     s[1] = 1.;

      for (i in 1:R) {
        space[i, i] = 0;
        for (j in 1:(i-1)) {
          s[2] = popmat[i,j] /pow(MIJ[i, j], k);
          space[i, j] = max(s);
          space[j, i] = space[i, j];
        }

      }
        space ./= rep_matrix(rep_row_vector(1, R) * space, R);

      "),
  x_gravity_interaction_x =  c("x_gravity_interaction_x" , "+  alpha_spat * space"),
  x_adjacency_interaction_x = c("x_adjacency_interaction_x", "+  alpha_adj * adjmat" ),
  x_con_mat_interaction_x = c("x_con_mat_interaction_x", "+  alpha_con * con_mat" ),

  x_overdispersion_x = c("x_overdispersion_x", "real<lower=0> beta;"),

  x_priors_grav_x = c("x_priors_grav_x",
  "gamma ~ gamma(1, 1);
   alpha_spat ~ gamma(1, 1);
   k ~ normal(1, 1);" ),



  x_priors_adj_x = c("x_priors_adj_x", "alpha_adj  ~ gamma(1, 1);"),
  x_priors_con_x = c("x_priors_con_x", "alpha_con  ~ gamma(1, 1);"),

  x_priors_over_x = c("x_priors_over_x", "beta ~ gamma(1, 1);")
  )


  all_options = c("x_distmat_x", "x_popmat_x", "x_spat_int_par_x", "x_power_for_dist_x",  "x_gravity_law_x",
                  "x_gravity_interaction_x", "x_adjmat_x", "x_adj_par_x", "x_adjacency_interaction_x",
                  "x_overdispersion_x", "x_priors_grav_x", "x_priors_adj_x", "x_priors_over_x", "x_con_mat_x", "x_con_par_x", "x_con_mat_interaction_x", "x_priors_con_x")

  to_include = list(
    c("x_distmat_x", "x_popmat_x", "x_spat_int_par_x", "x_power_for_dist_x",  "x_gravity_law_x", "x_gravity_interaction_x", "x_priors_grav_x"),
    c("x_distmat_x", "x_popmat_x", "x_spat_int_par_x", "x_power_for_dist_x",  "x_gravity_law_x", "x_gravity_interaction_x", "x_priors_grav_x"),
    c("x_distmat_x", "x_popmat_x", "x_spat_int_par_x", "x_power_for_dist_x",  "x_gravity_law_x", "x_gravity_interaction_x", "x_priors_grav_x"),
    c("x_adjmat_x", "x_adj_par_x", "x_adjacency_interaction_x", "x_priors_adj_x"),
    c("x_overdispersion_x", "x_priors_over_x"),
    c("x_con_mat_x", "x_con_par_x", "x_con_mat_interaction_x", "x_priors_con_x")
    )



   please_include = c()

   for (item in interactions){
     print(item)

     please_include = append(please_include, to_include[[item]])

     }

   please_include = unique(please_include)

   dont_include = all_options[!(all_options %in% please_include)]



   model = readChar(model_path, file.info(model_path)$size)

   for (item in please_include){

     model = gsub(includeparam[[item]][1], includeparam[[item]][2], model)

   }





   for (item in dont_include){
     model = gsub(item, ' ', model)

   }

   model
#
  }
