

construct_stan_model <- function(model_path, interactions=c(1,4))

{ 
  
  
  includeparam = list(
  x_distmat_x = c("x_distmat_x", "matrix[R, R] MIJ; "),
  x_popmat_x = c("x_popmat_x", "matrix[R, R] popmat;"),
  x_adjmat_x = c("x_adjmat_x","matrix[R, R] adjmat;") ,
  x_spat_int_par_x = c("x_spat_int_par_x", "real<lower=0> alpha_spat; "),
  x_adj_par_x = c("x_adj_par_x", "real<lower=0> alpha_adj; "),
  x_power_for_dist_x = c("x_power_for_dist_x", "real<lower=0> k; "), 
  
  
  
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
  
  x_overdispersion_x = c("x_overdispersion_x", "real<lower=0> beta;"), 
  
  x_priors_grav_x = c("x_priors_grav_x",
  "gamma ~ gamma(1, 1);
   alpha_spat ~ gamma(1, 1);
   k ~ normal(1, 1);" ), 
  
  
  
  x_priors_adj_x = c("x_priors_adj_x", "alpha_adj  ~ gamma(1, 1);"), 
  
  x_priors_over_x = c("x_priors_over_x", "beta ~ gamma(1, 1);")
  )
  
  
  all_options = c("x_distmat_x", "x_popmat_x", "x_spat_int_par_x", "x_power_for_dist_x",  "x_gravity_law_x", 
                  "x_gravity_interaction_x", "x_adjmat_x", "x_adj_par_x", "x_adjacency_interaction_x", 
                  "x_overdispersion_x", "x_priors_grav_x", "x_priors_adj_x", "x_priors_over_x")
  
  to_include = list( 
    c("x_distmat_x", "x_popmat_x", "x_spat_int_par_x", "x_power_for_dist_x",  "x_gravity_law_x", "x_gravity_interaction_x", "x_priors_grav_x"), 
    c("x_distmat_x", "x_popmat_x", "x_spat_int_par_x", "x_power_for_dist_x",  "x_gravity_law_x", "x_gravity_interaction_x", "x_priors_grav_x"),
    c("x_distmat_x", "x_popmat_x", "x_spat_int_par_x", "x_power_for_dist_x",  "x_gravity_law_x", "x_gravity_interaction_x", "x_priors_grav_x"), 
    c("x_adjmat_x", "x_adj_par_x", "x_adjacency_interaction_x", "x_priors_adj_x"),
    c("x_overdispersion_x", "x_priors_over_x")
    ) 

   
   
   please_include = c()
   
   for (item in interactions){ 
     
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