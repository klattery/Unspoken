// Kevin Lattery Nov 2021
// Conjoint Model in Stan
// con_sign has hard constraints for sign of coded parameter (pos/neg)
// paircon_matrix has pairwise constraints beyond con_sign

functions{
  matrix logistic_hinge(matrix x, matrix delta) {
    return delta .* log1p_exp(x ./ delta);
    // https://statmodeling.stat.columbia.edu/2017/05/19/continuous-hinge-function-bayesian-modeling/
  }
  
  real MNL_LL_par(array[] int array_slice,
                  int a_beg, int a_end, // Stan determines values of these for specific core
                  matrix beta_ind,  matrix X,  vector dep, array[] real wts,
                  array[] int start, array[] int end,
                  array[] int task_individual,
                  array[] int task_type,
                  vector timez_flip, real beta_scale, 
                  vector timescale_con, vector timescale_att,
                  array[] real wts_bound_att, array[] real wts_bound_conv,
                  vector thresh_base,
                  int love_exist, int hate_exist,
                  vector love, vector hate
  ) {
    real ll = 0;
    for (t in a_beg:a_end){
      if (task_type[t] == 1){ // Conversion 
        real wt_task = (wts_bound_conv[1] + (wts_bound_conv[2] - wts_bound_conv[1]) *
                          inv_logit(timez_flip[start[t]] * timescale_con[task_individual[t]])) *
                          wts[t];   
        ll+= wt_task * 
          dot_product(dep[start[t]:end[t]],
                      log_softmax(X[start[t]:end[t]] * col(beta_ind,task_individual[t]) * beta_scale)
          );       
      }
      if (task_type[t] == 2){ // Attraction task has multiple choices
        int num_items =  end[t] - start[t] + 1;
        vector[num_items] util =  X[start[t]:end[t]] * col(beta_ind,task_individual[t]);
        vector[num_items] wt_tasks = (wts_bound_att[1] + (wts_bound_att[2] - wts_bound_att[1]) *
                                        inv_logit(timez_flip[start[t]:end[t]] * timescale_att[task_individual[t]])) *
                                     wts[t]; // weight for each task
        
        vector[num_items] prob = inv_logit(util - thresh_base[task_individual[t]]); // applies scale factor to utilities
        ll+= dot_product(wt_tasks,
                         lmultiply(dep[start[t]:end[t]],     prob) +
                         lmultiply(1-dep[start[t]:end[t]], 1-prob));
        
        // Now Love > Not_Love        
        if (love_exist == 1){
          vector[num_items] util_exp = exp(util); # Used for love/hate  
          vector[num_items] love_d = love[start_t:end_t]; # love data for this task
          ll+= dot_product(
            love_d,
            util - log(util_exp + dot_product(util_exp, 1 - love_d))
          );            
        }
        // Now Hate < Others
        if (hate_exist == 1){
          vector[num_items] util_neg_exp = 1.0/util_exp; # Flip Util  
          vector[num_items] hate_d = hate[start_t:end_t]; # hate data for this task
          ll+= dot_product(
            hate_d,
            (-util) - log(util_neg_exp + dot_product(util_neg_exp, 1 - hate_d))
          );            
        }
    } // end for
    return ll;
  }
  
  int vec_not0(vector vec){
    int count = 0;
    for (i in 1:num_elements(vec)){
      if (abs(vec[i]) > .00000001) count += 1;
    }
    return count;
  }
} // End Functions

data {
  // Sizing Constants, Number of:
  int N; // Rows in data file
  int P; // Parameters = ncol(coded independent data)
  int I; // Individuals
  int T; // Total unique tasks (individuals x mean tasks per individual)
  
  // Main data
  vector<lower = 0, upper = 1>[N] dep; // Dep variable
  array[3] int<lower = 0> sizes; // 1=# columns coded, 2 = # cat vars w/levels, 3=# rows in code_master 
  matrix[N, sizes[1]] ind_coded; // Coded data mapped to ind 
  array[N, sizes[2]]int<lower = 0> ind_levels;
  
  // coding for each attribute combined in code_master, and code_blocks that define attributes in that
  matrix[sizes[3], P] code_master;
  int n_atts; // rows in code blocks = # atts
  array[n_atts, 5] int<lower = 0> code_blocks;
  
  // Upper level priors
  cov_matrix[P] prior_cov;  // Prior covariance matrix, recommend Lenk adjustment
  int<lower = 0> df; // df over P, recommend 5
  vector[P] prior_alpha;
  real<lower = 0> a_sig;  // alpha ~ N(prior_alpha, a_sig)
  
  // Constraints
  int<lower = 0, upper =1> con_use; // 0 = ignore constraints, 1 = use_constraints
  vector[P] con_sign; // Sign of constraints -1, +1 or 0
  int paircon_rows; // # rows in pairs constraint matrix, 0 for none
  matrix[paircon_rows, P] paircon_matrix; // Pair constraints: beta * paircon_matrix > 0
  array[2] real<lower = .01> con_factors; // multiplicative, bound
  
  // Covariates
  int P_cov; // Number of Covariate parameters (coded)
  matrix[I,P_cov] i_cov; // Resp covariates for each individual I
  
  // Misc: Weights, blocking, scale mult, multi-threading split
  array[T] real wts; // array of weights for each respondent
  matrix[P,P] cov_block; // Specifies blocks of covariance items
  real<lower = 0> prior_cov_scale;  // Multiply prior_cov by this:  Typically 1
  int<lower = 0> splitsize; // grainsize for parallel processing 
  
  // Ragged array matching, For each task in 1:T, Specify:
  array[T] int<lower = 1, upper = I> task_individual;
  array[T] int<lower = 1, upper = N> start;
  array[T] int<lower = 1, upper = N> end;
  
  // New Data for Unspoken not in Standard HB 
  array[T] int task_type; // 1 = Attraction, 2 = Conversion
  vector[N] timez_flip; // faster times will be bigger numbers
  
  // Bounds that we set
  array[2] real scale_bounds; // lb, ub for scale converting utilities: conversion to attraction
  array[2] real wts_bound_att;  // lb, ub for weights based on time (attraction)
  array[2] real wts_bound_conv; // lb, ub for weights based on time (conversion)
  
  // Attraction: Distingush Like vs Love, Hate vs Dislike
  vector[N] love; //  1 = Love, 0 = Everything else
  vector[N] hate;  // 1 = Hate, 0 = Everything else
  
  // Hyper parameters for the conversion from time to weights
 // real<lower = 0> timescale_con_hyper; // Default 2: resp scale/2 ~ beta(1,hyper);
//  real<lower = 0> timescale_att_hyper; // Default 2: resp scale/2 ~ beta(1,hyper);
}

transformed data{
  array[I] real beta_scale_1 = rep_array(1, I); // beta scale all 1
  int<lower = 0, upper = 1> love_exist = ((max(love) > 0) ? 1 : 0); // 0 if no love
  int<lower = 0, upper = 1> hate_exist = ((max(hate) > 0) ? 1 : 0); // 0 if hate

  int est_att = 0; // initial, revised later
  int est_conv = 0; // iniial, revised later
  int est_scale_util; 
  
  matrix[N, P] ind; // ind_coded and ind_levels will expand to ind (below)
  matrix[P_cov, I] i_covt = i_cov'; // transpose of i_cov is what we use
  int df_tot = P + df;
  matrix[P,P] L_IW = cholesky_decompose(prior_cov * (df -1)); // For Inverse Wishart

  int con_n = vec_not0(con_sign);
  array[con_n] int con_p;               // Array of which parameters are sign constrained
  matrix[con_n, I] con_delta;     // Con function scale for parameter and respondent
  int paircon_use = 0;
  vector[P] a_LB = rep_vector(-20, P); // LB for alpha (initial)
  vector[P] a_UB = rep_vector(20, P); // UB for alpha (initial)
  vector[P] prior_alpha_mod = prior_alpha; // Will modify alpha for constraints
  
  array[T] int array_slice;  // Parallel threads slice across tasks 1:T 
  int count = 1;
  int lev_col = 1;
  
  for (i in 1:T){
    array_slice[i] = i;
    if (task_type[i] == 1) est_conv = 1; // Whether we have attraction data
    if (task_type[i] == 2) est_att = 1; // Whether we have conversion data
  }
  est_scale_util = est_att * est_conv; // Binary indicator for whether we need utility scale conversion 
  
  for (i in 1:P){
    if (abs(con_sign[i]) > .00000001){
      con_p[count] = i;
      con_delta[count,1:I] = rep_row_vector(con_sign[i] * con_factors[1], I); // con_delta controls hinge function
      count += 1;
    }
  }
  if (con_use == 0) con_n = 0; // set number of ordinal constraints to 0 if no constraints
  if (con_n > 0){
    for (i in 1:P){
      if (con_sign[i] > 0){
        a_LB[i] = -con_factors[2]; // Pos Con LB alpha before hinge
        prior_alpha_mod[i] = -con_factors[2]; // lower bound will be mode
      } 
      if (con_sign[i] < 0){
        a_UB[i] = con_factors[2]; // Neg Con UB alpha before hinge
        prior_alpha_mod[i] = con_factors[2]; // upper bound will be mode
      } 
    }   
  }
  if ((paircon_rows > 0) && (con_use == 1)) paircon_use = 1; // use paircon with rows and flag set

  
  // Coding of attributes with levels
  for (i in 1:n_atts){
    if (code_blocks[i,5] > 0){ //already coded, just copy
      int col_beg = code_blocks[i, 1];
      int col_end = code_blocks[i, 2];
      int col_coded = code_blocks[i, 5]; // beginning of coded column
      ind[:,col_beg:col_end] = ind_coded[:,col_coded:(col_coded + col_end - col_beg)];
    } else {
      int col_beg = code_blocks[i,1];
      int col_end = code_blocks[i,2];
      int row_beg = code_blocks[i,3];
      int row_end = code_blocks[i,4];
      int nrows = row_end - row_beg +1;
      int ncols = col_end - col_beg + 1;
      matrix[nrows, ncols] code_att =  code_master[row_beg:row_end, col_beg:col_end]; // code for this attribute
      row_vector[ncols] zeroes = rep_row_vector(0, ncols);
      for (n in 1:N){
        int lev = ind_levels[n,lev_col];
        if (lev >= 1) ind[n,col_beg:col_end] = code_att[lev,];
        if (lev == 0) ind[n,col_beg:col_end] = zeroes;
      } // end 1:N loop coding
      lev_col = lev_col + 1;
    } // end if  
  } // end for (att coding)  

  print("Estimate Conversion (0/1) = ", est_conv);
  print("Estimate Attraction (0/1) = ", est_att);
  if (est_att == 1) print("Estimate Love = ", love_exist);
} // end transformed data

parameters {
  vector[P] alpha; // upper mean
  cholesky_factor_cov[P] cov_chol;
  matrix[P, I] z; // individual z scores N(0,1)
  matrix[P, P_cov] i_cov_load; // loadings for i_cov (resp covariates)
  
  // resp scale of utilities: conversion -> attraction
  // Use hyper parameters resp scale ~ N(mu, sigma)
  array[est_scale_util] real<lower = scale_bounds[1], upper = scale_bounds[2]> beta_scale;  // attraction to conversion
  vector[I * est_att] thresh_base;  // Attraction: threshold for swipe good or bad each respondent

  vector<lower = 0, upper = 2>[I] timescale_con; // conversion time to weight
  vector<lower = 0, upper = 2>[I] timescale_att; // attraction time to weight
}

transformed parameters {
  matrix[P, I] beta_ind = rep_matrix(alpha, I) + (i_cov_load * i_covt) + cov_chol * z; // Unconstrained;
  if (con_n > 0){
      beta_ind[con_p,1:I] = con_delta .* log1p_exp(beta_ind[con_p,1:I] ./ con_delta);
  } 
}

model {
  // priors on the parameters
  target += normal_lupdf(alpha|prior_alpha_mod, a_sig);
  target += -sum(log1p_exp(-100 * (a_UB - alpha))); // Penalty if alpha > UB
  target += -sum(log1p_exp(-100 * (alpha - a_LB))); // Penalty if alpha < LB
  
  target += inv_wishart_cholesky_lupdf(cov_chol|df_tot, L_IW);
  to_vector(z) ~ std_normal(); // log probabilities of each choice in the dataset
  if (P_cov > 0) to_vector(i_cov_load) ~ std_normal();
  if (paircon_use == 1) target += -sum(log1p_exp(-100 * (paircon_matrix * beta_ind))); // penalty for soft constraints

  if (est_scale_util == 1){ // both types exist --> scale factor
   beta_scale ~ normal(1,2);
  } 

  if (est_conv == 1) timescale_con ~ normal(0,5); // biasing toward 0 = flat weights
  if (est_att == 1){
   timescale_att ~ normal(0,5);
   thresh_base ~ normal(0, 10);
  }    

  target += reduce_sum(MNL_LL_par, array_slice, splitsize, 
                     beta_ind, ind, dep, wts, start, end, task_individual,
                     task_type, timez_flip,
                     (est_scale_util ? beta_scale[1] : 1.0),  // utility scale conv --> attr
                     timescale_con, timescale_att, 
                     wts_bound_att, wts_bound_conv,
                     thresh_base,
                     love_exist, hate_exist,
                     love, hate); 
} // End Model


