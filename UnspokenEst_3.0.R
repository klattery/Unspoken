#### FUNCTIONS #######################
clean_swipe <- function(att_data, att_cols, att_recode){
  dep_r <- att_recode[att_data[,att_cols$dep] + 1]
  # -2 = Hate, -1 = Bad,  0 = NA, 1 = Like, 2 = Love
  att_data[,att_cols$dep] <- dep_r
  
  att_data$swipeneg_dev <- c(-1,1,1,1,1)[dep_r + 3] # Hate loses to others
  att_data$swipepos_dev <- c(-1,-1,-1,-1,1)[dep_r + 3] # Love beats others
  att_data$dep_direction <- c(0,0,.5,1,1)[dep_r + 3] 
  
  ksplit <- split(att_data, att_data[,att_cols$id])
  scale_time <- function(kdata, att_cols){
    time <- log(kdata[,att_cols$time])
    time_bad <- (time > log(6000))
    time_lovehate <- kdata[,att_cols$dep] %in% c(-2,0,2) # -2 = Hate, 2 = LOVE
    time[time_bad | time_lovehate] <- NA
    time_z <- scale(time)
    if (sum(!is.na(time_z)) >=5) { # adjust for trend if we have 5+ obs
      task_z <- scale(1:length(time_z))
      time_z <- time_z - task_z * cor(task_z, time_z, use = "complete.obs")[1,1]
    }
    time_z_flip <- -1 * as.vector(time_z)
    time_z_flip[time_bad] <- 0 # mean value for bad times
    options(warn = -1)
    time_z_flip[time_lovehate] <- max(5, 2 * max(time_z_flip, na.rm =TRUE)) # very fast
    options(warn = 0)
    
    result <- data.frame(
      id = kdata[,att_cols$id],
      task = 9999,
      task_type = 2,
      time = kdata[,att_cols$time],
      dep_raw = kdata[,att_cols$dep],
      time_z_flip = time_z_flip,
      swipepos_dev = kdata$swipepos_dev,
      swipeneg_dev = kdata$swipeneg_dev,
      dep = kdata$dep_direction
    )
    result <- cbind(result, kdata[,att_cols$atts, drop = FALSE])
  }
  att_data_clean <- do.call(rbind, lapply(ksplit, function(x) scale_time(x, att_cols)))
  return(att_data_clean)  
}  

clean_conv <- function(conv_data, conv_cols){
  ksplit <- split(conv_data, conv_data[,conv_cols$id])
  scale_time <- function(kdata, conv_cols){
    time <- log(kdata[,conv_cols$time])
    time_bad <- (time > log(6000))
    time[time_bad] <- NA
    true_z <- scale(unique(time))
    time_z <- scale(time) 
    n_tasks <- length(unique(kdata[,conv_cols$task]))
    n_repeats <- length(time_z)/n_tasks
    n_validtimes <- sum(!is.na(time_z))/n_repeats # valid times per task
    if (n_validtimes > 1){# adjust scale for repeat time observations
      time_z <- time_z * sqrt((n_repeats * (n_validtimes - 1))/(n_repeats * n_validtimes - 1)) 
    } 
    if (n_validtimes >=5) { # adjust for trend if we have 5+ obs
      task_z <- scale(kdata[,conv_cols$task])
      time_z <- time_z - task_z * cor(task_z, time_z, use = "complete.obs")[1,1]
    }
    time_z[is.na(time_z)] <- 0
    time_z_flip <- -1 * as.vector(time_z)
    time_z_flip[time_bad] <- 0 # mean value for bad times
    
    result <- data.frame(
      id = kdata[,conv_cols$id],
      task = kdata[,conv_cols$task],
      task_type = 1,
      time = kdata[,conv_cols$time],
      dep_raw = kdata[,conv_cols$dep],
      time_z_flip = time_z_flip,
      swipepos_dev = rep(0, nrow(kdata)),
      swipeneg_dev = rep(0, nrow(kdata)),
      dep = kdata[,conv_cols$dep]
    )
    result <- cbind(result, kdata[,conv_cols$atts, drop = FALSE])
  }
  conv_data_clean <- do.call(rbind, lapply(ksplit, function(x) scale_time(x, conv_cols)))
  return(conv_data_clean)  
}  

process_utilities_unspoken <- function(data_stan, utilities, out_prefix, dir_work,
                                       task_type = NULL, inv_logit_thresh = NULL, beta_scale = 1){
  # Compute predictions
  row_weights <- data_stan$wts[data_stan$idtask_r] # convert task weights to row weights
  if (is.null(task_type)){ # Standard MNL
    pred_all <- do.call(rbind, lapply(1:data_stan$T,
                                      function(t){
                                        U <- exp(data_stan$ind[data_stan$start[t]:data_stan$end[t],] %*%
                                                   utilities[data_stan$task_individual[t],])
                                        pred <- U/sum(U)
                                        return(pred)
                                      }
    )) # lapply, rbind
  } else {
    pred_all <- do.call(rbind, lapply(1:data_stan$T,
                                      function(t){
                                        if (task_type[t] == 1){ # MNL
                                          U <- exp(data_stan$ind[data_stan$start[t]:data_stan$end[t],] %*%
                                                     utilities[data_stan$task_individual[t],] * beta_scale)
                                          pred <- U/sum(U)
                                        }
                                        if (task_type[t] == 2){ # Inv Logit
                                          U <- exp(data_stan$ind[data_stan$start[t]:data_stan$end[t],] %*%
                                                     utilities[data_stan$task_individual[t],])
                                          if (is.null(inv_logit_thresh)){
                                            pred <- U/(1+U)
                                          } else pred <- U/(exp(inv_logit_thresh[data_stan$task_individual[t]]) + U)
                                        }                                       
                                        return(pred)
                                      }
    )) # lapply, rbind                                   
  }
  utilities_r <- utilities %*% t(data_stan$code_master)
  if (!is.null(inv_logit_thresh)){
    utilities_r <- cbind(utilities_r, att_thresh = inv_logit_thresh)
  }
  
  util_name <- paste0(out_prefix,"_utilities_r.csv")
  pred_name <- paste0(out_prefix,"_preds_meanpt.csv")
  failcon_name <- paste0(out_prefix,"_utilities_failcon.csv")
  obs_vs_pred_name <- paste0(out_prefix,"_obs_vs_pred.csv")
  
  LL_id <- rowsum(log(pred_all) * data_stan$dep * row_weights, data_stan$match_id)
  sum_wts <- rowsum(data_stan$dep * row_weights, data_stan$match_id)
  sum_wts[sum_wts == 0] <- 1
  header <- data.frame(id = data_stan$resp_id, rlh = exp(LL_id/sum_wts))
  message(paste0(
    "Saving: \n",
    " Respondent mean utilities: ", util_name, "\n",
    " Predictions for data     : ", pred_name, "\n",
    " Observed vs Predicted    : ", obs_vs_pred_name 
  ))
  
  pred_all_export <- cbind(data_stan$idtask, wts = row_weights, dep = data_stan$dep, pred = pred_all)
  ind_wtask <- apply(data_stan$ind_levels, 2 , function(x){
    result <- paste(data_stan$other_data$task_type, formatC(x, width=3, flag="0"), sep="_")
  })
  colnames(ind_wtask) <- colnames(data_stan$ind_levels)
  obs_vs_pred <- obs_vs_pred(pred_all_export[,3:5], ind_wtask)
  
  write.table(cbind(header, utilities_r), file = file.path(dir_work, util_name), sep = ",", na = ".", row.names = FALSE)
  write.table(pred_all_export, file = file.path(dir_work, pred_name), sep = ",", na = ".", row.names = FALSE)
  write.table(obs_vs_pred, file = file.path(dir_work, obs_vs_pred_name), sep = ",", na = ".", row.names = FALSE)
  
  # Check if utilities meet constraints
  con_matrix <- diag(data_stan$con_sign)
  con_matrix <- rbind(con_matrix[rowSums(con_matrix !=0) > 0,,drop = FALSE], data_stan$paircon_matrix)
  bad_ids <- rowSums(((utilities %*% t(con_matrix)) < 0)) > 0
  if (sum(bad_ids) > 0){
    message(paste0(sum(bad_ids), " Respondents had reversals from constraints.\n",
                   "Reversals saved to: ", failcon_name))
    write.table(cbind(header, utilities_r)[bad_ids,], file = file.path(dir_work, failcon_name), sep = ",", na = ".", row.names = FALSE)
  } else message(" All respondent mean utilities obey constraints")
}

##########  END FUNCTIONS  ##################

# Clean data and merge with each other
if (!is.null(att_data)){
  est_att <- TRUE
  att_data_clean <- clean_swipe(att_data, att_cols, att_recode)
} else{
  est_att <- FALSE
  att_data_clean <- NULL
} 
if (!is.null(conv_data)){
  est_conv <- TRUE
  conv_data_clean <- clean_conv(conv_data, conv_cols)
} else{
  est_conv <- FALSE
  conv_data_clean <- NULL
} 

data_all <- rbind(conv_data_clean, att_data_clean)
korder <- order(data_all$id, data_all$task, 1:nrow(data_all))
data_all <- data_all[korder,] # Merged Data: Conversion (Source = 1) + Swipe (Source = 2) 

ind_names <- colnames(data_all)[-1:-9]
indcode_spec <- list()
for (i in seq_len(length(ind_names))){
  indcode_spec[[i]] <- catcode(data_all, ind_names[i], 3) # effects code  
}


################################################ 
## 3. RUN: CODE AND PREPARE DATA FOR STAN
#### CHECK: PRINTED OUTPUT LOOKS RIGHT
col_id_task_dep <- c(1,2,9) # Columns for id, task, dep

# 3) Code and Prepare data_stan =============
if (!is.null(indcode_spec)){
  indcode_list <- make_codefiles(indcode_spec) # Combine specifications above into one list
  # save_codemastercon(indcode_list,dir$work, out_prefix) # Save combined code_master and constraints 
  data_stan <- prep_file_stan(idtaskdep = data_all[,col_id_task_dep],
                              indcode_list = indcode_list,
                              data_cov = NULL, specs_cov_coding = NULL,
                              check_collinearity = FALSE,
                              other_data = data_all[,colnames(data_all) %in% 
                                                      c("task_type","time_z_flip","swipepos_dev","swipeneg_dev")],
                              force_sumtask_dep1 = FALSE)
} else cat("!! STOP NOW and FIX NAMES that do not match !!")


################################################### 
## 5. OPTIONALLY CHANGE DEFAULTS
data_model <- c(data_model,
  list(
    con_use = 1, # 0 = ignore constraints, 1 = use
    con_factors = c(mult = .25, bound = 1.5), # Default (.25, 1.5) must be >0
    threads_per_chain = data_stan$threads_rec,
    agg_model = NULL, tag = NULL, ind = NULL, resp_id = NULL, other_data = NULL,
    task_type = data_stan$other_data$task_type[data_stan$start],
    timez_flip = data_stan$other_data$time_z_flip,
    thresh_high_dev = data_stan$other_data$swipepos_dev,
    thresh_low_dev = data_stan$other_data$swipeneg_dev    
  )
)
data_model$splitsize <- round(.5 + data_stan$T/(4 * data_model$threads_per_chain))

##################################
# ESTIMATE MODEL =============
#####  Specify Stan Model 
stan_file <- "Unspoken_Final_v1.5.stan" # Name of stan model in dir$stanmodel
stan_outname <- paste0(out_prefix, "_StanOut_", 
                       format(Sys.time(), '%Y%m%d-%H%M%S')) # Base Name of Stan Output files  
rm(indcode_spec); rm(indcode_list); gc() # Clear RAM

#####  Compile and RunStan Model 
HB_model <- cmdstan_model(file.path(dir$stanmodel,stan_file), quiet = TRUE,
                          cpp_options = list(stan_threads = TRUE))
time_start <- format(Sys.time(), '%Y%m%d-%H%M%S')
stan_outname <- paste0(out_prefix, "_StanOut_", time_start)   
out_folder <- paste0(out_prefix, "_", time_start)
dir_run <- create_tempdir(dir, out_folder, save_specs = FALSE, code_master = data_stan$code_master) 
HB_fit <- HB_model$sample(modifyList(data_stan, data_model),
                          iter_warmup = data_model$iter_warmup,
                          iter_sampling = data_model$iter_sampling,
                          output_dir = dir$stanmodel,
                          output_basename = stan_outname, # set stan_outname above if changing
                          chains = 2,
                          parallel_chains = 2,
                          threads_per_chain = data_model$threads_per_chain,
                          save_warmup = TRUE,
                          refresh = 10,
                          adapt_delta = .8,
                          seed = 123,
                          init = .1,
                          show_messages = FALSE,
                          diagnostics = NULL)
saveRDS(HB_fit, file.path(dir_run, paste0(out_prefix, "_HB_fit.rds")))

#####  Check Convergence, Export Files  
if (min(HB_fit$return_codes() == 0)){
  cat("Reading draws from Stan csv output into R (large files take time)...")
  draws_beta <- read_cmdstan_csv(HB_fit$output_files(), variables = "beta_ind", format = "draws_list")
  draws_beta$warmup_draws <- NULL # to save space
  utilities <- get_mean_beta(draws_beta) # computes mean point estimates from draws
  checkconverge_export(draws_beta, vnames = colnames(data_stan$code_master), out_prefix, dir_run, data_model$export_draws)
  if (est_att){
    draws_thresh_base <- read_cmdstan_csv(HB_fit$output_files(), variables = "thresh_base", format = "draws_list")
    mychains <- lapply(draws_thresh_base$post_warmup_draws, function(x) do.call(cbind, x))
    fit_thresh_base <- plot_draws_df(mychains,  pdf_path = file.path(dir_run, paste0(out_prefix,"_thresh_base.pdf")))
    att_thresh <- fit_thresh_base$mean
  } else att_thresh <- NULL
  if (est_att & est_conv){
    draws_beta_scale <- read_cmdstan_csv(HB_fit$output_files(), variables = "beta_scale", format = "draws_list")
    mychains <- lapply(draws_beta_scale$post_warmup_draws, function(x) do.call(cbind, x))
    fit_beta_scale <- plot_draws_df(mychains, pdf_path = file.path(dir_run, paste0(out_prefix, "_beta_scale.pdf")))
    write.table(fit_beta_scale, file = file.path(dir_run, paste0(out_prefix,"_beta_scale.csv")), sep = ",", na = ".", row.names = FALSE)
    beta_scale <- fit_beta_scale$mean
  } else beta_scale <- 1
  process_utilities_unspoken(data_stan, utilities, out_prefix, dir_run,
                    task_type = data_stan$other_data$task_type[data_stan$start],
                    inv_logit_thresh = att_thresh,
                    beta_scale = beta_scale)
  if (data_stan$P_cov > 0) covariate_means <- get_covariates(HB_fit, data_stan, dir_run)
  zip_allout(dir_run, dir$work, paste0(out_folder,".zip"))
} else message("Stan Estimation Did not Finish")