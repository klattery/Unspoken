##================================================================##
#   Built using Estimation Ecosystem Functions                                        
#   Kevin Lattery Sep 2020
#
#   Combined Conversion and Attraction Module
#   Conversion is Master File, Attraction is Child
#    -> Respondents in Attraction but not Conversion will not be evaluated (bad_id)
#    -> Variables in Attraction must be subset of Conversion (otherwise FATAL ERROR)
#    -> Changed from July version: only memory management (no computations) 
##================================================================##

# Setup multi-threading
library(doParallel) # install.packages("doParallel")
cores_use <- max(detectCores() -1,1) # Number of cores to use

# Install functions/environments for estimation.  Unspoken uses MLEB within Estimation-Ecosystem
source("https://raw.githubusercontent.com/klattery/Estimation-Ecosystem/master/EE_Functions_Final.R")

##================================================================##
##       1.  SPECIFY DIRECTORIES, DATA, & COLUMNS                 ##
##       Action: Specify everything                               ##
dir_data <- "C:/Users/k.lattery/OneDrive - SKIM/Problems/Unspoken/PaulScale/F4269_Breeze"
dir_out <- "C:/Users/k.lattery/OneDrive - SKIM/Problems/Unspoken/PaulScale/F4269_Breeze/NewMethod"

# Choice/Conversion Data
choice_data <- read.csv(file.path(dir_data, "Conv.csv"), as.is=TRUE)
choice_cols <- list( # Specify the column numbers of choice_data for each of the inputs below
  id = 1,  # respondent id
  task = 3, # task number
  dep = 6, # 0 = Not swiped, 1 = swiped
  time = 7, # time (usually in milliseconds but will be standardized anyway)
  atts = 4 # ind vars that we want to estimate utilities for
)

# Swipe/Attraction Data used as Dual None
swipe_data <- read.csv(file.path(dir_data, file = "att.csv"), as.is=TRUE)
swipe_cols <- list(
  id = 1,
  task = 3,
  dep = 5, # 1 = bad swipe, 2 = positive swipe will be converted to binary below
  time = 6,  
  atts = 4
)

# change dep in swipe to 0/1
swipe_data[, swipe_cols$dep] <- swipe_data[, swipe_cols$dep] - 1 
n_atts <- length(choice_cols$atts) # define # of atts for convenience
# if there is more than 1 attribute column they must line up across the swipe and conversion.
#    1st attribute column of swipe = 1st attribute column of conversion, 2nd = 2nd, etc.
##======= End  1 ==================================================##

##================================================================##
##    2.  Merge the 2 Data Files to create stacked data_all       ##
##        TYpically just run all of this

merge_2 <- function(swipe_data, choice_data, swipe_cols, choice_cols){
  id_swipe <- unique(swipe_data[,swipe_cols$id])
  id_choice <- unique(choice_data[,choice_cols$id])
  bad_id <- id_swipe[!(id_swipe %in% id_choice)] #ids in swipe that are not in choice
  keep <- !(swipe_data[,swipe_cols$id] %in% bad_id)
  swipe_data <- swipe_data[keep,] # Remove ids from swipe that are not in conversion
  
  n_atts <- length(choice_cols$atts)
  vnames <- c("id", "task", "dep", "time", paste0("att", 1:n_atts), "None", "Source")
  att_cols <- 5:(4+ n_atts)
  
  # Align 2 Data Sets
  # Choice Data
  data_prep_c <- choice_data[,do.call(c, choice_cols)] # Choice Data
  data_prep_c$None <- 0
  data_prep_c$Source <- 1
  colnames(data_prep_c) <- vnames
  
  # Swipe data
  data_prep_s <- swipe_data[,do.call(c, swipe_cols)] 
  data_prep_s$None <- 0
  data_prep_s$Source <- 2
  colnames(data_prep_s) <- vnames
  key <- data_prep_s[,c(1,att_cols)]
  unique_items <- unique(key)
  kmatch <- match(data.frame(t(key)), data.frame(t(unique_items))) # unique id, items
  temp <- aggregate(data_prep_s, list(kmatch), mean)
  ksplit <- split(temp[,-1], temp$id)
  tasknew <- lapply(ksplit, function(x) {
    x$task <- 1:nrow(x)
    return(x)
  })
  data_prep_s <- do.call(rbind, tasknew)
  data_prep_s$task <- max(data_prep_c$task) + 100 + data_prep_s$task
  data_prep_s$concept <- 1
  
  # Now add synthetic none concept
  data_prep_s_2 <- data_prep_s
  data_prep_s_2[,att_cols] <- 0 # Attributes are 0
  data_prep_s_2$dep <- (1 - data_prep_s$dep)
  data_prep_s_2$None <- 1 
  data_prep_s_2$concept <- 2
  
  data_prep_s <- rbind(data_prep_s, data_prep_s_2)
  korder <- order(data_prep_s$id, data_prep_s$task, data_prep_s$concept)
  data_prep_s <- data_prep_s[korder,-ncol(data_prep_s)]
  data_all <- rbind(data_prep_c, data_prep_s)
  korder <- order(data_all$id, data_all$task, 1:nrow(data_all))
  data_all <- data_all[korder,] # Merged Data: Conversion (Source = 1) + Swipe (Source = 2) 

  indcode_ind <- do.call(cbind, lapply(c(1:n_atts), function(i) catcode(data_all, i+4, 1)$outcode)) # indicator coded
  check <- indcode_ind[data_all$Source == 1,]
  bad_vars <- colnames(indcode_ind)[colSums(check) == 0]
  if (length(bad_vars) > 0){
    message()
    message("!! FATAL DATA PROBLEM !!")
    message()
    message("Variable(s) found in Swipe but Not Conversion")
    message("YOU MUST REMOVE these variables from your Swipe .csv data file and start again")  
    cat(bad_vars)
  } else message("Variables matched")
  if (length(bad_id) > 0){
    message()
    message("Respondents found in Swipe that are Not in Conversion")
    message("These respondents (below) will not be part of the results")  
    cat(bad_id)
  } else message("Respondents and variables matched across Swipe and Conversion")
  .GlobalEnv$bad_vars <- bad_vars
  .GlobalEnv$bad_id <- bad_id
  return(data_all)  
}
data_all <- merge_2(swipe_data, choice_data, swipe_cols, choice_cols)
write.table(data_all, file = paste0(dir_data, "data_all.csv"), sep = ",", na = ".", row.names = FALSE)

##======= End  2 ==================================================##

##================================================================##
##       3.  SPECIFY CODING, RUN AGG MODEL                        ##
##           Default run is indicator code vs None = 0           ##
##           If so, just run all of this                          ##
indcode_specs <- lapply(c(1:n_atts),
          function(i) catcode(data_all, i+4, 1)) # indicator coded
  #==== 
  indcode_list <- make_codefiles(indcode_specs)

# Check code_master is right
# write.table(cbind(rownames(indcode_list$code_master), indcode_list$code_master), file = file.path(dir_out, "code_master.csv"), sep = ",", na = ".", row.names = FALSE)
data_conjoint <- prep_file(data_all[, 1:3], indcode_list,
                           train = TRUE) # TRUE = Use All data 
rm(indcode_specs, indcode_list) # Free RAM

# Base model for Agg Model Unconstrained
model_base <- list(
  func = list(
    pred = PredMNL,
    min = LL_Neg,
    gr = grad_MNL
  ),
  con = con_trivial(data_conjoint$P), # no constraints
  x0 =  rep(0, data_conjoint$P) # start at 0
)
agg_result <- agg_solve(data_conjoint, model_base) # Aggregate solution
agg_par <- agg_result$par
names(agg_par) <- colnames(data_conjoint$code_master)
agg_coded <- data_conjoint$code_master %*% agg_result$par
# View agg_coded to verify paramteres and results of coding

data_conjoint$agg_model$par_est <- agg_par
data_conjoint$agg_model$par_backcode <- agg_coded 
data_conjoint$agg_model$model_base <- model_base
save(data_conjoint, file = file.path(dir_out, "data_conjoint.RData"))
message(paste0("data_conjoint has " , data_conjoint$I, " respondents"))
message(paste0(" It has " , data_conjoint$P, " coded parameters"))
message(paste0(" It has ", data_conjoint$T, " total tasks, a mean of ", data_conjoint$T/data_conjoint$I, " tasks per respondent"))
  #====
##======= End  3 ==================================================##

##================================================================##
##       4.  Specify Constraints and initial x0                   ##
##           With no constraints, just run code
con_specs <- list(
  col_pos = NULL,
  col_neg = NULL,
  row_rel = NULL  # if not NULL you may have to manually specify x0
) 
# For list above, create constraints and x0 within constraints
# If you specify constraints but agg model conflicts,
# next code will try to find x0 within constraints 
  #====
# Run this to make constrain and x0
constraints <- make_con(con_specs = con_specs, code_master = data_conjoint$code_master,
                        x0_try = data_conjoint$agg_model$par_est)
constrain <- constraints$constrain 
x0 <- constraints$x0 
con_check <- check_beta(x0, constrain) # Check that x0 is within constraints
  #====
# If message says so, you must manually set x0 within constraints. Do this below: 
if (!con_check$good){
  x0 <- 0 # YOUR INPUT VECTOR THAT SATISFIES CONSTRAINTS
  con_check <- check_beta(x0, constrain)
}
##======= End  4 ==================================================##

##================================================================##
##       5.  Run MLEB                                             ##
  #====
agg_beta <- data_conjoint$agg_model$par_est
if (max(abs(x0 - agg_beta)) > .00001){ # If x0 is not unconstrained beta run again with constraints
  model_base_con <- data_conjoint$agg_model$model_base
  model_base_con$con <- as.matrix(constrain) # With constraints now
  model_base_con$x0 <- x0
  agg_result <- agg_solve(data_conjoint, model_base_con)
  agg_beta <- agg_result$par
  consistent <- FALSE
} else consistent <- TRUE

# Get Prior Covariance: cov_prior
# For big data set this may take like 30 minutes
JeffPriors <- JeffPrior(agg_beta, data_conjoint, data_conjoint$agg_model$model_base,
                        ObsFish = consistent, ExpFish = !consistent, score_n = 500)
JeffPriors[JeffPriors < .1] <- .1
if (consistent){ 
  cov_prior <- diag(sqrt(1/JeffPriors$ObsFish)) %*% cov2cor(data_conjoint$prior_cov) %*% diag(sqrt(1/JeffPriors$ObsFish))
} else cov_prior <- diag(sqrt(1/JeffPriors$ExpFish)) %*% cov2cor(data_conjoint$prior_cov) %*% diag(sqrt(1/JeffPriors$ExpFish))

# Specify Standard MLEB Model.
# x0 and alpha are normally the same and both are within your constraints
model_mleb <- list(
  func = list(
    pred = PredMNL,
    min = LL_wPriorPDF,
    gr = grad_MNLwMVN,
    logpdf = logpdf_mvnorm
  ),
  prior = list(
    alpha = agg_beta,
    cov_inv = solve(cov_prior), # Cov Inverse 
    upper_model = rep(TRUE, length(agg_beta)),
    scale = 1 # scale factor that we will solve for
  ),
  con = as.matrix(constrain), # must be matrix, 
  x0 = agg_beta # should be same as prior$alpha
)

# Setup cores and MLEB Control
mleb_control <- list(cal_resp = min(300,data_conjoint$I), cal_tasks = 1,
                     hold_resp = min(500, data_conjoint$I), hold_tasks = 2,
                     maxiter = 20, conv_n = 3,
                     tolerance = .1, hold_poss = (data_all$Source == 1),
                     dir_pdf = dir_out, solveonly = FALSE)
setup_cores(cores_use)
mleb(data_list = data_conjoint, model_list = model_mleb, mleb_control = mleb_control)

mleb_result <- model_env$mleb_result

kscale_hist <- mleb_result$kscale_hist
kiter <- 1:sum(!is.na(kscale_hist[,1]))
#plot(kiter, kscale_hist[kiter, 3])
best_fit <- max(kscale_hist[kiter, 3])
best_iter <- match(best_fit, kscale_hist[, 3]) # Picks iteration with best rlh

# Get betas, back code them, and export
best_betas <- mleb_result$eb_betas_hist[[best_iter]] #coded betas
betas_final_r <- cbind(best_betas[,1:2], best_betas[,-1:-2] %*% t(data_conjoint$code_master)) #back coded
write.table(betas_final_r, file = file.path(dir_out, "betas_notime.csv"), sep = ",", na = ".", row.names = FALSE)
write.table(mleb_result$prior_cov_hist[[best_iter]], file = file.path(dir_out, "best_cov.csv"), sep = ",", na = ".", row.names = FALSE)
write.table(mleb_result$prior_alpha_hist[[best_iter]], file = file.path(dir_out, "best_alpha.csv"), sep = ",", na = ".", row.names = FALSE)
  #====
##       Watch Console Convergence
##======= End  5 ==================================================##

##================================================================##
##       6.  Estimate Scale for Time and Solve Using Upper Level  ##
  #====
#    Add Time to Data 
# est z value of log(time) for each resp
# Make faster time bigger
if (!is.null(dev.list())) dev.off() # Reset plot area
time_clean <- function(time, data_list, too_long = 6000){
  time[(time > too_long)] <- NA # Set very long times to NA
  time <- log(time)
  # standardize times within id
  time_std_id <- do.call(rbind, tapply(X = time, INDEX = data_list$match_id , FUN = scale, simplify = FALSE))
  time_std_id[is.na(time_std_id)] <- 0  
  time_std_id <- time_std_id / -1
}
time_recode <- data_all$time
time_recode[data_all$Source == 2] <- NA
time_std_id <- time_clean(time_recode, data_conjoint)
time_std_id[time_std_id > 2,] <- 2
time_std_id[time_std_id < -2,] <- -2

# Use eb_result to create new model for time
model_time <- model_mleb
model_time$prior$alpha <- mleb_result$prior_alpha_hist[[best_iter]]
model_time$prior$cov_inv <- solve(mleb_result$prior_cov_hist[[best_iter]])
model_env <- list2env(model_time, parent = emptyenv()) # Writes over Env for EB

# Find scale of time for each respondnt
seq_test <- seq(0, 2, by = .2)
message("Estimating scale factor for each respondent's time (takes time)")
test_scale_id <- foreach(i = 1:length(data_conjoint$resp_id), .combine='cbind', .export = c('SolveID_TestWts','time_std_id')) %dopar% {
  SolveID_TestWts(i, data_conjoint, model_env, seq_test, mleb_control$hold_poss)
} 
# Above takes about 1 second per respondent with 3 cores
id_scale <- apply(test_scale_id, 2, function(coldata) best_x_fromy(seq_test, coldata))
hist(id_scale, breaks = 21)

# For each respondent apply the scale factor to get weights for each task
ksplit <- split(time_std_id, data_conjoint$match_id)
time_to_wts <- function(time_z) {
  cum_z <- pnorm(time_z)
  wts <- as.vector(cum_z / mean(cum_z))  # average to 1
  return(wts)
}
data_conjoint$wts <- do.call(c, lapply(1:length(ksplit), function(x) time_to_wts(ksplit[[x]] * id_scale[x])))
# Get final betas using weighted task
mleb(data_list = data_conjoint, model_list = model_time, mleb_control = list(solveonly = TRUE))
mleb_result2 <- model_env$mleb_result

# Get betas and back code
betas <- mleb_result2$eb_betas
betas_final_r <- as.data.frame(cbind(betas[,1:2], id_scale, betas[,-1:-2] %*% t(data_conjoint$code_master)))
betas_final_r$SynNone <- 0

# Check variables
if (length(bad_vars) > 0){
  message()
  message("Variable(s) found in Swipe but Not Conversion")
  message("Results Unreliable.  Remove these variables from Attraction and Run Again")  
  cat(bad_vars)
  betas_final_r$BAD <- "Error: Data Mismatch"
} else message("ANALYSIS DONE!  Check Output Directory for Results")
write.table(betas_final_r, file = file.path(dir_out, "betas_final_r.csv"), sep = ",", na = ".", row.names = FALSE)
  #====
##======= End  6 ==================================================##

