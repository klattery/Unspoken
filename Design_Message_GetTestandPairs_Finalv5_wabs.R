####      Generate Experimental Design of Pairwise Items      ###########
###       Written April 2016 Kevin Lattery, Update 8/2018     ##########
###       Inputs are Lines Immediately Below Between ##     ###########     
################################################################
# SPECIFY PARAMETERS HERE#
dir <- "C:/Users/k.lattery/OneDrive - SKIM/Projects/Google/MD_Design/"

set_test <- 1:17 # vector of item numbers to be used as test items
ntest_perver <- 7     # How many test items are used for each specific version

set_comp <- 18:25 # vector of item numbers to be used as comp items, always shown
                  # set_comp <- NULL if no items are always shown

show_eachitem <- 2.5 # How many times to show each item, can be decimal like 2.5

# If no restrictions on items shown in a versions then use:
must_haves <- as.data.frame(c(1:100)) # last number is # of versions
             #as.data.frame(c(1:number of versions))

# If restrictions on items shown in versions then use:
#must_haves <- read.csv(paste0(dir, "MustHave_Kellogs.csv"))
# .csv has version item1... itemn with item numbers needed
# Must have a row for each version

con_pairs <- NULL
#  list(c(12,53), c(20,21), c(26,27)) # constraints
# set con_pairs <- NULL if no constraints

#######################################################################
####   Just Run this Code to Create Defaults and Helper Variables  ####
#######################################################################
items_task <- 2 # How many items in a task (for Unspoken this is 2 and cannot be changed)

# Check that these are what you expect
n_ver <- nrow(must_haves) # number of versions to create
n_items <- ntest_perver + length(set_comp) # number of items shown in a version
n_tasks <-  round(n_items * show_eachitem/items_task,0) # number of tasks in a version, assume each task is a pair of items = 2
# n_tasks can be changed manually if needed

set_test_r <- 1:length(set_test) # recode test values to 1:n

#Create default target and target_wt for test items shown n x n where n = number of test items
target_freq <- (n_ver*ntest_perver)/length(set_test_r)
target_cov_num <- choose(ntest_perver, 2) # Pairs in version
target_cov_den <- (length(set_test_r)^2 - length(set_test_r))/2 # Total Pairs of target
target_cov <- n_ver * target_cov_num/target_cov_den
target_1 <- matrix(target_cov, nrow = length(set_test_r), ncol = length(set_test_r))
diag(target_1) <- target_freq
target_1_wt <- matrix(1, nrow = nrow(target_1), ncol=ncol(target_1))
diag(target_1_wt) <- (nrow(target_1)-1) * 1 

# Define Functions
ind_code <- function(vec_val, klength){
  vec0 <- rep(0, klength)
  vec0[vec_val] <- 1
  return(vec0)
}
gen_design <- function(dummy, design_in, npicks, target, target_wt){
  design <- design_in
  for (iter in 1:10){
    for (i in sample(nrow(design))){
      rowpick <- design[i,]
      musthave <- design_in[i,]
      design_minus <- design[-i,]
      show_minus <- t(design_minus) %*% design_minus
      need_minus <- target_wt * (target - show_minus)
      test <- do.call(rbind, lapply(1:2000, function(x) sample(x=ncol(design), size=npicks - sum(musthave), prob = 1-musthave)))
      test_ind <- t(apply(test,1,ind_code, ncol(design)))
      test_ind <- test_ind + do.call(rbind,lapply(1:nrow(test_ind), function(x) musthave))
      if (iter > 1){
        test_ind <- rbind(rowpick,test_ind)
      }
      kfit <- apply(test_ind %*% need_minus, 1, sum)
      fitbest <- which(kfit == max(kfit))
      pick <- sample(fitbest, 1)
      rowpick <- test_ind[pick,] # add new pick
      design[i,] <- rowpick
    }
  }
  return(design)
}  

compfit <- function (pdesign, target, target_wt){
  show_now <- t(pdesign) %*% pdesign
  return(sum(target_wt * abs(show_now - target)))
}
gen_design_pair <- function(task_code, n_tasks, target_en, target_wt){
  # pick 2 random tasks
  pick <- sample(1:nrow(task_code),2)
  design <- task_code[pick,]
  task_code <- task_code[-pick,] #drop pick
  
  # now loop through adding picks
  for (i in 1:(n_tasks - nrow(design))){
    show_now <- t(design) %*% design
    need <- target_wt * (target_en - show_now)
    fit_ver <- apply(task_code %*% need, 1, sum)
    vmax <- max(fit_ver)
    vbest <- (fit_ver >= (vmax - .00001)) #rows that are best (or nearly)
    pick <- sample(which(vbest),1)
    design <- rbind(design, task_code[pick,])
    task_code <- task_code[-pick,] #drop pick      
  }
  return(design)
}


fix_con <- function(x, con_pairs){
  # Fix constraints, if con_pairs is null does nothing
  check_bad <- function(x, con_pairs){
    bad <- duplicated(x)
    if (!is.null(con_pairs)){
      kmatch <- sapply(con_pairs, function(pair) rowSums(x[,pair]))
      bad <- bad | (rowSums(kmatch == 2) > 0)
    }
    return(bad)
  }
  bad <- check_bad(x, con_pairs)
  while (sum(bad) > 0){
    if (sum(bad) == nrow(x)) bad[1] <- FALSE # Can't have all bad
    x <- rbind(x[bad,], x[!bad,])
    for (i in 1:sum(bad)){
      badrow <- which(x[i,] == 1)
      subrow <- sample((i+1):nrow(x), 1)
      nextrow <-  which(x[subrow,] == 1)
      if (length(intersect(badrow, nextrow)) == 0) {
        x[i,badrow[2]] <- 0
        x[i,nextrow[2]] <- 1
        x[subrow,nextrow[2]] <- 0
        x[subrow,badrow[2]] <- 1
      }
    }
    bad <- check_bad(x, con_pairs)
  } # end while 
  return(x)
}


#########################################################
#                   STAGE 1                             #
#    Get Items to Be Tested: design_itemsall            #
#########################################################
# Check target (target 1-way and 2-way) and target_wt (how important to hit target)
# Adjust them if desired
diag(target_1_wt) <- (nrow(target_1)-1) * 2 # *2 or *5 for better one way
designs_stage1 <- 10 # 10 is reasonable, but can go higher if needed (30 takes time)

# Generate design_itemsall
design <- matrix(0, nrow = n_ver, ncol = length(set_test_r))
if ((ncol(must_haves)) > 1){
  must_haves_r <- must_haves[,-1] # drop first column
  must_haves_r <- apply(must_haves_r,2,match,set_test) # recode so this matches 1:n values of set_test_r 
  for (i in 1:nrow(design)){
    kvalue <- must_haves_r[i,]
    design[i,kvalue] <- 1
  } # design has initial design based on must_haves, but indicator coded
}
design_in <- design # Keep copy of initial design

# Next step TAKES TIME
multdesigns <- lapply(1:designs_stage1, gen_design,design_in, ntest_perver, target_1, target_1_wt) 
fit_all <- sapply(multdesigns, compfit, target_1, target_1_wt) # fit of each member above
design_testpick <- multdesigns[[which.min(fit_all)]]  # design of test items to show for each version (row)
summary_stat <- t(design_testpick) %*% design_testpick # shows frequency (diagnol) and cross tab of items shown together
design_itemsall <- cbind(design_testpick, matrix(1, nrow=nrow(design_testpick), ncol = length(set_comp)))
vnames <- paste0("Test", set_test)
if (!is.null(set_comp)) vnames <- c(vnames, paste0("Comp", set_comp))
colnames(design_itemsall) <- vnames 

# Export list of items
write.table(design_itemsall, file = paste0(dir,"design_itemsall.csv"), sep = ",", na = ".", row.names = FALSE) 

#########################################################################
#                   STAGE 2                                             #
#    Generate  pairwise tasks using fixed number of items n_items       #
#    Creates a list "designs" where each element is a design            #
#########################################################################
designs_stage2 <- 10000 # Number of Designs to create: 10,000; higher for more sparse

#a) Create task_code as tasks to choose from (indicator coding)
if (items_task > 2) { # Added 3/2019
  vec_0 <- rep(0, n_items)
  task_con <- do.call(rbind, lapply(1:100000, function(x){
    picks <- sample(1:n_items, items_task)
    result <- vec_0
    result[picks] <- 1
    return(result)
  }))
  task_code <- unique(task_con)
} else {
  task_pop <- expand.grid(1:n_items, 1:n_items) # Generate possible tasks - can be large sample
  tdrop <- (task_pop[,1] == task_pop[,2]) # picking same item twice is not allowed
  task_con <- task_pop[!tdrop,] # filter out bad tasks
  task_code <- diag(n_items)[task_con[,1],]  + diag(n_items)[task_con[,2],] # indicator coded
}

#b) Create target matrix for version
target_temp <- t(task_code) %*% task_code # empirical target sampling from
target_2 <- target_temp * ((n_tasks * items_task)/sum(diag(target_temp))) # natural target per version
target_2_wt <- matrix(1, nrow = nrow(target_2), ncol=ncol(target_2))
diag(target_2_wt) <- nrow(target_2) * 5 # Can be adjusted higher to fit one-way

#c) Find good n attribute designs (pairwise)
designs <- list() # designs will be stored here
stats_all <- NULL # summary stats for each design
if (nrow(task_code) > 10000) task_code <- task_code[sample(1:nrow(task_code), 10000),]
for (i in 1:designs_stage2){ 
  test <- gen_design_pair(task_code, n_tasks, target_2, target_2_wt) # Fixed 3/2019
  designs[[i]] <- test
  twoway <- t(test) %*% test
  kdet <- det(twoway)^(1/ncol(twoway))
  oneway <- diag(twoway)
  kstd <- sd(oneway)
  diag(twoway) <- 0
  max2way <- max(twoway)
  stats <- c(oneway, kstd, max2way, kdet)
  stats_all <- rbind(stats_all, stats)
}
#stats_all has one way freq, std of that, max 2 way, d-opt value
colnames(stats_all) <- c(paste0("OneWay", 1:n_items), "Std","max2way", "DOpt")
bad_dopt <- as.matrix(is.na(stats_all[,ncol(stats_all)]))
stats_export <- cbind(1:nrow(stats_all), stats_all)
colnames(stats_export)[1] <- "design_num"
write.table(stats_export, file = paste0(dir,"stats_all.csv"), sep = ",", na = ".", row.names = FALSE) 


###########################################################
#          STOP: MANUAL WORK  SECTION                     #
#    Creates "good_des" showing which designs are best    #
###########################################################

# Method 1
# Directory has .csv file "stats_all".
# Use Exel to keep desired rows and import as "good_des"
good_des <- read.csv(paste0(dir, "good_des.csv"))

# Method 2
# Create good_des in R  
# Some sample code below, but will depend upon your design
check_sd_one <- as.matrix(tapply(stats_all[,n_items+1], stats_all[,n_items+1], length))
check_max_two <- as.matrix(tapply(stats_all[,n_items+2], stats_all[,n_items+2], length))
#Use check files here to set criteria for bal_one and bal_two

  bal_one <- as.matrix(stats_all[,n_items+1] < .6) ### MANUAL ### sd of oneway criteria
  bal_two <- as.matrix(stats_all[,n_items+2] <= 1) ### MANUAL ### max of twoway criteria
  sum(bal_one * bal_two *!bad_dopt) # check count of designs meeting criteria above

check <- cbind(1:nrow(stats_all),stats_all)[bal_one & bal_two &!bad_dopt,] 
d_opt <- as.matrix(check[,ncol(check)])
khist <- hist(d_opt, breaks = 30) # histogram to check for outliers from d-optimal
khist <- hist(d_opt[(d_opt> .25)], breaks = 30) # Specify d_opt > x to get more refined estimates
# Use hist plot to set d-optimal criteria 
max(d_opt)
quantile(d_opt, .95)
good_des <- check[(check[,ncol(check)] > 1.6),] ### MANUAL ### cutoff from histogram
# first column of good_des has the design numbers we want to keep
# Other columns are not relevant

#######################################################
#         End of MANUAL WORK  SECTION                 #
#######################################################


#################################################################################
#                      STAGE 3                                                  #
#    Select Stage 2 Binary Tasks, Recoding the items to match Stage 1           #
#################################################################################
  
# Get new target based on test + competitive items
get_target <- function(x){
  design_new <- NULL
  for (i in 1:nrow(design_itemsall)){
    items_in <- (design_itemsall[i,]) 
    items_in2 <- items_in * cumsum(items_in)
    des_pick_n <- sample(good_des[,1],1) # random good design num
    des_pick <- designs[[des_pick_n]] # pairwise design chosen
    con_keep <- matrix(0, nrow = nrow(des_pick), ncol = ncol(design_itemsall))
    for (j in 1:nrow(des_pick)){
      pair_rel <- which(des_pick[j,] == 1)
      pair_abs <- items_in * (items_in2 %in% pair_rel)
      con_keep[j,] <- pair_abs
    }
    con_keep <- cbind(i, con_keep)
    design_new <- rbind(design_new, con_keep)
  }
  test <- design_new[,-1]
  check_test <- t(test) %*% test 
  target <- check_test
  return(target)
}
# Simple target without constraints
target_3 <- get_target(1)
sec1 <- 1:length(set_test)
diag(target_3)[sec1] <- mean(diag(target_3)[sec1])
if (length(set_comp) > 0) {
  sec2 <- (length(set_test) + 1):(length(set_test) + length(set_comp))
  diag(target_3)[sec2] <- mean(diag(target_3)[sec2])
}

keepdiag <- diag(target_3)
diag(target_3) <- NA
target_3[sec1, sec1] <- mean(target_3[sec1, sec1], na.rm = TRUE)
if (length(set_comp) > 0) {
  target_3[sec2, sec1] <- mean(target_3[sec2, sec1], na.rm = TRUE)
  target_3[sec1, sec2] <- mean(target_3[sec1, sec2], na.rm = TRUE)
  target_3[sec2, sec2] <- mean(target_3[sec2, sec2], na.rm = TRUE)
}  
diag(target_3) <- keepdiag

# More complex target with constraints/etc estimated empirically 
#emp_targets <- lapply(1:500, get_target) # Empirical execution
#target_3 <- Reduce("+", emp_targets)/length(emp_targets)

target_3_wt <- matrix(1, nrow =nrow(target_3), ncol=ncol(target_3))
diag(target_3_wt) <- nrow(target_3_wt) * 3 # can be different than all 1

x1 <- do.call(rbind, con_pairs) # No Constraints - skip 3 lines
target_3[x1] <- 0
target_3[cbind(x1[,2], x1[,1])] <- 0

# Assemble designs based on target
design_new <- NULL
for (i in 1:2){ 
  items_in <- (design_itemsall[i,]) 
  items_in2 <- items_in * cumsum(items_in)
  des_pick_n <- sample(good_des[,1],1) # random good design num
  des_pick <- designs[[des_pick_n]] # pairwise design chosen
  con_keep <- matrix(0, nrow = nrow(des_pick), ncol = ncol(design_itemsall))
  for (j in 1:nrow(des_pick)){
    pair_rel <- which(des_pick[j,] == 1)
    pair_abs <- items_in * (items_in2 %in% pair_rel)
    con_keep[j,] <- pair_abs
  }
  con_keep <- fix_con(con_keep, con_pairs)
  con_keep <- cbind(i, con_keep)
  design_new <- rbind(design_new, con_keep)
} #end of initial seeds

for (i in 3:nrow(design_itemsall)){
  items_in <- (design_itemsall[i,]) 
  items_in2 <- items_in * cumsum(items_in)
  kcounts <- t(design_new[,-1]) %*% design_new[,-1]
  diff <- target_3 - kcounts
  bestimp <- Inf
  for (ksamp in 1:nrow(good_des)){
    des_pick_n <- good_des[ksamp,1] # design number
    des_pick <- designs[[des_pick_n]] # pairwise design chosen
    con_keep <- matrix(0, nrow = nrow(des_pick), ncol = ncol(design_itemsall))
    for (j in 1:nrow(des_pick)){
      pair_rel <- which(des_pick[j,] == 1)
      pair_abs <- items_in * (items_in2 %in% pair_rel)
      con_keep[j,] <- pair_abs
    }
    con_keep <- fix_con(con_keep, con_pairs)
    con_counts <- t(con_keep) %*% con_keep
    con_diff <- (diff - con_counts) * target_3_wt
    con_fit <- sum(con_diff[con_diff > 0]) - .5 * sum(con_diff[con_diff < 0])  
    if (con_fit < bestimp) {
      bestimp <- con_fit
      winner <- con_keep
    }
  } #end ksamp
  design_new <- rbind(design_new, cbind(i,winner))
  print(paste0("Version " ,i))
 }
 # design_new contains initial design from one loop through
 tab_design <- t(design_new[,-1]) %*% design_new[,-1] # check one-way and two-way

colnames(design_new) <- c("Version", vnames)
colnames(tab_design) <- vnames

# Export Results
write.table(design_new, file = paste0(dir,"DesignNew.csv"), sep = ",", na = ".", row.names = FALSE)
write.table(tab_design, file = paste0(dir,"Design_Tab.csv"), sep = ",", na = ".", row.names = FALSE)
items <- t(sapply(1:nrow(design_new),function(i) which(1 == design_new[i,-1])))
write.table(items, file = paste0(dir,"DesignNew_Items.csv"), sep = ",", na = ".", row.names = FALSE)
write.table(target_3, file = paste0(dir,"target_3.csv"), sep = ",", na = ".", row.names = FALSE)

sort(unique(as.vector(items)))
swipe <- unique(items)
####################################################################
#       design_new has design                                      #
#       tab_design has the one-way and two-way freq                #
####################################################################

sum(duplicated(design_new)) # Check duplicates
