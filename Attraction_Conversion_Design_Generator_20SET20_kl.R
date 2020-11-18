conversion_function <- function(ntest, ntest_perver, ntest_comp, show_eachitem, show_eachitem_attraction, n_versions, restrictions_table, constraints_table, shiny = TRUE) {
  if(shiny) {
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = 'Calculating', value = 0)
  }
  ################################################################
  ####      Generate Experimental Design of Pairwise Items      ###########
  ###       Written April 2016 Kevin Lattery, Update 8/2018     ##########
  ###       Updated 9/20 by Eugenio Grant        ##########
  ###       Inputs are Lines Immediately Below Between ##     ########### 
  ###       Update 11/20 Kevin Lattery
  ################################################################
  
  #Download libraries necessary for the code 
  
  # SPECIFY PARAMETERS HERE#
  set_test <- 1:ntest # vector of item numbers to be used as test items
  ntest_perver <- ntest_perver # How many test items are used for each specific version
  
  if (ntest_comp == 0) {
    set_comp <- NULL
  } else {
    set_comp <- ((ntest+1):(ntest+ntest_comp))
  }

  show_eachitem <- show_eachitem # How many times to show each item in conversion, can be decimal like 2.5
  show_eachitem_attraction <- show_eachitem_attraction # How many times to show each item in Attraction, cannot be decimal
  
  if (is.null(restrictions_table)) {
    must_haves <- as.data.frame(c(1:n_versions))
  } else {
    must_haves <- as.data.frame(restrictions_table)
  }
  
  
  if (is.null(constraints_table)) {
    con_pairs <- NULL
  } else {
    df <- as.data.frame(constraints_table)
    newlist <- lapply(split(df, row.names(df)), unlist)
    con_pairs <- lapply(newlist, function(x) x[!is.na(x)])
  }
  
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
        test <- do.call(rbind, lapply(1:2000, function(x) sample(x=ncol(design), size=npicks, prob = 1-design_in[i,])))
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
  designs_stage1 <- 30 # 10 is reasonable, but can go higher if needed (30 takes time)
  
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
  if(shiny) progress$inc(.2, message = "1. Creating Design for Items Shown Each Version")
  if(!shiny) message("1. Creating Design for Items Shown Each Version")
  multdesigns <- lapply(1:designs_stage1, gen_design,design_in, ntest_perver, target_1, target_1_wt) 
  fit_all <- sapply(multdesigns, compfit, target_1, target_1_wt) # fit of each member above
  design_testpick <- multdesigns[[which.min(fit_all)]]  # design of test items to show for each version (row)
  summary_stat <- t(design_testpick) %*% design_testpick # shows frequency (diagnol) and cross tab of items shown together
  design_itemsall <- cbind(design_testpick, matrix(1, nrow=nrow(design_testpick), ncol = length(set_comp)))
  vnames <- paste0("Test", set_test)
  if (!is.null(set_comp)) vnames <- c(vnames, paste0("Comp", set_comp))
  colnames(design_itemsall) <- vnames 
  if(shiny) progress$inc(0.6, message = "2. Creating Designs")
  if(!shiny) message("2. Creating Designs")
  
  # Export list of items
  # write.table(design_itemsall, file = paste0(dir,"design_itemsall.csv"), sep = ",", na = ".", row.names = FALSE) 
  
  #########################################################################
  #                   STAGE 2                                             #
  #    Generate  pairwise tasks using fixed number of items n_items       #
  #    Creates a list "designs" where each element is a design            #
  #########################################################################
  designs_stage2 <- 10000 # Number of Designs to create: 10,000; higher for more sparse
  
  #a) Create task_code as tasks to choose from (indicator coding)
  task_pop <- expand.grid(1:n_items, 1:n_items) # Generate possible tasks - can be large sample
  tdrop <- (task_pop[,1] == task_pop[,2]) # picking same item twice is not allowed
  task_con <- task_pop[!tdrop,] # filter out bad tasks
  task_code <- diag(n_items)[task_con[,1],]  + diag(n_items)[task_con[,2],] # indicator coded
  
  #b) Create target matrix for version
  target_temp <- t(task_code) %*% task_code # empirical target sampling from
  target_2 <- target_temp * ((n_tasks * items_task)/sum(diag(target_temp))) # natural target per version
  target_2_wt <- matrix(1, nrow = nrow(target_2), ncol=ncol(target_2))
  diag(target_2_wt) <- nrow(target_2) * 5 # Can be adjusted higher to fit one-way
  
  #c) Find good n attribute designs (pairwise)
  designs <- list() # designs will be stored here
  stats_all <- NULL # summary stats for each design
  for (i in 1:designs_stage2){ 
    test <- gen_design_pair(task_code, n_tasks, target_2, target_2_wt)
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
  # write.table(stats_export, file = paste0(dir,"stats_all.csv"), sep = ",", na = ".", row.names = FALSE) 

  ###########################################################
  #          STOP: MANUAL WORK  SECTION                     #
  #    Creates "good_des" showing which designs are best    #
  ###########################################################
  
  # Method 1
  # Directory has .csv file "stats_all".
  # Use Exel to keep desired rows and import as "good_des"
  # good_des <- read.csv(paste0(dir, "good_des.csv"))
  
  # Method 2
  # Create good_des in R  
  # Some sample code below, but will depend upon your design
  check_sd_one <- as.matrix(tapply(stats_all[,n_items+1], stats_all[,n_items+1], length))
  check_max_two <- as.matrix(tapply(stats_all[,n_items+2], stats_all[,n_items+2], length))
  #Use check files here to set criteria for bal_one and bal_two
  
    bal_one <- as.matrix(stats_all[,n_items+1] <= quantile(stats_all[,n_items+1], .5)) ### MANUAL ### sd of oneway criteria
    bal_two <- as.matrix(stats_all[,n_items+2] <= quantile(stats_all[,n_items+2], .5)) ### MANUAL ### max of twoway criteria
    sum(bal_one * bal_two *!bad_dopt) # check count of designs meeting criteria above
    check <- cbind(1:nrow(stats_all),stats_all)[bal_one & bal_two &!bad_dopt,] 
    d_opt <- as.matrix(check[,ncol(check)])
   good_des <- check[(check[,ncol(check)] >= quantile(d_opt, .25)),] 
   if(shiny) progress$inc(.7, message = paste0("3. Selecting from ", nrow(good_des), " Designs for Versions"))
   if(!shiny) message(paste0("3. Selecting from ", nrow(good_des), " Designs for Versions"))   
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
  
  if (is.list(con_pairs)) {
    x1 <- do.call(rbind, con_pairs)
    target_3[x1] <- 0
    target_3[cbind(x1[,2], x1[,1])] <- 0
  }
  
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
  # write.table(design_new, file = paste0(dir,"DesignNew.csv"), sep = ",", na = ".", row.names = FALSE)
  # write.table(tab_design, file = paste0(dir,"Design_Tab.csv"), sep = ",", na = ".", row.names = FALSE)
  items <- t(sapply(1:nrow(design_new),function(i) which(1 == design_new[i,-1])))
  # write.table(items, file = paste0(dir,"DesignNew_Items.csv"), sep = ",", na = ".", row.names = FALSE)
  sort(unique(as.vector(items)))
  swipe <- unique(items)
  ####################################################################
  #       design_new has design                                      #
  #       tab_design has the one-way and two-way freq                #
  ####################################################################
  
  sum(duplicated(design_new)) # Check duplicates
  
  #####                                                                                 #####
  #####                          Unspoken design generation                             #####
  #####                                                                                 #####
  #####      Target output 1) Attraction design matrix 2) Conversion design matrix      #####
  
  
  #1. change the data frame format to edit the output from original Unspoken codes
  
  items <- as_tibble(items)
  design_new <- as_tibble(design_new)
  tab_design <- as_tibble(tab_design)
  design_itemsall <- as_tibble(design_itemsall)
  stats_export <- as_tibble(stats_export)
  
  
  #2. change conversion output format to designed format
  
  items_u1 <- items %>%
    rename('1'= 1, '2'=2)%>%    
    mutate(obs=1:nrow(items))%>%
    mutate(version= rep(1:(nrow(items)/(n_tasks)), each=n_tasks))
  
  # (Make a col for tasks)
  
  b=n_tasks
  c=nrow(items)/n_tasks
  a=1:b
  task <- as.data.frame(rep(a,c))
  
  items_u2 <- cbind(items_u1, task)
  colnames(items_u2)[5] <- "task"
  
  items_u3 <- items_u2 %>%
    gather('order', 'concept', 1:2)%>%
    
    arrange(obs)%>%
    transmute(version, task, order, concept)%>%
    select(version, task, concept, order)
  
  items_u3 <- items_u3 %>%
    mutate(product=0)
    
  items_u3 <- items_u3[,c(1,2,3,5,4)]
  
  
  # write.table(items_u3, file = paste0(dir,"ConversionDesign.csv"), sep = ",", na = ".", row.names = FALSE) 
  if(shiny) progress$inc(.8, message = "4. Finished Conversion")
  if(!shiny) message("4. Finished Conversion")
  return(items_u3)
} 




attraction_function_1 <- function(items_u3, ntest, ntest_comp, ntest_perver, show_eachitem_attraction) { 
  set_test <- 1:ntest
  
  if (ntest_comp == 0) {
    set_comp <- NULL
  } else {
    set_comp <- ((ntest+1):(ntest+ntest_comp))
  }
  
  #3. Change attraction output format to designed formats
  
  sumi3<- items_u3 %>%
    group_by(version)%>%
    summarize(ncon = n_distinct(concept))
  
  nconcept_a<- mean(sumi3$ncon) 
  
  items_u4 <- items_u3 %>%
    transmute(version, concept)%>%
    group_by(version)%>%
    distinct(concept)
  
  items_u4 <- as_tibble(items_u4)
  
  rand1=as.data.frame(runif(1:nrow(items_u4)))
  colnames(rand1)[1] <- "rand"    
  
  rand2=as.data.frame(runif(1:nrow(items_u4)))
  colnames(rand2)[1] <- "rand" 
  
  items_u5a <- cbind(items_u4, rand1)
  items_u5b <- cbind(items_u4, rand2)
  
  items_u6a <- items_u5a %>%
    arrange(version, rand)
  items_u6b <- items_u5b %>%
    arrange(version, rand)
  
  if(show_eachitem_attraction==3){
    rand3=as.data.frame(runif(1:nrow(items_u4)))
    colnames(rand3)[1] <- "rand"
    items_u5c <- cbind(items_u4, rand3)
    items_u6c <- items_u5c %>%
      arrange(version, rand)
  }
  
  # (make a new col for attention test question)
  
  nversion_a <- (max(items_u6a$version))
  nconcept_a <- (ntest_perver+length(set_comp))
  max_con <- max(set_test,set_comp)
  att_n_tasks = show_eachitem_attraction*nconcept_a
  
  attt1 <- tibble(version=(1:nversion_a), concept=max_con+1, rand=runif(1:nversion_a))
  attt2 <- tibble(version=(1:nversion_a), concept=max_con+2, rand=runif(1:nversion_a))
  attt3 <- tibble(version=(1:nversion_a), concept=max_con+1, rand=runif(1:nversion_a))
  
  # (merge)
  items_u7a <- rbind(items_u6a, attt1) %>%
    arrange(version, rand)
  items_u7b <- rbind(items_u6b, attt2) %>%
    arrange(version, rand)
  if(show_eachitem_attraction==3){
    items_u7c <- rbind(items_u6c, attt3) %>%
      arrange(version, rand)
  }
  
  # (add task)
  items_u8a <- items_u7a%>%
    group_by(version)%>%
    mutate(task=min_rank(rand))
  items_u8b <- items_u7b%>%
    group_by(version)%>%
    mutate(task=min_rank(rand))
  if(show_eachitem_attraction==3){
    items_u8c <- items_u7c%>%
      group_by(version)%>%
      mutate(task=min_rank(rand))
  }
  
  # (final output selection)
  if(show_eachitem_attraction==2){
    items_u9a <- as.data.frame(items_u8a)%>%
      transmute(version, task, concept)
    items_u9b <- as.data.frame(items_u8b)%>%
      transmute(version, task, concept)
    # write.table(items_u9a, file = paste0(dir,"AttractionDesign_1.csv"), sep = ",", na = ".", row.names = FALSE)
    # write.table(items_u9b, file = paste0(dir,"AttractionDesign_2.csv"), sep = ",", na = ".", row.names = FALSE)
    attraction_1 = items_u9a
  }
  
  if(show_eachitem_attraction==3){
    items_u9 <- rbind(items_u8a, items_u8b, items_u8c)
    items_u9 <- items_u9[,-4]
    items_u9 <- arrange(items_u9,version)
    
    v_tasks <- c(1:(att_n_tasks+3))
    v_tasks <- rep(v_tasks, times=nversion_a)
    items_u9 <- cbind(items_u9,task=v_tasks)
    
    half_task <- round(max(v_tasks)/2)
    
    items_u10a <- items_u9
    items_u10a <- items_u10a[items_u10a$task<(half_task+1),]
    
    items_u10b <- items_u9
    items_u10b <- items_u10b[items_u10b$task>half_task,]
    items_u10b <- items_u10b[,-4]
    v_tasks2 <- c(1:(att_n_tasks+3-half_task))
    v_tasks2 <- rep(v_tasks2, times=nversion_a)
    items_u10b <- cbind(items_u10b,task=v_tasks2)
    items_u10b <- rbind(items_u10b)
  
    items_u10a <- as.data.frame(items_u10a)%>%
      transmute(version, task, concept)
    items_u10b <- as.data.frame(items_u10b)%>%
      transmute(version, task, concept)
    # write.table(items_u10a, file = paste0(dir,"AttractionDesign_1.csv"), sep = ",", na = ".", row.names = FALSE)
    # write.table(items_u10b, file = paste0(dir,"AttractionDesign_2.csv"), sep = ",", na = ".", row.names = FALSE)
    attraction_1 = items_u10a
  }
  return(attraction_1)
}


attraction_function_2 <- function(items_u3, ntest, ntest_comp, ntest_perver, show_eachitem_attraction) { 
  set_test <- 1:ntest
  
  if (ntest_comp == 0) {
    set_comp <- NULL
  } else {
    set_comp <- ((ntest+1):(ntest+ntest_comp))
  }
  
  #3. Change attraction output format to designed formats
  
  sumi3<- items_u3 %>%
    group_by(version)%>%
    summarize(ncon = n_distinct(concept))
  
  nconcept_a<- mean(sumi3$ncon) 
  
  items_u4 <- items_u3 %>%
    transmute(version, concept)%>%
    group_by(version)%>%
    distinct(concept)
  
  items_u4 <- as_tibble(items_u4)
  
  rand1=as.data.frame(runif(1:nrow(items_u4)))
  colnames(rand1)[1] <- "rand"    
  
  rand2=as.data.frame(runif(1:nrow(items_u4)))
  colnames(rand2)[1] <- "rand" 
  
  items_u5a <- cbind(items_u4, rand1)
  items_u5b <- cbind(items_u4, rand2)
  
  items_u6a <- items_u5a %>%
    arrange(version, rand)
  items_u6b <- items_u5b %>%
    arrange(version, rand)
  
  if(show_eachitem_attraction==3){
    rand3=as.data.frame(runif(1:nrow(items_u4)))
    colnames(rand3)[1] <- "rand"
    items_u5c <- cbind(items_u4, rand3)
    items_u6c <- items_u5c %>%
      arrange(version, rand)
  }
  
  # (make a new col for attention test question)
  
  nversion_a <- (max(items_u6a$version))
  nconcept_a <- (ntest_perver+length(set_comp))
  max_con <- max(set_test,set_comp)
  att_n_tasks = show_eachitem_attraction*nconcept_a
  
  attt1 <- tibble(version=(1:nversion_a), concept=max_con+1, rand=runif(1:nversion_a))
  attt2 <- tibble(version=(1:nversion_a), concept=max_con+2, rand=runif(1:nversion_a))
  attt3 <- tibble(version=(1:nversion_a), concept=max_con+1, rand=runif(1:nversion_a))
  
  # (merge)
  items_u7a <- rbind(items_u6a, attt1) %>%
    arrange(version, rand)
  items_u7b <- rbind(items_u6b, attt2) %>%
    arrange(version, rand)
  if(show_eachitem_attraction==3){
    items_u7c <- rbind(items_u6c, attt3) %>%
      arrange(version, rand)
  }
  
  # (add task)
  items_u8a <- items_u7a%>%
    group_by(version)%>%
    mutate(task=min_rank(rand))
  items_u8b <- items_u7b%>%
    group_by(version)%>%
    mutate(task=min_rank(rand))
  if(show_eachitem_attraction==3){
    items_u8c <- items_u7c%>%
      group_by(version)%>%
      mutate(task=min_rank(rand))
  }
  
  # (final output selection)
  if(show_eachitem_attraction==2){
    items_u9a <- as.data.frame(items_u8a)%>%
      transmute(version, task, concept)
    items_u9b <- as.data.frame(items_u8b)%>%
      transmute(version, task, concept)
    # write.table(items_u9a, file = paste0(dir,"AttractionDesign_1.csv"), sep = ",", na = ".", row.names = FALSE)
    # write.table(items_u9b, file = paste0(dir,"AttractionDesign_2.csv"), sep = ",", na = ".", row.names = FALSE)
    attraction_2 = items_u9b
  }
  
  if(show_eachitem_attraction==3){
    items_u9 <- rbind(items_u8a, items_u8b, items_u8c)
    items_u9 <- items_u9[,-4]
    items_u9 <- arrange(items_u9,version)
    
    v_tasks <- c(1:(att_n_tasks+3))
    v_tasks <- rep(v_tasks, times=nversion_a)
    items_u9 <- cbind(items_u9,task=v_tasks)
    
    half_task <- round(max(v_tasks)/2)
    
    items_u10a <- items_u9
    items_u10a <- items_u10a[items_u10a$task<(half_task+1),]
    
    items_u10b <- items_u9
    items_u10b <- items_u10b[items_u10b$task>half_task,]
    items_u10b <- items_u10b[,-4]
    v_tasks2 <- c(1:(att_n_tasks+3-half_task))
    v_tasks2 <- rep(v_tasks2, times=nversion_a)
    items_u10b <- cbind(items_u10b,task=v_tasks2)
    items_u10b <- rbind(items_u10b)
    
    items_u10a <- as.data.frame(items_u10a)%>%
      transmute(version, task, concept)
    items_u10b <- as.data.frame(items_u10b)%>%
      transmute(version, task, concept)
    # write.table(items_u10a, file = paste0(dir,"AttractionDesign_1.csv"), sep = ",", na = ".", row.names = FALSE)
    # write.table(items_u10b, file = paste0(dir,"AttractionDesign_2.csv"), sep = ",", na = ".", row.names = FALSE)
    attraction_2 = items_u10b
  }
  return(attraction_2)
}


