##################################################################################################
### This script contains: 
###
### 1. Functions used to generate several datasets from the proposed model functions 
###     (found in 01_Fn_SimPropModel) using parallel processing 
###
### 2. Functions used to run program mark iteratively on a list of datasets generated 
###     from the proposed model
###
### 3. Functions used to plot and analyze the results from the simulation study 
### 
### NOTE: Code should not be run directly in here. Source() in this script using 
###       the code from 11_RunSimulationStudy.R  or 12_RunCompareChatEstimators.R 
###       as if you would do so with an R Package.
### NOTE: 01_Fn_SimPropModel needs to be sourced in as well otherwise the code here will not work
###
##################################################################################################

##################################################################################
# 1. Simulating Multiple CJS Datasets using the 01_Fn_SimPropoModel.R Functions ##
##################################################################################

#Generate Data
sim_cjs_dat <- function(parameter_list,iterations,model_choice="CJS_dgr",ncores = detectCores() - 1){
  
  #Assign Core count (will not use more than the system has minus one)
  cat("Initializing cluster for data generation...\n")
  numcores <- ncores
  
  #Set up Cluster
  cl <- makeCluster(numcores)
  
  # Upload Functions, Variables, and Load Libraries
  cat("Pushing objects to children...\n")
  
  #Get Script Path
  path2scripts <- getwd() %+% "/Scripts"
  
  #Add Libraries 
  clusterEvalQ(cl, library("tidyverse"))

  #Export Variables to clusters
  export <- list(
    "parameter_list", "iterations", "model_choice","path2scripts"
  )
  
  clusterExport(cl, export, envir = environment())
  
  #Add Custom Library for Simulation
  if(model_choice == "CJS"){
    clusterEvalQ(cl, source(paste0(path2scripts,"/01_Fn_SimPropModel.R")))
  } else if(model_choice == "CJS2"){
    clusterEvalQ(cl, source(paste0(path2scripts,"/01_Fn_SimPropModel.R")))
  } else if(model_choice == "CJS_dgr"){
    clusterEvalQ(cl, source(paste0(path2scripts,"/01_Fn_SimPropModel.R")))
  } else if(model_choice == "CJS3"){
    stop("Model under construction")
  } else {
    stop("Pick CJS, CJS2, or CJS3")
  }
  

  # Set Random Seeds
  clusterSetRNGStream(cl)
  
  # Dummy Function
  f <- function(i,parameter_list,model_choice) {
    if(model_choice == "CJS"){
      return(compile_cjs1(parameter_list))
    } else if(model_choice == "CJS2"){
      return(compile_cjs2(parameter_list,raw=TRUE))
    } else if(model_choice == "CJS_dgr"){
      return(compile_cjs_dgr(parameter_list,raw=TRUE))
    } else if(model_choice == "CJS3"){
      stop("Model under construction")
    } else {
      stop("Pick CJS, CJS2, or CJS3")
    }
  }
  
  #Simulate Data
  cat("Generating data...\n")
  
  cjs_dat_list <- parLapply(cl, 1:iterations, f,
                            parameter_list,
                            model_choice)
  
  cat("Data generation complete...\n")
  
  # Close Cluster
  stopCluster(cl)
  
  #Return Results
  return(cjs_dat_list)
}


############################################################################
# 2. Functions for calling program MARK and running large scale simulations#
############################################################################

#Run program MARK (ONE RUN)
model_cjs_data <- function(data,grouping,delete = TRUE){
  
  #Dummy Function
  run_mark <- function(i,data,grouping,delete){
    
    #Unlist Survival and Recapture Parameters
    phi.f <- data[[i]]$phi.f[1]
    phi.m <- data[[i]]$phi.m[1]
    phi <- mean(phi.f,phi.m)
    p.f <- data[[i]]$p.f[1]
    p.m <- data[[i]]$p.m[1]
    p <- mean(p.f,p.m)
    
    if(grouping == "N"){
      #Choose Appropriate CJS Model Settings
      phi.true <- phi
      p.true <- p
      param.true <- c(phi.true,p.true)
      phi.name <-c("phi")
      p.name <- c("p")
      param.names <- c(phi.name,p.name)
      
      #Process Mark Data (Extract Recapture/sex)
      dat.process <- data[[i]]$x %>% 
        as.data.frame() %>% 
        unite("ch",sep="") 
      
      mark.process <- process.data(dat.process,model="CJS") 
      
      #Design Data List
      mark.ddl <- make.design.data(mark.process) 
      
      #Generate Estimates
      mark_out <- mark(mark.process,
                       mark.ddl,
                       profile.int = FALSE,
                       invisible=TRUE,
                       brief = TRUE,
                       delete = delete,
                       output = FALSE)
      
    } else if(grouping == "B"){
      
      #Choose Appropriate CJS Model Settings
      phi.grp <- list(formula = ~sex)
      p.grp <- list(formula = ~sex)
      phi.true <- c(phi.f,phi.m)
      p.true <- c(p.f,p.m)
      param.true <- c(phi.true,p.true)
      phi.name <- c("phi.f","phi.m")
      p.name <- c("p.f","p.m")
      param.names <- c(phi.name,p.name)
      
      #Process Mark Data (Extract Recapture/sex)
      dat.process <- data[[i]]$x %>% 
        as.data.frame() %>% 
        unite("ch",sep="") %>% 
        mutate(sex = as.factor(data[[i]]$sex)) %>%
        arrange(sex) 
      
      mark.process <- process.data(dat.process,model="CJS",groups="sex") 
      
      #Design Data List
      mark.ddl <- make.design.data(mark.process) 
      
      #Generate Estimates
      mark_out <- mark(mark.process,
                       mark.ddl,
                       model.parameters=list(Phi=phi.grp,p=p.grp),
                       profile.int = FALSE,
                       invisible=TRUE,
                       brief = TRUE,
                       delete = delete,
                       output = FALSE)
      
    } else if(grouping == "S"){
      #Choose Appropriate CJS Model Settings
      phi.grp <- list(formula = ~sex)
      phi.true <- c(phi.f,phi.m)
      p.true <- c(p)
      param.true <- c(phi.true,p.true)
      phi.name <- c("phi.f","phi.m")
      p.name <- c("p")
      param.names <- c(phi.name,p.name)
      
      #Process Mark Data (Extract Recapture/sex)
      dat.process <- data[[i]]$x %>% 
        as.data.frame() %>% 
        unite("ch",sep="") %>% 
        mutate(sex = as.factor(data[[i]]$sex)) %>%
        arrange(sex) 
      
      mark.process <- process.data(dat.process,model="CJS",groups="sex") 
      
      #Design Data List
      mark.ddl <- make.design.data(mark.process) 
      
      #Generate Estimates
      mark_out <- mark(mark.process,
                       mark.ddl,
                       model.parameters=list(Phi=phi.grp),
                       profile.int = FALSE,
                       invisible=TRUE,
                       brief = TRUE,
                       delete = delete,
                       output = FALSE)
      
    } else if(grouping == "R"){
      
      #Choose Appropriate CJS Model Settings
      p.grp <- list(formula = ~sex)
      phi.true <- c(phi)
      p.true <- c(p.f,p.m)
      param.true <- c(phi.true,p.true)
      phi.name <- c("phi")
      p.name <- c("p.f","p.m")
      param.names <- c(phi.name,p.name)
      
      #Process Mark Data (Extract Recapture/sex)
      dat.process <- data[[i]]$x %>% 
        as.data.frame() %>% 
        unite("ch",sep="") %>% 
        mutate(sex = as.factor(data[[i]]$sex)) %>%
        arrange(sex) 
      
      mark.process <- process.data(dat.process,model="CJS",groups="sex") 
      
      #Design Data List
      mark.ddl <- make.design.data(mark.process) 
      
      #Generate Estimates
      mark_out <- mark(mark.process,
                       mark.ddl,
                       model.parameters=list(p=p.grp),
                       profile.int = FALSE,
                       invisible=TRUE,
                       brief = TRUE,
                       delete = delete,
                       output = FALSE)
      
    } else {
      stop("Choose appropriate grouping: \"N\", \"R\", \"S\", or \"B\"")
    }
    
    mark_out <- mark_out$results
    
    #Extract values of use
    gof.stats <- data.frame("iter"=i,
                            "lnl"=mark_out$lnl,
                            "npar" =mark_out$npar,
                            "deviance"=mark_out$deviance,
                            "deviance.df"=mark_out$deviance.df)  %>%
      rename(`-2lnl` = "lnl")
    
    mark.stats <- mark_out$real[,1:4] %>% t() %>% t() %>%
      unname() %>%
      data.frame() %>%
      rename("Est" = X1,"SE" = X2,"LB" = X3, "UB" = X4) %>%
      mutate(Range = UB - LB,Bias = Est - param.true,Parameter=param.names)
    
    
    mark_results <- cbind(mark.stats,gof.stats)
    #Return Output
    return(mark_results)
  }
  
  #List Apply 
  output <- lapply(1:length(data), run_mark,data,grouping,delete) %>% bind_rows()
  
  #Returns output
  return(output)
}



#Grid Generator (for multiple grid simulation)
generate_grid <- function(grid_base,to_vary){
  #Initialize Variables
  parameter_grid <- list()
  param_names <- names(to_vary)
  grid_entries <- bind_rows(to_vary[1])
  
  #Unpack Grid Entries 
  if(length(to_vary) > 1){
    for(i in 2:length(param_names)){
      param <- param_names[i]
      param_entries <- bind_rows(to_vary[param])
      grid_entries <- merge(grid_entries,param_entries,by=NULL)
    }
  }
  
  #Create Parameter Grid
  for(j in 1:nrow(grid_entries)){
    for(l in 1:length(param_names)){
      param <- param_names[l]
      dimension <- length(grid_base[param][[1]])
      grid_base[[param]] <- rep(as.vector(t(grid_entries[j,][param])),dimension)
    }
    parameter_grid[[j]] <- grid_base
  }
  #Return Results
  return(list(parameter_grid,grid_entries))
}


#Compute Monte Carlo Estimates for Simulation Grid
summarize_grid <- function(data,grid_entries){
  
  #For Storage 
  mc_list <- list()
  
  #Monte Carlo Estimates
  for(i in 1:length(data)){
    
    #True Parameters
    true_param <- grid_entries[i,] %>%
      t() %>%
      data.frame() %>%
      rownames_to_column(var = "Parameter") 
    colnames(true_param) <- c("Parameter", "Truth")
    
    #Perform Monte Carlo Calculation 
    mc_list[[i]] <- data[[i]] %>% left_join(true_param,by=c("Parameter")) %>%
      mutate(Coverage = 1*(Truth <= UB)*(Truth >= LB),
             Bias = Bias/Truth,
             Range = Range/Truth,
             c.hat = deviance/deviance.df) %>% 
      group_by(Parameter) %>% 
      summarize(Est = mean(Est),SE = mean(SE),Range=mean(Range),Bias=mean(Bias),
                Coverage=mean(Coverage),Chat=median(c.hat,na.rm = TRUE)) %>%
      ungroup() %>% 
      cbind(as.tibble(grid_entries[i,]))
  }
  #Combine Results and store in batch list 
  return(bind_rows(mc_list))
}

#Add covariance part for grid simulation within mathematical bounds
add_cov_grid <- function(to_vary,parameter_list,param,by=0.5){
  
  #Covariance Bounds
  grid_entries <- generate_grid(parameter_list,to_vary)[[2]]
  #Add Parameters
  if(param == "gamma"){
    #Gamma Grid
    gam_stats <- compute_jbin_cjs(grid_entries$phi.f,grid_entries$phi.m)
    gam_bounds <- cbind(gam_stats[[5]],gam_stats[[6]])
    gam_range <- cbind(max(gam_bounds[,1]),min(gam_bounds[,2]))
    gam <- seq(gam_range[1],gam_range[2],by=by) %>% round(3)
    to_vary[["gamma"]] <- gam
  } else if(param== "rho"){
    #Rho Grid
    rho_stats <- compute_jbin_cjs(grid_entries$p.f,grid_entries$p.m)
    rho_bounds <- cbind(rho_stats[[5]],rho_stats[[6]])
    rho_range <- cbind(max(rho_bounds[,1]),min(rho_bounds[,2]))
    vrho <- seq(rho_range[1],rho_range[2],by=by) %>% round(3)
    to_vary[["rho"]] <- vrho
  }
  #Return Results
  return(to_vary)
}


#Simulate Across Parameter Grid
simulate_grid <- function(grid_base,to_vary,model_choice,iterations,ncores = detectCores() - 1,
                          batches = NA,model_groups=c("B","N"),batch_save=T){
  
  #Create Parameter Grid
  grid_data <- generate_grid(grid_base,to_vary)
  parameter_grid <- grid_data[[1]]
  grid_entries <- grid_data[[2]]
  model_groups <- unique(model_groups)
  results_list <- list()
  mc_results_list <- list()
  
  #Generate batches automatically 
  if(is.na(batches)){
    batch_length <- 10
    batches <- round(length(parameter_grid)/batch_length,0)
    if(batches == 0){batches <- 1}
    cat("Fixing batches to #" %+% batches %+% " batches" %+% "...\n")
  }
  
  #Batch Process for Memory Purposes
  for(j in 1:batches){
    
    cat("Executing batch #" %+% j %+% " of " %+% batches %+% "...\n")
    
    #Length of Each Batch
    batch_length <- round(length(parameter_grid)/batches,0)
    if(batch_length == 0){batch_length <- 1}
    
    #Slice of the parameter_grid
    list_slice <- 1:batch_length + (j-1)*batch_length
    
    #Truncate Grid if it goes over the size of the parameter grid
    if(max(list_slice) >= length(parameter_grid)){list_slice <- list_slice[1]:length(parameter_grid)}
    
    #Report List Slice
    cat("Batch contains entries " %+% list_slice[1] %+% "-" %+% max(list_slice) %+% "...\n")  
    cat("Generating data for entries " %+% list_slice[1] %+% "-" %+% max(list_slice) %+% "...\n")  
    
    #Slice Parameter Grid
    parameter_grid_slice <- parameter_grid[list_slice]  
    
    #Generate Grid of data
    cjs_dat_grid <- lapply(parameter_grid_slice,sim_cjs_dat,
                           model_choice=model_choice,
                           iterations = iterations,
                           ncores=ncores)
    #Model Data
    cat("Modelling batch data with program MARK for entries " %+% list_slice[1] %+% "-" %+% max(list_slice) %+% "...\n")  
  
    #Compute Nested Models
    for(i in 1:length(model_groups)){
      grouping <- model_groups[i]
      if(j == 1){results_list[[i]] <- list()}
      results_list[[i]] <- c(results_list[[i]],lapply(cjs_dat_grid,model_cjs_data,grouping = grouping))
      #mc_results_list[[i]] <- c(mc_results_list[[i]],summarize_grid(results_list[[i]],grid_entries))
    }
    
    #Report Completion 
    cat("Results generated for the nested CJS Model:" %+% model_groups %+% "\n")
    cat("Modelling batch #" %+% j %+% " of " %+% batches %+% " complete...\n")

    
    #Save Batch Results in case of crash
    if(batch_save == T){
      names(results_list) <- model_groups
      saveRDS(list(results=results_list,
                   parameter_grid = parameter_grid,
                   grid_entries = grid_entries,
                   model_groups),
              path2out %+% "Grid_Batch#" %+% j %+% " of " %+% batches %+% "-Results.rds")
      saveRDS(cjs_dat_grid,path2out %+% "Grid_Batch#" %+% j %+% " of " %+% batches %+% "-Data.rds")
    }
    
    #Free Up Memeory
    rm(cjs_dat_grid)
    rm(parameter_grid_slice)
    
  }
  
  #Compute Monte Carlo Estimates
  cat("Monte Carlo estimates compiled and stored for all batches...\n")
  mc_results_list <- lapply(results_list,summarize_grid,grid_entries)
  
  #Assign Model Names
  names(results_list) <- model_groups
  names(mc_results_list) <- model_groups
  
  #Save Results if needed
  if(batch_save == T){
    saveRDS(list(results=results_list,mc_results=mc_results_list,
                 parameter_grid = parameter_grid,
                 grid_entries = grid_entries,
                 model_groups),file=path2out  %+% Sys.Date() %+% "-CJS2_Grid_Full.rds")
  }
  
  cat("Grid-batch modelling complete!...\n")
  #Return Dataframe 
  return(list(results=results_list,mc_results=mc_results_list,
              parameter_grid = parameter_grid,
              grid_entries = grid_entries,
              model_groups))
  
}

###########################################################################
# 4. Functions used to explore results in for one setting in the simulation. 
#
#
# Note the results of these functions do not show up the manuscript as they 
# allow for low-level testing and are difficult to present in a way that shows
# a broad outcome.However, they can be useful for exploration and providing
# insight into the CJS model. As such, I've left them here to be used for 
# further exploration.  
#
# See 13_Run_FurtherInvestigation.R for some examples 
############################################################################

#Duplicate Data observations in Simulated CJS Data (see how change impacts sd)
# Arg: Data refers to data generated with sim_cjs_dat
double_observed <- function(data){
  for(i in 1:length(data)){
    data[[i]]$first <- c(data[[i]]$first,data[[i]]$first)
    data[[i]]$sex <- c(data[[i]]$sex,data[[i]]$sex)
    data[[i]]$x <- rbind(data[[i]]$x, data[[i]]$x)
    data[[i]]$a <- rbind(data[[i]]$a, data[[i]]$a)
    data[[i]]$n <- data[[i]]$n*2 
  }  
  return(data)
}

#Randomize Gender (see how change impacts likelihood)
# Arg: Data refers to data generated with sim_cjs_dat
gender_randomize <- function(data){
  for(i in 1:length(data)){
    data[[i]]$sex <- ifelse(rbinom(data[[i]]$n,1,0.5)==1,"F","M")
  }  
  return(data)
}

#Reduce to one gender (see how change impacts c-hat)
gender_reduce <- function(data, gender_keep){
  for(i in 1:length(data)){
    mask <- data[[i]]$sex == gender_keep
    data[[i]]$n <- sum(mask)
    data[[i]]$x <- data[[i]]$x[mask,] 
    data[[i]]$a <- data[[i]]$a[mask,] 
    data[[i]]$sex <- data[[i]]$sex[mask] 
    data[[i]]$first <- data[[i]]$first[mask] 
  }
  return(data)
}

#Generate MARK CI for Probabilities
# Arg: Prob is value btwn [0, 1]
# Arg: Se is a standard error value
# Arg: Alpha is a confidence level
compute_mark_ci <- function(prob,se,alpha=0.05){
  var.logit <- (se^2)/((prob-1)^2*prob^2)
  se.logit <- sqrt(var.logit)
  est.logit <- log(prob/(1-prob),base=exp(1))
  ub.logit <- est.logit + round(qnorm(1-alpha/2),2)*se.logit
  lb.logit <- est.logit - round(qnorm(1-alpha/2),2)*se.logit  
  lb.real <- exp(lb.logit)/(1+exp(lb.logit))
  ub.real <- exp(ub.logit)/(1+exp(ub.logit))
  return(list(lb.real,ub.real))
} 


#Compute Median C-Hat using chi-sq approach
# Arg: Data refers to results generated with model_cjs_data()
compute_chat <- function(data){
  
  #Deviance C-hat
  chat <- data$deviance/data$deviance.df
  
  #Take Median
  med.chat <- median(chat)
  return(med.chat)
}

#Summarize Single Run with Monte Carlo Method
# Arg: Data refers to results generated with model_cjs_data()
mc_cjs_out <- function(data,parameter_list,grouping,chat=1){
  #List variables
  mc.list <- list()
  phi.name <- ifelse(grouping %in% c("B","S"),list(c("phi.f","phi.m")),c("phi"))[[1]]
  p.name <- ifelse(grouping %in% c("B","R"),list(c("p.f","p.m")),c("p"))[[1]]
  param.names <- c(phi.name,p.name)
  
  #Generate MC Results
  for(i in 1:length(param.names)){
    param.name <- param.names[i]
    param <- parameter_list[[param.name]][1]
    
    mc.list[[i]] <- data %>% 
      filter(Parameter == param.name) %>% 
      mutate(Coverage = 1*(param <= UB)*(param >= LB),
             Bias = (Bias )/param,
             Range = Range/param,
             SE.C = SE*sqrt(chat),
             UB.C = compute_mark_ci(Est,SE.C)[[2]],
             LB.C = compute_mark_ci(Est,SE.C)[[1]],
             Range.C = (UB.C-LB.C)/param,
             Coverage.C = 1*(param <= UB.C)*(param >= LB.C)) %>%
      dplyr::select(Est,SE,Range,Bias,Coverage,SE.C,Range.C,Coverage.C) %>% 
      colMeans() %>% round(5)
  }
  #Combine and Summarize
  mc.results <- bind_rows(!!!mc.list)
  mc.results$Parameter <- param.names
  return(mc.results)
}


#Plots for a Single Simulation (one parameter set)
# Arg: output refers to results generated with model_cjs_data()
# Arg: param refers to a true parameter value for phi or p
# Arg: refers to either Survival or Recapture
# Arg: gender refers to either Male or Female
model_plots <- function(output,param ,type,gender){
  
  #Survival or Recapture Filter
  if(type == "Survival"){
    if(gender == "Male"){
      dat <- filter(output,Parameter == "phi.m")
    } else if(gender == "Female"){
      dat <- filter(output,Parameter == "phi.f")
    } else if(gender == "Flat"){
      dat <- filter(output,Parameter == "phi")
    } else {
      stop("Gender must be Male or Female")
    }
  } else if(type == "Recapture") {
    if(gender == "Male"){
      dat <- filter(output,Parameter == "p.m")
    } else if(gender == "Female"){
      dat <- filter(output,Parameter == "p.f")
    } else if(gender == "Flat"){
      dat <- filter(output,Parameter == "p")
    } else {
      stop("Gender must be Male, Female, or Flat")
    }
  } else {
    stop("Pick either Survival or Recapture!")
  }
  
  if(missing(param)){
    stop("Pick the true parameter value for" %+% type)
  }
  
  #Estimates
  p1 <- ggplot(data = dat) + 
    geom_point(mapping = aes(y = Est, x = iter)) +
    geom_errorbar(aes(ymin=LB, ymax=UB,x=iter)) + 
    geom_hline(yintercept = param,col="blue") + 
    geom_hline(yintercept = colMeans(dat[1]) %>% as.double(),col="red") + 
    labs(title= paste0(gender," ",type, " Estimates"), y="Estimate +/- 95% CI", x = "Iteration")+
    theme_bw() + theme(plot.title = element_text(hjust = 0.5))
  
  #Standard Errors
  p2 <- ggplot(data = dat) + 
    geom_point(mapping = aes(y = SE, x = iter)) +
    geom_hline(yintercept = colMeans(dat[2])%>% as.double(),col="red") + 
    labs(title=paste0(gender," ",type, " Standard Error"), y="Standard Error", x = "Iteration")+
    theme_bw() + theme(plot.title = element_text(hjust = 0.5))
  
  #Interval Length
  p3 <- ggplot(data = dat) + 
    geom_point(mapping = aes(y = Range, x = iter)) +
    geom_hline(yintercept = colMeans(dat[5]) %>% as.double(),col="red") + 
    labs(title=paste0(gender," ",type," Interval Length"), y="Interval Length", x = "Iteration")+
    theme_bw() + theme(plot.title = element_text(hjust = 0.5))
  
  #Bias
  p4 <- ggplot(data = dat) + 
    geom_point(mapping = aes(y = abs(Bias), x = iter)) +
    geom_hline(yintercept = 0,col="blue") + 
    geom_hline(yintercept = colMeans(abs(dat[6])) %>% as.double(),col="red") + 
    labs(title=paste0(gender," ",type," Bias"), y="Bias", x = "Iteration")+
    theme_bw() + theme(plot.title = element_text(hjust = 0.5))
  
  #Return Results
  grid.arrange(p1,p2,p3,p4)
} 

#Likelihood Ratio Test (Red is chi-square distn, blue is observed from lrt)
# Arg: Null and alt are generated with model_cjs_data() and refer to null and alternative hypothesis of the LRT
lrt_plots <- function(null,alt){
  
  null <- null %>% dplyr::select(iter,`-2lnl`,npar) %>% distinct()
  alt <- alt %>% dplyr::select(iter,`-2lnl`,npar) %>% distinct()
  
  test <- null %>% inner_join(alt,by=c("iter")) %>% 
    rename(`-2lnl.null` = `-2lnl.x`,npar.null = npar.x,`-2lnl.alt` = `-2lnl.y`,npar.alt = npar.y) %>% 
    mutate(lrt.stat = `-2lnl.null`-`-2lnl.alt`, 
           df.lrt =  npar.alt-npar.null,
           chisq2.randint = rchisq(max(iter),df.lrt,0),
           p_value = pchisq(lrt.stat,df=df.lrt,lower.tail=FALSE),
           p_value_chisq2 = pchisq(chisq2.randint,df.lrt,lower.tail=FALSE)) 
  
  
  p1 <- test %>% ggplot(mapping = aes(x = iter, y = lrt.stat)) + 
    geom_point() + geom_line() + 
    labs(title="Deviance Statistic for Likelihood Ratio Test", 
         y= expression(-2*"Log("*Delta*")"), x = "Iteration")+
    theme_bw() + theme(plot.title = element_text(hjust = 0.5))
  p2 <- test %>% ggplot() + geom_density(mapping = aes(lrt.stat),fill="blue",alpha=0.2)  + 
    geom_density(mapping = aes(chisq2.randint),fill="red",alpha=0.2)  + 
    labs(title="Density of Deviance Statistic", 
         y="Density", x = expression(-2*"Log("*Delta*")"),
         fill = "")+
    theme_bw() + theme(plot.title = element_text(hjust = 0.5)) #+ 
  # theme(legend.position = "bottom") +
  # scale_fill_manual(values = c('red','blue'),
  #                   labels=c(expression(-2*"Log("*Delta*")"),expression(chi^2)))
  p3 <- test %>% ggplot(mapping = aes(x = iter, y = p_value)) + 
    geom_point() + geom_line() + geom_hline(yintercept = 0.05 ,col="red")  +
    labs(title="P-Values for Likelihood Ratio Test", 
         y=expression("P("*chi^2*">"*-2*"Log("*Delta*"))"), x = "Iteration")+
    theme_bw() + theme(plot.title = element_text(hjust = 0.5))
  p4 <- test %>% ggplot() + geom_density(mapping = aes(p_value),fill="blue",alpha=0.2)  + 
    geom_density(mapping = aes(p_value_chisq2),fill="red",alpha=0.2) +
    labs(title="Density of P-Values", 
         y="Density", x = expression("P("*chi^2*">"*-2*"Log("*Delta*"))"))+
    theme_bw() + theme(plot.title = element_text(hjust = 0.5))
  grid.arrange(p1,p2,p3,p4)
}

#LRT with quasi-correction applied (Red is F distn, blue is observed from q-lrt)
# Arg: Null and alt are generated with model_cjs_data() and refer to null and alternative hypothesis of the LRT
# chat: refers to the value of c-hat to be used for the quasi-likelihood correction
qlrt_plots <- function(null,alt,chat){
  
  null <- null %>% dplyr::select(iter,`-2lnl`,npar) %>% distinct()
  alt <- alt %>% dplyr::select(iter,`-2lnl`,npar,deviance,deviance.df) %>% distinct()
  
  if(missing(chat))  chat <- median(alt$deviance/alt$deviance.df)
  
  lrt.test <- null %>% inner_join(alt,by=c("iter")) %>%
    rename(`-2lnl.null` = `-2lnl.x`,npar.null = npar.x,`-2lnl.alt` = `-2lnl.y`,npar.alt = npar.y) %>%
    mutate(lrt.stat = `-2lnl.null`-`-2lnl.alt`,
           df.lrt =  npar.alt-npar.null,
           chisq2.randint = rchisq(max(iter),df.lrt,0),
           p_value = pchisq(lrt.stat,df=df.lrt,lower.tail=FALSE),
           p_value_chisq2 = 1-pchisq(chisq2.randint,df.lrt))
  
  f.test <- null %>% inner_join(alt,by=c("iter")) %>%
    rename(`-2lnl.null` = `-2lnl.x`,npar.null = npar.x,`-2lnl.alt` = `-2lnl.y`,npar.alt = npar.y) %>%
    mutate(f.stat = ((`-2lnl.null`-`-2lnl.alt`)/(chat*(npar.alt-npar.null))),
           df1 =  npar.alt-npar.null,
           df2 = npar.alt,
           f.randint = rf(max(iter),df1,df2),
           p_value = pf(f.stat,df1,df2,lower.tail=FALSE),
           p_value_f = pf(f.randint,df1,df2,lower.tail = FALSE))
  
  
  p1 <- lrt.test %>% ggplot() + geom_density(mapping = aes(lrt.stat),fill="blue",alpha=0.2)  + 
    geom_density(mapping = aes(chisq2.randint),fill="red",alpha=0.2) +
    labs(title="Density of Raw Deviance Statistic", 
         y="Density", x = expression(-2*"Log("*Delta*")"),
         fill = "")+
    theme_bw() + theme(plot.title = element_text(hjust = 0.5)) 
  
  p2 <- f.test %>% ggplot() + geom_density(mapping = aes(f.stat),fill="blue",alpha=0.2)  +  
    geom_density(mapping = aes(f.randint),fill="red",alpha=0.2) + 
    labs(title="Density of Corrected Deviance Statistic", 
         y="Density", x = expression(over(-2*"Log("*Delta*")/("*df[alt]-df[null]*")",hat(C))),
         fill = "")+
    theme_bw() + theme(plot.title = element_text(hjust = 0.5)) 
  
  p3 <- lrt.test %>% ggplot() + geom_density(mapping = aes(p_value),fill="blue",alpha=0.2)  + 
    geom_density(mapping = aes(p_value_chisq2),fill="red",alpha=0.2)  +
    labs(title="Density of P-Values for Raw Deviance", 
         y="Density", x = expression("P("*chi^2*">"*-2*"Log("*Delta*"))")) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5))
  
  p4 <- f.test %>% ggplot() + geom_density(mapping = aes(p_value),fill="blue",alpha=0.2)  + 
    geom_density(mapping = aes(p_value_f),fill="red",alpha=0.2) + 
    labs(title="Density of P-Values for Corrected Deviance", 
         y="Density", x = expression("P(F>"*over(-2*"Log("*Delta*")/("*df[alt]-df[null]*")",hat(C))*")")) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5)) 
  grid.arrange(p1,p2,p3,p4)
}


#C-hat estimate Plots
# Arg: Null and alt are generated with model_cjs_data() and refer to null and alternative hypothesis of the LRT
chat_plots <- function(null,alt){
  null <- null %>% 
    dplyr::select(iter,deviance ,deviance.df) %>%
    distinct() %>% 
    mutate(c.hat = deviance/deviance.df)
  
  alt <- alt %>% 
    dplyr::select(iter,deviance ,deviance.df) %>%
    distinct() %>% 
    mutate(c.hat = deviance/deviance.df)
  
  p1 <- null %>% ggplot(mapping = aes(x = iter, y = c.hat)) +
    geom_point() + 
    geom_line() + 
    geom_hline(yintercept = 1 ,col="blue") + 
    geom_hline(yintercept = 3 ,col="red") + 
    labs(title="Estimated Variance Inflation Factor for Null Model", 
         y=expression(hat(C)), x = "Iteration", fill = "") +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5)) 
  
  p2 <- null %>% ggplot() + 
    geom_density(mapping = aes(c.hat),fill="blue",alpha=0.2) + 
    geom_vline(xintercept = median(null$c.hat) ,col="red",linetype="dashed") + 
    labs(title="Density of Variance Inflation Factor for Null Model", 
         y="Density", x = expression(hat(C)), fill = "") +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5)) 
  
  p3 <- alt %>% ggplot(mapping = aes(x = iter, y = c.hat)) + 
    geom_point() + geom_line() + 
    geom_hline(yintercept = 1 ,col="blue")  +
    geom_hline(yintercept = 3 ,col="red")  + 
    labs(title="Estimated Variance Inflation Factor for Alternative Model", 
         y=expression(hat(C)), x = "Iteration", fill = "") +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5)) 
  
  
  p4 <- alt %>% ggplot() + geom_density(mapping = aes(c.hat),fill="blue",alpha=0.2) + 
    geom_vline(xintercept = median(alt$c.hat) ,col="red",linetype="dashed") + 
    labs(title="Density of Variance Inflation Factor for Alternative Model", 
         y="Density", x = expression(hat(C)), fill = "") +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5)) 
  
  
  grid.arrange(p1,p2,p3,p4)
}
