#############################################################################################################
### This script contains functions used to simulate a single dataset from the proposed model of the manuscript
### 
### NOTE: Code should not be run directly in here. Source() in this script using 
###       the code from 11_RunSimulationStudy.R  or 12_RunCompareChatEstimators.R 
###       as if you would do so with an R Package.
#############################################################################################################


#Extract joint binomial parameters
compute_jbin_cjs <- function(prob.f,prob.m){
  #Compute Components
  prob.prod <- prob.m * prob.f
  sig.prob.f <- sqrt(prob.f*(1-prob.f))
  sig.prob.m <- sqrt(prob.m*(1-prob.m))
  sig.prod <- sig.prob.f*sig.prob.m
  ##Upper and Lower bounds for cov of joint binomial
  lub <- pmin(1,(prob.f-prob.prod)/(sig.prod),(prob.m-prob.prod)/(sig.prod))
  glb <- pmax(-1,(-prob.prod)/(sig.prod),(prob.f+prob.m - prob.prod - 1)/(sig.prod))
  #Return values
  return(list(prob.prod,sig.prob.f,sig.prob.m,sig.prod,glb,lub))
}


#Generate Derived Survival Probabilities
param_surv_cjs_dgr <- function(prob.f,prob.m,cor){
  #Extract Parameters
  parameters <- compute_jbin_cjs(prob.f,prob.m)
  cor_upper_bound <- parameters[[6]]
  cor_lower_bound <- parameters[[5]]
  sig.prob.f <- parameters[[2]]
  sig.prob.m <- parameters[[3]]
  
  ###Joint Probability Distn for survivorship 
  prob.mf <- cor * sig.prob.m * sig.prob.f + (prob.f*prob.m)
  prob.f0 <- prob.f - prob.mf
  prob.m0 <- prob.m - prob.mf
  prob.00 <- 1 - prob.f0 - prob.m0 - prob.mf 
  #List of parameters
  param <- list(prob.mf,prob.f0,prob.m0,prob.00)
  #Return values
  return(param)
}


#Survival Function
survival_cjs_dgr <- function(previous_state,mated,j,phi.m,phi.f,phi.00,phi.f0,phi.m0,phi.mf){
  if(previous_state == 3){
    next_state <- 
      which(rmultinom(c(1,2,3,4),1,
                      c((1-phi.m[j]),0,phi.m[j],0)) == 1)
  } else if(previous_state == 2){
    next_state <- 
      which(rmultinom(c(1,2,3,4),1,
                      c((1-phi.f[j]),phi.f[j],0,0)) == 1)
  } else if(previous_state == 4) {
    if(mated == 1){
      next_state <-
        which(rmultinom(c(1,2,3,4),1,
                        c((1-phi.m[j])*(1-phi.f[j]),(1-phi.m[j])*phi.f[j],phi.m[j]*(1-phi.f[j]),phi.m[j]*phi.f[j])) == 1)
    } else if(mated ==2){
      next_state <-
        which(rmultinom(c(1,2,3,4),1,
                        c(phi.00[j],phi.f0[j],phi.m0[j],phi.mf[j])) == 1)
    }
  } else {
    next_state <- 1
  }
  return(next_state)
}

#Recapture Function
recapture_cjs_dgr <- function(current_state,mated,j,p.m,p.f,p.00,p.f0,p.m0,p.mf){
  if(current_state == 3){
    obs <- 
      which(rmultinom(c(1,2,3,4),1,
                      c((1-p.m[j]),0,p.m[j],0)) == 1)  
  } else if(current_state == 2){
    obs <- 
      which(rmultinom(c(1,2,3,4),1,
                      c((1-p.f[j]),p.f[j],0,0)) == 1)
  } else if(current_state == 4) {
    if(mated==1){
      obs <-
        which(rmultinom(c(1,2,3,4),1,
                        c((1-p.m[j])*(1-p.f[j]),(1-p.m[j])*p.f[j],p.m[j]*(1-p.f[j]),p.m[j]*p.f[j])) == 1)
    } else if(mated==2){
      obs <-
        which(rmultinom(c(1,2,3,4),1,
                        c(p.00[j],p.f0[j],p.m0[j],p.mf[j])) == 1)
    }
  } else {
    obs <- 1
  }
  return(obs)
}


#Divorce function
divorce_dgr <- function(previous_state,j,delta){
  if(previous_state == 4){
    mated <- rbinom(1,1,delta[j]) + 1
  } else {
    mated <- 1
  }
  return(mated)
}



init_pop_cjs_dgr <- function(n,prob.female,prob.partner=1,k){
  
  #Construct Data Frame with animals, genders, and initial entry
  population <- data_frame(animal_id = 1:n,
                           sex = ifelse(rbinom(n,1,prob.female)==1,"F","M"), #c(rep("F",n/2),rep("M",n/2)),
                           partner = ifelse(rbinom(n,1,prob.partner)==1,T,F)) 
  
  ##Split by mated animals and single animals
  mated_animals <- population %>% filter(partner==T)
  single_animals <- population %>% filter(partner==F)
  
  ###Assign single individual to female column or male column (indexed by postion in 1:n)
  ##0 implies that there is no mate (eg. Male - 177 and Female - 0 implies a single male)
  single_animals <- data.frame(female = ifelse(single_animals$sex == "F",
                                               single_animals$animal_id,0),
                               male = ifelse(single_animals$sex == "M",
                                             single_animals$animal_id,0)
  )
  
  ###Divide mated tables by gender 
  mated_females <- mated_animals %>% filter(sex == "F")
  mated_males <- mated_animals %>% filter(sex == "M") 
  
  #Assign initial pairings
  pairs <- mated_males %>% 
    mutate(partner_id = sample(c(mated_females$animal_id),replace=F)[1:nrow(mated_males)]) %>% 
    rename("male" = animal_id,"female"=partner_id) %>% dplyr::select("female","male") %>% 
    mutate(female = ifelse(is.na(female),0,female))
  
  #Set up data frame containing all relevant information
  entities <- bind_rows(pairs,single_animals,
                        right_join(filter(pairs,female>0),mated_females,by=c("female"="animal_id")) %>%
                          filter(is.na(male)) %>% 
                          select(female,male) %>% mutate(male = replace_na(male,0))) %>% 
    mutate(ID=1:n(),initial_entry=1)#as.integer(sample(c(1:round(k/2)),size=n(),replace=T)))
  
  entities2 <- merge(entities$ID,1:k) %>% 
    rename(ID=x,Time=y) %>% 
    inner_join(entities,by="ID") %>% 
    arrange(ID,Time) %>% 
    mutate(mated = ifelse(initial_entry==Time & pmin(female,male) != 0,2,ifelse(pmin(female,male)==0,1,NA)),
           SV = ifelse(initial_entry==Time,ifelse(female==0,3,ifelse(male==0,2,4)),NA),
           RC = ifelse(initial_entry==Time,ifelse(female==0,3,ifelse(male==0,2,4)),NA),
           a = ifelse(initial_entry==Time,ifelse(female==0,3,ifelse(male==0,2,4)),NA),
           d=mated)
  
  return(entities2)
}

sim_mated_cjs_dgr <- function(n,k,prob.female,prob.partner=1,phi.f,phi.m,gamma,p.f,p.m,rho,delta){
  
  #Generate Data Template
  entities <- init_pop_cjs_dgr(n,prob.female,prob.partner,k)
  
  #Probabilities
  s.param <- param_surv_cjs_dgr(phi.f,phi.m,gamma)
  r.param <- param_surv_cjs_dgr(p.f,p.m,rho)
  
  #Unpack survival probs
  phi.mf <- s.param[[1]]
  phi.f0 <- s.param[[2]]
  phi.m0 <- s.param[[3]]
  phi.00 <- s.param[[4]]
  
  #Unpack recapture Probs
  p.mf <- r.param[[1]]
  p.f0 <- r.param[[2]]
  p.m0 <- r.param[[3]]
  p.00 <- r.param[[4]]
  
  #Assign values sequentially
  for(i in unique(entities$ID)){
    init <- filter(entities,ID==i) %>% dplyr::select(initial_entry) %>% unique() %>% t() %>% as.vector()
    #Loop through unknown states/observations
    for(j in (init):k){
      if(j == init){
        next
      } else {
        #animal_id
        female <- entities[(i - 1) * k + j,3]
        male <-  entities[(i - 1) * k + j,4]
        #Previous survival state
        previous_state <- entities[(i - 1) * k + j - 1,7]
        #Simulated Mated 
        entities[(i - 1) * k + j,6] <- mated <- divorce_dgr(previous_state,j,delta)
        #Survival
        entities[(i - 1) * k + j,7] <- current_state <- survival_cjs_dgr(previous_state,mated,j,phi.m,phi.f,phi.00,phi.f0,phi.m0,phi.mf)
        #Recapture
        entities[(i - 1) * k + j,8] <- current_observation  <- recapture_cjs_dgr(current_state,mated,j,p.m,p.f,p.00,p.f0,p.m0,p.mf)
        #Confounded Survival 
        entities[(i - 1) * k + j,9] <- ifelse(current_observation==1,NA,current_observation)
        #Add confounded mated 
        entities[(i - 1) * k + j,10] <- ifelse(current_observation==4|female==0|male==0,mated,NA)
        
        rm(female,male,mated,current_state,current_observation)
      } 
    } 
  }
  
  #Confounded survival states at time i
  state_confounded <- entities %>% filter(ID == i) %>% dplyr::select(a) %>% t() %>% as.vector()
  
  #Positions 
  last.female.seen <- ifelse(abs(max(which(state_confounded %in%
                                             c(2))))==Inf,0,max(which(state_confounded %in% c(2))))
  last.male.seen <- ifelse(abs(max(which(state_confounded %in% 
                                           c(3))))==Inf,0,max(which(state_confounded %in% c(3))))
  last.both.seen <- ifelse(abs(max(which(state_confounded %in%
                                           c(4))))==Inf,0,max(which(state_confounded %in% c(4))))
  last.seen <- max(last.both.seen,last.male.seen,last.female.seen)
  
  #Update Confounded survival Information
  for(l in init:last.seen){
    if(l == init){
      next
    } else {
      state_confounded[l] <- ifelse(last.both.seen>=l|(last.male.seen>=l&last.female.seen>=l),4,
                                    ifelse(last.male.seen>=l,3,ifelse(last.female.seen>=l,2,NA)))
      entities[(i - 1) * k + l,9] <- state_confounded[l]
    }
  }
  return(entities)
}

#Format MRC data into std cjs model 
extract_cjs_dgr <- function(Data){
  
  #Split data apart for generic CJS model 
  females.seperate <- Data %>% 
    select(Time,female,initial_entry,RC,a) %>% 
    filter(female > 0) %>% 
    mutate(RC = ifelse(RC ==1,0,ifelse(RC==2|RC==4,1,0)),
           a = ifelse(a==4|a==2,1,ifelse(a==1|a==3,0,NA)),
           gender = "F") %>% 
    rename(ID = female)
  
  males.seperate <- Data %>% 
    select(Time,male,initial_entry,RC,a) %>% 
    filter(male > 0) %>% 
    mutate(RC = ifelse(RC ==1,0,ifelse(RC==3|RC==4,1,0)),
           a = ifelse(a==4|a==3,1,ifelse(a==1|a==2,0,NA)),
           gender = "M") %>% 
    rename(ID = male)
  
  #Table of individual results
  split.results <- bind_rows(females.seperate,males.seperate)
  
  #Set up results in list format
  a <- split.results %>% select(ID,Time,a) %>% tidyr::spread(Time,a) %>% select(-ID)
  x <- split.results %>% select(ID,Time,RC) %>% tidyr::spread(Time,RC) %>% select(-ID)
  first <- split.results %>% 
    select(ID,initial_entry) %>% distinct() %>% 
    arrange(ID) %>% select(-ID) %>% t() %>% as.vector()
  sex <- split.results %>% 
    select(ID,gender) %>% distinct() %>% 
    arrange(ID) %>% select(-ID) %>% t() %>% as.vector()
  k <- ncol(a)
  n <- nrow(a) 
  
  # Prior generating function
  PGF <- function(n) {
    phi <- rbeta(n, 40, 10)
    p <- rbeta(n, 40, 10) 
    return(list("phi" = phi, "p" = p))
  }
  
  DGF <- function(n, phi, p, log = TRUE){
    dbeta(phi, 40, 10, log = log) +
      dbeta(p, 40, 10, log = log)
  }
  
  # Param Names
  param.names <- c("Phi", "P")
  
  #Results
  cjs_dat <- list(
    "PGF" = PGF,
    "DGF" = DGF,
    "k" = k,
    "n" = n,
    "a" = a,
    "x" = x,
    "first" = first,
    "param" = param.names,
    "sex" = sex
  )
  
  return(cjs_dat)
}

compile_cjs_dgr <- function(parameter_list,raw= TRUE){
  #Parameters
  n <- parameter_list[["n"]] ##Sample Size 
  k <- parameter_list[["k"]] ##Sampling Occasions
  delta <- parameter_list[["delta"]] #Breeding Probability
  phi.f <- parameter_list[["phi.f"]] #Marginal Female Survival from j to j+1
  phi.m <- parameter_list[["phi.m"]] #Marginal Male Survival from j to j+1
  gamma <- parameter_list[["gamma"]] #Correlation between males and female survival
  p.f <- parameter_list[["p.f"]] #Marginal Female Recapture at j
  p.m <- parameter_list[["p.m"]] #Marginal Male Recapture at j
  rho <- parameter_list[["rho"]] #Correlation between Female and Male Recapture
  prob.female <- parameter_list[["prob.female"]] #Proportion of females in the population
  
  #Simulate Results
  sim <- sim_mated_cjs_dgr(n,k,prob.female,prob.partner=1,phi.f,phi.m,gamma,p.f,p.m,rho,delta)

  #Filter Down to Raw CJS data
  cjs_dat <- extract_cjs_dgr(sim)
  
  #Add Survival/Recapture Probabilities
  cjs_dat$phi.m <- phi.m
  cjs_dat$p.m <- p.m
  cjs_dat$phi.f <- phi.f
  cjs_dat$p.f <- p.f
  
  #Which Data do you want?
  if(raw == TRUE){
    return(cjs_dat)
  } else {
    return(sim)
  }
}
