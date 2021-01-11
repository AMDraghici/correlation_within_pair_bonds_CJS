###################################
# Set up libraries and directories#
###################################
rm(list = ls()) # Clear memory
options(scipen = 999)
options(digits = 3)
## Libraries
libs <- c("tidyverse", "parallel","RMark","gridExtra")

## If package is not installed then do so
for (pkg in libs) {
  if (!pkg %in% rownames(installed.packages())) {
    install.packages(pkg)
  }
}

## Attach our libraries
lapply(libs, require, character.only = TRUE)
## Inline object concatenation (like python)
`%+%` <- function(a, b) paste0(a, b)

# ## Set Working Directory and reference paths
path2dat <- getwd() %+% "/Data/"
path2out <- getwd() %+% "/Output/"
path2scripts <- getwd() %+% "/Scripts/"

#####################################
# Call in code from relevant scripts#
#####################################
source(path2scripts %+% "01_Fn_SimPropModel.R")
source(path2scripts %+% "02_Fn_ModelCode.R")

#####################################
# Generate Mated CJS Simulation Data#
#####################################

#Parameters
n <- 200 ##Sample Size 
k <- 4 ##Sampling Occasions
delta <- rep(1,k) #Breeding Probability
phi.f <- rep(0.5,k) #Marginal Female Survival from j to j+1
phi.m <- rep(0.5,k) #Marginal Male Survival from j to j+1
phi <- rep(0.5,k) #Marginal survival for any animal
gamma <- rep(1,k) #Correlation between males and female survival
p.f <- rep(0.5,k) #Marginal Female Recapture at j
p.m <- rep(0.5,k) #Marginal Male Recapture at j
p <- rep(0.5,k) #Marginal recapture for any animal
rho <- rep(1,k) #Correlation between Female and Male Recapture
prob.female <- 0.5 #Proportion of females in the population

#Parameters
parameter_list <- list(n = n,
                       k = k,
                       delta = delta,
                       phi.f = phi.f,
                       phi.m = phi.m,
                       gamma = gamma,
                       p.f = p.f,
                       p.m = p.m,
                       rho = rho,
                       prob.female = prob.female,
                       phi = phi,
                       p = p)

#####################
#Simulate with Rmark#
#####################

###################
#Single Simulation#
###################
# #Generate Data
cjs_dat_list <- sim_cjs_dat(parameter_list,iterations = 100,model_choice = "CJS_dgr")
cjs_dat_list2 <- double_observed(cjs_dat_list)
cjs_dat_list3 <- sim_cjs_dat(parameter_list,iterations = 100,model_choice = "CJS_dgr") #set correlations high
cjs_dat_list4 <- gender_randomize(cjs_dat_list3)
cjs_dat_list5 <- double_observed(cjs_dat_list3)

#Generate Results
fit.sr <- model_cjs_data(data = cjs_dat_list, grouping = "B",delete = TRUE)
fit.s1 <-  model_cjs_data(data = cjs_dat_list, grouping = "S",delete = TRUE)
fit.1r <- model_cjs_data(data = cjs_dat_list, grouping = "R",delete = TRUE)
fit.11 <- model_cjs_data(data = cjs_dat_list, grouping = "N",delete = TRUE)
fit.all <- list(fit.sr,fit.s1,fit.1r,fit.11,parameter_list,cjs_dat_list)

#LRT Gof Test Plot
lrt_plots(null = fit.s1,alt = fit.sr)
lrt_plots(null = fit.1r,alt = fit.sr)
lrt_plots(null = fit.11, alt = fit.s1)
lrt_plots(null = fit.11, alt = fit.1r)
lrt_plots(null = fit.11, alt = fit.sr)

#C-Hat Distribution
chat_plots(null = fit.s1,alt = fit.sr)
chat_plots(null = fit.1r,alt = fit.sr)
chat_plots(null = fit.11, alt = fit.s1)
chat_plots(null = fit.11, alt = fit.1r)
chat_plots(null = fit.11, alt = fit.sr)

#F-Corrected (C-Hat) Gof Test Plot
qlrt_plots(null = fit.1r,alt = fit.sr,compute_chat(fit.sr))
qlrt_plots(null = fit.11, alt = fit.s1,compute_chat(fit.s1))

#MC Result
mc_cjs_out(fit.s1,parameter_list,"S",compute_chat(fit.s1))
mc_cjs_out(fit.11,parameter_list,"N",1)

#View Survival Results
model_plots(fit.s1,param=phi.f,type="Survival",gender="Female")
model_plots(fit.s1,param=phi.m,type="Survival",gender="Male")
model_plots(fit.11,param=phi,type="Survival",gender="Flat")

#View Recapture Results
model_plots(fit.1r,param=p.m,type="Recapture",gender="Male")
model_plots(fit.1r,param=p.f,type="Recapture",gender="Female")
model_plots(fit.11,param=p,type="Recapture",gender="Flat")

#Save Results
saveRDS(fit.all,file=path2out %+% Sys.Date() %+% "fit_all.rds")