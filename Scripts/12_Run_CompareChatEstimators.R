###################################################################################################
### This script is used to run the simulation study found in Appendix B.3 Part Three. Namely, the 
### study which compares the deviance, Pearson, and Fletcher C-hat estimators. 
###
### To our knowledge, as of writing this, the R library Rmark does not extract c-hat values from the 
### mark###.out files produced by program Mark. This means we need to do so using this script
###
### This script will generate the simulation data and run the models L122-131 (WARNING: SLOW)
### using the settings in L87-L120. L15-L85 calls in the necessary libs, and custom functions
###
### L133-L149 will extract the c-hat results from the mark###.out files and return a tidy dataframe
### L151 - L163 will produce the plots found in Figure 7 (subject to expected statistical variation) 
####################################################################################################

###################################
# Set up libraries and directories#
###################################
rm(list = ls()) # Clear memory
options(scipen = 999)
options(digits = 3)
## Libraries
libs <- c("tidyverse", "parallel","RMark","gridExtra", "ggpubr")

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

######################################
# Extract Functions only working here#
######################################

# Gather C-Hat results from mark output
extract_manual_chat2 <- function(path2dat){
  data  <- read_table(path2dat) %>% unlist() %>% unname() %>% data.frame(stringsAsFactors = FALSE)
  fletcher <- str_split(data$.[grep("Fletcher",data$.)],"= ")[[1]][2] %>% as.numeric()
  pearson <- str_split(data$.[grep("Pearson chat",data$.)],"= ")[[1]][2] %>% as.numeric()
  deviance <- str_split(data$.[grep("c-hat ",data$.)],"= ")[[1]][2] %>% as.numeric()
  model <-   str_split(data$.[grep("c-hat ",data$.)],"= ")[[1]][1]
  return(data.frame(model,deviance,pearson,fletcher,stringsAsFactors = FALSE))
}

# Compare Estimators  
chat_compare <- function(data, model = "N",
                         main = expression("Model ("*phi*","*p*")"),
                         label = c("Deviance","Pearson","Fletcher")){
  
  # Convert data into long format
  data <- pivot_longer(data,
                       cols = c("deviance","pearson","fletcher"),
                       names_to = "estimator")
  
  #filter rename
  model_id <- model
  
  # Plot Results
  p1 <- data %>% 
    dplyr::filter(model == model_id) %>% 
    ggplot(aes(color=factor(estimator))) + 
    geom_density(mapping = aes(value)) +
    geom_vline(xintercept = 1,linetype="dashed") +
    labs(title = main,y="Density",x=expression(hat(c))) +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme(legend.position="bottom") +
    scale_color_manual(name = expression("Estimator of "*hat(c)),
                       values = c("#F8766D","#00BA38","#619CFF"),labels = label)
  return(p1)
}

#####################################
# Generate Mated CJS Simulation Data#
#####################################
#Parameters
n <- 200 ##Sample Size 
k <- 4 ##Sampling Occasions
delta <- rep(1,k) #Breeding Probability
phi.f <- rep(0.7,k) #Marginal Female Survival from j to j+1
phi.m <- rep(0.7,k) #Marginal Male Survival from j to j+1
phi <- rep(0.7,k) #Marginal survival for any animal
gamma <- rep(1,k) #Correlation between males and female survival
p.f <- rep(0.8,k) #Marginal Female Recapture at j
p.m <- rep(0.8,k) #Marginal Male Recapture at j
p <- rep(0.8,k) #Marginal recapture for any animal
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

######################################
#Simulation for one set of parameters#
######################################
# #Generate Data
cjs_dat_list <- sim_cjs_dat(parameter_list,iterations = 10,model_choice = "CJS_dgr")
#Generate Results
fit.s1 <-  model_cjs_data(data = cjs_dat_list, grouping = "S",delete=FALSE)
fit.1r <- model_cjs_data(data = cjs_dat_list, grouping = "R",delete=FALSE)
fit.sr <- model_cjs_data(data = cjs_dat_list, grouping = "B",delete=FALSE)
fit.11 <- model_cjs_data(data = cjs_dat_list, grouping = "N",delete=FALSE)

#########################################
# Extract Results from mark###.out files#
#########################################

chat_dat <- data.frame(model = "",deviance=0,pearson=0,fletcher=0,stringsAsFactors = FALSE)

for(i in 1:4000){
  if(i < 1000){index <- substr("00" %+% i,nchar("00" %+% i)-2,nchar("00" %+% i))} else {index <- i}
  path2dat <- getwd() %+% "/mark" %+% index %+% ".out"
  chat_dat[i,1] <- suppressMessages(extract_manual_chat2(path2dat)[1,1])
  chat_dat[i,2:4] <- suppressMessages(extract_manual_chat2(path2dat)[1,2:4])
}

# Organize Results for plotting
naming <- data.frame(model = chat_dat$model %>% unique(),short = c("S","R","B","N"),stringsAsFactors = FALSE)
chat_dat <- inner_join(chat_dat,naming) %>% select(-model) %>% rename(`model` = short) %>% 
  select(model,deviance,pearson,fletcher)

#########################################
# Plot Summary of Results for all models#
#########################################

p1 <- chat_compare(chat_dat,"N",main = expression("Model ("*phi*","*p*")")) + ylim(0, 1)
p2 <- chat_compare(chat_dat,"R",main = expression("Model ("*phi*","*p^G*")")) + ylim(0, 1)
p3 <- chat_compare(chat_dat,"S",main = expression("Model ("*phi^G*","*p*")")) + ylim(0, 1)
p4 <- chat_compare(chat_dat,"B",main = expression("Model ("*phi^G*","*p^G*")")) + ylim(0, 1)

# Put into one 4-panel plot 
ggarrange(p1, p2, p3, p4, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
# Save Data
saveRDS(chat_dat,file=path2out %+% "chat_dat_resim.rds")
