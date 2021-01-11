############################################################################################################
### This script can be used to produce the plots and tables in the main body and appendix of the manuscript#
### 
### Before running this, you will need the output from the three simulation studies done in the paper. 
###
### 1. The main one discussed in the methods section
### 2. The sample size LRT one in Appendix B.2 Part Three
### 3. The C-hat estimator one in Appendix B.3 Part Three
###
###
### In this script I have assumed that they are called:
### 
###  1. "2020-04-26-CJS_Grid_Full_Resim.rds"
###  2. "2020-05-14-CJS_LRT_N_Resim.rds"
###  3. "chat_dat_resim.rds"
###
### and that they are stored in a folder called output in your working directory
###
### L24-L281 will load libraries, custom functions, and functions (in this document) for plotting the results
### 
### Remaining lines will generate the Figure or Table in their corresponding headers
############################################################################################################

##################################
# Set up libraries and directories#
##################################

## Libraries
libs <- c("tidyverse", "parallel","RMark","gridExtra","ggpubr","knitr","kableExtra")

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
path2scripts <- getwd() %+% "/Scripts/"
path2out <- getwd() %+% "/Output/"

#####################################
# Call in code from relevant scripts#
#####################################
source(path2scripts %+% "01_Fn_SimPropModel.R")
source(path2scripts %+% "02_Fn_ModelCode.R")

#Call in Results
data <- readRDS(path2out %+% "2020-04-26-CJS_Grid_Full_Resim.rds")
lrt_dat <- readRDS(path2out %+% "2020-05-14-CJS_LRT_N_Resim.rds")
chat_dat <- readRDS(path2out %+% "chat_dat_resim.rds")

##########################################
# Functions to Produce Plots in Manuscript
##########################################

sim_surv_plots <- function(data,model,param,x){
  
  #Filter Down Data to Parameter of Interest and Extract True Value
  plot_data <- data$mc_results[[model]] %>% 
    filter(Parameter == param,rho == x) %>% 
    mutate(lb = Est - 1.96*(SE/sqrt(1000)),
           ub = Est + 1.96*(SE/sqrt(1000)))
  
  Truth <- plot_data[,c(param)][1]
  
  Lower.Est <- min(Truth - 0.01,plot_data$Est)
  Upper.Est <- max(Truth + 0.01,plot_data$Est)
  
  p1 <- plot_data %>% 
    ggplot(aes(y = Est, x = gamma)) + 
    geom_point() +
    geom_line() + 
    geom_hline(yintercept = Truth,col="red") +
    labs(title="Estimates of Survival Probability", 
         x=expression(gamma), y = expression(phi))+
    theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + 
    ylim(Lower.Est,Upper.Est) +
    geom_errorbar(aes(ymin=lb, ymax=ub), width=.05) 
  
  Lower.Range <- max(plot_data$Range - 0.005,0)
  Upper.Range <- min(plot_data$Range + 0.005,1)
  
  p2 <- plot_data %>%  
    ggplot(aes(y = Range, x = gamma)) + 
    geom_point() + 
    geom_line() +
    labs(title="95% CI Interval Range on Survival", 
         x=expression(gamma), y = "Relative Interval Range (R)")+
    theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + 
    ylim(Lower.Range,Upper.Range)
  
  p3 <- plot_data %>%  
    ggplot(aes(y = Coverage, x = gamma)) + 
    geom_point() +
    geom_line() +
    geom_hline(yintercept = 0.95,col="red") +
    labs(title="Coverage Percentage of 95% CI on Survival", 
         x=expression(gamma), y = "Coverage Percentage (C)")+
    theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + 
    ylim(0.8,1)
  
  Lower.Bias <- max(plot_data$Bias - 0.01)
  Upper.Bias <- min(plot_data$Bias + 0.01)
  
  p4 <- plot_data %>%  
    ggplot(aes(y = Bias, x = gamma)) + 
    geom_point() +
    geom_line() +   
    geom_hline(yintercept = 0,col="red") +
    labs(title="Relative Bias for Survival", 
         x=expression(gamma), y = "Relative Bias (B)")+
    theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + 
    ylim(Lower.Bias,Upper.Bias) 
  
  grid.arrange(p1,p2,p3,p4)
}

sim_surv_ci_plots <- function(data,x){
  
  #Filter Down Data to Parameter of Interest and Extract True Value
  plot_data_N <- data$mc_results$N %>% filter(Parameter == "phi",rho == x)
  plot_data_R <- data$mc_results$R %>% filter(Parameter == "phi",rho == x)
  plot_data_S <- data$mc_results$S %>% filter(Parameter %in% c("phi.f","phi.m"),rho == x)
  plot_data_B <- data$mc_results$B %>% filter(Parameter %in% c("phi.f","phi.m"),rho == x)
  
  p1 <- plot_data_N %>%  
    ggplot(aes(y = Coverage, x = gamma)) + 
    geom_point() +
    geom_line() +
    geom_hline(yintercept = 0.95,col="red") +
    labs(title=expression("Model ("*phi*","*p*")"), 
         x=expression(gamma), y = "Coverage Percentage (C)")+
    theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + 
    ylim(0.8,1)
  
  p2 <- plot_data_R %>%  
    ggplot(aes(y = Coverage, x = gamma)) + 
    geom_point() +
    geom_line() +
    geom_hline(yintercept = 0.95,col="red") +
    labs(title=expression("Model ("*phi*","*p^G*")"), 
         x=expression(gamma), y = "Coverage Percentage (C)")+
    theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + 
    ylim(0.8,1)
  
  p3 <- plot_data_S %>%  
    ggplot(aes(y = Coverage, x = gamma,color= Parameter)) + 
    geom_point() +
    geom_line() +
    geom_hline(yintercept = 0.95,col="red") +
    labs(title=expression("Model ("*phi^G*","*p*")"), 
         x=expression(gamma), y = "Coverage Percentage (C)")+
    theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + 
    ylim(0.8,1) + #+ guides(color = FALSE) 
    scale_color_discrete(name = "",labels = c(expression(phi^F),expression(phi^M)))
  
  p4 <- plot_data_B %>%  
    ggplot(aes(y = Coverage, x = gamma,color= Parameter)) + 
    geom_point() +
    geom_line() +
    geom_hline(yintercept = 0.95,col="red") +
    labs(title=expression("Model ("*phi^G*","*p^G*")"),  
         x=expression(gamma), y = "Coverage Percentage (C)")+
    theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + 
    ylim(0.8,1) + scale_color_discrete(name = "",
                                       labels = c(expression(phi^F),expression(phi^M)))
  
  # grid.arrange(p1,p2,p3,p4)
  return(ggarrange(p1,p2,p3,p4,common.legend = TRUE,legend = "bottom"))
}


grid_lrt_plots <- function(data,null="N",alt="S",slice,case = "gamma",name = expression("Survival Correlation"*(gamma)),chi_include = FALSE,xmax1 =1, ymax1 = 1,ymax2 = 1){
  
  null.list <- data$results[[null]]
  alt.list <- data$results[[alt]]
  test.list <- list()
  
  for(i in slice){
    
    null.lrt <- null.list[[i]] %>% dplyr::select(iter,`-2lnl`,npar) %>% distinct()
    alt.lrt <- alt.list[[i]] %>% dplyr::select(iter,`-2lnl`,npar) %>% distinct()
    test.list[[i]] <- null.lrt %>% inner_join(alt.lrt,by=c("iter")) %>% 
      rename(`-2lnl.null` = `-2lnl.x`,npar.null = npar.x,
             `-2lnl.alt` = `-2lnl.y`,npar.alt = npar.y) %>% 
      mutate(lrt.stat = `-2lnl.null`-`-2lnl.alt`, 
             df.lrt =  npar.alt-npar.null,
             p_value = pchisq(lrt.stat,df=df.lrt,lower.tail=FALSE),
             case=data$grid_entries[i,case]) 
  }
  
  test <- bind_rows(test.list) 
  
  
  #If we want to include samples from the theoretical distn
  if(chi_include == TRUE){
    n <- max(test$iter)
    x <- rchisq(n,unique(test$df.lrt)[1])
    p <- pchisq(x,unique(test$df.lrt)[1],lower.tail = FALSE)
    chi_sim <- data.frame("iter" = 1:n, "-2lnl.null" = NA, "npar.null" = NA,  "-2lnl.alt" = NA, "npar.alt" = NA,
                          "lrt.stat" = x,   "df.lrt" = unique(test$df.lrt)[1],  "p_value" = p, 
                          "case" = 1,stringsAsFactors = FALSE)
    rm(n,x,p)
  }
  
  p1 <- test %>% ggplot(aes(color=factor(case))) + 
    geom_density(mapping = aes(lrt.stat),key_glyph = "path") +
    labs(title = "Density of the Deviance Statistic",y="Density",x=expression("G"^2)) + 
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5),legend.position="none") +
    scale_color_manual(values = c("black","blue","purple","orange","firebrick4","red")) + 
    coord_cartesian(xlim = c(0.001,xmax1),ylim = c(-0.001,ymax1),expand = c(0))
  
  p2 <- test %>% ggplot(aes(color=factor(case))) + 
    geom_density(mapping = aes(p_value),key_glyph = "path") +
    geom_vline(xintercept = 0.05,linetype="dashed") +
    labs(title = "Density of the p-Value for the Likelihood Ratio Test",y="Density",x=expression("P(X"[1]^2*">G"^2*")")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme(legend.position="bottom") + 
    coord_cartesian(xlim = c(0.001,1.001),ylim = c(-0.001,ymax2),expand = c(0))
  
  
  #Add densities of theoretical distn to the plots
  if(chi_include == TRUE){
    p1 <- p1  + geom_density(data = chi_sim, aes(lrt.stat), col = "forestgreen",linetype = "dotted")  
    p2 <- p2  + geom_density(data = chi_sim, aes(p_value), col = "forestgreen",linetype = "dotted") + 
      scale_color_manual(name = name, values = c("black","blue","purple","orange","firebrick4","red")) 
  } else {
    p2 <- p2 +
      scale_color_manual(name = name,values = c("black","blue","purple","orange","firebrick4","red")) 
  }
  
  return(grid.arrange(p1,p2))
  
}

grid_chat_plots <- function(data,model="B",slice,main,labs=c(0.0,0.3,0.6,0.9,1.0)){
  
  data.list <- data$results[[model]]
  chat.list <- list()
  
  for(i in slice){
    
    chat.list[[i]] <- data.list[[i]] %>% dplyr::select(iter,deviance ,deviance.df) %>%
      distinct() %>% 
      mutate(c.hat = deviance/deviance.df,
             case = i)
  }
  
  chat.dat <- bind_rows(chat.list)
  
  p1 <- chat.dat %>% ggplot(aes(color=factor(case))) + 
    geom_density(mapping = aes(c.hat),key_glyph = "path") +
    geom_vline(xintercept = 1,linetype="dashed") +
    labs(title = main,y="Density",x=expression(hat(c))) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + theme(legend.position="bottom") +
    scale_color_manual(name = expression("Survival Correlation"*(gamma)),values = c("black","blue","purple","orange","firebrick4","red"),labels = labs)
  return(p1)
  
}

chat_compare <- function(data, mod = "N",main = expression("Model ("*phi*","*p*")"),label = c("Deviance","Pearson","Fletcher")){
  
  data <- pivot_longer(data,cols = c("deviance","pearson","fletcher"),names_to = "estimator")
  
  p1 <- data %>% filter(model == mod) %>% ggplot(aes(color=factor(estimator))) + 
    geom_density(mapping = aes(value),key_glyph = "path") +
    geom_vline(xintercept = 1,linetype="dashed") +
    labs(title = main,y="Density",x=expression(hat(c))) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + theme(legend.position="bottom") +
    scale_color_manual(name = expression("Estimator of "*hat(c)),
                       values = c("#F8766D","#00BA38","#619CFF"),labels = label)
  return(p1)
}

###########
# FIGURE 1#
###########
sim_surv_plots(data,"N","phi",0)

###########
# FIGURE 2#
###########
sim_surv_ci_plots(data,0)

###########
# FIGURE 3#
###########
grid_lrt_plots(data,"N","S",c(20,23,26,29,30),ymax1 = 5,ymax2 = 2.8,xmax1 = 7.5)

###########
# FIGURE 4#
###########
grid_lrt_plots(data,"N","R",c(20,23,26,29,30),ymax1 = 1.0,ymax2 = 1.25,xmax1 = 7.5)

###########
# FIGURE 5#
###########
slice <- c(80,83,86,89,90)

p1 <- grid_chat_plots(data,"N",slice,main = expression("Model ("*phi*","*p*")"))  + 
  coord_cartesian(xlim = c(0.05,9),ylim = c(-0.001,1),expand = c(0))
p2 <- grid_chat_plots(data,"R",slice,main = expression("Model ("*phi*","*p^G*")")) + 
  coord_cartesian(xlim = c(0.05,9),ylim = c(-0.001,1),expand = c(0))
p3 <- grid_chat_plots(data,"S",slice,main = expression("Model ("*phi^G*","*p*")"))  + 
  coord_cartesian(xlim = c(0.05,9),ylim = c(-0.001,1),expand = c(0))
p4 <- grid_chat_plots(data,"B",slice,main = expression("Model ("*phi^G*","*p^G*")"))  + 
  coord_cartesian(xlim = c(0.05,9),ylim = c(-0.001,1),expand = c(0))
ggarrange(p1, p2, p3, p4, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")

###########
# FIGURE 6#
###########

grid_lrt_plots(lrt_dat,"N","S",c(1,2),case = "n",name="Sample Size",ymax1 = 0.85,ymax2 = 1.15,xmax1 = 7.5)

###########
# FIGURE 7#
###########

p1 <- chat_compare(chat_dat,mod = "N",main = expression("Model ("*phi*","*p*")")) + 
  coord_cartesian(xlim = c(0.05,9),ylim = c(0,1),expand = c(0))
p2 <- chat_compare(chat_dat,mod = "R",main = expression("Model ("*phi*","*p^G*")"))  + 
  coord_cartesian(xlim = c(0.05,9),ylim = c(0,1),expand = c(0))
p3 <- chat_compare(chat_dat,mod = "S",main = expression("Model ("*phi^G*","*p*")"))  + 
  coord_cartesian(xlim = c(0.05,9),ylim = c(0,1),expand = c(0))
p4 <- chat_compare(chat_dat,mod = "B",main = expression("Model ("*phi^G*","*p^G*")"))  + 
  coord_cartesian(xlim = c(0.05,9),ylim = c(0,1),expand = c(0))
ggarrange(p1, p2, p3, p4, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")


##########
# Table 1#
##########

gamma <- c("$\\gamma=0.0$","$\\gamma=0.3$","$\\gamma=0.6$","$\\gamma=0.9$","$\\gamma=1.0$")
model <- c("$(\\phi,p)$","$(\\phi,p^G)$","$(\\phi^G,p)$","$(\\phi^G,p^G)$")
slice <- c(80,83,86,89,90)
mchat.results <- matrix(NA,ncol=5,nrow=4)
mods <- c("N","R","S","B")


compute_chat_mean <- function(data){
  chat <- data$deviance/data$deviance.df
  med.chat <- mean(chat)
  return(med.chat)
}

for(i in 1:nrow(mchat.results)){
  mod <- mods[i]
  for(j in 1:ncol(mchat.results)){
    entry <- slice[j]
    mchat.results[i,j] <- compute_chat(data$results[[mod]][[entry]])
  }
}

mchat.results <- data.frame(round(mchat.results,2)) %>% 
  mutate(Model = model) %>% 
  select(Model,X1,X2,X3,X4,X5)


kable(mchat.results,col.names = c("Model",gamma),escape = FALSE,align="c",caption = "Median($\\hat{c}$) for varying levels of ($\\gamma$) across all models") %>% 
  collapse_rows(columns = 1:5, latex_hline = "none") %>% 
  kable_styling(position="center", bootstrap_options = c("bordered"),
                font_size=11) %>%
  add_header_above(c(" " = 1, "Survival Correlation" = 5)) %>%
  kable_styling(latex_options = "HOLD_position") 

##########
# Table 2#
##########

phi <- 0.7
p <- 0.8

outcomes <- data.frame(
  Histories = c("1000","1011","1101","1110","1100","1010","1001","1111"),
  Probability = c((1-phi) + phi*(1-p)*(1-phi) + phi^2*(1-p)^2*((1-phi)+phi*(1-p)),
                  phi^3*(1-p)*p^2,
                  phi^3*(1-p)*p^2,
                  phi^2*p^2*((1-phi)+phi*(1-p)),
                  phi*p*((1-phi)+phi*(1-p)*((1-phi)+phi*(1-p))),
                  phi^2*(1-p)*p*((1-phi)+phi*(1-p)),
                  phi^3*(1-p)^2*p,
                  (phi*p)^3),
  stringsAsFactors = FALSE) %>% 
  mutate(`Expected(n=100)` = round(Probability*100,1),
         `Expected(n=200)` = round(Probability*200,1))

kable(outcomes,escape = FALSE,align="c",digits = 3, caption = "Recapture history cell probabilities and expected number of observed histories (for populations with n=100 and n=200 individuals) used in simulation study") %>% 
  collapse_rows(columns = 1, latex_hline = "none") %>% 
  kable_styling(latex_options = "HOLD_position")

##########
# Table 3#
##########

estimator <- c("Deviance","Pearson","Fletcher")
model <- c("$(\\phi,p)$","$(\\phi,p^G)$","$(\\phi^G,p)$","$(\\phi^G,p^G)$")

mchat.results <- matrix(NA,ncol=3,nrow=4)
mods <- c("N","R","S","B")

for(i in 1:nrow(mchat.results)){
  mod <- mods[i]
  for(j in 1:ncol(mchat.results)){
    mchat.results[i,j] <- median(chat_dat[chat_dat$model == mod,j+1])
  }
}

mchat.results <- data.frame(round(mchat.results,2)) %>% 
  mutate(Model = model) %>% 
  select(Model,X1,X2,X3)


kable(mchat.results,col.names = c("Model",estimator),escape = FALSE,align="c",caption = "Median($\\hat{c}$) for common estimators across all models") %>% 
  collapse_rows(columns = 1, latex_hline = "none") %>% 
  kable_styling(position="center",bootstrap_options = c("bordered"),
                font_size=11) %>%
  add_header_above(c(" " = 1, "Estimator" = 3)) %>%
  kable_styling(latex_options = "HOLD_position")
