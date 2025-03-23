## Load libraries
library(mrgsolve)    # Needed for Loading mrgsolve code into r via mcode from the 'mrgsolve' pckage
library(magrittr)    # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(dplyr)       # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(ggplot2)     # Needed for plot
library(FME)         # Package for MCMC simulation and model fitting
library(minpack.lm)  # Package for model fitting
library(reshape)     # Package for melt function to reshape the table
library(truncnorm)   # Package for the truncated normal distribution function   
library(EnvStats)    # Package for Environmental Statistics, Including US EPA Guidance
library(invgamma)    # Package for inverse gamma distribution function
library(foreach)     # Package for parallel computing
library(doParallel)  # Package for parallel computing
library(bayesplot)   # Package for MCMC traceplot
library(gridExtra)   # Package for combine ggplot
library(ggExtra)     # R-package for extend the ggplot2
library(reshape2)    # melt function to reshape the table
library(tidyverse)   # Needed for the pipe %>% operator
library(coda)        # for convergence diagnosis
library(DescTools)   # descriptive statistics and exploratory data analysis 
library(scales)      # for plotting the figure
library(grid)        # for plotting the figure
library(lattice)     # for plotting the figure

## Input mrgsolve-based PBPK Model
source (file = "8：2 RMod.R")

## Build mrgsolve-based PBPK Model
modby <- mcode ("RatPBPKby", RatPBPKby.code)

## input data set for model calibration/ oral
Data_A1    <- read.csv(file = "A1（ng ml）.csv")

# Model calibration for PBPK model based on the data of TK study                                                                                        #
# Dose regimen: oral(A), IV(B): dose to 0.01 mg/kg
# Abreviation: venous plasma (VP), liver, Kidney, Lung
# A1. : SD rat, matrix: VP, Sampling time: 2,4,8,12,24,36,48,72,120h
# A2. : SD rat, matrix: Liver, Sampling time: 2,4,8,12,24,48,72,168h  
# A3. : SD rat, matrix: Kidney, Sampling time: 2,4,8,12,24,48,72,168h    
# A4. : SD rat, matrix: Lung, Sampling time: 2,4,8,12,24,48,72,168h    
# B1. : SD rat, matrix: VP, Sampling time: 0.5,1,2,4,8,12,24,48,72h
# B2. : SD rat, matrix: Liver, Sampling time: 2,4,8,12,24,48,72,168h  
# B3. : SD rat, matrix: Kidney, Sampling time: 2,4,8,12,24,48,72,168h    
# B4. : SD rat, matrix: Lung, Sampling time: 2,4,8,12,24,48,72,168h                            
#===================================================================================================

## Read these datasets and later used in model calibration
OBS.A1  <- Data_A1 %>% filter(Study == 1 & Sample == "VP" & Dose == 0.01) %>% select(Time = "Time", CVP = "Conc")
OBS.A2  <- Data_A1 %>% filter(Study == 1 & Sample == "Liver" & Dose == 0.01) %>% select(Time = "Time", CL= "Conc")
OBS.A3  <- Data_A1 %>% filter(Study == 1 & Sample == "Kidney" & Dose == 0.01) %>% select(Time = "Time", CK = "Conc")
OBS.A4  <- Data_A1 %>% filter(Study == 1 & Sample == "Lung" & Dose == 0.01) %>% select(Time = "Time", CLu   = "Conc")

OBS.B1  <- Data_A1 %>% filter(Study == 2 & Sample == "VP" & Dose == 0.01) %>% select(Time = "Time", CVP = "Conc")
OBS.B2  <- Data_A1 %>% filter(Study == 2 & Sample == "Liver" & Dose == 0.01) %>% select(Time = "Time", CL= "Conc")
OBS.B3  <- Data_A1 %>% filter(Study == 2 & Sample == "Kidney" & Dose == 0.01) %>% select(Time = "Time", CK = "Conc")
OBS.B4  <- Data_A1 %>% filter(Study == 2 & Sample == "Lung" & Dose == 0.01) %>% select(Time = "Time", CLu   = "Conc")

## Define the prediction function
pred.rat <- function(pars) {
  
  ## Get out of log domain
  pars %<>% lapply(exp)                 ## return a list of exp (parameters) from log domain
  
  ## Define the exposure scenario
  BW          = 0.3                                       # Rat body weight
  tinterval   = 250                                       # Time interval
  TDoses      = 1                                         # Total dosing/Dose times
  PDOSEoral.A = 0.01                                      # Single oral dose from In-house experience
  DOSEoral.A  = PDOSEoral.A*BW                            # Amount of oral dose
  
  ## To create a data set of 1 subject receiving DOSE for 0,01 mg/kg
  ex.oral.A<- ev  (ID = 1,               ## One individual
                   amt  = DOSEoral.A,     ## Amount of dose 
                   ii   = tinterval,     ## Time interval
                   addl = TDoses-1,      ## Addtional doseing 
                   cmt  = "AST",         ## The dosing comaprtment: AST Stomach  
                   replicate = FALSE)    ## No replicate
  
  ex.IV.B  <- ev  (ID = 1,               ## One individual
                   amt  = DOSEoral.A,     ## Amount of dose 
                   ii   = tinterval,     ## Time interval
                   addl = TDoses-1,      ## Addtional doseing 
                   cmt  = "AVPlas_free", ## The dosing comaprtment: AST Stomach  
                   replicate = FALSE)    ## No replicate
  
  ## set up the exposure time
  tsamp=tgrid(0,tinterval*(TDoses-1)+24*10,0.5)
  
  ## Get a prediction
  out.A <- 
    modby %>%                                             # model object
    param(pars) %>%                                       # to update the parameters in the model subject
    update(atol = 1E-70, maxsteps=50000) %>%              # solver setting, atol: Absolute tolerance parameter
    mrgsim_d(data = ex.oral.A, tgrid=tsamp)               # Set up the simulation run
  out.A<-cbind.data.frame(Time=out.A$time, 
                          CVP=out.A$VP*1000,              # change μg/ml to ng/ml
                          CPlas=out.A$Plas*1000,
                          CL=out.A$Liver*1000,
                          CK=out.A$Kidney*1000,
                          CLu=out.A$Lung*1000,
                          CF=out.A$Fat*1000,
                          CR=out.A$Rest*1000,
                          BAL=out.A$Balance)
  
  out.B <- 
    modby %>%                                             # model object
    param(pars) %>%                                       # to update the parameters in the model subject
    update(atol = 1E-70, maxsteps=50000) %>%              # solver setting, atol: Absolute tolerance parameter
    mrgsim_d(data = ex.IV.B, tgrid=tsamp)                 # Set up the simulation run
  out.B<-cbind.data.frame(Time=out.B$time, 
                          CVP=out.B$VP*1000, 
                          CPlas=out.B$Plas*1000,
                          CL=out.B$Liver*1000,
                          CK=out.B$Kidney*1000,
                          CLu=out.B$Lung*1000,
                          CF=out.B$Fat*1000,
                          CR=out.B$Rest*1000,
                          BAL=out.B$Balance)
  
  out.A <- out.A %>% filter (Time > 0)                    # filter the value at time = 0
  out.B <- out.B %>% filter (Time > 0)                                 
  
  return (list("out.A"=out.A,
               "out.B"=out.B))                            # Return Dataframe
}

## initial parmaeters
theta.int <- log(c(
  KbileC                         = 0.00074,                      ## Biliary elimiatnion rate
  KurineC                        = 1.73e-5,                      ## Urinary elimination rate
  Free                           = 0.01,                       ## Free fraction in plasma
  PL                             = 8.78,                        ## Liver/plasma partition coefficient (PC)
  PK                             = 0.58,                        ## Kidney/plasma PC
  PLu                            = 0.93,                        ## Lung/plasma PC
  PF                             = 0.028,                        ## Fat/plasma PC
  PRest                          = 0.64,                        ## Rest of body/plasma PC
  K0C                            = 1,                           ## Rate of absorption of 8:2 Cl-PFESA  in the stomach
  Kabsc                          = 0.33,                        ## Rate of absorption of 8:2 Cl-PFESA  in the small intestines
  KunabsC                        = 1.85e-4,                     ## Rate of unabsorbed dose to appear in feces
  Tm                             = 84,                          ## Transport maximum
  Kt                             = 30                           ## Transport affinity constant
))
result <- pred.rat(theta.int)

plot(result$out.A$Time,result$out.A$CVP,type="l",lwd=2,xlab="Time(hour)",ylab="8:2 Cl-PFESA concentration in plasma")
plot(result$out.A$Time,result$out.A$CPlas,type="l",lwd=2,xlab="Time(hour)",ylab="8:2 Cl-PFESA concentration in plasma")
plot(result$out.A$Time,result$out.A$CL,type="l",lwd=2,xlab="Time(hour)",ylab="8:2 Cl-PFESA concentration in liver")
plot(result$out.B$Time,result$out.B$CVP,type="l",lwd=2,xlab="Time(hour)",ylab="8:2 Cl-PFESA concentration in plasma")
plot(result$out.B$Time,result$out.B$CPlas,type="l",lwd=2,xlab="Time(hour)",ylab="8:2 Cl-PFESA concentration in plasma")
plot(result$out.B$Time,result$out.B$CK,type="l",lwd=2,xlab="Time(hour)",ylab="8:2 Cl-PFESAconcentration in Kidney")
plot(result$out.A$Time,result$out.A$BAL,type="l",lwd=2,xlab="Time(hour)",ylab="Mass Balance")

## Cost fuction (FME) 
## Estimate the model residual by modCost function
MCcost<-function (pars){
  outdf <- pred.rat (pars)
  cost<- modCost  (model = outdf$out.A, obs = OBS.A1, x ="Time" ,weight = "mean")
  cost<- modCost  (model = outdf$out.A, obs = OBS.A2, x ="Time" ,weight = "mean",cost = cost)
  cost<- modCost  (model = outdf$out.A, obs = OBS.A3, x ="Time" ,weight = "mean",cost = cost)
  cost<- modCost  (model = outdf$out.A, obs = OBS.A4, x ="Time" ,weight = "mean",cost = cost)
  cost<- modCost  (model = outdf$out.B, obs = OBS.B1, x ="Time" ,weight = "mean",cost = cost)
  cost<- modCost  (model = outdf$out.B, obs = OBS.B2, x ="Time" ,weight = "mean",cost = cost)
  cost<- modCost  (model = outdf$out.B, obs = OBS.B3, x ="Time" ,weight = "mean",cost = cost)
  cost<- modCost  (model = outdf$out.B, obs = OBS.B4, x ="Time" ,weight = "mean",cost = cost)
  return(cost)
}

## Senstivity function (FME) 
## Check the senstive parameters in the model
SnsPlasma <- sensFun(func = MCcost, parms = theta.int, varscale = 1)
Sen       <- summary(SnsPlasma)
plot(Sen)

## Selected senstive parameters
theta <- theta.int[abs(Sen$Mean) > 1.2*mean(abs(Sen$Mean))]
theta 

## initial parmaeters
theta.int <- log(c(
  #KbileC                         = 0.00074,                      ## Biliary elimiatnion rate
  #KurineC                        = 1.73e-5,                      ## Urinary elimination rate
  Free                           = 0.01,                       ## Free fraction in plasma
  #PL                             = 8.78,                        ## Liver/plasma partition coefficient (PC)
  PK                             = 0.58,                        ## Kidney/plasma PC
  PLu                            = 0.93,                        ## Lung/plasma PC
  PF                             = 0.028,                        ## Fat/plasma PC
  PRest                          = 0.64,                         ## Rest of body/plasma PC
  K0C                            = 1,                          ## Rate of absorption of 8:2 Cl-PFESA  in the stomach
  Kabsc                          = 0.33                        ## Rate of absorption of 8:2 Cl-PFESA  in the small intestines
  #KunabsC                        = 1.85e-4                     ## Rate of unabobded dose to appear in feces
  #Tm                             = 84,                          ## Transport maximum
  #Kt                             = 30                           ## Transport affinity constant
))

## PBPK model fitting 
## least squres fit using Nelder-Mead algorithm

Fit<- modFit(f=MCcost, p=theta.int, method ="Nelder-Mead",
             control = nls.lm.control(nprint=1))

summary(Fit)                                 ## Summary of fit 
exp(Fit$par)                                 ## Get the arithmetic value out of the log domain
Cost <- MCcost(Fit$par)

## Save the fitting results to RDS files
saveRDS(Fit, file = "Fit_C10.rds") 

FDataA <- cbind.data.frame (name= Cost$residuals$name,
                            OBS = Cost$residuals$obs, 
                            PRE = Cost$residuals$mod)

## Transformed the predicted and obseved values using log10-sacle to do the plot
FDataA %<>% mutate (Log.OBS = log(OBS,10), Log.PRE = log(PRE,10), Species = "Rat")

## Estimating the R-squared and goodness-of-fit using linear regression model
fit <- lm(Log.OBS ~Log.PRE, data = FDataA)
summary(fit)

FDataA %<>% mutate(res = residuals(fit), 
                   prediction = predict(fit), 
                   OPR = PRE/OBS,            ## OPR: the ratio of prediction value and observed data
                   log.OPR =  log(OPR,10),
                   ID = 1:nrow(FDataA),
                   Label = "C10") 
p <- 
  ggplot(FDataA, aes(Log.OBS, Log.PRE)) +    ## using log-sacle axis
  geom_abline (intercept = 0, 
               slope     = 1,
               color     ="black", linetype = "dashed",linewidth = 1, alpha = 0.8) +
  geom_point  (aes(shape   = as.factor(name)),size = 4)  +
  annotation_logticks() +
  scale_y_continuous(limits = c(0,4), labels = scales::math_format(10^.x))+
  scale_x_continuous(limits = c(0,4),labels = scales::math_format(10^.x))

## Set up your theme and font
windowsFonts("Times" = windowsFont("Times New Roman"))

p1 <- p + 
  theme (
    plot.background         = element_rect (fill="White"),
    text                    = element_text (family = "Times"),   # text front (Time new roman)
    panel.border            = element_rect (colour = "black", fill=NA, linewidth =2),
    panel.background        = element_rect (fill="White"),
    panel.grid.major        = element_blank(),
    panel.grid.minor        = element_blank(), 
    axis.text               = element_text (size   = 15, colour = "black", face = "bold"),    # tick labels along axes 
    axis.title              = element_text (size   = 18, colour = "black", face = "bold"),   # label of axes
    legend.position         ='none') +
  labs (x = "",  y = "")

p2 <-
  ggplot(FDataA, aes(Log.PRE, log.OPR)) +
  geom_hline(yintercept = log10(2),linetype = 3, color   = "black", size =1) +
  geom_hline(yintercept = log10(0.5),linetype = 3, color   = "black", size =1) +
  geom_hline(yintercept = log10(3),linetype = 3, color   = "grey", size =1) +
  geom_hline(yintercept = log10(0.33),linetype = 3, color   = "grey", size =1) +
  geom_point(color   = "steelblue4",size = 3) +
  annotation_logticks() +
  scale_y_continuous(limits = c(-2,2), labels = scales::math_format(10^.x))+
  scale_x_continuous(limits = c(-1,3),labels = scales::math_format(10^.x))

p2 <- p2 +
  theme (
    plot.background         = element_rect (fill ="White"),
    text                    = element_text (family = "Times"),   # text front (Time new roman)
    panel.border            = element_rect (colour = "black", fill=NA, linewidth =2),
    panel.background        = element_rect (fill ="White"),
    panel.grid.major        = element_blank(),
    panel.grid.minor        = element_blank(), 
    axis.text               = element_text (size   = 15, colour = "black", face = "bold"),    # tick labels along axes 
    axis.title              = element_text (size   = 18, colour = "black", face = "bold"),   # label of axes
    legend.position='none') +
  labs (x = "", y = "")

grid.arrange(p1,p2,ncol = 2)

FDataA$ID <- 1:nrow(FDataA)
p3 <- ggplot(FDataA, aes(ID, log.OPR)) + geom_point() + theme_bw()+
  scale_y_continuous(limits = c(-2, 2), labels = scales::math_format(10^.x)) +
  annotation_logticks()

p3 <- ggMarginal(p3, type = "histogram", margins = "y", size = 4,
                 color="black", alpha = 0.4, bins = 50,
                 fill = "steelblue4", position="identity")
print(p3)

write.csv(FDataA_C10,"FDataA_C10.csv")

## Model calibration plot using ggplot2 
Sim.fit.A = pred.rat(Fit$par)
df.sim.A1 <- cbind.data.frame(Time = Sim.fit.A$out.A$Time,
                              CVP = Sim.fit.A$out.A$CVP,
                              CL = Sim.fit.A$out.A$CL,
                              CK = Sim.fit.A$out.A$CK,
                              CLu = Sim.fit.A$out.A$CLu)

df.sim.B1 <- cbind.data.frame(Time = Sim.fit.A$out.B$Time,
                              CVP = Sim.fit.A$out.B$CVP,
                              CL = Sim.fit.A$out.B$CL,
                              CK = Sim.fit.A$out.B$CK,
                              CLu = Sim.fit.A$out.B$CLu)

## Setting an initial value
df.sim.A1 <- rbind(data.frame(Time = 0, CVP = 1e-6, CL = 1e-6, CK = 1e-6, CLu = 1e-6),df.sim.A1)
df.sim.B1 <- rbind(data.frame(Time = 0, CVP = 1e-6, CL = 1e-6, CK = 1e-6, CLu = 1e-6),df.sim.B1)

OBS.A1.1  <- Data_A1 %>% filter(Study == 1 & Sample == "VP" & Dose == 0.01) %>% select(Time = "Time", CVP = "Conc", SD = "SD")
OBS.A2.1  <- Data_A1 %>% filter(Study == 1 & Sample == "Liver" & Dose == 0.01) %>% select(Time = "Time", CL= "Conc", SD = "SD")
OBS.A3.1  <- Data_A1 %>% filter(Study == 1 & Sample == "Kidney" & Dose == 0.01) %>% select(Time = "Time", CK = "Conc", SD = "SD")
OBS.A4.1  <- Data_A1 %>% filter(Study == 1 & Sample == "Lung" & Dose == 0.01) %>% select(Time = "Time", CLu   = "Conc", SD = "SD")

OBS.B1.1  <- Data_A1 %>% filter(Study == 2 & Sample == "VP" & Dose == 0.01) %>% select(Time = "Time", CVP = "Conc", SD = "SD")
OBS.B2.1  <- Data_A1 %>% filter(Study == 2 & Sample == "Liver" & Dose == 0.01) %>% select(Time = "Time", CL= "Conc", SD = "SD")
OBS.B3.1  <- Data_A1 %>% filter(Study == 2 & Sample == "Kidney" & Dose == 0.01) %>% select(Time = "Time", CK = "Conc", SD = "SD")
OBS.B4.1  <- Data_A1 %>% filter(Study == 2 & Sample == "Lung" & Dose == 0.01) %>% select(Time = "Time", CLu   = "Conc", SD = "SD")

plot.A1 <- ggplot() +
  geom_line(data = df.sim.A1, aes(Time, CVP), col = "#5A5A5A", lwd = 0.8) +
  geom_point(data = OBS.A1, aes(Time, CVP), col = "#e3716e", size = 2.5) +
  geom_errorbar(data = OBS.A1.1, aes(x = Time, ymin = CVP - SD, ymax = CVP + SD), col = "#e3716e", width = 2, size = 0.8) +
  ylab("Concentration in the plasma (ng/mL)")+
  xlab("Time (h)") +
  xlim(c(0, 125)) +
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2))

plot.A2 <- ggplot() +
  geom_line(data = df.sim.A1, aes(Time, CL), col = "#5A5A5A", lwd = 0.8) +
  geom_point(data = OBS.A2, aes(Time, CL), col = "#e3716e", size = 2.5) + 
  geom_errorbar(data = OBS.A2.1, aes(x = Time, ymin = CL - SD, ymax = CL + SD), col = "#e3716e", width = 2, size = 0.8) +
  ylab("Concentration in the liver (ng/g)")+
  xlab("Time (h)") + 
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2))

plot.A3 <- ggplot() +
  geom_line(data = df.sim.A1, aes(Time, CK), col = "#5A5A5A", lwd = 0.8) +
  geom_point(data = OBS.A3, aes(Time, CK), col = "#e3716e", size = 2.5) + 
  geom_errorbar(data = OBS.A3.1, aes(x = Time, ymin = CK - SD, ymax = CK + SD), col = "#e3716e", width = 2, size = 0.8) +
  ylab("Concentration in the kidney (ng/g)")+
  xlab("Time (h)") + 
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2))

plot.A4 <- ggplot() +
  geom_line(data = df.sim.A1, aes(Time, CLu), col = "#5A5A5A", lwd = 0.8) +
  geom_point(data = OBS.A4, aes(Time, CLu), col = "#e3716e", size = 2.5) + 
  geom_errorbar(data = OBS.A4.1, aes(x = Time, ymin = CLu - SD, ymax = CLu + SD), col = "#e3716e", width = 2, size = 0.8) +
  ylab("Concentration in the lung (ng/g)")+
  xlab("Time (h)") + 
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2))

plot.B1 <- ggplot() +
  geom_line(data = df.sim.B1, aes(Time, CVP), col = "#5A5A5A", lwd = 0.8) +
  geom_point(data = OBS.B1, aes(Time, CVP), shape = 17, col = "#4A90E2", size = 2.5) +
  geom_errorbar(data = OBS.B1.1, aes(x = Time, ymin = CVP - SD, ymax = CVP + SD), col = "#4A90E2", width = 2, size = 0.8) +
  ylab("Concentration in the plasma (ng/mL)")+
  xlab("Time (h)") +
  xlim(c(0, 80)) +
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2))

plot.B2 <- ggplot() +
  geom_line(data = df.sim.B1, aes(Time, CL), col = "#5A5A5A", lwd = 0.8) +
  geom_point(data = OBS.B2, aes(Time, CL), shape = 17, col = "#4A90E2", size = 2.5) + 
  geom_errorbar(data = OBS.B2.1, aes(x = Time, ymin = CL - SD, ymax = CL + SD), col = "#4A90E2", width = 2, size = 0.8) +
  ylab("Concentration in the liver (ng/g)")+
  xlab("Time (h)") + 
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2))

plot.B3 <- ggplot() +
  geom_line(data = df.sim.B1, aes(Time, CK), col = "#5A5A5A", lwd = 0.8) +
  geom_point(data = OBS.B3, aes(Time, CK), shape = 17, col = "#4A90E2", size = 2.5) + 
  geom_errorbar(data = OBS.B3.1, aes(x = Time, ymin = CK - SD, ymax = CK + SD), col = "#4A90E2", width = 2, size = 0.8) +
  ylab("Concentration in the kidney (ng/g)")+
  xlab("Time (h)") + 
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2))

plot.B4 <- ggplot() +
  geom_line(data = df.sim.B1, aes(Time, CLu), col = "#5A5A5A", lwd = 0.8) +
  geom_point(data = OBS.B4, aes(Time, CLu), shape = 17, col = "#4A90E2", size = 2.5) + 
  geom_errorbar(data = OBS.B4.1, aes(x = Time, ymin = CLu - SD, ymax = CLu + SD), col = "#4A90E2", width = 2, size = 0.8) +
  ylab("Concentration in the lung (ng/g)")+
  xlab("Time (h)") + 
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2))

grid.arrange(plot.A1, plot.B1, plot.A2, plot.B2, plot.A3, plot.B3, plot.A4, plot.B4, ncol = 2)

############################################# Optimization with MCMC Analaysis ####################################################
## Fixed the physiological parameters;
## Input a new initial parameters
## Population mean and model error (sig2)
theta.MCMC<-log(c(
  KbileC                         = 0.00074,                      ## Biliary elimiatnion rate
  KurineC                        = 1.73e-5,                      ## Urinary elimination rate
  Free                           = 0.009,                       ## Free fraction in plasma
  PL                             = 8.78,                        ## Liver/plasma partition coefficient (PC)
  PK                             = 0.72,                        ## Kidney/plasma PC
  PLu                            = 1.49,                        ## Lung/plasma PC
  PF                             = 0.064,                        ## Fat/plasma PC
  PRest                          = 1.90,                         ## Rest of body/plasma PC
  K0C                            = 0.35,                          ## Rate of absorption of 8:2 Cl-PFESA  in the stomach
  Kabsc                          = 0.077,                        ## Rate of absorption of 8:2 Cl-PFESA  in the small intestines
  KunabsC                        = 1.85e-4,                     ## Rate of unabsorbed dose to appear in feces
  Tm                             = 84,                          ## Transport maximum
  Kt                             = 30,                           ## Transport affinity constant
  sig2                           = 0.5,                         ## Model error (residuals); mostly between 0.3 and 0.5 (corresponding to coefficients of variation of about 30-50%); Bois et al. (1998, 2000)
  
  ## population variance; equal to the CV of parameters (this study assumed the cv of most parameters is 0.3, some parameters is 0.5 due to possible variation)
  sig_KbileC                     = 0.5,                         ## Default value of 0.3 and 0.5 was used to represent a moderate level of variation (Hack et al., 2006; Chiu et al., 2009)
  sig_KurineC                    = 0.5,
  sig_TmC                        = 0.5,
  sig_KtC                        = 0.5,
  sig_Free                       = 0.3,
  sig_PL                         = 0.3,
  sig_PK                         = 0.3,
  sig_PLu                        = 0.3,
  sig_PF                         = 0.3,
  sig_PRest                      = 0.3,
  sig_K0C                        = 0.3,
  sig_Kabsc                      = 0.3,
  sig_KunabsC                    = 0.3
))
which_sig <- grep("sig", names(theta.MCMC))           ## Get a list of parameters with "sig" character

## Maximum likelihood estimation (MLE) fuction for MCMC
mcmc.fun <- function (pars, pred=FALSE){
  
  ## Get out of log domain
  pars.data <- pars [-which_sig] %>% lapply(exp)      ## Return a list of exp (parameters) from log scale

  BW               = 0.30                              ## kg, Rat body weight                                      
  tinterval        = 240                               ## hr, Time interval                                 
  TDoses           = 1                                ## Dose times                                    

  PDOSEoral.A        = 0.01                                 ## mg/kg-d, Single Oral dose                               
  DOSEoral.A         = PDOSEoral.A*BW                      ## mg, amount of oral dose              
  ex.oral.A          <- ev(ID =1, 
                           amt   = DOSEoral.A, 
                           ii    =tinterval, 
                           addl  =TDoses-1, 
                           cmt   ="AST", 
                           replicate = FALSE)
  
  PDOSEiv.B        = 0.01                                 ## mg/kg-d, Single iv dose                               
  DOSEiv.B         = PDOSEiv.B*BW                      ## mg, amount of iv dose              
  ex.iv.B          <- ev(ID =1, 
                         amt   = DOSEiv.B, 
                         ii    =tinterval, 
                         addl  =TDoses-1, 
                         cmt   ="AVPlas_free", 
                         replicate = FALSE)
  
  # set up the exposure time
  tsamp = tgrid (0,tinterval*(TDoses-1)+24*10,1)       ## simulation 24*10 hr (10 days)
  
  ## Simulation of exposure scenaior A (oral single dose to 0.01 mg/kg)
  out.A <- 
    modby %>% 
    param (pars.data) %>%
    update(atol = 1E-8, maxsteps = 5000) %>%          ## atol: Absolute tolerance parameter
    mrgsim_d (data = ex.oral.A, tgrid=tsamp) %>%
    filter(time!=0)
  
  out.A1 = cbind.data.frame(Time=out.A$time, 
                            CVP=out.A$VP*1000)
  out.A2 = cbind.data.frame(Time=out.A$time, 
                            CL=out.A$Liver*1000,
                            CK=out.A$Kidney*1000,
                            CLu=out.A$Lung*1000)
  
  ## Simulation of exposure scenaior B (iv single dose to 0.01 mg/kg)
  out.B <- 
    modby %>% 
    param (pars.data) %>%
    update(atol = 1E-8, maxsteps = 5000) %>%          ## atol: Absolute tolerance parameter
    mrgsim_d (data = ex.iv.B, tgrid=tsamp) %>%
    filter(time!=0)
  
  out.B1 = cbind.data.frame(Time=out.B$time, 
                            CVP=out.B$VP*1000)
  out.B2 = cbind.data.frame(Time=out.B$time, 
                            CL=out.B$Liver*1000,
                            CK=out.B$Kidney*1000,
                            CLu=out.B$Lung*1000) 
  
  if (pred) return(list("out.A" = out.A,     
                        "out.B" = out.B))   
  
  ## Get the same time-course profiles with experiment data from the simulation
  out.A1 = out.A1[which(out.A1$Time %in% OBS.A1$Time),]    
  out.A2 = out.A2[which(out.A2$Time %in% OBS.A2$Time),]    
  out.B1 = out.B1[which(out.B1$Time %in% OBS.B1$Time),]
  out.B2 = out.B2[which(out.B2$Time %in% OBS.B2$Time),]
  
  ## log.Predition 
  log.yhat.oral1.CV     <-log(out.A1$CVP)
  log.yhat.oral1.CL     <-log(out.A2$CL)
  log.yhat.oral1.CK     <-log(out.A2$CK)
  log.yhat.oral1.CLu     <-log(out.A2$CLu)
  
  log.yhat.iv.CV     <-log(out.B1$CVP)
  log.yhat.iv.CL     <-log(out.B2$CL)
  log.yhat.iv.CK     <-log(out.B2$CK)
  log.yhat.iv.CLu     <-log(out.B2$CLu)
  
  ## log.Observed data
  
  log.y.oral1.CV     <-log(OBS.A1$CVP)
  log.y.oral1.CL     <-log(OBS.A2$CL)
  log.y.oral1.CK     <-log(OBS.A3$CK)
  log.y.oral1.CLu     <-log(OBS.A4$CLu)
  
  log.y.iv.CV     <-log(OBS.B1$CVP)
  log.y.iv.CL     <-log(OBS.B2$CL)
  log.y.iv.CK     <-log(OBS.B3$CK)
  log.y.iv.CLu     <-log(OBS.B4$CLu)
  
  ## The method of Maximum likelihood
  log.yhat        <- c(log.yhat.oral1.CV, log.yhat.oral1.CL, log.yhat.oral1.CK, log.yhat.oral1.CLu, log.yhat.iv.CV, log.yhat.iv.CL,
                       log.yhat.iv.CK, log.yhat.iv.CLu)
  log.y           <- c(log.y.oral1.CV,log.y.oral1.CL, log.y.oral1.CK, log.y.oral1.CLu,log.y.iv.CV,log.y.iv.CL,
                       log.y.iv.CK,log.y.iv.CLu)
  
  sig2            <- as.numeric((exp(pars[which_sig][1]))) # Get the parameter of sig2 from the parameters
  
  log_likelihood  <- -2*sum ((dnorm (log.y,
                                     mean =log.yhat,
                                     sd = sqrt(sig2),
                                     log = TRUE)))
  return(log_likelihood)
  
}

log_likelihood <- mcmc.fun(theta.MCMC)

## Define the Prior distributions: either normal or lognormal distribution
## nomral distribution
Prior <- function(pars) {
  
  ## Population level
  # The likelihood for population mean (parameters)
  pars.data = exp(pars[-which_sig])
  sig  <- as.numeric (exp(pars[which_sig][2:14]))                 # Coefficient of variation from population variance; sigmal0
  sig2 <- as.numeric (exp(pars[which_sig][1]))                    # error variances from model residual
  
  mean           = exp(theta.MCMC[-which_sig])
  CV             = 0.5                                            # Coefficient of variation; Default value of 0.5 in all parameters (Bois,2000; Bois et al., 1996)
  sd             = mean*CV
  
  # Calculate likelihoods of each parameters; P(u|M,S)
  prior_pars     = dtruncnorm(pars.data, 
                              a = qnorm(0.025, mean = mean, sd = sd), 
                              b = qnorm(0.975, mean = mean, sd = sd), 
                              mean = mean, sd = sd ) 
  
  # The likelihood for population variance; P(sigmal^2|sigmal0^2)
  CU             = 1                                              # Coefficient of uncertainty (CU) (Hack et al., 2006)
  CV.sig         = exp(theta.MCMC[which_sig])[2:14]               # Singmal0
  alpha          = (2+1)/(CU^2)                                   # Shape parametrer of gamma distribution; Appendix Table A-7 from EPA (2011) :EPA/635/R-09/011F
  beta           = (alpha-1)*CV.sig^2                             # Scale parameter  of gamma distribution; Appendix Table A-7 from EPA (2011) :EPA/635/R-09/011F
  
  # Calculate likelihoods of model error (sig2) and population variance (sig) parameters
  prior_sig      = dinvgamma (sig, shape = alpha , scale = beta)  # Prior distribution for population vraicne; sigma2
  prior_sig2     = dunif (sig2, min = 0.01, max = 3.3)            # Error variances, Lower and upper boundary from Chiu et al., 2009; Chiu et al., 2014)   
  
  ## individual level; P(theta|u,sigmal^2)
  mean_i         = prior_pars
  sd_i           = sqrt(prior_sig)
  prior_pars_i   = dtruncnorm (prior_pars, 
                               a = qnorm(0.025, mean = mean_i, sd = sd_i), 
                               b = qnorm(0.975, mean = mean_i, sd = sd_i), 
                               mean = mean_i, sd = sd_i) 
  
  # log-transformed (log-likelihoods of each parameters)
  log.pri.pars   = log (prior_pars)
  log.pri.sig    = log (prior_sig)
  log.pri.pars.i = log (prior_pars_i)
  log.pri.sig2   = log (prior_sig2)
  
  # maximau likelihood estimation (MLE): negative log-likelihood function, (-2 times sum of log-likelihoods)
  MLE =  -2*sum(log.pri.pars, log.pri.sig , log.pri.pars.i,log.pri.sig2)  
  
  return(MLE)
}

#################### MCMC simulation with parallel computing ############################
detectCores()                                ## check the cores
cl<- makeCluster(detectCores())              ## use all cores in our system     
registerDoParallel(cl)                       ## registers a cluster of all the cores on our system

# start time
strt<-Sys.time()

# parallel
#need a Supercomputing Center 
system.time(
  MCMC <- foreach( i = 1:4, .packages = c('mrgsolve','magrittr','FME','truncnorm','EnvStats','invgamma','dplyr')) %dopar% {
    modby <- mcode ("RatPBPKby", RatPBPKby.code)
    modMCMC(f             = mcmc.fun, 
            p             = theta.MCMC, 
            niter         = 500000,           ## iteration number 
            jump          = 0.01,             ## jump function generation new parameters distribution using covariate matrix
            prior         = Prior,            ## prior function
            updatecov     = 50,               ## Adaptative Metropolis
            var0          = NULL,             ## initial model variance;
            wvar0         = 0.01,             ## "Weight" for the initial model variance
            ntrydr        = 2,                ## Delayed Rejection
            burninlength  = 250000,           ## number of initial iterations to be removed from output.
            outputlength  = 50000)            ## number of output iterations           
  }
)

#end time
print(Sys.time()-strt)

stopCluster(cl)   

## Performance four chains to check the convergences
MC.rat.1 = as.mcmc (MCMC[[1]]$pars)   # first  chain
combinedchains = mcmc.list(MC.rat.1,MC.rat.2,MC.rat.3,MC.rat.4) ## combine all chains
diag<-gelman.diag (combinedchains)          # Gelman convergence diagnosis
diagpsrf <- diag$psrf

## Save the Posterior parameters (95% CI)
## Median (2.5%, 97.5%) for Posterior parameters of population mean
quan.rat = exp(summary(MC.rat.1)$quantiles)  

## Trace plot using bayesplot
## Covergences plot
mcmc_trace (
  combinedchains,
  pars       = names(theta.MCMC[1:19]),
  size       = 0.5,
  facet_args = list(nrow = 4)) +
  ggplot2::scale_color_brewer()

## output the MCMC results with .csv and .rds
write.csv(quan.rat,file="rat.summary_pos.csv")
#write.csv(MC.rat.1,file="rat.pos.csv")
saveRDS(MCMC[[1]],file ='rat.MCMC.rds')
saveRDS(names(theta.MCMC),file ="theta.names.rds")

rat.MCMC <- readRDS(file = "rat.MCMC.rds")
theta.names <- readRDS(file = "theta.names.rds")
which_sig   <- grep("sig", theta.names)
# Median (2.5%, 97.5%) of Prior distribution and Posterior distribution

## Prior parameters
## M value: Mean of population mean u
## S vaelu: SD of population mean u: set the CV equal to 50%
mean    = (exp(theta.MCMC[-which_sig])) # M-value for rat
sd      = mean*0.5                       # S-value for Rat

## 50%, 2.5% and 97.5% for prior population mean u

a.rat    = qnorm(0.025, mean = mean, sd = sd)
b.rat    = qnorm(0.975, mean = mean, sd = sd)
c.rat    = qnorm(0.5, mean = mean, sd = sd)

prior <-  cbind(theta.names[-which_sig], c.rat, a.rat,b.rat)
colnames(prior) <- c("parameter", "50%", "2.5%", "97.5%")

## Evaluation
pred.ratE <- function (pars.rat) {
  
  ## Exposure scenario for rat oral daily exposue to 0.01 mg/kg
  
  BW.ratA                = 0.3                                      ## Rat body weight
  tinterval.rat          = 240                                        ## Time interval
  TDoses.rat             = 1                                        ## times, Dose times
  
  PDOSE.rat.a            = 0.01                                      ## mg/kg, BW Oral dose
  DOSEoral.rat.a         = PDOSE.rat.a*BW.ratA                                    ## mg, amount of oral dose
  DOSEIV.rat.b           = PDOSE.rat.a*BW.ratA                                    ## mg, amount of oral dose
  ex.rat.a               <- ev(ID=1, amt= DOSEoral.rat.a,
                               ii=tinterval.rat, addl=TDoses.rat-1,cmt="AST",replicate = FALSE)
  ex.rat.b               <- ev(ID=1, amt= DOSEIV.rat.b,
                               ii=tinterval.rat, addl=TDoses.rat-1,cmt="AVPlas_free",replicate = FALSE)
  tsamp.rat              = tgrid(0,tinterval.rat*(TDoses.rat-1)+24*10,1)  
  
  ## Rat ouput
  pars.rat %<>% lapply(exp)
  names(pars.rat) <- names(pars.rat)
  pars.rat        <- pars.rat[-which_sig ]
  
  outdf.rat.a <- 
    modby %>% 
    param(pars.rat) %>%
    update(atol = 1E-5,rtol= 1e-5, maxsteps = 2000) %>% # 2000 is big enough to ensure stable simulations. rtol is the key factor that determines the speed. The smaller of rtol, the longer of the simulation time.
    mrgsim_d (data = ex.rat.a, tgrid = tsamp.rat)
  
  outdf.rat.b <- 
    modby %>% 
    param(pars.rat) %>%
    update(atol = 1E-5,rtol= 1e-5, maxsteps = 2000) %>%
    mrgsim_d (data = ex.rat.b, tgrid = tsamp.rat)
  
  outdf.rat.a <- cbind.data.frame(Time   = outdf.rat.a$time, 
                                  CR     = outdf.rat.a$Rest*1000)      
  
  outdf.rat.b <- cbind.data.frame(Time   = outdf.rat.b$time, 
                                  CR     = outdf.rat.b$Rest*1000)      
  
  return (list("outdf.rat.a"   = outdf.rat.a, 
               "outdf.rat.b"   = outdf.rat.b))
}

Newtime.r   = pred.ratE(theta.MCMC[-which_sig])$outdf.rat.a$Time  # this is the new time variable.
nrwo.r = length (Newtime.r)

# Create the matrix 
MC.rat.a.CR    = matrix(nrow = nrwo.r, ncol = 5000) 
MC.rat.b.CR    = matrix(nrow = nrwo.r, ncol = 5000)

## Input paramters
for(i in 1:500){
  
  j = i *10  
  pars.rat             = rat.MCMC$pars    [j,]     # sample parameter set once every ten sets, so you will have 5000 sets from 50000 total sets
  
  MCdata               = pred.ratE (pars.rat)
  MC.rat.a.CR    [,i]  = MCdata $outdf.rat.a$CR    
  MC.rat.b.CR    [,i]  = MCdata $outdf.rat.b$CR
  cat("iteration = ", i , "\n") 

# get the simulation of Rest
M.rat.a.CR  = MC.rat.a.CR %>% apply(1,mean, na.rm = TRUE)
SD.rat.a.CR = MC.rat.a.CR %>% apply(1,sd, na.rm = TRUE)
M.rat.b.CR  = MC.rat.b.CR %>% apply(1,mean, na.rm = TRUE)
SD.rat.b.CR = MC.rat.b.CR %>% apply(1,sd, na.rm = TRUE)

NewMC.rat.a.CR = cbind.data.frame(
  Time = Newtime.r, Mean = M.rat.a.CR, SD = SD.rat.a.CR) # Create this object to plot the histogram. This dataset only contains predicted values. The observed values are below.

NewMC.rat.b.CR = cbind.data.frame(
  Time = Newtime.r, Mean = M.rat.b.CR, SD = SD.rat.b.CR)

NewMC.rat.a.CR %<>% filter(Time == 168) # only obtain the predicted value at 168 days for the histogram.
NewMC.rat.a.CR$Dose = c("0.01")
NewMC.rat.b.CR %<>% filter(Time == 168)
NewMC.rat.b.CR$Dose = c("0.011")

Rat.CR <- rbind (NewMC.rat.a.CR, NewMC.rat.b.CR) # combine the data points.
Rat.CR$Matrix <- c("Pre.Rest") # create a matrix column as the predicted Rest values.

## Observed data from In house TK data
Obs.Time      = c(168,168)                          # days,Treatment 98 days
Obs.Dose      = c("0.01","0.011")                 # double quotation marks are not necessary, but it is OK to use it as it could be used as the x axis label.
Obs.rat.CR.M  = c(2.08,1.07)                  # ug/g, mean of PFOS concentration in Plasma for male rat after dosing to 14 weeks; Obtained from tale 2 (Seacat et al., 2003)
Obs.rat.CR.SD = c(1.15,0.62)                # ug/g, SD of PFOS concentration in Plasma for male rat after dosing to 14 weeks; Obtained from tale 2 (Seacat et al., 2003)
Obs.rat.CR    = cbind.data.frame (Time = Obs.Time, Mean = Obs.rat.CR.M, SD = Obs.rat.CR.SD, Dose = Obs.Dose, Matrix = c("Obs.Rest")) # Create a matrix of obs.plasma just to categorize the data.
Comb.Rat = rbind (Rat.CR,Obs.rat.CR)

## Histogram plot
p1.E<- 
  ggplot(Comb.Rat, aes(x=as.factor(Dose), y = Mean, fill= as.factor(Matrix))) + # do the histogram plot by the factor of matrix, so you have a bar for pre.plasma, a bar for obs.plasma, so on.
  geom_bar(stat="identity", color="black", size = 1.2, 
           position=position_dodge()) + 
  scale_fill_manual(values=c("#999999","#E69F00"))+
  geom_errorbar(aes(ymin=Mean, ymax=Mean+SD), width=.2,
                position=position_dodge(.9))+
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2))

###################################NCA##############################################
### Sensitivity analysis
Rat.theta <- log(c(
  # Physiological parameters
  BW                             = 0.3,
  Htc                            = 0.46,
  QCC                            = 14.1,
  QLC                            = 0.183,
  QKC                            = 0.141,
  QLuC                           = 1,   
  QFC                            = 0.07,
  VLC                            = 0.037,
  VKC                            = 0.0073,
  VLuC                           = 0.005, 
  VFC                            = 0.07,
  VPlasC                         = 0.0312,
  # Chemical-specific parameters (FME-fitting values)
  KbileC                         = 0.00074,                     
  KurineC                        = 1.73e-5,                      
  Free                           = 0.008,                             ### fitting parameters
  PL                             = 8.78,                       
  PK                             = 0.72,                              ### fitting parameters
  PLu                            = 1.49,                              ### fitting parameters
  PF                             = 0.06,                              ### fitting parameters
  PRest                          = 1.90,                              ### fitting parameters
  K0C                            = 0.35,                              ### fitting parameters
  Kabsc                          = 0.077,                             ### fitting parameters
  KunabsC                        = 1.85e-4,                    
  Tm                             = 84,                        
  Kt                             = 30                         
))

pred.NSC <- function(pars) {
  
  ## Get out of log domain
  pars <- exp(pars)                   ## Return a list of exp parametrs from log scale
  
  ## Exposure scenario for gestational exposure
  BW          = 0.3                                       # Rat body weight
  tinterval   = 240                                       # Time interval
  TDoses      = 1                                         # Total dosing/Dose times
  PDOSEoral.A = 0.01                                      # Single oral dose from TK experience
  DOSEoral.A  = PDOSEoral.A*BW                            # Amount of oral dose
  
  # To create a exposure scenario
  ex.oral <- ev (ID   = 1,             ## One individual
                 amt  = DOSEoral.A,    ## Amount of dose 
                 ii   = tinterval,     ## Time interval
                 addl = TDoses-1,      ## Addtional doseing 
                 cmt  = "AST",         ## The dosing comaprtment: AST Stomach  
                 replicate = FALSE)    ## No replicate
  
  ## set up the exposure time
  tsamp=tgrid(0,tinterval*(TDoses-1)+24*7,1)
  
  ## Simulation of exposure scenaior (oral dose to 10 μg/kg)
  out <- 
    modby %>%
    param (pars) %>%
    update(atol = 1E-6, maxsteps = 500000) %>%          
    mrgsim_d (data = ex.oral, tgrid = tsamp)
  
  outdf = cbind.data.frame (Time   = out$time, 
                            CPlas = out$Plas, 
                            AUCCV = out$AUC_CV,
                            AUCCL = out$AUC_CL,
                            AUCCK = out$AUC_CK,
                            AUCCLu = out$AUC_CLu,
                            AUCCF = out$AUC_CF)
  
  outdf <- outdf %>% filter (Time == 7*24) ## The model output on last days
  
  
  return (list("R" = outdf))
}

result  <-  pred.NSC(Fit$par)

NSC_func <- function (pars, Pred) {
  n <- length(pars)
  NSC_R     = matrix(NA, nrow = length(pars) , ncol = 5)
  
  for (i in 1:n) {
    pars.new       <- pars %>% replace(i, log(exp((pars[i]))*1.01))
    Mnew.G         <- Pred(pars.new)
    M.G            <- Pred(pars)
    delta.pars     <- exp(pars[i])/(exp(pars[i])*0.01)
    
    ## Estimated the AUC
    AUC.Plas.new       =  Mnew.G$R %>% select (AUCCV)
    AUC.Plas.ori       =  M.G$R    %>% select (AUCCV)
    AUC.L.new          =  Mnew.G$R %>% select (AUCCL)
    AUC.L.ori          =  M.G$R    %>% select (AUCCL)
    AUC.K.new          =  Mnew.G$R %>% select (AUCCK)
    AUC.K.ori          =  M.G$R    %>% select (AUCCK)
    AUC.Lu.new         =  Mnew.G$R %>% select (AUCCLu)
    AUC.Lu.ori         =  M.G$R    %>% select (AUCCLu)
    AUC.F.new          =  Mnew.G$R %>% select (AUCCF)
    AUC.F.ori          =  M.G$R    %>% select (AUCCF)
    
    delta.AUC.Plas     =  AUC.Plas.new - AUC.Plas.ori
    delta.AUC.L        =  AUC.L.new - AUC.L.ori
    delta.AUC.K        =  AUC.K.new -  AUC.K.ori
    delta.AUC.Lu       =  AUC.Lu.new - AUC.Lu.ori
    delta.AUC.F        =  AUC.F.new - AUC.F.ori
    
    NSC_R     [i, 1]   <- as.numeric((delta.AUC.Plas/AUC.Plas.ori) * delta.pars)
    NSC_R     [i, 2]   <- as.numeric((delta.AUC.L /AUC.L.ori) * delta.pars)
    NSC_R     [i, 3]   <- as.numeric((delta.AUC.K  /AUC.K.ori) * delta.pars)
    NSC_R     [i, 4]   <- as.numeric((delta.AUC.Lu  /AUC.Lu.ori) * delta.pars)
    NSC_R     [i, 5]   <- as.numeric((delta.AUC.F /AUC.F.ori) * delta.pars)
  }
  return (NSC_R)
}

# Model results
A <- NSC_func (Rat.theta, pred.NSC)

rownames (A)  = names(Rat.theta)
colnames (A)  = c("NSC_CPlas", "NSC_CL", "NSC_CK", "NSC_CLu", "NSC_CF")
NSC <- data.frame(A)

NSC_Clean <- NSC %>%
  mutate_all(~ifelse(abs(.) < 1e-5, "<1E-5", 
                     ifelse(abs(.) < 0.1, formatC(., format = "e", digits = 2), 
                            round(., 2))))
write.csv(NSC_Clean, file = "NSC.csv", row.names = FALSE)

##################################### Circle barplot function ###############################################
## plot modifed from "R graph gallery: https://www.r-graph-gallery.com/297-circular-barplot-with-groups/ "  #
#############################################################################################################
melt.Plas        = melt(NSC[,1]) 
melt.Plas$group  = c("Plasma") 
melt.Liver       = melt(NSC[,2])
melt.Liver$group = c("Liver")
melt.Kidney      = melt(NSC[,3])
melt.Kidney$group = c("Kidney")
melt.Lung       = melt(NSC[,4])
melt.Lung$group = c("Lung")
melt.Fat    = melt(NSC[,5])
melt.Fat$group = c("Fat")

melt.data         = rbind (melt.Plas,melt.Liver,melt.Kidney,melt.Lung,melt.Fat)
melt.data$par     = rep(rownames(NSC),5) 
data.NSC     = melt.data%>%filter(abs(value)>=0.2)
data.NSC$group <- factor(data.NSC$group)

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 3
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data.NSC$group), ncol(data.NSC)))
colnames(to_add) <- colnames(data.NSC)
to_add$group <- rep(levels(data.NSC$group), each=empty_bar)
data.NSC  <- rbind(data.NSC, to_add)
data.NSC  <- data.NSC %>% arrange(group)
data.NSC$id <- seq(1, nrow(data.NSC))

# Get the name and the y position of each label
label_data <- data.NSC 
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse(angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

#prepare a data frame for base lines
base_data <- data.NSC%>% 
  group_by(group)%>% 
  dplyr::summarize(start=min(id)-0.2, end=max(id) - empty_bar+0.2) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

windowsFonts(Times=windowsFont("Times New Roman"))

# Make the plot
p1 <- ggplot(data.NSC , aes(x=as.factor(id), y=value, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  
  # Add a val=90/60/30 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 80, xend = start, yend = 80), colour = "grey", alpha=1, size=0.5 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 50, xend = start, yend = 50), colour = "grey", alpha=1, size=0.5 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "grey", alpha=1, size=0.5 , inherit.aes = FALSE ) +

  annotate("text", x = rep(max(data.NSC$id),3), y = c(20, 50, 80), label = c("20%", "50%", "80%")  , color="red", size=4, angle=0, fontface="bold", hjust=1) +
  geom_bar(aes(x=as.factor(id), y=abs(value*100), fill=group), stat="identity", alpha=0.5) +
  ylim(-120,120) +
  theme_minimal() +
  theme(
    legend.position = "none",
    text= element_text (family = "Times"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id,y=abs(value*100)+10, label=par, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE) +
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE) 
