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
library(ggbreak)     # for plotting the figure

## Input mrgsolve-based PBPK Model
source (file = "6：2 Hmod.R")

## Build mrgsolve-based PBPK Model
mod <- mcode ("HumanPBPK.code", HumanPBPK.code)

theta.int <- log(c(
  Vmax_baso_invitro              = 479,                      
  Km_baso                        = 64.4,                        
  Vmax_apical_invitro            = 51803,                        
  Km_apical                      = 20.1,                        
  RAFapi                         = 0.001,                       
  RAFbaso                        = 1,                        
  KeffluxC                       = 0.74,                        
  KbileC                         = 3.55e-4,                       
  KurineC                        = 5.68e-6,                        
  Free                           = 4.33e-3,                        
  PL                             = 3.73,                        
  PK                             = 0.75,                         
  PLu                            = 0.63,
  PF                             = 0.028,
  PRest                          = 0.64,                         
  K0C                            = 0.13,                           
  Kabsc                          = 0.015,                        
  Kdif                           = 0.001,                       
  KunabsC                        = 1.49e-3                     
))

saveRDS(theta.int, file = "int_H.rds") 

Pred.child <- function(pars,DOSE_A) {
  
  ## Get out of log domain
  pars <- exp(pars)
  
  tinterval <- 24  
  TDoses    <- 365 * 17 
  ex <- tibble(
    ID   = rep(1, TDoses - 365 + 1), 
    time = seq(from = 24 * 365, to = tinterval * TDoses, by = tinterval)  
  ) %>%
    mutate(
      DAY   = time / 24,  
      YEAR  = DAY / 365, 
      BW    = if_else(YEAR <= 17, 
                      true  = (9.86 + 0.370*YEAR)/(1 - 0.0789*YEAR + 0.00205*YEAR^2),
                      ifelse(YEAR > 17, 63, 63)),           ##Equation from Chang 2022
      amt   = DOSE_A * BW,  
      cmt   = "AST",        
      ii    = tinterval,    
      evid  = 1)
  
  tsamp <- tgrid(start = 24 * 365, end = tinterval * TDoses, delta = tinterval)
  
  out <- mod %>%
    param (pars) %>%
    update(atol = 1E-8, maxsteps = 5000) %>%
    mrgsim_d(data = ex, tgrid = tsamp)%>%
    filter (time > 0) %>%
    mutate(Time = time / (24 * 365)) %>%  
    select(Time, everything())  
  
  out <- out %>% distinct(time, .keep_all = TRUE)
  
  return(out)
}

## The assumed exposure dose were estiamted from previous literatures
# Init value: Chinese population ;Dose: 2.75 ng/kg/week (0.393 ng/kg/day) from the Sixth China Total Diet Study;
DOSE_A <- 0.393e-6 

out <- Pred.child(theta.int,DOSE_A)

Cchild <- out %>%  select(Time, Plas) %>% mutate(Plas = Plas * 1000) 

Init <- Pred.child(theta.int,DOSE_A) %>% filter (row_number()== n()) %>% select(-c("Plas"))

## Exposure sceinario 
pred.adult <- function(pars, DOSE,Init) {
  
  ## Get out of log domain
  pars <- lapply(pars, exp)## Return a list of exp parametrs from log scale
  
  ## Exposure scenario 
  GBW          = 63                 ## Body weight during gestation 
  tinterval    = 24                 ## Time interval; 
  GTDOSE       = 365 * 33           ## Total dosing/Dose times; 
  GDOSE        = DOSE               ## Input oral dose  
  GDOSEoral    = GDOSE*GBW          ## Amount of oral dose
  
  # To create exposure scenario
  Gex.oral <- ev (ID   = 1, 
                  time = 0,             ## Dosing start time
                  amt  = GDOSEoral,     ## Amount of dose 
                  ii   = tinterval,     ## Time interval
                  addl = GTDOSE - 1,    ## Additional doseing 
                  cmt  = "AST",         ## dosing: AST Stomach  
                  replicate = FALSE)    ## No replicate
  
  Gtsamp  = tgrid(0, tinterval*(GTDOSE - 1) + 24*1, 24) ## Simulation time 
  
  ## Simulation of exposure scenario
  Gout <- 
    mod %>%
    init(AVPlas_free = Init$AVPlas_free, AAPlas_free = Init$AAPlas_free, APTC = Init$APTC, AFil = Init$AFil, AKb = Init$AKb,
         ARest = Init$ARest,AL = Init$AL, ALu = Init$ALu, AF = Init$AF, A_baso = Init$A_baso, A_apical = Init$A_apical, 
         Adif = Init$Adif,Aefflux = Init$Aefflux, ACI = Init$ACI, Aurine = Init$Aurine, AST = Init$AST, ASI= Init$ASI, 
         Afeces= Init$Afeces, Ameta= Init$Ameta, ADOSE = Init$AT+GDOSEoral) %>% 
    param (pars) %>%
    update(atol = 1E-3, maxsteps = 500000) %>%          
    mrgsim_d (data = Gex.oral, tgrid = Gtsamp)
  
  Goutdf = cbind.data.frame(Time   = Gout$time/(24 * 365), 
                            Plas  = Gout$Plas*1000)
  
  return (Goutdf)
}

## Estimate the model residual with experimental data by modCost function (from FME package)
outdf.A <- pred.adult (pars = theta.int, DOSE = DOSE_A, Init = Init)

Cadult <- outdf.A %>%  mutate(Time = Time + 17) 

Plas <- rbind(Cchild,Cadult)

## NSC
Human.theta <- log(c(
  # Physiological parameters
  BW                             = 63,
  Htc                            = 0.44,
  QCC                            = 15,
  QLC                            = 0.255, 
  QKC                            = 0.19,
  QLuC                           = 1,
  QFC                            = 0.05,
  VLC                            = 0.022,
  VKC                            = 0.0046,
  VFC                            = 0.143,
  VLuC                           = 0.02,
  VPlasC                         = 0.0312,
  VFilC                          = 4.6e-4,  
  VPTCC                          = 1.35e-4, 
  # Chemical-specific parameters (final mean values)
  Vmax_baso_invitro              = 479,                      
  Km_baso                        = 64.4,                        
  Vmax_apical_invitro            = 51803,                        
  Km_apical                      = 20.1,                        
  RAFapi                         = 0.001,                       
  RAFbaso                        = 1,                        
  KeffluxC                       = 0.74,                        
  KbileC                         = 3.55e-4,   # fitting parameter                    
  KurineC                        = 5.68e-6,   # fitting parameter                     
  Free                           = 4.33e-3,   # fitting parameter                     
  PL                             = 3.73,      # fitting parameter                   
  PK                             = 0.75,      # fitting parameter                  
  PLu                            = 0.63,      # fitting parameter
  PF                             = 0.028,     # fitting parameter
  PRest                          = 0.64,      # fitting parameter                   
  K0C                            = 0.13,      # fitting parameter                    
  Kabsc                          = 0.015,     # fitting parameter                   
  Kdif                           = 0.001,                       
  KunabsC                        = 1.49e-3    # fitting parameter
))

pred.Human <- function(pars, DOSE) {
  
  ## Get out of log domain
  pars <- lapply(pars, exp)## Return a list of exp parametrs from log scale
  
  ## Exposure scenario
  GBW          = 63                 ## Body weight during gestation 
  tinterval    = 24                 ## Time interval; 
  GTDOSE       = 365 * 30           ## Total dosing/Dose times; 
  GDOSE        = DOSE               ## Input oral dose  
  GDOSEoral    = GDOSE*GBW          ## Amount of oral dose
  
  # To create exposure scenario
  Gex.oral <- ev (ID   = rep(1, 24*365*30+1), 
                  time = 0,             ## Dosing start time
                  amt  = GDOSEoral,     ## Amount of dose 
                  ii   = tinterval,     ## Time interval
                  addl = GTDOSE - 1,    ## Additional doseing 
                  cmt  = "AST",         ## dosing: AST Stomach  
                  replicate = FALSE)    ## No replicate
  
  Gtsamp  = tgrid(0, tinterval*(GTDOSE - 1) + 24*1, 1) ## Simulation time 
  
  ## Simulation of exposure scenario
  Gout <- 
    mod %>%
    param (pars) %>%
    update(atol = 1E-3, maxsteps = 500000) %>%          
    mrgsim_d (data = Gex.oral, tgrid = Gtsamp)
  
  Goutdf = cbind.data.frame(Time   = Gout$time/(24 * 365), 
                            CPlas = Gout$Plas, 
                            AUCCV = Gout$AUC_CV,
                            AUCCL = Gout$AUC_CL,
                            AUCCK = Gout$AUC_CK,
                            AUCCLu= Gout$AUC_CLu,
                            AUCCF = Gout$AUC_CF)
  
  Goutdf <- Goutdf %>% filter (Time == 30) ## common sampling time in biomonitoring studies
  return (list("G" = Goutdf))
  
}

NSC_func <- function (pars, Pred, DOSE) {
  nG <- length(pars)
  NSC_H    = matrix(NA, nrow = length(pars) ,ncol = 5)
  
  for (i in 1:nG) {
    pars.new      <- pars %>% replace(i, log(exp((pars[i]))*1.01))
    Rnew.G         <- Pred(pars.new, DOSE)
    R.G            <- Pred(pars, DOSE)
    delta.pars    <- exp(pars[i])/(exp(pars[i])*0.01)
    
    ## Estimated the AUC
    AUC.Plas.new       =  Rnew.G$G %>% select (AUCCV)
    AUC.Plas.ori       =  R.G$G    %>% select (AUCCV)
    AUC.L.new          =  Rnew.G$G %>% select (AUCCL)
    AUC.L.ori          =  R.G$G    %>% select (AUCCL)
    AUC.K.new          =  Rnew.G$G    %>% select (AUCCK)
    AUC.K.ori          =  R.G$G    %>% select (AUCCK)
    AUC.Lu.new         =  Rnew.G$G    %>% select (AUCCLu)
    AUC.Lu.ori         =  R.G$G    %>% select (AUCCLu)
    AUC.F.new          =  Rnew.G$G    %>% select (AUCCF)
    AUC.F.ori          =  R.G$G    %>% select (AUCCF)
    
    delta.AUC.Plas     =  AUC.Plas.new - AUC.Plas.ori
    delta.AUC.L        =  AUC.L.new - AUC.L.ori
    delta.AUC.K        =  AUC.K.new -  AUC.K.ori
    delta.AUC.Lu       =  AUC.Lu.new - AUC.Lu.ori
    delta.AUC.F        =  AUC.F.new - AUC.F.ori
    
    NSC_H     [i, 1]   <- as.numeric((delta.AUC.Plas/AUC.Plas.ori) * delta.pars)
    NSC_H     [i, 2]   <- as.numeric((delta.AUC.L /AUC.L.ori) * delta.pars)
    NSC_H     [i, 3]   <- as.numeric((delta.AUC.K  /AUC.K.ori) * delta.pars)
    NSC_H     [i, 4]   <- as.numeric((delta.AUC.Lu  /AUC.Lu.ori) * delta.pars)
    NSC_H     [i, 5]   <- as.numeric((delta.AUC.F /AUC.F.ori) * delta.pars)
  }
  return (NSC_H)
}

## Gather all the data for plotting
B <- NSC_func (Human.theta, pred.Human, 0.393e-6)

rownames (B)  = names(Human.theta)
colnames (B)  = c("NSC_CPlas", "NSC_CL", "NSC_CK", "NSC_CLu", "NSC_CF")

NSC_H <- data.frame(B)

NSC_Clean <- NSC_H %>%
  mutate_all(~ifelse(. == 0, "<1e-5", 
                     ifelse(. > 0.1, round(., 2), 
                            formatC(., format = "e", digits = 2))))
write.csv(NSC_Clean, file = "NSC_H.csv", row.names = FALSE)

##################################### Circle barplot function ###############################################
## plot modifed from "R graph gallery: https://www.r-graph-gallery.com/297-circular-barplot-with-groups/ "  #
#############################################################################################################
melt.Plas        = melt(NSC_H[,1]) 
melt.Plas$group  = c("Plasma") 
melt.Liver       = melt(NSC_H[,2])
melt.Liver$group = c("Liver")
melt.Kidney      = melt(NSC_H[,3])
melt.Kidney$group = c("Kidney")
melt.Lung       = melt(NSC_H[,4])
melt.Lung$group = c("Lung")
melt.Fat    = melt(NSC_H[,5])
melt.Fat$group = c("Fat")

melt.data         = rbind (melt.Plas,melt.Liver,melt.Kidney,melt.Lung,melt.Fat)
melt.data$par     = rep(rownames(NSC_H),5) 
data.NSC_H     = melt.data%>%filter(abs(value)>=0.2)
data.NSC_H$group <- factor(data.NSC_H$group)

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 3
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data.NSC_H$group), ncol(data.NSC_H)))
colnames(to_add) <- colnames(data.NSC_H)
to_add$group <- rep(levels(data.NSC_H$group), each=empty_bar)
data.NSC_H  <- rbind(data.NSC_H, to_add)
data.NSC_H  <- data.NSC_H %>% arrange(group)
data.NSC_H$id <- seq(1, nrow(data.NSC_H))

# Get the name and the y position of each label
label_data <- data.NSC_H 
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse(angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

#prepare a data frame for base lines
base_data <- data.NSC_H%>% 
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
p1 <- ggplot(data.NSC_H , aes(x=as.factor(id), y=value, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  
  # Add a val=90/60/30 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 80, xend = start, yend = 80), colour = "grey", alpha=1, size=0.5 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 50, xend = start, yend = 50), colour = "grey", alpha=1, size=0.5 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "grey", alpha=1, size=0.5 , inherit.aes = FALSE ) +
 
  annotate("text", x = rep(max(data.NSC_H$id),3), y = c(20, 50, 80), label = c("20%", "50%", "80%")  , color="red", size=4, angle=0, fontface="bold", hjust=1) +
  
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

##############################  MCMC  ##############################################
PBPK_H_pop <- function(pars, DOSE_A, N) {
  
  Init <- Pred.child (theta.int,DOSE_A) %>% filter (row_number()== n()) %>% select(-c("Plas"))
  
  ## Get out of log domain
  pars <- exp(pars)
  
  ## Time parameters
  tinterval <- 24  
  TDoses    <- 365 * 33     # Modelled exposure duration of 50 years
  
  ## Generate individual-specific parameters
  idata <- tibble(ID = 1:N) %>% 
    mutate(
      BW        = rnormTrunc(N, min = 25.96, max = 100.04, mean = 63, sd = 18.9),
      QluC      = rnormTrunc(N, min = 0.41, max = 1.59, mean = 1, sd = 0.3),
      VLC       = rnormTrunc(N, min = 0.009, max = 0.034, mean = 0.022, sd = 6.60e-3),
      PL        = rlnormTrunc(N, min = 2.48, max = 5.39, meanlog = 1.3, sdlog = 0.20),
      Free      = rlnormTrunc(N, min = 0.0023, max = 0.0074, meanlog = -5.49, sdlog = 0.29),
      PRest     = rlnormTrunc(N, min = 0.43, max = 0.93, meanlog = -0.47, sdlog = 0.20),
      KunabsC   = rlnormTrunc(N, min = 0.0008, max = 0.00254, meanlog = -6.55, sdlog = 0.29),
      KabsC     = rlnormTrunc(N, min = 0.0081, max = 0.0255, meanlog = -4.24, sdlog = 0.29),
      KbileC    = 3.55e-4,   # Fitting parameter                    
      PK        = 0.75,      # Fitting parameter                  
      PLu       = 0.63,      # Fitting parameter
      PF        = 0.028,     # Fitting parameter
      K0C       = 0.13,      # Fitting parameter                    
      DOSEoral.A = DOSE_A * BW
    )
  
  ## Define time grid for simulation
  tsamp <- tgrid(start = 24 * 365, end = tinterval * TDoses, delta = 24 * 365)
  
  ## Define dosing events using idata
  Gex.oral_1 <- ev(
    ID   = 1:N,             # One individual
    time = 0,               # Dose start time
    amt  = idata$DOSEoral.A, # Individual-specific dose
    ii   = tinterval,       # Time interval
    addl = TDoses - 1,      # Additional dosing
    cmt  = "AST",           # The dosing compartment: AST Stomach  
    replicate = FALSE       # No replicate
  )
  
  ## Run the simulation
  out <- mod %>%
    init(AVPlas_free = Init$AVPlas_free, AAPlas_free = Init$AAPlas_free, APTC = Init$APTC, AFil = Init$AFil, AKb = Init$AKb,
         ARest = Init$ARest,AL = Init$AL, ALu = Init$ALu, AF = Init$AF, A_baso = Init$A_baso, A_apical = Init$A_apical, 
         Adif = Init$Adif,Aefflux = Init$Aefflux, ACI = Init$ACI, Aurine = Init$Aurine, AST = Init$AST, ASI= Init$ASI, 
         Afeces= Init$Afeces, Ameta= Init$Ameta) %>% 
    param(pars) %>%
    data_set(Gex.oral_1) %>%
    idata_set(idata) %>%
    update(atol = 1E-8, maxsteps = 5000) %>%
    mrgsim(obsonly = TRUE, tgrid = tsamp)
  
  ## Format output
  out <- cbind.data.frame(
    ID    = out$ID,
    Time  = out$time / (24 * 365),  # Convert time to years
    CVP   = out$VP * 1000)          # Convert concentration to ng/mL
  
  out <- out %>% filter(Time == 33) 
  
}

N = 1000
PlotDat_A1   <- PBPK_H_pop (theta.int,DOSE_A,N = N)
H_CPlas <- PlotDat_A1  %>% summarize (median = quantile (CVP, probs = 0.5), 
                                      ci_05 = quantile (CVP, probs = 0.025),
                                      ci_95 = quantile (CVP, probs = 0.975))

# box plot
Human <- read.csv(file ="Human Evaluation.csv")
table(Human$Category)
Human$Category <- factor(Human$Category, 
                         levels = c("observed", "40Simulate", "50Simulate"))

ggplot(Human, aes(Category, Median, color = Category, fill = Category)) + 
  geom_jitter(data = subset(Human, Category == "observed"), width = 0.15, size = 2) +
  geom_boxplot(width = 0.5, alpha = 0.7, outlier.shape = NA) +  
  scale_color_manual(values = c("#BF5960", "#6F99AD","#F9A363")) + 
  scale_fill_manual(values = c("#BF5960", "#6F99AD","#F9A363")) +  
  theme_test() +
  labs(y = "6:2 Cl-PFESA Concerntration in adult plasma (ng/mL)", x = "") +
  scale_y_continuous(expand = c(0.1, 0.1)) +
  scale_x_discrete(limits = c("observed", "40Simulate", "50Simulate")) +  
  theme(axis.text = element_text(color = "black", size = 10),
        axis.line = element_blank())

##############################  HED Derivation  ##############################################
HED <- function(pars, DOSE_A, N) {
  
  Init <- Pred.child (theta.int,DOSE_A) %>% filter (row_number()== n()) %>% select(-c("Plas"))
  
  ## Get out of log domain
  pars <- exp(pars)
  
  ## Time parameters
  tinterval <- 24  
  TDoses    <- 365 * 33 
  
  ## Generate individual-specific parameters
  idata <- tibble(ID = 1:N) %>% 
    mutate(
      BW        = rnormTrunc(N, min = 25.96, max = 100.04, mean = 63, sd = 18.9),
      QluC      = rnormTrunc(N, min = 0.41, max = 1.59, mean = 1, sd = 0.3),
      VLC       = rnormTrunc(N, min = 0.009, max = 0.034, mean = 0.022, sd = 6.60e-3),
      PL        = rlnormTrunc(N, min = 2.48, max = 5.39, meanlog = 1.3, sdlog = 0.20),
      Free      = rlnormTrunc(N, min = 0.0023, max = 0.0074, meanlog = -5.49, sdlog = 0.29),
      PRest     = rlnormTrunc(N, min = 0.43, max = 0.93, meanlog = -0.47, sdlog = 0.20),
      KunabsC   = rlnormTrunc(N, min = 0.0008, max = 0.00254, meanlog = -6.55, sdlog = 0.29),
      KabsC     = rlnormTrunc(N, min = 0.0081, max = 0.0255, meanlog = -4.24, sdlog = 0.29),
      KbileC    = 3.55e-4,   # Fitting parameter                    
      PK        = 0.75,      # Fitting parameter                  
      PLu       = 0.63,      # Fitting parameter
      PF        = 0.028,     # Fitting parameter
      K0C       = 0.13,      # Fitting parameter                    
      DOSEoral.A = DOSE_A * BW
    )
  
  ## Define time grid for simulation
  tsamp <- tgrid(start = 24 * 365, end = tinterval * TDoses, delta = 24 * 365)
  
  ## Define dosing events using idata
  Gex.oral_1 <- ev(
    ID   = 1:N,             # One individual
    time = 0,               # Dose start time
    amt  = idata$DOSEoral.A, # Individual-specific dose
    ii   = tinterval,       # Time interval
    addl = TDoses - 1,      # Additional dosing
    cmt  = "AST",           # The dosing compartment: AST Stomach  
    replicate = FALSE       # No replicate
  )
  
  ## Run the simulation
  out <- mod %>%
    init(AVPlas_free = Init$AVPlas_free, AAPlas_free = Init$AAPlas_free, APTC = Init$APTC, AFil = Init$AFil, AKb = Init$AKb,
         ARest = Init$ARest,AL = Init$AL, ALu = Init$ALu, AF = Init$AF, A_baso = Init$A_baso, A_apical = Init$A_apical, 
         Adif = Init$Adif,Aefflux = Init$Aefflux, ACI = Init$ACI, Aurine = Init$Aurine, AST = Init$AST, ASI= Init$ASI, 
         Afeces= Init$Afeces, Ameta= Init$Ameta) %>% 
    param(pars) %>%
    data_set(Gex.oral_1) %>%
    idata_set(idata) %>%
    update(atol = 1E-8, maxsteps = 5000) %>%
    mrgsim(obsonly = TRUE, tgrid = tsamp)
  
  ## Format output
  out <- cbind.data.frame(
    ID    = out$ID,
    Time  = out$time / (24 * 365),  # Convert time to years
    CVP   = out$VP,
    AUCCV = out$AUC_CV,
    AUCCL = out$AUC_CL)   
  
  out <- out %>% filter(Time == 33) }

N = 1000

DOSE_B = 100 #mg/kg/day，same as rat toxic experience

HAUC    <- HED (theta.int,DOSE_B,N = N) %>% select(c("ID","Time","AUCCV","AUCCL"))

HAUC$Avg.CV = HAUC$AUCCV/(33*365*24) 
HAUC$Avg.CL = HAUC$AUCCL/(33*365*24)
RAUC <- read.csv(file ="RAUC.csv")

##Calculate HED values from rat toxicity test
BMDL_values <- c(28.51, 52.4, 41.38, 364.76, 91.55, 24.18, 53.68, 73.82)
BMD_names <- c("ALT", "AST", "ALP", "ALB", "TC", "TG", "LDL", "BUN")

HED_results <- list()

for (i in seq_along(BMDL_values)) {
  BMDL <- BMDL_values[i]
  
  HED_results[[paste0(BMD_names[i], ".Plas")]] <- ((RAUC$Avg.CV) / (HAUC$AUCCV)) * BMDL
  HED_results[[paste0(BMD_names[i], ".Liver")]] <- ((RAUC$Avg.CL) / (HAUC$AUCCV)) * BMDL
}

for (i in seq_along(BMDL_values)) {
  HAUC[[paste0("HED.Plas.", BMD_names[i])]] <- HED_results[[paste0(BMD_names[i], ".Plas")]]
  HAUC[[paste0("HED.Liver.", BMD_names[i])]] <- HED_results[[paste0(BMD_names[i], ".Liver")]]
}

HED.human.range <- do.call(cbind.data.frame, lapply(seq_along(BMDL_values), function(i) {
  data.frame(
    setNames(
      list(
        quantile(HAUC[[paste0("HED.Plas.", BMD_names[i])]], probs = 0.5, names = FALSE, na.rm = TRUE),
        quantile(HAUC[[paste0("HED.Plas.", BMD_names[i])]], probs = 0.025, names = FALSE, na.rm = TRUE),
        quantile(HAUC[[paste0("HED.Plas.", BMD_names[i])]], probs = 0.975, names = FALSE, na.rm = TRUE),
        quantile(HAUC[[paste0("HED.Liver.", BMD_names[i])]], probs = 0.5, names = FALSE, na.rm = TRUE),
        quantile(HAUC[[paste0("HED.Liver.", BMD_names[i])]], probs = 0.025, names = FALSE, na.rm = TRUE),
        quantile(HAUC[[paste0("HED.Liver.", BMD_names[i])]], probs = 0.975, names = FALSE, na.rm = TRUE)
      ),
      c(
        paste0(BMD_names[i], ".Median"),
        paste0(BMD_names[i], ".lower"),
        paste0(BMD_names[i], ".upper"),
        paste0(BMD_names[i], ".Median.CL"),
        paste0(BMD_names[i], ".lower.CL"),
        paste0(BMD_names[i], ".upper.CL"))))
  }))

##Change the BMDL(POD) unit from mg/kg/d to ng/kg/d
HED.Rat.range <- HED.Rat.range %>%
  mutate(across(where(is.numeric), ~ . * 10^6))

##Divided by UFs(3.16*10)
HBGV.Rat.range <- HED.Rat.range %>%
  mutate(across(where(is.numeric), ~ . /31.6))

##Calculate HED values from population epidemiology data
DOSE_C =   1 #mg/kg/day, to calculate stable concentration

Css     <- HED (theta.int,DOSE_C,N = N) %>% select(c("ID","Time","CVP"))

HBMDL_values <- c(2.28, 5.11, 14.82, 6.39, 2.32, 3.04, 0.63, 11.63, 2.59)
HBMD_names <- c("ALT", "AST", "ALP", "ALB","TC", "TG", "LDL", "HDL", "Uric Acid")

Css$inverse_Css <- 1 / (Css$CVP/10^3)      ##Change the CVP unit from mg/L to mg/mL 
 
# Calculate quantiles（2.5%, 50%, 97.5%）
quantiles <- quantile(Css$inverse_Css, probs = c(0.025, 0.5, 0.975))

HED.Human.range <- data.frame(Name = character(), HBMDL = numeric(), P2.5 = numeric(), P50 = numeric(), P97.5 = numeric())

for (i in seq_along(HBMDL_values)) {
  HBMDL <- HBMDL_values[i]
  HBMD_name <- HBMD_names[i]
  
  P2.5 <- quantiles[1] * HBMDL
  P50 <- quantiles[2] * HBMDL
  P97.5 <- quantiles[3] * HBMDL
  
  HED.Human.range <- rbind(HED.Human.range, data.frame(Name = HBMD_name, HBMDL = HBMDL, P2.5 = P2.5, P50 = P50, P97.5 = P97.5))
}

##Divided by UFs(10)
HBGV.Human.range <- HED.Human.range%>%
  mutate(across(c(P2.5, P50, P97.5), ~ . / 10))

write.csv(HBGV.Human.range,"HBGV.Human.range.csv") 

##HQ for Cl-PFESAs 
HQ <- read.csv("HI.csv")
HQ_UF10 <- HQ %>%
  pivot_longer(cols = c(C8HQUF10, C10HQUF10),
               names_to = "HQ_type",
               values_to = "HQ_value") %>%
  mutate(category = case_when(
      HQ_value > 1 ~ "High concern",
      HQ_value <= 1 ~ "Low concern"),
  shape_type = case_when( HQ_type == "C8HQUF10" ~ "C8 (Triangle)",
      HQ_type == "C10HQUF10" ~ "C10 (Square)"))

HQ_UF10$Province <- factor(HQ_UF10$Province, levels = HQ$Province)

ggplot(HQ_UF10, aes(x = HQ_value, y = Province, color = category, shape = shape_type)) +
  geom_rect(aes(xmin = 0.0001, xmax = 1, ymin = -Inf, ymax = Inf),
            fill = "lightblue", alpha = 0.2, inherit.aes = FALSE) +
  geom_rect(aes(xmin = 1, xmax = 100, ymin = -Inf, ymax = Inf),
            fill = "#FD9EB2", alpha = 0.2, inherit.aes = FALSE) +
  geom_point(size = 3) +
  scale_x_log10(breaks = c(0.0001, 0.01, 1, 100),
    labels = c("0.0001", "0.01", "1", "100")) +
  scale_color_manual(
    values = c("High concern" = "darkred", "Low concern" = "blue"),
    name = "HQ Category") +
  scale_shape_manual(
    values = c("C8 (Triangle)" = 17, "C10 (Square)" = 15),
    name = "HQ Type") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.line.x = element_line(color = "black", linewidth = 0.8),
    axis.line.y = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 12, color = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.spacing = unit(0.5, "cm")) +
  labs(x = "Hazard Quotient (HQ) [log scale]") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray")

HQ_UF1 <- HQ %>%
  pivot_longer(cols = c(C8HQUF1, C10HQUF1),
               names_to = "HQ_type",
               values_to = "HQ_value") %>%
  mutate(category = case_when(
    HQ_value > 1 ~ "High concern",
    HQ_value <= 1 ~ "Low concern"),
    shape_type = case_when( HQ_type == "C8HQUF1" ~ "C8 (Triangle)",
                            HQ_type == "C10HQUF1" ~ "C10 (Square)"))

HQ_UF1$Province <- factor(HQ_UF1$Province, levels = HQ$Province)

ggplot(HQ_UF1, aes(x = HQ_value, y = Province, color = category, shape = shape_type)) +
  geom_rect(aes(xmin = 0.0001, xmax = 1, ymin = -Inf, ymax = Inf),
            fill = "lightblue", alpha = 0.2, inherit.aes = FALSE) +
  geom_rect(aes(xmin = 1, xmax = 100, ymin = -Inf, ymax = Inf),
            fill = "#FD9EB2", alpha = 0.2, inherit.aes = FALSE) +
  geom_point(size = 3) +
  scale_x_log10(breaks = c(0.0001, 0.01, 1, 100),
                labels = c("0.0001", "0.01", "1", "100")) +
  scale_color_manual(
    values = c("High concern" = "darkred", "Low concern" = "blue"),
    name = "HQ Category") +
  scale_shape_manual(
    values = c("C8 (Triangle)" = 17, "C10 (Square)" = 15),
    name = "HQ Type") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.line.x = element_line(color = "black", linewidth = 0.8),
    axis.line.y = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 12, color = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.spacing = unit(0.5, "cm")) +
  labs(x = "Hazard Quotient (HQ) [log scale]") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray")


###HI for Cl-PFESAs
data <- read.csv("HI.csv") %>% 
  select(Province, C8HQUF10, C10HQUF10, C8HQUF1, C10HQUF1) %>%
  pivot_longer(
    cols = -Province,
    names_to = "Compound",
    values_to = "Value") %>%
  mutate(
    Group = case_when(
    Compound %in% c("C8HQUF10", "C10HQUF10") ~ "HIUF10",
    Compound %in% c("C8HQUF1", "C10HQUF1") ~ "HIUF1"),
    Compound = factor(Compound, 
                      levels = c("C8HQUF10", "C10HQUF10", "C8HQUF1", "C10HQUF1")))
color_palette <- c(
  "C8HQUF10" = "#2A6F97","C10HQUF10" = "#E69F00","C8HQUF1" = "#009E73", "C10HQUF1" = "#D55E00")

#UF=10
ggplot(data %>% filter(Group == "HIUF10"), 
        aes(x = Province, y = Value, fill = Compound)) +
  geom_col(position = "stack", width = 0.7) + 
  scale_fill_manual(values = color_palette) + 
  scale_y_break(breaks = c(7, 26), scales = c(0.2, 0.8)) +  
  scale_y_continuous(
    breaks = c(seq(0, 7, 0.5), seq(26, 30, 1)), 
    labels = c(seq(0, 7, 0.5), seq(26, 30, 1)),  
    expand = expansion(mult = c(0, 0.05)) ) +
  labs(
    x = "Province", 
    y = "Hazard Index (HIUF10)",
    title = "Hazard Index Composition (HIUF10) by Province"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.ticks.x = element_line(color = "black", linewidth = 1),
    axis.ticks.y = element_line(color = "black", linewidth = 0.5), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold")  )
                     
#UF=1
ggplot(data %>% filter(Group == "HIUF1"), 
       aes(x = Province, y = Value, fill = Compound)) +
  geom_col(position = "stack", width = 0.7) + 
  scale_fill_manual(values = color_palette) + 
  scale_y_break(breaks = c(0.65, 2.6), scales = c(0.2, 0.8)) +  
  scale_y_continuous(
    breaks = c(seq(0, 0.65, 0.05), seq(2.6, 3, 0.1)), 
    labels = c(seq(0, 0.65, 0.05), seq(2.6, 3, 0.1)),  
    expand = expansion(mult = c(0, 0.05)) ) +
  labs(
    x = "Province", 
    y = "Hazard Index (HIUF1)",
    title = "Hazard Index Composition (HIUF1) by Province"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.ticks.x = element_line(color = "black", linewidth = 1),
    axis.ticks.y = element_line(color = "black", linewidth = 0.5), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold")  )
