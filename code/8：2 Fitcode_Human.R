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
source (file = "8：2 Hmod.R")

## Build mrgsolve-based PBPK Model
mod <- mcode ("HumanPBPK.code", HumanPBPK.code)

theta.int <- log(c(
  KbileC                         = 0.000263,                    ## Biliary elimiatnion rate
  KurineC                        = 4.54e-6,                     ## Urinary elimination rate
  Free                           = 0.009,                       ## Free fraction in plasma
  PL                             = 8.78,                        ## Liver/plasma partition coefficient (PC)
  PK                             = 0.72,                        ## Kidney/plasma PC
  PLu                            = 1.49,                        ## Lung/plasma PC
  PF                             = 0.064,                       ## Fat/plasma PC
  PRest                          = 1.90,                        ## Rest of body/plasma PC
  K0C                            = 0.14,                        ## Rate of absorption of 8:2 Cl-PFESA  in the stomach
  Kabsc                          = 0.010,                       ## Rate of absorption of 8:2 Cl-PFESA  in the small intestines
  KunabsC                        = 7.93e-5,                     ## Rate of unabobded dose to appear in feces
  Tm                             = 2032,                        ## Transport maximum
  Kt                             = 30                           ## Transport affinity constant                 
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
  
  tsamp <- tgrid(start = 0, end = tinterval * TDoses, delta = tinterval)
  
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
# Init value: Chinese population;  Dose: 5.70 pg/kg/week (0.815 pg/kg/day) from the Sixth China Total Diet Study;
DOSE_A <- 0.815e-9

out <- Pred.child(theta.int,DOSE_A)

Cchild <- out %>%  select(Time, Plas) %>% mutate(Plas = Plas * 1000) 

Init <- out %>% filter (row_number()== n()) %>% select(-c("Plas"))

## Exposure sceinario 
pred.adult <- function(pars, DOSE,Init) {
  
  ## Get out of log domain
  pars <- lapply(pars, exp)## Return a list of exp (parametrs for gestational model) from log scale
  
  ## Exposure scenario for gestational exposure
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
    init(AVPlas_free = Init$AVPlas_free, AAPlas_free = Init$AAPlas_free, AFil = Init$AFil, AKb = Init$AKb,
         ARest = Init$ARest, AL = Init$AL, ALu = Init$ALu, AF = Init$AF, Aurine = Init$Aurine, AST = Init$AST, ASI= Init$ASI, 
         Afeces= Init$Afeces,ADOSE = Init$AT+GDOSEoral) %>% 
    param (pars) %>%
    update(atol = 1E-3, maxsteps = 500000) %>%          
    mrgsim_d (data = Gex.oral, tgrid = Gtsamp)
  
  Goutdf = cbind.data.frame(Time   = Gout$time/(24 * 365), 
                            Plas  = Gout$Plas*1000)
  
  Goutdf <- Goutdf#
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
   # Chemical-specific parameters (final mean values)
  KbileC                         = 2.63e-4,   # fitting parameter                    
  KurineC                        = 4.54e-6,   # fitting parameter                     
  Free                           = 0.009,     # fitting parameter                     
  PL                             = 8.78,                       
  PK                             = 0.72,      # fitting parameter                  
  PLu                            = 1.49,      # fitting parameter
  PF                             = 0.064,     # fitting parameter
  PRest                          = 1.90,      # fitting parameter                   
  K0C                            = 0.14,      # fitting parameter                    
  Kabsc                          = 0.010,     # fitting parameter                   
  KunabsC                        = 7.93e-5,   # fitting parameter
  Tm                             = 2032,      
  Kt                             = 30))

pred.Human <- function(pars, DOSE) {
  
  ## Get out of log domain
  pars <- lapply(pars, exp)## Return a list of exp (parametrs for gestational model) from log scale
  
  ## Exposure scenario for gestational exposure
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
  
  Goutdf <- Goutdf %>% filter (Time == 30) 
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
B <- NSC_func (Human.theta, pred.Human, DOSE_A)

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
  TDoses    <- 365 * 33    # Modelled exposure duration of 50 years
  
  ## Generate individual-specific parameters
  idata <- tibble(ID = 1:N) %>% 
    mutate(
      BW        = rnormTrunc(N, min = 25.96, max = 100.04, mean = 63, sd = 18.9),
      QluC      = rnormTrunc(N, min = 0.41, max = 1.59, mean = 1, sd = 0.3),
      Free      = rlnormTrunc(N, min = 0.0048, max = 0.015, meanlog = -4.75, sdlog = 0.29),
      PRest     = rlnormTrunc(N, min = 1.26, max = 2.75, meanlog = 0.62, sdlog = 0.20),
      KabsC     = 0.010,     # Fitting parameter   
      KunabsC   = 7.93e-5,   # Fitting parameter 
      KurineC   = 4.54e-6,   # Fitting parameter
      PK        = 0.72,      # Fitting parameter                  
      PLu       = 1.49,      # Fitting parameter
      PF        = 0.064,     # Fitting parameter
      K0C       = 0.14,      # Fitting parameter                    
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
    init(AVPlas_free = Init$AVPlas_free, AAPlas_free = Init$AAPlas_free, AFil = Init$AFil, AKb = Init$AKb,
         ARest = Init$ARest, AL = Init$AL, ALu = Init$ALu, AF = Init$AF, Aurine = Init$Aurine, AST = Init$AST, ASI= Init$ASI, 
         Afeces= Init$Afeces)  %>%
    param(pars) %>%
    data_set(Gex.oral_1) %>%
    idata_set(idata) %>%
    update(atol = 1E-8, maxsteps = 500000) %>%
    mrgsim(obsonly = TRUE, tgrid = tsamp)
  
  ## Format output
  out <- cbind.data.frame(
    ID    = out$ID,
    Time  = out$time / (24 * 365),  # Convert time to years
    CVP   = out$VP * 1000000)          # Convert concentration to pg/mL
  
  out <- out %>% filter(Time == 33) 
  
}

N = 1000
PlotDat_A1   <- PBPK_H_pop (theta.int,DOSE_A,N = N)
H_CPlas <- PlotDat_A1  %>% summarize (median = quantile (CVP, probs = 0.5), 
                                      ci_05 = quantile (CVP, probs = 0.025),
                                      ci_95 = quantile (CVP, probs = 0.975))

# Read xlsx
Human <- read.csv(file = "Human Evaluation.csv")
Human$Category <- factor(Human$Category, 
                         levels = c("observed", "50Simulate", "60Simulate"))

ggplot(Human, aes(x = Category, y = Median, color = Category, fill = Category)) + 
  geom_boxplot(width = 0.5, alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(
    data = subset(Human, Category == "observed"),
    width = 0.15, size = 2) +
  scale_color_manual(values = c("#F9A363","#BF5960", "#6F99AD" )) +
  scale_fill_manual(values = c("#F9A363","#BF5960", "#6F99AD")) +
  scale_y_break( breaks = c(0.1, 0.2), 
    scales = c(0.3, 0.7)) +  
  scale_y_continuous(
    breaks = c(seq(0, 0.1, 0.01), seq(0.2, 0.5, 0.1)), 
    labels = c(
      sprintf("%.2f", seq(0, 0.1, 0.01)),  
      sprintf("%.1f", seq(0.2, 0.5, 0.1))),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_x_discrete(
    limits = c("observed", "50Simulate", "60Simulate"),
    labels = c("Observed", "50% Simulate", "60% Simulate")
  ) +
  theme_test() +
  labs(y = "8:2 Cl-PFESA Concerntration in adult plasma (ng/mL)", x = "") +
  theme(
    axis.text = element_text(color = "black", size = 10),
    axis.line = element_blank()  )

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
      Free      = rlnormTrunc(N, min = 0.0048, max = 0.015, meanlog = -4.75, sdlog = 0.29),
      PRest     = rlnormTrunc(N, min = 1.26, max = 2.75, meanlog = 0.62, sdlog = 0.20),
      KabsC     = 0.010,     # Fitting parameter   
      KunabsC   = 7.93e-5,   # Fitting parameter 
      KurineC   = 4.54e-6,   # Fitting parameter
      PK        = 0.72,      # Fitting parameter                  
      PLu       = 1.49,      # Fitting parameter
      PF        = 0.064,     # Fitting parameter
      K0C       = 0.14,      # Fitting parameter                    
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
    init(AVPlas_free = Init$AVPlas_free, AAPlas_free = Init$AAPlas_free, AFil = Init$AFil, AKb = Init$AKb,
         ARest = Init$ARest, AL = Init$AL, ALu = Init$ALu, AF = Init$AF, Aurine = Init$Aurine, AST = Init$AST, ASI= Init$ASI, 
         Afeces= Init$Afeces)  %>%
    param(pars) %>%
    data_set(Gex.oral_1) %>%
    idata_set(idata) %>%
    update(atol = 1E-8, maxsteps = 500000) %>%
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

##Calculate HED values from population epidemiology data
DOSE_B =   1 #mg/kg/day, to calculate stable concentration

Css     <- HED (theta.int,DOSE_B,N = N) %>% select(c("ID","Time","CVP"))

HBMDL_values <- c(0.086, 0.078, 0.104, 0.026, 0.086)
HBMD_names <- c("AST", "ALB","TC", "LDL", "HDL")

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
