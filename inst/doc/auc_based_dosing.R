## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  collapse=TRUE,
  comment="#>"
)
set.seed(1)
library(rxode2)
setRxThreads(2L)  # limit the number of threads

## ----setup--------------------------------------------------------------------
library(posologyr)

## ----model--------------------------------------------------------------------
mod_vancomycin_Goti2018 <- list(
  ppk_model   = rxode2::rxode({
    centr(0) = 0;
    TVCl  = THETA_Cl*(CLCREAT/120)^0.8*(0.7^DIAL);
    TVVc  = THETA_Vc*(WT/70)          *(0.5^DIAL);
    TVVp  = THETA_Vp;
    TVQ   = THETA_Q;
    Cl    = TVCl*exp(ETA_Cl);
    Vc    = TVVc*exp(ETA_Vc);
    Vp    = TVVp*exp(ETA_Vp);
    Q     = TVQ;
    ke    = Cl/Vc;
    k12   = Q/Vc;
    k21   = Q/Vp;
    Cc    = centr/Vc;
    d/dt(centr)  = - ke*centr - k12*centr + k21*periph;
    d/dt(periph) =            + k12*centr - k21*periph;
    d/dt(AUC)    =   Cc;
  }),
  error_model = function(f,sigma){
    g <- sigma[1] + sigma[2]*f
    return(g)
  },
  theta = c(THETA_Cl=4.5, THETA_Vc=58.4, THETA_Vp=38.4,THETA_Q=6.5),
  omega = lotri::lotri({ETA_Cl + ETA_Vc + ETA_Vp + ETA_Q ~
      c(0.147,
        0     ,   0.510,
        0     ,       0,   0.282,
        0     ,       0,       0,    0)}),
  covariates  = c("CLCREAT","WT","DIAL"),
  sigma       = c(additive_a = 3.4, proportional_b = 0.227))

## ----tdm_patientB-------------------------------------------------------------
df_patientB <- data.frame(ID=1,TIME=c(0.0,13.0,24.2,48),
                          DV=c(NA,12,NA,9.5),
                          AMT=c(2000,0,1000,0),
                          DUR=c(2,NA,2,NA),
                          EVID=c(1,0,1,0),
                          CLCREAT=65,WT=70,DIAL=0)
df_patientB

## ----estim_map----------------------------------------------------------------
patB_map <- poso_estim_map(dat=df_patientB,
                           prior_model=mod_vancomycin_Goti2018)

## -----------------------------------------------------------------------------
plot(patB_map$model,Cc)

## ----map_plot_tdm-------------------------------------------------------------
#Get the observations from the patient record
indiv_obs             <- df_patientB[,c("DV","TIME")]
names(indiv_obs)      <- c("value","time")

#Overlay the MAP profile and the observations
plot(patB_map$model,Cc) +
  ggplot2::ylab("Central concentration") +
  ggplot2::geom_point(data=indiv_obs, size= 3, na.rm=TRUE)


## ----AUC_map_dose-------------------------------------------------------------
#AUC 0_24
AUC_map_first_dose <- patB_map$model$AUC[which(patB_map$model$time == 24)]
AUC_map_first_dose

#AUC 24_48
AUC_map_second_dose <- patB_map$model$AUC[which(patB_map$model$time == 48)] - AUC_map_first_dose
AUC_map_second_dose

## ----optim_next_dose----------------------------------------------------------
poso_dose_auc(dat=df_patientB,
              prior_model=mod_vancomycin_Goti2018,
              tdm=TRUE,
              time_auc=24,            #AUC24
              time_dose = 48,         #48 h: immediately following the last observation
              duration=2,             #infused over 2 h
              target_auc=400)

## ----optim_maintenance_dose---------------------------------------------------
poso_dose_auc(dat=df_patientB,
              prior_model=mod_vancomycin_Goti2018,
              time_auc=24,
              starting_time=24*9,
              interdose_interval=24,
              add_dose=10,
              duration=2,
              target_auc=400)

## ----continuous_infusion------------------------------------------------------
poso_dose_auc(dat=df_patientB,
              prior_model=mod_vancomycin_Goti2018,
              time_auc=24,
              starting_time=24*9,
              interdose_interval=24,
              add_dose=10,
              duration=24,
              target_auc=400)

