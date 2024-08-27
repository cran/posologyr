## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
set.seed(1)
library(rxode2)
setRxThreads(2L)  # limit the number of threads

## ----setup--------------------------------------------------------------------
library(posologyr)

## ----model--------------------------------------------------------------------
mod_amikacin_Burdet2015 <- function() {
    ini({
      THETA_Cl=4.3
      THETA_Vc=15.9
      THETA_Vp=21.4
      THETA_Q=12.1
      ETA_Cl + ETA_Vc + ETA_Vp + ETA_Q ~
        c(0.1,
          0.01     ,   0.05 ,
          0.01     ,   0.02 ,   0.2  ,
          -0.06    ,   0.004,   0.003,    0.08)
      add_sd <- 0.2
      prop_sd <- 0.1
    })
    model({
      TVCl  = THETA_Cl*(CLCREAT4H/82)^0.7
      TVVc  = THETA_Vc*(TBW/78)^0.9*(PoverF/169)^0.4
      TVVp  = THETA_Vp
      TVQ   = THETA_Q
      Cl    = TVCl*exp(ETA_Cl)
      Vc    = TVVc*exp(ETA_Vc)
      Vp    = TVVp*exp(ETA_Vp)
      Q     = TVQ *exp(ETA_Q)
      ke    = Cl/Vc
      k12   = Q/Vc
      k21   = Q/Vp
      Cp    = centr/Vc
      d/dt(centr)  = - ke*centr - k12*centr + k21*periph
      d/dt(periph) =            + k12*centr - k21*periph

      Cp ~ add(add_sd) + prop(prop_sd) + combined1()
    })
  }

## ----tdm_patientA-------------------------------------------------------------
df_patientA <- data.frame(ID=1,TIME=c(0,1,6),
                              DV=c(NA,58,14),
                              EVID=c(1,0,0),
                              AMT=c(2000,0,0),
                              DUR=c(0.5,NA,NA),
                              CLCREAT4H=50,TBW=62,PoverF=169)
df_patientA

## ----estim_map----------------------------------------------------------------
patA_map <- poso_estim_map(dat=df_patientA,
                           prior_model=mod_amikacin_Burdet2015)

## ----map_plot_tdm-------------------------------------------------------------
plot(patA_map$model,Cc)

## ----map_Cmin_priordose_patientA----------------------------------------------
poso_time_cmin(dat=df_patientA,
               prior_model=mod_amikacin_Burdet2015,
               tdm = TRUE,
               target_cmin = 2.5)

## ----map_Cmax_optim_patientA--------------------------------------------------
map_dose <- poso_dose_conc(dat=df_patientA,
                           prior_model=mod_amikacin_Burdet2015,
                           tdm=TRUE,
                           time_c = 35,               #target concentration at t = 35 h
                           time_dose = 34,            #dosing at t = 34 h
                           duration = 0.5,
                           target_conc = 80)
map_dose

## ----map_Cmin_optim_patientA--------------------------------------------------
map_interval <- poso_inter_cmin(dat=df_patientA,
                                prior_model=mod_amikacin_Burdet2015,
                                dose = map_dose$dose,
                                duration = 0.5,
                                target_cmin = 2.5)
map_interval

