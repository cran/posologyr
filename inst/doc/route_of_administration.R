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
mod_ganciclovir_Caldes_2009 <- function() {
  ini({
    THETA_cl  <- 7.49
    THETA_v1  <- 31.90
    THETA_cld <- 10.20
    THETA_v2  <- 32.0
    THETA_ka  <- 0.895
    THETA_baf <- 0.825
    ETA_cl ~ 0.107
    ETA_v1 ~ 0.227
    ETA_ka ~ 0.464
    ETA_baf ~ 0.049
    add.sd <- 0.465
    prop.sd <- 0.143
  })
  model({
    TVcl  = THETA_cl*(ClCr/57);
    TVv1  = THETA_v1;
    TVcld = THETA_cld;
    TVv2  = THETA_v2;
    TVka  = THETA_ka;
    TVbaf = THETA_baf;

    cl  = TVcl*exp(ETA_cl);
    v1  = TVv1*exp(ETA_v1);
    cld = TVcld;
    v2  = TVv2;
    ka  = TVka*exp(ETA_ka);
    baf = TVbaf*exp(ETA_baf);

    k10 = cl/v1;
    k12 = cld / v1;
    k21 = cld / v2;
    Cc = centr/v1;

    d/dt(depot)  = -ka*depot
    d/dt(centr)  =  ka*depot - k10*centr - k12*centr + k21*periph;
    d/dt(periph) =                         k12*centr - k21*periph;
    d/dt(AUC)    = Cc;

    f(depot)=baf;
    alag(depot)=0.382;

    Cc ~ add(add.sd) + prop(prop.sd) + combined1()
  })
}

## ----tdm_patientA-------------------------------------------------------------
patient <- data.frame(ID=1,TIME=c(0,121,122,126,144),
                      DV=c(NA,10.8,5.8,3.3,NA),
                      ADDL=c(5,0,0,0,0),
                      II=c(24,0,0,0,0),
                      EVID=c(1,0,0,0,1),
                      CMT=c("centr",NA,NA,NA,"centr"),
                      AMT=c(250,0,0,0,0),
                      DUR=c(0.5,NA,NA,NA,NA),
                      ClCr=25)
patient

## ----estim_map----------------------------------------------------------------
map_patient <- poso_estim_map(patient,mod_ganciclovir_Caldes_2009)

## ----plot_map-----------------------------------------------------------------
plot(map_patient$model,Cc)

## ----auc_24_ganci-------------------------------------------------------------
library(data.table)
data.table(map_patient$model)[time==144,AUC] - 
  data.table(map_patient$model)[time==120,AUC]

## ----optimal_intravenous_ganci------------------------------------------------
poso_dose_auc(patient,mod_ganciclovir_Caldes_2009,tdm=TRUE,
              time_dose = 145,
              duration = 1,
              time_auc = 24,
              target_auc = 50,
              cmt_dose = "centr")

## ----optimal_oral_valganci----------------------------------------------------
poso_dose_auc(patient,mod_ganciclovir_Caldes_2009,tdm=TRUE,
              time_dose = 145,
              time_auc = 24,
              target_auc = 50,
              cmt_dose = "depot")

