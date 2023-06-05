## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----minimal_patient_df-------------------------------------------------------
data.frame(ID=1,
           TIME=c(0.0,3),
           DV=c(NA,60.0),
           AMT=c(1000,0),
           EVID=c(101,0))

## ----before_dosing_df---------------------------------------------------------
data.frame(ID=1,
           TIME=0,
           DV=0,
           AMT=0,
           EVID=0,
           COVAR1=c("X"),
           COVAR2=c("Y"))

## ----oral_IV_bolus_df---------------------------------------------------------
data.frame(ID=1,
           TIME=c(0.0,3),
           DV=c(NA,60.0),
           AMT=c(1000,0),
           EVID=c(1,0),
           COVAR1=c("X"),
           COVAR2=c("Y"))

## ----intermittent_infusion_df-------------------------------------------------
data.frame(ID=1,
           TIME=c(0.0,1.0,14.0),
           DV=c(NA,25.0,5.5),
           AMT=c(1000,0,0),
           DUR=c(0.5,NA,NA),
           EVID=c(1,0,0),
           COVAR1=c("X"),
           COVAR2=c("Y"))

## ----iov_df-------------------------------------------------------------------
data.frame(ID=1,
           TIME=c(0.0,1.0,14.0,24.0,25.0,36.0),
           DV=c(NA,25.0,5.5,NA,30.0,6.0),
           AMT=c(1000,0,0,1000,0,0),
           DUR=c(0.5,NA,NA,0.5,NA,NA),
           EVID=c(1,0,0,1,0,0),
           OCC=c(1,1,1,2,2,2),
           COVAR1=c("X"),
           COVAR2=c("Y"))

## ----endpoints----------------------------------------------------------------
data.frame(ID=1,
           TIME=c(0.0,1.0,14.0,24.0,25.0,36.0),
           DV=c(NA,20.0,80,35.5,60.0,40.0),
           AMT=c(1000,0,0,0,0,0),
           EVID=c(1,0,0,0,0,0),
           DVID=c("parent","parent","metabolite","parent","metabolite","metabolite"),
           COVAR1=c("X"),
           COVAR2=c("Y"))

