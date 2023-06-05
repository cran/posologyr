## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
set.seed(1)

## ----setup--------------------------------------------------------------------
library(posologyr)

## ----model--------------------------------------------------------------------
mod_amikacin_Burdet2015 <- list(
  ppk_model   = rxode2::rxode({
    centr(0) = 0;
    TVCl  = THETA_Cl*(CLCREAT4H/82)^0.7;
    TVVc  = THETA_Vc*(TBW/78)^0.9*(PoverF/169)^0.4;
    TVVp  = THETA_Vp;
    TVQ   = THETA_Q;
    Cl    = TVCl*exp(ETA_Cl);
    Vc    = TVVc*exp(ETA_Vc);
    Vp    = TVVp*exp(ETA_Vp);
    Q     = TVQ *exp(ETA_Q);
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
  theta = c(THETA_Cl=4.3, THETA_Vc=15.9, THETA_Vp=21.4,THETA_Q=12.1),
  omega = lotri::lotri({ETA_Cl + ETA_Vc + ETA_Vp + ETA_Q ~
      c(0.1,
        0.01     ,   0.05 ,
        0.01     ,   0.02 ,   0.2  ,
        -0.06    ,   0.004,   0.003,    0.08)}),
  covariates  = c("CLCREAT4H","TBW","PoverF"),
  sigma       = c(additive_a = 0.2, proportional_b = 0.1))

## ----prior_patientA-----------------------------------------------------------
df_patientA <- data.frame(ID=1,TIME=0,
                                DV=0,
                                EVID=0,
                                AMT=0,
                                CLCREAT4H=50,TBW=62,PoverF=169)
df_patientA

## ----prior_Cmax_optim_patientA------------------------------------------------
prior_dose <- poso_dose_conc(dat=df_patientA,
                             prior_model=mod_amikacin_Burdet2015,
                             time_c = 1,                        #30 min after a  
                             duration = 0.5,                    #30 min infusion
                             target_conc = 80)
prior_dose

## ----prior_Cmin_optim_patientA------------------------------------------------
poso_time_cmin(dat=df_patientA,
               prior_model=mod_amikacin_Burdet2015,
               dose = prior_dose$dose,
               duration = 0.5,                                  #30 min infusion
               target_cmin = 2.5)

## ----prior_plot_model---------------------------------------------------------
# generate a model using the individual covariates 
simu_patA      <- poso_simu_pop(dat=df_patientA,
                                prior_model=mod_amikacin_Burdet2015,
                                n_simul = 0)

## -----------------------------------------------------------------------------
simu_patA$model$time <- seq(0,20,b=0.1)
simu_patA$model$add.dosing(dose=prior_dose$dose,rate=prior_dose$dose/0.5)

## -----------------------------------------------------------------------------
plot(simu_patA$model,Cc)

## ----prior_plot---------------------------------------------------------------
plot(simu_patA$model,Cc) + 
  ggplot2::ylab("Central concentration") +
  ggplot2::geom_vline(xintercept=1, linetype="dashed") +
  ggplot2::geom_ribbon(ggplot2::aes(ymin=60, ymax=80),
                       fill="seagreen",show.legend = FALSE, alpha=0.15)

