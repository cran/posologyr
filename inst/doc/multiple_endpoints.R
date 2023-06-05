## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(posologyr)

## -----------------------------------------------------------------------------
mod_warfarin_nlmixr <- list(
  ppk_model   = rxode2::rxode({
    ktr <- exp(THETA_ktr + ETA_ktr)
    ka <- exp(THETA_ka + ETA_ka)
    cl <- exp(THETA_cl + ETA_cl)
    v <- exp(THETA_v + ETA_v)
    emax = expit(THETA_emax + ETA_emax)
    ec50 =  exp(THETA_ec50 + ETA_ec50)
    kout = exp(THETA_kout + ETA_kout)
    e0 = exp(THETA_e0 + ETA_e0)
    ##
    DCP = center/v
    PD=1-emax*DCP/(ec50+DCP)
    ##
    effect(0) = e0
    kin = e0*kout
    ##
    d/dt(depot) = -ktr * depot
    d/dt(gut) =  ktr * depot -ka * gut
    d/dt(center) =  ka * gut - cl / v * center
    d/dt(effect) = kin*PD -kout*effect
    ##
    cp = center / v
    pca = effect
  }),
  error_model = list(
    cp = function(f,sigma){
      g <- sigma[1]^2 + (sigma[2]^2)*(f^2)
      return(sqrt(g))
    },
    pca = function(f,sigma){
      g <- sigma[1]^2 + (sigma[2]^2)*(f^2)
      return(sqrt(g))
    }
  ),
  theta = c(THETA_ktr=0.106,
            THETA_ka=-0.087,
            THETA_cl=-2.03,
            THETA_v=2.07,
            THETA_emax=3.4,
            THETA_ec50=0.00724,
            THETA_kout=-2.9,
            THETA_e0=4.57),
  omega = lotri::lotri({ETA_ktr + ETA_ka + ETA_cl + ETA_v + ETA_emax + ETA_ec50 +
      ETA_kout + ETA_e0 ~
      c(1.024695,
        0.00     ,   0.9518403 ,
        0.00     ,   0.00 ,   0.5300943  ,
        0.00    ,    0.00,   0.00,   0.4785394,
        0.00    ,    0.00,   0.00,   0.00,   0.7134424,
        0.00    ,    0.00,   0.00,   0.00,   0.00,   0.7204165,
        0.00    ,    0.00,   0.00,   0.00,   0.00,   0.00,   0.3563706,
        0.00    ,    0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.2660827)}),
  sigma       = list(
    cp=c(additive_a = 0.144, proportional_b = 0.15),
    pca=c(additive_a = 3.91, proportional_b = 0.0)
    )
  )

## -----------------------------------------------------------------------------
warf_01 <- data.frame(ID=1,
                      TIME=c(0.0,0.5,1.0,2.0,3.0,6.0,9.0,12.0,24.0,24.0,36.0,
                             36.0,48.0,48.0,72.0,72.0,96.0,120.0,144.0),
                      DV=c(0.0,0.0,1.9,3.3,6.6,9.1,10.8,8.6,5.6,44.0,4.0,27.0,
                           2.7,28.0,0.8,31.0,60.0,65.0,71.0),
                      DVID=c("cp","cp","cp","cp","cp","cp","cp","cp","cp","pca",
                              "cp","pca","cp","pca","cp","pca","pca","pca","pca"),
                      EVID=c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                      AMT=c(100,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
warf_01

## -----------------------------------------------------------------------------
map_warf_01  <- poso_estim_map(warf_01,mod_warfarin_nlmixr)
map_warf_01

## -----------------------------------------------------------------------------
plot(map_warf_01$model,"cp")

## -----------------------------------------------------------------------------
plot(map_warf_01$model,"pca")

