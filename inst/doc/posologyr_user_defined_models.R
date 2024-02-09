## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----ppk_model----------------------------------------------------------------
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
  })

## ----error_model--------------------------------------------------------------
error_model <- function(f,sigma){     #additive model if sigma[2] == 0
  g <- sigma[1] + sigma[2]*f          #proportional model if sigma[1] == 0
  return(g)
}

## ----error_model_nonmem-------------------------------------------------------
error_model <- function(f,sigma){
  dv <- cbind(f,1)
  g  <- diag(dv%*%sigma%*%t(dv))     #sigma is the square matrix of the residual
  return(sqrt(g))                    #errors
}

## ----theta--------------------------------------------------------------------
theta = c(THETA_Cl=4.5, THETA_Vc=58.4, THETA_Vp=38.4, THETA_Q=6.5)

## ----omega--------------------------------------------------------------------
omega = lotri::lotri({ETA_Cl + ETA_Vc + ETA_Vp + ETA_Q ~
                          c(0.147,
                            0    ,  0.510 ,
                            0    ,  0     ,   0.282,
                            0    ,  0     ,   0    ,    0)})


## ----sigma--------------------------------------------------------------------
sigma       = c(additive_a = 3.4, proportional_b = 0.227)

## ----matrix-------------------------------------------------------------------
sigma       = lotri::lotri({prop + add ~ c(0.227,0.0,3.4)})

## ----named_list---------------------------------------------------------------
sigma       = list(
    cp=c(additive_a = 0.144, proportional_b = 0.15),
    pca=c(additive_a = 3.91, proportional_b = 0.0)
    )

## ----pi_matrix----------------------------------------------------------------
pi_matrix = lotri::lotri({KAPPA_Cl + KAPPA_Vc ~
      c(0.1934626,
        0.00     ,  0.05783106)})


## ----covariates---------------------------------------------------------------
covariates  = c("CLCREAT","WT","DIAL")

## ----vancomycin_pososlogyr_model_list-----------------------------------------
mod_vancomyin_Goti2018 <- list(
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
                            0    ,  0.510 ,
                            0    ,  0     ,   0.282,
                            0    ,  0     ,   0    ,    0)}),
  sigma       = c(additive_a = 3.4, proportional_b = 0.227),
  covariates  = c("CLCREAT","WT","DIAL"))

## ----vancomycin_pososlogyr_model----------------------------------------------
mod_vancomyin_Goti2018

