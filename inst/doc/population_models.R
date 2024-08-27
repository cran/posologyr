## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----gentamicin---------------------------------------------------------------
mod_gentamicin_Xuan2003 <- function() {
  ini({
    #Fixed effects: population estimates 
    THETA_Cl  = 0.047
    THETA_V   = 0.28
    THETA_k12 = 0.092
    THETA_k21 = 0.071
    
    #Random effects: inter-individual variability
    ETA_Cl  ~ 0.084
    ETA_V   ~ 0.003
    ETA_k12 ~ 0.398
    ETA_k21 ~ 0.342
    
    #Unexplained residual variability
    add_sd  <- 0.230
    prop_sd <- 0.237
  })
  model({
    #Individual model and covariates
    TVl   = THETA_Cl*ClCr
    TVV   = THETA_V*WT
    TVk12 = THETA_k12
    TVk21 = THETA_k21
    Cl    = TVl*exp(ETA_Cl)
    V     = TVV*exp(ETA_V)
    k12   = TVk12*exp(ETA_k12)
    k21   = TVk21*exp(ETA_k21)
    
    #Structural model defined using ordinary differential equations (ODE)
    ke    = Cl/V
    Cp    = centr/V

    d/dt(centr)  = - ke*centr - k12*centr + k21*periph
    d/dt(periph) =            + k12*centr - k21*periph

    #Model for unexplained residual variability
    Cp ~ add(add_sd) + prop(prop_sd) + combined1()
  })
}

## ----amikacin-----------------------------------------------------------------
mod_amikacin_Burdet2015 <- function() {
    ini({
      #Fixed effects: population estimates 
      THETA_Cl=4.3
      THETA_Vc=15.9
      THETA_Vp=21.4
      THETA_Q=12.1
      
      #Random effects: inter-individual variability
      ETA_Cl + ETA_Vc + ETA_Vp + ETA_Q ~
        c(0.1,
          0.01     ,    0.05 ,
          0.01     ,    0.02 ,   0.2  ,
         -0.06     ,    0.004,   0.003,    0.08)
          
      #Unexplained residual variability
      add_sd <- 0.2
      prop_sd <- 0.1
    })
    model({
      #Individual model and covariates
      TVCl  = THETA_Cl*(CLCREAT4H/82)^0.7
      TVVc  = THETA_Vc*(TBW/78)^0.9*(PoverF/169)^0.4
      TVVp  = THETA_Vp
      TVQ   = THETA_Q
      Cl    = TVCl*exp(ETA_Cl)
      Vc    = TVVc*exp(ETA_Vc)
      Vp    = TVVp*exp(ETA_Vp)
      Q     = TVQ *exp(ETA_Q)
      
      #Structural model defined using ordinary differential equations (ODE)
      ke    = Cl/Vc
      k12   = Q/Vc
      k21   = Q/Vp
      Cp    = centr/Vc      
      
      d/dt(centr)  = - ke*centr - k12*centr + k21*periph
      d/dt(periph) =            + k12*centr - k21*periph

      #Model for unexplained residual variability
      Cp ~ add(add_sd) + prop(prop_sd) + combined1()
    })
  }

## ----ganciclovir--------------------------------------------------------------
mod_ganciclovir_Caldes2009 <- function() {
  ini({
    #Fixed effects: population estimates 
    THETA_cl  <- 7.49
    THETA_v1  <- 31.90
    THETA_cld <- 10.20
    THETA_v2  <- 32.0
    THETA_ka  <- 0.895
    THETA_baf <- 0.825
      
    #Random effects: inter-individual variability
    ETA_cl ~ 0.107
    ETA_v1 ~ 0.227
    ETA_ka ~ 0.464
    ETA_baf ~ 0.049
    
    #Unexplained residual variability
    add.sd <- 0.465
    prop.sd <- 0.143
  })
  model({
    #Individual model and covariates
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

    #Structural model defined using ordinary differential equations (ODE)
    k10 = cl/v1;
    k12 = cld / v1;
    k21 = cld / v2;
    Cc = centr/v1;

    d/dt(depot)  = -ka*depot
    d/dt(centr)  =  ka*depot - k10*centr - k12*centr + k21*periph;
    d/dt(periph) =                         k12*centr - k21*periph;
    d/dt(AUC)    = Cc;

    #Special model event changes
    f(depot)=baf;
    alag(depot)=0.382;
    
    #Model for unexplained residual variability
    Cc ~ add(add.sd) + prop(prop.sd) + combined1()
  })
}

