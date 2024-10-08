#-------------------------------------------------------------------------
# posologyr: individual dose optimization using population PK
# Copyright (C) Cyril Leven
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as
#  published by the Free Software Foundation, either version 3 of the
#  License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
#  Adapted from: http://shiny.webpopix.org/mcmc/bayes1/
#  Copyright (C) Marc Lavielle, Inria Saclay, CeCILL-B
#
#  Modifications:
#   - interfacing with rxode2
#   - deletion of shiny-specific parts
#   - variable names changed to snake_case
#   - square matrix taken as input, not diagonal
#   - functions return values for both etas and theta
#   - inter-occasion variability (IOV)
#-------------------------------------------------------------------------

#' Estimate the prior distribution of population parameters
#'
#' Estimates the prior distribution of population parameters by Monte Carlo
#' simulations
#'
#' @param dat Dataframe. An individual subject dataset following the
#'     structure of NONMEM/rxode2 event records.
#' @param prior_model A \code{posologyr} prior population pharmacokinetics
#'    model, a list of six objects.
#' @param n_simul An integer, the number of simulations to be run. For `n_simul
#'   =0`, all ETAs are set to 0.
#' @param return_model A boolean. Returns a rxode2 model using the simulated
#'    ETAs if set to `TRUE`.
#'
#' @return If `return_model` is set to `FALSE`, a list of one element: a
#' dataframe `$eta` of the individual values of ETA.
#' If `return_model` is set to `TRUE`, a list of the dataframe of the
#' individual values of ETA, and a rxode2 model using the simulated ETAs.
#'
#' @examples
#' # model
#' mod_run001 <- function() {
#'   ini({
#'     THETA_Cl <- 4.0
#'     THETA_Vc <- 70.0
#'     THETA_Ka <- 1.0
#'     ETA_Cl ~ 0.2
#'     ETA_Vc ~ 0.2
#'     ETA_Ka ~ 0.2
#'     prop.sd <- sqrt(0.05)
#'   })
#'   model({
#'     TVCl <- THETA_Cl
#'     TVVc <- THETA_Vc
#'     TVKa <- THETA_Ka
#'
#'     Cl <- TVCl*exp(ETA_Cl)
#'     Vc <- TVVc*exp(ETA_Vc)
#'     Ka <- TVKa*exp(ETA_Ka)
#'
#'     K20 <- Cl/Vc
#'     Cc <- centr/Vc
#'
#'     d/dt(depot) = -Ka*depot
#'     d/dt(centr) = Ka*depot - K20*centr
#'     Cc ~ prop(prop.sd)
#'   })
#' }
#' # df_patient01: event table for Patient01, following a 30 minutes intravenous
#' # infusion
#' df_patient01 <- data.frame(ID=1,
#'                         TIME=c(0.0,1.0,14.0),
#'                         DV=c(NA,25.0,5.5),
#'                         AMT=c(2000,0,0),
#'                         EVID=c(1,0,0),
#'                         DUR=c(0.5,NA,NA))
#' # estimate the prior distribution of population parameters
#' poso_simu_pop(dat=df_patient01,prior_model=mod_run001,n_simul=100)
#'
#' @export
poso_simu_pop <- function(dat=NULL,prior_model=NULL,n_simul=1000,
                          return_model=TRUE){
  prior_model <- get_prior_model(prior_model)
  object <- posologyr(prior_model,dat)
  no_covariates <- is.null(object$covariates)

  omega      <- object$omega
  ind_eta    <- which(diag(omega)>0)          # only parameters with IIV
  omega_eta  <- omega[ind_eta,ind_eta,drop=FALSE]
  eta_mat    <- matrix(0,nrow=1,ncol=ncol(omega))

  if (n_simul > 0) {
    eta_mat <- matrix(0,nrow=n_simul,ncol=ncol(omega))
    eta_sim <- mvtnorm::rmvnorm(n_simul,mean=rep(0,ncol(omega_eta)),
                                sigma=omega_eta)
    eta_mat[,ind_eta] <- eta_sim
  }

  eta_df             <- data.frame(eta_mat)
  names(eta_df)      <- attr(omega,"dimnames")[[1]]

  eta_pop            <- list(eta=eta_df)

  # outputs
  if(return_model){
    model_pop         <- object$solved_ppk_model
    theta             <- rbind(object$theta)

    if(no_covariates){
      params <- cbind(theta,eta_df,row.names=NULL)
    } else {
      covar             <- as.data.frame(object$tdm_data[1,object$covariates])
      names(covar)      <- object$covariates
      params <- cbind(theta,eta_df,covar,row.names=NULL)
    }

    if (!is.null(object$pi_matrix)){
      kappa_mat         <- matrix(0,nrow=1,ncol=ncol(omega))
      kappa_df          <- data.frame(kappa_mat)
      names(kappa_df)   <- attr(object$pi_matrix,"dimnames")[[1]]
      params            <- cbind(params,kappa_df)
    }

    model_pop$params  <- params
    eta_pop$model     <- model_pop
  }

  return(eta_pop)
}

#' Estimate the Maximum A Posteriori individual parameters
#'
#' Estimates the Maximum A Posteriori (MAP) individual parameters,
#' also known as Empirical Bayes Estimates (EBE).
#'
#' @param dat Dataframe. An individual subject dataset following the
#'     structure of NONMEM/rxode2 event records.
#' @param prior_model A \code{posologyr} prior population pharmacokinetics
#'    model, a list of six objects.
#' @param return_model A boolean. Returns a rxode2 model using the estimated
#'    ETAs if set to `TRUE`.
#' @param return_ofv A boolean. Returns a the Objective Function Value (OFV)
#'    if set to `TRUE`.
#' @param nocb A boolean. for time-varying covariates: the next observation
#'     carried backward (nocb) interpolation style, similar to NONMEM.  If
#'     `FALSE`, the last observation carried forward (locf) style will be used.
#'     Defaults to `FALSE`.
#'
#' @return A named list consisting of one or more of the following elements
#' depending on the input parameters of the function: `$eta` a named vector
#' of the MAP estimates of the individual values of ETA, `$model` an rxode2
#' model using the estimated ETAs, `$event` the `data.table` used to solve the
#' returned rxode2 model.
#'
#' @import data.table
#'
#' @importFrom stats setNames
#'
#' @examples
#' rxode2::setRxThreads(1) # limit the number of threads
#'
#' # model
#' mod_run001 <- function() {
#'   ini({
#'     THETA_Cl <- 4.0
#'     THETA_Vc <- 70.0
#'     THETA_Ka <- 1.0
#'     ETA_Cl ~ 0.2
#'     ETA_Vc ~ 0.2
#'     ETA_Ka ~ 0.2
#'     prop.sd <- sqrt(0.05)
#'   })
#'   model({
#'     TVCl <- THETA_Cl
#'     TVVc <- THETA_Vc
#'     TVKa <- THETA_Ka
#'
#'     Cl <- TVCl*exp(ETA_Cl)
#'     Vc <- TVVc*exp(ETA_Vc)
#'     Ka <- TVKa*exp(ETA_Ka)
#'
#'     K20 <- Cl/Vc
#'     Cc <- centr/Vc
#'
#'     d/dt(depot) = -Ka*depot
#'     d/dt(centr) = Ka*depot - K20*centr
#'     Cc ~ prop(prop.sd)
#'   })
#' }
#' # df_patient01: event table for Patient01, following a 30 minutes intravenous
#' # infusion
#' df_patient01 <- data.frame(ID=1,
#'                         TIME=c(0.0,1.0,14.0),
#'                         DV=c(NA,25.0,5.5),
#'                         AMT=c(2000,0,0),
#'                         EVID=c(1,0,0),
#'                         DUR=c(0.5,NA,NA))
#' # estimate the Maximum A Posteriori individual parameters
#' poso_estim_map(dat=df_patient01,prior_model=mod_run001)
#'
#' @export
poso_estim_map <- function(dat=NULL,prior_model=NULL,return_model=TRUE,
                           return_ofv=FALSE,nocb=FALSE){
  prior_model <- get_prior_model(prior_model)
  object <- posologyr(prior_model,dat,nocb)
  endpoints <- get_endpoints(object)
  estim_with_iov <- check_for_iov(object)
  no_covariates  <- is.null(object$covariates)

  dat           <- object$tdm_data
  solved_model  <- object$solved_ppk_model
  omega         <- object$omega
  theta         <- object$theta
  sigma         <- object$sigma
  error_model   <- object$error_model
  interpolation <- object$interpolation

  ind_eta      <- which(diag(omega)>0)               # only parameters with IIV
  omega_eta    <- omega[ind_eta,ind_eta,drop=FALSE]  # only variances > 0
  solve_omega  <- try(solve(omega_eta))              # inverse of omega_eta

  eta_map      <- diag(omega)*0

  diag_varcovar_matrix <- diag(omega_eta)

  # initialize the list of outputs
  estim_map    <- list(eta=eta_map)

  #standard MAP estimation

    # avoid empty (NULL) arguments for stats::optim()
    omega_dim <- 0
    iov_col   <- 0
    pimat     <- 0
    eta_df    <- 0

    if (estim_with_iov){
      data_iov     <- dat
      pimat        <- object$pi_matrix

      ind_kappa    <- which(diag(pimat)>0)

      omega_dim    <- ncol(omega_eta)

      if(length(ind_kappa)==1){
        pimat_kappa  <- pimat
        pimat_dim    <- 1
      } else{
        pimat_kappa  <- pimat[ind_kappa,ind_kappa]
        pimat_dim    <- ncol(pimat_kappa)
      }

      iov_col      <- init_iov_col(dat=dat,pimat=pimat)
      all_the_mat  <- merge_covar_matrices(omega_eta=omega_eta,
                                           omega_dim=omega_dim,
                                           pimat_dim=pimat_dim,
                                           pimat_kappa=pimat_kappa,
                                           dat=dat)

      solve_omega   <- try(solve(all_the_mat))
      diag_varcovar_matrix <- diag(all_the_mat)
    }

    start_eta        <- init_eta(object,estim_with_iov,
                                 omega_iov=all_the_mat,endpoints=endpoints)

    ifelse(setequal(endpoints,"Cc"),
           y_obs <- data.frame(DV=dat[dat$EVID==0,"DV"],DVID="Cc"),
           y_obs <- dat[dat$EVID==0,c("DV","DVID")])

    # initial bounds for the optimization
    bfgs_bounds       <- stats::qnorm(25e-3,0,sqrt(diag_varcovar_matrix),
                                      lower.tail = F)

    optim_attempt     <- 1
    max_attempt       <- 40
    one_more_time     <- TRUE

    # create a table to log the estimates after each attempt
    optim_attempt_log        <- matrix(Inf,nrow=max_attempt,
                                       ncol=(1+length(start_eta)))
    optim_attempt_log        <- data.table::data.table(optim_attempt_log)

    data.table::setnames(optim_attempt_log,
                         1:(length(start_eta)+1),
                         c("OFV",names(start_eta)),
                         skip_absent=TRUE)

    while(one_more_time & optim_attempt <= max_attempt){
      r <- try(stats::optim(start_eta,errpred,
                        run_model=run_model,
                        y_obs=y_obs,
                        endpoints=endpoints,
                        theta=theta,
                        ind_eta=ind_eta,
                        sigma=sigma,
                        solve_omega=solve_omega,
                        omega=omega,
                        omega_dim=omega_dim,
                        iov_col=iov_col,
                        pimat=pimat,
                        dat=dat,
                        solved_model=solved_model,
                        error_model=error_model,
                        estim_with_iov=estim_with_iov,
                        interpolation=interpolation,
                        method="L-BFGS-B",
                        upper=bfgs_bounds,
                        lower=-bfgs_bounds),
               silent=TRUE)

      if(!inherits(r,'try-error')){

        # Objection Function Value: OFV
        OFV_current       <- r$value

        # log the detection of the anomalous estimates, OFV, and ETA estimates
        optim_attempt_log[optim_attempt,] <- data.table("OFV"=OFV_current,
                                                        rbind(r$par))
        best_attempt_ofv  <- min(optim_attempt_log$OFV)
        second_best_ofv   <- min(sort(optim_attempt_log$OFV)[-1])

        # detection of anomalous estimates calling for a new attempt
        stuck_on_bound    <- TRUE %in% (abs(r$par) >= bfgs_bounds)
        all_eta_are_zero  <- !(FALSE %in% (r$par == 0))
        identical_abs_eta <- isTRUE(length(unique(abs(r$par))) < length(r$par))
        sky_high_ofv      <- isTRUE(OFV_current >= 1e10)
        not_the_best      <- isTRUE(OFV_current - best_attempt_ofv > 1e-7)
        far_from_2nd_best <- isTRUE(abs(second_best_ofv -
                                          best_attempt_ofv) > 1e-5)

        need_a_new_start  <- isTRUE(optim_attempt == 1|
                                      all_eta_are_zero|
                                      identical_abs_eta|
                                      sky_high_ofv|
                                      not_the_best|
                                      far_from_2nd_best)

        if(optim_attempt < max_attempt){

          one_more_time <- FALSE

          # check conditions calling for a new attempt at minimizing the OFV

          if(stuck_on_bound){
            one_more_time <- TRUE
            bfgs_bounds   <- bfgs_bounds + 1

          } else if(need_a_new_start){
            one_more_time <- TRUE
            start_eta     <- init_eta(object,estim_with_iov,
                                      omega_iov=all_the_mat,
                                      endpoints=endpoints)
          }

        } else if(optim_attempt == max_attempt){

          # if all fails, the "less bad" solution is probably
          # the estimation with the lowest OFV

          OFV   <- NULL    # avoid undefined global variables

          r$par <- unlist(optim_attempt_log[OFV==min(OFV),
                                              2:(length(start_eta)+1)][1,])
          r$value <- min(optim_attempt_log[,OFV])

        }
      } else{ # class(r) == 'try-error'
        one_more_time <- TRUE
        start_eta     <- init_eta(object,estim_with_iov,
                                  omega_iov=all_the_mat,
                                  endpoints=endpoints)
      }

      optim_attempt     <- optim_attempt + 1
    }

    if (estim_with_iov){
      eta_map[ind_eta] <- r$par[1:omega_dim]
    } else {
      eta_map[ind_eta] <- r$par
    }

    if(!no_covariates){
      covar            <- t(dat[1,object$covariates]) #results in a matrix
      names(covar)     <- object$covariates
    }

  # list of all outputs
  estim_map$eta      <- eta_map

  if(return_model){
    et_poso <- rxode2::as.et(object$tdm_data)
    et_poso$clearSampling()
    et_poso$add.sampling(seq(dat$TIME[1],
                             dat$TIME[length(dat$TIME)]+1,
                             by=.1))

    #model_map        <- solved_model
    if(estim_with_iov){
      iov_kappa <- attr(pimat,"dimnames")[[1]]
      iov_col <- iov_proposition_as_cols(iov_col=iov_col,dat=dat,pimat=pimat,
                                         omega_dim=omega_dim,
                                         eta_estim=r$par)
      data_iov <- data.frame(dat,iov_col)
      iov_kappa_mat <- sapply(iov_kappa,FUN=extrapol_iov,dat=data_iov,
                              iov_kappa=iov_kappa,
                              event_table=et_poso)

      et_poso <- cbind(et_poso,iov_kappa_mat)
    }
    if(!no_covariates){
      covar_mat <- sapply(object$covariates,FUN=extrapol_cov,dat=dat,
                          covar=object$covariates,
                          interpol_approx="constant",
                          f=ifelse(object$interpolation == "nocb",1,0),
                          event_table=et_poso)

      et_poso <- cbind(et_poso,covar_mat)
    }
    estim_map$model <- rxode2::rxSolve(object$ppk_model,et_poso,
                                       c(object$theta,estim_map$eta),
                                       covsInterpolation = interpolation)
    estim_map$event <- data.table::data.table(et_poso)
  }

  if(return_ofv){
    estim_map$ofv <- r$value
  }

  return(estim_map)
}

# poso_estim_mcmc: distribution of the individual parameters
#
# Adapted from the saemix R package
# Copyright (C) Emmanuelle Comets <emmanuelle.comets@inserm.fr>,
# Audrey Lavenu, Marc Lavielle (authors of the saemix R package)
#
# The saemix package is free software; licensed under the terms of the GNU
# General Public License as published by the Free Software Foundation;
# either version 2 of the License, or (at your option) any later version.
#
# (GPLv2+)

#' Estimate the posterior distribution of individual parameters by MCMC
#'
#' Estimates the posterior distribution of individual parameters by Markov
#' Chain Monte Carlo (using a Metropolis-Hastings algorithm)
#'
#' @param dat Dataframe. An individual subject dataset following the
#'     structure of NONMEM/rxode2 event records.
#' @param prior_model A \code{posologyr} prior population pharmacokinetics
#'    model, a list of six objects.
#' @param return_model A boolean. Returns a rxode2 model using the estimated
#'    ETAs if set to `TRUE`.
#' @param burn_in Number of burn-in iterations for the Metropolis-Hastings
#'    algorithm.
#' @param n_iter Total number of iterations (following the burn-in iterations)
#'  for each Markov chain of the Metropolis-Hastings algorithm.
#' @param n_chains Number of Markov chains
#' @param nocb A boolean. for time-varying covariates: the next observation
#'     carried backward (nocb) interpolation style, similar to NONMEM.  If
#'     `FALSE`, the last observation carried forward (locf) style will be used.
#'     Defaults to `FALSE`.
#' @param control A list of parameters controlling the Metropolis-Hastings
#' algorithm.
#'
#' @return If `return_model` is set to `FALSE`, a list of one element: a
#' dataframe `$eta` of ETAs from the posterior distribution, estimated by
#' Markov Chain Monte Carlo.
#' If `return_model` is set to `TRUE`, a list of the dataframe of the posterior
#' distribution of ETA, and a rxode2 model using the estimated distributions of
#' ETAs.
#'
#' @author Emmanuelle Comets, Audrey Lavenu, Marc Lavielle, Cyril Leven
#'
#' @references Comets  E, Lavenu A, Lavielle M. Parameter estimation in
#' nonlinear mixed effect models using saemix, an R implementation of the SAEM
#' algorithm. Journal of Statistical Software 80, 3 (2017), 1-41.
#'
#' @examples
#' # model
#' mod_run001 <- function() {
#'   ini({
#'     THETA_Cl <- 4.0
#'     THETA_Vc <- 70.0
#'     THETA_Ka <- 1.0
#'     ETA_Cl ~ 0.2
#'     ETA_Vc ~ 0.2
#'     ETA_Ka ~ 0.2
#'     prop.sd <- sqrt(0.05)
#'   })
#'   model({
#'     TVCl <- THETA_Cl
#'     TVVc <- THETA_Vc
#'     TVKa <- THETA_Ka
#'
#'     Cl <- TVCl*exp(ETA_Cl)
#'     Vc <- TVVc*exp(ETA_Vc)
#'     Ka <- TVKa*exp(ETA_Ka)
#'
#'     K20 <- Cl/Vc
#'     Cc <- centr/Vc
#'
#'     d/dt(depot) = -Ka*depot
#'     d/dt(centr) = Ka*depot - K20*centr
#'     Cc ~ prop(prop.sd)
#'   })
#' }
#' # df_patient01: event table for Patient01, following a 30 minutes intravenous
#' # infusion
#' df_patient01 <- data.frame(ID=1,
#'                         TIME=c(0.0,1.0,14.0),
#'                         DV=c(NA,25.0,5.5),
#'                         AMT=c(2000,0,0),
#'                         EVID=c(1,0,0),
#'                         DUR=c(0.5,NA,NA))
#' # estimate the posterior distribution of population parameters
#' \donttest{poso_estim_mcmc(dat=df_patient01,prior_model=mod_run001,
#' n_iter=50,n_chains=2)}
#'
#' @export
poso_estim_mcmc <- function(dat=NULL,prior_model=NULL,return_model=TRUE,
                            burn_in=50,n_iter=1000,n_chains=4,nocb=FALSE,
                            control=list(n_kernel=c(2,2,2),
                            stepsize_rw=0.4,proba_mcmc=0.3,nb_max=3)){
  prior_model <- get_prior_model(prior_model)
  object <- posologyr(prior_model,dat,nocb)
  if(check_for_iov(object)) stop("IOV is not supported, poso_estim_sir() can be
                                 used instead to estimate the posterior
                                 distributions")
  endpoints <- get_endpoints(object)
  no_covariates <- is.null(object$covariates)

  dat          <- object$tdm_data
  solved_model <- object$solved_ppk_model
  omega        <- object$omega
  theta        <- object$theta

  one_mcmc_please <- function(chains,object,burn_in,n_iter,control){

    dat          <- object$tdm_data
    solved_model <- object$solved_ppk_model
    omega        <- object$omega
    sigma        <- object$sigma
    error_model  <- object$error_model

    ifelse(setequal(endpoints,"Cc"),
           y_obs <- data.frame(DV=dat[dat$EVID==0,"DV"],DVID="Cc"),
           y_obs <- dat[dat$EVID==0,c("DV","DVID")])
    ind_eta      <- which(diag(omega)>0)      # only parameters with IIV
    nb_etas      <- length(ind_eta)
    omega_eta    <- omega[ind_eta,ind_eta,drop=FALSE]    # only variances > 0
    solve_omega  <- try(solve(omega_eta))     # inverse of omega_eta
    chol_omega   <- chol(omega_eta)
    rw_init      <- 0.5                  #initial variance parameter for kernels
    d_omega      <- diag(omega_eta)*rw_init
    VK           <- rep(c(1:nb_etas),2)
    n_iter       <- n_iter + burn_in

    # Metropolis-Hastings algorithm------------------------------------------
    theta    <- object$theta
    eta      <- diag(omega_eta)*0
    f_all_endpoints <- do.call(run_model,list(c(theta,eta),
                                              solved_model=solved_model,
                                              estim_with_iov=FALSE,
                                              endpoints=endpoints))

    obs_res <- residual_error_all_endpoints(f_all_endpoints=f_all_endpoints,
                                            y_obs=y_obs,
                                            error_model=error_model,
                                            sigma=sigma,
                                            endpoints=endpoints)
    f <- obs_res$f_all_endpoints$f
    g <- obs_res$g_all_endpoints$g
    g[which(g == 0)] <- 1 # avoid NaN, idem issue #28

    U_y      <- 0.5*sum(((y_obs[,"DV"] - f)/g)^2 + log(g^2))
    U_eta    <- 0.5 * eta %*% solve_omega %*% eta

    eta_mat     <- matrix(0,nrow=n_iter+1,ncol=ncol(omega))
    eta_mat[1,] <- diag(omega)*0

    for (k_iter in 1:n_iter)
    {
      if (control$n_kernel[1] > 0)
      {
        for (u in 1:control$n_kernel[1])
        {
          etac          <- as.vector(chol_omega%*%stats::rnorm(nb_etas))
          names(etac)   <- attr(omega_eta,"dimnames")[[1]]
          f_all_endpoints <- do.call(run_model,list(c(theta,etac),
                                                    solved_model=solved_model,
                                                    estim_with_iov=FALSE,
                                                    endpoints=endpoints))

          obs_res <- residual_error_all_endpoints(f_all_endpoints=f_all_endpoints,
                                                  y_obs=y_obs,
                                                  error_model=error_model,
                                                  sigma=sigma,
                                                  endpoints=endpoints)
          f <- obs_res$f_all_endpoints$f
          g <- obs_res$g_all_endpoints$g
          g[which(g == 0)] <- 1 # avoid NaN, idem issue #28

          Uc_y          <- 0.5*sum(((y_obs[,"DV"] - f)/g)^2 + log(g^2))
          deltu         <- Uc_y - U_y
          if(deltu < (-1) * log(stats::runif(1)))
          {
            eta        <- etac
            U_y        <- Uc_y
          }
        }
      }
      if (control$n_kernel[2] > 0)
      {
        nb_max        <- min(nb_etas,control$nb_max)
        nbc2          <- nt2     <- replicate(nb_etas,0)
        U_eta         <- 0.5 * eta %*% solve_omega %*% eta
        for (u in 1:control$n_kernel[2])
        {
          for (nrs2 in 1:nb_max)
          {
            for (j in 1:nb_etas)
            {
              jr            <- sample(c(1:nb_etas), nrs2)
              jr            <- jr -jr[1] + j
              vk2           <- jr%%nb_etas + 1
              etac          <- eta
              etac[vk2]     <- eta[vk2] + stats::rnorm(nrs2)*d_omega[vk2]
              f_all_endpoints <- do.call(run_model,
                                         list(c(theta,etac),
                                              solved_model=solved_model,
                                              estim_with_iov=FALSE,
                                              endpoints=endpoints))

              obs_res <- residual_error_all_endpoints(f_all_endpoints=f_all_endpoints,
                                                      y_obs=y_obs,
                                                      error_model=error_model,
                                                      sigma=sigma,
                                                      endpoints=endpoints)
              f <- obs_res$f_all_endpoints$f
              g <- obs_res$g_all_endpoints$g
              g[which(g == 0)] <- 1 # avoid NaN, idem issue #28

              Uc_y          <- 0.5*sum(((y_obs[,"DV"] - f)/g)^2 + log(g^2))
              Uc_eta        <- 0.5 * etac %*% solve_omega %*% etac
              deltu         <- Uc_y - U_y + Uc_eta - U_eta
              if(deltu < (-1) * log(stats::runif(1)))
              {
                eta         <- etac
                U_y         <- Uc_y
                U_eta       <- Uc_eta
                nbc2[vk2]   <- nbc2[vk2]+1
              }
              nt2[vk2]      <- nt2[vk2] + 1
            }
          }
        }
        d_omega <- d_omega*(1 + control$stepsize_rw*(nbc2/nt2 - control$proba_mcmc))
      }
      if(control$n_kernel[3]>0) {
        nt2          <- nbc2     <- matrix(data=0,nrow=nb_etas,ncol=1)
        nrs2         <- k_iter%%(nb_etas-1)+2
        for (u in 1:control$n_kernel[3]) {
          if(nrs2<nb_etas) {
            vk        <- c(0,sample(c(1:(nb_etas-1)),nrs2-1))
            nb_iter2  <- nb_etas
          } else {
            vk        <- 0:(nb_etas-1)
            nb_iter2  <- 1
          }
          for(k2 in 1:nb_iter2) {
            vk2             <- VK[k2+vk]
            etac            <- eta
            etac[vk2]       <- eta[vk2]+matrix(stats::rnorm(nrs2),
                                               ncol=nrs2)%*%diag(d_omega[vk2])
            f_all_endpoints <- do.call(run_model,list(c(theta,etac),
                                                      solved_model=solved_model,
                                                      estim_with_iov=FALSE,
                                                      endpoints=endpoints))

            obs_res <- residual_error_all_endpoints(f_all_endpoints=f_all_endpoints,
                                                    y_obs=y_obs,
                                                    error_model=error_model,
                                                    sigma=sigma,
                                                    endpoints=endpoints)
            f <- obs_res$f_all_endpoints$f
            g <- obs_res$g_all_endpoints$g
            g[which(g == 0)] <- 1 # avoid NaN, idem issue #28

            Uc_y            <- 0.5*sum(((y_obs[,"DV"] - f)/g)^2 + log(g^2))
            Uc_eta          <- 0.5*rowSums(etac*(etac%*%solve(omega_eta)))
            deltu           <- Uc_y-U_y+Uc_eta-U_eta
            ind             <- which(deltu<(-log(stats::runif(1))))
            eta[ind]        <- etac[ind]
            U_y[ind]        <- Uc_y[ind]
            U_eta[ind]      <- Uc_eta[ind]
            nbc2[vk2]       <- nbc2[vk2]+length(ind)
            nt2[vk2]        <- nt2[vk2]+1
          }
        }
        d_omega <- d_omega*(1+control$stepsize_rw * (nbc2/nt2-control$proba_mcmc))
      }
      eta_mat[k_iter+1,ind_eta]   <- eta
    }
    return(eta_mat)
  }

  split_iter_id_burn_in <- function(chain,n_iter,burn_in){
    iter_to_keep <- chain*((burn_in+1):(n_iter+burn_in))
    return(iter_to_keep)
  }

  eta_mat_list <- lapply((1:n_chains),one_mcmc_please,object=object,
                         burn_in=burn_in,n_iter=n_iter,control=control)

  iter_to_keep <- sapply((1:n_chains),FUN=split_iter_id_burn_in,
                         n_iter=n_iter,burn_in=burn_in)

  eta_mat                <- do.call(rbind,eta_mat_list)
  eta_df_mcmc            <- data.frame(eta_mat[iter_to_keep,])
  names(eta_df_mcmc)     <- attr(omega,"dimnames")[[1]]

  estim_mcmc             <- list(eta=eta_df_mcmc)

  if(return_model){
    model_mcmc        <- solved_model
    theta_return      <- rbind(theta)
    if(no_covariates){
      model_mcmc$params <- cbind(theta_return,eta_df_mcmc,row.names=NULL)
    } else {
      covar             <- as.data.frame(dat[1,object$covariates])
      names(covar)      <- object$covariates
      model_mcmc$params <- cbind(theta_return,eta_df_mcmc,covar,row.names=NULL)
    }
    estim_mcmc$model  <- model_mcmc
  }

  return(estim_mcmc)
}

#' Estimate the posterior distribution of individual parameters by SIR
#'
#' Estimates the posterior distribution of individual parameters by
#' Sequential Importance Resampling (SIR)
#'
#'
#' @param dat Dataframe. An individual subject dataset following the
#'     structure of NONMEM/rxode2 event records.
#' @param prior_model A \code{posologyr} prior population pharmacokinetics
#'    model, a list of six objects.
#' @param n_sample Number of samples from the S-step
#' @param n_resample Number of samples from the R-step
#' @param return_model A boolean. Returns a rxode2 model using the estimated
#'    ETAs if set to `TRUE`.
#' @param nocb A boolean. for time-varying covariates: the next observation
#'     carried backward (nocb) interpolation style, similar to NONMEM.  If
#'     `FALSE`, the last observation carried forward (locf) style will be used.
#'     Defaults to `FALSE`.
#'
#' @return If `return_model` is set to `FALSE`, a list of one element: a
#' dataframe `$eta` of ETAs from the posterior distribution, estimated by
#' Sequential Importance Resampling.
#' If `return_model` is set to `TRUE`, a list of the dataframe of the posterior
#' distribution of ETA, and a rxode2 model using the estimated distributions of
#' ETAs.
#'
#' @import data.table
#' @examples
#' # model
#' mod_run001 <- function() {
#'   ini({
#'     THETA_Cl <- 4.0
#'     THETA_Vc <- 70.0
#'     THETA_Ka <- 1.0
#'     ETA_Cl ~ 0.2
#'     ETA_Vc ~ 0.2
#'     ETA_Ka ~ 0.2
#'     prop.sd <- sqrt(0.05)
#'   })
#'   model({
#'     TVCl <- THETA_Cl
#'     TVVc <- THETA_Vc
#'     TVKa <- THETA_Ka
#'
#'     Cl <- TVCl*exp(ETA_Cl)
#'     Vc <- TVVc*exp(ETA_Vc)
#'     Ka <- TVKa*exp(ETA_Ka)
#'
#'     K20 <- Cl/Vc
#'     Cc <- centr/Vc
#'
#'     d/dt(depot) = -Ka*depot
#'     d/dt(centr) = Ka*depot - K20*centr
#'     Cc ~ prop(prop.sd)
#'   })
#' }
#' # df_patient01: event table for Patient01, following a 30 minutes intravenous
#' # infusion
#' df_patient01 <- data.frame(ID=1,
#'                         TIME=c(0.0,1.0,14.0),
#'                         DV=c(NA,25.0,5.5),
#'                         AMT=c(2000,0,0),
#'                         EVID=c(1,0,0),
#'                         DUR=c(0.5,NA,NA))
#' # estimate the posterior distribution of population parameters
#' poso_estim_sir(dat=df_patient01,prior_model=mod_run001,
#' n_sample=1e3,n_resample=1e2)
#'
#' @export
poso_estim_sir <- function(dat=NULL,prior_model=NULL,n_sample=1e4,
                           n_resample=1e3,return_model=TRUE,nocb=FALSE){

  prior_model <- get_prior_model(prior_model)

  object <- posologyr(prior_model,dat,nocb)
  endpoints <- get_endpoints(object)
  estim_with_iov <- check_for_iov(object)
  no_covariates  <- is.null(object$covariates)

  dat           <- object$tdm_data
  solved_model  <- object$solved_ppk_model
  omega         <- object$omega
  sigma         <- object$sigma
  error_model   <- object$error_model
  interpolation <- object$interpolation

  ifelse(endpoints=="Cc",
         y_obs <- data.frame(DV=dat[dat$EVID==0,"DV"],DVID="Cc"),
         y_obs <- dat[dat$EVID==0,c("DV","DVID")])
  ind_eta      <- which(diag(omega)>0)      # only parameters with IIV
  nb_etas      <- length(ind_eta)
  omega_eta    <- omega[ind_eta,ind_eta,drop=FALSE]    # only variances > 0
  omega_sim    <- omega_eta
  omega_dim    <- ncol(omega_eta)
  solve_omega  <- try(solve(omega_eta))     # inverse of omega_eta

  theta        <- rbind(object$theta)

  if (estim_with_iov){
    pimat        <- object$pi_matrix

    ind_kappa    <- which(diag(pimat)>0)

    if(length(ind_kappa)==1){
      pimat_kappa  <- pimat
      pimat_names  <- attr(pimat_kappa,"dimnames")[[1]]
      pimat_dim    <- 1
    } else{
      pimat_kappa  <- pimat[ind_kappa,ind_kappa]
      pimat_names  <- attr(pimat_kappa,"dimnames")[[1]]
      pimat_dim    <- ncol(pimat_kappa)
    }

    n_occ        <- length(unique(dat$OCC))

    all_the_mat  <- merge_covar_matrices(omega_eta=omega_eta,
                                         omega_dim=omega_dim,
                                         pimat_dim=pimat_dim,
                                         pimat_kappa=pimat_kappa,
                                         dat=dat)
    omega_sim   <- all_the_mat        # needed for I-Step
    solve_omega <- solve(all_the_mat) # needed for LL_func in I-Step
  }

  #SIR algorithm
  # doi: 10.1002/psp4.12492; doi: 10.1007/s10928-016-9487-8

  #S-step
  eta_sim       <- mvtnorm::rmvnorm(n_sample,mean=rep(0,ncol(omega_sim)),
                                    sigma=omega_sim)

  if (estim_with_iov){
    # matrix large enough for omega + pi
    eta_mat           <- matrix(0,nrow=n_sample,
                                ncol=ncol(all_the_mat)-
                                  ncol(omega[ind_eta,ind_eta,drop=FALSE])+
                                  ncol(omega))

    # IIV
    eta_mat[,ind_eta] <-
      eta_sim[,1:length(ind_eta)]

    # IOV
    eta_mat[,(ncol(omega)+1):ncol(eta_mat)] <-
      eta_sim[,(length(ind_eta)+1):ncol(eta_sim)]

  } else{
    eta_mat           <- matrix(0,nrow=n_sample,ncol=ncol(omega))
    eta_mat[,ind_eta] <- eta_sim
  }

  eta_df            <- data.frame(eta_mat)

  names(eta_df) <- attr(omega,"dimnames")[[1]]
  eta_dt        <- data.table::data.table(eta_df)

  param_cols    <- attr(omega,"dimnames")[[1]]
  params        <- cbind(ID=1,eta_dt[,param_cols,with=F],theta)

  if (estim_with_iov){
    eta_dt[,ID:=(1:n_sample)]                           # 1:n_samples ID

    dat_dt <- data.table::data.table(dat)
    dat_dt <- dat_dt[rep(dat_dt[,.I],n_sample)]         # one table per sample
    dat_dt[,ID:=rep(1:n_sample,1,each=nrow(dat))]       # 1:n_samples IDs

    # bind random effects to patient data.table
    ID     <- NULL    # avoid undefined global variables
    dat_dt <- dat_dt[eta_dt,on = list(ID = ID), roll = TRUE]

    # everything but ID, OCC and kappas
    tdm_dt            <- data.table::data.table(dat)
    names_tdm_dt_drop <- names(tdm_dt[,!c("ID","OCC")])
    names_tdm_dt_full <- names(tdm_dt)
    drop_cols         <- c(names_tdm_dt_drop,attr(omega,"dimnames")[[1]])

    # reshape IOV to colums
    iov_col <- t(apply(dat_dt[,!drop_cols,with=F],
                       MARGIN=1,
                       FUN=link_kappa_to_occ,
                       pimat_dim=pimat_dim,
                       pimat_names=pimat_names))
    iov_col_dt <- data.table::data.table(iov_col)

    data_iov   <- cbind(dat_dt[,names_tdm_dt_full,with=F],
                        iov_col_dt[,!c("ID","OCC")])

    param_cols <- c("ID",attr(omega,"dimnames")[[1]])
    params     <- cbind(eta_dt[,param_cols,with=F],theta)
  }

  #I-step
  if(estim_with_iov){
    group_index   <- data.frame(cbind(c(1:10),10,n_sample,nrow(dat)))

    loads_omodels <- apply(group_index,
                          pkmodel=object$solved_ppk_model,
                          params=params,
                          dat=data_iov,
                          interpolation=interpolation,
                          MARGIN=1,FUN=solve_by_groups)
    f_all_endpoints  <- do.call(rbind,loads_omodels)

    data.table::setnames(f_all_endpoints,"id","sim.id")
  } else {
    f_all_endpoints <- rxode2::rxSolve(solved_model,
                                       cbind(theta,eta_dt,row.names=NULL),
                                       dat,covsInterpolation=interpolation,
                                       returnType="data.table")
  }

  obs_res <- residual_error_all_endpoints(f_all_endpoints=f_all_endpoints,
                                          y_obs=y_obs,
                                          error_model=error_model,
                                          sigma=sigma,
                                          endpoints=endpoints)

  f_all_sim <- dcast(obs_res$f_all_endpoints, formula = sim.id ~ rowid(sim.id),
                     value.var = "f")
  g_all_sim <- dcast(obs_res$g_all_endpoints, formula = sim.id ~ rowid(sim.id),
                     value.var = "g")

  LL_func  <- function(simu_obs){ #doi: 10.4196/kjpp.2012.16.2.97
    eta_id   <- simu_obs[1]
    eta      <- eta_sim[eta_id,]
    f        <- simu_obs[-1]
    g        <- g_all_sim[eta_id,-1]
    minus_LL <- 0.5*objective_function(y_obs=y_obs$DV,f=f,g=g,eta=eta,
                                       solve_omega=solve_omega)
    return(-minus_LL)
  }

  lf        <- apply(f_all_sim,MARGIN=1,FUN=LL_func)
  lp        <- mvtnorm::dmvnorm(eta_sim,mean=rep(0,ncol(omega_sim)),
                                sigma=omega_sim,log=TRUE)
  md        <- max(lf - lp)
  wt        <- exp(lf - lp - md)
  probs     <- wt/sum(wt)

  #R-step
  indices   <- sample(1:n_sample, size = n_resample, prob = probs,
                      replace = TRUE)
  if (nb_etas > 1) {
    eta_sim <- eta_sim[indices, ]
  }
  else {
    eta_sim <- eta_sim[indices]
  }

  eta_mat           <- matrix(0,nrow=n_resample,ncol=ncol(omega))
  eta_mat[,ind_eta] <- eta_sim[,1:omega_dim]
  eta_df            <- data.frame(eta_mat)
  names(eta_df)     <- attr(omega,"dimnames")[[1]]

  estim_sir         <- list(eta=eta_df)

  if(return_model){
    if(estim_with_iov){
      params_resample   <- cbind(data.frame(ID=1:n_resample),eta_df,theta)

      # return a list of data.tables of resampled IDs
      loads_otables     <- lapply(indices,
                                  data_iov,
                                  FUN=function(indices,dat)
                                    {dat[ID == indices,]})

      # bind the data.tables nicely
      dat_resample      <- do.call(rbind,loads_otables)

      # overwrite IDs to avoid duplicates, and solve the model once again
      dat_resample[,ID:=rep(1:n_resample,1,each=nrow(dat))]
      estim_sir$model   <- rxode2::rxSolve(object$solved_ppk_model,
                                        params_resample,
                                        dat_resample,
                                        covsInterpolation=interpolation)
    } else {
      params_resample <- cbind(eta_df,theta,row.names=NULL)

      if(no_covariates){
        model_sir     <- rxode2::rxSolve(object$solved_ppk_model,
                                        params_resample,
                                        dat,covsInterpolation=interpolation)
      } else {
        covar         <- as.data.frame(dat[1,object$covariates])
        names(covar)  <- object$covariates
        model_sir     <- rxode2::rxSolve(object$solved_ppk_model,
                                        cbind(params_resample,covar,
                                              row.names=NULL),
                                        dat,covsInterpolation=interpolation)
      }
      estim_sir$model   <- model_sir
    }
  }
  return(estim_sir)
}
