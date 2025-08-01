# ==============================================================================
# Sample Size Calculation for Single‑Arm HIV Trial with Counterfactual Comparator
# ------------------------------------------------------------------------------
# Last Modified: 01 August 2025
#
# PURPOSE
#   Compute sample size, power, and variance components for trials that compare
#   an active PrEP arm against a *counterfactual* placebo HIV‑incidence estimate
#   derived from recency‑assay data.  The statistical framework follows
#   Gao et al. (2021), and equation numbers below refer to that paper.
#
# ==============================================================================

################################################################################
## Input parameters from PWID India study 
################################################################################
#  p            HIV prevalence at baseline (screening for trial)
#  OmegaT       MDRI (years)
#  betaT        False recent rate (FRR)
#  R1 / R0      Alternative / null incidence rate ratios
#  tau          Clinical follow-up time (years)
#  Time         Cut-off that defines “recent” (years)
#  sigmaOmegaT  MDRI SE (years)
#  sigmabetaT   FRR SE 
#  r            Proportion of HIV-negative subjects enrolled to the trial
################################################################################

# --------------------------------------------------------------------------
# Numerically safe utilities
eps <- 1e-8
clip01 <- function(x) pmin(pmax(x, eps), 1 - eps)


################################################################################
## Core statistical functions
################################################################################

#### Sample-size calculation based on Z (Gao eqn 4) ###########################
samplesize_CF <- function(para, R1, R0 = 1, alpha = 0.05, power = 0.9) {
  gamma0   <- rec_gamma(para)                       # counterfactual component
  gamma1   <- trial_gamma(para, para$lambda0 * R1)  # active-arm component

  gamma00  <- gamma0$gamma00
  gamma01  <- gamma0$gamma01
  gamma1v  <- gamma1$gamma1

  vZ       <- varZ(para, R1, R0)
  Zsum     <- qnorm(1 - alpha / 2) + qnorm(power, 0, sqrt(vZ))

  denom <- ((log(R1) - log(R0)) / Zsum)^2 - gamma01
  if (denom <= 0) return(Inf)                       ## FIX: infeasibility guard

  ## Equation (4): required screening size in presence of VR1 uncertainty
  ceiling((gamma00 + gamma1v) / denom)
}

#### Power calculation (expected power via Gao eqn 3) #########################
power_CF <- function(para, R1, R0 = 1, alpha = 0.05, N) {
  gamma0  <- rec_gamma(para, N)
  gamma1  <- trial_gamma(para, para$lambda0 * R1, N)
  gamma00 <- gamma0$gamma00
  gamma01 <- gamma0$gamma01
  gamma1v <- gamma1$gamma1

  qbeta   <- abs((log(R1) - log(R0)) /
                 sqrt((gamma00 + gamma1v) / N + gamma01)) -
             qnorm(1 - alpha / 2)

  vZ      <- varZ(para, R1, R0)
  ## probability(|Z|>z_{1-α/2}) under H1
  pnorm(qbeta, 0, sqrt(vZ))
}

################################################################################
## Data Generation (uses φ for recency-test coverage) 
################################################################################
generatedata <- function(N, para, R = 1) {
  p  <- para$p
  r  <- para$r
  lambda0 <- para$lambda0
  lambda1 <- lambda0 * R
  tau     <- para$tau

  OmegaT       <- para$OmegaT
  sigmaOmegaT  <- para$sigmaOmegaT
  betaT        <- para$betaT
  sigmabetaT   <- para$sigmabetaT
  Time         <- para$Time

  Npos <- rbinom(1, N, p)
  Nneg <- N - Npos

  Npos_test  <- rbinom(1, Npos, para$q)     
  Nneg_trial <- rbinom(1, Nneg, r)

  PR <- betaT + lambda0 * (1 - p) / p * (OmegaT - betaT * Time)
  PR <- clip01(PR)
  NR <- rbinom(1, Npos_test, PR)

  wbetaT  <- rnorm(1, betaT, sigmabetaT)
  wOmegaT <- rnorm(1, OmegaT, sigmaOmegaT)

  Nevent <- rpois(1, lambda1 * tau * Nneg_trial)

  list(Npos        = Npos,
       Nneg        = Nneg,
       Npos_test   = Npos_test,
       Nneg_trial  = Nneg_trial,
       NR          = NR,
       Nevent      = Nevent,
       PR          = PR,
       wbetaT      = wbetaT,
       wOmegaT     = wOmegaT)
}

################################################################################
## λ0 estimation (Gao eqn 1; uses empirical coverage q̂ )
################################################################################
rec_est <- function(N, Npos, Npos_test, NR, para) {
  pest <- Npos / N
  qhat <- Npos_test / Npos
  Nneg <- N - Npos        
  NR_adj    <- NR / qhat

  betaT  <- para$betaT
  OmegaT <- para$OmegaT
  Time   <- para$Time

  ## Equation (1) estimator of counterfactual incidence
  lambda_est <- (NR_adj - betaT * Npos) /
                (Nneg * (OmegaT - betaT * Time))

  para_est <- list(lambda0      = lambda_est,
                   p            = pest,
                   q            = qhat,
                   OmegaT       = OmegaT,
                   sigmaOmegaT  = para$sigmaOmegaT,
                   betaT        = betaT,
                   sigmabetaT   = para$sigmabetaT,
                   Time         = Time)

  varlog <- rec_gamma(para_est, N)$varlog
  SE     <- sqrt(varlog) * lambda_est
  if (qhat == 0) return(list(Est = NA, SE = NA, CI = NA, CI_log = NA))

  list(Est   = lambda_est,
       SE    = SE,
       CI    = lambda_est + c(-1, 1) * qnorm(0.975) * SE,
       CI_log = lambda_est * exp(c(-1, 1) * qnorm(0.975) *
                                 sqrt(varlog)))
}

################################################################################
## λ1 estimation (Gao eqn 2)
################################################################################
trial_est <- function(N, Nneg, Nneg_trial, Nevent, tau) {
  pest <- 1 - Nneg / N
  rest <- Nneg_trial / Nneg

  ## Equation (2) estimator for active-arm incidence
  lambda_est <- Nevent / (Nneg_trial * tau)

  para_est <- list(p = pest, r = rest, tau = tau)
  varlog   <- trial_gamma(para_est, lambda_est, N)$varlog
  SE       <- sqrt(varlog) * lambda_est
  if (!is.finite(lambda_est)) return(list(Est = NA, SE = NA, CI = NA, CI_log = NA))

  list(Est   = lambda_est,
       SE    = SE,
       CI    = lambda_est + c(-1, 1) * qnorm(0.975) * SE,
       CI_log = lambda_est * exp(c(-1, 1) * qnorm(0.975) *
                                 sqrt(varlog)))
}

################################################################################
## Hypothesis-test utilities
################################################################################
Rejection <- function(lambda0_est, lambda1_est,
                      R0 = 1, Zcut = qnorm(0.975)) {
  if (lambda0_est$Est <= 0) return(FALSE)
  if (lambda1_est$Est == 0) {
    Zstat <- (log(lambda0_est$Est)) /
             (lambda0_est$SE / lambda0_est$Est)
    return(abs(Zstat) > Zcut)
  }
  Zstat <- efficacy_test(lambda0_est, lambda1_est, R0)
  abs(Zstat) > Zcut
}

efficacy_test <- function(lambda0_est, lambda1_est, R0 = 1) {
  lambda0 <- lambda0_est$Est
  V0      <- (lambda0_est$SE / lambda0)^2
  lambda1 <- lambda1_est$Est
  V1      <- (lambda1_est$SE / lambda1)^2
  wR      <- lambda1 / lambda0
  (log(wR) - log(R0)) / sqrt(V0 + V1)
}

################################################################################
## Expected number helpers
################################################################################
rec_num <- function(N, para) {
  lambda0 <- para$lambda0
  p       <- para$p
  q       <- para$q
  OmegaT  <- para$OmegaT
  betaT   <- para$betaT
  Time    <- para$Time
  PR      <- betaT + lambda0 * (1 - p) / p * (OmegaT - betaT * Time)
  N * p * PR * q
}

inc_num <- function(N, lambda1, para) {
  p   <- para$p
  r   <- para$r
  tau <- para$tau
  N * tau * (1 - p) * r * lambda1
}

summarize_expectnum <- function(N, para, R) {
  p     <- para$p
  q     <- para$q
  MDRI  <- para$OmegaT
  FRR   <- para$betaT
  Time  <- para$Time
  r     <- para$r
  tau   <- para$tau
  lambda0 <- para$lambda0
  c(N,
    N * p * q,
    rec_num(N, para),
    N * (1 - p) * r,
    inc_num(N, lambda0 * R, para))
}

################################################################################
## Variance components (γ₀₀ and γ₀₁ from Gao Appx A)
################################################################################
rec_gamma <- function(para, N = NULL) {
  lambda0 <- para$lambda0
  p       <- para$p
  q       <- para$q       
  OmegaT  <- para$OmegaT
  betaT   <- para$betaT
  Time    <- para$Time
  sigma2Omega <- para$sigmaOmegaT^2
  sigma2beta  <- para$sigmabetaT^2

  PR <- betaT + lambda0 * (1 - p) / p * (OmegaT - betaT * Time)
  PR <- clip01(PR)
  den1 <- PR - betaT
  if (abs(den1) < eps)
    return(list(gamma00 = Inf, gamma01 = Inf, varlog = Inf))

  gamma00 <- ( PR * (1 - PR) / den1^2 / q +                
               1 / (1 - p) +
               (1 - p * q) * sigma2beta / den1^2 / q ) / p 

  gamma01 <- sigma2Omega / (OmegaT - betaT * Time)^2 +
             sigma2beta *
             ((OmegaT - PR * Time) / den1 / (OmegaT - betaT * Time))^2

  varlog  <- if (is.null(N)) NA else gamma00 / N + gamma01
  list(gamma00 = gamma00, gamma01 = gamma01, varlog = varlog)
}

trial_gamma <- function(para, lambda1, N = NULL) {
  p   <- para$p
  r   <- para$r
  tau <- para$tau
  gamma1 <- 1 / (lambda1 * (1 - p) * r * tau)
  if (!is.finite(gamma1)) gamma1 <- Inf
  varlog <- if (is.null(N)) NA else gamma1 / N
  list(gamma1 = gamma1, varlog = varlog)
}

#### Covariance matrix varW  #####################
varW <- function(p, q, PR, OmegaT, betaT, r, lambda1, tau) {
  variance <- matrix(0, 7, 7)
  variance[1,1] <- p * q * PR * (1 - PR) + p * q * (1 - p * q) * (PR - betaT)^2
  variance[2,2] <- p * q * (1 - p * q)
  variance[3,3] <- p * (1 - p)
  variance[5,5] <- (1 - p) * r * lambda1 * tau *
                   (1 + lambda1 * tau * (1 - r) + lambda1 * tau * r * p)
  variance[6,6] <- (1 - p) * r * (1 - r + p * r)
  variance[7,7] <- p * q * PR * (1 - PR * q * p)

  variance[1,2] <- variance[2,1] <- p * q * (1 - p * q) * (PR - betaT)
  variance[1,3] <- variance[3,1] <- p * (1 - p) * q * (PR - betaT)
  variance[1,5] <- variance[5,1] <- -p * (1 - p) * (PR - betaT) * q * r * lambda1 * tau
  variance[1,6] <- variance[6,1] <- -p * (1 - p) * (PR - betaT) * q * r
  variance[1,7] <- variance[7,1] <- p * q * PR * (1 - PR) +
                                    p * q * (1 - p * q) * (PR - betaT) * PR
  variance[2,3] <- variance[3,2] <- p * (1 - p) * q
  variance[2,5] <- variance[5,2] <- -p * (1 - p) * q * r * lambda1 * tau
  variance[2,6] <- variance[6,2] <- -p * (1 - p) * q * r
  variance[2,7] <- variance[7,2] <- p * q * (1 - p * q) * PR
  variance[3,5] <- variance[5,3] <- -p * (1 - p) * r * lambda1 * tau
  variance[3,6] <- variance[6,3] <- -p * (1 - p) * r
  variance[3,7] <- variance[7,3] <- p * (1 - p) * q * PR
  variance[5,6] <- variance[6,5] <- (1 - p) * r * (1 - r + p * r) * lambda1 * tau
  variance[5,7] <- variance[7,5] <- -p * (1 - p) * PR * q * r * lambda1 * tau
  variance[6,7] <- variance[7,6] <- -p * (1 - p) * PR * q * r
  variance
}

#### Gradient matrix var_AB – Appendix B derivatives ##########################
var_AB <- function(p, q, PR, OmegaT, betaT, r, lambda1, tau) {
  den1 <- PR - betaT
  if (abs(den1) < eps) return(matrix(Inf, 2, 2))
  a <- matrix(0, 7, 2)

  a[1,1] <- -1 / (p * q * den1)
  a[2,1] <-  1 / (p * q)
  a[3,1] <- -1 / (p * (1 - p))
  a[4,1] <-  1 /  OmegaT
  a[5,1] <-  1 / ((1 - p) * r * lambda1 * tau)
  a[6,1] <- -1 / ((1 - p) * r)
  
  a[1,2] <- -2 * PR * (1 - PR) / (p * q)^2 / den1^3
  a[2,2] <- (PR^2) / (p * q * den1)^2
  a[3,2] <- -(p^2 + (1 - p)^2) / (p * (1 - p))^2
  a[4,2] <-  0
  a[5,2] <- -1 / ((1 - p) * r * lambda1 * tau)^2
  a[6,2] <-  0
  a[7,2] <- (1 - 2 * PR) / (p * q * den1)^2

  cov_mat <- varW(p, q, PR, OmegaT, betaT, r, lambda1, tau)
  t(a) %*% cov_mat %*% a
}

#### var_AB_sigmabeta – Appendix B with assay SE ##############################
var_AB_sigmabeta <- function(p, q, PR, OmegaT, betaT, r, lambda1,
                              tau, sigmabetaT) {
  den1 <- PR - betaT
  if (abs(den1) < eps) return(matrix(Inf, 2, 2))
  a <- matrix(0, 7, 2)

  a[1,1] <- -1 / (p * q * den1)
  a[2,1] <-  1 / (p * q)
  a[3,1] <- -1 / (p * (1 - p))
  a[4,1] <-  1 /  OmegaT
  a[5,1] <-  1 / ((1 - p) * r * lambda1 * tau)
  a[6,1] <- -1 / ((1 - p) * r)
  
  a[1,2] <- -2 * PR * (1 - PR) / (p * q)^2 / den1^3 -
             2 * sigmabetaT^2 * (1 - p) / (p * q)^2 / den1^3
  a[2,2] <- (PR^2 + sigmabetaT^2 * (1 - 2 * p)) / (p * q * den1)^2
  a[3,2] <- -(p^2 + (1 - p)^2) / (p * (1 - p))^2
  a[4,2] <-  0
  a[5,2] <- -1 / ((1 - p) * r * lambda1 * tau)^2
  a[6,2] <-  0
  a[7,2] <- (1 - 2 * PR) / (p * q * den1)^2

  cov_mat <- varW(p, q, PR, OmegaT, betaT, r, lambda1, tau)
  t(a) %*% cov_mat %*% a
}

#### varZ uses corrected var_AB (Appendix B variance of Z) ####################
varZ <- function(para, R1, R0 = 1) {
  p       <- para$p
  q       <- para$q
  OmegaT  <- para$OmegaT
  betaT   <- para$betaT
  r       <- para$r
  tau     <- para$tau
  Time    <- para$Time
  lambda0 <- para$lambda0
  lambda1 <- lambda0 * R1

  PR <- betaT + lambda0 * (1 - p) / p * (OmegaT - betaT * Time)

  varAB <- var_AB(p, q, PR, OmegaT, betaT, r, lambda1, tau)

  tA <- log(R1) - log(R0)
  tB <- PR * (1 - PR) / (p * q * (PR - betaT)^2) +
        1 / p + 1 / (1 - p) + 1 / ((1 - p) * r * lambda1 * tau)

  cvec <- c(1 / sqrt(tB), -tA / (2 * tB^(3/2)))
  drop(t(cvec) %*% varAB %*% cvec)
}

################################################################################
##  5×5 and 6×6 covariance kernels + σΩ / σβ propagation 
################################################################################

## ---- 5 × 5 covariance (no assay SE terms; Appendix B) ----------------------
varW5 <- function(p, PR, OmegaT, betaT, r, lambda1, tau) {
  V <- matrix(0, 5, 5)
  V[1,1] <- p*PR*(1-PR) + p*(1-p)*(PR-betaT)^2
  V[2,2] <- p*(1-p)
  V[3,3] <- (1-p)*r*lambda1*tau *
            (1 + lambda1*tau*(1-r) + lambda1*tau*r*p)
  V[4,4] <- (1-p)*r*(1-r + p*r)
  V[5,5] <- p*PR*(1 - PR*p)

  V[1,2] <- V[2,1] <- p*(1-p)*(PR-betaT)
  V[1,3] <- V[3,1] <- -p*(1-p)*(PR-betaT)*r*lambda1*tau
  V[1,4] <- V[4,1] <- -p*(1-p)*(PR-betaT)*r
  V[1,5] <- V[5,1] <- p*PR*(1-PR) + p*(1-p)*(PR-betaT)*PR
  V[2,3] <- V[3,2] <- -p*(1-p)*r*lambda1*tau
  V[2,4] <- V[4,2] <- -p*(1-p)*r
  V[2,5] <- V[5,2] <- p*(1-p)*PR
  V[3,4] <- V[4,3] <- (1-p)*r*(1-r+p*r)*lambda1*tau
  V[3,5] <- V[5,3] <- -p*(1-p)*PR*r*lambda1*tau
  V[4,5] <- V[5,4] <- -p*(1-p)*PR*r
  V
}

## ---- 6 × 6 covariance (adds assay SE block; Appendix B) --------------------
varW6 <- function(p, PR, OmegaT, betaT, r, lambda1, tau,
                  sigmaOmegaT, sigmabetaT, Time, N) {

  V <- varW5(p, PR, OmegaT, betaT, r, lambda1, tau)
  V <- rbind(cbind(V, rep(0, 5)), c(rep(0, 5), 0))  # expand to 6 × 6
  V[1,6] <- V[6,1] <- p*sigmabetaT^2 * Time         # cross-cov term
  V <- V * N                                        # scale ∝ N
  V[6,6] <- sigmaOmegaT^2 + sigmabetaT^2 * Time^2
  V
}

################################################################################
## var_Z_sigmaOmega and var_Z_both (propagate assay SE; Appx B)
################################################################################

var_Z_sigmaOmega <- function(p, q, PR, OmegaT, betaT, r, lambda1,
                             tau, sigmaOmegaT, Time, N,
                             R1 = 1, R0 = 1) {
  V <- varW5(p, PR, OmegaT, betaT, r, lambda1, tau) * N
  a <- matrix(0, 5, 2)
  a[1,1] <- -1 / (p * (PR - betaT)) / N
  a[2,1] <- -1 / (1 - p) / N
  a[3,1] <-  1 / ((1 - p) * r * lambda1 * tau) / N
  a[4,1] <- -1 / ((1 - p) * r) / N

  a[1,2] <- -2*PR*(1-PR)/p^2/(PR-betaT)^3/N^2
  a[2,2] <- (PR^2)/(p*(PR-betaT))^2 /N^2 -
            (p^2+(1-p)^2)/(p*(1-p))^2 /N^2
  a[3,2] <- -1/((1-p)*r*lambda1*tau)^2 /N^2
  a[5,2] <- (1-2*PR)/(p*(PR-betaT))^2 /N^2

  mat <- t(a) %*% V %*% a
  last <- c(1/(OmegaT-betaT*Time),
            -2*sigmaOmegaT^2/(OmegaT-betaT*Time)^3)
  mat  <- mat + last %*% t(last) * sigmaOmegaT^2

  tA <- log(R1) - log(R0)
  tB <- (PR*(1-PR)/(p*q*(PR-betaT)^2) + 1/p + 1/(1-p) + 1/((1-p)*r*lambda1*tau))/N + sigmaOmegaT^2/(OmegaT - betaT*Time)^2
  cvec <- c(1/sqrt(tB), -tA/(2*tB^(3/2)))
  drop(t(cvec) %*% mat %*% cvec)
}

var_Z_both <- function(p, q, PR, OmegaT, betaT, r, lambda1, tau,
                       sigmaOmegaT, sigmabetaT, Time, N,
                       R1 = 1, R0 = 1) {
  V <- varW6(p, PR, OmegaT, betaT, r, lambda1, tau,
             sigmaOmegaT, sigmabetaT, Time, N)
  a <- matrix(0, 6, 2)

  a[1,1] <- -1/(p*(PR - betaT)) / N
  a[2,1] <- -1/(1 - p) / N
  a[3,1] <-  1/((1 - p)*r*lambda1*tau) / N
  a[4,1] <- -1/((1 - p)*r) / N
  a[6,1] <-  1/(OmegaT - betaT*Time)

  a[1,2] <- -2*(PR*(1-PR)+sigmabetaT^2*(1-p)) /
            p^2 / (PR - betaT)^3 / N^2 -
            2*sigmabetaT^2 * (OmegaT-PR*Time) /
            p / (PR-betaT)^3 / (OmegaT-betaT*Time) / N
  a[2,2] <- ((PR^2 + sigmabetaT^2*(1-2*p))
            /(p*(PR-betaT))^2 -
             (p^2+(1-p)^2)/(p*(1-p))^2) / N^2 +
            sigmabetaT^2 * (OmegaT-PR*Time) /
            p / (PR-betaT)^2 / (OmegaT-betaT*Time) / N
  a[3,2] <- -1 / ((1 - p)*r*lambda1*tau)^2 / N^2
  a[5,2] <-  (1 - 2*PR)/(p*(PR - betaT))^2 / N^2
  a[6,2] <- -2*sigmaOmegaT^2 / (OmegaT-betaT*Time)^2 -
            2*sigmabetaT^2*Time*(OmegaT-PR*Time) /
            (PR-betaT)/(OmegaT-betaT*Time)^3

  mat <- t(a) %*% V %*% a
  tA  <- log(R1) - log(R0)
  tB  <- ((PR*(1-PR) + sigmabetaT^2*(1-p))
         /(p*q*(PR-betaT)^2) +
         1/p + 1/(1-p) +
         1/((1-p)*r*lambda1*tau))/N +
        sigmaOmegaT^2/(OmegaT - betaT*Time)^2 +
        sigmabetaT^2*((OmegaT-PR*Time)/
                      (PR-betaT)/(OmegaT-betaT*Time))^2
  cvec <- c(1/sqrt(tB), -tA/(2*tB^(3/2)))
  drop(t(cvec) %*% mat %*% cvec)
}

################################################################################
## Wrappers for special variance propagation 
################################################################################

varZ_sigmabeta <- function(para, R1, R0 = 1, N) {
  p <- para$p; q <- para$q; OmegaT <- para$OmegaT; betaT <- para$betaT
  r <- para$r; tau <- para$tau; Time <- para$Time
  sigmabetaT <- para$sigmabetaT
  lambda0 <- para$lambda0; lambda1 <- lambda0 * R1
  PR <- betaT + lambda0*(1-p)/p*(OmegaT-betaT*Time)

  varAB <- var_AB_sigmabeta(p, q, PR, OmegaT, betaT, r, lambda1,
                            tau, sigmabetaT)
  tA <- log(R1) - log(R0)
  tB <- (PR*(1-PR)+sigmabetaT^2*(1-p))/(p*q*(PR-betaT)^2) +
      1/p + 1/(1-p) + 1/((1-p)*r*lambda1*tau)
  cvec <- c(1/sqrt(tB), -tA/(2*tB^(3/2)))
  drop(t(cvec) %*% varAB %*% cvec)
}

varZ_sigmaOmega <- function(para, R1, R0 = 1, N) {
  p <- para$p; OmegaT <- para$OmegaT; betaT <- para$betaT
  r <- para$r; tau <- para$tau; Time <- para$Time; q <- para$q 
  sigmaOmegaT <- para$sigmaOmegaT
  lambda0 <- para$lambda0; lambda1 <- lambda0 * R1
  PR <- betaT + lambda0*(1-p)/p*(OmegaT-betaT*Time)

  var_Z_sigmaOmega(p, q, PR, OmegaT, betaT, r, lambda1,
                   tau, sigmaOmegaT, Time, N, R1, R0)
}

varZ_both <- function(para, R1, R0 = 1, N) {
  p <- para$p; OmegaT <- para$OmegaT; betaT <- para$betaT
  r <- para$r; tau <- para$tau; Time <- para$Time; q <- para$q 
  sigmaOmegaT <- para$sigmaOmegaT; sigmabetaT <- para$sigmabetaT
  lambda0 <- para$lambda0; lambda1 <- lambda0 * R1
  PR <- betaT + lambda0*(1-p)/p*(OmegaT-betaT*Time)

  var_Z_both(p, q, PR, OmegaT, betaT, r, lambda1, tau,
             sigmaOmegaT, sigmabetaT, Time, N, R1, R0)
}

###############################################################################
## Run Scenarios
###############################################################################
run_scenario <- function(p, R1, lambda0, OmegaT, betaT, r, tau, alpha, power,
                         sigmaOmegaT, sigmabetaT, Time, outfile, R0 = 1) {
  ## Utility to recreate scenario tables

  q <- 1    # Assay coverage; set to 1 for published scenarios

  Results <- as.data.frame(matrix(0, nrow = length(lambda0), ncol = 16))
  names_vec <- c("lambda0","p","q","r","OmegaT","sigmaOmegaT","betaT",
                 "sigmabetaT","Time","tau","NToScreen","NHIVpos","NHIVneg",
                 "Power","NHIVenrol","R1")
  names(Results) <- names_vec
  for (i in seq_along(lambda0)) {
    para <- data.frame(lambda0 = lambda0[i], p = p, q = q, r = r,
                       OmegaT = OmegaT, sigmaOmegaT = sigmaOmegaT,
                       betaT = betaT, sigmabetaT = sigmabetaT,
                       Time = Time, tau = tau)
    NToScreen <- samplesize_CF(para, R1, R0, alpha, power)
    Results$NToScreen[i] <- NToScreen
    Results$NHIVpos[i]   <- p * NToScreen
    Results$NHIVneg[i]   <- (1 - p) * NToScreen
    Results$Power[i]     <- 100 * power
    Results$NHIVenrol[i] <- r * (1 - p) * NToScreen
    Results$R1[i]        <- R1
    Results[i, 1:10]     <- para
  }
  write.csv(Results, outfile, row.names = FALSE)
  Results
}

# --- Run scenarios --------------------
scenario1 <- run_scenario(p = 0.4, R1 = 0.3, lambda0 = c(0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10),
                          OmegaT = 163/365.25, betaT = 0.014, r = 0.9, tau = 1.5,
                          alpha = 0.05, power = 0.8,
                          sigmaOmegaT = (163/365.25)*0.083, sigmabetaT = 0.014*1.003,
                          Time = 2, 
                          outfile = "SampleSize_Results_prev40_eff70.csv",
                          R0 = 1)

scenario2 <- run_scenario(p = 0.3, R1 = 0.3, lambda0 = c(0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10),
                          OmegaT = 163/365.25, betaT = 0.014, r = 0.9, tau = 1.5,
                          alpha = 0.05, power = 0.8,
                          sigmaOmegaT = (163/365.25)*0.083, sigmabetaT = 0.014*1.003,
                          Time = 2,
                          outfile = "SampleSize_Results_prev30_eff70.csv",
                          R0 = 1)

scenario3 <- run_scenario(p = 0.4, R1 = 0.5, lambda0 = c(0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10),
                          OmegaT = 163/365.25, betaT = 0.014, r = 0.9, tau = 1.5,
                          alpha = 0.05, power = 0.8,
                          sigmaOmegaT = (163/365.25)*0.083, sigmabetaT = 0.014*1.003,
                          Time = 2,
                          outfile = "SampleSize_Results_prev40_eff50.csv",
                          R0 = 1)

scenario4 <- run_scenario(p = 0.3, R1 = 0.5, lambda0 = c(0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10),
                          OmegaT = 163/365.25, betaT = 0.014, r = 0.9, tau = 1.5,
                          alpha = 0.05, power = 0.8,
                          sigmaOmegaT = (163/365.25)*0.083, sigmabetaT = 0.014*1.003,
                          Time = 2,
                          outfile = "SampleSize_Results_prev30_eff50.csv",
                          R0 = 1)
