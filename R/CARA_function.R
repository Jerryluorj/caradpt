###############################################################################
####################   Covariate Adjusted Doubly Biased Coin Design   #########
###############################################################################
#' Allocation Function of Covariate Adjusted Doubly Biased Coin Design for Binary and Continuous Response
#' @importFrom stats binomial coef glm lm pnorm predict qchisq rnorm runif vcov rexp
#' @param ptsb.cov a \code{n x k} covariate matrix of previous patients.
#' @param ptsb.t a treatment vector of previous patients with length \code{n}.
#' @param ptsb.Y a response vector of previous patients with length \code{n}.
#' @param ptnow.cov a covariate vector of the incoming patient with length \code{k}.
#' @param v a non-negative integer that controls the randomness of CADBCD design. The default value is 2.
#' @param response the type of the response. Options are \code{"Binary"} or \code{"Cont"}.
#' @param target the type of optimal allocation target. Options are \code{"Neyman"} or \code{"RSIHR"}.
#'
#' @description Calculating the probability of assigning the upcoming patient to treatment A based on the patient's covariates and the previous patients' covariates and responses for Covariate Adjusted Doubly Biased Coin procedure proposed by Zhang and Hu.
#' @return \item{prob}{Probability of assigning the upcoming patient to treatment A.}
#' @details
#' To start, we let \eqn{\boldsymbol{\theta}_0} be an initial estimate of \eqn{\boldsymbol{\theta}}, and assign \eqn{m_0} subjects to each treatment using a restricted randomization.
#' Assume that \eqn{m} (with \eqn{m \geq 2 m_0}) subjects have been assigned to treatments. Their responses \eqn{\{\boldsymbol{Y}_j, j = 1, \ldots, m\}} and the corresponding covariates \eqn{\{\boldsymbol{\xi}_j, j = 1, \ldots, m\}} are observed.
#' We let \eqn{\widehat{\boldsymbol{\theta}}_m = (\widehat{\boldsymbol{\theta}}_{m, 1}, \widehat{\boldsymbol{\theta}}_{m, 2})} be an estimate of \eqn{\boldsymbol{\theta} = (\boldsymbol{\theta}_1, \boldsymbol{\theta}_2)}.
#' For each \eqn{k = 1, 2}, the estimator \eqn{\widehat{\boldsymbol{\theta}}_{m, k} = \widehat{\boldsymbol{\theta}}_{m, k}(Y_{j,k}, \boldsymbol{\xi}_j : X_{j,k} = 1, j = 1, \ldots, m)} is based on the observed sample of size \eqn{N_{m,k}}, that is, \eqn{\{(Y_{j,k}, \boldsymbol{\xi}_j) : X_{j,k} = 1, j = 1, \ldots, m\}}.
#'
#' Define \eqn{\widehat{\rho}_m = \frac{1}{m} \sum_{i=1}^m \pi_1(\widehat{\boldsymbol{\theta}}_m, \boldsymbol{\xi}_i)} and \eqn{\widehat{\pi}_m = \pi_1(\widehat{\boldsymbol{\theta}}_m, \boldsymbol{\xi}_{m+1})}.
#' When the \eqn{(m+1)}-th subject is ready for randomization and the corresponding covariate \eqn{\boldsymbol{\xi}_{m+1}} is recorded, we assign the patient to treatment 1 with the probability:
#'
#' \deqn{
#' \psi_{m+1,1} = \frac{\widehat{\pi}_m \left( \frac{\widehat{\rho}_m}{N_{m,1}/m} \right)^v}
#' {\widehat{\pi}_m \left( \frac{\widehat{\rho}_m}{N_{m,1}/m} \right)^v+ (1 - \widehat{\pi}_m) \left( \frac{1 - \widehat{\rho}_m}{1 - N_{m,1}/m} \right)^v}
#' }
#'
#' and to treatment 2 with probability \eqn{\psi_{m+1,2} = 1 - \psi_{m+1,1}}, where \eqn{v \geq 0} is a constant controlling the degree of randomness—from the most random when \eqn{v = 0} to the most deterministic when \eqn{v \rightarrow \infty}. See Zhang and Hu(2009) for more details.
#' @references
#' Zhang, L. X., Hu, F., Cheung, S. H., & Chan, W. S. (2007). Asymptotic properties of covariate-adjusted response-adaptive designs.
#' \emph{Annals of Statistics}, 35(3), 1166–1182.
#'
#' Zhang, L. X., & Hu, F. F. (2009). A new family of covariate-adjusted response adaptive designs and their properties.
#' \emph{Applied Mathematics—A Journal of Chinese Universities}, 24(1), 1–13.
#' @concept CADBCD Design
#' @examples
#' set.seed(123)
#'
#' n_prev = 40
#' covariates = cbind(Z1 = rnorm(n_prev), Z2 = rnorm(n_prev))
#' treatment = sample(c(0, 1), n_prev, replace = TRUE)
#' response = rbinom(n_prev, size = 1, prob = 0.6)
#'
#' # Simulate new incoming patient
#' new_patient_cov = c(Z1 = rnorm(1), Z2 = rnorm(1))
#'
#' # Run allocation function
#' result = CADBCD_Alloc(
#'   ptsb.cov = covariates,
#'   ptsb.t = treatment,
#'   ptsb.Y = response,
#'   ptnow.cov = new_patient_cov,
#'   response = "Binary",
#'   target = "Neyman"
#' )
#' print(result$prob)

#' @export
CADBCD_Alloc = function(ptsb.cov,
                        ptsb.t,
                        ptsb.Y,
                        ptnow.cov,
                        v = 2,
                        response,
                        target) {
  # --- Step 1: Input validation ---
  if (length(ptsb.t) != length(ptsb.Y) ||
      nrow(ptsb.cov) != length(ptsb.Y)
      || length(ptsb.t) != nrow(ptsb.cov)) {
    stop(
      "Covariate matrix, treatment vector, and response vector must have the same number of observations."
    )
  }
  if (!is.matrix(ptsb.cov) && !is.data.frame(ptsb.cov)) stop("ptsb.cov must be a matrix or data.frame.")
  if (!is.vector(ptnow.cov)) stop("ptnow.cov must be a numeric vector.")
  if (!all(response %in% c("Binary", "Cont"))) stop("response must be 'Binary' or 'Cont'.")
  if (!all(target %in% c("Neyman", "RSIHR"))) stop("target must be 'Neyman' or 'RSIHR'.")

  # --- Step 2: Model fitting for group A and B ---
  ptsb = data.frame(ptsb.cov, Treat = ptsb.t, Y = ptsb.Y)
  ptsb.A = ptsb[ptsb$Treat == 1, ]
  ptsb.B = ptsb[ptsb$Treat == 0, ]

  if (response == "Binary") {
    fitA = glm(Y ~ ., data = ptsb.A[, !(names(ptsb.A) %in% "Treat")], family = binomial)
    fitB = glm(Y ~ ., data = ptsb.B[, !(names(ptsb.B) %in% "Treat")], family = binomial)
  }
  else if (response == "Cont") {
    fitA = lm(Y ~ ., data = ptsb.A[, !(names(ptsb.A) %in% "Treat")])
    fitB = lm(Y ~ ., data = ptsb.B[, !(names(ptsb.B) %in% "Treat")])
  }

  # --- Step 3: Estimate individual and average optimal allocation probabilities ---
  pi.m = pi.patient(
    pt.cov = ptnow.cov,
    fitA = fitA,
    fitB = fitB,
    response = response,
    target = target
  )
  rho.m = mean(
    apply(
      ptsb.cov,
      1,
      pi.patient,
      fitA = fitA,
      fitB = fitB
      ,
      response = response,
      target = target
    )
  )

# --- Step 4: Compute allocation probability via CADBCD ---
  prob = CADBCD.prob(
    pi.m = pi.m,
    rho.m = rho.m,
    m = nrow(ptsb.cov) + 1,
    v = v,
    pts = ptsb
  )

  return(list(prob = prob))
}


#' Allocation Function of Covariate Adjusted Doubly Biased Coin Design for Survival Response
#'
#' @importFrom stats binomial coef glm lm pnorm predict qchisq rnorm runif rexp vcov as.formula
#' @importFrom survival Surv coxph survreg
#' @param ptsb.cov a \code{n x k} covariate matrix of previous patients.
#' @param ptsb.t a treatment vector of previous patients with length \code{n}.
#' @param ptsb.Y a response vector of previous patients with length \code{n}.
#' @param ptsb.E a censoring indicator vector (1 = event observed, 0 = censored)with length \code{n}.
#' @param ptnow.cov a covariate vector of the incoming patient with length \code{k}.
#' @param v a non-negative integer that controls the randomness of CADBCD design. The default value is 2.
#' @param target the type of optimal allocation target. Options are \code{"Neyman"} or \code{"RSIHR"}.
#' @description Calculating the probability of assigning the upcoming patient to treatment A based on the patient's covariates and the previous patients' covariates and responses
#' for Covariate Adjusted Doubly Biased Coin procedure for survival trial.
#' @details
#' For the first \eqn{m = 2m_0} patients, a restricted randomization procedure is used to allocate them equally to treatments A and B.
#'
#' After the initial enrollment of \eqn{2m_0} patients, adaptive treatment assignment begins from patient \eqn{2m_0 + 1} onward. Importantly, this adaptive assignment does \emph{not} require complete observation of outcomes from earlier patients. Instead, at each new patient arrival, the treatment effect parameters \eqn{(\tilde{\boldsymbol{\beta}}_{A,m}, \tilde{\boldsymbol{\beta}}_{B,m})} are re-estimated via partial likelihood using the accrued survival data subject to \emph{dynamic administrative censoring}. Specifically, each previously enrolled patient is censored at the current arrival time if their event has not yet occurred. These updated estimates are then used to guide covariate-adjusted allocation.
#'
#' When the \eqn{(m+1)}-th patient enters the trial with covariate vector \eqn{\boldsymbol{Z}_{m+1}}, the probability of assigning this patient to treatment A is computed as:
#' \deqn{
#' \phi_{m+1} = g(\tilde{\boldsymbol{\beta}}_{A,m}, \tilde{\boldsymbol{\beta}}_{B,m}, \boldsymbol{Z}_{m+1}),
#' }
#' where \eqn{0 \leq g(\cdot) \leq 1} is an appropriately chosen allocation function that favors the better-performing treatment arm.
#'
#' Following Zhang and Hu's CADBCD procedure, let \eqn{N_A(m)} and \eqn{N_B(m) = m - N_A(m)} denote the numbers of patients allocated to treatments A and B after \eqn{m} assignments. Let \eqn{\tilde{\pi}_m = \pi_A(\tilde{\boldsymbol{\beta}}_{A,m}, \tilde{\boldsymbol{\beta}}_{B,m}, \boldsymbol{Z}_{m+1})} be the target allocation probability for the incoming patient. The average target allocation proportion among the first \eqn{m} patients is defined as:
#' \deqn{
#' \tilde{\rho}_m = \frac{1}{m} \sum_{i=1}^m \pi_A(\tilde{\boldsymbol{\beta}}_{A,m}, \tilde{\boldsymbol{\beta}}_{B,m}, \boldsymbol{Z}_i).
#' }
#'
#' Under the CADBCD scheme, the probability of assigning treatment A to the \eqn{(m+1)}-th patient is given by:
#' \deqn{
#' \phi_{m+1} = \frac{\tilde{\pi}_m \left( \frac{\tilde{\rho}_m}{N_A(m)/m} \right)^v}
#' {\tilde{\pi}_m \left( \frac{\tilde{\rho}_m}{N_A(m)/m} \right)^v + (1 - \tilde{\pi}_m) \left( \frac{1 - \tilde{\rho}_m}{N_B(m)/m} \right)^v},
#' }
#' and to treatment 2 with probability \eqn{\psi_{m+1,2} = 1 - \psi_{m+1,1}}, where \eqn{v \geq 0} is a constant controlling the degree of randomness—from the most random when \eqn{v = 0} to the most deterministic when \eqn{v \rightarrow \infty}.
#' similiar procedure for survival responses are used in all other designs.
#' @returns \item{prob}{Probability of assigning the upcoming patient to treatment A.}
#' @references
#' Mukherjee, A., Jana, S., & Coad, S. (2024). Covariate-adjusted response-adaptive designs for semiparametric survival models.
#' \emph{Statistical Methods in Medical Research}, 09622802241287704.
#' @export
#' @concept CADBCD Design
#' @examples set.seed(123)
#'n = 40
#'covariates = cbind(rexp(40),rexp(40))
#'treatment = sample(c(0, 1), n, replace = TRUE)
#'survival_time = rexp(n, rate = 1)
#'censoring = runif(n)
#'event = as.numeric(survival_time < censoring)

# Simulated new patient
#'new_patient_cov = c(Z1 = 1, Z2 = 0.5)

#'result = CADBCD_Alloc_Surv(
#'  ptsb.cov = covariates,
#'  ptsb.t = treatment,
#'  ptsb.Y = survival_time,
#'  ptsb.E = event,
#'  ptnow.cov = new_patient_cov,
#'  v = 2,
#'  target = "Neyman"
#')
#'print(result$prob)
CADBCD_Alloc_Surv = function(ptsb.cov,
                             ptsb.t,
                             ptsb.Y,
                             ptsb.E,
                             ptnow.cov,
                             v = 2,
                             target) {
  # if (nrow(ptsb.cov)<40){
  #   warning("Sample size is too small.Estimations may be biased.")
  # }
  if (length(ptsb.t) != length(ptsb.Y) ||
      nrow(ptsb.cov) != length(ptsb.Y)
      || length(ptsb.t) != nrow(ptsb.cov)) {
    stop(
      "Covariate matrix, treatment vector, and response vector must have the same number of observations."
    )
  }
  colnames(ptsb.cov) = paste0("Z", 1:ncol(ptsb.cov))
  ptsb = data.frame(ptsb.cov,
                    Treat = ptsb.t,
                    Y = ptsb.Y,
                    E = ptsb.E)
  ptsb.A = ptsb[ptsb$Treat == 1, ]
  ptsb.B = ptsb[ptsb$Treat == 0, ]

  # if (model == "Cox") {
    formula_str = paste("Surv(Y, E) ~", paste(colnames(ptsb.cov), collapse = " + "))
    fitA = coxph(as.formula(formula_str), data = ptsb.A)
    fitB = coxph(as.formula(formula_str), data = ptsb.B)
#
#   }

  pi.m = pi.patient.surv(
    pt.cov = ptnow.cov,
    fitA = fitA,
    fitB = fitB,
    target = target
  )
  rho.m = mean(
    apply(
      ptsb.cov,
      1,
      pi.patient.surv,
      fitA = fitA,
      fitB = fitB,
      target = target
    )
  )

  prob = CADBCD.prob(
    pi.m = pi.m,
    rho.m = rho.m,
    m = nrow(ptsb.cov) + 1,
    v = v,
    pts = ptsb
  )

  return(list(prob = prob))
}



#' Simulation Function of Covariate Adjusted Doubly Biased Coin Design for Binary and Continuous Response
#' @importFrom stats binomial coef glm lm pnorm predict qchisq rnorm runif vcov rexp
#' @param n a positive integer. The value specifies the total number of participants involved in each round of the simulation.
#' @param thetaA a vector of length \code{k+1}. The true coefficient parameter value for treatment A.
#' @param thetaB a vector of length \code{k+1}. The true coefficient parameter value for treatment B.
#' @param m0 a positive integer. The number of first 2m0 patients will be allocated equally for estimation. The default value is 40.
#' @param pts.cov a \code{n x k} matrix. The simulated covariate matrix for patients.
#' @param v a non-negative integer that controls the randomness of CADBCD design. The default value is 2.
#' @param response the type of the response. Options are \code{"Binary"} or \code{"Cont"}.
#' @param target the type of optimal allocation target. Options are \code{"Neyman"} or \code{"RSIHR"}.
#' @concept CADBCD Design
#' @description
#' This function simulates a clinical trial using the Covariate Adjusted Doubly Biased Coin Design (CADBCD) with Binary or Continuous Responses.
#' @return A list with the following elements:
#' \item{method}{The name of the procedure.}
#' \item{sampleSize}{Total number of patients.}
#' \item{parameter}{Estimated parameter values.}
#' \item{assignment}{Treatment assignment vector.}
#' \item{proportion}{Proportion of patients allocated to treatment A.}
#' \item{responses}{Simulated response values.}
#' \item{failureRate}{Proportion of treatment failures (if \code{response = "Binary"}).}
#' \item{meanResponse}{Mean response value (if \code{response = "Cont"}).}
#' \item{rejectNull}{Logical. Indicates whether the treatment effect is statistically significant based on a Wald test.}
#'
#' @examples
#' set.seed(123)
#' results = CADBCD_Sim(n = 400,
#'                       pts.cov = cbind(rnorm(400), rnorm(400)),
#'                       thetaA = c(-1, 1, 1),
#'                       thetaB = c(3, 1, 1),
#'                       response = "Binary",
#'                       target = "Neyman")
## view the outputs
#' results
#'
#' ## view the settings
#' results$method
#' results$sampleSize
#'
#' ## view the simulation results
#' results$parameter
#' results$assignments
#' results$proportion
#' results$responses
#' results$failureRate
#'
#' @export
CADBCD_Sim = function(n,
                      thetaA,
                      thetaB,
                      m0 = 40,
                      pts.cov,
                      v = 2,
                      response,
                      target) {
  if (length(thetaA) != ncol(pts.cov) + 1 ||
      length(thetaB) != ncol(pts.cov) + 1) {
    stop("Length of true parameter vector must be equal to the number of covariates+1.")
  }

  # Check sample size
  if (!is.numeric(n) || length(n) != 1 || n <= 0 || n != floor(n)) {
    stop("`n` must be a positive integer.")
  }

  # Check m0
  if (!is.numeric(m0) || length(m0) != 1 || m0 <= 0 || m0 != floor(m0)) {
    stop("`m0` must be a positive integer.")
  }

  # Check pts.cov
  if (!is.matrix(pts.cov) || nrow(pts.cov) != n) {
    stop("`pts.cov` must be a numeric matrix with `n` rows.")
  }

  k = ncol(pts.cov)

  # Check thetaA and thetaB
  if (length(thetaA) != k + 1) {
    stop("`thetaA` must be of length k + 1, where k is the number of covariates.")
  }
  if (length(thetaB) != k + 1) {
    stop("`thetaB` must be of length k + 1, where k is the number of covariates.")
  }

  # Check v
  if (!is.numeric(v) || length(v) != 1 || v < 0 ) {
    stop("`v` must be a non-negative value.")
  }

  # Check response
  if (!response %in% c("Binary", "Cont")) {
    stop("`response` must be either 'Binary' or 'Cont'.")
  }

  # Check target
  if (!target %in% c("Neyman", "RSIHR")) {
    stop("`target` must be either 'Neyman' or 'RSIHR'.")
  }

  pts = CADBCD_generate_data(
    n = n,
    pts.cov = pts.cov,
    m0 = m0,
    thetaA = thetaA,
    thetaB = thetaB,
    response = response
  )
  for (m in ((2 * m0 + 1):nrow(pts))) {
    ptsb = pts[1:(m - 1), ]

    ptsb.A = data.frame(ptsb[ptsb[, "Treat"] == 1, ])
    ptsb.B = data.frame(ptsb[ptsb[, "Treat"] == 0, ])

    if (response == "Binary") {
      fitA = glm(Y ~ ., data = ptsb.A[, !(names(ptsb.A) %in% "Treat")], family = binomial)
      fitB = glm(Y ~ ., data = ptsb.B[, !(names(ptsb.B) %in% "Treat")], family = binomial)
    }
    else if (response == "Cont") {
      fitA = lm(Y ~ ., data = ptsb.A[, !(names(ptsb.A) %in% "Treat")])
      fitB = lm(Y ~ ., data = ptsb.B[, !(names(ptsb.B) %in% "Treat")])
    }
    pi.m = pi.patient(
      pt.cov = pts[m, 1:ncol(pts.cov)],
      fitA = fitA,
      fitB = fitB,
      response = response,
      target = target
    )
    rho.m = mean(
      apply(
        ptsb[, 1:ncol(pts.cov)],
        1,
        pi.patient,
        fitA = fitA,
        fitB = fitB
        ,
        response = response,
        target = target
      )
    )

    prob = CADBCD.prob(
      pi.m = pi.m,
      rho.m = rho.m,
      m = m,
      v = v,
      pts = pts
    )

    if (anyNA(prob)) {
      pts[m, "Treat"] = sample(c(1, 0), 1, replace = TRUE, c(0.5, 0.5))
    }
    else{
      pts[m, "Treat"] = sample(c(1, 0), 1, replace = TRUE, c(prob, 1 - prob))
    }

    z_m = c(1, pts[m, 1:ncol(pts.cov)])
    if (response == "Binary") {
      exp_m = exp(-(z_m %*% thetaA * pts[m, "Treat"] + z_m %*% thetaB * (1 - pts[m, "Treat"])))
      pts[m, "Y"] = as.numeric(runif(1) < 1 / (1 + exp_m))
    }
    else if (response == "Cont") {
      pts[m, "Y"] = z_m %*% thetaA * pts[m, "Treat"] + z_m %*% thetaB * (1 - pts[m, "Treat"]) +
        rnorm(1)
    }
  }

  return(
    CARA_Output(
      name = "Covariate Adjusted Doubly Biased Coin Design",
      pts = pts,
      response = response
    )
  )
}


#' Simulation For Covariate Adjusted Doubly Biased Coin Design for Survival Response
#' @description
#' This function simulates a clinical trial with time-to-event (survival) outcomes using
#' the Covariate Adjusted Doubly Biased Coin Design (CADBCD).
#' Patient responses are generated under the Cox proportional hazards model,
#' assuming the proportional hazards assumption holds.
#' @importFrom stats binomial coef glm lm pnorm predict qchisq rnorm runif vcov rexp
#' @importFrom survival Surv coxph survreg
#' @param n a positive integer. The value specifies the total number of participants involved in each round of the simulation.
#' @param thetaA a vector of length \code{k}. The true coefficient parameter value for treatment A.
#' @param thetaB a vector of length \code{k}. The true coefficient parameter value for treatment B.
#' @param m0 a positive integer. The number of first 2m0 patients will be allocated equally for estimation. The default value is 40.
#' @param pts.cov a \code{n x k} matrix. The simulated covariate matrix for patients.
#' @param v a non-negative integer that controls the randomness of CADBCD design. The default value is 2.
#' @param target the type of optimal allocation target. Options are \code{"Neyman"} or \code{"RSIHR"}.
#' @param censor.time a positive value. The upper bound to the simulated uniform censor time.
#' @param arrival.rate a positive value. The rate of simulated exponential arrival time.
#' @return A list with the following elements:
#' \item{method}{The name of procedure.}
#' \item{sampleSize}{Sample size of the trial.}
#' \item{parameter}{Estimated parameters used to do the simulations.}
#' \item{N.events}{Total number of events of the trial.}
#' \item{assignment}{The randomization sequence.}
#' \item{proportion}{Average allocation proportion for treatment A.}
#' \item{responses}{The simulated observed survival responses of patients.}
#' \item{events}{Whether events are observed for patients(1=event,0=censored).}
#' \item{rejectNull}{Whether the study to detect a significant difference of treatment effect using Wald test.}
#' @concept CADBCD Design
#' @export
#' @examples
#' set.seed(123)
#'
#' ## Run CADBCD simulation with survival response
#' results = CADBCD_Sim_Surv(
#'   thetaA = c(0.1, 0.1),
#'   thetaB = c(-1, 0.1),
#'   n = 400,
#'   pts.cov = cbind(sample(c(1, 0), 400, replace = TRUE), rnorm(400)),
#'   target = "RSIHR",
#'   censor.time = 2,
#'   arrival.rate = 150
#' )

CADBCD_Sim_Surv = function(n,
                           thetaA,
                           thetaB,
                           m0 = 40,
                           pts.cov,
                           v = 2,
                           target,
                           censor.time,
                           arrival.rate) {
  if (!is.numeric(n) || length(n) != 1 || n <= 0 || n != floor(n)) {
    stop("`n` must be a single positive integer.")
  }

  if (!is.numeric(m0) || length(m0) != 1 || m0 <= 0 || m0 != floor(m0)) {
    stop("`m0` must be a single positive integer.")
  }

  if (!is.numeric(v) || length(v) != 1 || v < 0 ) {
    stop("`v` must be a non-negative number.")
  }

  if (!is.matrix(pts.cov)) {
    stop("`pts.cov` must be a numeric matrix.")
  }

  if (nrow(pts.cov) != n) {
    stop("`pts.cov` must have exactly `n` rows.")
  }

  k = ncol(pts.cov)

  if (length(thetaA) != k) {
    stop("`thetaA` must be a numeric vector of length k, where k is the number of covariates.")
  }

  if (length(thetaB) != k) {
    stop("`thetaB` must be a numeric vector of length k, where k is the number of covariates.")
  }

  if (!is.character(target) || !(target %in% c("Neyman", "RSIHR"))) {
    stop("`target` must be either 'Neyman' or 'RSIHR'.")
  }

  if (!is.numeric(censor.time) || length(censor.time) != 1 || censor.time <= 0) {
    stop("`censor.time` must be a single positive number.")
  }

  if (!is.numeric(arrival.rate) || length(arrival.rate) != 1 || arrival.rate <= 0) {
    stop("`arrival.rate` must be a single positive number.")
  }
  pts = CADBCD_generate_data_surv(
    n = n,
    pts.cov = pts.cov,
    m0 = m0,
    thetaA = thetaA,
    thetaB = thetaB,
    censor.time = censor.time,
    arrival.rate = arrival.rate
  )
  colnames(pts.cov) = paste0("Z", 1:ncol(pts.cov))
  for (m in ((2 * m0 + 1):nrow(pts))) {
    current_time = pts[m, "A"]
    ptsb = pts[1:(m - 1), ]
    ptsb.A = ptsb[ptsb[, "Treat"] == 1, ]
    ptsb.B = ptsb[ptsb[, "Treat"] == 0, ]

    Ai_A = current_time - ptsb.A[, "A"]
    Ai_B = current_time - ptsb.B[, "A"]

    E_star_A = as.numeric(ptsb.A[, "Y"] < Ai_A)
    E_star_B = as.numeric(ptsb.B[, "Y"] < Ai_B)

    ptsb.A_star = cbind(ptsb.A, E_star_A)
    ptsb.B_star = cbind(ptsb.B, E_star_B)

    ptsb.A_star[, "E_star_A"][ptsb.A_star[, "Y"] == ptsb.A_star[, "C"]] = 0
    ptsb.B_star[, "E_star_B"][ptsb.B_star[, "Y"] == ptsb.B_star[, "C"]] = 0

    Y_star_A = pmin(ptsb.A_star[, "Y"], Ai_A)
    Y_star_B = pmin(ptsb.B_star[, "Y"], Ai_B)

    ptsb.A_star = cbind(ptsb.A_star, Y_star_A)
    ptsb.B_star = cbind(ptsb.B_star, Y_star_B)

    formula_strA = paste("Surv(Y_star_A, E_star_A) ~",
                         paste(colnames(pts.cov), collapse = " + "))
    formula_strB = paste("Surv(Y_star_B, E_star_B) ~",
                         paste(colnames(pts.cov), collapse = " + "))
    fitA = coxph(as.formula(formula_strA), data = data.frame(ptsb.A_star))
    fitB = coxph(as.formula(formula_strB), data = data.frame(ptsb.B_star))


    pi.m = pi.patient.surv(
      pt.cov = pts[m, 1:ncol(pts.cov)],
      fitA = fitA,
      fitB = fitB,
      target = target
    )
    rho.m = mean(
      apply(
        ptsb[, 1:ncol(pts.cov)],
        1,
        pi.patient.surv,
        fitA = fitA,
        fitB = fitB,
        target = target
      )
    )

    prob = CADBCD.prob(
      pi.m = pi.m,
      rho.m = rho.m,
      m = m,
      v = v,
      pts = ptsb
    )
    if (anyNA(prob)) {
      pts[m, "Treat"] = sample(c(1, 0), 1, replace = TRUE, c(0.5, 0.5))
    }
    else{
      pts[m, "Treat"] = sample(c(1, 0), 1, replace = TRUE, c(prob, 1 - prob))
    }

    z_m = pts[m, 1:ncol(pts.cov)]
    linpred = z_m %*% thetaA * pts[m, "Treat"] + z_m %*% thetaB * (1 - pts[m, "Treat"])
    pts[m, "S"] = rexp(1, rate = 1 / exp(linpred))
    pts[m, "E"] = as.numeric(pts[m, "S"] < pts[m, "C"])
    pts[m, "Y"] = pmin(pts[m, "S"], pts[m, "C"])
  }

  return(CARA_Output_Surv(name = "Covariate Adjusted Doubly Biased Coin Design", pts =
                            pts))
}

###############################################################################
####################   CARA Designs Based on Efficiency and Ethics   ##########
###############################################################################
#' Allocation Function of CARA Designs Based on Efficiency and Ethics for Binary and Continuous Response.
#' @importFrom stats binomial coef glm lm pnorm predict qchisq rnorm runif vcov rexp
#' @param ptsb.cov a \code{n x k} covariate matrix of previous patients.
#' @param ptsb.t a treatment vector of previous patients with length \code{n}.
#' @param ptsb.Y a response vector of previous patients with length \code{n}.
#' @param ptnow.cov a covariate vector of the incoming patient with length \code{k}.
#' @param response the type of the response. Options are \code{"Binary"} or \code{"Cont"}.
#' @param gamma a non-negative number. A tuning parameter that reflects the importance of the efficiency component compared to the ethics component.
#' @concept CARAEE Design
#' @description Calculating the probability of assigning the upcoming patient to treatment A based on the patient's covariates and the previous patients' covariates and responses for CARAEE procedure.
#' @return \item{prob}{Probability of assigning the upcoming patient to treatment A.}
#' @details
#' Covariate-Adjusted Response-Adaptive with Ethics and Efficiency (CARAEE) Design:
#' The CARAEE procedure balances both ethical considerations and statistical efficiency when assigning subjects to treatments.
#'
#' As a start-up rule, \eqn{m_0} subjects are assigned to each treatment using a balanced randomization scheme.
#'
#' Assume that \eqn{m \geq 2m_0} subjects have been assigned, and their responses \eqn{\{\boldsymbol{X}_i, i = 1, \ldots, m\}} and covariates \eqn{\{\boldsymbol{Z}_i, i = 1, \ldots, m\}} are observed. Let \eqn{\hat{\boldsymbol{\theta}}(m) = \left( \hat{\theta}_1(m), \hat{\theta}_2(m) \right)}, where \eqn{\hat{\theta}_k(m)} is the maximum likelihood estimate of the treatment-specific parameter \eqn{\theta_k} based on the data for treatment group \eqn{k}.
#'
#' For the incoming subject \eqn{(m+1)} with covariates \eqn{\boldsymbol{Z}_{m+1}}, we define the efficiency and ethics measures for each treatment as:
#' \eqn{
#' \boldsymbol{d}(\boldsymbol{Z}, \boldsymbol{\theta}) = \left(d_1(\boldsymbol{Z}, \theta), d_2(\boldsymbol{Z}, \theta)\right), \quad
#' \boldsymbol{e}(\boldsymbol{Z}, \boldsymbol{\theta}) = \left(e_1(\boldsymbol{Z}, \theta), e_2(\boldsymbol{Z}, \theta)\right).
#' }
#'
#' The allocation probability of assigning subject \eqn{(m+1)} to treatment 1 is given by:
#' \deqn{
#' \phi_{m+1}(\boldsymbol{Z}_{m+1}, \hat{\boldsymbol{\theta}}(m)) =
#' \frac{e_1(\boldsymbol{Z}_{m+1}, \hat{\boldsymbol{\theta}}(m)) \cdot d_1^\gamma(\boldsymbol{Z}_{m+1}, \hat{\boldsymbol{\theta}}(m))}
#' {e_1(\boldsymbol{Z}_{m+1}, \hat{\boldsymbol{\theta}}(m)) \cdot d_1^\gamma(\boldsymbol{Z}_{m+1}, \hat{\boldsymbol{\theta}}(m)) + e_2(\boldsymbol{Z}_{m+1}, \hat{\boldsymbol{\theta}}(m)) \cdot d_2^\gamma(\boldsymbol{Z}_{m+1}, \hat{\boldsymbol{\theta}}(m))}.
#' }
#'
#' This allocation rule is scale-invariant in both efficiency and ethics components due to its ratio-based form. The tuning parameter \eqn{\gamma \geq 0} controls the trade-off between the two: when \eqn{\gamma = 0}, the assignment is based purely on ethical considerations; larger values of \eqn{\gamma} increase the emphasis on statistical efficiency. More details can be found in Hu, Zhu & Zhang(2015).
#' @references
#' Hu, J., Zhu, H., & Hu, F. (2015). A unified family of covariate-adjusted response-adaptive designs based on efficiency and ethics.
#' \emph{Journal of the American Statistical Association}, 110(509), 357–367.
#' @export
#'
#' @examples
#' set.seed(123)
#'
#'n_prev = 40
#'covariates = cbind(Z1 = rnorm(n_prev), Z2 = rnorm(n_prev))
#'treatment = sample(c(0, 1), n_prev, replace = TRUE)
#'response = rbinom(n_prev, size = 1, prob = 0.6)
#'
#'# Simulate new incoming patient
#'new_patient_cov = c(Z1 = rnorm(1), Z2 = rnorm(1))
#'
#'# Run allocation function
#'result = CARAEE_Alloc(
#'  ptsb.cov = covariates,
#'  ptsb.t = treatment,
#'  ptsb.Y = response,
#'  ptnow.cov = new_patient_cov,
#'  response = "Binary",
#'  gamma=1
#')
#'print(result$prob)
#'
CARAEE_Alloc = function(ptsb.cov,
                        ptsb.t,
                        ptsb.Y,
                        ptnow.cov,
                        gamma,
                        response) {
  # --- Step 1: Input validation ---
  if (length(ptsb.t) != length(ptsb.Y) ||
      nrow(ptsb.cov) != length(ptsb.Y)
      || length(ptsb.t) != nrow(ptsb.cov)) {
    stop(
      "Covariate matrix, treatment vector, and response vector must have the same number of observations."
    )
  }
  if (!is.matrix(ptsb.cov) && !is.data.frame(ptsb.cov)) stop("ptsb.cov must be a matrix or data.frame.")
  # if (!is.vector(ptnow.cov)) stop("ptnow.cov must be a numeric vector.")
  if (!all(response %in% c("Binary", "Cont"))) stop("response must be 'Binary' or 'Cont'.")
  if (!is.numeric(gamma) || length(gamma) != 1 || gamma < 0) {
    stop("`gamma` must be a non-negative number.")
  }
  # --- Step 2: Model fitting for group A and B ---
  ptsb = data.frame(ptsb.cov, Treat = ptsb.t, Y = ptsb.Y)
  ptsb.A = ptsb[ptsb$Treat == 1, ]
  ptsb.B = ptsb[ptsb$Treat == 0, ]

  ptnow.cov=as.data.frame(t(ptnow.cov))
  colnames(ptnow.cov)=colnames(ptsb.cov)

  if (response == "Binary") {
    fitA = glm(Y ~ ., data = ptsb.A[, !(names(ptsb.A) %in% "Treat")], family = binomial)
    fitB = glm(Y ~ ., data = ptsb.B[, !(names(ptsb.B) %in% "Treat")], family = binomial)


    pi1 = predict(fitA, newdata = ptnow.cov, type = "response")
    pi2 = predict(fitB, newdata = ptnow.cov, type = "response")

    # Ethics and Efficiency
    e1 = pi1
    e2 = pi2
    d1 = 1 / (pi1 * (1 - pi1))
    d2 = 1 / (pi2 * (1 - pi2))
  }
  else if (response == "Cont") {
    fitA = lm(Y ~ ., data = ptsb.A[, !(names(ptsb.A) %in% "Treat")])
    fitB = lm(Y ~ ., data = ptsb.B[, !(names(ptsb.B) %in% "Treat")])

    ptnow.cov=as.data.frame(ptnow.cov)

    e1 = predict(fitA, newdata = ptnow.cov)
    e2 = predict(fitB, newdata = ptnow.cov)

    # Efficiency component
    X1 = model.matrix(fitA)
    X2 = model.matrix(fitB)
    sigma1_sq = summary(fitA)$sigma^2
    sigma2_sq = summary(fitB)$sigma^2

    invXtX1 = solve(t(X1) %*% X1)
    invXtX2 = solve(t(X2) %*% X2)

    z_vec = c(1, as.numeric(ptnow.cov))  # intercept + covariates
    d1 = sigma1_sq * t(z_vec) %*% invXtX1 %*% z_vec
    d2 = sigma2_sq * t(z_vec) %*% invXtX2 %*% z_vec
  }

  numerator = e1 * d1^gamma
  denominator = numerator + e2 * d2^gamma
  prob = numerator / denominator

  return(list(prob = prob))
}


#' Allocation Function of CARA Designs Based on Efficiency and Ethics for Survival Response
#'
#' @importFrom stats binomial coef glm lm pnorm predict qchisq rnorm runif rexp vcov as.formula
#' @importFrom survival Surv coxph survreg
#' @param ptsb.cov a \code{n x k} covariate matrix of previous patients.
#' @param ptsb.t a treatment vector of previous patients with length \code{n}.
#' @param ptsb.Y a response vector of previous patients with length \code{n}.
#' @param ptsb.E a censoring indicator vector (1 = event observed, 0 = censored)with length \code{n}.
#' @param ptnow.cov a covariate vector of the incoming patient with length \code{k}.
#' @param gamma a non-negative number. A tuning parameter that reflects the importance of the efficiency component compared to the ethics component.
#' @param event.prob a vector with length 2. The probability of events of upcoming patient for two treatments.
#' @description Calculating the probability of assigning the upcoming patient to treatment A based on the patient's covariates and the previous patients' covariates and responses
#' using CARA Designs Based on Efficiency and Ethics for survival trial.
#' @concept CARAEE Design
#' @returns \item{prob}{Probability of assigning the upcoming patient to treatment A.}
#' @export
#'
#' @examples set.seed(123)
#'n = 40
#'covariates = cbind(rexp(40),rexp(40))
#'treatment = sample(c(0, 1), n, replace = TRUE)
#'survival_time = rexp(n, rate = 1)
#'censoring = runif(n)
#'event = as.numeric(survival_time < censoring)
#'
#'new_patient_cov = c(Z1 = 1, Z2 = 0.5)
#'
#'result = CARAEE_Alloc_Surv(
#'  ptsb.cov = covariates,
#'  ptsb.t = treatment,
#'  ptsb.Y = survival_time,
#'  ptsb.E = event,
#'  ptnow.cov = new_patient_cov,
#'  gamma=1,
#'  event.prob = c(0.5,0.7)
#')
#'print(result$prob)
CARAEE_Alloc_Surv = function(ptsb.cov,
                             ptsb.t,
                             ptsb.Y,
                             ptsb.E,
                             ptnow.cov,
                             gamma,
                             event.prob) {
  # if (nrow(ptsb.cov)<40){
  #   warning("Sample size is too small.Estimations may be biased.")
  # }
  if (length(ptsb.t) != length(ptsb.Y) ||
      nrow(ptsb.cov) != length(ptsb.Y)
      || length(ptsb.t) != nrow(ptsb.cov)) {
    stop(
      "Covariate matrix, treatment vector, and response vector must have the same number of observations."
    )
  }
  colnames(ptsb.cov) = paste0("Z", 1:ncol(ptsb.cov))

  ptsb = data.frame(ptsb.cov,
                    Treat = ptsb.t,
                    Y = ptsb.Y,
                    E = ptsb.E)
  ptsb.A = ptsb[ptsb$Treat == 1, ]
  ptsb.B = ptsb[ptsb$Treat == 0, ]

  formula_str = paste("Surv(Y, E) ~", paste(colnames(ptsb.cov), collapse = " + "))
  fitA = coxph(as.formula(formula_str), data = ptsb.A)
  fitB = coxph(as.formula(formula_str), data = ptsb.B)

  lp1 = sum(ptnow.cov * coef(fitA))
  lp2 = sum(ptnow.cov * coef(fitB))

  e1 = exp(lp2)
  e2 = exp(lp1)

  # efficiency: ε × z^T J^{-1} z
  Jinv1 = vcov(fitA)
  Jinv2 = vcov(fitB)

  d1 = as.numeric(event.prob[1] * t(ptnow.cov) %*% Jinv1 %*% ptnow.cov)
  d2 = as.numeric(event.prob[2] * t(ptnow.cov) %*% Jinv2 %*% ptnow.cov)

  # CARAEE allocation probability
  numer = e1 * d1^gamma
  denom = numer + e2 * d2^gamma
  prob = numer / denom

  return(list(prob = prob))
}


#' Simulation Function of of CARA Designs Based on Efficiency and Ethics for Binary and Continuous Response.
#' @importFrom stats binomial coef glm lm pnorm predict qchisq rnorm runif vcov rexp model.matrix
#' @param n a positive integer. The value specifies the total number of participants involved in each round of the simulation.
#' @param thetaA a vector of length \code{k+1}. The true coefficient parameter value for treatment A.
#' @param thetaB a vector of length \code{k+1}. The true coefficient parameter value for treatment B.
#' @param m0 a positive integer. The number of first 2m0 patients will be allocated equally for estimation. The default value is 40.
#' @param pts.cov a \code{n x k} matrix. The simulated covariate matrix for patients.
#' @param response the type of the response. Options are \code{"Binary"} or \code{"Cont"}.
#' @param gamma a non-negative number. A tuning parameter that reflects the importance of the efficiency component compared to the ethics component.
#' @concept CARAEE Design
#' @description
#' This function simulates a clinical trial using CARA Designs Based on Efficiency and Ethics (CARAEE) with Binary or Continuous Responses.
#' @return A list with the following elements:
#' \item{method}{The name of the procedure.}
#' \item{sampleSize}{Total number of patients.}
#' \item{parameter}{Estimated parameter values.}
#' \item{assignment}{Treatment assignment vector.}
#' \item{proportion}{Proportion of patients allocated to treatment A.}
#' \item{responses}{Simulated response values.}
#' \item{failureRate}{Proportion of treatment failures (if \code{response = "Binary"}).}
#' \item{meanResponse}{Mean response value (if \code{response = "Cont"}).}
#' \item{rejectNull}{Logical. Indicates whether the treatment effect is statistically significant based on a Wald test.}
#'
#' @examples
#'set.seed(123)
#'results = CARAEE_Sim(n = 400,
#'                      pts.cov = cbind(rnorm(400), rnorm(400)),
#'                      thetaA = c(-1, 1, 1),
#'                      thetaB = c(3, 1, 1),
#'                      response = "Binary",
#'                      gamma=1)
#' @export
CARAEE_Sim = function(n,
                      thetaA,
                      thetaB,
                      m0 = 40,
                      pts.cov,
                      response,
                      gamma) {
  if (length(thetaA) != ncol(pts.cov) + 1 ||
      length(thetaB) != ncol(pts.cov) + 1) {
    stop("Length of true parameter vector must be equal to the number of covariates+1.")
  }

  # Check sample size
  if (!is.numeric(n) || length(n) != 1 || n <= 0 || n != floor(n)) {
    stop("`n` must be a positive integer.")
  }

  # Check m0
  if (!is.numeric(m0) || length(m0) != 1 || m0 <= 0 || m0 != floor(m0)) {
    stop("`m0` must be a positive integer.")
  }

  # Check pts.cov
  if (!is.matrix(pts.cov) || nrow(pts.cov) != n) {
    stop("`pts.cov` must be a numeric matrix with `n` rows.")
  }

  k = ncol(pts.cov)

  # Check thetaA and thetaB
  if (length(thetaA) != k + 1) {
    stop("`thetaA` must be of length k + 1, where k is the number of covariates.")
  }
  if (length(thetaB) != k + 1) {
    stop("`thetaB` must be of length k + 1, where k is the number of covariates.")
  }

  # Check v
  if (!is.numeric(gamma) || length(gamma) != 1 || gamma < 0) {
    stop("`gamma` must be a non-negative number.")
  }

  # Check response
  if (!response %in% c("Binary", "Cont")) {
    stop("`response` must be either 'Binary' or 'Cont'.")
  }



  pts = CADBCD_generate_data(
    n = n,
    pts.cov = pts.cov,
    m0 = m0,
    thetaA = thetaA,
    thetaB = thetaB,
    response = response
  )
  for (m in ((2 * m0 + 1):nrow(pts))) {
    ptsb = pts[1:(m - 1), ]
    ptnow.cov=pts[m, 1:ncol(pts.cov)]
    ptnow.cov=data.frame(t(ptnow.cov))
    colnames(ptnow.cov)=colnames(pts)[1:ncol(pts.cov)]

    ptsb.A = data.frame(ptsb[ptsb[, "Treat"] == 1, ])
    ptsb.B = data.frame(ptsb[ptsb[, "Treat"] == 0, ])

    if (response == "Binary") {
      fitA = glm(Y ~ ., data = ptsb.A[, !(names(ptsb.A) %in% "Treat")], family = binomial)
      fitB = glm(Y ~ ., data = ptsb.B[, !(names(ptsb.B) %in% "Treat")], family = binomial)


      pi1 = predict(fitA, newdata = ptnow.cov, type = "response")
      pi2 = predict(fitB, newdata = ptnow.cov, type = "response")

      # Ethics and Efficiency
      e1 = pi1
      e2 = pi2
      d1 = 1 / (pi1 * (1 - pi1))
      d2 = 1 / (pi2 * (1 - pi2))
      }
    else if (response == "Cont") {
      fitA = lm(Y ~ ., data = ptsb.A[, !(names(ptsb.A) %in% "Treat")])
      fitB = lm(Y ~ ., data = ptsb.B[, !(names(ptsb.B) %in% "Treat")])


      ptnow.cov=as.data.frame(ptnow.cov)

      e1 = predict(fitA, newdata = ptnow.cov)
      e2 = predict(fitB, newdata = ptnow.cov)

      # Efficiency component
      X1 = model.matrix(fitA)
      X2 = model.matrix(fitB)
      sigma1_sq = summary(fitA)$sigma^2
      sigma2_sq = summary(fitB)$sigma^2

      invXtX1 = solve(t(X1) %*% X1)
      invXtX2 = solve(t(X2) %*% X2)

      z_vec = c(1, as.numeric(ptnow.cov))  # intercept + covariates
      d1 = sigma1_sq * t(z_vec) %*% invXtX1 %*% z_vec
      d2 = sigma2_sq * t(z_vec) %*% invXtX2 %*% z_vec
    }

    numerator = e1 * d1^gamma
    denominator = numerator + e2 * d2^gamma
    prob = numerator / denominator


    if (anyNA(prob)) {
      pts[m, "Treat"] = sample(c(1, 0), 1, replace = TRUE, c(0.5, 0.5))
    }
    else{
      pts[m, "Treat"] = sample(c(1, 0), 1, replace = TRUE, c(prob, 1 - prob))
    }

    z_m = c(1, pts[m, 1:ncol(pts.cov)])
    if (response == "Binary") {
      exp_m = exp(-(z_m %*% thetaA * pts[m, "Treat"] + z_m %*% thetaB * (1 - pts[m, "Treat"])))
      pts[m, "Y"] = as.numeric(runif(1) < 1 / (1 + exp_m))
    }
    else if (response == "Cont") {
      pts[m, "Y"] = z_m %*% thetaA * pts[m, "Treat"] + z_m %*% thetaB * (1 - pts[m, "Treat"]) +
        rnorm(1)
    }
  }

  return(
    CARA_Output(
      name = "CARA Designs Based on Efficiency and Ethics",
      pts = pts,
      response = response
    )
  )
}

#' Simulation Function of of CARA Designs Based on Efficiency and Ethics for Survival Response.
#' @description
#' This function simulates a clinical trial with time-to-event (survival) outcomes using
#' the CARA Designs Based on Efficiency and Ethics for Survival Response(CARAEE).
#' Patient responses are generated under the Cox proportional hazards model,
#' assuming the proportional hazards assumption holds.
#' @importFrom stats binomial coef glm lm pnorm predict qchisq rnorm runif vcov rexp
#' @importFrom survival Surv coxph survreg
#' @param n a positive integer. The value specifies the total number of participants involved in each round of the simulation.
#' @param thetaA a vector of length \code{k}. The true coefficient parameter value for treatment A.
#' @param thetaB a vector of length \code{k}. The true coefficient parameter value for treatment B.
#' @param m0 a positive integer. The number of first 2m0 patients will be allocated equally for estimation. The default value is 40.
#' @param pts.cov a \code{n x k} matrix. The simulated covariate matrix for patients.
#' @param gamma a non-negative number. A tuning parameter that reflects the importance of the efficiency component compared to the ethics component.
#' @param censor.time a positive value. The upper bound to the simulated uniform censor time.
#' @param arrival.rate a positive value. The rate of simulated exponential arrival time.
#' @return A list with the following elements:
#' \item{method}{The name of procedure.}
#' \item{sampleSize}{Sample size of the trial.}
#' \item{parameter}{Estimated parameters used to do the simulations.}
#' \item{N.events}{Total number of events of the trial.}
#' \item{assignment}{The randomization sequence.}
#' \item{proportion}{Average allocation proportion for treatment A.}
#' \item{responses}{The simulated observed survival responses of patients.}
#' \item{events}{Whether events are observed for patients(1=event,0=censored).}
#' \item{rejectNull}{Logical. Indicates whether the treatment effect is statistically significant based on a Wald test.}
#' @concept CARAEE Design
#' @export
#' @examples
#'set.seed(123)
#'
#'results = CARAEE_Sim_Surv(
#'  thetaA = c(0.1, 0.1),
#'  thetaB = c(-1, 0.1),
#'  n = 400,
#'  pts.cov = cbind(sample(c(1, 0), 400, replace = TRUE), rnorm(400)),
#'  gamma=1,
#'  censor.time = 2,
#'  arrival.rate = 150
#')

CARAEE_Sim_Surv = function(n,
                           thetaA,
                           thetaB,
                           m0 = 40,
                           pts.cov,
                           gamma,
                           censor.time,
                           arrival.rate) {
  if (!is.numeric(n) || length(n) != 1 || n <= 0 || n != floor(n)) {
    stop("`n` must be a single positive integer.")
  }

  if (!is.numeric(m0) || length(m0) != 1 || m0 <= 0 || m0 != floor(m0)) {
    stop("`m0` must be a single positive integer.")
  }

  if (!is.numeric(gamma) || length(gamma) != 1 || gamma < 0 ) {
    stop("`gamma` must be a non-negative number.")
  }

  if (!is.matrix(pts.cov)) {
    stop("`pts.cov` must be a numeric matrix.")
  }

  if (nrow(pts.cov) != n) {
    stop("`pts.cov` must have exactly `n` rows.")
  }

  k = ncol(pts.cov)

  if (length(thetaA) != k) {
    stop("`thetaA` must be a numeric vector of length k, where k is the number of covariates.")
  }

  if (length(thetaB) != k) {
    stop("`thetaB` must be a numeric vector of length k, where k is the number of covariates.")
  }

  if (!is.numeric(censor.time) || length(censor.time) != 1 || censor.time <= 0) {
    stop("`censor.time` must be a single positive number.")
  }

  if (!is.numeric(arrival.rate) || length(arrival.rate) != 1 || arrival.rate <= 0) {
    stop("`arrival.rate` must be a single positive number.")
  }
  pts = CADBCD_generate_data_surv(
    n = n,
    pts.cov = pts.cov,
    m0 = m0,
    thetaA = thetaA,
    thetaB = thetaB,
    censor.time = censor.time,
    arrival.rate = arrival.rate
  )
  colnames(pts.cov) = paste0("Z", 1:ncol(pts.cov))
  for (m in ((2 * m0 + 1):nrow(pts))) {
    current_time = pts[m, "A"]
    ptnow.cov=pts[m,1:ncol(pts.cov)]
    ptsb = pts[1:(m - 1), ]
    ptsb.A = ptsb[ptsb[, "Treat"] == 1, ]
    ptsb.B = ptsb[ptsb[, "Treat"] == 0, ]

    Ai_A = current_time - ptsb.A[, "A"]
    Ai_B = current_time - ptsb.B[, "A"]

    E_star_A = as.numeric(ptsb.A[, "Y"] < Ai_A)
    E_star_B = as.numeric(ptsb.B[, "Y"] < Ai_B)

    ptsb.A_star = cbind(ptsb.A, E_star_A)
    ptsb.B_star = cbind(ptsb.B, E_star_B)

    ptsb.A_star[, "E_star_A"][ptsb.A_star[, "Y"] == ptsb.A_star[, "C"]] = 0
    ptsb.B_star[, "E_star_B"][ptsb.B_star[, "Y"] == ptsb.B_star[, "C"]] = 0

    Y_star_A = pmin(ptsb.A_star[, "Y"], Ai_A)
    Y_star_B = pmin(ptsb.B_star[, "Y"], Ai_B)

    ptsb.A_star = cbind(ptsb.A_star, Y_star_A)
    ptsb.B_star = cbind(ptsb.B_star, Y_star_B)

    formula_strA = paste("Surv(Y_star_A, E_star_A) ~",
                         paste(colnames(pts.cov), collapse = " + "))
    formula_strB = paste("Surv(Y_star_B, E_star_B) ~",
                         paste(colnames(pts.cov), collapse = " + "))
    fitA = coxph(as.formula(formula_strA), data = data.frame(ptsb.A_star))
    fitB = coxph(as.formula(formula_strB), data = data.frame(ptsb.B_star))

    lp1 = sum(ptnow.cov * -coef(fitA))
    lp2 = sum(ptnow.cov * -coef(fitB))

    e1 = exp(lp1)
    e2 = exp(lp2)

    event.probA=get_event_prob_from_cox(fit=fitA,z_new=pts[m,1:ncol(pts.cov)],c_max=censor.time)
    event.probB=get_event_prob_from_cox(fit=fitB,z_new=pts[m,1:ncol(pts.cov)],c_max=censor.time)


    # efficiency: ε × z^T J^{-1} z
    Jinv1 = vcov(fitA)
    Jinv2 = vcov(fitB)

    d1 = as.numeric(event.probA * t(ptnow.cov) %*% Jinv1 %*% ptnow.cov)
    d2 = as.numeric(event.probB * t(ptnow.cov) %*% Jinv2 %*% ptnow.cov)

    # CARAEE allocation probability
    numer = e1 * d1^gamma
    denom = numer + e2 * d2^gamma
    prob = numer / denom

    if (anyNA(prob)) {
      pts[m, "Treat"] = sample(c(1, 0), 1, replace = TRUE, c(0.5, 0.5))
    }
    else{
      pts[m, "Treat"] = sample(c(1, 0), 1, replace = TRUE, c(prob, 1 - prob))
    }

    z_m = pts[m, 1:ncol(pts.cov)]
    linpred = z_m %*% thetaA * pts[m, "Treat"] + z_m %*% thetaB * (1 - pts[m, "Treat"])
    pts[m, "S"] = rexp(1, rate = 1 / exp(linpred))
    pts[m, "E"] = as.numeric(pts[m, "S"] < pts[m, "C"])
    pts[m, "Y"] = pmin(pts[m, "S"], pts[m, "C"])
  }

  return(CARA_Output_Surv(name = "Covariate Adjusted Doubly Biased Coin Design", pts =
                            pts))
}



###############################################################################
####################   Zhao's New Design   ####################################
###############################################################################
#' Allocation Function of Zhao's New Design for Binary and Continuous Response
#' @importFrom stats binomial coef glm lm pnorm predict qchisq rnorm runif vcov rexp model.matrix
#' @param ptsb.X a vector of length \code{n} of the predictive covariates of previous patients. Must be binary.
#' @param ptsb.Z a \code{n x k}  of the prognostic covariates of previous patients. Must be binary.
#' @param ptsb.t a vector of length \code{n} of the treatment allocation of previous patients.
#' @param ptsb.Y a vector of length \code{n} of the responses of previous patients.
#' @param ptnow.X a binary value of the predictive covariate of the present patient.
#' @param ptnow.Z a vector of length \code{k} of the prognostic covariate of the present patient.
#' @param response the type of the response. Options are \code{"Binary"} or \code{"Cont"}.
#' @param omega a vector of length \code{2+k}. The weight of imbalance.
#' @param p a positive value between 0.75 and 0.95. The probability parameter of Efron's biased coin design.
#' @description Calculating the probability of assigning the upcoming patient to treatment A based on the patient's predictive covariates and the previous patients' predictive covariates and responses for Zhao's New procedure.
#' @details
#' This function implements a stratified covariate-adjusted response-adaptive design that balances treatment allocation within and across strata defined by prognostic covariates.
#'
#' The first \eqn{2K} patients are randomized using a restricted randomization procedure, with \eqn{K} patients assigned to each treatment. For patient \eqn{n > 2K}, let \eqn{\mathbf{X}_n} denote predictive covariates, and \eqn{\mathbf{Z}_n} denote stratification covariates, placing the patient into stratum \eqn{(k_1^*, \ldots, k_I^*)}.
#'
#' Based on the covariate profiles and responses of the first \eqn{n-1} patients, we estimate the target allocation probability \eqn{\widehat{\rho}(\mathbf{X}_n)} for the current patient.
#'
#' If assigned to treatment 1, we compute imbalance measures between actual and target allocation at three levels:
#' \itemize{
#'   \item Overall imbalance: \eqn{D_n^{(1)}(\mathbf{X}_n)},
#'   \item Marginal imbalance: \eqn{D_n^{(1)}(i; k_i^*; \mathbf{X}_n)},
#'   \item Within-stratum imbalance: \eqn{D_n^{(1)}(k_1^*, \ldots, k_I^*; \mathbf{X}_n)}.
#' }
#'
#' These are combined into a weighted imbalance function:
#' \deqn{
#' \operatorname{Imb}_n^{(1)}(\mathbf{X}_n) = w_o [D_n^{(1)}(\mathbf{X}_n)]^2 + \sum_{i=1}^I w_{m,i} [D_n^{(1)}(i; k_i^*; \mathbf{X}_n)]^2 + w_s [D_n^{(1)}(k_1^*, \ldots, k_I^*; \mathbf{X}_n)]^2.
#' }
#'
#' A similar imbalance \eqn{\operatorname{Imb}_n^{(2)}(\mathbf{X}_n)} is defined for treatment 2. The patient is then assigned to treatment 1 with probability:
#' \deqn{
#' \phi_n = g\left( \operatorname{Imb}_n^{(1)}(\mathbf{X}_n) - \operatorname{Imb}_n^{(2)}(\mathbf{X}_n) \right),
#' }
#' where \eqn{g(x)} is a biasing function satisfying \eqn{g(-x) = 1 - g(x)} and \eqn{g(x) < 0.5} for \eqn{x \geq 0}. One common choice is Efron's biased coin function:
#' \deqn{
#' g(x) =
#' \begin{cases}
#' q, & \text{if } x > 0 \\
#' 0.5, & \text{if } x = 0 \\
#' p, & \text{if } x < 0
#' \end{cases}
#' }
#' with \eqn{p > 0.5} and \eqn{q = 1 - p}.
#'
#' This design unifies covariate-adjusted response-adaptive randomization and marginal/stratified balance. It reduces to Hu & Hu's design when \eqn{\mathbf{X}_n} is excluded, and to CARA designs when \eqn{\mathbf{Z}_n} is ignored. More detail can be found in Zhao et al.(2022).
#' @references
#' #' Zhao, W., Ma, W., Wang, F., & Hu, F. (2022). Incorporating covariates information in adaptive clinical trials for precision medicine.
#' \emph{Pharmaceutical Statistics}, 21(1), 176–195.
#' @return \item{prob}{Probability of assigning the upcoming patient to treatment A.}
#' @examples
#'set.seed(123)
#'ptsb.X = sample(c(1, -1), 400, replace = TRUE)
#'ptsb.Z = cbind(
#'  sample(c(1, -1), 400, replace = TRUE),
#'  sample(c(1, -1), 400, replace = TRUE)
#')
#'ptsb.Y = sample(c(1, 0), 400, replace = TRUE)
#'ptsb.t = sample(c(1, 0), 400, replace = TRUE)
#'
#'## Incoming patient
#'ptnow.X = 1
#'ptnow.Z = c(1, -1)
#'
#'## Run allocation probability calculation
#'prob = ZhaoNew_Alloc(
#'  ptsb.X = ptsb.X,
#'  ptsb.Z = ptsb.Z,
#'  ptsb.Y = ptsb.Y,
#'  ptsb.t = ptsb.t,
#'  ptnow.X = ptnow.X,
#'  ptnow.Z = ptnow.Z,
#'  response = "Binary",
#'  omega = rep(0.25, 4),
#'  p = 0.8
#')
#'
#'## View allocation probability for treatment A
#'prob
#' @export
#' @concept Zhao's New Design

ZhaoNew_Alloc = function(ptsb.X,
                         ptsb.Z,
                         ptsb.t,
                         ptsb.Y,
                         ptnow.X,
                         ptnow.Z,
                         response,
                         omega,
                         p = 0.8) {
  # Check ptsb.X is binary
  if (length(unique(ptsb.X)) != 2) {
    stop("ptsb.X must be binary (only two unique values).")
  }

  # Check each column of ptsb.Z is binary
  if (!all(apply(ptsb.Z, 2, function(col) length(unique(col)) == 2))) {
    stop("Each column of ptsb.Z must be binary (only two unique values).")
  }

  # Check ptnow.X in the same domain as ptsb.X
  if (!(ptnow.X %in% unique(ptsb.X))) {
    stop("ptnow.X must be one of the two unique values in ptsb.X.")
  }

  # Check ptnow.Z in the same domain as ptsb.Z
  for (j in seq_along(ptnow.Z)) {
    if (!(ptnow.Z[j] %in% unique(ptsb.Z[, j]))) {
      stop(paste0("ptnow.Z[", j, "] must be one of the two unique values in ptsb.Z[,", j, "]."))
    }
  }

  # Check response type
  if (!(response %in% c("Binary", "Cont"))) {
    stop('response must be either "Binary" or "Cont".')
  }

  # Check p
  if (!is.numeric(p) || p < 0.75 || p > 0.95) {
    stop("p must be a numeric value between 0.75 and 0.95.")
  }

  # Check omega
  if (!is.numeric(omega) || length(omega) != (2 + ncol(ptsb.Z))) {
    stop(paste0("omega must be a numeric vector of length ", 2 + ncol(ptsb.Z), "."))
  }

  # Check matching dimensions
  if (length(ptsb.X) != length(ptsb.t) || length(ptsb.X) != length(ptsb.Y)) {
    stop("Length of ptsb.X, ptsb.t, and ptsb.Y must be equal.")
  }

  if (nrow(ptsb.Z) != length(ptsb.X)) {
    stop("Number of rows in ptsb.Z must match the length of ptsb.X.")
  }

  if (length(ptnow.Z) != ncol(ptsb.Z)) {
    stop("Length of ptnow.Z must match the number of columns in ptsb.Z.")
  }

  ptsb = cbind.data.frame(ptsb.X, ptsb.Z, ptsb.t, ptsb.Y)
  colnames(ptsb) = c("X", paste0("Z", 1:ncol(ptsb.Z)), "Treat", "Y")

  ptsb.new1 = rbind(ptsb, c(ptnow.X, ptnow.Z, 1, NA))
  ptsb.new0 = rbind(ptsb, c(ptnow.X, ptnow.Z, 0, NA))

  if (response == "Binary") {
    p1hat = sum(ptsb[, "Y"] == 1 &
                  ptsb[, "X"] == ptnow.X &
                  ptsb[, "Treat"] == 1) / sum(ptsb[, "X"] == ptnow.X &
                                                ptsb[, "Treat"] == 1)
    p0hat = sum(ptsb[, "Y"] == 1 &
                  ptsb[, "X"] == ptnow.X &
                  ptsb[, "Treat"] == 0) / sum(ptsb[, "X"] == ptnow.X  &
                                                ptsb[, "Treat"] == 0)
    rho = sqrt(p1hat) / (sqrt(p1hat) + sqrt(p0hat))
  }

  else if (response == "Cont") {
    ptsb.Xn = ptsb[ptsb[, "X"] == ptnow.X, ]
    fit.Xn = lm(Y ~ Treat - 1, data = ptsb.Xn)
    predict.new1 = data.frame(ptsb.new1[nrow(ptsb.new1), ])
    predict.new0 = data.frame(ptsb.new0[nrow(ptsb.new0), ])
    Y1 = predict(fit.Xn, predict.new1)
    Y0 = predict(fit.Xn, predict.new0)
    rho = pnorm(Y1) / (pnorm(Y1) + pnorm(Y0))
  }

  Dn1 = sum(ptsb.new1[, "Treat"] == 1 &
              ptsb.new1[, "X"] == ptnow.X) - rho * sum(ptsb.new1[, "X"] == ptnow.X)
  Dn2 = sum(ptsb.new0[, "Treat"] == 1 &
              ptsb.new0[, "X"] == ptnow.X) - rho * sum(ptsb.new0[, "X"] == ptnow.X)

  Dn1k = numeric()
  Dn2k = numeric()
  k = ncol(ptsb.Z)
  for (cov in 1:k) {
    Dn1k = append(
      Dn1k,
      sum(ptsb.new1[, "Treat"] == 1 &
            ptsb.new1[, "X"] == ptnow.X &
            ptsb.new1[, paste0("Z", cov)] == ptnow.Z[cov]) -
        rho * sum(ptsb.new1[, "X"] == ptnow.X &
                    ptsb.new1[, paste0("Z", cov)] == ptnow.Z[cov])
    )

    Dn2k = append(
      Dn2k,
      sum(ptsb.new0[, "Treat"] == 1 &
            ptsb.new0[, "X"] == ptnow.X &
            ptsb.new0[, paste0("Z", cov)] == ptnow.Z[cov]) -
        rho * sum(ptsb.new0[, "X"] == ptnow.X &
                    ptsb.new0[, paste0("Z", cov)] == ptnow.Z[cov])
    )
  }

  Dn1k1k2 = sum(ptsb.new1[, "Treat"] == 1 &
                  ptsb.new1[, "X"] == ptnow.X &
                  setequal(ptsb.new1[, paste0("Z", 1:k)], ptnow.Z)) -
    rho * sum(ptsb.new1[, "X"] == ptnow.X &
                setequal(ptsb.new1[, paste0("Z", 1:k)], ptnow.Z))

  Dn2k1k2 = sum(ptsb.new0[, "Treat"] == 1 &
                  ptsb.new0[, "X"] == ptnow.X &
                  setequal(ptsb.new0[, paste0("Z", 1:k)], ptnow.Z)) -
    rho * sum(ptsb.new0[, "X"] == ptnow.X &
                setequal(ptsb.new0[, paste0("Z", 1:k)], ptnow.Z))

  Imbn1 = omega[1] * Dn1^2 + sum(omega[2:(length(omega) - 1)] * Dn1k^2) + omega[length(omega)] * Dn1k1k2^2
  Imbn2 = omega[1] * Dn2^2 + sum(omega[2:(length(omega) - 1)] * Dn2k^2) + omega[length(omega)] * Dn2k1k2^2
  Imb = c(Imbn1, Imbn2)
  names(Imb) = c("Imbalance 1", "Imbalance0")
  if (is.na(Imbn1 > Imbn2) |
      is.na(Imbn1 < Imbn2) |
      is.na(Imbn1 == Imbn2)) {
    return(list(prob = 0.5))
  }
  else if (Imbn1 > Imbn2) {
    return(list(prob = 1 - p))
  }
  else if (Imbn1 == Imbn2) {
    return(list(prob = 0.5))
  }
  else if (Imbn1 < Imbn2) {
    return(list(prob = p))
  }


}


#' Allocation Function of Zhao's Design for Survival Response
#' @importFrom stats binomial coef glm lm pnorm predict qchisq rnorm runif vcov rexp
#' @importFrom survival Surv coxph survreg
#' @param ptsb.X a vector of length \code{n} of the predictive covariates of previous patients. Must be binary.
#' @param ptsb.Z a \code{n x k}  of the prognostic covariates of previous patients. Must be binary.
#' @param ptsb.t a vector of length \code{n} of the treatment allocation of previous patients.
#' @param ptsb.Y a vector of length \code{n} of the responses of previous patients.
#' @param ptsb.E a vector of length \code{n} with value 1 or 0 of the status of event and censoring.
#' @param ptnow.X a binary value of the predictive covariate of the present patient.
#' @param ptnow.Z a vector of length \code{k} of the binary prognostic covariate of the present patient.
#' @param omega a vector of length \code{2+k}. The weight of imbalance.
#' @param p a positive value between 0.75 and 0.95. The probability parameter of Efron's biased coin design.
#' @description Calculating the probability of assigning the upcoming patient to treatment A based on the patient's covariates and the previous patients' covariates and responses for Zhao's New procedure for survival trials.
#' @return \item{prob}{Probability of assigning the upcoming patient to treatment A.}
#' @export
#' @concept Zhao's New Design
#' @examples
#'set.seed(123)
#'
#'# Generate historical data for 400 patients
#'ptsb.X = sample(c(1, -1), 400, replace = TRUE)  # predictive covariate
#'ptsb.Z = cbind(
#'  sample(c(1, -1), 400, replace = TRUE),         # prognostic covariate 1
#'  sample(c(1, -1), 400, replace = TRUE)          # prognostic covariate 2
#')
#'ptsb.Y = rexp(400, rate = 1)                    # survival time (response)
#'ptsb.E = sample(c(1, 0), 400, replace = TRUE)   # event indicator (1 = event, 0 = censored)
#'ptsb.t = sample(c(1, 0), 400, replace = TRUE)   # treatment assignment
#'
#'# Incoming patient covariates
#'ptnow.X = 1
#'ptnow.Z = c(1, -1)
#'
#'# Allocation probability calculation
#'prob = ZhaoNew_Alloc_Surv(
#'  ptsb.X = ptsb.X,
#'  ptsb.Z = ptsb.Z,
#'  ptsb.Y = ptsb.Y,
#'  ptsb.E = ptsb.E,
#'  ptsb.t = ptsb.t,
#'  ptnow.X = ptnow.X,
#'  ptnow.Z = ptnow.Z,
#'  omega = rep(0.25, 4),
#'  p = 0.8
#')
#'
#'# View the allocation probability for treatment A
#'prob

ZhaoNew_Alloc_Surv = function(ptsb.X,
                         ptsb.Z,
                         ptsb.t,
                         ptsb.Y,
                         ptsb.E,
                         ptnow.X,
                         ptnow.Z,
                         omega,
                         p = 0.8) {
  if (length(unique(ptsb.X)) != 2) {
    stop("ptsb.X must be binary (only two unique values).")
  }

  if (!all(apply(ptsb.Z, 2, function(col) length(unique(col)) == 2))) {
    stop("Each column of ptsb.Z must be binary (only two unique values).")
  }

  if (!(ptnow.X %in% unique(ptsb.X))) {
    stop("ptnow.X must match one of the unique values in ptsb.X.")
  }

  for (j in seq_along(ptnow.Z)) {
    if (!(ptnow.Z[j] %in% unique(ptsb.Z[, j]))) {
      stop(paste0("ptnow.Z[", j, "] must match values in ptsb.Z[,", j, "]."))
    }
  }

  if (!is.numeric(p) || p < 0.75 || p > 0.95) {
    stop("p must be a numeric value between 0.75 and 0.95.")
  }

  if (!is.numeric(omega) || length(omega) != (2 + ncol(ptsb.Z))) {
    stop(paste0("omega must be a numeric vector of length ", 2 + ncol(ptsb.Z), "."))
  }

  n = length(ptsb.X)
  if (!all(lengths(list(ptsb.t, ptsb.Y, ptsb.E)) == n)) {
    stop("ptsb.X, ptsb.t, ptsb.Y, and ptsb.E must have the same length.")
  }

  if (nrow(ptsb.Z) != n) {
    stop("Number of rows in ptsb.Z must match the length of ptsb.X.")
  }

  if (length(ptnow.Z) != ncol(ptsb.Z)) {
    stop("Length of ptnow.Z must match number of columns in ptsb.Z.")
  }

  if (any(is.na(c(ptsb.X, ptsb.t, ptsb.Y, ptsb.E, ptsb.Z)))) {
    stop("Input data contains missing values. Please remove or impute them before calling the function.")
  }

  ptsb = cbind.data.frame(ptsb.X, ptsb.Z, ptsb.t, ptsb.Y,ptsb.E)
  colnames(ptsb) = c("X", paste0("Z", 1:ncol(ptsb.Z)), "Treat", "Y","E")
  X_m=ptnow.X

  ptsb.new1 = rbind(ptsb, c(ptnow.X, ptnow.Z, 1, NA,NA))
  ptsb.new0 = rbind(ptsb, c(ptnow.X, ptnow.Z, 0, NA,NA))

  Z_terms = paste0("Z", 1:ncol(ptsb.Z))
  rhs = paste(c("Treat", "X", "X:Treat", Z_terms), collapse = " + ")
  form = as.formula(paste("Surv(Y, E) ~", rhs))

  # 然后拟合模型
  fit = coxph(formula = form, data = data.frame(ptsb))

  b_T  = -coef(fit)["Treat"]
  b_X  = -coef(fit)["X"]
  b_TX = -coef(fit)["Treat:X"]

  # 定义 target function

  logh_A = exp(b_T + (b_X + b_TX) * X_m)
  logh_B = exp(b_X * X_m)

  rho=logh_A/(logh_A+logh_B)

  Dn1 = sum(ptsb.new1[, "Treat"] == 1 &
              ptsb.new1[, "X"] == ptnow.X) - rho * sum(ptsb.new1[, "X"] == ptnow.X)
  Dn2 = sum(ptsb.new0[, "Treat"] == 1 &
              ptsb.new0[, "X"] == ptnow.X) - rho * sum(ptsb.new0[, "X"] == ptnow.X)

  Dn1k = numeric()
  Dn2k = numeric()
  k = ncol(ptsb.Z)
  for (cov in 1:k) {
    Dn1k = append(
      Dn1k,
      sum(ptsb.new1[, "Treat"] == 1 &
            ptsb.new1[, "X"] == ptnow.X &
            ptsb.new1[, paste0("Z", cov)] == ptnow.Z[cov]) -
        rho * sum(ptsb.new1[, "X"] == ptnow.X &
                    ptsb.new1[, paste0("Z", cov)] == ptnow.Z[cov])
    )

    Dn2k = append(
      Dn2k,
      sum(ptsb.new0[, "Treat"] == 1 &
            ptsb.new0[, "X"] == ptnow.X &
            ptsb.new0[, paste0("Z", cov)] == ptnow.Z[cov]) -
        rho * sum(ptsb.new0[, "X"] == ptnow.X &
                    ptsb.new0[, paste0("Z", cov)] == ptnow.Z[cov])
    )
  }

  Dn1k1k2 = sum(ptsb.new1[, "Treat"] == 1 &
                  ptsb.new1[, "X"] == ptnow.X &
                  setequal(ptsb.new1[, paste0("Z", 1:k)], ptnow.Z)) -
    rho * sum(ptsb.new1[, "X"] == ptnow.X &
                setequal(ptsb.new1[, paste0("Z", 1:k)], ptnow.Z))

  Dn2k1k2 = sum(ptsb.new0[, "Treat"] == 1 &
                  ptsb.new0[, "X"] == ptnow.X &
                  setequal(ptsb.new0[, paste0("Z", 1:k)], ptnow.Z)) -
    rho * sum(ptsb.new0[, "X"] == ptnow.X &
                setequal(ptsb.new0[, paste0("Z", 1:k)], ptnow.Z))

  Imbn1 = omega[1] * Dn1^2 + sum(omega[2:(length(omega) - 1)] * Dn1k^2) + omega[length(omega)] * Dn1k1k2^2
  Imbn2 = omega[1] * Dn2^2 + sum(omega[2:(length(omega) - 1)] * Dn2k^2) + omega[length(omega)] * Dn2k1k2^2
  Imb = c(Imbn1, Imbn2)


  if (is.na(Imbn1 > Imbn2) |
      is.na(Imbn1 < Imbn2) |
      is.na(Imbn1 == Imbn2)) {
    return(list(prob = 0.5))
  }
  else if (Imbn1 > Imbn2) {
    return(list(prob = 1 - p))
  }
  else if (Imbn1 == Imbn2) {
    return(list(prob = 0.5))
  }
  else if (Imbn1 < Imbn2) {
    return(list(prob = p))
  }
}



#' Simulation Function of Zhao's New Design for Binary and Continuous Response
#' @importFrom stats binomial coef glm lm pnorm predict qchisq rnorm runif vcov rexp
#' @param n a positive integer. The sample size of the simulated data.
#' @param mu a vector of length 2. The true parameters of treatment.
#' @param beta a vector of length 2. The true parameters of predictive covariate and interaction with treatment.
#' @param gamma a vector of length k. The true parameters of prognostic covariates.
#' @param m0 a positive integer. The number of first 2m0 patients will be allocated equally to both treatments.
#' @param pts.X a vector of length n. The vector of patients' binary predictive covariates.
#' @param pts.Z a matrix of \code{n x k}. The matrix of patients' binary prognostic covariates.
#' @param response the type of the response. Options are \code{"Binary"} or \code{"Cont"}.
#' @param omega a vector of length \code{2+k}. The weight of imbalance.
#' @param p a positive value between 0.75 and 0.95. The probability parameter of Efron's biased coin design.
#' @description This function simulates a trial using Zhao's new design for binary and continuous responses.
#' @return
#' \item{method}{The name of procedure.}
#' \item{sampleSize}{The sample size of the trial.}
#' \item{assignment}{The randomization sequence.}
#' \item{X1proportion}{Average allocation proportion for treatment A when predictive covariate equals the smaller value.}
#' \item{X2proportion}{Average allocation proportion for treatment A when predictive covariate equals the larger value.}
#' \item{proportion}{Average allocation proportion for treatment A.}
#' \item{failureRate}{Proportion of treatment failures (if \code{response = "Binary"}).}
#' \item{meanResponse}{Mean response value (if \code{response = "Cont"}).}
#' \item{rejectNull}{Logical. Indicates whether the treatment effect is statistically significant based on a Wald test.}
#' @export
#' @concept Zhao's New Design
#' @examples
#' set.seed(123)
#'
#' # Simulation settings
#' n = 400                    # total number of patients
#' m0 = 40                    # initial burn-in sample size
#' mu = c(0.5, 0.8)           # potential means (for continuous or logistic link)
#' beta = c(1, 1)             # treatment effect and predictive covariate effect
#' gamma = c(0.1, 0.5)        # prognostic covariate effects
#' omega = rep(0.25, 4)       # imbalance weights
#' p = 0.8                    # biased coin probability
#'
#' # Generate patient covariates
#' pts.X = sample(c(1, -1), n, replace = TRUE)  # predictive covariate
#' pts.Z = cbind(
#'   sample(c(1, -1), n, replace = TRUE),        # prognostic Z1
#'   sample(c(1, -1), n, replace = TRUE)         # prognostic Z2
#' )
#'
#' # Run the simulation (binary response setting)
#' result = ZhaoNew_Sim(
#'   n = n,
#'   mu = mu,
#'   beta = beta,
#'   gamma = gamma,
#'   m0 = m0,
#'   pts.X = pts.X,
#'   pts.Z = pts.Z,
#'   response = "Binary",
#'   omega = omega,
#'   p = p
#' )
#'

ZhaoNew_Sim = function(n,
                       mu,
                       beta,
                       gamma,
                       m0 = 40,
                       pts.X,
                       pts.Z,
                       response,
                       omega,
                       p = 0.8) {
  if (!is.numeric(n) || length(n) != 1 || n <= 0 || n != as.integer(n)) {
    stop("n must be a positive integer.")
  }

  if (!is.numeric(m0) || length(m0) != 1 || m0 <= 0 || m0 != as.integer(m0)) {
    stop("m0 must be a positive integer.")
  }

  if (2 * m0 >= n) {
    stop("2 * m0 must be less than n.")
  }

  if (!is.numeric(mu) || length(mu) != 2) {
    stop("mu must be a numeric vector of length 2.")
  }

  if (!is.numeric(beta) || length(beta) != 2) {
    stop("beta must be a numeric vector of length 2.")
  }

  if (!is.matrix(pts.Z)) {
    stop("pts.Z must be a matrix.")
  }

  k = ncol(pts.Z)

  if (!is.numeric(gamma) || length(gamma) != k) {
    stop(paste0("gamma must be a numeric vector of length ", k, "."))
  }

  if (!is.numeric(omega) || length(omega) != (2 + k)) {
    stop(paste0("omega must be a numeric vector of length ", 2 + k, "."))
  }

  if (!is.numeric(p) || p < 0.75 || p > 0.95) {
    stop("p must be a numeric value between 0.75 and 0.95.")
  }

  if (!is.vector(pts.X) || length(pts.X) != n) {
    stop("pts.X must be a vector of length n.")
  }

  if (nrow(pts.Z) != n) {
    stop("pts.Z must have n rows.")
  }

  if (!all(apply(pts.Z, 2, function(col) length(unique(col)) == 2))) {
    stop("Each column in pts.Z must be binary (i.e., only two unique values).")
  }

  if (length(unique(pts.X)) != 2) {
    stop("pts.X must be binary (i.e., contain exactly two unique values).")
  }

  if (!(response %in% c("Binary", "Cont"))) {
    stop('response must be either "Binary" or "Cont".')
  }

  pts = ZhaoNew_generate_data(
    n = n,
    pts.X = pts.X,
    pts.Z = pts.Z,
    mu = mu,
    beta = beta,
    gamma = gamma,
    m0 = m0,
    response = response
  )
  for (m in (2 * m0 + 1):nrow(pts)) {
    ptsb = pts[1:(m - 1), ]
    ptsb.new1 = pts[1:m, ]
    ptsb.new1[m, "Treat"] = 1
    ptsb.new0 = pts[1:m, ]
    ptsb.new0[m, "Treat"] = 0
    ptnow.X = pts[m, "X"]
    ptnow.Z = pts.Z[m, ]

    if (response == "Binary") {
      p1hat = sum(ptsb[, "Y"] == 1 &
                    ptsb[, "X"] == ptnow.X &
                    ptsb[, "Treat"] == 1) / sum(ptsb[, "X"] == ptnow.X &
                                                  ptsb[, "Treat"] == 1)
      p0hat = sum(ptsb[, "Y"] == 1 &
                    ptsb[, "X"] == ptnow.X &
                    ptsb[, "Treat"] == 0) / sum(ptsb[, "X"] == ptnow.X  &
                                                  ptsb[, "Treat"] == 0)
      rho = sqrt(p1hat) / (sqrt(p1hat) + sqrt(p0hat))
    }

    else if (response == "Cont") {
      ptsb.Xn = ptsb[ptsb[, "X"] == ptnow.X, ]
      fit.Xn = lm(Y ~ Treat - 1, data = ptsb.Xn)
      predict.new1 = data.frame(ptsb.new1[nrow(ptsb.new1), ])
      predict.new0 = data.frame(ptsb.new0[nrow(ptsb.new0), ])
      Y1 = predict(fit.Xn, predict.new1)
      Y0 = predict(fit.Xn, predict.new0)
      rho = pnorm(Y1) / (pnorm(Y1) + pnorm(Y0))
    }

    Dn1 = sum(ptsb.new1[, "Treat"] == 1 &
                ptsb.new1[, "X"] == ptnow.X) - rho * sum(ptsb.new1[, "X"] == ptnow.X)
    Dn2 = sum(ptsb.new0[, "Treat"] == 1 &
                ptsb.new0[, "X"] == ptnow.X) - rho * sum(ptsb.new0[, "X"] == ptnow.X)

    Dn1k = numeric()
    Dn2k = numeric()
    k = ncol(pts.Z)
    for (cov in 1:k) {
      Dn1k = append(
        Dn1k,
        sum(
          ptsb.new1[, "Treat"] == 1 &
            ptsb.new1[, "X"] == ptnow.X &
            ptsb.new1[, paste0("Z", cov)] == ptnow.Z[cov]
        ) -
          rho * sum(ptsb.new1[, "X"] == ptnow.X &
                      ptsb.new1[, paste0("Z", cov)] == ptnow.Z[cov])
      )

      Dn2k = append(
        Dn2k,
        sum(
          ptsb.new0[, "Treat"] == 1 &
            ptsb.new0[, "X"] == ptnow.X &
            ptsb.new0[, paste0("Z", cov)] == ptnow.Z[cov]
        ) -
          rho * sum(ptsb.new0[, "X"] == ptnow.X &
                      ptsb.new0[, paste0("Z", cov)] == ptnow.Z[cov])
      )
    }

    Dn1k1k2 = sum(ptsb.new1[, "Treat"] == 1 &
                    ptsb.new1[, "X"] == ptnow.X &
                    setequal(ptsb.new1[, paste0("Z", 1:k)], ptnow.Z)) -
      rho * sum(ptsb.new1[, "X"] == ptnow.X &
                  setequal(ptsb.new1[, paste0("Z", 1:k)], ptnow.Z))

    Dn2k1k2 = sum(ptsb.new0[, "Treat"] == 1 &
                    ptsb.new0[, "X"] == ptnow.X &
                    setequal(ptsb.new0[, paste0("Z", 1:k)], ptnow.Z)) -
      rho * sum(ptsb.new0[, "X"] == ptnow.X &
                  setequal(ptsb.new0[, paste0("Z", 1:k)], ptnow.Z))

    Imbn1 = omega[1] * Dn1^2 + sum(omega[2:(length(omega) - 1)] * Dn1k^2) + omega[length(omega)] * Dn1k1k2^2
    Imbn2 = omega[1] * Dn2^2 + sum(omega[2:(length(omega) - 1)] * Dn2k^2) + omega[length(omega)] * Dn2k1k2^2
    Imb = c(Imbn1, Imbn2)
    names(Imb) = c("Imbalance 1", "Imbalance0")
    if (is.na(Imbn1 > Imbn2) |
        is.na(Imbn1 < Imbn2) |
        is.na(Imbn1 == Imbn2)) {
      pts[m, "Treat"] = sample(c(1, 0), 1, prob = c(0.5, 0.5))
    }
    else if (Imbn1 > Imbn2) {
      pts[m, "Treat"] = sample(c(1, 0), 1, prob = c(1 - p, p))
    }
    else if (Imbn1 == Imbn2) {
      pts[m, "Treat"] = sample(c(1, 0), 1, prob = c(0.5, 0.5))
    }
    else if (Imbn1 < Imbn2) {
      pts[m, "Treat"] = sample(c(1, 0), 1, prob = c(p, 1 - p))
    }

    if (response == "Binary") {
      term.m = c(pts[m, "Treat"], 1 - pts[m, "Treat"]) %*% mu +
        c(pts[m, "X"], pts[m, "X"] * pts[m, "Treat"]) %*% beta +
        pts.Z[m, ] %*% gamma
      prob.Y0.m = 1 / (1 + exp(-term.m))
      pts[m, "Y"] = as.numeric(runif(1) < prob.Y0.m)
    }
    else if (response == "Cont") {
      pts[m, "Y"] = c(pts[m, "Treat"], 1 - pts[m, "Treat"]) %*% mu +
        c(pts[m, "X"], pts[m, "X"] * pts[m, "Treat"]) %*% beta +
        pts.Z[m, ] %*% gamma + rnorm(1)
    }
  }
  return(ZhaoNew_Output(
    name = "Zhao's New Design",
    pts = pts,
    response = response
  ))
}


#' Simulation Function for Zhao's New Design for Survival Trial
#' @importFrom stats binomial coef glm lm pnorm predict qchisq rnorm runif vcov rexp
#' @importFrom survival Surv coxph survreg
#' @param n a positive integer. The sample size of the simulated data.
#' @param mu a number. The true parameters of treatment.
#' @param beta a vector of length 2. The true parameters of predictive covariate and interaction with treatment.
#' @param gamma a vector of length k. The true parameters of prognostic covariates.
#' @param m0 a positive integer. The number of first 2m0 patients will be allocated equally to both treatments.
#' @param pts.X a vector of length n. The vector of patients' binary predictive covariates.Must be binary.
#' @param pts.Z a matrix of \code{n x k}. The matrix of patients' binary prognostic covariates.Must be binary.
#' @param censor.time a positive number. The upper bound of the uniform censor time in year.
#' @param arrival.rate a positive integer. The arrival rate of patients each year.
#' @param omega a vector of length \code{2+k}. The weight of imbalance.
#' @param p a positive value between 0.75 and 0.95. The probability parameter of Efron's biased coin design.
#' @concept Zhao's New Design
#' @description
#' This function simulates a clinical trial using Zhao's New design for survival responses.
#' @return A list with the following elements:
#' \item{method}{The name of procedure.}
#' \item{sampleSize}{The sample size of the trial.}
#' \item{assignment}{The randomization sequence.}
#' \item{X1proportion}{Average allocation proportion for treatment A when predictive covariate equals the smaller value.}
#' \item{X2proportion}{Average allocation proportion for treatment A when predictive covariate equals the larger value.}
#' \item{proportion}{Average allocation proportion for treatment A.}
#' \item{N.events}{Total number of events occured of the trial.}
#' \item{responses}{Observed survival responses of patients.}
#' \item{events}{Survival status vector of patients(1=event,0=censored)}
#' \item{rejectNull}{Logical. Indicates whether the treatment effect is statistically significant based on a Wald test.}

#' @export
#' @examples

#' set.seed(123)
#'
#' # Simulation settings
#' n = 400                    # total number of patients
#' m0 = 40                    # initial burn-in sample size
#' mu = 0.5           # potential means (for continuous or logistic link)
#' beta = c(1, 1)             # treatment effect and predictive covariate effect
#' gamma = c(0.1, 0.5)        # prognostic covariate effects
#' omega = rep(0.25, 4)       # imbalance weights
#' p = 0.8                    # biased coin probability
#' censor.time = 2            # maximum censoring time
#' arrival.rate = 150         # arrival rate of patients
#' # Generate patient covariates
#' pts.X = sample(c(1, -1), n, replace = TRUE)  # predictive covariate
#' pts.Z = cbind(
#'   sample(c(1, -1), n, replace = TRUE),        # prognostic Z1
#'   sample(c(1, -1), n, replace = TRUE)         # prognostic Z2
#' )
#'
#' # Run the simulation (binary response setting)
#' result = ZhaoNew_Sim_Surv(
#'   n = n,
#'   mu = mu,
#'   beta = beta,
#'   gamma = gamma,
#'   m0 = m0,
#'   pts.X = pts.X,
#'   pts.Z = pts.Z,
#'   omega = omega,
#'   p = p,
#'   censor.time = censor.time,
#'   arrival.rate = arrival.rate
#' )
#'
ZhaoNew_Sim_Surv = function(n,
                       mu,
                       beta,
                       gamma,
                       m0 = 40,
                       pts.X,
                       pts.Z,
                       censor.time,
                       arrival.rate,
                       omega,
                       p = 0.8) {
  # Sample size
  if (!is.numeric(n) || length(n) != 1 || n <= 0 || n != as.integer(n)) {
    stop("n must be a positive integer.")
  }

  # Burn-in
  if (!is.numeric(m0) || length(m0) != 1 || m0 <= 0 || m0 != as.integer(m0)) {
    stop("m0 must be a positive integer.")
  }

  if (2 * m0 >= n) {
    stop("2 * m0 must be less than n.")
  }

  # mu: scalar
  if (!is.numeric(mu) || length(mu) != 1) {
    stop("mu must be a single numeric value (scalar).")
  }

  # beta
  if (!is.numeric(beta) || length(beta) != 2) {
    stop("beta must be a numeric vector of length 2.")
  }

  # pts.X
  if (!is.vector(pts.X) || length(pts.X) != n) {
    stop("pts.X must be a vector of length n.")
  }

  if (length(unique(pts.X)) != 2) {
    stop("pts.X must be binary (i.e., contain exactly two unique values).")
  }

  # pts.Z
  if (!is.matrix(pts.Z)) {
    stop("pts.Z must be a matrix.")
  }

  if (nrow(pts.Z) != n) {
    stop("pts.Z must have n rows.")
  }

  k = ncol(pts.Z)

  if (!all(apply(pts.Z, 2, function(col) length(unique(col)) == 2))) {
    stop("Each column in pts.Z must be binary (i.e., only two unique values).")
  }

  # gamma
  if (!is.numeric(gamma) || length(gamma) != k) {
    stop(paste0("gamma must be a numeric vector of length ", k, "."))
  }

  # omega
  if (!is.numeric(omega) || length(omega) != (2 + k)) {
    stop(paste0("omega must be a numeric vector of length ", 2 + k, "."))
  }

  # p
  if (!is.numeric(p) || length(p) != 1 || p < 0.75 || p > 0.95) {
    stop("p must be a numeric value between 0.75 and 0.95.")
  }

  # censor.time
  if (!is.numeric(censor.time) || length(censor.time) != 1 || censor.time <= 0) {
    stop("censor.time must be a positive number.")
  }

  # arrival.rate
  if (!is.numeric(arrival.rate) || length(arrival.rate) != 1 || arrival.rate <= 0) {
    stop("arrival.rate must be a positive number.")
  }

  pts = ZhaoNew_generate_data_Surv(
    n = n,
    pts.X = pts.X,
    pts.Z = pts.Z,
    mu = mu,
    beta = beta,
    gamma = gamma,
    m0 = m0,
    censor.time=censor.time,
    arrival.rate=arrival.rate
  )
  for (m in (2 * m0 + 1):nrow(pts)) {
    current_time=pts[m,"A"]
    pts_before = pts[1:(m - 1), ]

    Ai=current_time-pts_before[,"A"]

    E_star=as.numeric(pts_before[,"Y"]<Ai)

    pts_before_star=cbind(pts_before,E_star)

    pts_before_star[,"E_star"][pts_before_star[,"Y"]==pts_before_star[,"C"]]=0

    Y_star=pmin(pts_before_star[,"Y"],Ai)

    pts_before_star=cbind(pts_before_star,Y_star)

    X_m=pts[m,"X"]

    pts_before_Xm=pts_before_star[pts_before_star[,"X"]==X_m,]

    pts_new1=pts[m,]
    pts_new0=pts[m,]

    pts_new1["Treat"]=1
    pts_new0["Treat"]=0

    Z1=pts_new1["Z1"]
    Z2=pts_new1["Z2"]

    data_new1=rbind(pts_before,pts_new1)
    data_new0=rbind(pts_before,pts_new0)

    Z_terms = grep("^Z", colnames(pts_before), value = TRUE)
    rhs = paste(c("Treat", "X", "X:Treat", Z_terms), collapse = " + ")
    form = as.formula(paste("Surv(Y, E) ~", rhs))

    fit = coxph(formula = form, data = data.frame(pts_before))

    b_T  = -coef(fit)["Treat"]
    b_X  = -coef(fit)["X"]
    b_TX = -coef(fit)["Treat:X"]

    logh_A = exp(b_T + (b_X + b_TX) * X_m)
    logh_B = exp(b_X * X_m)

    rho=logh_A/(logh_A+logh_B)

    Dn1 = sum(data_new1[, "Treat"] == 1 &
                data_new1[, "X"] == X_m) - rho * sum(data_new1[, "X"] == X_m)
    Dn2 = sum(data_new0[, "Treat"] ==1 &
                data_new0[, "X"] == X_m) - rho * sum(data_new0[, "X"] == X_m)

    Dn1k = numeric()
    Dn2k = numeric()

    Dn1k=append(
      Dn1k,
      sum(
        data_new1[, "Treat"] == 1 &
          data_new1[, "X"] == X_m &
          data_new1[, "Z1"] == Z1
      ) -
        rho * sum(data_new1[, "X"] == X_m &
                    data_new1[, "Z1"] == Z1)
    )

    Dn1k=append(
      Dn1k,
      sum(
        data_new1[, "Treat"] == 1 &
          data_new1[, "X"] == X_m &
          data_new1[, "Z2"] == Z2
      ) -
        rho * sum(data_new1[, "X"] == X_m &
                    data_new1[, "Z2"] == Z2)
    )

    Dn2k=append(
      Dn2k,
      sum(
        data_new0[, "Treat"] == 1 &
          data_new0[, "X"] == X_m &
          data_new0[, "Z1"] == Z1
      ) -
        rho * sum(data_new0[, "X"] == X_m &
                    data_new0[, "Z1"] == Z1)
    )

    Dn2k=append(
      Dn2k,
      sum(
        data_new0[, "Treat"] == 1 &
          data_new0[, "X"] == X_m &
          data_new0[, "Z2"] == Z2
      ) -
        rho * sum(data_new0[, "X"] == X_m &
                    data_new0[, "Z2"] == Z2)
    )


    Dn1k1k2 = sum(data_new1[, "Treat"] == 1 &
                    data_new1[, "X"] == X_m &
                    data_new1[, "Z1"] == Z1 &
                    data_new1[, "Z2"] == Z2) -
      rho * sum(data_new1[, "X"] == X_m &
                  data_new1[, "Z1"] == Z1&
                  data_new1[, "Z2"] == Z2)

    Dn2k1k2 = sum(data_new0[, "Treat"] == 1 &
                    data_new0[, "X"] == X_m &
                    data_new0[, "Z1"] == Z1 &
                    data_new0[, "Z2"] == Z2) -
      rho * sum(data_new0[, "X"] == X_m &
                  data_new0[, "Z1"] == Z1&
                  data_new0[, "Z2"] == Z2)

    Imbn1 = omega[1] * Dn1 ^ 2 + sum(omega[2:3] * Dn1k ^
                                       2) + omega[4] * Dn1k1k2 ^ 2
    Imbn2 = omega[1] * Dn2 ^ 2 + sum(omega[2:3] * Dn2k ^
                                       2) + omega[4] * Dn2k1k2 ^ 2


    if (is.na(Imbn1 > Imbn2) |
        is.na(Imbn1 < Imbn2) |
        is.na(Imbn1 == Imbn2)) {
      pts[m, "Treat"] = sample(c(1, 0), 1, prob = c(0.5, 0.5))
    }
    else if (Imbn1 > Imbn2) {
      pts[m, "Treat"] = sample(c(1, 0), 1, prob = c(1 - p, p))
    }
    else if (Imbn1 == Imbn2) {
      pts[m, "Treat"] = sample(c(1, 0), 1, prob = c(0.5, 0.5))
    }
    else if (Imbn1 < Imbn2) {
      pts[m, "Treat"] = sample(c(1, 0), 1, prob = c(p, 1 - p))
    }
    # pts[m,"S"]=-log(runif(1))/lambda_m
    surv.rate.m=exp(c(pts[m, "Treat"]) * mu +
      c(pts[m, "X"], pts[m, "X"] * pts[m, "Treat"]) %*% beta +
      pts.Z[m, ] %*% gamma)

    pts[m,"S"]=rexp(1,rate=1/surv.rate.m)
    pts[m,"E"]=as.numeric(pts[m,"S"]<pts[m,"C"])
    pts[m,"Y"]=min(pts[m,"S"],pts[m,"C"])
  }

  return(ZhaoNew_Output_Surv(name="Zhao's New Design",pts))

}


###############################################################################
####################   Weighted Balance Ratio Design   ########################
###############################################################################
#' Allocation Function of Weighted Balance Ratio Design for Binary and Continuous Response
#' @importFrom stats binomial coef glm lm pnorm predict qchisq rnorm runif vcov rexp
#' @param ptsb.X a vector of length \code{n} of the predictive covariates of previous patients. Must be binary.
#' @param ptsb.Z a \code{n x k}  of the prognostic covariates of previous patients. Must be binary.
#' @param ptsb.t a vector of length \code{n} of the treatment allocation of previous patients.
#' @param ptsb.Y a vector of length \code{n} of the responses of previous patients.
#' @param ptnow.X a binary value of the predictive covariate of the present patient.
#' @param ptnow.Z a vector of length \code{k} of the binary prognostic covariate of the present patient.
#' @param weight a vector of length \code{2+k}. The weight of balance ratio in overall,margin and stratum levels.
#' @param response the type of the response. Options are \code{"Binary"} or \code{"Cont"}.
#' @param v a positive value that controls the randomness of allocation probability function.
#' @description Calculating the probability of assigning the upcoming patient to treatment A based on the patient's covariates and the previous patients' covariates and responses for Zhao's New procedure.
#' @details
#' This function implements a covariate-adjusted response-adaptive design using a weighted balancing rule (WBR) combined with a CARA-type allocation function.
#'
#' The first two steps follow Zhao et al. (2022): the first \eqn{2K} patients are randomized using restricted randomization with \eqn{K} patients in each treatment group. For patient \eqn{n > 2K}, suppose the prognostic covariates \eqn{\boldsymbol{Z}_n} fall into stratum \eqn{(k_1^*, \ldots, k_J^*)}.
#'
#' Define the weighted imbalance ratio \eqn{\boldsymbol{r}_{n-1}(\boldsymbol{X}_n)} for patient \eqn{n} as:
#' \deqn{
#' \boldsymbol{r}_{n-1}(\boldsymbol{X}_n) =
#' w_o \frac{\boldsymbol{N}_{n-1}(\boldsymbol{X}_n)}{N_{n-1}(\boldsymbol{X}_n)} +
#' \sum_{j=1}^{J} w_{m,j} \frac{\boldsymbol{N}_{(j; k_j^*), n-1}(\boldsymbol{X}_n)}{N_{(j; k_j^*), n-1}(\boldsymbol{X}_n)} +
#' w_s \frac{\boldsymbol{N}_{(k_1^*, \ldots, k_J^*), n-1}(\boldsymbol{X}_n)}{N_{(k_1^*, \ldots, k_J^*), n-1}(\boldsymbol{X}_n)},
#' }
#' where \eqn{w_o}, \eqn{w_{m,1}}, \ldots, \eqn{w_{m,J}}, \eqn{w_s} are non-negative weights that sum to 1.
#'
#' Based on the patient’s covariates, the target allocation probability \eqn{\hat{\rho}_{n-1}(\boldsymbol{X}_n)} is estimated from previous data.
#'
#' The final probability of assigning patient \eqn{n} to treatment \eqn{k} is computed using a CARA-type allocation function from Hu and Zhang (2009):
#' \deqn{
#' \phi_{n,k}(\boldsymbol{X}_n) = g_k(\boldsymbol{r}_{n-1}(\boldsymbol{X}_n), \hat{\rho}_{n-1}(\boldsymbol{X}_n)) =
#' \frac{ \hat{\rho}_{n-1} \left( \frac{\hat{\rho}_{n-1}}{r_{n-1}} \right)^v }
#' { \hat{\rho}_{n-1} \left( \frac{\hat{\rho}_{n-1}}{r_{n-1}} \right)^v + (1 - \hat{\rho}_{n-1}) \left( \frac{1 - \hat{\rho}_{n-1}}{1 - r_{n-1}} \right)^v },
#' }
#' where \eqn{v \geq 0} is a tuning parameter controlling the degree of randomness. A value of \eqn{v = 2} is commonly recommended.
#'
#' This approach combines stratified covariate balancing with covariate-adjusted optimal targeting, and supports both discrete and continuous covariate inputs. More details can be found in Yu(2025).
#' @references
#' Yu, J. (2025). A New Family of Covariate-Adjusted Response-Adaptive Randomization Procedures for Precision Medicine (Doctoral dissertation, The George Washington University).
#' @return \item{prob}{Probability of assigning the upcoming patient to treatment A.}
#' @export
#' @concept Weighted Balance Ratio Design
#' @examples
#' set.seed(123)
#'
#' # Generate historical data for 400 patients
#' ptsb.X = sample(c(1, -1), 400, replace = TRUE)  # predictive covariate
#' ptsb.Z = cbind(
#'   sample(c(1, -1), 400, replace = TRUE),         # prognostic covariate 1
#'   sample(c(1, -1), 400, replace = TRUE)          # prognostic covariate 2
#' )
#' ptsb.Y = sample(c(1, 0), 400, replace = TRUE)   # binary responses
#' ptsb.t = sample(c(1, 0), 400, replace = TRUE)   # treatment assignments
#'
#' # Incoming patient covariates
#' ptnow.X = 1
#' ptnow.Z = c(1, -1)
#'
#' # Calculate allocation probability using WBR method
#' result = WBR_Alloc(
#'   ptsb.X = ptsb.X,
#'   ptsb.Z = ptsb.Z,
#'   ptsb.Y = ptsb.Y,
#'   ptsb.t = ptsb.t,
#'   ptnow.X = ptnow.X,
#'   ptnow.Z = ptnow.Z,
#'   weight = rep(0.25, 4),
#'   response = "Binary"
#' )
#'
#' # View probability of assigning to treatment A
#' result$prob
WBR_Alloc = function(ptsb.X,
                        ptsb.Z,
                        ptsb.Y,
                        ptsb.t,
                        ptnow.X,
                        ptnow.Z,
                        weight,
                        v=2,
                        response) {
  # ptsb.X
  if (!is.vector(ptsb.X) || !is.numeric(ptsb.X)) {
    stop("ptsb.X must be a numeric vector.")
  }

  if (length(unique(ptsb.X)) != 2) {
    stop("ptsb.X must be binary (i.e., contain exactly two unique values).")
  }

  n = length(ptsb.X)

  # ptsb.Z
  if (!is.matrix(ptsb.Z)) {
    stop("ptsb.Z must be a matrix.")
  }

  if (nrow(ptsb.Z) != n) {
    stop("ptsb.Z must have the same number of rows as ptsb.X.")
  }

  if (!all(apply(ptsb.Z, 2, function(col) length(unique(col)) == 2))) {
    stop("Each column in ptsb.Z must be binary (i.e., only two unique values).")
  }

  k = ncol(ptsb.Z)

  # ptsb.Y, ptsb.t
  if (!is.numeric(ptsb.Y) || length(ptsb.Y) != n) {
    stop("ptsb.Y must be a numeric vector of length equal to ptsb.X.")
  }

  if (!is.numeric(ptsb.t) || length(ptsb.t) != n) {
    stop("ptsb.t must be a numeric vector of length equal to ptsb.X.")
  }

  # ptnow.X
  if (!(ptnow.X %in% unique(ptsb.X))) {
    stop("ptnow.X must be one of the values in ptsb.X.")
  }

  # ptnow.Z
  if (!is.vector(ptnow.Z) || length(ptnow.Z) != k) {
    stop(paste0("ptnow.Z must be a vector of length ", k, "."))
  }

  for (j in seq_len(k)) {
    if (!(ptnow.Z[j] %in% unique(ptsb.Z[, j]))) {
      stop(paste0("ptnow.Z[", j, "] must be one of the values in ptsb.Z[,", j, "]."))
    }
  }

  # weight
  if (!is.numeric(weight) || length(weight) != (2 + k)) {
    stop(paste0("weight must be a numeric vector of length ", 2 + k, "."))
  }

  # v
  if (!is.numeric(v) || v <= 0) {
    stop("v must be a positive numeric value.")
  }

  # response
  if (!(response %in% c("Binary", "Cont"))) {
    stop('response must be either "Binary" or "Cont".')
  }

  ptsb = cbind.data.frame(ptsb.X, ptsb.Z, ptsb.t, ptsb.Y)
  colnames(ptsb) = c("X", paste0("Z", 1:ncol(ptsb.Z)), "Treat", "Y")

  Z_terms = paste0("Z", 1:ncol(ptsb.Z))
  rhs = paste(c("Treat", "X", "X:Treat", Z_terms), collapse = " + ")
  form = as.formula(paste("Y ~", rhs))

  if (response == "Binary") {
    p1hat = sum(ptsb[, "Y"] == 1 &
                  ptsb[, "X"] == ptnow.X &
                  ptsb[, "Treat"] == 1) / sum(ptsb[, "X"] == ptnow.X &
                                                ptsb[, "Treat"] == 1)
    p0hat = sum(ptsb[, "Y"] == 1 &
                  ptsb[, "X"] == ptnow.X &
                  ptsb[, "Treat"] == 0) / sum(ptsb[, "X"] == ptnow.X  &
                                                ptsb[, "Treat"] == 0)
    pi_m = sqrt(p1hat) / (sqrt(p1hat) + sqrt(p0hat))
  }

  else if (response == "Cont") {
    fit = lm(formula = form, data = data.frame(ptsb))
    b_I = coef(fit)["(Intercept)"]
    b_T  = coef(fit)["Treat"]
    b_X  = coef(fit)["X"]
    b_TX = coef(fit)["Treat:X"]

    # 定义 target function

    logh_A = pnorm(b_I + b_T + (b_X + b_TX) * X_m)
    logh_B = pnorm(b_I + b_X * X_m)
    pi_m = logh_A / (logh_A + logh_B)
  }

  X_m = ptnow.X
  Z_m = ptnow.Z
  pts_before = ptsb
  pts_beforeX = pts_before[pts_before[, "X"] == X_m, ]

  Overall_NA_m = sum(pts_beforeX[, "Treat"] == 1) / nrow(pts_beforeX)

  Z_names = paste0("Z",1:ncol(ptsb.Z))

  # 筛选符合条件的样本
  pts_before_Z = pts_beforeX[Z_names, ]

  Z_NA=get_allocation_by_margin_and_stratum(data=pts_beforeX,Z_cols = Z_names,Z_values = Z_m)
  Z_NA_Alloc=Z_NA$Allocation

  NA_m = sum(c(Overall_NA_m , Z_NA_Alloc) * weight)
  # NA_m=sum(pts_before[,"Treat"]==1)/(m-1)
  NB_m = 1 - NA_m
  term1 = (pi_m / NA_m)^v
  term2 = ((1 - pi_m) / NB_m)^v
  # Compute phi_m+1
  phi_m1 = (pi_m * term1) / (pi_m * term1 + (1 - pi_m) * term2)
  return(list(prob = phi_m1))

}




#' Allocation Function of Weighted Balance Ratio Design for Survival Response
#' @importFrom stats binomial coef glm lm pnorm predict qchisq rnorm runif vcov rexp
#' @importFrom survival Surv coxph survreg
#' @param ptsb.X a vector of length \code{n} of the predictive covariates of previous patients. Must be binary.
#' @param ptsb.Z a \code{n x k}  of the prognostic covariates of previous patients. Must be binary.
#' @param ptsb.t a vector of length \code{n} of the treatment allocation of previous patients.
#' @param ptsb.Y a vector of length \code{n} of the responses of previous patients.
#' @param ptsb.E a vector of length \code{n} with value 1 or 0 of the status of event and censoring.
#' @param ptnow.X a binary value of the predictive covariate of the present patient.
#' @param ptnow.Z a vector of length \code{k} of the binary prognostic covariate of the present patient.
#' @param weight a vector of length \code{2+k}. The weight of balance ratio in overall,margin and stratum levels.
#' @param v a positive value that controls the randomness of allocation probability function.
#' @description Calculating the probability of assigning the upcoming patient to treatment A based on the patient's covariates and the previous patients' covariates and responses for Zhao's New procedure.
#' @return \item{prob}{Probability of assigning the upcoming patient to treatment A.}
#' @export
#' @concept Weighted Balance Ratio Design
#' @examples
#' set.seed(123)
#'
#' # Generate historical data for 400 patients
#' ptsb.X = sample(c(1, -1), 400, replace = TRUE)  # predictive covariate
#' ptsb.Z = cbind(
#'   sample(c(1, -1), 400, replace = TRUE),         # prognostic covariate 1
#'   sample(c(1, -1), 400, replace = TRUE)          # prognostic covariate 2
#' )
#' ptsb.Y = rexp(400, rate = 1)                    # observed time
#' ptsb.E = sample(c(1, 0), 400, replace = TRUE)   # event indicator (1=event, 0=censored)
#' ptsb.t = sample(c(1, 0), 400, replace = TRUE)   # treatment assignments
#'
#' # Incoming patient covariates
#' ptnow.X = 1
#' ptnow.Z = c(1, -1)
#'
#' # Calculate allocation probability using WBR for survival response
#' result = WBR_Alloc_Surv(
#'   ptsb.X = ptsb.X,
#'   ptsb.Z = ptsb.Z,
#'   ptsb.Y = ptsb.Y,
#'   ptsb.E = ptsb.E,
#'   ptsb.t = ptsb.t,
#'   ptnow.X = ptnow.X,
#'   ptnow.Z = ptnow.Z,
#'   weight = rep(0.25, 4)
#' )
#'
#' # View probability of assigning to treatment A
#' result$prob
WBR_Alloc_Surv = function(
                        ptsb.X,
                        ptsb.Z,
                        ptsb.Y,
                        ptsb.t,
                        ptsb.E,
                        ptnow.X,
                        ptnow.Z,
                        v=2,
                        weight) {
  # ptsb.X
  if (!is.numeric(ptsb.X) || !is.vector(ptsb.X)) {
    stop("ptsb.X must be a numeric vector.")
  }
  if (length(unique(ptsb.X)) != 2) {
    stop("ptsb.X must be binary (i.e., contain exactly two unique values).")
  }

  n = length(ptsb.X)

  # ptsb.Z
  if (!is.matrix(ptsb.Z)) {
    stop("ptsb.Z must be a matrix.")
  }
  if (nrow(ptsb.Z) != n) {
    stop("ptsb.Z must have the same number of rows as ptsb.X.")
  }
  k = ncol(ptsb.Z)
  if (!all(apply(ptsb.Z, 2, function(col) length(unique(col)) == 2))) {
    stop("Each column of ptsb.Z must be binary (i.e., exactly two unique values).")
  }

  # ptsb.Y
  if (!is.numeric(ptsb.Y) || length(ptsb.Y) != n) {
    stop("ptsb.Y must be a numeric vector of the same length as ptsb.X.")
  }

  # ptsb.t
  if (!is.numeric(ptsb.t) || length(ptsb.t) != n) {
    stop("ptsb.t must be a numeric vector of the same length as ptsb.X.")
  }

  # ptsb.E
  if (!is.numeric(ptsb.E) || length(ptsb.E) != n) {
    stop("ptsb.E must be a numeric vector of the same length as ptsb.X.")
  }
  if (!all(ptsb.E %in% c(0, 1))) {
    stop("ptsb.E must be binary (values 0 or 1).")
  }

  # ptnow.X
  if (!(ptnow.X %in% unique(ptsb.X))) {
    stop("ptnow.X must be one of the values in ptsb.X.")
  }

  # ptnow.Z
  if (!is.numeric(ptnow.Z) || length(ptnow.Z) != k) {
    stop(paste0("ptnow.Z must be a numeric vector of length ", k, "."))
  }
  for (j in seq_len(k)) {
    if (!(ptnow.Z[j] %in% unique(ptsb.Z[, j]))) {
      stop(paste0("ptnow.Z[", j, "] must match one of the values in ptsb.Z[,", j, "]."))
    }
  }

  # weight
  if (!is.numeric(weight) || length(weight) != (2 + k)) {
    stop(paste0("weight must be a numeric vector of length ", 2 + k, "."))
  }

  # v
  if (!is.numeric(v) || length(v) != 1 || v <= 0) {
    stop("v must be a single positive numeric value.")
  }

  ptsb = cbind.data.frame(ptsb.X, ptsb.Z, ptsb.t, ptsb.Y, ptsb.E)
  colnames(ptsb) = c("X", paste0("Z", 1:ncol(ptsb.Z)), "Treat", "Y","E")

  Z_terms = paste0("Z", 1:ncol(ptsb.Z))
  X_m = ptnow.X
  rhs = paste(c("Treat", "X", "X:Treat", Z_terms), collapse = " + ")
  form = as.formula(paste("Surv(Y,E) ~", rhs))

    fit = coxph(formula = form, data = data.frame(ptsb))

    b_T  = coef(fit)["Treat"]
    b_X  = coef(fit)["X"]
    b_TX = coef(fit)["Treat:X"]

    # 定义 target function

    logh_A = exp(b_T + (b_X + b_TX) * X_m)
    logh_B = exp(b_X * X_m)
    pi_m = logh_A / (logh_A + logh_B)



  Z_m = ptnow.Z
  pts_before = ptsb
  pts_beforeX = pts_before[pts_before[, "X"] == X_m, ]

  Overall_NA_m = sum(pts_beforeX[, "Treat"] == 1) / nrow(pts_beforeX)

  Z_names = paste0("Z",1:ncol(ptsb.Z))

  # 筛选符合条件的样本
  pts_before_Z = pts_beforeX[Z_names, ]

  Z_NA=get_allocation_by_margin_and_stratum(data=pts_beforeX,Z_cols = Z_names,Z_values = Z_m)
  Z_NA_Alloc=Z_NA$Allocation

  NA_m = sum(c(Overall_NA_m , Z_NA_Alloc) * weight)
  # NA_m=sum(pts_before[,"Treat"]==1)/(m-1)
  NB_m = 1 - NA_m
  term1 = (pi_m / NA_m)^v
  term2 = ((1 - pi_m) / NB_m)^v
  # Compute phi_m+1
  phi_m1 = (pi_m * term1) / (pi_m * term1 + (1 - pi_m) * term2)
  return(list(prob = as.numeric(phi_m1)))
}



#' Simulation Function of Weighted Balance Ratio Design for Binary and Continuous Response
#' @importFrom stats binomial coef glm lm pnorm predict qchisq rnorm runif vcov rexp
#' @param n a positive integer. The sample size of the simulated data.
#' @param mu a vector of length 2. The true parameters of treatment.
#' @param beta a vector of length 2. The true parameters of predictive covariate and interaction with treatment.
#' @param gamma a vector of length k. The true parameters of prognostic covariates.
#' @param m0 a positive integer. The number of first 2m0 patients will be allocated equally to both treatments.
#' @param pts.X a vector of length n. The vector of patients' binary predictive covariates.
#' @param pts.Z a matrix of \code{n x k}. The matrix of patients' binary prognostic covariates.
#' @param response the type of the response. Options are \code{"Binary"} or \code{"Cont"}.
#' @param weight a vector of length \code{2+k}. The weight of balance ratio in overall,margin and stratum levels.
#' @param v a positive value that controls the randomness of allocation probability function.
#' @description This function simulates a trial using Weighted Balance Ratio design for binary and continuous responses.
#' @return
#' \item{method}{The name of procedure.}
#' \item{sampleSize}{The sample size of the trial.}
#' \item{assignment}{The randomization sequence.}
#' \item{X1proportion}{Average allocation proportion for treatment A when predictive covariate equals the smaller value.}
#' \item{X2proportion}{Average allocation proportion for treatment A when predictive covariate equals the larger value.}
#' \item{proportion}{Average allocation proportion for treatment A.}
#' \item{failureRate}{Proportion of treatment failures (if \code{response = "Binary"}).}
#' \item{meanResponse}{Mean response value (if \code{response = "Cont"}).}#'
#' \item{rejectNull}{Logical. Indicates whether the treatment effect is statistically significant based on a Wald test.}

#' @export
#' @concept Weighted Balance Ratio Design
#' @examples
#' set.seed(123)
#'
#' # Simulation settings
#' n = 400                              # total number of patients
#' mu = c(0.8, 0.8)                     # treatment effects (muA, muB)
#' beta = c(0.8, -0.8)                  # predictive effect and interaction
#' gamma = c(0.8, 0.8)                  # prognostic covariate effects
#' weight = rep(0.25, 4)               # weights for imbalance components
#'
#' # Generate patient covariates
#' pts.X = sample(c(1, -1), n, replace = TRUE)  # predictive covariate
#' pts.Z = cbind(
#'   sample(c(1, -1), n, replace = TRUE),        # prognostic Z1
#'   sample(c(1, -1), n, replace = TRUE)         # prognostic Z2
#' )
#'
#' # Run simulation for continuous response
#' result = WBR_Sim(
#'   n = n,
#'   mu = mu,
#'   beta = beta,
#'   gamma = gamma,
#'   pts.X = pts.X,
#'   pts.Z = pts.Z,
#'   response = "Cont",
#'   weight = weight
#' )
WBR_Sim = function(n,
                        mu,
                        beta,
                        gamma,
                        m0 = 40,
                        pts.X,
                        pts.Z,
                        response,
                        weight,v=2) {
  pts=ZhaoNew_generate_data(
    n = n,
    pts.X = pts.X,
    pts.Z = pts.Z,
    mu = mu,
    beta = beta,
    gamma = gamma,
    m0 = m0,
    response = response
  )
  for (m in (2 * m0 + 1):nrow(pts)) {
    # n 和 m0
    if (!is.numeric(n) || length(n) != 1 || n <= 0 || n != as.integer(n)) {
      stop("n must be a positive integer.")
    }
    if (!is.numeric(m0) || length(m0) != 1 || m0 <= 0 || m0 != as.integer(m0)) {
      stop("m0 must be a positive integer.")
    }
    if (2 * m0 >= n) {
      stop("2 * m0 must be less than n.")
    }

    # mu
    if (!is.numeric(mu) || length(mu) != 2) {
      stop("mu must be a numeric vector of length 2.")
    }

    # beta
    if (!is.numeric(beta) || length(beta) != 2) {
      stop("beta must be a numeric vector of length 2.")
    }

    # pts.X
    if (!is.numeric(pts.X) || !is.vector(pts.X) || length(pts.X) != n) {
      stop("pts.X must be a numeric vector of length n.")
    }
    if (length(unique(pts.X)) != 2) {
      stop("pts.X must be binary (i.e., contain exactly two unique values).")
    }

    # pts.Z
    if (!is.matrix(pts.Z)) {
      stop("pts.Z must be a matrix.")
    }
    if (nrow(pts.Z) != n) {
      stop("pts.Z must have n rows.")
    }
    k = ncol(pts.Z)
    if (!all(apply(pts.Z, 2, function(col) length(unique(col)) == 2))) {
      stop("Each column in pts.Z must be binary (i.e., contain exactly two unique values).")
    }

    # gamma
    if (!is.numeric(gamma) || length(gamma) != k) {
      stop(paste0("gamma must be a numeric vector of length ", k, "."))
    }

    # weight
    if (!is.numeric(weight) || length(weight) != (2 + k)) {
      stop(paste0("weight must be a numeric vector of length ", 2 + k, "."))
    }

    # v
    if (!is.numeric(v) || length(v) != 1 || v <= 0) {
      stop("v must be a single positive numeric value.")
    }

    # response
    if (!(response %in% c("Binary", "Cont"))) {
      stop('response must be either "Binary" or "Cont".')
    }

  ptsb=pts[1:(m-1),]
  X_m = as.numeric(pts[m,"X"])
  colnames(pts) = c("X", paste0("Z", 1:ncol(pts.Z)), "Treat", "Y")

  Z_terms = paste0("Z", 1:ncol(pts.Z))
  rhs = paste(c("Treat", "X", "X:Treat", Z_terms), collapse = " + ")
  form = as.formula(paste("Y ~", rhs))

  if (response == "Binary") {
    p1hat = sum(ptsb[, "Y"] == 1 &
                  ptsb[, "X"] == X_m &
                  ptsb[, "Treat"] == 1) / sum(ptsb[, "X"] == X_m &
                                                ptsb[, "Treat"] == 1)
    p0hat = sum(ptsb[, "Y"] == 1 &
                  ptsb[, "X"] == X_m &
                  ptsb[, "Treat"] == 0) / sum(ptsb[, "X"] == X_m  &
                                                ptsb[, "Treat"] == 0)
    pi_m = sqrt(p1hat) / (sqrt(p1hat) + sqrt(p0hat))
  }

  else if (response == "Cont") {
    fit = lm(formula = form, data = data.frame(ptsb))
    b_I = coef(fit)["(Intercept)"]
    b_T  = coef(fit)["Treat"]
    b_X  = coef(fit)["X"]
    b_TX = coef(fit)["Treat:X"]


    logh_A = pnorm(b_I + b_T + (b_X + b_TX) * X_m)
    logh_B = pnorm(b_I + b_X * X_m)
    pi_m = logh_A / (logh_A + logh_B)
  }


  Z_names = paste0("Z",1:ncol(pts.Z))

  Z_m = as.numeric(pts[m, Z_names])
  pts_before = ptsb
  pts_beforeX = pts_before[pts_before[, "X"] == X_m, ]

  Overall_NA_m = sum(pts_beforeX[, "Treat"] == 1) / nrow(pts_beforeX)

  pts_before_Z = pts_beforeX[,Z_names]

  Z_NA=get_allocation_by_margin_and_stratum(data=pts_beforeX,Z_cols = Z_names,Z_values = Z_m)
  Z_NA_Alloc=Z_NA$Allocation

  NA_m = sum(c(Overall_NA_m , Z_NA_Alloc) * weight)
  # NA_m=sum(pts_before[,"Treat"]==1)/(m-1)
  NB_m = 1 - NA_m
  term1 = (pi_m / NA_m)^v
  term2 = ((1 - pi_m) / NB_m)^v
  # Compute phi_m+1
  phi_m1 = (pi_m * term1) / (pi_m * term1 + (1 - pi_m) * term2)

  if (is.na(phi_m1)){pts[m,"Treat"]=sample(c(1,0),1,prob=c(0.5,0.5))
  }
  else{
  pts[m,"Treat"]=sample(c(1,0),1,prob=c(phi_m1,1-phi_m1))
  }

  if (response == "Binary") {
    term.m = c(1,pts[m, "Treat"]) %*% mu +
      c(pts[m, "X"], pts[m, "X"] * pts[m, "Treat"]) %*% beta +
      pts.Z[m, ] %*% gamma
    prob.Y0.m = 1 / (1 + exp(-term.m))
    pts[m, "Y"] = as.numeric(runif(1) < prob.Y0.m)
  }
  else if (response == "Cont") {
    pts[m, "Y"] = c(1,pts[m, "Treat"]) %*% mu +
      c(pts[m, "X"], pts[m, "X"] * pts[m, "Treat"]) %*% beta +
      pts.Z[m, ] %*% gamma + rnorm(1)
  }
  }
  return(ZhaoNew_Output(
    name = "Weighted Balace Ratio Design",
    pts = pts,
    response = response
  ))
}


#' Simulation Function of Weighted Balance Ratio Design for Survival Response
#' @importFrom stats binomial coef glm lm pnorm predict qchisq rnorm runif vcov rexp
#' @importFrom survival Surv coxph survreg
#' @param n a number. The sample size of the simulated data.
#' @param mu a number. The true parameters of treatment effect.
#' @param beta a vector of length 2. The true parameters of predictive covariate and interaction with treatment.
#' @param gamma a vector of length k. The true parameters of prognostic covariates.
#' @param m0 a positive integer. The number of first 2m0 patients will be allocated equally to both treatments.
#' @param pts.X a vector of length n. The vector of patients' binary predictive covariates.
#' @param pts.Z a matrix of \code{n x k}. The matrix of patients' binary prognostic covariates.
#' @param censor.time a positive number. The upper bound of the uniform censor time in year.
#' @param arrival.rate a positive integer. The arrival rate of patients each year.
#' @param weight a vector of length \code{2+k}. The weight of balance ratio in overall,margin and stratum levels.
#' @param v a positive value that controls the randomness of allocation probability function.
#' @description This function simulates a trial using Weighted Balance Ratio design for survival responses.
#' @concept Weighted Balance Ratio Design
#' @return A list with the following elements:
#' \item{method}{The name of procedure.}
#' \item{sampleSize}{The sample size of the trial.}
#' \item{assignment}{The randomization sequence.}
#' \item{X1proportion}{Average allocation proportion for treatment A when predictive covariate equals the smaller value.}
#' \item{X2proportion}{Average allocation proportion for treatment A when predictive covariate equals the larger value.}
#' \item{proportion}{Average allocation proportion for treatment A.}
#' \item{N.events}{Total number of events occured of the trial.}
#' \item{responses}{Observed survival responses of patients.}
#' \item{events}{Survival status vector of patients(1=event,0=censored)}
#' \item{rejectNull}{Logical. Indicates whether the treatment effect is statistically significant based on a Wald test.}

#' @export
#' @examples
#' set.seed(123)
#'
#' # Simulation settings
#' n = 400                            # total number of patients
#' mu = 0.5                           # treatment effect (log hazard ratio)
#' beta = c(0.5, -0.5)                # predictive effect and interaction
#' gamma = c(0.5, 0.5)                # prognostic covariate effects
#' censor.time = 2                   # maximum censoring time (years)
#' arrival.rate = 1.5                # arrival rate per year
#' weight = rep(0.25, 4)             # imbalance weights for overall, margins, and stratum
#'
#' # Generate patient covariates
#' pts.X = sample(c(1, -1), n, replace = TRUE)  # predictive covariate
#' pts.Z = cbind(
#'   sample(c(1, -1), n, replace = TRUE),        # prognostic Z1
#'   sample(c(1, -1), n, replace = TRUE)         # prognostic Z2
#' )
#'
#' # Run simulation for survival outcome
#' result = WBR_Sim_Surv(
#'   n = n,
#'   mu = mu,
#'   beta = beta,
#'   gamma = gamma,
#'   pts.X = pts.X,
#'   pts.Z = pts.Z,
#'   censor.time = censor.time,
#'   arrival.rate = arrival.rate,
#'   weight = weight
#'
#' )
#'
WBR_Sim_Surv = function(n,
                   mu,
                   beta,
                   gamma,
                   m0 = 40,
                   pts.X,
                   pts.Z,
                   censor.time,
                   arrival.rate,
                   weight,v=2) {
  # ---- Input Checks ----

  # n 和 m0
  if (!is.numeric(n) || length(n) != 1 || n <= 0 || n != as.integer(n)) {
    stop("n must be a positive integer.")
  }
  if (!is.numeric(m0) || length(m0) != 1 || m0 <= 0 || m0 != as.integer(m0)) {
    stop("m0 must be a positive integer.")
  }
  if (2 * m0 >= n) {
    stop("2 * m0 must be less than n.")
  }

  # mu
  if (!is.numeric(mu) || length(mu) != 1) {
    stop("mu must be a single numeric value.")
  }

  # beta
  if (!is.numeric(beta) || length(beta) != 2) {
    stop("beta must be a numeric vector of length 2.")
  }

  # pts.X
  if (!is.numeric(pts.X) || length(pts.X) != n) {
    stop("pts.X must be a numeric vector of length n.")
  }
  if (length(unique(pts.X)) != 2) {
    stop("pts.X must be binary (i.e., contain exactly two unique values).")
  }

  # pts.Z
  if (!is.matrix(pts.Z)) {
    stop("pts.Z must be a matrix.")
  }
  if (nrow(pts.Z) != n) {
    stop("pts.Z must have n rows.")
  }
  k = ncol(pts.Z)
  if (!all(apply(pts.Z, 2, function(col) length(unique(col)) == 2))) {
    stop("Each column in pts.Z must be binary (i.e., exactly two unique values).")
  }

  # gamma
  if (!is.numeric(gamma) || length(gamma) != k) {
    stop(paste0("gamma must be a numeric vector of length ", k, "."))
  }

  # weight
  if (!is.numeric(weight) || length(weight) != (2 + k)) {
    stop(paste0("weight must be a numeric vector of length ", 2 + k, "."))
  }

  # v
  if (!is.numeric(v) || length(v) != 1 || v <= 0) {
    stop("v must be a single positive numeric value.")
  }

  # censor.time
  if (!is.numeric(censor.time) || length(censor.time) != 1 || censor.time <= 0) {
    stop("censor.time must be a single positive number.")
  }

  # arrival.rate
  if (!is.numeric(arrival.rate) || length(arrival.rate) != 1 || arrival.rate <= 0) {
    stop("arrival.rate must be a single positive number.")
  }


  pts=ZhaoNew_generate_data_Surv(
    n = n,
    pts.X = pts.X,
    pts.Z = pts.Z,
    mu = mu,
    beta = beta,
    gamma = gamma,
    m0 = m0,
    censor.time=censor.time,
    arrival.rate=arrival.rate
  )
  for (m in (2 * m0 + 1):nrow(pts)) {
    ptsb=pts[1:(m-1),]
    X_m = as.numeric(pts[m,"X"])

    current_time=pts[m,"A"]
    pts_before = pts[1:(m - 1), ]

    Ai=current_time-pts_before[,"A"]

    E_star=as.numeric(pts_before[,"Y"]<Ai)

    pts_before_star=cbind(pts_before,E_star)

    pts_before_star[,"E_star"][pts_before_star[,"Y"]==pts_before_star[,"C"]]=0

    Y_star=pmin(pts_before_star[,"Y"],Ai)

    pts_before_star=cbind(pts_before_star,Y_star)

    pts_before_Xm=pts_before_star[pts_before_star[,"X"]==X_m,]

    Z_terms = grep("^Z", colnames(pts_before), value = TRUE)
    rhs = paste(c("Treat", "X", "X:Treat", Z_terms), collapse = " + ")
    form = as.formula(paste("Surv(Y_star, E_star) ~", rhs))

    fit = coxph(formula = form, data = data.frame(pts_before_star))

    b_T  = -coef(fit)["Treat"]
    b_X  = -coef(fit)["X"]
    b_TX = -coef(fit)["Treat:X"]

    logh_A = exp(b_T + (b_X + b_TX) * X_m)
    logh_B = exp(b_X * X_m)

    pi_m=logh_A/(logh_A+logh_B)

    Z_names = paste0("Z",1:ncol(pts.Z))

    Z_m = as.numeric(pts[m, Z_names])
    pts_before = ptsb
    pts_beforeX = pts_before[pts_before[, "X"] == X_m, ]

    Overall_NA_m = sum(pts_beforeX[, "Treat"] == 1) / nrow(pts_beforeX)

    pts_before_Z = pts_beforeX[,Z_names]

    Z_NA=get_allocation_by_margin_and_stratum(data=pts_beforeX,Z_cols = Z_names,Z_values = Z_m)
    Z_NA_Alloc=Z_NA$Allocation

    NA_m = sum(c(Overall_NA_m , Z_NA_Alloc) * weight)

    NB_m = 1 - NA_m
    term1 = (pi_m / NA_m)^v
    term2 = ((1 - pi_m) / NB_m)^v

    phi_m1 = (pi_m * term1) / (pi_m * term1 + (1 - pi_m) * term2)

    if (is.na(phi_m1)){pts[m,"Treat"]=sample(c(1,0),1,prob=c(0.5,0.5))
    }
    else{
      pts[m,"Treat"]=sample(c(1,0),1,prob=c(phi_m1,1-phi_m1))
    }

    surv.rate.m=exp(c(pts[m, "Treat"]) * mu +
                      c(pts[m, "X"], pts[m, "X"] * pts[m, "Treat"]) %*% beta +
                      pts.Z[m, ] %*% gamma)

    pts[m,"S"]=rexp(1,rate=1/surv.rate.m)
    pts[m,"E"]=as.numeric(pts[m,"S"]<pts[m,"C"])
    pts[m,"Y"]=min(pts[m,"S"],pts[m,"C"])
  }
  return(ZhaoNew_Output_Surv(
    name = "Weighted Balace Ratio Design",
    pts = pts
  ))
}
