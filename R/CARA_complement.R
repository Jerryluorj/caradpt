#Function for equally allocate first 2m0 patients to each treatment.
permuted_block = function(n, block_size = 4) {
  assignments = numeric(n)
  ts = c(1, 0)

  num_full_blocks = n %/% block_size
  remaining_patients = n %% block_size

  if (remaining_patients == 1 || remaining_patients > 2) {
    stop("Remaining patients must be even or zero.")
  }

  for (i in 1:num_full_blocks) {
    block = rep(ts, each = block_size / length(ts))
    assignments[((i - 1) * block_size + 1):(i * block_size)] = sample(block)
  }


  if (remaining_patients == 2) {
    last_block = sample(ts)
    assignments[(n - remaining_patients + 1):n] = last_block
  }

  return(as.numeric(assignments))
}

#Calculate the estimated target for a patient.
pi.patient = function(pt.cov, fitA, fitB, response, target) {
  if ((target != "Neyman") & (target != "RSIHR")) {
    stop("Your target allocation is not in the list")
  }
  if (response == "Binary") {
    theta_hatA = coef(fitA)
    theta_hatB = coef(fitB)
    tempA = exp(-sum(theta_hatA * c(1, pt.cov)))
    tempB = exp(-sum(theta_hatB * c(1, pt.cov)))

    pA = 1 / (1 + tempA)
    pB = 1 / (1 + tempB)

    if (target == "Neyman") {
      rho = sqrt(pA * (1 - pA)) / (sqrt(pA * (1 - pA)) + sqrt(pB * (1 - pB)))
    } else if (target == "RSIHR") {
      rho = sqrt(pA) / (sqrt(pA) + sqrt(pB))
    }

  }
  else if (response == "Cont") {
    theta_hatA = coef(fitA)
    theta_hatB = coef(fitB)

    ztilde = c(1, pt.cov)

    muA = sum(theta_hatA * ztilde)
    muB = sum(theta_hatB * ztilde)

    sigmaA = summary(fitA)$sigma
    sigmaB = summary(fitB)$sigma

    if (target == "Neyman") {
      # rho= sigmaA / (sigmaA + sigmaB)
      rho <- (sigmaA * sqrt(pnorm(muB))) / (sigmaA * sqrt(pnorm(muB)) + sigmaB * sqrt(pnorm(muA)))
    } else if (target == "RSIHR") {
      # if (muA <= 0 || muB <= 0) {
      #   warning("Predicted muA or muB is non-positive. Using Neyman allocation instead.")
      #   rho <- sigmaA / (sigmaA + sigmaB)
      # } else {
      rho <- (sigmaA * sqrt(pnorm(muB))) / (sigmaA * sqrt(pnorm(muB)) + sigmaB * sqrt(pnorm(muA)))
    }
  }


  return(rho)
}

#Calculate the estimated target for a patient with survival response.
pi.patient.surv = function(pt.cov, fitA, fitB, target) {
  if ((target != "Neyman") & (target != "RSIHR")) {
    stop("Your target allocation is not in the list")
  }
  # if (model=="Cox"){
  theta_hatA = -coef(fitA)
  theta_hatB = -coef(fitB)
  varA = fitA$var
  varB = fitB$var
  if (target == "Neyman") {
    tempA = exp(pt.cov %*% theta_hatA) * sqrt(t(pt.cov) %*% pt.cov * (t(pt.cov) %*%
                                                                        varA %*% pt.cov))
    tempB = exp(pt.cov %*% theta_hatB) * sqrt(t(pt.cov) %*% pt.cov * (t(pt.cov) %*%
                                                                        varB %*% pt.cov))
    rho = tempA / (tempA + tempB)
  }
  else if (target == "RSIHR") {
    tempA = exp(pt.cov %*% theta_hatA) * sqrt(t(pt.cov) %*% pt.cov * (t(pt.cov) %*%
                                                                        varA %*% pt.cov) * exp(pt.cov %*% theta_hatB))
    tempB = exp(pt.cov %*% theta_hatB) * sqrt(t(pt.cov) %*% pt.cov * (t(pt.cov) %*%
                                                                        varB %*% pt.cov) * exp(pt.cov %*% theta_hatA))
    rho = tempA / (tempA + tempB)
  }
  return(rho)
  # }
}
#Probability function for CADBCD design.
CADBCD.prob <- function(pi.m, rho.m, m, v, pts) {
  NA.m = sum(pts[1:(m - 1), "Treat"] == 1) / (m - 1)
  NB.m = 1 - NA.m
  term1 = (rho.m / NA.m)^v
  term2 = ((1 - rho.m) / NB.m)^v
  # Compute phi_m+1
  phi.m = (pi.m * term1) / (pi.m * term1 + (1 - pi.m) * term2)
  return(phi.m)
}

#Generate patients' profiles for simulation.
CADBCD_generate_data = function(n, pts.cov, m0 = 40, thetaA, thetaB, response) {
  if (length(thetaA) != ncol(pts.cov) + 1 ||
      length(thetaB) != ncol(pts.cov) + 1) {
    stop("Length of true parameter vector must be equal to the number of covariates+1.")
  }
  Treat = c(permuted_block(2 * m0), rep(NA, n - 2 * m0))
  if (response == "Binary") {
    exp.terms = exp(-Treat[1:(2 * m0)] * cbind(1, pts.cov[1:(2 * m0), , drop = FALSE]) %*% thetaA -
                      (1 - Treat[1:(2 * m0)]) * cbind(1, pts.cov[1:(2 * m0), , drop = FALSE]) %*% thetaB)
    prob.Y0 = 1 / (1 + exp.terms)
    Y = c(as.numeric(runif(2 * m0) < prob.Y0), rep(NA, n - 2 * m0))
  }
  else if (response == "Cont") {
    Y = c(
      Treat[1:(2 * m0)] * cbind(1, pts.cov[1:(2 * m0), , drop = FALSE]) %*% thetaA +
        (1 - Treat[1:(2 * m0)]) * cbind(1, pts.cov[1:(2 * m0), , drop = FALSE]) %*% thetaB +
        rnorm(2 * m0),
      rep(NA, n - 2 * m0)
    )
  }
  colnames(pts.cov) = paste0("Z", 1:ncol(pts.cov))
  patients = cbind(pts.cov, Treat, Y)
  colnames(patients) = c(paste0("Z", 1:ncol(pts.cov)), "Treat", "Y")
  return(patients)
}


CADBCD_generate_data_surv = function(n,
                                     pts.cov,
                                     m0 = 40,
                                     thetaA,
                                     thetaB,
                                     model,
                                     censor.time,
                                     arrival.rate) {
  # if (length(thetaA)!=ncol(pts.cov)+1 || length(thetaB)!=ncol(pts.cov)+1){
  #   stop("Length of true parameter vector must be equal to the number of covariates+1.")
  # }
  Treat = c(permuted_block(2 * m0), rep(NA, n - 2 * m0))
  colnames(pts.cov) = paste0("Z", 1:ncol(pts.cov))
  # U= runif(2*m0)
  linpred = exp(Treat[1:(2 * m0)] * pts.cov[1:(2 * m0), ] %*% thetaA +
                  (1 - Treat)[1:(2 * m0)] * pts.cov[1:(2 * m0), ] %*% thetaB)
  S = c(rexp(2 * m0, rate = 1 / linpred), rep(NA, n - 2 * m0))
  inter_arrival = rexp(n, rate = arrival.rate)
  A = cumsum(inter_arrival)
  C = runif(n, 0, censor.time)
  E = as.numeric(S < C)
  Y = pmin(S, C)
  patients = cbind(pts.cov, Treat, S, C, A, Y, E)
  return(patients)
}

ZhaoNew_generate_data = function(n,
                                 pts.X,
                                 pts.Z,
                                 mu,
                                 beta,
                                 gamma,
                                 m0 = 40,
                                 response) {
  Treat = c(permuted_block(2 * m0), rep(NA, n - 2 * m0))
  if (response == "Binary") {
    term = cbind(Treat[1:(2 * m0)], (1 - Treat[1:(2 * m0)])) %*% mu +
      cbind(pts.X[1:(2 * m0)], pts.X[1:(2 * m0)] * Treat[1:(2 * m0)]) %*%
      beta +
      pts.Z[1:(2 * m0), ] %*% gamma
    prob.Y0 = 1 / (1 + exp(-term))
    Y = c(as.numeric(runif(2 * m0) < prob.Y0), rep(NA, n - 2 * m0))
  }
  else if (response == "Cont") {
    Y = c(
      cbind(Treat[1:(2 * m0)], (1 - Treat[1:(2 * m0)])) %*% mu +
        cbind(pts.X[1:(2 * m0)], pts.X[1:(2 * m0)] * Treat[1:(2 * m0)]) %*%
        beta +
        pts.Z[1:(2 * m0), ] %*% gamma + rnorm(2 * m0),
      rep(NA, n - 2 * m0)
    )
  }
  colnames(pts.Z) = paste0("Z", 1:ncol(pts.Z))
  pts = cbind.data.frame(X = pts.X, pts.Z, Treat, Y)

  return(pts)
}


ZhaoNew_generate_data_Surv = function(n,
                                      pts.X,
                                      pts.Z,
                                      mu,
                                      beta,
                                      gamma,
                                      m0 = 40,
                                      censor.time,
                                      arrival.rate) {
  Treat = c(permuted_block(2 * m0), rep(NA, n - 2 * m0))
  surv.rate = exp(Treat[1:(2 * m0)] * mu +
    cbind(pts.X[1:(2 * m0)], pts.X[1:(2 * m0)] * Treat[1:(2 * m0)]) %*%
    beta + pts.Z[1:(2 * m0), ] %*% gamma)
  inter_arrival = rexp(n, rate = arrival.rate)
  A = cumsum(inter_arrival)
  S = c(rexp(2 * m0, rate = 1 / surv.rate), rep(NA, n - 2 * m0))
  C = runif(n, 0, censor.time)
  E = c(as.numeric(S[1:(2 * m0)] <= C[1:(2 * m0)]), rep(NA, n - 2 * m0))
  Y = pmin(S, C)
  patients = cbind(pts.X, pts.Z, Treat, S, C, Y, A, E)
  colnames(patients)=c("X",paste0("Z",1:ncol(pts.Z)),"Treat","S","C","Y","A","E")
  return(patients)
}

get_event_prob_from_cox <- function(fit, z_new, c_max) {
  linpred <- sum(z_new * coef(fit))
  lambda <- exp(-linpred)  # because mean = exp(z %*% beta), so hazard = 1/mean = exp(-z^T beta)
  p_event <- 1 - (1 - exp(-c_max * lambda)) / (c_max * lambda)
  return(p_event)
}



#Estimate parameters for simulated data.
para_est = function(pts, response) {
  ptsA = data.frame(pts[pts[, "Treat"] == 1, ])
  ptsB = data.frame(pts[pts[, "Treat"] == 0, ])
  if (response == "Binary") {
    fit_paraA = glm(Y ~ ., data = ptsA[, !(names(ptsA) %in% "Treat")], family = binomial)
    fit_paraB = glm(Y ~ ., data = ptsB[, !(names(ptsB) %in% "Treat")], family = binomial)
  }
  if (response == "Cont") {
    fit_paraA = lm(Y ~ ., data = ptsA[, !(names(ptsA) %in% "Treat")])
    fit_paraB = lm(Y ~ ., data = ptsB[, !(names(ptsB) %in% "Treat")])
  }
  return(list(
    "Theta A" = coef(fit_paraA),
    "Theta B" = coef(fit_paraB)
  ))
}

#Calculate if the treatment effect difference is deteteced for simulated data.
power_est = function(pts, response, alpha = 0.05) {
  ptsA = data.frame(pts[pts[, "Treat"] == 1, ])
  ptsB = data.frame(pts[pts[, "Treat"] == 0, ])
  if (response == "Binary") {
    fit_paraA = glm(Y ~ ., data = ptsA[, !(names(ptsA) %in% "Treat")], family = binomial)
    fit_paraB = glm(Y ~ ., data = ptsB[, !(names(ptsB) %in% "Treat")], family = binomial)
  }
  else if (response == "Cont") {
    fit_paraA = lm(Y ~ ., data = ptsA[, !(names(ptsA) %in% "Treat")])
    fit_paraB = lm(Y ~ ., data = ptsB[, !(names(ptsB) %in% "Treat")])
  }
  wald.stat = t(coef(fit_paraA) - coef(fit_paraB)) %*% solve(vcov(fit_paraA) +
                                                               vcov(fit_paraB)) %*% (coef(fit_paraA) - coef(fit_paraB))
  return(wald.stat > qchisq(1 - alpha, df = length(coef(fit_paraA))))
}




para_est_surv = function(pts) {
  ptsA = data.frame(pts[pts[, "Treat"] == 1, ])
  ptsB = data.frame(pts[pts[, "Treat"] == 0, ])

  Z_vars <- grep("^Z", colnames(pts), value = TRUE)
  formula_str <- paste("Surv(Y, E) ~", paste(Z_vars, collapse = " + "))

  fit_paraA = coxph(as.formula(formula_str), data = ptsA)
  fit_paraB = coxph(as.formula(formula_str), data = ptsB)

  return(list(
    "Theta A" = -coef(fit_paraA),
    "Theta B" = -coef(fit_paraB)
  ))
}

#Calculate if the treatment effect difference is deteteced for simulated data.
power_est_surv = function(pts, alpha = 0.05) {
  ptsA = data.frame(pts[pts[, "Treat"] == 1, ])
  ptsB = data.frame(pts[pts[, "Treat"] == 0, ])

  Z_vars <- grep("^Z", colnames(pts), value = TRUE)
  formula_str <- paste("Surv(Y, E) ~", paste(Z_vars, collapse = " + "))

  fit_paraA = coxph(as.formula(formula_str), data = ptsA)
  fit_paraB = coxph(as.formula(formula_str), data = ptsB)

  wald.stat = t(coef(fit_paraA) - coef(fit_paraB)) %*% solve(fit_paraA$var +
                                                               fit_paraB$var) %*% (coef(fit_paraA) - coef(fit_paraB))
  return(wald.stat > qchisq(1 - alpha, df = length(coef(fit_paraA))))
}

#Output function
CARA_Output = function(name, pts, response, alpha = 0.05) {
  # statistics
  Mean.proportion = mean(pts[, "Treat"])

  para = para_est(pts = pts, response = response)
  reject = power_est(pts = pts,
                     response = response,
                     alpha = alpha)

  # parameter
  # paraNum = paste0('theta',  LETTERS[1:k])
  # names(parameter) = paraNum

  # proportion
  # trt = stringr::str_c("treatment ", LETTERS[1:k])
  # proportion = as.data.frame(proportion)
  # colnames(proportion) = trt


  Mean.response = mean(pts[, "Y"])
  # output
  if (response == "Binary") {
    outputList = list(
      method = name,
      "sampleSize" = nrow(pts),
      "parameter" = para,
      "proportion" = Mean.proportion,
      "assignments" = pts[, "Treat"],


      # "sd of proportion" = sdprop,
      "responses" = pts[, "Y"],
      "failureRate" = Mean.response,
      # "sd of failure rate" = sdfrt,
      'rejectNull' = reject
      # "data: failureRate" = failRate,
      # "data: test" = pwCalc,
      # "data: assignment" = assignment,
      # "data: proportion" = proportion)
    )
  }
  else if (response == "Cont") {
    outputList = list(
      method = name,
      "sampleSize" = nrow(pts),
      "parameter" = para,
      "proportion" = Mean.proportion,
      "assignments" = pts[, "Treat"],

      # "sd of proportion" = sdprop,
      "responses" = pts[, "Y"],
      "meanResponse" = Mean.response,
      # "sd of failure rate" = sdfrt,
      'rejectNull' = reject
      # "data: failureRate" = failRate,
      # "data: test" = pwCalc,
      # "data: assignment" = assignment,
      # "data: proportion" = proportion
    )
  }
  # cat(stringr::str_c(trt, round(Mean.proportion, 3), sep = ": ", collapse = "   "), "\n",
  #     sprintf('SD of allocation = %g', round(sdprop, 3)), "\n","\n")
  # cat(sprintf('Failure Rate = %g', round(Mean.failRate, 3)), "\n",
  #     sprintf('SD of Failure Rate = %g', round(sdfrt, 3)), "\n", "\n")
  # cat(sprintf('Power = %g', round(power, 3)), "\n")
  return(invisible(outputList))
}

CARA_Output_Surv = function(name, pts, alpha = 0.05) {
  # statistics
  Mean.proportion = mean(pts[, "Treat"])

  para = para_est_surv(pts = pts)
  reject = power_est_surv(pts = pts, alpha = alpha)

  # parameter
  # paraNum = paste0('theta',  LETTERS[1:k])
  # names(parameter) = paraNum

  # proportion
  # trt = stringr::str_c("treatment ", LETTERS[1:k])
  # proportion = as.data.frame(proportion)
  # colnames(proportion) = trt

  ptsA = pts[pts[, "Treat"] == 1, ]
  ptsB = pts[pts[, "Treat"] == 0, ]
  Total.EventsA = mean(ptsA[, "E"])
  Total.EventsB = mean(ptsB[, "E"])
  Total.Events = sum(pts[, "E"])

  outputList = list(
    method = name,
    "sampleSize" = nrow(pts),
    "parameter" = para,
    "N.events" = Total.Events,
    "assignments" = pts[, "Treat"],
    "proportion" = Mean.proportion,

    # "sd of proportion" = sdprop,
    "responses" = pts[, "Y"],
    "events" = pts[, "E"],
    # "sd of failure rate" = sdfrt,
    'rejectNull' = reject
    # "data: failureRate" = failRate,
    # "data: test" = pwCalc,
    # "data: assignment" = assignment,
    # "data: proportion" = proportion)
    # cat(stringr::str_c(trt, round(Mean.proportion, 3), sep = ": ", collapse = "   "), "\n",
    #     sprintf('SD of allocation = %g', round(sdprop, 3)), "\n","\n")
    # cat(sprintf('Failure Rate = %g', round(Mean.failRate, 3)), "\n",
    #     sprintf('SD of Failure Rate = %g', round(sdfrt, 3)), "\n", "\n")
    # cat(sprintf('Power = %g', round(power, 3)), "\n")
  )
  return(invisible(outputList))
}




ZhaoNew_Output = function(name, pts, response, alpha = 0.05) {
  # statistics
  ptsX1 = pts[pts[, "X"] == sort(unique(pts[, "X"]))[1], ]
  ptsXN1 = pts[pts[, "X"] == sort(unique(pts[, "X"]))[2], ]
  X1.proportion=mean(ptsX1[,"Treat"])
  XN1.proportion=mean(ptsXN1[,"Treat"])
  Mean.proportion = mean(pts[, "Treat"])

  # para=para_est(pts=pts,response=response)
  reject = power_est(pts = pts,
                     response = response,
                     alpha = alpha)

  # parameter
  # paraNum = paste0('theta',  LETTERS[1:k])
  # names(parameter) = paraNum

  # proportion
  # trt = stringr::str_c("treatment ", LETTERS[1:k])
  # proportion = as.data.frame(proportion)
  # colnames(proportion) = trt


  Mean.response = mean(pts[, "Y"])
  # output
  if (response == "Binary") {
    outputList = list(
      "method" = name,
      "sampleSize" = nrow(pts),
      # "parameter" = para,
      "X1proportion"=mean(ptsX1[,"Treat"]),
      "X2proportion"=mean(ptsXN1[,"Treat"]),
      "proportion" = Mean.proportion,
      "assignments" = pts[, "Treat"],


      # "sd of proportion" = sdprop,
      "responses" = pts[, "Y"],
      "failureRate" = Mean.response
      # "sd of failure rate" = sdfrt,
      # 'reject null' = reject
      # "data: failureRate" = failRate,
      # "data: test" = pwCalc,
      # "data: assignment" = assignment,
      # "data: proportion" = proportion)
    )
  }
  else if (response == "Cont") {
    outputList = list(
      "method" = name,
      "sampleSize" = nrow(pts),
      # "parameter" = para,
      "X1proportion"=mean(ptsX1[,"Treat"]),
      "X2proportion"=mean(ptsXN1[,"Treat"]),
      "proportion" = Mean.proportion,
      "assignments" = pts[, "Treat"],

      # "sd of proportion" = sdprop,
      "responses" = pts[, "Y"],
      "meanResponse" = Mean.response
      # "sd of failure rate" = sdfrt,
      # 'reject null' = reject
      # "data: failureRate" = failRate,
      # "data: test" = pwCalc,
      # "data: assignment" = assignment,
      # "data: proportion" = proportion
    )
  }
  # cat(stringr::str_c(trt, round(Mean.proportion, 3), sep = ": ", collapse = "   "), "\n",
  #     sprintf('SD of allocation = %g', round(sdprop, 3)), "\n","\n")
  # cat(sprintf('Failure Rate = %g', round(Mean.failRate, 3)), "\n",
  #     sprintf('SD of Failure Rate = %g', round(sdfrt, 3)), "\n", "\n")
  # cat(sprintf('Power = %g', round(power, 3)), "\n")
  return(invisible(outputList))
}


ZhaoNew_Output_Surv = function(name, pts, alpha = 0.05) {
  # statistics
  ptsX1 = pts[pts[, "X"] == sort(unique(pts[, "X"]))[1], ]
  ptsXN1 = pts[pts[, "X"] == sort(unique(pts[, "X"]))[2], ]
  X1.proportion=mean(ptsX1[,"Treat"])
  XN1.proportion=mean(ptsXN1[,"Treat"])
  Mean.proportion = mean(pts[, "Treat"])
  #
  # para = para_est_surv(pts = pts)
  # reject = power_est_surv(pts = pts, alpha = alpha)

  # parameter
  # paraNum = paste0('theta',  LETTERS[1:k])
  # names(parameter) = paraNum

  # proportion
  # trt = stringr::str_c("treatment ", LETTERS[1:k])
  # proportion = as.data.frame(proportion)
  # colnames(proportion) = trt

  ptsA = pts[pts[, "Treat"] == 1, ]
  ptsB = pts[pts[, "Treat"] == 0, ]
  # Total.EventsA = mean(ptsA[, "E"])
  # Total.EventsB = mean(ptsB[, "E"])
  Total.Events = sum(pts[, "E"])

  outputList = list(
    "method" = name,
    "sampleSize" = nrow(pts),
    # "parameter" = para,
    # "eventsA" = Total.EventsA,
    # "eventsB" = Total.EventsB,

    "assignments" = pts[, "Treat"],
    "X1proportion"=X1.proportion,
    "X2proportion"=XN1.proportion,
    "proportion" = Mean.proportion,
    "N.events" = Total.Events,

    # "sd of proportion" = sdprop,
    "responses" = pts[, "Y"],
    "events" = pts[, "E"]
    # "sd of failure rate" = sdfrt,
    # 'reject null' = reject
    # "data: failureRate" = failRate,
    # "data: test" = pwCalc,
    # "data: assignment" = assignment,
    # "data: proportion" = proportion)
    # cat(stringr::str_c(trt, round(Mean.proportion, 3), sep = ": ", collapse = "   "), "\n",
    #     sprintf('SD of allocation = %g', round(sdprop, 3)), "\n","\n")
    # cat(sprintf('Failure Rate = %g', round(Mean.failRate, 3)), "\n",
    #     sprintf('SD of Failure Rate = %g', round(sdfrt, 3)), "\n", "\n")
    # cat(sprintf('Power = %g', round(power, 3)), "\n")
  )
  return(invisible(outputList))
}


get_allocation_by_margin_and_stratum <- function(data, Z_cols, Z_values) {
  result <- list()
  n <- length(Z_cols)

  # margin
  for (i in seq_len(n)) {
    col <- Z_cols[i]
    value <- Z_values[i]
    subset_data <- data[data[, col] == value, ]
    alloc <- if (nrow(subset_data) == 0) NA else sum(subset_data[ , "Treat"] == 1) / nrow(subset_data)
    result[[col]] <- alloc
  }

  # stratum
  condition <- rep(TRUE, nrow(data))
  for (i in seq_len(n)) {
    condition <- condition & (data[, Z_cols[i]] == Z_values[i])
  }
  subset_data <- data[condition, ]
  stratum_name <- paste(Z_cols, collapse = "")
  alloc <- if (nrow(subset_data) == 0) NA else sum(subset_data[ , "Treat"] == 1) / nrow(subset_data)
  result[[stratum_name]] <- alloc

  # 输出为 data.frame
  data.frame(
    Term = names(result),
    Allocation = round(unlist(result), 3),
    row.names = NULL
  )
}
