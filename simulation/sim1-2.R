####################
# Simulation design 1
# Outcome type: continuous
# Estimand: ATE on the difference scale
# Randomization design: simple, rerandomization, or stratified rerandomization with pi = 0.5
# Outcome missing: Yes
# sample size: 400 in total
# Estimators (missing outcomes under MAR): DR-WLS, DML
# Measures: bias, ese (empirical se), ase* (averge of se estimators assuming simple randomization is used), 
#           CP_N (coverage prob based on normal distribution), CP_true (coverage prob based on true asymptotic distribution)
###################
rm(list = ls())
set.seed(123)
library(tidyverse)
library(SuperLearner)
library(foreach)
library(doSNOW)
cl <- makeCluster(4)
registerDoSNOW(cl)
package_list <- c("tidyverse", "SuperLearner")
SLmethods <- c("SL.glm", "SL.rpart", "SL.nnet")
expit <- function(x){1/(1+exp(-x))}

n <- 400
# Sample size - varied
# n <- c(30, 100, 250, 500, 1000)
pi <- 0.5
n_sim <- 100
# Get threshold to vary
# qchisq(c(0.4, 0.2, 0.1, 0.01), df = 2)
tt <- 1
# tt <- qchisq(c(0.4, 0.2, 0.1, 0.01), df = 2)

tictoc::tic()
results <- foreach(iter = 1:n_sim, .combine = cbind, .packages = package_list) %dopar% {
  
  # Generate X1
  X1 <- rnorm(n, 1, 1)
  # Generate Strata (Based on X1, when X1 is present, probability of assignment to S_1 strata is higher)
  S <- rbinom(n, size = 1, prob = 0.4 + 0.2 * (X1 < 1))
  # Generate X2
  X2 <- rnorm(n, 0, 1)
  # Outcomes under Treatment (Y1) and Control (Y0) 
  # Double check in the paper that treatment effect is 4 * X2^2 * S
  Y1 <- 2 * exp(X1) + 4 * X2^2 * S + abs(X2) + rnorm(n, 0, 1)
  Y0 <- 2 * exp(X1) + abs(X2) + rnorm(n, 0, 1) 
  # Is the data missing under Treatment (R1) and Control (R0) - Informative missing?
  R1 <- rbinom(n, size = 1, prob = mean(expit(1.2 + X2 + S)))
  R0 <- rbinom(n, size = 1, prob = mean(expit(0.6 + X2 + S)))
  # R1 <- rbinom(n, size = 1, prob = 0.7)
  # R0 <- rbinom(n, size = 1, prob = 0.7)
  Xr <- cbind(X1, X2) # rerandomization variables
  
  tryCatch({
  
  # simple randomization----------------
  A_s <- rbinom(n, 1, pi)
  d_s <- data.frame(RY=Y1 * R1 * A_s + Y0 * R0 * (1-A_s), R=R1 * A_s + R0 * (1-A_s), A=A_s, X1, X2, S) %>%
    mutate(Y = ifelse(R==0, NA, RY))
  
  # DR-WLS
  propensity.fit <- glm(R ~ A+X1+X2+S, data = d_s, family = "binomial")
  propensity_score <- predict(propensity.fit, type = "response")
  ancova.fit <- lm(Y~ A+X1+X2+S, data = d_s, weights = 1/propensity_score)
  p1 <- predict(ancova.fit, newdata=mutate(d_s, A = 1), type = "response")
  p0 <- predict(ancova.fit, newdata=mutate(d_s, A = 0), type = "response")
  est_dr_s <- mean(p1-p0)
  pA <- predict(ancova.fit, newdata = d_s, type = "response")
  Xmat <- model.matrix(propensity.fit)
  dot_h1 <- dot_h0 <- dot_h <- Xmat
  dot_h1[,"A"] <- 1
  dot_h0[,"A"] <- 0
  dot_e <- diag(propensity_score - propensity_score^2) %*% Xmat
  c_1 <- colMeans(diag((d_s$A - pi)/pi/(1-pi) * (d_s$R/propensity_score - 1)) %*% dot_h) %*%
    solve(t(dot_h) %*% diag(d_s$R/propensity_score) %*% Xmat/n)
  c_2 <- colMeans(dot_h1 - dot_h0) %*%
    solve(t(dot_h) %*% diag(d_s$R/propensity_score) %*% Xmat/n) %*%
    (t(dot_e) %*% diag(d_s$R/propensity_score^2 * (d_s$RY - pA)) %*% Xmat/n) %*%
    solve(t(dot_e) %*% Xmat/n)
  IF_dr_s <- ((d_s$A-pi)/pi/(1-pi) - Xmat %*% t(c_1)) * (d_s$R/propensity_score * (d_s$RY - pA)) - 
    Xmat %*% t(c_2) * (d_s$R - propensity_score) + p1 - p0 - est_dr_s
  V_dr_s <- mean(IF_dr_s^2)
  
  # DML
  K <- 5
  sample_split_index <- sample(rep(1:K, each = n/K), size = n)
  sample_split <- map(1:K, ~which(sample_split_index==.))
  temp <- map(1:K, function(k){
    validation_index <- sample_split[[k]]
    training_index1 <- intersect(setdiff(1:n, validation_index), which(d_s$A==1))
    training_index0 <- intersect(setdiff(1:n, validation_index), which(d_s$A==0))
    training_index1_r1 <- intersect(training_index1, which(!is.na(d_s$Y)))
    training_index0_r1 <- intersect(training_index0, which(!is.na(d_s$Y)))
    kappa_fit1 <- SuperLearner(d_s$R[training_index1], X = d_s[training_index1,c("X1","X2","S")], family = "binomial", SL.library = SLmethods)
    kappa_fit0 <- SuperLearner(d_s$R[training_index0], X = d_s[training_index0,c("X1","X2","S")], family = "binomial", SL.library = SLmethods)
    kappa.pred1 <- predict(kappa_fit1, newdata = d_s[validation_index,c("X1","X2","S")])$pred
    kappa.pred0 <- predict(kappa_fit0, newdata = d_s[validation_index,c("X1","X2","S")])$pred
    eta_fit1 <- SuperLearner(d_s$Y[training_index1_r1], X = d_s[training_index1_r1,c("X1","X2","S")], family = "gaussian", SL.library = SLmethods)
    eta_fit0 <- SuperLearner(d_s$Y[training_index0_r1], X = d_s[training_index0_r1,c("X1","X2","S")], family = "gaussian", SL.library = SLmethods)
    eta.pred1 <- predict(eta_fit1, newdata = d_s[validation_index,c("X1","X2","S")])$pred
    eta.pred0 <- predict(eta_fit0, newdata = d_s[validation_index,c("X1","X2","S")])$pred
    f1 <- d_s$A[validation_index]/pi * d_s$R[validation_index]/kappa.pred1 * (d_s$RY[validation_index] - eta.pred1) + eta.pred1
    f0 <- (1-d_s$A[validation_index])/(1-pi) * d_s$R[validation_index]/kappa.pred0 * (d_s$RY[validation_index] - eta.pred0) + eta.pred0
    data.frame(pred1=f1, pred0=f0, fold=k)
  }) %>% Reduce(rbind,.)
  est_dml_s <- mean(temp$pred1 - temp$pred0)
  V_dml_s <- map_dbl(1:K, function(k){
    var(temp$pred1[temp$fold==k]-temp$pred0[temp$fold==k])
  }) %>% mean
  
  # rerandomization----------------
  Q <- 100
  while (Q >= tt){
    A_star <- rbinom(n, 1, pi)
    I <- colMeans(Xr[A_star == 1,]) - colMeans(Xr[A_star == 0,])
    VI <- n/sum(A_star)/sum(1-A_star) * var(Xr)
    Q <- as.numeric(t(I) %*% solve(VI) %*% I)
  }
  A_re <- A_star
  d_re <- data.frame(RY=Y1 * R1 * A_re + Y0 * R0 * (1-A_re), R=R1 * A_re + R0 * (1-A_re), A=A_re, X1, X2, S) %>%
    mutate(Y = ifelse(R==0, NA, RY))
  
  # DR-WLS
  propensity.fit <- glm(R ~ A+X1+X2+S, data = d_re, family = "binomial")
  propensity_score <- predict(propensity.fit, type = "response")
  ancova.fit <- lm(Y~ A+X1+X2+S, data = d_re, weights = 1/propensity_score)
  p1 <- predict(ancova.fit, newdata=mutate(d_re, A = 1), type = "response")
  p0 <- predict(ancova.fit, newdata=mutate(d_re, A = 0), type = "response")
  est_dr_re <- mean(p1-p0)
  pA <- predict(ancova.fit, d_re, type = "response")
  Xmat <- model.matrix(propensity.fit)
  dot_h1 <- dot_h0 <- dot_h <- Xmat
  dot_h1[,"A"] <- 1
  dot_h0[,"A"] <- 0
  dot_e <- diag(propensity_score - propensity_score^2) %*% Xmat
  c_1 <- colMeans(diag((d_re$A - pi)/pi/(1-pi) * (d_re$R/propensity_score - 1)) %*% dot_h) %*%
    solve(t(Xmat) %*% diag(d_re$R/propensity_score) %*% Xmat/n)
  c_2 <- colMeans(dot_h1 - dot_h0) %*%
    solve(t(Xmat) %*% diag(d_re$R/propensity_score) %*% Xmat/n) %*%
    (t(dot_e) %*% diag(d_re$R/propensity_score^2 * (d_re$RY - pA)) %*% Xmat/n) %*%
    solve(t(dot_e) %*% Xmat/n)
  IF_dr_re <- ((d_re$A-pi)/pi/(1-pi) - Xmat %*% t(c_1)) * (d_re$R/propensity_score * (d_re$RY - pA)) - 
    Xmat %*% t(c_2) * (d_re$R - propensity_score) + p1 - p0 - est_dr_re
  V_dr_re <- mean(IF_dr_re^2)
  VI <- n/sum(A_re)/sum(1-A_re) * var(Xr)
  cIIF <- colMeans(diag((d_re$A-pi)/pi/(1-pi) * as.vector(IF_dr_re)) %*% (Xr - rep(1, n) %*% t(colMeans(Xr))))
  R2_dr_re <- t(cIIF) %*% solve(n * VI) %*% cIIF /V_dr_re
  
  # DML
  K <- 5
  sample_split_index <- sample(rep(1:K, each = n/K), size = n)
  sample_split <- map(1:K, ~which(sample_split_index==.))
  temp <- map(1:K, function(k){
    validation_index <- sample_split[[k]]
    training_index1 <- intersect(setdiff(1:n, validation_index), which(d_re$A==1))
    training_index0 <- intersect(setdiff(1:n, validation_index), which(d_re$A==0))
    training_index1_r1 <- intersect(training_index1, which(!is.na(d_re$Y)))
    training_index0_r1 <- intersect(training_index0, which(!is.na(d_re$Y)))
    kappa_fit1 <- SuperLearner(d_re$R[training_index1], X = d_re[training_index1,c("X1","X2","S")], family = "binomial", SL.library = SLmethods)
    kappa_fit0 <- SuperLearner(d_re$R[training_index0], X = d_re[training_index0,c("X1","X2","S")], family = "binomial", SL.library = SLmethods)
    kappa.pred1 <- predict(kappa_fit1, newdata = d_re[validation_index,c("X1","X2","S")])$pred
    kappa.pred0 <- predict(kappa_fit0, newdata = d_re[validation_index,c("X1","X2","S")])$pred
    eta_fit1 <- SuperLearner(d_re$Y[training_index1_r1], X = d_re[training_index1_r1,c("X1","X2","S")], family = "gaussian", SL.library = SLmethods)
    eta_fit0 <- SuperLearner(d_re$Y[training_index0_r1], X = d_re[training_index0_r1,c("X1","X2","S")], family = "gaussian", SL.library = SLmethods)
    eta.pred1 <- predict(eta_fit1, newdata = d_re[validation_index,c("X1","X2","S")])$pred
    eta.pred0 <- predict(eta_fit0, newdata = d_re[validation_index,c("X1","X2","S")])$pred
    f1 <- d_re$A[validation_index]/pi * d_re$R[validation_index]/kappa.pred1 * (d_re$RY[validation_index] - eta.pred1) + eta.pred1
    f0 <- (1-d_re$A[validation_index])/(1-pi) * d_re$R[validation_index]/kappa.pred0 * (d_re$RY[validation_index] - eta.pred0) + eta.pred0
    mutate(d_re[validation_index,], pred1=f1, pred0=f0, fold=k)
  }) %>% Reduce(rbind,.)
  est_dml_re <- mean(temp$pred1 - temp$pred0)
  V_dml_re <- map_dbl(1:K, function(k){
    var(temp$pred1[temp$fold==k]-temp$pred0[temp$fold==k])
  }) %>% mean
  IF <- as.vector(temp$pred1-temp$pred0-mean(temp$pred1-temp$pred0))
  cIIF <- colMeans(diag((temp$A-pi)/pi/(1-pi) * IF) %*% as.matrix(temp[,c("X1","X2")] - rep(1, n) %*% t(colMeans(temp[,c("X1","X2")]))))
  VI <- n/sum(A_re)/sum(1-A_re) * var(Xr)
  R2_dml_re <- t(cIIF) %*% solve(n * VI) %*% cIIF /V_dml_re

  # stratified rerandomization----------------
  Q <- 100
  while(Q >= tt){
    A_tilde <- rep(0, n)
    for(s in unique(S)){
      A_tilde[sample(which(S==s), size = pi * sum(S==s))] <- 1
    }
    I <- colMeans(Xr[A_tilde == 1,]) - colMeans(Xr[A_tilde == 0,])
    VI <- n/sum(A_tilde)/sum(1-A_tilde) * map(unique(S), function(s){var(Xr[S==s,]) * mean(S==s)}) %>% Reduce(`+`,.)
    Q <- as.numeric(t(I) %*% solve(VI) %*% I)
  }
  A_stre <- A_tilde
  d_stre <- data.frame(RY=Y1 * R1 * A_stre + Y0 * R0 * (1-A_stre), R=R1 * A_stre + R0 * (1-A_stre), A=A_stre, X1, X2, S) %>%
    mutate(Y = ifelse(R==0, NA, RY))
  
  # DR-WLS
  propensity.fit <- glm(R ~ A+X1+X2+S, data = d_stre, family = "binomial")
  propensity_score <- predict(propensity.fit, type = "response")
  ancova.fit <- lm(Y~ A+X1+X2+S, data = d_stre, weights = 1/propensity_score)
  p1 <- predict(ancova.fit, newdata=mutate(d_stre, A = 1), type = "response")
  p0 <- predict(ancova.fit, newdata=mutate(d_stre, A = 0), type = "response")
  est_dr_stre <- mean(p1-p0)
  pA <- predict(ancova.fit, d_stre, type = "response")
  Xmat <- model.matrix(propensity.fit)
  dot_h1 <- dot_h0 <- dot_h <- Xmat
  dot_h1[,"A"] <- 1
  dot_h0[,"A"] <- 0
  dot_e <- diag(propensity_score - propensity_score^2) %*% Xmat
  c_1 <- colMeans(diag((d_stre$A - pi)/pi/(1-pi) * (d_stre$R/propensity_score - 1)) %*% dot_h) %*%
    solve(t(Xmat) %*% diag(d_stre$R/propensity_score) %*% Xmat/n)
  c_2 <- colMeans(dot_h1 - dot_h0) %*%
    solve(t(Xmat) %*% diag(d_stre$R/propensity_score) %*% Xmat/n) %*%
    (t(dot_e) %*% diag(d_stre$R/propensity_score^2 * (d_stre$RY - pA)) %*% Xmat/n) %*%
    solve(t(dot_e) %*% Xmat/n)
  IF_dr_stre <- ((d_stre$A-pi)/pi/(1-pi) - Xmat %*% t(c_1)) * (d_stre$R/propensity_score * (d_stre$RY - pA)) - 
    Xmat %*% t(c_2) * (d_stre$R - propensity_score) + p1 - p0 - est_dr_stre
  V_dr_stre <- mean(IF_dr_stre^2)
  VI <- n/sum(A_stre)/sum(1-A_stre) * var(Xr)
  cIIF <- colMeans(diag((d_stre$A-pi)/pi/(1-pi) * as.vector(IF_dr_stre)) %*% (Xr - rep(1, n) %*% t(colMeans(Xr))))
  R2_dr_stre <- t(cIIF) %*% solve(n * VI) %*% cIIF /V_dr_stre
  
  # DML
  K <- 5
  temp <- map(unique(d_stre$S), function(s){
    d_stre_s <- d_stre[d_stre$S==s,]
    sample_split_index <- rep(0, nrow(d_stre_s))
    for(a in 0:1){
      nn <- sum(d_stre_s$A == a)
      sample_split_index[d_stre_s$A == a] <- sample(rep(1:K, each = ceiling(nn/K)), size = nn)
    }
    sample_split <- map(1:K, ~which(sample_split_index==.))
    temp1 <- map(1:K, function(k){
      validation_index <- sample_split[[k]]
      training_index1 <- intersect(setdiff(1:nrow(d_stre_s), validation_index), which(d_stre_s$A==1))
      training_index0 <- intersect(setdiff(1:nrow(d_stre_s), validation_index), which(d_stre_s$A==0))
      training_index1_r1 <- intersect(training_index1, which(!is.na(d_stre_s$Y)))
      training_index0_r1 <- intersect(training_index0, which(!is.na(d_stre_s$Y)))
      kappa_fit1 <- SuperLearner(d_stre_s$R[training_index1], X = d_stre_s[training_index1,c("X1","X2")], family = "binomial", SL.library = SLmethods)
      kappa_fit0 <- SuperLearner(d_stre_s$R[training_index0], X = d_stre_s[training_index0,c("X1","X2")], family = "binomial", SL.library = SLmethods)
      kappa.pred1 <- predict(kappa_fit1, newdata = d_stre_s[validation_index,c("X1","X2")])$pred
      kappa.pred0 <- predict(kappa_fit0, newdata = d_stre_s[validation_index,c("X1","X2")])$pred
      eta_fit1 <- SuperLearner(d_stre_s$Y[training_index1_r1], X = d_stre_s[training_index1_r1,c("X1","X2")], family = "gaussian", SL.library = SLmethods)
      eta_fit0 <- SuperLearner(d_stre_s$Y[training_index0_r1], X = d_stre_s[training_index0_r1,c("X1","X2")], family = "gaussian", SL.library = SLmethods)
      eta.pred1 <- predict(eta_fit1, newdata = d_stre_s[validation_index,c("X1","X2")])$pred
      eta.pred0 <- predict(eta_fit0, newdata = d_stre_s[validation_index,c("X1","X2")])$pred
      f1 <- d_stre_s$A[validation_index]/pi * d_stre_s$R[validation_index]/kappa.pred1 * (d_stre_s$RY[validation_index] - eta.pred1) + eta.pred1
      f0 <- (1-d_stre_s$A[validation_index])/(1-pi) *  d_stre_s$R[validation_index]/kappa.pred0 * (d_stre_s$RY[validation_index] - eta.pred0) + eta.pred0
      mutate(d_stre_s[validation_index,], pred1=f1, pred0=f0, fold=k)
    }) %>% Reduce(rbind,.)
  }) %>% Reduce(rbind,.)
  est_dml_stre <- mean(temp$pred1 - temp$pred0)
  V_dml_stre <- map_dbl(1:K, function(k){
    map_dbl(unique(d_stre$S), function(s){
      var(temp$pred1[temp$fold==k & temp$S==s]-temp$pred0[temp$fold==k & temp$S==s]) * mean(d_stre$S==s)
    }) %>% sum
  }) %>% mean
  VI <- n/sum(A_stre)/sum(1-A_stre) * map(unique(S), function(s){var(Xr[S==s,]) * mean(S==s)}) %>% Reduce(`+`,.)
  IF <- as.vector(temp$pred1-temp$pred0-mean(temp$pred1-temp$pred0))
  cIIF <- colMeans(diag((temp$A-pi)/pi/(1-pi) * IF) %*% as.matrix(temp[,c("X1","X2")] - rep(1, n) %*% t(colMeans(temp[,c("X1","X2")]))))
  R2_dml_stre <- t(cIIF) %*% solve(n * VI) %*% cIIF /V_dml_stre

 
  # result - 18 metrics per simulation
  c(est_dr_s, V_dr_s, 0, est_dml_s, V_dml_s, 0, 
    est_dr_re, V_dr_re, R2_dr_re, est_dml_re, V_dml_re, R2_dml_re, 
    est_dr_stre, V_dr_stre, R2_dr_stre, est_dml_stre, V_dml_stre, R2_dml_stre)
  }, error = function(e) { return(rep(NA, 18))})
}

tictoc::toc()
stopCluster(cl)

saveRDS(results, file = "simulation/sim1-2.rds")

# summarizing results

results <- readRDS(file = "simulation/sim1-2.rds")
n_sim <- ncol(results)
tt <- 1
n <- 400

est_index <- seq(1, 18, by = 3)
ase_index <- est_index + 1
r2_index <- est_index + 2
delta <- 2
z <- rnorm(10000)
d1 <- rnorm(10000)
d2 <- rnorm(10000)
index <- which(d1^2 + d2^2 < tt)
z <- z[index]
d1 <- d1[index]

results[is.infinite(results)] <- NA
for(iter in 1:n_sim){
  if(max(abs(results[,iter]), na.rm = T) > 1e6){results[,iter] <- NA}
}

cp_true <- map_dbl(1:6, function(j){
  map_dbl(1: n_sim, function(i){
    emp <- (sqrt(1-results[r2_index[j], i]) * z + sqrt(results[r2_index[j], i]) * d1) * sqrt(results[ase_index[j],i]/n) + delta
    (results[est_index[j], i] >= quantile(emp, 0.025, na.rm=T)) & (results[est_index[j], i] <= quantile(emp, 0.975, na.rm=T))
  }) %>% mean(na.rm = T)
})



summary_results <- cbind(bias = apply(results[est_index,], 1, mean, na.rm = T) - delta,
                         ese = apply(results[est_index,], 1, sd, na.rm = T), 
                         ase = apply(results[ase_index,], 1, function(x){mean(sqrt(x/n),na.rm=T)}),
                         cp_normal = apply(abs(results[est_index,] - delta)/sqrt(results[ase_index,]/n) <= qnorm(0.975), 1, mean, na.rm = T),
                         cp_true = cp_true
)

rownames(summary_results) <- rep(c("DR-WLS", "DML"), 3)

round(summary_results, 2)
xtable::xtable(summary_results)

