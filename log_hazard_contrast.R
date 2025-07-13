
print("# Packages et donnees  ----------------------------------------------------")
library(dplyr)
library(tidyverse)
library(tmle3)
library(sl3)
library(kableExtra)
library(data.table)
library(ggplot2)
library(parallel)
library(future)
library(delayed)
library(msm)
library(profvis)
library(Hmisc)
library(tidyverse)
library(ranger)
library(hrbrthemes)
library(tidyr)
library(viridis)
library(glmnet)
library(doParallel)
library(foreach)
library(mgcv)
library(MatchIt)
library(leaps)

print("# Parallel ----------------------------------------------------------------")
detectCores()
availableCores()
cpus_physical = detectCores()
supportsMulticore()
Sys.setenv(R_FUTURE_FORK_ENABLE = TRUE)
supportsMulticore()
plan(multisession)

print("# Parallel ----------------------------------------------------------------")
simulate_data <- function(n_sim = 2e2) {
  library(simcausal)
  D <- DAG.empty()
  D <- D +
    node("W1", distr = "rbinom", size = 1, prob = .5) +
    node("W2", distr = "runif", min = 0, max = 1.5) +
    node("W3", distr = "rnorm", mean = 0, sd = 1) +
    node("W4", distr = "rnorm", mean = 0, sd = 1) +

    node("W", distr = "runif", min = 0, max = 1.5) +
    node("A", distr = "rbinom", size = 1, prob = plogis(.15 + .5 * as.numeric(W > .75) - 0.75*W2 + W3 - 0.3*W2*W3 + 1.5*W4 - 0.2*W4^2)) +
    node("Trexp", distr = "rexp", rate = 1 + .7 * W^2 - 0.8*A) +
    node("Cweib", distr = "rweibull", shape = 1 + .5 * W, scale = 75) +
    node("T", distr = "rconst", const = round(Trexp * 2)) +
    node("C", distr = "rconst", const = round(Cweib * 2)) +
    # Observed random variable (follow-up time):
    node("T.tilde", distr = "rconst", const = ifelse(T <= C, T, C)) +
    # Observed random variable (censoring indicator, 1 - failure event, 0 - censored):
    node("Delta", distr = "rconst", const = ifelse(T <= C, 1, 0))
  setD <- set.DAG(D)
  dat <- sim(setD, n = n_sim)
  # only grab ID, W's, A, T.tilde, Delta
  Wname <- grep("W", colnames(dat), value = TRUE)
  dat <- dat[, c("ID", Wname, "A", "T.tilde", "Delta")]
  # input: scalar q, W vector. computes for all W, the S(q|A,W)
  true_surv_one <- function(q, W, A = 1) sapply(W, function(w) {
    1 - pexp(q, rate = 1 + .7 * w^2 - .8 * A)
  })
  # input: vector q. mean(S(q|A,W)|A), average out W. loop over q
  true_surv <- function(q_grid, surv_fn, A) {
    W_grid <- seq(0, 1.5, .01)
    survout <- numeric()
    for (q in q_grid) survout <- c(survout, mean(surv_fn(q = q / 2, W = W_grid, A = A)))
    return(survout)
  }
  truth_surv <- function(q) true_surv(q_grid = q, surv_fn = true_surv_one, A = 1)
  truth_surv0 <- function(q) true_surv(q_grid = q, surv_fn = true_surv_one, A = 0)
  return(list(dat = dat, true_surv1 = truth_surv, true_surv0 = truth_surv0))
}

sim <- simulate_data(n = 500)
d <- sim$dat

# max(d$T.tilde)
# table(d$Delta)
# table(d$A)

# sim$true_surv1(1:36)
# sim$true_surv0(1:36)
# psi_test
# psi_control
# initial_psi_test
# initial_psi_control
# 
# exp(mean(log(log_psi_test/log_psi_control)))
# exp(mean(log(log_psi_test/log_psi_control)))
# 
# exp(mean(log(log(initial_psi_test)/log(initial_psi_control))))
# exp(mean(log(log(sim$true_surv1(1:36))/log(sim$true_surv0(1:36)))))



data <- d

#data <- data %>% select(-c("BMI", "ETHNI_5CL"))
str(data)
data <- data %>% mutate(Delta = factor(Delta))
data <- data %>% mutate(A = factor(A))
data <- data %>% mutate(W1 = factor(W1))

str(data)


# print("# Sampling  ----------------------------------------------------------")
# tab0 <- subset(data,TREAT==0) #sélectionne les lignes de data qui ont un "0" comme équivalent dans tab0
# tab1 <- subset(data,TREAT==1) #sélectionne les lignes de data qui ont un "1" comme équivalent dans tab1
# newtab0 <- tab0[sample(1:nrow(tab0),2*3256),] #tirage aléatoire de n0 lignes
# newtab1 <- tab1[sample(1:nrow(tab1),3256),] #tirage aléatoire de n1 lignes
# data <- rbind(newtab0,newtab1)


# print("#Censurer les patients dont la durée du suivi est > 180 jours (6 mois environ) -------------------")
# cens <- quantile(data$TIME, 0.75) 
# data <- data %>%
#   mutate(OUTCOME = replace(OUTCOME, TIME >= cens, "0"))
# 
# data <- data %>%
#   mutate(TIME = replace(TIME, TIME >= cens, cens))


print("# Create the time point grid ----------------------------------------------")
k_grid <- 1:max(data$T.tilde)


print("# TMLE : Specify right-censored survival data arguments ---------------")
vars <- colnames(data)
Covars <- vars[! vars %in% c("ID", "A", "Delta", "T.tilde")]

var_types <- list(
  T_tilde = Variable_Type$new("continuous"),
  t = Variable_Type$new("continuous"),
  Delta = Variable_Type$new("binomial"))

survival_spec <- tmle_survival(
  treatment_level = 1, control_level = 0,
  target_times =  k_grid,
  variable_types = var_types)

node_list <- list(
  W = Covars,
  A = "A", T_tilde = "T.tilde", Delta = "Delta", id = "ID")


print("# transform data ------------------------------------------------------")
long_data_tuple <- survival_spec$transform_data(data, node_list)
df_long <- long_data_tuple$long_data
long_node_list <- long_data_tuple$long_node_list


print("# TL  task  ----------------------------------------------------------")
tmle_task <- survival_spec$make_tmle_task(df_long, long_node_list)


print("# Obtain initial estimates with SuperLearner sl3 R package ------------------------------------------------")
sl3_list_learners(c("binomial"))
lrnr_glm <- make_learner(Lrnr_glm)                       
lrnr_rf <- make_learner(Lrnr_randomForest)
lrnr_gam <- make_learner(Lrnr_gam, deg.gam = 2)

print("# preparation du super learner  ------------------------------------------")
sl_N <- Lrnr_sl$new(learners = list(lrnr_glm, lrnr_rf, lrnr_gam))
sl_C <- Lrnr_sl$new(learners = list(lrnr_glm, lrnr_rf, lrnr_gam))
sl_A <- Lrnr_sl$new(learners = list(lrnr_glm, lrnr_rf, lrnr_gam))
learner_list <- list(A = sl_A, N = sl_N, A_c = sl_C)


print("# Initial likelihood fit --------------------------------------------------")
start_time <- Sys.time()
initial_likelihood <- survival_spec$make_initial_likelihood(tmle_task, learner_list)
# save(initial_likelihood, file = "initial_likelihood.RData")
end_time <- Sys.time()
end_time - start_time


print("# Initial fit --------------------------------------------------")
initial_fit <- initial_likelihood$get_likelihoods(tmle_task)
head(initial_fit)


print(" PS --------------------------------------------------")
ps <- initial_fit$A
summary(ps)

data_ps <- cbind(data, ps = ps)
p6 <- ggplot(data = data_ps, aes(x = ps, group = A, fill = A)) +
  geom_density(adjust = 1.5, alpha = 0.4) +
  theme_ipsum() +
  xlab("pscore")+
  ggtitle("")


# Learner and ps plot  -----------------------------------------------------
N_learner <- initial_likelihood$factor_list[["N"]]
N_learner$learner


Ac_learner <- initial_likelihood$factor_list[["A_c"]]
Ac_learner$learner


ps_learner <- initial_likelihood$factor_list[["A"]]
ps_learner$learner

# update -----------------------------------------

up <- tmle3_Update_survival$new(
  maxit = 5000,
  cvtmle = FALSE,
  convergence_type = "sample_size",
  delta_epsilon = 1e-2,
  fit_method = "l2",
  constrain_step = TRUE,
  verbose = TRUE
)

print("# Perform TMLE adjustment of the initial conditional survival estimate --------")

targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood,
                                               updater = up)

tmle_params <- survival_spec$make_params(tmle_task,
                                         targeted_likelihood)

tmle_fit_manual <- fit_tmle3(tmle_task,
                             targeted_likelihood, tmle_params,
                             targeted_likelihood$updater)

initial_psi_test <- tmle_fit_manual$initial_psi
estimates_test <- tmle_fit_manual$estimates
psi_test <- estimates_test[[1]]$psi
log_psi_test <- log(psi_test)


print("# New spec, control arm ----------------------------------------------------------------")
vars <- colnames(data)
Covars <- vars[! vars %in% c("ID", "A", "Delta", "T.tilde")]
var_types <- list(
  T_tilde = Variable_Type$new("continuous"),
  t = Variable_Type$new("continuous"),
  Delta = Variable_Type$new("binomial"))

survival_spec <- tmle_survival(
  treatment_level = 0, control_level = 1,
  target_times = k_grid,
  variable_types = var_types
)

node_list <- list(
  W = Covars,
  A = "A", T_tilde = "T.tilde", Delta = "Delta", id = "ID"
)

print("# transform data ------------------------------------------------------")
long_data_tuple <- survival_spec$transform_data(data, node_list)
df_long <- long_data_tuple$long_data
long_node_list <- long_data_tuple$long_node_list


print("# TL  task  ----------------------------------------------------------")
tmle_task <- survival_spec$make_tmle_task(df_long, long_node_list)


print("# Perform TMLE adjustment of the initial conditional survival esti --------")
up <- tmle3_Update_survival$new(
  maxit = 5000,
  cvtmle = FALSE,
  convergence_type = "sample_size",
  delta_epsilon = 1e-2,
  fit_method = "l2",
  constrain_step = TRUE,
  verbose = TRUE
)

targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater = up)
tmle_params <- survival_spec$make_params(tmle_task, targeted_likelihood)

tmle_fit_manual <- fit_tmle3(tmle_task,
                             targeted_likelihood, tmle_params,
                             targeted_likelihood$updater)

initial_psi_control <- tmle_fit_manual$initial_psi
estimates_control <- tmle_fit_manual$estimates
psi_control <- estimates_control[[1]]$psi
log_psi_control <- log(psi_control)


print("# HRm -----------------------------------------------------------")
res <- data.frame(True = exp(mean(log(log(sim$true_surv1(1:max(data$T.tilde)))/log(sim$true_surv0(1:max(data$T.tilde)))))),  HRm = exp(mean(log(log_psi_test/log_psi_control))), HRm_init = exp(mean(log(log(initial_psi_test)/log(initial_psi_control)))))


print(" Return res ------------------------------------------")
write.csv(res, file =  "res.csv")

print(" End ------------------------------------------")
