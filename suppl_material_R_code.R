###############################################################################################################################################
##                                                                                                    
##    Supplementary R code for                                                              
##    "Handling time-dependent exposures and confounders when estimating attributable fractions       
##     -- bridging the gap between multistate and counterfactual modeling"                          
##                                                                                                    
##    Author: Johan Steen                                                                             
##    Date: March 30, 2022           
##
##    ------
##
##    The datasets.RData file contains three data frames
##
##     events: wide format data frame which includes the following event-related variables (as defined in section 2 of the main text)
##         - id: unique patient identifier
##         - T: time from ICU admission to either ICU death or discharge, whichever occurs first (expressed in days)
##         - epsilon: event type indicator (1 = ICU death, 2 = ICU discharge)
##         - C: time from ICU admission to infection onset (expressed in days)
##         - delta: infection status at ICU death or discharge (1 = infected, 0 = uninfected)
##
##     baseline: wide format data frame which includes the following baseline patient characteristics
##         - id: unique patient identifier
##         - sex: F = female, M = male
##         - age: age (in years) at ICU admission
##         - admissionyr: year of ICU admission (coded as factor)
##         - admissioncat: admission category ("Medicine", "Emergency surgery" or "Scheduled surgery")
##         - apache2: APACHE II score at ICU admission
##         - SOFA_total_0: total SOFA score at ICU admission
##         - CCI: Charlson Comorbidity Index
##
##     timedep: long format data frame which includes the following time-dependent patient information
##         - id: unique patient identifier
##         - tstop: end time of the corresponding 24-h time interval (with follow-up starting at the time of ICU admission)
##         - SOFA_total: daily total SOFA score
##                                                                                                   
###############################################################################################################################################

rm(list=ls())

### load required packages ----
library(dplyr) # for convenient data wrangling tools
library(tidyr) # for LOCF functionality
library(magrittr) # for piping
library(survival) # for survival analysis functions
library(ipw) # for inverse probability of censoring weighting
library(splines) # for fitting splines

# optional packages (for plots)
library(prodlim)
library(ggplot2)
library(ggfortify)


### load data ----
load("/.../datasets.RData")


### transform the event time data to counting process format ----

## calculate tilde_T and tilde_epsilon as defined in the main text
events %<>% mutate(tilde_T = pmin(T,C), 
                   tilde_epsilon = delta*epsilon)

## calculate discrete time variables (under both competing risk models 1 and 2 in the main text) 
## to enable transformation to counting process format later on
## and recode event type indicators as labelled factors for convenience and clarity
events %<>% mutate(time1 = ceiling(T),
                   time2 = ceiling(tilde_T),
                   timeC = ifelse(delta == 1, ceiling(tilde_T), floor(tilde_T)),
                   event1 = factor(epsilon, levels = 0:2, labels = c("censor", "death", "discharge")),
                   event2 = factor(tilde_epsilon, levels = -1:2, labels = c("censor", "HAI onset", "HAI-free death", "HAI-free discharge")),
                   eventC = factor(tilde_epsilon, levels = 0:2, labels = c("censor", "death", "discharge")),
                   A_tau = 1-delta)
#   Note 1: these new time variables encode the end of the considered 24-h time intervals (with follow-up starting at the time of ICU admission)
# except for timeC, which can be considered a copy of time1 and time2, but with (artificial) censoring due to HAI onset at the start rather than the end of the 24-h time interval
# to enforce the temporal ordering assumption that censoring precedes ICU death or discharge (as described in section 2 of the main text)    
#   Note 2: it is important that the first factor level in the event type variables always corresponds to 'cernsoring'!
# Even if no censoring occurs, a factor level for censoring should be added and explicitly coded as first factor level.
# If not, the Surv() function may not provide the intended result when dealing with competing events.

# check
# with(events, table(event1, epsilon, useNA = "ifany"))
# with(events, table(event2, tilde_epsilon, useNA = "ifany"))
# with(events, table(eventC, tilde_epsilon, useNA = "ifany"))

## use the survival::survSplit function to transform to counting process format, setting the end of follow-up intervals (day level) as pre-specified cut times  
CPdata1 <- survSplit(Surv(time1, event1) ~ ., data = events %>% select(id, time1, event1, delta, A_tau), cut = 0:60, end = "time1") %>% rename(tstop = time1, event = event1)
CPdata2 <- survSplit(Surv(time2, event2) ~ ., data = events %>% select(id, time2, event2), cut = 0:60, end = "time2") %>% rename(tstop = time2, event = event2)
CPdata <- full_join(CPdata1 %>% rename(event1 = event), CPdata2 %>% rename(event2 = event) %>% select(-tstart), by = c("id", "tstop"))
rm(CPdata1, CPdata2)
#   Note: to enable calculation of all the discrete time counting process variables defined in section 2 of the main text, 
# we create a counting process format dataframe for the event times and types of both competing risk models 1 and 2,
# which are then joined

## calculate the discrete time counting process variables defined in section 2 of the main text
CPdata %<>% mutate(k = tstop,
                   A_k = as.numeric(event2 == "HAI onset"),
                   D_k = as.numeric(event1 == "discharge"),
                   Y_k = as.numeric(event1 == "death"))
CPdata %<>% group_by(id) %>% fill(A_k) 
#   Note: ensure that A_k is coded 1 at each follow-up interval following the interval at which HAI onset occurs (LOCF for missing A_k)
CPdata %<>% mutate(tilde_D_k = (1-A_k)*D_k,
                   tilde_Y_k = (1-A_k)*Y_k)
# calculate lag variables
CPdata %<>% group_by(id) %>% mutate(A_k_lag = lag(A_k, default = 0),
                                    D_k_lag = lag(D_k, default = 0),
                                    Y_k_lag = lag(Y_k, default = 0))


### merge dataframes ----
CPdata %<>% full_join(baseline, by = "id")
CPdata %<>% left_join(timedep, by = c("id", "tstop"))

## create lagged SOFA_total variable
CPdata %<>% group_by(id) %>% mutate(SOFA_total_lag2 = lag(SOFA_total, 2, default = 0))


### estimate observable/factual cumulative incidences ----

# calculate total sample size
n <- length(unique(CPdata$id))

## cumulative incidence of ICU death (competing risk model 1)

# Aalen-Johansen estimator via survfit function applied on original (wide) data format
mod1 <- survfit(Surv(time1, event1) ~ 1, data = events)
plot(mod1["death"], xlim = c(0, 30), fun = "event", xlab = "Days from ICU admission", ylab = "Cumulative risk of ICU death")
#ggplot2::autoplot.survfit(mod1["death"], xlim = c(0, 30), xlab = "Days from ICU admission", ylab = "Cumulative risk of ICU death")
#prodlim::prodlim(prodlim::Hist(time1, event1) ~ 1, data = events) %>% plot(cause = "death", xlim = c(0, 30), ylim = c(0, 0.2))

# Aalen-Johansen estimator via survfit function applied on counting process data format
mod1CP <- survfit(Surv(tstart, tstop, event1) ~ 1, data = CPdata, id = id)
lines(mod1CP["death"], fun = "event", col = "red")

# store cumulative risk estimates based on AJ estimator (eq2 in main text, eq15 in Appendix A.1) 
E_Y_AJ <- with(mod1CP["death"], data.frame(k = time, E_Y_AJ = pstate))

# compare with empirical cumulative distribution function (eq16 in Appendix A.1)
E_Y_ecdf <- CPdata %>% group_by(k) %>% summarise(tmp = sum(Y_k*(1-D_k)*(1-Y_k_lag))/n) %>% 
  mutate(E_Y_ecdf = cumsum(tmp)) %>% select(k, E_Y_ecdf)
E_Y <- full_join(E_Y_AJ, E_Y_ecdf, by = "k") %>% arrange(k)
E_Y


## cumulative incidence of HAI-free ICU death (competing risk model 2)

# Aalen-Johansen estimator via survfit function applied on original (wide) data format
mod2 <- survfit(Surv(time2, event2) ~ 1, data = events)
plot(mod2["HAI-free death"], fun = "event", xlim = c(0, 30), xlab = "Days from ICU admission", ylab = "Cumulative risk of HAI-free ICU death")

# Aalen-Johansen estimator via survfit function applied on counting process data format
mod2CP <- survfit(Surv(tstart, tstop, event2) ~ 1, data = CPdata, id = id)
lines(mod2CP["HAI-free death"], fun = "event", col = "red")

# store cumulative risk estimates based on AJ estimator (eq3 in main text, eq17 in Appendix A.2)
E_tilde_Y_AJ <- with(mod2CP["HAI-free death"], data.frame(k = time, E_tilde_Y_AJ = pstate))

# compare with empirical cumulative distribution function (eq4 in main text, eq18 in Appendix A.2)
E_tilde_Y_ecdf <- CPdata %>% group_by(k) %>% summarise(tmp = sum(Y_k*(1-D_k)*(1-A_k)*(1-Y_k_lag))/n) %>% 
  mutate(E_tilde_Y_ecdf = cumsum(tmp)) %>% select(k, E_tilde_Y_ecdf)
E_tilde_Y <- full_join(E_tilde_Y_AJ, E_tilde_Y_ecdf, by = "k") %>% arrange(k)
E_tilde_Y



### treating exposure as an exclusion criterion (naive approach to estimate hypothetical/counterfactual cumulative incidence) ----

## Aalen-Johansen estimator via survfit function applied on original (wide) data format
mod_naive1 <- survfit(Surv(time1, event1) ~ 1, data = events %>% filter(A_tau == 0))
mod_naive2 <- survfit(Surv(time2, event2) ~ 1, data = events %>% filter(A_tau == 0))
mod_naive3 <- survfit(Surv(time1, event1) ~ A_tau, data = events)
mod_naive4 <- survfit(Surv(time2, event2) ~ A_tau, data = events)
mod_naive5 <- survfit(Surv(time1, event1) ~ 1, data = events, weights = 1-A_tau)
mod_naive6 <- survfit(Surv(time2, event2) ~ 1, data = events, weights = 1-A_tau)
plot(mod_naive1["death"], fun = "event", xlim = c(0, 30), main = "Treating HAI as\n an exclusion criterion", xlab = "Days from ICU admission", ylab = "Counterfactual cumulative risk of ICU death")
lines(mod_naive2["HAI-free death"], fun = "event", col = "red")
lines(mod_naive3["A_tau=0", "death"], fun = "event", col = "green")
lines(mod_naive4["A_tau=0", "HAI-free death"], fun = "event", col = "purple")
lines(mod_naive5["death"], fun = "event", col = "red")
lines(mod_naive6["HAI-free death"], fun = "event", col = "green")
#   Note that all the above survfit commands provide the same result because
# 1/ when restricting the data to patients who don't develop an infection during hospitalization,
# for each k, we have Y_k = tilde_Y_k;
# 2/ restricting the data to patients who don't develop an infection is equivalent to stratification 
# and to setting weights of infected patients to zero.

## Aalen-Johansen estimator via survfit function applied on counting process data format
mod_naiveCP1 <- survfit(Surv(tstart, tstop, event1) ~ 1, data = CPdata %>% filter(A_tau == 0), id = id)
mod_naiveCP2 <- survfit(Surv(tstart, tstop, event1) ~ A_tau, data = CPdata, id = id)
mod_naiveCP3 <- survfit(Surv(tstart, tstop, event1) ~ 1, data = CPdata, id = id, weights = 1-A_tau) 
lines(mod_naiveCP1["death"], fun = "event", col = "purple")
lines(mod_naiveCP2["A_tau=0", "death"], fun = "event", col = "red")
lines(mod_naiveCP3["death"], fun = "event", col = "green")

## store cumulative risk estimates based on AJ estimator (eq24 in Appendix B.2)
E_Y0_naive_AJ <- with(mod_naiveCP1["death"], data.frame(k = time, E_Y0_naive_AJ = pstate))

## compare with weighted empirical cumulative distribution function (eq10 in main text, eq25 in Appendix B.2)
#fit_naive <- glm(A_tau ~ 1, data = CPdata %>% group_by(id) %>% slice(n()))
#CPdata$W_naive <- 1/(1-predict(fit_naive, type = "response", newdata = CPdata))
class(CPdata) <- "data.frame"
#   Note: recode CPdata as a data.frame object to avoid error in ipwpoint function below
ipw_naive <- ipwpoint(exposure = A_tau, family = "binomial", link = "logit", denominator = ~ 1, data = events)
events$W_naive <- ipw_naive$ipw.weights
CPdata %<>% left_join(events %>% select(id, W_naive))
#   Note: first calculate naive weights (see lines above). 
# As this corresponds to (inappropriately) considering HAI onset as a point exposure at baseline, 
# we can simply use the ipwpoint function from the ipw package for this purpose, specifying an empty covariate set.
# However, in order to obtain the correct weights, this function needs to applied to the original (wide format) dataset, 
# and the resulting weights then need to be added to the counting process format dataset.
E_Y0_naive_ecdf <- CPdata %>% group_by(k) %>% summarise(tmp = sum(Y_k*(1-D_k)*(1-A_k)*(1-Y_k_lag)*W_naive)/n) %>% 
  mutate(E_Y0_naive_ecdf = cumsum(tmp)) %>% select(k, E_Y0_naive_ecdf)
E_Y0_naive <- full_join(E_Y0_naive_AJ, E_Y0_naive_ecdf, by = "k") %>% arrange(k)
E_Y0_naive

## estimate population attributable fraction
PAF_naive <- full_join(E_Y_AJ, E_Y0_naive_AJ, by = "k") %>% arrange(k)
PAF_naive %<>% mutate(PAF_naive = (E_Y_AJ - E_Y0_naive_AJ)/E_Y_AJ)
PAF_naive


### treating exposure onset as a time-dependent exclusion criterion (Schumacher et al 2007 approach to estimate hypothetical/counterfactual cumulative incidence) ----

#   Note: because the time-dependent exclusion criterion (time-dependent weights) is (are) indexed by subsequent landmark intervals K at which the cumulative incidence is evaluated
# rather than at follow-up intervals k over which the (weighted) proportion of deaths are summed, 
# this time-dependent exclusion criterion cannot be directly applied in the survfit function, neither by exclusion, stratification nor weighting (as before)
# Instead, we apply the original estimator as proposed by Schumacher et al 2007 as a functional of the estimated cumulative incidence of ICU death and HAI onset, respectively,
# as detailed in Appendix B3.

Pprime01 <- E_tilde_Y_AJ %>% rename(Pprime01 = E_tilde_Y_AJ)
Pprime03 <- with(mod2CP["HAI onset"], data.frame(k = time, Pprime03 = pstate))

## store cumulative risk estimates based on AJ estimators (see Appendix B.2)
E_Y0_schu_AJ <- full_join(Pprime01, Pprime03, by = "k")
E_Y0_schu_AJ %<>% mutate(E_Y0_schu_AJ = Pprime01/(1-Pprime03)) %>% select(k, E_Y0_schu_AJ)

plot(E_Y0_schu_AJ ~ k, data = E_Y0_schu_AJ, type = "s", xlim = c(0, 30), main = "Treating HAI onset as\n a time-dependent exclusion criterion", xlab = "Days from ICU admission", ylab = "Counterfactual cumulative risk of ICU death")

## compare with weighted empirical cumulative distribution function (eq11 in main text, eq28 in Appendix B.3)
#   Note: calculating the weights defined in section 4.3 in the main text (and in Appendix B.3)
# requires to extend the counting process format data to the maximal follow-up interval tau (i.e. beyond the observed event time) for each patient 
tau <- 60
CPdata %<>% filter(k <= tau)
CPdata_ext <- full_join(expand.grid(k = 1:tau, id = unique(CPdata$id)), CPdata, by = c("id", "k"))
CPdata_ext %<>% fill(A_k, D_k, Y_k, tilde_D_k, tilde_Y_k)
# recalculate lag variables
CPdata_ext %<>% group_by(id) %>% mutate(A_k_lag = lag(A_k, default = 0),
                                        D_k_lag = lag(D_k, default = 0),
                                        Y_k_lag = lag(Y_k, default = 0))
#fit_schu <- glm(A_k ~ factor(k), data = CPdata_ext)
#CPdata_ext$W_schu <- 1/(1-predict(fit_schu, type = "response", newdata = CPdata_ext))
class(CPdata_ext) <- "data.frame"
#   Note: recode CPdata as a data.frame object to avoid error in ipwpoint function below
W_schu <- data.frame(k = 1:60, W_schu = NA)
for(K in 1:60) {
  ipw_schu <- ipwpoint(exposure = A_k, family = "binomial", link = "logit", denominator = ~ 1, data = CPdata_ext %>% filter(k == K))
  W_schu[K, "W_schu"] <- ipw_schu$ipw.weights[1] 
  # Note: make sure to pick the id/index of a subject that is uninfected during the entire hospitalization to extract the weight!
}
CPdata %<>% left_join(W_schu, by = "k")
CPdata_ext %<>% left_join(W_schu, by = "k")

#   Note: to calculate the weighted ecdf first sum over intervals k per subject, then multiply by weight indexed by K, and then average over subjects for each k 
# (rather than first averaging over subjects at each k and then taking cumulative sum as before)
E_Y0_schu_ecdf <- CPdata_ext %>% arrange(id, k) %>% group_by(id) %>% mutate(tmp = cumsum(Y_k*(1-D_k)*(1-A_k)*(1-Y_k_lag))*W_schu) %>% 
  group_by(k) %>% summarise(E_Y0_schu_ecdf = sum(tmp)/n) %>% select(k, E_Y0_schu_ecdf)
E_Y0_schu <- full_join(E_Y0_schu_AJ, E_Y0_schu_ecdf, by = "k") %>% arrange(k)
E_Y0_schu

## estimate population attributable fraction
PAF_schu <- full_join(E_Y_AJ, E_Y0_schu_AJ, by = "k") %>% arrange(k)
PAF_schu %<>% mutate(PAF_schu = (E_Y_AJ - E_Y0_schu_AJ)/E_Y_AJ)
PAF_schu


### treating exposure onset as a (non-informative) censoring event (Keiding et al 2001 approach to estimate hypothetical/counterfactual cumulative incidence) ----

## Aalen-Johansen estimator via survfit function applied on original (wide) data format
mod_keiding <- survfit(Surv(timeC, eventC) ~ 1, data = events)
plot(mod_keiding["death"], fun = "event", xlim = c(0, 30), main = "Treating HAI onset as\n non-informative censoring", xlab = "Days from ICU admission", ylab = "Counterfactual cumulative risk of ICU death")

## Aalen-Johansen estimator via survfit function applied on counting process data format
mod_keidingCP1 <- survfit(Surv(tstart, tstop, event1) ~ 1, data = CPdata %>% filter(A_k == 0), id = id)
mod_keidingCP2 <- survfit(Surv(tstart, tstop, event1) ~ A_k, data = CPdata, id = 1:nrow(CPdata))
mod_keidingCP3 <- survfit(Surv(tstart, tstop, event1) ~ 1, data = CPdata, id = 1:nrow(CPdata), weights = 1-A_k) 
lines(mod_keidingCP1["death"], fun = "event", col = "red")
lines(mod_keidingCP2["A_k=0", "death"], fun = "event", col = "green")
lines(mod_keidingCP3["death"], fun = "event", col = "purple")
#   Note that, to guarantee identical results, when stratifying on the time-dependent exposure status A_k or applying time-dependent weights 1-A_k,
# the id argument needs to be specified differently, treating each row in the dataset as a 'pseudo-patient'.

## store cumulative risk estimates based on AJ estimator (eq6 in the main text, eq22 in Appendix B.1)
E_Y0_keiding_AJ <- with(mod_keidingCP1["death"], data.frame(k = time, E_Y0_keiding_AJ = pstate))

## compare with weighted empirical cumulative distribution function (eq7 in main text, eq23 in Appendix B.1)
#fit_naive <- glm(A_tau ~ 1, data = CPdata %>% group_by(id) %>% slice(n()))
#CPdata$W_naive <- 1/(1-predict(fit_naive, type = "response", newdata = CPdata))
class(CPdata) <- "data.frame"
#   Note: recode CPdata as a data.frame object to avoid error in ipwtm function below
ipw_keiding <- ipwtm(exposure = A_k, family = "survival",
                     denominator = ~ 1, id = id,
                     tstart = tstart, timevar = tstop,
                     type = "first", data = CPdata)
CPdata$W_keiding <- ipw_keiding$ipw.weights
#   Note: first calculate IPC weights based on time-dependent propensity score models conditional on empty covariate sets (see lines above)
E_Y0_keiding_ecdf <- CPdata %>% group_by(k) %>% summarise(tmp = sum(Y_k*(1-D_k)*(1-A_k)*(1-Y_k_lag)*W_keiding)/n) %>% 
  mutate(E_Y0_keiding_ecdf = cumsum(tmp)) %>% select(k, E_Y0_keiding_ecdf)
E_Y0_keiding <- full_join(E_Y0_keiding_AJ, E_Y0_keiding_ecdf, by = "k") %>% arrange(k)
E_Y0_keiding
#   Note: although estimators are algebraically equivalent, estimates are not identical because of the way weights are calculated using the ipwtm function.
# Identical estimates can be obtained by self-calculated IPC weights (e.g. based on a pooled logit model or the KM estimator). 

## estimate population attributable fraction
PAF_keiding <- full_join(E_Y_AJ, E_Y0_keiding_AJ, by = "k") %>% arrange(k)
PAF_keiding %<>% mutate(PAF_keiding = (E_Y_AJ - E_Y0_keiding_AJ)/E_Y_AJ)
PAF_keiding


### treating exposure onset as an informative censoring event (IPCW approach to estimate hypothetical/counterfactual cumulative incidence) ----

## first calculate IPC weights that account for confounder history up to each time
class(CPdata) <- "data.frame"
#   Note: recode CPdata as a data.frame object to avoid error in ipwtm function below
ipw_ipcw <- ipwtm(exposure = A_k, family = "survival",
                  denominator = ~ sex + as.numeric(admissioncat != "Medicine") + admissionyr + ns(age, df = 4) + ns(CCI, df = 2) + ns(SOFA_total_0, df = 4) + ns(SOFA_total_lag2, df = 4), id = id,
                  tstart = tstart, timevar = tstop,
                  type = "first", data = CPdata)
CPdata$W_ipcw <- ipw_ipcw$ipw.weights

## IPC weighted Aalen-Johansen estimator cannot be obtained via survfit function applied on original (wide) data format, because this does not allow accounting for time-dependent confounding

## IPC weighted Aalen-Johansen estimator via survfit function applied on counting process data format
mod_ipcwCP1 <- survfit(Surv(tstart, tstop, event1) ~ 1, data = CPdata %>% filter(A_k == 0), id = 1:nrow(CPdata %>% filter(A_k == 0)), weights = W_ipcw)
mod_ipcwCP2 <- survfit(Surv(tstart, tstop, event1) ~ A_k, data = CPdata, id = 1:nrow(CPdata), weights = W_ipcw)
mod_ipcwCP3 <- survfit(Surv(tstart, tstop, event1) ~ 1, data = CPdata, id = 1:nrow(CPdata), weights = (1-A_k)*W_ipcw) 
plot(mod_ipcwCP1["death"], fun = "event", xlim = c(0, 30), main = "Treating HAI onset as\n informative censoring", xlab = "Days from ICU admission", ylab = "Counterfactual cumulative risk of ICU death")
lines(mod_ipcwCP2["A_k=0", "death"], fun = "event", col = "red")
lines(mod_ipcwCP3["death"], fun = "event", col = "green")
#   Note: to guarantee appropriate estimation as intended,
# the id argument needs to be specified such that each row in the dataset is treated as an independent 'pseudo-patient'.
# As a result of this (and the uncertainty in estimation of the weights), confidence intervals can no longer be expected to have nominal coverage.
# Instead, we recommend using the non-parametric bootstrap.

## store cumulative risk estimates based on AJ estimator (eq12 in the main text, eq34 in Appendix B.4)
E_Y0_ipcw_AJ <- with(mod_ipcwCP1["death"], data.frame(k = time, E_Y0_ipcw_AJ = pstate))

## compare with weighted empirical cumulative distribution function (eq14 in main text, eq35 in Appendix B.4)
E_Y0_ipcw_ecdf <- CPdata %>% group_by(k) %>% summarise(tmp = sum(Y_k*(1-D_k)*(1-A_k)*(1-Y_k_lag)*W_ipcw)/n) %>% 
  mutate(E_Y0_ipcw_ecdf = cumsum(tmp)) %>% select(k, E_Y0_ipcw_ecdf)
E_Y0_ipcw <- full_join(E_Y0_ipcw_AJ, E_Y0_ipcw_ecdf, by = "k") %>% arrange(k)
E_Y0_ipcw

## estimate population attributable fraction
PAF_ipcw <- full_join(E_Y_AJ, E_Y0_ipcw_AJ, by = "k") %>% arrange(k)
PAF_ipcw %<>% mutate(PAF_ipcw = (E_Y_AJ - E_Y0_ipcw_AJ)/E_Y_AJ)


### compare estimated PAFs ----
PAF_naive %>% select(k, PAF_naive) %>%
  full_join(PAF_schu %>% select(k, PAF_schu), by = "k") %>%
  full_join(PAF_keiding %>% select(k, PAF_keiding), by = "k") %>%
  full_join(PAF_ipcw %>% select(k, PAF_ipcw), by = "k") %>% arrange(k)
