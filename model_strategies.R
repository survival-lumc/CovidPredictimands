#------------------------------------------------------------------------------#
# Load packages and data
#------------------------------------------------------------------------------#
library(tidyverse)
library(survival)
library(mstate) # requires version 0.3.2
library(ggpubr)

load("Data/datahosp.rda")

dfr <- data_hosp

# Variables description:
# * Age_cat: Patient age, categorized (less or equal to 50, 50 - 59, 60 - 69, 70 - 79, 80 - 89, over 90)
# * Sex_m: Patient gender (0 = Male, 1 = Female)
# * N_comobidities: Number of medical conditions, capped at 3
# * DaysToDeath: Days from hospitalization and positive test to death, censored at 28 days
# * Death: Death status indicator related to DaysToDeath (0 = alive, 1 = dead)
# * DaysToICU: Days from hospitalization and positive test to death, censored at 28 days
# * Death: ICU treatment status indicator related to DaysToICU (0 = untreated, 1 = treated)

#------------------------------------------------------------------------------#
# Analysis 
#------------------------------------------------------------------------------#
### 1) Ignore Treatment --------------------------------------------------------
### ----------------------------------------------------------------------------
# Fit Cox proportional hazards model using death as event
IT_hosp <- coxph(Surv(DaysToDeath, Death) ~ Age_cat + Sex_m + N_Comorbidities + WaveII, 
                 data = dfr)

### 2) Composite Strategy ------------------------------------------------------
### ----------------------------------------------------------------------------
# Fit Cox proportional hazards model using the composite outcome as event
dfr$ICUOrDeath <- if_else(dfr$ICU == 1 | dfr$Death == 1, 1, 0)
dfr$DaysToICUOrDeath <- pmin(dfr$DaysToICU, dfr$DaysToDeath, na.rm = TRUE)
CS_hosp <- coxph(Surv(DaysToICUOrDeath, ICUOrDeath) ~ Age_cat + Sex_m + N_Comorbidities + WaveII, 
                 data=dfr)

### 3) While untreated ---------------------------------------------------------
### ----------------------------------------------------------------------------
# Fit two cause-specific hazard models using the mstate package.
# Explanation of the code here below can be found on the tutorial on competing risks 
# by Putter, Fiocco and Geskus at https://onlinelibrary.wiley.com/doi/10.1002/sim.2712
tmat <- trans.comprisk(2, c("Death","ICU")) # Define transition matrix 
dfr$MstateICUDeath <-  case_when(dfr$ICU == 1 ~ 2, dfr$Death == 1 ~ 1, TRUE ~ 0) 
dfr$stat1 <- as.numeric(dfr$MstateICUDeath == 1) # Define status variable for untreated death
dfr$stat2 <- as.numeric(dfr$MstateICUDeath == 2) # Define status variable for ICU admittance

# Prepare data for multistate model
data_mstate <- msprep(time = c(NA, "DaysToICUOrDeath", "DaysToICUOrDeath"), 
                      status = c(NA, "stat1", "stat2"),
                      data = dfr,
                      keep = c("Age_cat", "Sex_m", "N_Comorbidities", "WaveII"),
                      trans = tmat)
data_mstate <- expand.covs(data_mstate, c("Age_cat", "Sex_m", "N_Comorbidities", "WaveII"))

# Fit model
pred_varnames <- c(names(data_mstate)[c(grep(".1",names(data_mstate)), grep(".2",names(data_mstate)))])
fmla <- formula(paste("Surv(time, status) ~ ", paste(pred_varnames, collapse = "+"),
                      " + strata(trans)"))
WU_hosp <- coxph(fmla, data = data_mstate, method = "breslow")

### 4) Hypothetical ------------------------------------------------------------
### ----------------------------------------------------------------------------
# Turn dfr into a long dataframe where we the variable ICU is now a treatment indicator
# that is 0 before a patient receives treatment and 1 after treatment (if a patient
# is never treated, it will only take the value 0). In this new long dataset
# patients that are never admitted to the ICU will have only have line,
# while patients who are, at a certain point, admitted to the ICU will have two lines.
dfr_long <- dfr %>%
  mutate(id = 1:nrow(.)) %>%
  pivot_longer(
    cols = c(DaysToICU, DaysToDeath),
    names_to = "DaysTo",
    values_to = "tstop"
  ) %>%
  filter(tstop != 0) %>%
  mutate(
    ICU = if_else(DaysTo == "DaysToICU", 0, ICU),
    Death = if_else(DaysTo == "DaysToICU", 0, Death),
    tstart = 0
  ) %>%
  select(id, tstart, tstop, Age_cat, Sex_m, N_Comorbidities, WaveII, ICU, Death) 
dfr_long$tstart[diff(c(0, dfr_long$id)) == 0] <- dfr_long$tstop[diff(c(dfr_long$id, 0)) == 0]

# Fit model
Hyp_hosp <- coxph(Surv(tstart, tstop, Death) ~ ICU * (Age_cat + Sex_m + N_Comorbidities + WaveII), 
                  data=dfr_long) 


### 5) Pictures ----------------------------------------------------------------
### ----------------------------------------------------------------------------

### Figure 1: scatterplots -----------------------------------------------------
# Choose covariate values for the predictions
scatter_tpred <- 28
scatter_age <- levels(dfr$Age_cat)
scatter_sex <- c(0,1)
scatter_comorb <- 0:3
scatter_therapy <- 0
scatter_wave <- 1

# Create dataset for prediction with the chosen values
scatt_hosp <- expand.grid(Age_cat = scatter_age,
                          Sex_m = scatter_sex, 
                          N_Comorbidities = scatter_comorb,
                          ICU = scatter_therapy,
                          WaveII = scatter_wave,
                          DaysToDeath = scatter_tpred,
                          Death = 0,
                          DaysToICUOrDeath = scatter_tpred,
                          tstart = 0,
                          tstop = scatter_tpred,
                          ICUOrDeath = 0,
                          stat1 = 0,
                          stat2 = 0)

scatt_IT_hosp <- 1 - predict(IT_hosp, newdata = scatt_hosp, type = "survival")   # predict IT strategy 
scatt_CS_hosp <- 1 - predict(CS_hosp, newdata = scatt_hosp, type = "survival")   # predict CS strategy 
scatt_Hyp_hosp <- 1 - predict(Hyp_hosp, newdata = scatt_hosp, type = "survival") # predict Hyp strategy

# Expand the prediction dataset for mstate prediction
data_mstate_scatt <- msprep(time = c(NA, "DaysToICUOrDeath", "DaysToICUOrDeath"),
                            status = c(NA, "stat1", "stat2"),
                            data = scatt_hosp,
                            keep = c("Age_cat", "Sex_m", "N_Comorbidities", "WaveII"),
                            trans = tmat)
data_mstate_scatt <- expand.covs(data_mstate_scatt, c("Age_cat", "Sex_m", "N_Comorbidities", "WaveII"))
data_mstate_scatt$strata <- data_mstate_scatt$trans

scatt_WU_hosp <- rep(NA_real_, nrow(data_mstate_scatt)/2)
for (i in 1:(nrow(data_mstate_scatt)/2)) { 
  msf_scatt <- msfit(WU_hosp, data_mstate_scatt[c(2*i-1, 2*i),], trans = tmat) # predict from multistate model
  pti <- probtrans(msf_scatt, 0)[[1]] # WU prediction at all times
  scatt_WU_hosp[i] <- pti$pstate2[pti$time == scatter_tpred] # WU prediction for patient i at the wanted time
} 

# Prepare data frame for ggplot figure
scatt_hosp <- scatt_hosp %>%
  mutate(IT = scatt_IT_hosp, CS = scatt_CS_hosp, WU = scatt_WU_hosp, 
         Hyp = scatt_Hyp_hosp, Sex = if_else(Sex_m == 1, "Male", "Female"))
scatt_hosp$N_Comorbidities <- factor(scatt_hosp$N_Comorbidities, levels = c(0,1,2,3),
                                     labels = c("0", "1", "2", "\u2265 3"))

# The 6 scatterplots of Figure 1:
p1 <- ggplot(scatt_hosp, aes(IT, WU, color = Age_cat, shape = Sex, size = N_Comorbidities)) +
  geom_point(alpha = 0.7)+
  scale_size_discrete(range = c(1,2)) +
  scale_color_viridis_d(direction =  -1) +
  labs(x = "Ignore Treatment", y = "While Untreated", color = "Age",         
       size = "Medical conditions", shape = "Sex") +
  geom_abline(intercept = 0, slope = 1, color = "lightgrey") +
  xlim(c(0, 0.8)) +
  ylim(c(0, 0.8)) +
  theme_bw()
p2 <- ggplot(scatt_hosp, aes(IT, CS, color = Age_cat, shape = Sex, size = N_Comorbidities)) +
  geom_point(alpha = 0.7)+
  scale_size_discrete(range = c(1,2)) +
  scale_color_viridis_d(direction =  -1) +
  labs(x = "Ignore Treatment", y = "Composite Strategy", color = "Age",         
       size = "Medical conditions", shape = "Sex") +
  geom_abline(intercept = 0, slope = 1, color = "lightgrey") +
  xlim(c(0, 0.8)) +
  ylim(c(0, 0.8)) +
  theme_bw()
p3 <- ggplot(scatt_hosp, aes(IT, Hyp, color = Age_cat, shape = Sex, size = N_Comorbidities)) +
  geom_point(alpha = 0.7)+
  scale_size_discrete(range = c(1,2)) +
  scale_color_viridis_d(direction =  -1) +
  labs(x = "Ignore Treatment", y = "Hypothetical", color = "Age",         
       size = "Medical conditions", shape = "Sex") +
  geom_abline(intercept = 0, slope = 1, color = "lightgrey") +
  xlim(c(0, 0.8)) +
  ylim(c(0, 0.8)) +
  theme_bw()
p4 <- ggplot(scatt_hosp, aes(WU, CS, color = Age_cat, shape = Sex, size = N_Comorbidities)) +
  geom_point(alpha = 0.7)+
  scale_size_discrete(range = c(1,2)) +
  scale_color_viridis_d(direction =  -1) +
  labs(x = "While Untreated", y = "Composite Strategy", color = "Age",         
       size = "Medical conditions", shape = "Sex") +
  geom_abline(intercept = 0, slope = 1, color = "lightgrey") +
  xlim(c(0, 0.8)) +
  ylim(c(0, 0.8)) +
  theme_bw()
p5 <- ggplot(scatt_hosp, aes(Hyp, CS, color = Age_cat, shape = Sex, size = N_Comorbidities)) +
  geom_point(alpha = 0.7)+
  scale_size_discrete(range = c(1,2)) +
  scale_color_viridis_d(direction =  -1) +
  labs(x = "Hypothetical", y = "Composite Strategy", color = "Age",         
       size = "Medical conditions", shape = "Sex") +
  geom_abline(intercept = 0, slope = 1, color = "lightgrey") +
  xlim(c(0, 0.8)) +
  ylim(c(0, 0.8)) +
  theme_bw()
p6 <- ggplot(scatt_hosp, aes(Hyp, WU, color = Age_cat, shape = Sex, size = N_Comorbidities)) +
  geom_point(alpha = 0.7)+
  scale_size_discrete(range = c(1,2)) +
  scale_color_viridis_d(direction =  -1) +
  labs(x = "Hypothetical", y = "While Untreated", color = "Age",         
       size = "Medical conditions", shape = "Sex") +
  geom_abline(intercept = 0, slope = 1, color = "lightgrey") +
  xlim(c(0, 0.8)) +
  ylim(c(0, 0.8)) +
  theme_bw()

# Figure 1:
ggpubr::ggarrange(p1, p2, p3, p4, p6, p5, ncol=2, nrow=3, common.legend = TRUE, legend="right")


### Figure 2 -------------------------------------------------------------------
# Choose covariate values for the predictions
agepredict <- c("70-79")                # patient aged 70-79 
sexpredict <- 1                         # Male patient
wavepredict <- 1                        # Second wave
compredict <- 2                         # 2 comorbidities 
therapypredict <- 0                     # non-treated patients
timepredict <- seq(0.01,28, by=0.01) # Days since start of hospitalization
l <- length(timepredict)

# Create dataset for prediction with the chosen values
covs_hosp <- data.frame(Age_cat = factor(rep(agepredict, l), levels = levels(dfr$Age_cat)),
                        Sex_m = sexpredict, 
                        N_Comorbidities = compredict,
                        WaveII = wavepredict,
                        ICU = therapypredict,
                        DaysToDeath = timepredict,
                        Death = rep(0, l),
                        DaysToICUOrDeath = timepredict,
                        tstart = 0,
                        tstop = timepredict,
                        ICUOrDeath = rep(0, l),
                        stat1 = 0,
                        stat2 = 0)
preds_IT_hosp <- 1 - predict(IT_hosp, newdata = covs_hosp, type = "survival")   # predict IT strategy
preds_CS_hosp <- 1 - predict(CS_hosp, newdata = covs_hosp, type = "survival")   # predict CS strategy
preds_Hyp_hosp <- 1 - predict(Hyp_hosp, newdata = covs_hosp, type = "survival") # predict hyp strategy

# Expand the prediction dataset for mstate prediction
data_mstate_plot <- msprep(time = c(NA, "DaysToICUOrDeath", "DaysToICUOrDeath"),
                           status = c(NA, "stat1", "stat2"),
                           data = covs_hosp[nrow(covs_hosp), ],
                           keep = c("Age_cat", "Sex_m", "N_Comorbidities", "WaveII"),
                           trans = tmat)
data_mstate_plot <- expand.covs(data_mstate_plot, c("Age_cat", "Sex_m", "N_Comorbidities", "WaveII"))
data_mstate_plot$strata <- data_mstate_plot$trans

msf <- msfit(WU_hosp, data_mstate_plot, trans = tmat) # predict from multistate model
pt <- probtrans(msf, 0)[[1]] # WU prediction 
# WU predictions are only reported at the times the mortality risk changes. We 
# here expand the predictions to get one predicted risk per desired timepoint 
# (i.e. tiems in "time predict")
WU_hosp_df <- data.frame(time = pt$time[pt$time <= 28], risk = pt$pstate2[pt$time <= 28]) %>%
  merge(data.frame(time = timepredict), all = TRUE) %>% 
  arrange(time) %>%
  fill(risk) %>%
  filter(time %in% timepredict)

# Prepare data frame for ggplot figure
strategies <- c("Ignore Treatment", "Composite", "While Untreated", "Hypothetical")
df_pred_hosp <- data.frame(
  Time = c(rep(timepredict, 4)),
  Strategy  = factor(rep(strategies, each = l), levels = strategies),
  Prediction = c(preds_IT_hosp, preds_CS_hosp, WU_hosp_df$risk, preds_Hyp_hosp)
)

# Figure 2:
ggplot(df_pred_hosp, aes(x = Time, y = Prediction, color = Strategy)) +
  geom_line(aes(linetype = Strategy), size = 1.05) +
  scale_linetype_manual(values=c(4,5,1,6))+
  theme_bw() +
  labs(title = "75 year old male patient with two medical conditions", 
       x = "Days since hospitalization", y = "Mortality risk")





