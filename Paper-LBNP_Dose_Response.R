##########################################
### Preamble
##########################################

rm(list = ls())
graphics.off()
cat("\014")

library(tidyverse)
library(ggpubr)
library(BayesFactor)
library(ggpubr)
library(brms)
library(emmeans)
library(ggeffects)
library(rstan)
library(lme4)
library(lmerTest)
library(ggResidpanel)
library(bayesplot)
library(tidybayes)
library(bayestestR)
library(corrplot)
library(ggraph)
library(igraph)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
options(buildtools.check = function(action) TRUE)

#####################
### Simulate Generative Model
#####################

# sub_n  <- 100 # number of subjects in this simulation
# b_Face <- -0.2
# b_Sex.c <- -0.2
# b_Index <- 0.1
# 
# 
# sim <- tibble(
#   ID = rep(c(1:sub_n),each = 12),
#   sub_i  = rep(rnorm(sub_n, 0, 0.5), each = 12), # random intercept
#   Face = rep(c(0,1), each = 6, times = sub_n),
#   Index = rep(seq(0,5,1), 2*sub_n),
#   Sex.c = rep(c(-0.5,0.5), each = 6*sub_n)
# ) %>%
#   mutate(mu = sub_i + b_Face*Face +b_Sex.c*Sex.c + b_Index*Index) %>%
#   mutate(Y = rstudent_t(n(),5,mu,0.2))
# 
# 
# priors.sim <- c(prior(normal(0, 1), class = b), # slope prior
#                prior(normal(0, 1), class = b, coef = "Intercept"), # intercept prior
#                prior(cauchy(0,1), class = sigma), # population variance
#                prior(cauchy(0,1), class = sd), # group variance
#                prior(gamma(2, .1), class = nu) # degrees of freedom
# )
# 
# model.sim.p <- brm(Y ~ 0 + Intercept + Index + Sex.c + Face + (1|ID),
#                  data = sim,
#                  prior = priors.sim,
#                  sample_prior = "only",
#                  family = student(),
#                  warmup = 1000, # burn-in
#                  iter = 5000, # number of iterations
#                  chains = 4,  # number of MCMC chains
#                  control = list(adapt_delta = 0.95),# advanced MC settings
#                  save_pars=save_pars(all=TRUE))
# 
# pp_check(model.sim.p, ndraws = 50) + coord_cartesian(xlim = c(-10,10))
# 
# model.sim <- brm(Y ~ 0 + Intercept + Index + Sex.c + Face + (1|ID),
#              data = sim,
#              prior = priors.sim,
#              sample_prior = "yes",
#              family = student(),
#              warmup = 1000, # burn-in
#              iter = 5000, # number of iterations
#              chains = 4,  # number of MCMC chains
#              control = list(adapt_delta = 0.95),# advanced MC settings
#              save_pars=save_pars(all=TRUE))
# 
# summary(model.sim)
# pp_check(model.sim, ndraws = 50)

#####################
### Subject Characteristics (Table 1)
#####################

# Read Data
subjects <- read.csv("LBNP-Characteristics.csv")

# Summarize by Sex
summary.stats <- subjects %>% dplyr::group_by(Sex) %>%
  summarise(n = n(),
            mean.Age = mean(Age),
            sd.Age = sd(Age),
            mean.Height = mean(Height),
            sd.Height = sd(Height),
            mean.Weight = mean(Weight),
            sd.Weight = sd(Weight),
            mean.BMI = mean(BMI),
            sd.BMI = sd(BMI),
            mean.HR = mean(HR),
            sd.HR = sd(HR),
            mean.MAP = mean(MAP),
            sd.MAP = sd(MAP))

# Calculate BF10 ([1] = BF(-0.1<d<0.1/d=0),[2] = BF(!(-0.1<d<0.1)/d=0))
Age.bf <- ttestBF(formula = Age ~ Sex, data = subjects, nullInterval = c(-0.1,0.1))
Height.bf <- ttestBF(formula = Height ~ Sex, data = subjects, nullInterval = c(-0.1,0.1))
Weight.bf <- ttestBF(formula = Weight ~ Sex, data = subjects, nullInterval = c(-0.1,0.1))
BMI.bf <- ttestBF(formula = BMI ~ Sex, data = subjects, nullInterval = c(-0.1,0.1))
HR.bf <- ttestBF(formula = HR ~ Sex, data = subjects, nullInterval = c(-0.1,0.1))
MAP.bf <- ttestBF(formula = MAP ~ Sex, data = subjects, nullInterval = c(-0.1,0.1))

# Calculate BF10ROPE: BF(!(-0.1<d<0.1)/(-0.1<d<0.1))
Age.bf10 <- extractBF(Age.bf[2]/Age.bf[1])$bf
Height.bf10 <- extractBF(Height.bf[2]/Height.bf[1])$bf
Weight.bf10 <- extractBF(Weight.bf[2]/Weight.bf[1])$bf
BMI.bf10 <- extractBF(BMI.bf[2]/BMI.bf[1])$bf
HR.bf10 <- extractBF(HR.bf[2]/HR.bf[1])$bf
MAP.bf10 <- extractBF(MAP.bf[2]/MAP.bf[1])$bf

# Bind BFs to Summary Stats
summary.stats <- rbind(summary.stats, list("BF.10.ROPE",NA,Age.bf10,NA,Height.bf10,NA,Weight.bf10,NA,BMI.bf10,NA,
                       HR.bf10,NA,MAP.bf10,NA))

# Print Summary Table
#tibble(summary.stats)

# Extract Latent Variables

subjects.latent <- subjects %>% select(ID,Age,BMI)

#####################
### Variables Considered
#####################

## Systemic:
# Heart Rate *
# Stroke Volume *
# Cardiac Output *
# Oxygen Consumption *
# Systolic Blood Pressure *
# Diastolic Blood Pressure *
# Rate Pressure Product *
# Myocardial Oxygen Supply:Demand Index (DPTI/SPTI) *
# Total Peripheral Resistance *

## Systemic Indices:
# Cardiac Index *
# Stroke Index *

## Cardiovascular Control
# SDNN *
# HRVTi *
# RMSDD *
# BRS *
# LF (Norm) *
# HF (Norm) *
# LF/HF *

## Head/Neck
# Intraocular Pressure *
# Ocular Perfusion Pressure *
# IJV Area
# IJVP *
# IJV Pattern *
# CCA Area
# CCA Peak Systolic Velocity *
# CCA End Diastolic Velocity *


#####################
### Load Data
#####################

full <- read.csv("LBNP_Use.csv")
df <- full %>% filter(Face != "Seated") %>%
  mutate(Sex.c = case_when(Sex == "Male" ~ 0.5,
                                      Sex == "Female" ~ -0.5)) %>%
  mutate(Index = Index - 1)
  #mutate(IJVF.D = IJVF.D - 1,
  #       IJVF.S = IJVF.S - 1) %>%
  #mutate(IJVF.D = factor(IJVF.D, levels = c("1","2","3"), ordered = TRUE),
  #       IJVF.S = factor(IJVF.S, levels = c("1","2","3"), ordered = TRUE))

Pressure_seq <- c("0 mmHg","-10 mmHg","-20 mmHg","-30 mmHg","-40 mmHg","-50 mmHg")
pd <- position_dodge(width = 0.4)

face.labs <- c("0° Supine", "15° HDT")
names(face.labs) <- c("L0", "L15")

side.labs <- c("Left","Right")
names(side.labs) <- c("-0.5","0.5")

int.labs <- c("F, 0° Supine", "M, 0° Supine", "F, 15° HDT", "M, 15° HDT")
names(int.labs) <- c("-0.5.L0","0.5.L0","-0.5.L15","0.5.L15")

nd <- tidyr::crossing(Index = seq(from = 0, to = 5, length.out = 51),
                      Sex.c = -0.5:0.5,
                      Face = c("L0","L15"))

nd2 <- tidyr::crossing(Index = seq(from = 0, to = 5, length.out = 51),
                       Sex.c = -0.5:0.5,
                       Side = -0.5:0.5,
                       Face = c("L0","L15"))

# Scale function to mean = 0, sd = 0.5
scalebayes <- function(x) {
  scale(x, center = TRUE, scale = 2*sd(x,na.rm = TRUE))
}

df.s <- df %>% 
  merge(subjects.latent, by = "ID") %>%
  pivot_longer(IOP.D:IJV.S, names_to = c(".value", "Side"),
               names_sep = "\\.") %>%
  mutate(across(Time.F:PredCO, ~ ifelse(Side == "S",NA,.))) %>%
  mutate(IJVF = factor(IJVF, levels = c("1","2","3"), ordered = TRUE)) %>%
  arrange(Side) %>%
  mutate(Side = case_when(Side == "D" ~ 0.5,
                          Side == "S" ~ -0.5)) %>%
  mutate(HR.s = scalebayes(HR.ECG),
         SV.s = scalebayes(SV),
         CO.s = scalebayes(CO),
         VO2.s = scalebayes(VO2),
         SBP.s = scalebayes(reSYS),
         DBP.s = scalebayes(reDIA),
         RPP.s = scalebayes(RPP),
         MO.s = scalebayes(DPTI.SPTI),
         TPR.calc = reMAP/CO*0.06,
         TPR.s = scalebayes(TPR.calc),
         CI.s = scalebayes(CI),
         SI.s = scalebayes(SI),
         Age.s = scalebayes(Age.y),
         Height.s = scalebayes(Height),
         Weight.s = scalebayes(Weight),
         BMI.s = scalebayes(BMI),
         SDNN.s = scalebayes(SDNN),
         HRVTI.s = scalebayes(HRVTI),
         RMSDD.s = scalebayes(RMSDD),
         BRS.s = scalebayes(BRS),
         LFNorm.s = scalebayes(LFNorm),
         HFNorm.s = scalebayes(HFNorm),
         LFHF.s = scalebayes(LF.HF),
         IOP.s = scalebayes(IOP),
         VPC.s = scalebayes(VPC),
         VP.s = scalebayes(VP),
         PSV.s = scalebayes(PSV),
         EDV.s = scalebayes(EDV),
         OPP.s = scalebayes(OPP),
         CCA.s = scalebayes(CCA),
         IJV.s = scalebayes(IJV))


#### Summarise
summarystats1 <- df.s %>% reframe(across(c(HR.ECG,SV,CO,VO2,reSYS,reDIA,RPP,DPTI.SPTI,TPR.calc,SI,CI),
                                list(mean = \(x) mean(x, na.rm = TRUE),
                                     se =\(x) sd(x, na.rm = TRUE)/sqrt(sum(!is.na(x))))),
                         .by = c(Index, Sex, Face)) %>%
  arrange(desc(Sex),Face,Index)
summarystats1t <- df.s %>% reframe(across(c(HR.ECG,SV,CO,VO2,reSYS,reDIA,RPP,DPTI.SPTI,TPR.calc,SI,CI),
                                         list(mean = \(x) mean(x, na.rm = TRUE),
                                              se =\(x) sd(x, na.rm = TRUE)/sqrt(sum(!is.na(x))))),
                                  .by = c(Index)) %>%
  arrange(Index)
summarystats2 <- df.s %>% reframe(across(c(SDNN,HRVTI,RMSDD,BRS,LFNorm,HFNorm,LF.HF),
                                         list(mean = \(x) mean(x, na.rm = TRUE),
                                              se =\(x) sd(x, na.rm = TRUE)/sqrt(sum(!is.na(x))))),
                                  .by = c(Index, Sex, Face)) %>%
  arrange(desc(Sex),Face,Index)
summarystats2t <- df.s %>% reframe(across(c(SDNN,HRVTI,RMSDD,BRS,LFNorm,HFNorm,LF.HF),
                                         list(mean = \(x) mean(x, na.rm = TRUE),
                                              se =\(x) sd(x, na.rm = TRUE)/sqrt(sum(!is.na(x))))),
                                  .by = c(Index)) %>%
  arrange(Index)
test <- df.s %>% reframe(across(c(SDNN,HRVTI,RMSDD,BRS,LFNorm,HFNorm,LF.HF),
                                         list(mean = \(x) mean(x, na.rm = TRUE),
                                              se =\(x) sd(x, na.rm = TRUE)/sqrt(sum(!is.na(x))))),
                                  .by = c(Index)) %>%
  arrange(Index)
summarystats3 <- df.s %>% reframe(across(c(IOP,OPP,IJV,VPC,CCA,PSV,EDV),
                                         list(mean = \(x) mean(x, na.rm = TRUE),
                                              se =\(x) sd(x, na.rm = TRUE)/sqrt(sum(!is.na(x))))),
                                  .by = c(Index, Sex, Face, Side)) %>%
  arrange(desc(Sex),Face,Index,desc(Side))

test <- df.s %>% reframe(across(c(IOP,OPP,IJV,VPC,CCA,PSV,EDV),
                                         list(mean = \(x) mean(x, na.rm = TRUE),
                                              se =\(x) sd(x, na.rm = TRUE)/sqrt(sum(!is.na(x))))),
                                  .by = c(Index, Face)) %>%
  arrange(Index, Face)

test <- df.s %>% reframe(across(c(IOP,OPP,IJV,VPC,CCA,PSV,EDV),
                                list(mean = \(x) mean(x, na.rm = TRUE),
                                     se =\(x) sd(x, na.rm = TRUE)/sqrt(sum(!is.na(x))))),
                         .by = c(Face)) %>%
  arrange(Face)

##########################################
### Systemic Hemodynamics
##########################################
#####################
### Heart Rate
#####################
# Find scales
scale.hr <- attributes(df.s$HR.s)$`scaled:scale`
center.hr <- attributes(df.s$HR.s)$`scaled:center`
sd.hr <- sd(df$HR.ECG[!is.na(df$HR.ECG)])

# Plot Data
data.hr <- ggplot(data = df, aes(x = Index, color = as.factor(Sex), linetype = Face, y = HR.ECG)) +
  #geom_jitter(alpha = 0.2) +
  stat_summary(aes(group = interaction(Sex,Face)), fun=mean, geom="line", position = pd) +
  stat_summary(aes(group = interaction(Sex,Face)), fun.data = mean_se, geom="errorbar", width = 0.4, position = pd, linetype = "solid") +
  stat_summary(aes(group = interaction(Sex,Face)), fun=mean, geom="point", position = pd) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
  scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq) +
  scale_y_continuous(name = "Heart Rate (bpm)") +
  scale_color_manual(name = "Sex:", values = c("Female" = "hotpink", "Male" = "royalblue"),
                     breaks = c("Male","Female")) +
  scale_linetype_manual(name = "Position:", values = c("L0" = "solid", "L15" = "dashed"),
                        labels = c("L0" = "0° Supine", "L15" = "15° HDT")) +
  theme(legend.position = "right")
# data.hr

# Define Priors
# priors.hr <- c(prior(normal(0, 1), class = b), # slope prior
#                prior(normal(0, 1), class = b, coef = "Intercept"), # intercept prior
#                prior(cauchy(0,1), class = sigma), # population variance
#                prior(cauchy(0,1), class = sd), # group variance
#                prior(gamma(2, .1), class = nu) # degrees of freedom
# )

# # Fit Model
# model.hr <- brm(HR.s ~ 0 + Intercept + Index + Sex.c + Face + (1|ID),
#              data = df.s,
#              prior = priors.hr,
#              sample_prior = "yes",
#              family = student(),
#              warmup = 1000, # burn-in
#              iter = 5000, # number of iterations
#              chains = 4,  # number of MCMC chains
#              control = list(adapt_delta = 0.95),# advanced MC settings
#              save_pars=save_pars(all=TRUE))
# saveRDS(model.hr,file = "BRMSFits/model.hr.rds")
# model.hr <- readRDS("BRMSFits/model.hr.rds")

# Posterior Checks
# summary(model.hr)
# plot(model.hr, N = 7)
# pp_check(model.hr, ndraws = 50)
# model.hr.loo <- loo(model.hr, moment_match = TRUE, save_psis = TRUE)
# model.hr.loo
# model.hr.w <- weights(model.hr.loo$psis_object)
# ppc_loo_pit_overlay(df.s$HR.s[!is.na(df.s$HR.s)], 
#                     posterior_predict(model.hr), 
#                     lw = model.hr.w)

# Describe Posterior
# describe_posterior(model.hr,
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.05, 0.05), rope_ci = 1,
#                    diagnostic = NULL, effects = "fixed", component = "all")
# 
# # Rescale Posterior
# rescale.hr <- model.hr %>% spread_draws(`b_.*`, sigma, `sd_.*`, nu, regex = TRUE) %>%
#   mutate(b_Intercept = (b_Intercept * scale.hr) + center.hr,
#          b_Index = b_Index * scale.hr,
#          b_Sex.c = b_Sex.c * scale.hr,
#          b_FaceL15 = b_FaceL15 * scale.hr,
#          sigma = sigma * scale.hr,
#          sd_ID__Intercept = sd_ID__Intercept * scale.hr)
# describe_posterior(rescale.hr[,-c(1:3)],
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.1*sd.hr, 0.1*sd.hr), rope_ci = 1,
#                    diagnostic = NULL)

# Plot Model
# model.hr.dr <- ggemmeans(model.hr, terms = c("Index","Sex.c","Face"), typical = "mode", ppd = TRUE, ci.lvl = 0.89) %>%
#   mutate(predicted = (predicted * scale.hr) + center.hr,
#          conf.low = (conf.low * scale.hr) + center.hr,
#          conf.high = (conf.high * scale.hr) + center.hr)
# 
# dose.hr <-ggplot(data = model.hr.dr, aes(x = x, fill = group)) +
#   geom_line(aes(color = group, y = predicted)) +
#   geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
#   facet_wrap(~ facet, labeller = labeller(facet = face.labs), scales = "free") +
#   theme_classic2() +
#   theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
#   scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq, expand = c(0,0)) +
#   scale_y_continuous(name = "Heart Rate (bpm)") +
#   scale_color_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                      labels = c("-0.5" = "Female", "0.5" = "Male"),
#                      breaks = c("0.5","-0.5")) +
#   scale_fill_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                     labels = c("-0.5" = "Female", "0.5" = "Male"),
#                     breaks = c("0.5","-0.5")) +
#   theme(legend.position = "bottom", panel.spacing = unit(1,"lines")) +
#   coord_cartesian(ylim = c(50,100), expand = FALSE)
# dose.hr

# model.hr.dr <- fitted(model.hr, newdata = nd, re_formula = NA, probs = c(.055,.945)) %>%
#   data.frame() %>%
#   bind_cols(nd) %>%
#   mutate(Estimate.r = (Estimate * scale.hr) + center.hr,
#          Q5.5.r = (Q5.5 * scale.hr) + center.hr,
#          Q94.5.r = (Q94.5 * scale.hr) + center.hr)
# 
# dose.hr <- ggplot(model.hr.dr, aes(x = Index, color = factor(Sex.c), fill = factor(Sex.c))) +
#   geom_smooth(aes(y = Estimate.r, ymin = Q5.5.r, ymax = Q94.5.r),
#               stat = "identity",
#               alpha = 0.2) +
#   facet_wrap(~ Face, labeller = labeller(Face = face.labs), scales = "free") +
#   theme_classic2() +
#   theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
#   scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq, expand = c(0,0)) +
#   scale_y_continuous(name = "Heart Rate (bpm)") +
#   scale_color_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                      labels = c("-0.5" = "Female", "0.5" = "Male"),
#                      breaks = c("0.5","-0.5")) +
#   scale_fill_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                     labels = c("-0.5" = "Female", "0.5" = "Male"),
#                     breaks = c("0.5","-0.5")) +
#   theme(legend.position = "bottom", panel.spacing = unit(1,"lines")) +
#   coord_cartesian(ylim = c(50,100), expand = FALSE)
# dose.hr

#####################
### Stroke Volume
#####################

# Find scales
scale.sv <- attributes(df.s$SV.s)$`scaled:scale`
center.sv <- attributes(df.s$SV.s)$`scaled:center`
sd.sv <- sd(df$SV[!is.na(df$SV)])

# Plot Data
data.sv <- ggplot(data = df, aes(x = Index, color = as.factor(Sex), linetype = Face, y = SV)) +
  #geom_jitter(alpha = 0.2) +
  stat_summary(aes(group = interaction(Sex,Face)), fun=mean, geom="line", position = pd) +
  stat_summary(aes(group = interaction(Sex,Face)), fun.data = mean_se, geom="errorbar", width = 0.4, position = pd, linetype = "solid") +
  stat_summary(aes(group = interaction(Sex,Face)), fun=mean, geom="point", position = pd) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
  scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq) +
  scale_y_continuous(name = "Stroke Volume (ml)") +
  scale_color_manual(name = "Sex:", values = c("Female" = "hotpink", "Male" = "royalblue"),
                     breaks = c("Male","Female")) +
  scale_linetype_manual(name = "Position:", values = c("L0" = "solid", "L15" = "dashed"),
                        labels = c("L0" = "0° Supine", "L15" = "15° HDT")) +
  theme(legend.position = "right")
# data.sv

# # Define Priors
# priors.sv <- c(prior(normal(0, 1), class = b), # slope prior
#                prior(normal(0, 1), class = b, coef = "Intercept"), # intercept prior
#                prior(cauchy(0,1), class = sigma), # population variance
#                prior(cauchy(0,1), class = sd), # group variance
#                prior(gamma(2, .1), class = nu) # degrees of freedom
# )
# 
# # # Fit Model
# # model.sv <- brm(SV.s ~ 0 + Intercept + Index + Sex.c + Face + Index:Sex.c + (1|ID),
# #              data = df.s,
# #              prior = priors.sv,
# #              sample_prior = "yes",
# #              family = student(),
# #              warmup = 1000, # burn-in
# #              iter = 5000, # number of iterations
# #              chains = 4,  # number of MCMC chains
# #              control = list(adapt_delta = 0.95),# advanced MC settings
# #              save_pars=save_pars(all=TRUE))
# # saveRDS(model.sv,file = "BRMSFits/model.sv.rds")
# model.sv <- readRDS("BRMSFits/model.sv.rds")
# 
# # Posterior Checks
# summary(model.sv)
# plot(model.sv, N = 8)
# pp_check(model.sv, ndraws = 50)
# model.sv.loo <- loo(model.sv, moment_match = TRUE, save_psis = TRUE)
# model.sv.loo
# model.sv.w <- weights(model.sv.loo$psis_object)
# ppc_loo_pit_overlay(df.s$SV.s[!is.na(df.s$SV.s)], 
#                     posterior_predict(model.sv), 
#                     lw = model.sv.w)
# 
# # Describe Posterior
# describe_posterior(model.sv,
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.05, 0.05), rope_ci = 1,
#                    diagnostic = NULL, effects = "fixed", component = "all")
# 
# # Rescale Posterior
# rescale.sv <- model.sv %>% spread_draws(`b_.*`, sigma, `sd_.*`, nu, regex = TRUE) %>%
#   mutate(b_Intercept = (b_Intercept * scale.sv) + center.sv,
#          b_Index = b_Index * scale.sv,
#          b_Sex.c = b_Sex.c * scale.sv,
#          b_FaceL15 = b_FaceL15 * scale.sv,
#          `b_Index:Sex.c` = `b_Index:Sex.c` * scale.sv,
#          sigma = sigma * scale.sv,
#          sd_ID__Intercept = sd_ID__Intercept * scale.sv)
# describe_posterior(rescale.sv[,-c(1:3)],
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.1*sd.sv, 0.1*sd.sv), rope_ci = 1,
#                    diagnostic = NULL)
# 
# # Plot Model
# model.sv.dr <- fitted(model.sv, newdata = nd, re_formula = NA, probs = c(.055,.945)) %>%
#   data.frame() %>%
#   bind_cols(nd) %>%
#   mutate(Estimate.r = (Estimate * scale.sv) + center.sv,
#          Q5.5.r = (Q5.5 * scale.sv) + center.sv,
#          Q94.5.r = (Q94.5 * scale.sv) + center.sv)
# 
# dose.sv <- ggplot(model.sv.dr, aes(x = Index, color = factor(Sex.c), fill = factor(Sex.c))) +
#   geom_smooth(aes(y = Estimate.r, ymin = Q5.5.r, ymax = Q94.5.r),
#               stat = "identity",
#               alpha = 0.2) +
#   facet_wrap(~ Face, labeller = labeller(Face = face.labs), scales = "free") +
#   theme_classic2() +
#   theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
#   scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq, expand = c(0,0)) +
#   scale_y_continuous(name = "Stroke Volume (ml)") +
#   scale_color_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                      labels = c("-0.5" = "Female", "0.5" = "Male"),
#                      breaks = c("0.5","-0.5")) +
#   scale_fill_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                     labels = c("-0.5" = "Female", "0.5" = "Male"),
#                     breaks = c("0.5","-0.5")) +
#   theme(legend.position = "bottom", panel.spacing = unit(1,"lines")) +
#   coord_cartesian(ylim = c(0,120), expand = FALSE)
# dose.sv

#####################
### Cardiac Output
#####################

# Find scales
scale.co <- attributes(df.s$CO.s)$`scaled:scale`
center.co <- attributes(df.s$CO.s)$`scaled:center`
sd.co <- sd(df$CO[!is.na(df$CO)])

# Plot Data
data.co <- ggplot(data = df, aes(x = Index, color = as.factor(Sex), linetype = Face, y = CO)) +
  #geom_jitter(alpha = 0.2) +
  stat_summary(aes(group = interaction(Sex,Face)), fun=mean, geom="line", position = pd) +
  stat_summary(aes(group = interaction(Sex,Face)), fun.data = mean_se, geom="errorbar", width = 0.4, position = pd, linetype = "solid") +
  stat_summary(aes(group = interaction(Sex,Face)), fun=mean, geom="point", position = pd) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
  scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq) +
  scale_y_continuous(name = "Cardiac Output (l/min)") +
  scale_color_manual(name = "Sex:", values = c("Female" = "hotpink", "Male" = "royalblue"),
                     breaks = c("Male","Female")) +
  scale_linetype_manual(name = "Position:", values = c("L0" = "solid", "L15" = "dashed"),
                          labels = c("L0" = "0° Supine", "L15" = "15° HDT")) +
  theme(legend.position = "right")
# data.co

# # Define Priors
# priors.co <- c(prior(normal(0, 1), class = b), # slope prior
#             prior(normal(0, 1), class = b, coef = "Intercept"), # intercept prior
#             prior(cauchy(0,1), class = sigma), # population variance
#             prior(cauchy(0,1), class = sd), # group variance
#             prior(gamma(2, .1), class = nu) # degrees of freedom
# )
# 
# # # Fit Model
# # model.co <- brm(CO.s ~ 0 + Intercept + Index + Sex.c + Face + (1|ID),
# #              data = df.s,
# #              prior = priors.co,
# #              sample_prior = "yes",
# #              family = student(),
# #              warmup = 1000, # burn-in
# #              iter = 5000, # number of iterations
# #              chains = 4,  # number of MCMC chains
# #              control = list(adapt_delta = 0.95),# advanced MC settings
# #              save_pars=save_pars(all=TRUE))
# # model.co <- brm(CO.s ~ 0 + Intercept + Index + Sex.c + Face + (1|ID),
# #              data = df.s,
# #              prior = priors.co,
# #              sample_prior = "yes",
# #              family = student(),
# #              warmup = 1000, # burn-in
# #              iter = 5000, # number of iterations
# #              chains = 4,  # number of MCMC chains
# #              control = list(adapt_delta = 0.95),# advanced MC settings
# #              save_pars=save_pars(all=TRUE))
# # saveRDS(model.co,file = "BRMSFits/model.co.rds")
# model.co <- readRDS("BRMSFits/model.co.rds")
# 
# # Posterior Checks
# summary(model.co)
# plot(model.co, N = 7)
# pp_check(model.co, ndraws = 50)
# model.co.loo <- loo(model.co, moment_match = TRUE, save_psis = TRUE)
# model.co.loo
# model.co.w <- weights(model.co.loo$psis_object)
# ppc_loo_pit_overlay(df.s$CO.s[!is.na(df.s$CO.s)], 
#                     posterior_predict(model.co), 
#                     lw = model.co.w)
# 
# # Describe Posterior
# describe_posterior(model.co,
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.05, 0.05), rope_ci = 1,
#                    diagnostic = NULL, effects = "fixed", component = "all")
# 
# # Rescale Posterior
# rescale.co <- model.co %>% spread_draws(`b_.*`, sigma, `sd_.*`, nu, regex = TRUE) %>%
#   mutate(b_Intercept = (b_Intercept * scale.co) + center.co,
#          b_Index = b_Index * scale.co,
#          b_Sex.c = b_Sex.c * scale.co,
#          b_FaceL15 = b_FaceL15 * scale.co,
#          sigma = sigma * scale.co,
#          sd_ID__Intercept = sd_ID__Intercept * scale.co)
# describe_posterior(rescale.co[,-c(1:3)],
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.1*sd.co, 0.1*sd.co), rope_ci = 1,
#                    diagnostic = NULL)
# 
# # Plot Model
# model.co.dr <- fitted(model.co, newdata = nd, re_formula = NA, probs = c(.055,.945)) %>%
#   data.frame() %>%
#   bind_cols(nd) %>%
#   mutate(Estimate.r = (Estimate * scale.co) + center.co,
#          Q5.5.r = (Q5.5 * scale.co) + center.co,
#          Q94.5.r = (Q94.5 * scale.co) + center.co)
# 
# dose.co <- ggplot(model.co.dr, aes(x = Index, color = factor(Sex.c), fill = factor(Sex.c))) +
#   geom_smooth(aes(y = Estimate.r, ymin = Q5.5.r, ymax = Q94.5.r),
#               stat = "identity",
#               alpha = 0.2) +
#   facet_wrap(~ Face, labeller = labeller(Face = face.labs), scales = "free") +
#   theme_classic2() +
#   theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
#   scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq, expand = c(0,0)) +
#   scale_y_continuous(name = "Cardiac Output (l/min)") +
#   scale_color_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                      labels = c("-0.5" = "Female", "0.5" = "Male"),
#                      breaks = c("0.5","-0.5")) +
#   scale_fill_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                     labels = c("-0.5" = "Female", "0.5" = "Male"),
#                     breaks = c("0.5","-0.5")) +
#   theme(legend.position = "bottom", panel.spacing = unit(1,"lines")) +
#   coord_cartesian(ylim = c(0,8), expand = FALSE)
# dose.co

#####################
### Oxygen Consumption
#####################

# Find scales
scale.vo2 <- attributes(df.s$VO2.s)$`scaled:scale`
center.vo2 <- attributes(df.s$VO2.s)$`scaled:center`
sd.vo2 <- sd(df$VO2[!is.na(df$VO2)])

# Plot Data
data.vo2 <- ggplot(data = df, aes(x = Index, color = as.factor(Sex), linetype = Face, y = VO2)) +
  #geom_jitter(alpha = 0.2) +
  stat_summary(aes(group = interaction(Sex,Face)), fun=mean, geom="line", position = pd) +
  stat_summary(aes(group = interaction(Sex,Face)), fun.data = mean_se, geom="errorbar", width = 0.4, position = pd, linetype = "solid") +
  stat_summary(aes(group = interaction(Sex,Face)), fun=mean, geom="point", position = pd) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
  scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq) +
  scale_y_continuous(name = "Oxygen Consumption (l/min)") +
  scale_color_manual(name = "Sex:", values = c("Female" = "hotpink", "Male" = "royalblue"),
                     breaks = c("Male","Female")) +
  scale_linetype_manual(name = "Position:", values = c("L0" = "solid", "L15" = "dashed"),
                        labels = c("L0" = "0° Supine", "L15" = "15° HDT")) +
  theme(legend.position = "right")
# data.vo2

# # Define Priors
# priors.vo2 <- c(prior(normal(0, 1), class = b), # slope prior
#                prior(normal(0, 1), class = b, coef = "Intercept"), # intercept prior
#                prior(cauchy(0,1), class = sigma), # population variance
#                prior(cauchy(0,1), class = sd), # group variance
#                prior(gamma(2, .1), class = nu) # degrees of freedom
# )
# 
# model.vo2 <- lmer(VO2.s ~ Index + Sex.c + Face +
#                     Index:Sex.c + Sex.c:Face +
#                     (1|ID), data = df.s, REML = TRUE)
# 
# # # Fit Model
# # model.vo2 <- brm(VO2.s ~ 0 + Intercept + Index + Sex.c + Face + Sex.c:Face + (1|ID),
# #              data = df.s,
# #              prior = priors.vo2,
# #              sample_prior = "yes",
# #              family = student(),
# #              warmup = 1000, # burn-in
# #              iter = 5000, # number of iterations
# #              chains = 4,  # number of MCMC chains
# #              control = list(adapt_delta = 0.95),# advanced MC settings
# #              save_pars=save_pars(all=TRUE))
# # saveRDS(model.vo2,file = "BRMSFits/model.vo2.rds")
# model.vo2 <- readRDS("BRMSFits/model.vo2.rds")
# 
# # Posterior Checks
# summary(model.vo2)
# plot(model.vo2, N = 8)
# pp_check(model.vo2, ndraws = 50)
# model.vo2.loo <- loo(model.vo2, moment_match = TRUE, save_psis = TRUE)
# model.vo2.loo
# model.vo2.w <- weights(model.vo2.loo$psis_object)
# ppc_loo_pit_overlay(df.s$VO2.s[!is.na(df.s$VO2.s)], 
#                     posterior_predict(model.vo2), 
#                     lw = model.vo2.w)
# 
# # Describe Posterior
# describe_posterior(model.vo2,
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.05, 0.05), rope_ci = 1,
#                    diagnostic = NULL, effects = "fixed", component = "all")
# 
# # Rescale Posterior
# rescale.vo2 <- model.vo2 %>% spread_draws(`b_.*`, sigma, `sd_.*`, nu, regex = TRUE) %>%
#   mutate(b_Intercept = (b_Intercept * scale.vo2) + center.vo2,
#          b_Index = b_Index * scale.vo2,
#          b_Sex.c = b_Sex.c * scale.vo2,
#          b_FaceL15 = b_FaceL15 * scale.vo2,
#          `b_Sex.c:FaceL15` = `b_Sex.c:FaceL15` * scale.vo2,
#          sigma = sigma * scale.vo2,
#          sd_ID__Intercept = sd_ID__Intercept * scale.vo2)
# describe_posterior(rescale.vo2[,-c(1:3)],
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.1*sd.vo2, 0.1*sd.vo2), rope_ci = 1,
#                    diagnostic = NULL)
# 
# # Plot Model
# model.vo2.dr <- fitted(model.vo2, newdata = nd, re_formula = NA, probs = c(.055,.945)) %>%
#   data.frame() %>%
#   bind_cols(nd) %>%
#   mutate(Estimate.r = (Estimate * scale.vo2) + center.vo2,
#          Q5.5.r = (Q5.5 * scale.vo2) + center.vo2,
#          Q94.5.r = (Q94.5 * scale.vo2) + center.vo2)
# 
# dose.vo2 <- ggplot(model.vo2.dr, aes(x = Index, color = factor(Sex.c), fill = factor(Sex.c))) +
#   geom_smooth(aes(y = Estimate.r, ymin = Q5.5.r, ymax = Q94.5.r),
#               stat = "identity",
#               alpha = 0.2) +
#   facet_wrap(~ Face, labeller = labeller(Face = face.labs), scales = "free") +
#   theme_classic2() +
#   theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
#   scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq, expand = c(0,0)) +
#   scale_y_continuous(name = "Oxygen Consumption (l/min)") +
#   scale_color_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                      labels = c("-0.5" = "Female", "0.5" = "Male"),
#                      breaks = c("0.5","-0.5")) +
#   scale_fill_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                     labels = c("-0.5" = "Female", "0.5" = "Male"),
#                     breaks = c("0.5","-0.5")) +
#   theme(legend.position = "bottom", panel.spacing = unit(1,"lines")) +
#   coord_cartesian(ylim = c(0,0.3), expand = FALSE)
# dose.vo2

#####################
### Systolic Blood Pressure
#####################

# Find scales
scale.sbp <- attributes(df.s$SBP.s)$`scaled:scale`
center.sbp <- attributes(df.s$SBP.s)$`scaled:center`
sd.sbp <- sd(df$reSYS[!is.na(df$reSYS)])

# Plot Data
data.sbp <- ggplot(data = df, aes(x = Index, color = as.factor(Sex), linetype = Face, y = reSYS)) +
  #geom_jitter(alpha = 0.2) +
  stat_summary(aes(group = interaction(Sex,Face)), fun=mean, geom="line", position = pd) +
  stat_summary(aes(group = interaction(Sex,Face)), fun.data = mean_se, geom="errorbar", width = 0.4, position = pd, linetype = "solid") +
  stat_summary(aes(group = interaction(Sex,Face)), fun=mean, geom="point", position = pd) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
  scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq) +
  scale_y_continuous(name = "Systolic Blood Pressure (mmHg)") +
  scale_color_manual(name = "Sex:", values = c("Female" = "hotpink", "Male" = "royalblue"),
                     breaks = c("Male","Female")) +
  scale_linetype_manual(name = "Position:", values = c("L0" = "solid", "L15" = "dashed"),
                        labels = c("L0" = "0° Supine", "L15" = "15° HDT")) +
  theme(legend.position = "right")
# data.sbp

# # Define Priors
# priors.sbp <- c(prior(normal(0, 1), class = b), # slope prior
#                 prior(normal(0, 1), class = b, coef = "Intercept"), # intercept prior
#                 prior(cauchy(0,1), class = sigma), # population variance
#                 prior(cauchy(0,1), class = sd), # group variance
#                 prior(gamma(2, .1), class = nu) # degrees of freedom
# )
# 
# # # Fit Model
# # model.sbp <- brm(SBP.s ~ 0 + Intercept + Index + Sex.c + Face + (1|ID),
# #              data = df.s,
# #              prior = priors.sbp,
# #              sample_prior = "yes",
# #              family = student(),
# #              warmup = 1000, # burn-in
# #              iter = 5000, # number of iterations
# #              chains = 4,  # number of MCMC chains
# #              control = list(adapt_delta = 0.95),# advanced MC settings
# #              save_pars=save_pars(all=TRUE))
# # saveRDS(model.sbp,file = "BRMSFits/model.sbp.rds")
# model.sbp <- readRDS("BRMSFits/model.sbp.rds")
# 
# # Posterior Checks
# summary(model.sbp)
# plot(model.sbp, N = 7)
# pp_check(model.sbp, ndraws = 50)
# model.sbp.loo <- loo(model.sbp, moment_match = TRUE, save_psis = TRUE)
# model.sbp.loo
# model.sbp.w <- weights(model.sbp.loo$psis_object)
# ppc_loo_pit_overlay(df.s$SBP.s[!is.na(df.s$SBP.s)], 
#                     posterior_predict(model.sbp), 
#                     lw = model.sbp.w)
# 
# # Describe Posterior
# describe_posterior(model.sbp,
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.05, 0.05), rope_ci = 1,
#                    diagnostic = NULL, effects = "fixed", component = "all")
# 
# # Rescale Posterior
# rescale.sbp <- model.sbp %>% spread_draws(`b_.*`, sigma, `sd_.*`, nu, regex = TRUE) %>%
#   mutate(b_Intercept = (b_Intercept * scale.sbp) + center.sbp,
#          b_Index = b_Index * scale.sbp,
#          b_Sex.c = b_Sex.c * scale.sbp,
#          b_FaceL15 = b_FaceL15 * scale.sbp,
#          sigma = sigma * scale.sbp,
#          sd_ID__Intercept = sd_ID__Intercept * scale.sbp)
# describe_posterior(rescale.sbp[,-c(1:3)],
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.1*sd.sbp, 0.1*sd.sbp), rope_ci = 1,
#                    diagnostic = NULL)
# 
# # Plot Model
# model.sbp.dr <- fitted(model.sbp, newdata = nd, re_formula = NA, probs = c(.055,.945)) %>%
#   data.frame() %>%
#   bind_cols(nd) %>%
#   mutate(Estimate.r = (Estimate * scale.sbp) + center.sbp,
#          Q5.5.r = (Q5.5 * scale.sbp) + center.sbp,
#          Q94.5.r = (Q94.5 * scale.sbp) + center.sbp)
# 
# dose.sbp <- ggplot(model.sbp.dr, aes(x = Index, color = factor(Sex.c), fill = factor(Sex.c))) +
#   geom_smooth(aes(y = Estimate.r, ymin = Q5.5.r, ymax = Q94.5.r),
#               stat = "identity",
#               alpha = 0.2) +
#   facet_wrap(~ Face, labeller = labeller(Face = face.labs), scales = "free") +
#   theme_classic2() +
#   theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
#   scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq, expand = c(0,0)) +
#   scale_y_continuous(name = "Systolic Blood Pressure (mmHg)") +
#   scale_color_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                      labels = c("-0.5" = "Female", "0.5" = "Male"),
#                      breaks = c("0.5","-0.5")) +
#   scale_fill_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                     labels = c("-0.5" = "Female", "0.5" = "Male"),
#                     breaks = c("0.5","-0.5")) +
#   theme(legend.position = "bottom", panel.spacing = unit(1,"lines")) +
#   coord_cartesian(ylim = c(100,150), expand = FALSE)
# dose.sbp

#####################
### Diastolic Blood Pressure
#####################

# Find scales
scale.dbp <- attributes(df.s$DBP.s)$`scaled:scale`
center.dbp <- attributes(df.s$DBP.s)$`scaled:center`
sd.dbp <- sd(df$reDIA[!is.na(df$reDIA)])

# Plot Data
data.dbp <- ggplot(data = df, aes(x = Index, color = as.factor(Sex), linetype = Face, y = reDIA)) +
  #geom_jitter(alpha = 0.2) +
  stat_summary(aes(group = interaction(Sex,Face)), fun=mean, geom="line", position = pd) +
  stat_summary(aes(group = interaction(Sex,Face)), fun.data = mean_se, geom="errorbar", width = 0.4, position = pd, linetype = "solid") +
  stat_summary(aes(group = interaction(Sex,Face)), fun=mean, geom="point", position = pd) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
  scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq) +
  scale_y_continuous(name = "Diastolic Blood Pressure (mmHg)") +
  scale_color_manual(name = "Sex:", values = c("Female" = "hotpink", "Male" = "royalblue"),
                     breaks = c("Male","Female")) +
  scale_linetype_manual(name = "Position:", values = c("L0" = "solid", "L15" = "dashed"),
                        labels = c("L0" = "0° Supine", "L15" = "15° HDT")) +
  theme(legend.position = "right")
# data.dbp

# # Define Priors
# priors.dbp <- c(prior(normal(0, 1), class = b), # slope prior
#                 prior(normal(0, 1), class = b, coef = "Intercept"), # intercept prior
#                 prior(cauchy(0,1), class = sigma), # population variance
#                 prior(cauchy(0,1), class = sd), # group variance
#                 prior(gamma(2, .1), class = nu) # degrees of freedom
# )
# 
# # # Fit Model
# # model.dbp <- brm(DBP.s ~ 0 + Intercept + Index + Sex.c + Face + (1|ID),
# #              data = df.s,
# #              prior = priors.dbp,
# #              sample_prior = "yes",
# #              family = student(),
# #              warmup = 1000, # burn-in
# #              iter = 5000, # number of iterations
# #              chains = 4,  # number of MCMC chains
# #              control = list(adapt_delta = 0.95),# advanced MC settings
# #              save_pars=save_pars(all=TRUE))
# # saveRDS(model.dbp,file = "BRMSFits/model.dbp.rds")
# model.dbp <- readRDS("BRMSFits/model.dbp.rds")
# 
# # Posterior Checks
# summary(model.dbp)
# plot(model.dbp, N = 7)
# pp_check(model.dbp, ndraws = 50)
# model.dbp.loo <- loo(model.dbp, moment_match = TRUE, save_psis = TRUE)
# model.dbp.loo
# model.dbp.w <- weights(model.dbp.loo$psis_object)
# ppc_loo_pit_overlay(df.s$DBP.s[!is.na(df.s$DBP.s)], 
#                     posterior_predict(model.dbp), 
#                     lw = model.dbp.w)
# 
# # Describe Posterior
# describe_posterior(model.dbp,
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.05, 0.05), rope_ci = 1,
#                    diagnostic = NULL, effects = "fixed", component = "all")
# 
# # Rescale Posterior
# rescale.dbp <- model.dbp %>% spread_draws(`b_.*`, sigma, `sd_.*`, nu, regex = TRUE) %>%
#   mutate(b_Intercept = (b_Intercept * scale.dbp) + center.dbp,
#          b_Index = b_Index * scale.dbp,
#          b_Sex.c = b_Sex.c * scale.dbp,
#          b_FaceL15 = b_FaceL15 * scale.dbp,
#          sigma = sigma * scale.dbp,
#          sd_ID__Intercept = sd_ID__Intercept * scale.dbp)
# describe_posterior(rescale.dbp[,-c(1:3)],
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.1*sd.dbp, 0.1*sd.dbp), rope_ci = 1,
#                    diagnostic = NULL)
# 
# # Plot Model
# model.dbp.dr <- fitted(model.dbp, newdata = nd, re_formula = NA, probs = c(.055,.945)) %>%
#   data.frame() %>%
#   bind_cols(nd) %>%
#   mutate(Estimate.r = (Estimate * scale.dbp) + center.dbp,
#          Q5.5.r = (Q5.5 * scale.dbp) + center.dbp,
#          Q94.5.r = (Q94.5 * scale.dbp) + center.dbp)
# 
# dose.dbp <- ggplot(model.dbp.dr, aes(x = Index, color = factor(Sex.c), fill = factor(Sex.c))) +
#   geom_smooth(aes(y = Estimate.r, ymin = Q5.5.r, ymax = Q94.5.r),
#               stat = "identity",
#               alpha = 0.2) +
#   facet_wrap(~ Face, labeller = labeller(Face = face.labs), scales = "free") +
#   theme_classic2() +
#   theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
#   scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq, expand = c(0,0)) +
#   scale_y_continuous(name = "Diastolic Blood Pressure (mmHg)") +
#   scale_color_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                      labels = c("-0.5" = "Female", "0.5" = "Male"),
#                      breaks = c("0.5","-0.5")) +
#   scale_fill_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                     labels = c("-0.5" = "Female", "0.5" = "Male"),
#                     breaks = c("0.5","-0.5")) +
#   theme(legend.position = "bottom", panel.spacing = unit(1,"lines")) +
#   coord_cartesian(ylim = c(50,100), expand = FALSE)
# dose.dbp

#####################
### Rate Pressure Product
#####################

# Find scales
scale.rpp <- attributes(df.s$RPP.s)$`scaled:scale`
center.rpp <- attributes(df.s$RPP.s)$`scaled:center`
sd.rpp <- sd(df$RPP[!is.na(df$RPP)])

# Plot Data
data.rpp <- ggplot(data = df, aes(x = Index, color = as.factor(Sex), linetype = Face, y = RPP)) +
  #geom_jitter(alpha = 0.2) +
  stat_summary(aes(group = interaction(Sex,Face)), fun=mean, geom="line", position = pd) +
  stat_summary(aes(group = interaction(Sex,Face)), fun.data = mean_se, geom="errorbar", width = 0.4, position = pd, linetype = "solid") +
  stat_summary(aes(group = interaction(Sex,Face)), fun=mean, geom="point", position = pd) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
  scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq) +
  scale_y_continuous(name = "Rate Pressure Product (mmHg/min)") +
  scale_color_manual(name = "Sex:", values = c("Female" = "hotpink", "Male" = "royalblue"),
                     breaks = c("Male","Female")) +
  scale_linetype_manual(name = "Position:", values = c("L0" = "solid", "L15" = "dashed"),
                        labels = c("L0" = "0° Supine", "L15" = "15° HDT")) +
  theme(legend.position = "right")
# data.rpp
# 
# # Define Priors
# priors.rpp <- c(prior(normal(0, 1), class = b), # slope prior
#                 prior(normal(0, 1), class = b, coef = "Intercept"), # intercept prior
#                 prior(cauchy(0,1), class = sigma), # population variance
#                 prior(cauchy(0,1), class = sd) # group variance
# )
# 
# # # Fit Model
# # model.rpp <- brm(RPP.s ~ 0 + Intercept + Index + Sex.c + Face + (1|ID),
# #              data = df.s,
# #              prior = priors.rpp,
# #              sample_prior = "yes",
# #              family = gaussian(),
# #              warmup = 1000, # burn-in
# #              iter = 5000, # number of iterations
# #              chains = 4,  # number of MCMC chains
# #              control = list(adapt_delta = 0.95),# advanced MC settings
# #              save_pars=save_pars(all=TRUE))
# # saveRDS(model.rpp,file = "BRMSFits/model.rpp.rds")
# model.rpp <- readRDS("BRMSFits/model.rpp.rds")
# 
# # Posterior Checks
# summary(model.rpp)
# plot(model.rpp, N = 7)
# pp_check(model.rpp, ndraws = 50)
# model.rpp.loo <- loo(model.rpp, moment_match = TRUE, save_psis = TRUE)
# model.rpp.loo
# model.rpp.w <- weights(model.rpp.loo$psis_object)
# ppc_loo_pit_overlay(df.s$RPP.s[!is.na(df.s$RPP.s)], 
#                     posterior_predict(model.rpp), 
#                     lw = model.rpp.w)
# 
# # Describe Posterior
# describe_posterior(model.rpp,
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.05, 0.05), rope_ci = 1,
#                    diagnostic = NULL, effects = "fixed", component = "all")
# 
# # Rescale Posterior
# rescale.rpp <- model.rpp %>% spread_draws(`b_.*`, sigma, `sd_.*`, regex = TRUE) %>%
#   mutate(b_Intercept = (b_Intercept * scale.rpp) + center.rpp,
#          b_Index = b_Index * scale.rpp,
#          b_Sex.c = b_Sex.c * scale.rpp,
#          b_FaceL15 = b_FaceL15 * scale.rpp,
#          sigma = sigma * scale.rpp,
#          sd_ID__Intercept = sd_ID__Intercept * scale.rpp)
# describe_posterior(rescale.rpp[,-c(1:3)],
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.1*sd.rpp, 0.1*sd.rpp), rope_ci = 1,
#                    diagnostic = NULL)
# 
# # Plot Model
# model.rpp.dr <- fitted(model.rpp, newdata = nd, re_formula = NA, probs = c(.055,.945)) %>%
#   data.frame() %>%
#   bind_cols(nd) %>%
#   mutate(Estimate.r = (Estimate * scale.rpp) + center.rpp,
#          Q5.5.r = (Q5.5 * scale.rpp) + center.rpp,
#          Q94.5.r = (Q94.5 * scale.rpp) + center.rpp)
# 
# dose.rpp <- ggplot(model.rpp.dr, aes(x = Index, color = factor(Sex.c), fill = factor(Sex.c))) +
#   geom_smooth(aes(y = Estimate.r, ymin = Q5.5.r, ymax = Q94.5.r),
#               stat = "identity",
#               alpha = 0.2) +
#   facet_wrap(~ Face, labeller = labeller(Face = face.labs), scales = "free") +
#   theme_classic2() +
#   theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
#   scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq, expand = c(0,0)) +
#   scale_y_continuous(name = "Rate Pressure Product (mmHg/min)") +
#   scale_color_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                      labels = c("-0.5" = "Female", "0.5" = "Male"),
#                      breaks = c("0.5","-0.5")) +
#   scale_fill_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                     labels = c("-0.5" = "Female", "0.5" = "Male"),
#                     breaks = c("0.5","-0.5")) +
#   theme(legend.position = "bottom", panel.spacing = unit(1,"lines")) +
#   coord_cartesian(ylim = c(6000,11000), expand = FALSE)
# dose.rpp

#####################
### Myocardial Oxygen Supply:Demand Index (DPTI/SPTI)
#####################

# Find scales
scale.mo <- attributes(df.s$MO.s)$`scaled:scale`
center.mo <- attributes(df.s$MO.s)$`scaled:center`
sd.mo <- sd(df$DPTI.SPTI[!is.na(df$DPTI.SPTI)])

# Plot Data
data.mo <- ggplot(data = df, aes(x = Index, color = as.factor(Sex), linetype = Face, y = DPTI.SPTI)) +
  #geom_jitter(alpha = 0.2) +
  stat_summary(aes(group = interaction(Sex,Face)), fun=mean, geom="line", position = pd) +
  stat_summary(aes(group = interaction(Sex,Face)), fun.data = mean_se, geom="errorbar", width = 0.4, position = pd, linetype = "solid") +
  stat_summary(aes(group = interaction(Sex,Face)), fun=mean, geom="point", position = pd) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
  scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq) +
  scale_y_continuous(name = "Myocardial Oxygen Supply:Demand Index") +
  scale_color_manual(name = "Sex:", values = c("Female" = "hotpink", "Male" = "royalblue"),
                     breaks = c("Male","Female")) +
  scale_linetype_manual(name = "Position:", values = c("L0" = "solid", "L15" = "dashed"),
                        labels = c("L0" = "0° Supine", "L15" = "15° HDT")) +
  theme(legend.position = "right")
# data.mo
# 
# # Define Priors
# priors.mo <- c(prior(normal(0, 1), class = b), # slope prior
#                 prior(normal(0, 1), class = b, coef = "Intercept"), # intercept prior
#                 prior(cauchy(0,1), class = sigma), # population variance
#                 prior(cauchy(0,1), class = sd), # group variance
#                 prior(gamma(2, .1), class = nu) # degrees of freedom
# )
# 
# # # Fit Model
# # model.mo <- brm(MO.s ~ 0 + Intercept + Index + Sex.c + Face + (1|ID),
# #              data = df.s,
# #              prior = priors.mo,
# #              sample_prior = "yes",
# #              family = student(),
# #              warmup = 1000, # burn-in
# #              iter = 5000, # number of iterations
# #              chains = 4,  # number of MCMC chains
# #              control = list(adapt_delta = 0.95),# advanced MC settings
# #              save_pars=save_pars(all=TRUE))
# # saveRDS(model.mo,file = "BRMSFits/model.mo.rds")
# model.mo <- readRDS("BRMSFits/model.mo.rds")
# 
# # Posterior Checks
# summary(model.mo)
# plot(model.mo, N = 7)
# pp_check(model.mo, ndraws = 50)
# model.mo.loo <- loo(model.mo, moment_match = TRUE, save_psis = TRUE)
# model.mo.loo
# model.mo.w <- weights(model.mo.loo$psis_object)
# ppc_loo_pit_overlay(df.s$MO.s[!is.na(df.s$MO.s)], 
#                     posterior_predict(model.mo), 
#                     lw = model.mo.w)
# 
# # Describe Posterior
# describe_posterior(model.mo,
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.05, 0.05), rope_ci = 1,
#                    diagnostic = NULL, effects = "fixed", component = "all")
# 
# # Rescale Posterior
# rescale.mo <- model.mo %>% spread_draws(`b_.*`, sigma, `sd_.*`, nu, regex = TRUE) %>%
#   mutate(b_Intercept = (b_Intercept * scale.mo) + center.mo,
#          b_Index = b_Index * scale.mo,
#          b_Sex.c = b_Sex.c * scale.mo,
#          b_FaceL15 = b_FaceL15 * scale.mo,
#          sigma = sigma * scale.mo,
#          sd_ID__Intercept = sd_ID__Intercept * scale.mo)
# describe_posterior(rescale.mo[,-c(1:3)],
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.1*sd.mo, 0.1*sd.mo), rope_ci = 1,
#                    diagnostic = NULL)
# 
# # Plot Model
# model.mo.dr <- fitted(model.mo, newdata = nd, re_formula = NA, probs = c(.055,.945)) %>%
#   data.frame() %>%
#   bind_cols(nd) %>%
#   mutate(Estimate.r = (Estimate * scale.mo) + center.mo,
#          Q5.5.r = (Q5.5 * scale.mo) + center.mo,
#          Q94.5.r = (Q94.5 * scale.mo) + center.mo)
# 
# dose.mo <- ggplot(model.mo.dr, aes(x = Index, color = factor(Sex.c), fill = factor(Sex.c))) +
#   geom_smooth(aes(y = Estimate.r, ymin = Q5.5.r, ymax = Q94.5.r),
#               stat = "identity",
#               alpha = 0.2) +
#   facet_wrap(~ Face, labeller = labeller(Face = face.labs), scales = "free") +
#   theme_classic2() +
#   theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
#   scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq, expand = c(0,0)) +
#   scale_y_continuous(name = "Myocardial Oxygen Supply:Demand Index") +
#   scale_color_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                      labels = c("-0.5" = "Female", "0.5" = "Male"),
#                      breaks = c("0.5","-0.5")) +
#   scale_fill_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                     labels = c("-0.5" = "Female", "0.5" = "Male"),
#                     breaks = c("0.5","-0.5")) +
#   theme(legend.position = "bottom", panel.spacing = unit(1,"lines")) +
#   coord_cartesian(ylim = c(0.3,0.7), expand = FALSE)
# dose.mo

#####################
### Total Peripheral Resistance
#####################

# Find scales
scale.tpr <- attributes(df.s$TPR.s)$`scaled:scale`
center.tpr <- attributes(df.s$TPR.s)$`scaled:center`
sd.tpr <- sd(df.s$TPR.calc[!is.na(df.s$TPR.calc)])

# Plot Data
data.tpr <- ggplot(data = df.s, aes(x = Index, color = as.factor(Sex), linetype = Face, y = TPR.calc)) +
  #geom_jitter(alpha = 0.2) +
  stat_summary(aes(group = interaction(Sex,Face)), fun=mean, geom="line", position = pd) +
  stat_summary(aes(group = interaction(Sex,Face)), fun.data = mean_se, geom="errorbar", width = 0.4, position = pd, linetype = "solid") +
  stat_summary(aes(group = interaction(Sex,Face)), fun=mean, geom="point", position = pd) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
  scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq) +
  scale_y_continuous(name = "Total Peripheral Resistance (mmHg.s/ml)") +
  scale_color_manual(name = "Sex:", values = c("Female" = "hotpink", "Male" = "royalblue"),
                     breaks = c("Male","Female")) +
  scale_linetype_manual(name = "Position:", values = c("L0" = "solid", "L15" = "dashed"),
                        labels = c("L0" = "0° Supine", "L15" = "15° HDT")) +
  theme(legend.position = "right")
# data.tpr
# 
# # Define Priors
# priors.tpr <- c(prior(normal(0, 1), class = b), # slope prior
#                prior(normal(0, 1), class = b, coef = "Intercept"), # intercept prior
#                prior(cauchy(0,1), class = sigma), # population variance
#                prior(cauchy(0,1), class = sd), # group variance
#                prior(gamma(2, .1), class = nu) # degrees of freedom
# )
# 
# # # Fit Model
# # model.tpr <- brm(TPR.s ~ 0 + Intercept + Index + Sex.c + Face + (1|ID),
# #              data = df.s,
# #              prior = priors.tpr,
# #              sample_prior = "yes",
# #              family = student(),
# #              warmup = 1000, # burn-in
# #              iter = 5000, # number of iterations
# #              chains = 4,  # number of MCMC chains
# #              control = list(adapt_delta = 0.95),# advanced MC settings
# #              save_pars=save_pars(all=TRUE))
# # saveRDS(model.tpr,file = "BRMSFits/model.tpr.rds")
# model.tpr <- readRDS("BRMSFits/model.tpr.rds")
# 
# # Posterior Checks
# summary(model.tpr)
# plot(model.tpr, N = 7)
# pp_check(model.tpr, ndraws = 50)
# model.tpr.loo <- loo(model.tpr, moment_match = TRUE, save_psis = TRUE)
# model.tpr.loo
# model.tpr.w <- weights(model.tpr.loo$psis_object)
# ppc_loo_pit_overlay(df.s$TPR.s[!is.na(df.s$TPR.s)], 
#                     posterior_predict(model.tpr), 
#                     lw = model.tpr.w)
# 
# # Describe Posterior
# describe_posterior(model.tpr,
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.05, 0.05), rope_ci = 1,
#                    diagnostic = NULL, effects = "fixed", component = "all")
# 
# # Rescale Posterior
# rescale.tpr <- model.tpr %>% spread_draws(`b_.*`, sigma, `sd_.*`, nu, regex = TRUE) %>%
#   mutate(b_Intercept = (b_Intercept * scale.tpr) + center.tpr,
#          b_Index = b_Index * scale.tpr,
#          b_Sex.c = b_Sex.c * scale.tpr,
#          b_FaceL15 = b_FaceL15 * scale.tpr,
#          sigma = sigma * scale.tpr,
#          sd_ID__Intercept = sd_ID__Intercept * scale.tpr)
# describe_posterior(rescale.tpr[,-c(1:3)],
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.1*sd.tpr, 0.1*sd.tpr), rope_ci = 1,
#                    diagnostic = NULL)
# 
# # Plot Model
# model.tpr.dr <- fitted(model.tpr, newdata = nd, re_formula = NA, probs = c(.055,.945)) %>%
#   data.frame() %>%
#   bind_cols(nd) %>%
#   mutate(Estimate.r = (Estimate * scale.tpr) + center.tpr,
#          Q5.5.r = (Q5.5 * scale.tpr) + center.tpr,
#          Q94.5.r = (Q94.5 * scale.tpr) + center.tpr)
# 
# dose.tpr <- ggplot(model.tpr.dr, aes(x = Index, color = factor(Sex.c), fill = factor(Sex.c))) +
#   geom_smooth(aes(y = Estimate.r, ymin = Q5.5.r, ymax = Q94.5.r),
#               stat = "identity",
#               alpha = 0.2) +
#   facet_wrap(~ Face, labeller = labeller(Face = face.labs), scales = "free") +
#   theme_classic2() +
#   theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
#   scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq, expand = c(0,0)) +
#   scale_y_continuous(name = "Total Peripheral Resistance (mmHg.s/ml)") +
#   scale_color_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                      labels = c("-0.5" = "Female", "0.5" = "Male"),
#                      breaks = c("0.5","-0.5")) +
#   scale_fill_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                     labels = c("-0.5" = "Female", "0.5" = "Male"),
#                     breaks = c("0.5","-0.5")) +
#   theme(legend.position = "bottom", panel.spacing = unit(1,"lines")) +
#   coord_cartesian(ylim = c(0.7,2.1), expand = FALSE)
# dose.tpr

##########################################
### Systemic Indices
##########################################
#####################
### Cardiac Index
#####################

# Find scales
scale.ci <- attributes(df.s$CI.s)$`scaled:scale`
center.ci <- attributes(df.s$CI.s)$`scaled:center`
sd.ci <- sd(df.s$CI[!is.na(df.s$CI)])

# Plot Data
data.ci <- ggplot(data = df.s, aes(x = Index, color = as.factor(Sex), linetype = Face, y = CI)) +
  #geom_jitter(alpha = 0.2) +
  stat_summary(aes(group = interaction(Sex,Face)), fun=mean, geom="line", position = pd) +
  stat_summary(aes(group = interaction(Sex,Face)), fun.data = mean_se, geom="errorbar", width = 0.4, position = pd, linetype = "solid") +
  stat_summary(aes(group = interaction(Sex,Face)), fun=mean, geom="point", position = pd) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
  scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq) +
  scale_y_continuous(name = expression(paste("Cardiac Index (l/min/", m^{2},")"))) +
  scale_color_manual(name = "Sex:", values = c("Female" = "hotpink", "Male" = "royalblue"),
                     breaks = c("Male","Female")) +
  scale_linetype_manual(name = "Position:", values = c("L0" = "solid", "L15" = "dashed"),
                        labels = c("L0" = "0° Supine", "L15" = "15° HDT")) +
  theme(legend.position = "right")
# data.ci
# 
# # Define Priors
# priors.ci <- c(prior(normal(0, 1), class = b), # slope prior
#                 prior(normal(0, 1), class = b, coef = "Intercept"), # intercept prior
#                 prior(cauchy(0,1), class = sigma), # population variance
#                 prior(cauchy(0,1), class = sd), # group variance
#                 prior(gamma(2, .1), class = nu) # degrees of freedom
# )
# 
# # # Fit Model
# # model.ci <- brm(CI.s ~ 0 + Intercept + Index + Sex.c + Face + (1|ID),
# #              data = df.s,
# #              prior = priors.ci,
# #              sample_prior = "yes",
# #              family = student(),
# #              warmup = 1000, # burn-in
# #              iter = 5000, # number of iterations
# #              chains = 4,  # number of MCMC chains
# #              control = list(adapt_delta = 0.95),# advanced MC settings
# #              save_pars=save_pars(all=TRUE))
# # saveRDS(model.ci,file = "BRMSFits/model.ci.rds")
# model.ci <- readRDS("BRMSFits/model.ci.rds")
# 
# # Posterior Checks
# summary(model.ci)
# plot(model.ci, N = 7)
# pp_check(model.ci, ndraws = 50)
# model.ci.loo <- loo(model.ci, moment_match = TRUE, save_psis = TRUE)
# model.ci.loo
# model.ci.w <- weights(model.ci.loo$psis_object)
# ppc_loo_pit_overlay(df.s$CI.s[!is.na(df.s$CI.s)], 
#                     posterior_predict(model.ci), 
#                     lw = model.ci.w)
# 
# # Describe Posterior
# describe_posterior(model.ci,
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.05, 0.05), rope_ci = 1,
#                    diagnostic = NULL, effects = "fixed", component = "all")
# 
# # Rescale Posterior
# rescale.ci <- model.ci %>% spread_draws(`b_.*`, sigma, `sd_.*`, nu, regex = TRUE) %>%
#   mutate(b_Intercept = (b_Intercept * scale.ci) + center.ci,
#          b_Index = b_Index * scale.ci,
#          b_Sex.c = b_Sex.c * scale.ci,
#          b_FaceL15 = b_FaceL15 * scale.ci,
#          sigma = sigma * scale.ci,
#          sd_ID__Intercept = sd_ID__Intercept * scale.ci)
# describe_posterior(rescale.ci[,-c(1:3)],
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.1*sd.ci, 0.1*sd.ci), rope_ci = 1,
#                    diagnostic = NULL)
# 
# # Plot Model
# model.ci.dr <- fitted(model.ci, newdata = nd, re_formula = NA, probs = c(.055,.945)) %>%
#   data.frame() %>%
#   bind_cols(nd) %>%
#   mutate(Estimate.r = (Estimate * scale.ci) + center.ci,
#          Q5.5.r = (Q5.5 * scale.ci) + center.ci,
#          Q94.5.r = (Q94.5 * scale.ci) + center.ci)
# 
# dose.ci <- ggplot(model.ci.dr, aes(x = Index, color = factor(Sex.c), fill = factor(Sex.c))) +
#   geom_smooth(aes(y = Estimate.r, ymin = Q5.5.r, ymax = Q94.5.r),
#               stat = "identity",
#               alpha = 0.2) +
#   facet_wrap(~ Face, labeller = labeller(Face = face.labs), scales = "free") +
#   theme_classic2() +
#   theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
#   scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq, expand = c(0,0)) +
#   scale_y_continuous(name = expression(paste("Cardiac Index (l/min/", m^{2},")"))) +
#   scale_color_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                      labels = c("-0.5" = "Female", "0.5" = "Male"),
#                      breaks = c("0.5","-0.5")) +
#   scale_fill_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                     labels = c("-0.5" = "Female", "0.5" = "Male"),
#                     breaks = c("0.5","-0.5")) +
#   theme(legend.position = "bottom", panel.spacing = unit(1,"lines")) +
#   coord_cartesian(ylim = c(1,3.5), expand = FALSE)
# dose.ci

#####################
### Stroke Index
#####################

# Find scales
scale.si <- attributes(df.s$SI.s)$`scaled:scale`
center.si <- attributes(df.s$SI.s)$`scaled:center`
sd.si <- sd(df.s$SI[!is.na(df.s$SI)])

# Plot Data
data.si <- ggplot(data = df.s, aes(x = Index, color = as.factor(Sex), linetype = Face, y = SI)) +
  #geom_jitter(alpha = 0.2) +
  stat_summary(aes(group = interaction(Sex,Face)), fun=mean, geom="line", position = pd) +
  stat_summary(aes(group = interaction(Sex,Face)), fun.data = mean_se, geom="errorbar", width = 0.4, position = pd, linetype = "solid") +
  stat_summary(aes(group = interaction(Sex,Face)), fun=mean, geom="point", position = pd) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
  scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq) +
  scale_y_continuous(name = expression(paste("Stroke Index (ml/", m^{2},")"))) +
  scale_color_manual(name = "Sex:", values = c("Female" = "hotpink", "Male" = "royalblue"),
                     breaks = c("Male","Female")) +
  scale_linetype_manual(name = "Position:", values = c("L0" = "solid", "L15" = "dashed"),
                        labels = c("L0" = "0° Supine", "L15" = "15° HDT")) +
  theme(legend.position = "right")
# data.si
# 
# # Define Priors
# priors.si <- c(prior(normal(0, 1), class = b), # slope prior
#                prior(normal(0, 1), class = b, coef = "Intercept"), # intercept prior
#                prior(cauchy(0,1), class = sigma), # population variance
#                prior(cauchy(0,1), class = sd), # group variance
#                prior(gamma(2, .1), class = nu) # degrees of freedom
# )
# 
# # # Fit Model
# # model.si <- brm(SI.s ~ 0 + Intercept + Index + Sex.c + Face + (1|ID),
# #              data = df.s,
# #              prior = priors.si,
# #              sample_prior = "yes",
# #              family = student(),
# #              warmup = 1000, # burn-in
# #              iter = 5000, # number of iterations
# #              chains = 4,  # number of MCMC chains
# #              control = list(adapt_delta = 0.95),# advanced MC settings
# #              save_pars=save_pars(all=TRUE))
# # saveRDS(model.si,file = "BRMSFits/model.si.rds")
# model.si <- readRDS("BRMSFits/model.si.rds")
# 
# # Posterior Checks
# summary(model.si)
# plot(model.si, N = 7)
# pp_check(model.si, ndraws = 50)
# model.si.loo <- loo(model.si, moment_match = TRUE, save_psis = TRUE)
# model.si.loo
# model.si.w <- weights(model.si.loo$psis_object)
# ppc_loo_pit_overlay(df.s$SI.s[!is.na(df.s$SI.s)], 
#                     posterior_predict(model.si), 
#                     lw = model.si.w)
# 
# # Describe Posterior
# describe_posterior(model.si,
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.05, 0.05), rope_ci = 1,
#                    diagnostic = NULL, effects = "fixed", component = "all")
# 
# # Rescale Posterior
# rescale.si <- model.si %>% spread_draws(`b_.*`, sigma, `sd_.*`, nu, regex = TRUE) %>%
#   mutate(b_Intercept = (b_Intercept * scale.si) + center.si,
#          b_Index = b_Index * scale.si,
#          b_Sex.c = b_Sex.c * scale.si,
#          b_FaceL15 = b_FaceL15 * scale.si,
#          sigma = sigma * scale.si,
#          sd_ID__Intercept = sd_ID__Intercept * scale.si)
# describe_posterior(rescale.si[,-c(1:3)],
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.1*sd.si, 0.1*sd.si), rope_ci = 1,
#                    diagnostic = NULL)
# 
# # Plot Model
# model.si.dr <- fitted(model.si, newdata = nd, re_formula = NA, probs = c(.055,.945)) %>%
#   data.frame() %>%
#   bind_cols(nd) %>%
#   mutate(Estimate.r = (Estimate * scale.si) + center.si,
#          Q5.5.r = (Q5.5 * scale.si) + center.si,
#          Q94.5.r = (Q94.5 * scale.si) + center.si)
# 
# dose.si <- ggplot(model.si.dr, aes(x = Index, color = factor(Sex.c), fill = factor(Sex.c))) +
#   geom_smooth(aes(y = Estimate.r, ymin = Q5.5.r, ymax = Q94.5.r),
#               stat = "identity",
#               alpha = 0.2) +
#   facet_wrap(~ Face, labeller = labeller(Face = face.labs), scales = "free") +
#   theme_classic2() +
#   theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
#   scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq, expand = c(0,0)) +
#   scale_y_continuous(name = expression(paste("Stroke Index (ml/", m^{2},")"))) +
#   scale_color_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                      labels = c("-0.5" = "Female", "0.5" = "Male"),
#                      breaks = c("0.5","-0.5")) +
#   scale_fill_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                     labels = c("-0.5" = "Female", "0.5" = "Male"),
#                     breaks = c("0.5","-0.5")) +
#   theme(legend.position = "bottom", panel.spacing = unit(1,"lines")) +
#   coord_cartesian(ylim = c(10,60), expand = FALSE)
# dose.si

##########################################
### Autonomic Indices
##########################################
#####################
### SDNN
#####################

# Find scales
scale.sdnn <- attributes(df.s$SDNN.s)$`scaled:scale`
center.sdnn <- attributes(df.s$SDNN.s)$`scaled:center`
sd.sdnn <- sd(df.s$SDNN[!is.na(df.s$SDNN)])

# Plot Data
data.sdnn <- ggplot(data = df.s, aes(x = Index, color = as.factor(Sex), linetype = Face, y = SDNN)) +
  #geom_jitter(alpha = 0.2) +
  stat_summary(aes(group = interaction(Sex,Face)), fun=mean, geom="line", position = pd) +
  stat_summary(aes(group = interaction(Sex,Face)), fun.data = mean_se, geom="errorbar", width = 0.4, position = pd, linetype = "solid") +
  stat_summary(aes(group = interaction(Sex,Face)), fun=mean, geom="point", position = pd) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
  scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq) +
  scale_y_continuous(name = "SDNN (ms)") +
  scale_color_manual(name = "Sex:", values = c("Female" = "hotpink", "Male" = "royalblue"),
                     breaks = c("Male","Female")) +
  scale_linetype_manual(name = "Position:", values = c("L0" = "solid", "L15" = "dashed"),
                        labels = c("L0" = "0° Supine", "L15" = "15° HDT")) +
  theme(legend.position = "right")
# data.sdnn
# 
# # Define Priors
# priors.sdnn <- c(prior(normal(0, 1), class = b), # slope prior
#                prior(normal(0, 1), class = b, coef = "Intercept"), # intercept prior
#                prior(cauchy(0,1), class = sigma), # population variance
#                prior(cauchy(0,1), class = sd), # group variance
#                prior(gamma(2, .1), class = nu) # degrees of freedom
# )
# 
# # # Fit Model
# # model.sdnn <- brm(SDNN.s ~ 0 + Intercept + Index + Sex.c + Face + (1|ID),
# #              data = df.s,
# #              prior = priors.sdnn,
# #              sample_prior = "yes",
# #              family = student(),
# #              warmup = 1000, # burn-in
# #              iter = 5000, # number of iterations
# #              chains = 4,  # number of MCMC chains
# #              control = list(adapt_delta = 0.95),# advanced MC settings
# #              save_pars=save_pars(all=TRUE))
# # saveRDS(model.sdnn,file = "BRMSFits/model.sdnn.rds")
# model.sdnn <- readRDS("BRMSFits/model.sdnn.rds")
# 
# # Posterior Checks
# summary(model.sdnn)
# plot(model.sdnn, N = 7)
# pp_check(model.sdnn, ndraws = 50)
# model.sdnn.loo <- loo(model.sdnn, moment_match = TRUE, save_psis = TRUE)
# model.sdnn.loo
# model.sdnn.w <- weights(model.sdnn.loo$psis_object)
# ppc_loo_pit_overlay(df.s$SDNN.s[!is.na(df.s$SDNN.s)], 
#                     posterior_predict(model.sdnn), 
#                     lw = model.sdnn.w)
# 
# # Describe Posterior
# describe_posterior(model.sdnn,
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.05, 0.05), rope_ci = 1,
#                    diagnostic = NULL, effects = "fixed", component = "all")
# 
# # Rescale Posterior
# rescale.sdnn <- model.sdnn %>% spread_draws(`b_.*`, sigma, `sd_.*`, nu, regex = TRUE) %>%
#   mutate(b_Intercept = (b_Intercept * scale.sdnn) + center.sdnn,
#          b_Index = b_Index * scale.sdnn,
#          b_Sex.c = b_Sex.c * scale.sdnn,
#          b_FaceL15 = b_FaceL15 * scale.sdnn,
#          sigma = sigma * scale.sdnn,
#          sd_ID__Intercept = sd_ID__Intercept * scale.sdnn)
# describe_posterior(rescale.sdnn[,-c(1:3)],
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.1*sd.sdnn, 0.1*sd.sdnn), rope_ci = 1,
#                    diagnostic = NULL)
# 
# # Plot Model
# model.sdnn.dr <- fitted(model.sdnn, newdata = nd, re_formula = NA, probs = c(.055,.945)) %>%
#   data.frame() %>%
#   bind_cols(nd) %>%
#   mutate(Estimate.r = (Estimate * scale.sdnn) + center.sdnn,
#          Q5.5.r = (Q5.5 * scale.sdnn) + center.sdnn,
#          Q94.5.r = (Q94.5 * scale.sdnn) + center.sdnn)
# 
# dose.sdnn <- ggplot(model.sdnn.dr, aes(x = Index, color = factor(Sex.c), fill = factor(Sex.c))) +
#   geom_smooth(aes(y = Estimate.r, ymin = Q5.5.r, ymax = Q94.5.r),
#               stat = "identity",
#               alpha = 0.2) +
#   facet_wrap(~ Face, labeller = labeller(Face = face.labs), scales = "free") +
#   theme_classic2() +
#   theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
#   scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq, expand = c(0,0)) +
#   scale_y_continuous(name = "SDNN (ms)") +
#   scale_color_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                      labels = c("-0.5" = "Female", "0.5" = "Male"),
#                      breaks = c("0.5","-0.5")) +
#   scale_fill_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                     labels = c("-0.5" = "Female", "0.5" = "Male"),
#                     breaks = c("0.5","-0.5")) +
#   theme(legend.position = "bottom", panel.spacing = unit(1,"lines")) +
#   coord_cartesian(ylim = c(35,85), expand = FALSE)
# dose.sdnn

#####################
### HRVTI
#####################

# Find scales
scale.hrvti <- attributes(df.s$HRVTI.s)$`scaled:scale`
center.hrvti <- attributes(df.s$HRVTI.s)$`scaled:center`
sd.hrvti <- sd(df.s$HRVTI[!is.na(df.s$HRVTI)])

# Plot Data
data.hrvti <- ggplot(data = df.s, aes(x = Index, color = as.factor(Sex), linetype = Face, y = HRVTI)) +
  #geom_jitter(alpha = 0.2) +
  stat_summary(aes(group = interaction(Sex,Face)), fun=mean, geom="line", position = pd) +
  stat_summary(aes(group = interaction(Sex,Face)), fun.data = mean_se, geom="errorbar", width = 0.4, position = pd, linetype = "solid") +
  stat_summary(aes(group = interaction(Sex,Face)), fun=mean, geom="point", position = pd) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
  scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq) +
  scale_y_continuous(name = "HRVTi") +
  scale_color_manual(name = "Sex:", values = c("Female" = "hotpink", "Male" = "royalblue"),
                     breaks = c("Male","Female")) +
  scale_linetype_manual(name = "Position:", values = c("L0" = "solid", "L15" = "dashed"),
                        labels = c("L0" = "0° Supine", "L15" = "15° HDT")) +
  theme(legend.position = "right")
# data.hrvti
# 
# # Define Priors
# priors.hrvti <- c(prior(normal(0, 1), class = b), # slope prior
#                  prior(normal(0, 1), class = b, coef = "Intercept"), # intercept prior
#                  prior(cauchy(0,1), class = sigma), # population variance
#                  prior(cauchy(0,1), class = sd), # group variance
#                  prior(gamma(2, .1), class = nu) # degrees of freedom
# )
# 
# # # Fit Model
# # model.hrvti <- brm(HRVTI.s ~ 0 + Intercept + Index + Sex.c + Face + (1|ID),
# #              data = df.s,
# #              prior = priors.hrvti,
# #              sample_prior = "yes",
# #              family = student(),
# #              warmup = 1000, # burn-in
# #              iter = 5000, # number of iterations
# #              chains = 4,  # number of MCMC chains
# #              control = list(adapt_delta = 0.95),# advanced MC settings
# #              save_pars=save_pars(all=TRUE))
# # saveRDS(model.hrvti,file = "BRMSFits/model.hrvti.rds")
# model.hrvti <- readRDS("BRMSFits/model.hrvti.rds")
# 
# # Posterior Checks
# summary(model.hrvti)
# plot(model.hrvti, N = 7)
# pp_check(model.hrvti, ndraws = 50)
# model.hrvti.loo <- loo(model.hrvti, moment_match = TRUE, save_psis = TRUE)
# model.hrvti.loo
# model.hrvti.w <- weights(model.hrvti.loo$psis_object)
# ppc_loo_pit_overlay(df.s$HRVTI.s[!is.na(df.s$HRVTI.s)], 
#                     posterior_predict(model.hrvti), 
#                     lw = model.hrvti.w)
# 
# # Describe Posterior
# describe_posterior(model.hrvti,
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.05, 0.05), rope_ci = 1,
#                    diagnostic = NULL, effects = "fixed", component = "all")
# 
# # Rescale Posterior
# rescale.hrvti <- model.hrvti %>% spread_draws(`b_.*`, sigma, `sd_.*`, nu, regex = TRUE) %>%
#   mutate(b_Intercept = (b_Intercept * scale.hrvti) + center.hrvti,
#          b_Index = b_Index * scale.hrvti,
#          b_Sex.c = b_Sex.c * scale.hrvti,
#          b_FaceL15 = b_FaceL15 * scale.hrvti,
#          sigma = sigma * scale.hrvti,
#          sd_ID__Intercept = sd_ID__Intercept * scale.hrvti)
# describe_posterior(rescale.hrvti[,-c(1:3)],
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.1*sd.hrvti, 0.1*sd.hrvti), rope_ci = 1,
#                    diagnostic = NULL)
# 
# # Plot Model
# model.hrvti.dr <- fitted(model.hrvti, newdata = nd, re_formula = NA, probs = c(.055,.945)) %>%
#   data.frame() %>%
#   bind_cols(nd) %>%
#   mutate(Estimate.r = (Estimate * scale.hrvti) + center.hrvti,
#          Q5.5.r = (Q5.5 * scale.hrvti) + center.hrvti,
#          Q94.5.r = (Q94.5 * scale.hrvti) + center.hrvti)
# 
# dose.hrvti <- ggplot(model.hrvti.dr, aes(x = Index, color = factor(Sex.c), fill = factor(Sex.c))) +
#   geom_smooth(aes(y = Estimate.r, ymin = Q5.5.r, ymax = Q94.5.r),
#               stat = "identity",
#               alpha = 0.2) +
#   facet_wrap(~ Face, labeller = labeller(Face = face.labs), scales = "free") +
#   theme_classic2() +
#   theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
#   scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq, expand = c(0,0)) +
#   scale_y_continuous(name = "HRVTI") +
#   scale_color_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                      labels = c("-0.5" = "Female", "0.5" = "Male"),
#                      breaks = c("0.5","-0.5")) +
#   scale_fill_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                     labels = c("-0.5" = "Female", "0.5" = "Male"),
#                     breaks = c("0.5","-0.5")) +
#   theme(legend.position = "bottom", panel.spacing = unit(1,"lines")) +
#   coord_cartesian(ylim = c(9,16), expand = FALSE)
# dose.hrvti

#####################
### RMSDD
#####################

# Find scales
scale.rmsdd <- attributes(df.s$RMSDD.s)$`scaled:scale`
center.rmsdd <- attributes(df.s$RMSDD.s)$`scaled:center`
sd.rmsdd <- sd(df.s$RMSDD[!is.na(df.s$RMSDD)])

# Plot Data
data.rmsdd <- ggplot(data = df.s, aes(x = Index, color = as.factor(Sex), linetype = Face, y = RMSDD)) +
  #geom_jitter(alpha = 0.2) +
  stat_summary(aes(group = interaction(Sex,Face)), fun=mean, geom="line", position = pd) +
  stat_summary(aes(group = interaction(Sex,Face)), fun.data = mean_se, geom="errorbar", width = 0.4, position = pd, linetype = "solid") +
  stat_summary(aes(group = interaction(Sex,Face)), fun=mean, geom="point", position = pd) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
  scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq) +
  scale_y_continuous(name = "RMSDD (ms)") +
  scale_color_manual(name = "Sex:", values = c("Female" = "hotpink", "Male" = "royalblue"),
                     breaks = c("Male","Female")) +
  scale_linetype_manual(name = "Position:", values = c("L0" = "solid", "L15" = "dashed"),
                        labels = c("L0" = "0° Supine", "L15" = "15° HDT")) +
  theme(legend.position = "right")
# data.rmsdd
# 
# # Define Priors
# priors.rmsdd <- c(prior(normal(0, 1), class = b), # slope prior
#                   prior(normal(0, 1), class = b, coef = "Intercept"), # intercept prior
#                   prior(cauchy(0,1), class = sigma), # population variance
#                   prior(cauchy(0,1), class = sd), # group variance
#                   prior(gamma(2, .1), class = nu) # degrees of freedom
# )
# 
# # # Fit Model
# # model.rmsdd <- brm(RMSDD.s ~ 0 + Intercept + Index + Sex.c + Face + (1|ID),
# #              data = df.s,
# #              prior = priors.rmsdd,
# #              sample_prior = "yes",
# #              family = student(),
# #              warmup = 1000, # burn-in
# #              iter = 5000, # number of iterations
# #              chains = 4,  # number of MCMC chains
# #              control = list(adapt_delta = 0.95),# advanced MC settings
# #              save_pars=save_pars(all=TRUE))
# # saveRDS(model.rmsdd,file = "BRMSFits/model.rmsdd.rds")
# model.rmsdd <- readRDS("BRMSFits/model.rmsdd.rds")
# 
# # Posterior Checks
# summary(model.rmsdd)
# plot(model.rmsdd, N = 7)
# pp_check(model.rmsdd, ndraws = 50)
# model.rmsdd.loo <- loo(model.rmsdd, moment_match = TRUE, save_psis = TRUE)
# model.rmsdd.loo
# model.rmsdd.w <- weights(model.rmsdd.loo$psis_object)
# ppc_loo_pit_overlay(df.s$RMSDD.s[!is.na(df.s$RMSDD.s)], 
#                     posterior_predict(model.rmsdd), 
#                     lw = model.rmsdd.w)
# 
# # Describe Posterior
# describe_posterior(model.rmsdd,
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.05, 0.05), rope_ci = 1,
#                    diagnostic = NULL, effects = "fixed", component = "all")
# 
# # Rescale Posterior
# rescale.rmsdd <- model.rmsdd %>% spread_draws(`b_.*`, sigma, `sd_.*`, nu, regex = TRUE) %>%
#   mutate(b_Intercept = (b_Intercept * scale.rmsdd) + center.rmsdd,
#          b_Index = b_Index * scale.rmsdd,
#          b_Sex.c = b_Sex.c * scale.rmsdd,
#          b_FaceL15 = b_FaceL15 * scale.rmsdd,
#          sigma = sigma * scale.rmsdd,
#          sd_ID__Intercept = sd_ID__Intercept * scale.rmsdd)
# describe_posterior(rescale.rmsdd[,-c(1:3)],
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.1*sd.rmsdd, 0.1*sd.rmsdd), rope_ci = 1,
#                    diagnostic = NULL)
# 
# # Plot Model
# model.rmsdd.dr <- fitted(model.rmsdd, newdata = nd, re_formula = NA, probs = c(.055,.945)) %>%
#   data.frame() %>%
#   bind_cols(nd) %>%
#   mutate(Estimate.r = (Estimate * scale.rmsdd) + center.rmsdd,
#          Q5.5.r = (Q5.5 * scale.rmsdd) + center.rmsdd,
#          Q94.5.r = (Q94.5 * scale.rmsdd) + center.rmsdd)
# 
# dose.rmsdd <- ggplot(model.rmsdd.dr, aes(x = Index, color = factor(Sex.c), fill = factor(Sex.c))) +
#   geom_smooth(aes(y = Estimate.r, ymin = Q5.5.r, ymax = Q94.5.r),
#               stat = "identity",
#               alpha = 0.2) +
#   facet_wrap(~ Face, labeller = labeller(Face = face.labs), scales = "free") +
#   theme_classic2() +
#   theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
#   scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq, expand = c(0,0)) +
#   scale_y_continuous(name = "RMSDD (ms)") +
#   scale_color_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                      labels = c("-0.5" = "Female", "0.5" = "Male"),
#                      breaks = c("0.5","-0.5")) +
#   scale_fill_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                     labels = c("-0.5" = "Female", "0.5" = "Male"),
#                     breaks = c("0.5","-0.5")) +
#   theme(legend.position = "bottom", panel.spacing = unit(1,"lines")) +
#   coord_cartesian(ylim = c(5,55), expand = FALSE)
# dose.rmsdd

#####################
### Baroreflex Sensitivity
#####################

# Find scales
scale.brs <- attributes(df.s$BRS.s)$`scaled:scale`
center.brs <- attributes(df.s$BRS.s)$`scaled:center`
sd.brs <- sd(df.s$BRS[!is.na(df.s$BRS)])

# Plot Data
data.brs <- ggplot(data = df.s, aes(x = Index, color = as.factor(Sex), linetype = Face, y = BRS)) +
  #geom_jitter(alpha = 0.2) +
  stat_summary(aes(group = interaction(Sex,Face)), fun=mean, geom="line", position = pd) +
  stat_summary(aes(group = interaction(Sex,Face)), fun.data = mean_se, geom="errorbar", width = 0.4, position = pd, linetype = "solid") +
  stat_summary(aes(group = interaction(Sex,Face)), fun=mean, geom="point", position = pd) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
  scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq) +
  scale_y_continuous(name = "Baroreflex Sensitivity (ms/mmHg)") +
  scale_color_manual(name = "Sex:", values = c("Female" = "hotpink", "Male" = "royalblue"),
                     breaks = c("Male","Female")) +
  scale_linetype_manual(name = "Position:", values = c("L0" = "solid", "L15" = "dashed"),
                        labels = c("L0" = "0° Supine", "L15" = "15° HDT")) +
  theme(legend.position = "right")
# data.brs
# 
# # Define Priors
# priors.brs <- c(prior(normal(0, 1), class = b), # slope prior
#                   prior(normal(0, 1), class = b, coef = "Intercept"), # intercept prior
#                   prior(cauchy(0,1), class = sigma), # population variance
#                   prior(cauchy(0,1), class = sd), # group variance
#                   prior(gamma(2, .1), class = nu) # degrees of freedom
# )
# 
# # # Fit Model
# # model.brs <- brm(BRS.s ~ 0 + Intercept + Index + Sex.c + Face + (1|ID),
# #              data = df.s,
# #              prior = priors.brs,
# #              sample_prior = "yes",
# #              family = student(),
# #              warmup = 1000, # burn-in
# #              iter = 5000, # number of iterations
# #              chains = 4,  # number of MCMC chains
# #              control = list(adapt_delta = 0.95),# advanced MC settings
# #              save_pars=save_pars(all=TRUE))
# # saveRDS(model.brs,file = "BRMSFits/model.brs.rds")
# model.brs <- readRDS("BRMSFits/model.brs.rds")
# 
# # Posterior Checks
# summary(model.brs)
# plot(model.brs, N = 7)
# pp_check(model.brs, ndraws = 50)
# model.brs.loo <- loo(model.brs, moment_match = TRUE, save_psis = TRUE)
# model.brs.loo
# model.brs.w <- weights(model.brs.loo$psis_object)
# ppc_loo_pit_overlay(df.s$BRS.s[!is.na(df.s$BRS.s)], 
#                     posterior_predict(model.brs), 
#                     lw = model.brs.w)
# 
# # Describe Posterior
# describe_posterior(model.brs,
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.05, 0.05), rope_ci = 1,
#                    diagnostic = NULL, effects = "fixed", component = "all")
# 
# # Rescale Posterior
# rescale.brs <- model.brs %>% spread_draws(`b_.*`, sigma, `sd_.*`, nu, regex = TRUE) %>%
#   mutate(b_Intercept = (b_Intercept * scale.brs) + center.brs,
#          b_Index = b_Index * scale.brs,
#          b_Sex.c = b_Sex.c * scale.brs,
#          b_FaceL15 = b_FaceL15 * scale.brs,
#          sigma = sigma * scale.brs,
#          sd_ID__Intercept = sd_ID__Intercept * scale.brs)
# describe_posterior(rescale.brs[,-c(1:3)],
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.1*sd.brs, 0.1*sd.brs), rope_ci = 1,
#                    diagnostic = NULL)
# 
# # Plot Model
# model.brs.dr <- fitted(model.brs, newdata = nd, re_formula = NA, probs = c(.055,.945)) %>%
#   data.frame() %>%
#   bind_cols(nd) %>%
#   mutate(Estimate.r = (Estimate * scale.brs) + center.brs,
#          Q5.5.r = (Q5.5 * scale.brs) + center.brs,
#          Q94.5.r = (Q94.5 * scale.brs) + center.brs)
# 
# dose.brs <- ggplot(model.brs.dr, aes(x = Index, color = factor(Sex.c), fill = factor(Sex.c))) +
#   geom_smooth(aes(y = Estimate.r, ymin = Q5.5.r, ymax = Q94.5.r),
#               stat = "identity",
#               alpha = 0.2) +
#   facet_wrap(~ Face, labeller = labeller(Face = face.labs), scales = "free") +
#   theme_classic2() +
#   theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
#   scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq, expand = c(0,0)) +
#   scale_y_continuous(name = "Baroreflex Sensitivity (ms/mmHg)") +
#   scale_color_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                      labels = c("-0.5" = "Female", "0.5" = "Male"),
#                      breaks = c("0.5","-0.5")) +
#   scale_fill_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                     labels = c("-0.5" = "Female", "0.5" = "Male"),
#                     breaks = c("0.5","-0.5")) +
#   theme(legend.position = "bottom", panel.spacing = unit(1,"lines")) +
#   coord_cartesian(ylim = c(5,20), expand = FALSE)
# dose.brs

#####################
### LF Norm
#####################

# Find scales
scale.lfn <- attributes(df.s$LFNorm.s)$`scaled:scale`
center.lfn <- attributes(df.s$LFNorm.s)$`scaled:center`
sd.lfn <- sd(df.s$LFNorm[!is.na(df.s$LFNorm)])

# Plot Data
data.lfn <- ggplot(data = df.s, aes(x = Index, color = as.factor(Sex), linetype = Face, y = LFNorm)) +
  #geom_jitter(alpha = 0.2) +
  stat_summary(aes(group = interaction(Sex,Face)), fun=mean, geom="line", position = pd) +
  stat_summary(aes(group = interaction(Sex,Face)), fun.data = mean_se, geom="errorbar", width = 0.4, position = pd, linetype = "solid") +
  stat_summary(aes(group = interaction(Sex,Face)), fun=mean, geom="point", position = pd) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
  scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq) +
  scale_y_continuous(name = "Low Frequency (Normalized)") +
  scale_color_manual(name = "Sex:", values = c("Female" = "hotpink", "Male" = "royalblue"),
                     breaks = c("Male","Female")) +
  scale_linetype_manual(name = "Position:", values = c("L0" = "solid", "L15" = "dashed"),
                        labels = c("L0" = "0° Supine", "L15" = "15° HDT")) +
  theme(legend.position = "right")
# data.lfn
# 
# # Define Priors
# priors.lfn <- c(prior(normal(0, 1), class = b), # slope prior
#                 prior(normal(0, 1), class = b, coef = "Intercept"), # intercept prior
#                 prior(cauchy(0,1), class = sigma), # population variance
#                 prior(cauchy(0,1), class = sd), # group variance
#                 prior(gamma(2, .1), class = nu) # degrees of freedom
# )
# 
# # # Fit Model
# # model.lfn <- brm(LFNorm.s ~ 0 + Intercept + Index + Sex.c + Face + (1|ID),
# #              data = df.s,
# #              prior = priors.lfn,
# #              sample_prior = "yes",
# #              family = student(),
# #              warmup = 1000, # burn-in
# #              iter = 5000, # number of iterations
# #              chains = 4,  # number of MCMC chains
# #              control = list(adapt_delta = 0.95),# advanced MC settings
# #              save_pars=save_pars(all=TRUE))
# # saveRDS(model.lfn,file = "BRMSFits/model.lfn.rds")
# model.lfn <- readRDS("BRMSFits/model.lfn.rds")
# 
# # Posterior Checks
# summary(model.lfn)
# plot(model.lfn, N = 7)
# pp_check(model.lfn, ndraws = 50)
# model.lfn.loo <- loo(model.lfn, moment_match = TRUE, save_psis = TRUE)
# model.lfn.loo
# model.lfn.w <- weights(model.lfn.loo$psis_object)
# ppc_loo_pit_overlay(df.s$LFNorm.s[!is.na(df.s$LFNorm.s)], 
#                     posterior_predict(model.lfn), 
#                     lw = model.lfn.w)
# 
# # Describe Posterior
# describe_posterior(model.lfn,
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.05, 0.05), rope_ci = 1,
#                    diagnostic = NULL, effects = "fixed", component = "all")
# 
# # Rescale Posterior
# rescale.lfn <- model.lfn %>% spread_draws(`b_.*`, sigma, `sd_.*`, nu, regex = TRUE) %>%
#   mutate(b_Intercept = (b_Intercept * scale.lfn) + center.lfn,
#          b_Index = b_Index * scale.lfn,
#          b_Sex.c = b_Sex.c * scale.lfn,
#          b_FaceL15 = b_FaceL15 * scale.lfn,
#          sigma = sigma * scale.lfn,
#          sd_ID__Intercept = sd_ID__Intercept * scale.lfn)
# describe_posterior(rescale.lfn[,-c(1:3)],
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.1*sd.lfn, 0.1*sd.lfn), rope_ci = 1,
#                    diagnostic = NULL)
# 
# # Plot Model
# model.lfn.dr <- fitted(model.lfn, newdata = nd, re_formula = NA, probs = c(.055,.945)) %>%
#   data.frame() %>%
#   bind_cols(nd) %>%
#   mutate(Estimate.r = (Estimate * scale.lfn) + center.lfn,
#          Q5.5.r = (Q5.5 * scale.lfn) + center.lfn,
#          Q94.5.r = (Q94.5 * scale.lfn) + center.lfn)
# 
# dose.lfn <- ggplot(model.lfn.dr, aes(x = Index, color = factor(Sex.c), fill = factor(Sex.c))) +
#   geom_smooth(aes(y = Estimate.r, ymin = Q5.5.r, ymax = Q94.5.r),
#               stat = "identity",
#               alpha = 0.2) +
#   facet_wrap(~ Face, labeller = labeller(Face = face.labs), scales = "free") +
#   theme_classic2() +
#   theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
#   scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq, expand = c(0,0)) +
#   scale_y_continuous(name = "Low Frequency (Normalized)") +
#   scale_color_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                      labels = c("-0.5" = "Female", "0.5" = "Male"),
#                      breaks = c("0.5","-0.5")) +
#   scale_fill_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                     labels = c("-0.5" = "Female", "0.5" = "Male"),
#                     breaks = c("0.5","-0.5")) +
#   theme(legend.position = "bottom", panel.spacing = unit(1,"lines")) +
#   coord_cartesian(ylim = c(50,100), expand = FALSE)
# dose.lfn

#####################
### HF Norm
#####################

# Find scales
scale.hfn <- attributes(df.s$HFNorm.s)$`scaled:scale`
center.hfn <- attributes(df.s$HFNorm.s)$`scaled:center`
sd.hfn <- sd(df.s$HFNorm[!is.na(df.s$HFNorm)])

# Plot Data
data.hfn <- ggplot(data = df.s, aes(x = Index, color = as.factor(Sex), linetype = Face, y = HFNorm)) +
  #geom_jitter(alpha = 0.2) +
  stat_summary(aes(group = interaction(Sex,Face)), fun=mean, geom="line", position = pd) +
  stat_summary(aes(group = interaction(Sex,Face)), fun.data = mean_se, geom="errorbar", width = 0.4, position = pd, linetype = "solid") +
  stat_summary(aes(group = interaction(Sex,Face)), fun=mean, geom="point", position = pd) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
  scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq) +
  scale_y_continuous(name = "High Frequency (Normalized)") +
  scale_color_manual(name = "Sex:", values = c("Female" = "hotpink", "Male" = "royalblue"),
                     breaks = c("Male","Female")) +
  scale_linetype_manual(name = "Position:", values = c("L0" = "solid", "L15" = "dashed"),
                        labels = c("L0" = "0° Supine", "L15" = "15° HDT")) +
  theme(legend.position = "right")
# data.hfn
# 
# # Define Priors
# priors.hfn <- c(prior(normal(0, 1), class = b), # slope prior
#                 prior(normal(0, 1), class = b, coef = "Intercept"), # intercept prior
#                 prior(cauchy(0,1), class = sigma), # population variance
#                 prior(cauchy(0,1), class = sd), # group variance
#                 prior(gamma(2, .1), class = nu) # degrees of freedom
# )
# 
# # # Fit Model
# # model.hfn <- brm(HFNorm.s ~ 0 + Intercept + Index + Sex.c + Face + (1|ID),
# #              data = df.s,
# #              prior = priors.hfn,
# #              sample_prior = "yes",
# #              family = student(),
# #              warmup = 1000, # burn-in
# #              iter = 5000, # number of iterations
# #              chains = 4,  # number of MCMC chains
# #              control = list(adapt_delta = 0.95),# advanced MC settings
# #              save_pars=save_pars(all=TRUE))
# # saveRDS(model.hfn,file = "BRMSFits/model.hfn.rds")
# model.hfn <- readRDS("BRMSFits/model.hfn.rds")
# 
# # Posterior Checks
# summary(model.hfn)
# plot(model.hfn, N = 7)
# pp_check(model.hfn, ndraws = 50)
# model.hfn.loo <- loo(model.hfn, moment_match = TRUE, save_psis = TRUE)
# model.hfn.loo
# model.hfn.w <- weights(model.hfn.loo$psis_object)
# ppc_loo_pit_overlay(df.s$HFNorm.s[!is.na(df.s$HFNorm.s)], 
#                     posterior_predict(model.hfn), 
#                     lw = model.hfn.w)
# 
# # Describe Posterior
# describe_posterior(model.hfn,
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.05, 0.05), rope_ci = 1,
#                    diagnostic = NULL, effects = "fixed", component = "all")
# 
# # Rescale Posterior
# rescale.hfn <- model.hfn %>% spread_draws(`b_.*`, sigma, `sd_.*`, nu, regex = TRUE) %>%
#   mutate(b_Intercept = (b_Intercept * scale.hfn) + center.hfn,
#          b_Index = b_Index * scale.hfn,
#          b_Sex.c = b_Sex.c * scale.hfn,
#          b_FaceL15 = b_FaceL15 * scale.hfn,
#          sigma = sigma * scale.hfn,
#          sd_ID__Intercept = sd_ID__Intercept * scale.hfn)
# describe_posterior(rescale.hfn[,-c(1:3)],
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.1*sd.hfn, 0.1*sd.hfn), rope_ci = 1,
#                    diagnostic = NULL)
# 
# # Plot Model
# model.hfn.dr <- fitted(model.hfn, newdata = nd, re_formula = NA, probs = c(.055,.945)) %>%
#   data.frame() %>%
#   bind_cols(nd) %>%
#   mutate(Estimate.r = (Estimate * scale.hfn) + center.hfn,
#          Q5.5.r = (Q5.5 * scale.hfn) + center.hfn,
#          Q94.5.r = (Q94.5 * scale.hfn) + center.hfn)
# 
# dose.hfn <- ggplot(model.hfn.dr, aes(x = Index, color = factor(Sex.c), fill = factor(Sex.c))) +
#   geom_smooth(aes(y = Estimate.r, ymin = Q5.5.r, ymax = Q94.5.r),
#               stat = "identity",
#               alpha = 0.2) +
#   facet_wrap(~ Face, labeller = labeller(Face = face.labs), scales = "free") +
#   theme_classic2() +
#   theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
#   scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq, expand = c(0,0)) +
#   scale_y_continuous(name = "High Frequency (Normalized)") +
#   scale_color_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                      labels = c("-0.5" = "Female", "0.5" = "Male"),
#                      breaks = c("0.5","-0.5")) +
#   scale_fill_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                     labels = c("-0.5" = "Female", "0.5" = "Male"),
#                     breaks = c("0.5","-0.5")) +
#   theme(legend.position = "bottom", panel.spacing = unit(1,"lines")) +
#   coord_cartesian(ylim = c(0,50), expand = FALSE)
# dose.hfn

#####################
### LF/HF
#####################

# Find scales
scale.lfhf <- attributes(df.s$LFHF.s)$`scaled:scale`
center.lfhf <- attributes(df.s$LFHF.s)$`scaled:center`
sd.lfhf <- sd(df.s$LF.HF[!is.na(df.s$LF.HF)])

# Plot Data
data.lfhf <- ggplot(data = df.s, aes(x = Index, color = as.factor(Sex), linetype = Face, y = LF.HF)) +
  #geom_jitter(alpha = 0.2) +
  stat_summary(aes(group = interaction(Sex,Face)), fun=mean, geom="line", position = pd) +
  stat_summary(aes(group = interaction(Sex,Face)), fun.data = mean_se, geom="errorbar", width = 0.4, position = pd, linetype = "solid") +
  stat_summary(aes(group = interaction(Sex,Face)), fun=mean, geom="point", position = pd) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
  scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq) +
  scale_y_continuous(name = "Low Frequency/High Frequency Ratio") +
  scale_color_manual(name = "Sex:", values = c("Female" = "hotpink", "Male" = "royalblue"),
                     breaks = c("Male","Female")) +
  scale_linetype_manual(name = "Position:", values = c("L0" = "solid", "L15" = "dashed"),
                        labels = c("L0" = "0° Supine", "L15" = "15° HDT")) +
  theme(legend.position = "right")
# data.lfhf
# 
# # Define Priors
# priors.lfhf <- c(prior(normal(0, 1), class = b), # slope prior
#                 prior(normal(0, 1), class = b, coef = "Intercept"), # intercept prior
#                 prior(cauchy(0,1), class = sigma), # population variance
#                 prior(cauchy(0,1), class = sd), # group variance
#                 prior(gamma(2, .1), class = nu) # degrees of freedom
# )
# 
# # # Fit Model
# # model.lfhf <- brm(LFHF.s ~ 0 + Intercept + Index + Sex.c + Face + 
# #                     Index:Sex.c + (1|ID),
# #              data = df.s,
# #              prior = priors.lfhf,
# #              sample_prior = "yes",
# #              family = student(),
# #              warmup = 1000, # burn-in
# #              iter = 5000, # number of iterations
# #              chains = 4,  # number of MCMC chains
# #              control = list(adapt_delta = 0.95),# advanced MC settings
# #              save_pars=save_pars(all=TRUE))
# # saveRDS(model.lfhf,file = "BRMSFits/model.lfhf.rds")
# model.lfhf <- readRDS("BRMSFits/model.lfhf.rds")
# 
# # Posterior Checks
# summary(model.lfhf)
# plot(model.lfhf, N = 8)
# pp_check(model.lfhf, ndraws = 50)
# model.lfhf.loo <- loo(model.lfhf, moment_match = TRUE, save_psis = TRUE)
# model.lfhf.loo
# model.lfhf.w <- weights(model.lfhf.loo$psis_object)
# ppc_loo_pit_overlay(df.s$LFHF.s[!is.na(df.s$LFHF.s)], 
#                     posterior_predict(model.lfhf), 
#                     lw = model.lfhf.w)
# 
# # Describe Posterior
# describe_posterior(model.lfhf,
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.05, 0.05), rope_ci = 1,
#                    diagnostic = NULL, effects = "fixed", component = "all")
# 
# # Rescale Posterior
# rescale.lfhf <- model.lfhf %>% spread_draws(`b_.*`, sigma, `sd_.*`, nu, regex = TRUE) %>%
#   mutate(b_Intercept = (b_Intercept * scale.lfhf) + center.lfhf,
#          b_Index = b_Index * scale.lfhf,
#          b_Sex.c = b_Sex.c * scale.lfhf,
#          b_FaceL15 = b_FaceL15 * scale.lfhf,
#          `b_Index:Sex.c` = `b_Index:Sex.c` * scale.lfhf,
#          sigma = sigma * scale.lfhf,
#          sd_ID__Intercept = sd_ID__Intercept * scale.lfhf)
# describe_posterior(rescale.lfhf[,-c(1:3)],
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.1*sd.lfhf, 0.1*sd.lfhf), rope_ci = 1,
#                    diagnostic = NULL)
# 
# # Plot Model
# model.lfhf.dr <- fitted(model.lfhf, newdata = nd, re_formula = NA, probs = c(.055,.945)) %>%
#   data.frame() %>%
#   bind_cols(nd) %>%
#   mutate(Estimate.r = (Estimate * scale.lfhf) + center.lfhf,
#          Q5.5.r = (Q5.5 * scale.lfhf) + center.lfhf,
#          Q94.5.r = (Q94.5 * scale.lfhf) + center.lfhf)
# 
# dose.lfhf <- ggplot(model.lfhf.dr, aes(x = Index, color = factor(Sex.c), fill = factor(Sex.c))) +
#   geom_smooth(aes(y = Estimate.r, ymin = Q5.5.r, ymax = Q94.5.r),
#               stat = "identity",
#               alpha = 0.2) +
#   facet_wrap(~ Face, labeller = labeller(Face = face.labs), scales = "free") +
#   theme_classic2() +
#   theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
#   scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq, expand = c(0,0)) +
#   scale_y_continuous(name = "Low Frequency/High Frequency Ratio") +
#   scale_color_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                      labels = c("-0.5" = "Female", "0.5" = "Male"),
#                      breaks = c("0.5","-0.5")) +
#   scale_fill_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                     labels = c("-0.5" = "Female", "0.5" = "Male"),
#                     breaks = c("0.5","-0.5")) +
#   theme(legend.position = "bottom", panel.spacing = unit(1,"lines")) +
#   coord_cartesian(ylim = c(0,10), expand = FALSE)
# dose.lfhf

##########################################
### Head/Neck
##########################################
#####################
### Intraocular Pressure
#####################

# Find scales
scale.iop <- attributes(df.s$IOP.s)$`scaled:scale`
center.iop <- attributes(df.s$IOP.s)$`scaled:center`
sd.iop <- sd(df.s$IOP[!is.na(df.s$IOP)])

# Plot Data
data.iop <- ggplot(data = df.s, aes(x = Index, color = as.factor(Sex), linetype = Face, y = IOP, size = as.factor(Side))) +
  #geom_jitter(alpha = 0.2) +
  stat_summary(aes(group = interaction(Sex,Face,Side)), fun=mean, geom="line", position = pd) +
  stat_summary(aes(group = interaction(Sex,Face,Side)), fun.data = mean_se, geom="errorbar", width = 1, position = pd, linetype = "solid") +
  stat_summary(aes(group = interaction(Sex,Face,Side)), fun=mean, geom="point", position = pd) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
  scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq) +
  scale_y_continuous(name = "Intraocular Pressure (mmHg)") +
  scale_color_manual(name = "Sex:", values = c("Female" = "hotpink", "Male" = "royalblue"),
                     breaks = c("Male","Female")) +
  scale_linetype_manual(name = "Position:", values = c("L0" = "solid", "L15" = "dashed"),
                        labels = c("L0" = "0° Supine", "L15" = "15° HDT")) +
  scale_size_manual("Side:", values = c("0.5" = 0.75, "-0.5" = 0.25), labels = c("0.5" = "Right", "-0.5" = "Left"),
                    breaks = c("0.5","-0.5")) +
  theme(legend.position = "right")
# data.iop
# 
# # Define Priors
# priors.iop <- c(prior(normal(0, 1), class = b), # slope prior
#                 prior(normal(0, 1), class = b, coef = "Intercept"), # intercept prior
#                 prior(cauchy(0,1), class = sigma), # population variance
#                 prior(cauchy(0,1), class = sd), # group variance
#                 prior(gamma(2, .1), class = nu) # degrees of freedom
# )
# 
# # # Fit Model
# # model.iop <- brm(IOP.s ~ 0 + Intercept + Index + Sex.c + Face + Side + (1|ID),
# #              data = df.s,
# #              prior = priors.iop,
# #              sample_prior = "yes",
# #              family = student(),
# #              warmup = 1000, # burn-in
# #              iter = 5000, # number of iterations
# #              chains = 4,  # number of MCMC chains
# #              control = list(adapt_delta = 0.95),# advanced MC settings
# #              save_pars=save_pars(all=TRUE))
# # saveRDS(model.iop,file = "BRMSFits/model.iop.rds")
# model.iop <- readRDS("BRMSFits/model.iop.rds")
# 
# # Posterior Checks
# summary(model.iop)
# plot(model.iop, N = 8)
# pp_check(model.iop, ndraws = 50)
# model.iop.loo <- loo(model.iop, moment_match = TRUE, save_psis = TRUE)
# model.iop.loo
# model.iop.w <- weights(model.iop.loo$psis_object)
# ppc_loo_pit_overlay(df.s$IOP.s[!is.na(df.s$IOP.s)], 
#                     posterior_predict(model.iop), 
#                     lw = model.iop.w)
# 
# # Describe Posterior
# describe_posterior(model.iop,
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.05, 0.05), rope_ci = 1,
#                    diagnostic = NULL, effects = "fixed", component = "all")
# 
# # Rescale Posterior
# rescale.iop <- model.iop %>% spread_draws(`b_.*`, sigma, `sd_.*`, nu, regex = TRUE) %>%
#   mutate(b_Intercept = (b_Intercept * scale.iop) + center.iop,
#          b_Index = b_Index * scale.iop,
#          b_Sex.c = b_Sex.c * scale.iop,
#          b_FaceL15 = b_FaceL15 * scale.iop,
#          b_Side = b_Side * scale.iop,
#          sigma = sigma * scale.iop,
#          sd_ID__Intercept = sd_ID__Intercept * scale.iop)
# describe_posterior(rescale.iop[,-c(1:3)],
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.1*sd.iop, 0.1*sd.iop), rope_ci = 1,
#                    diagnostic = NULL)
# 
# # Plot Model
# model.iop.dr <- fitted(model.iop, newdata = nd2, re_formula = NA, probs = c(.055,.945)) %>%
#   data.frame() %>%
#   bind_cols(nd2) %>%
#   mutate(Estimate.r = (Estimate * scale.iop) + center.iop,
#          Q5.5.r = (Q5.5 * scale.iop) + center.iop,
#          Q94.5.r = (Q94.5 * scale.iop) + center.iop)
# 
# dose.iop <- ggplot(model.iop.dr, aes(x = Index, color = factor(Sex.c), fill = factor(Sex.c))) +
#   geom_smooth(aes(y = Estimate.r, ymin = Q5.5.r, ymax = Q94.5.r),
#               stat = "identity",
#               alpha = 0.2) +
#   facet_grid(Side ~ Face, labeller = labeller(Face = face.labs, Side = side.labs), scales = "free") +
#   theme_classic2() +
#   theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
#   scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq, expand = c(0,0)) +
#   scale_y_continuous(name = "Intraocular Pressure (mmHg)") +
#   scale_color_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                      labels = c("-0.5" = "Female", "0.5" = "Male"),
#                      breaks = c("0.5","-0.5")) +
#   scale_fill_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                     labels = c("-0.5" = "Female", "0.5" = "Male"),
#                     breaks = c("0.5","-0.5")) +
#   theme(legend.position = "bottom", panel.spacing = unit(1,"lines")) +
#   coord_cartesian(ylim = c(10,25), expand = FALSE)
# dose.iop

#####################
### IJVP (Close)
#####################

# Find scales
scale.vpc <- attributes(df.s$VPC.s)$`scaled:scale`
center.vpc <- attributes(df.s$VPC.s)$`scaled:center`
sd.vpc <- sd(df.s$VPC[!is.na(df.s$VPC)])

# Plot Data
data.vpc <- ggplot(data = df.s, aes(x = Index, color = as.factor(Sex), linetype = Face, y = VPC, size = as.factor(Side))) +
  #geom_jitter(alpha = 0.2) +
  stat_summary(aes(group = interaction(Sex,Face,Side)), fun=mean, geom="line", position = pd) +
  stat_summary(aes(group = interaction(Sex,Face,Side)), fun.data = mean_se, geom="errorbar", width = 1, position = pd, linetype = "solid") +
  stat_summary(aes(group = interaction(Sex,Face,Side)), fun=mean, geom="point", position = pd) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
  scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq) +
  scale_y_continuous(name = "IJV Pressure (mmHg)") +
  scale_color_manual(name = "Sex:", values = c("Female" = "hotpink", "Male" = "royalblue"),
                     breaks = c("Male","Female")) +
  scale_linetype_manual(name = "Position:", values = c("L0" = "solid", "L15" = "dashed"),
                        labels = c("L0" = "0° Supine", "L15" = "15° HDT")) +
  scale_size_manual("Side:", values = c("0.5" = 0.75, "-0.5" = 0.25), labels = c("0.5" = "Right", "-0.5" = "Left"),
                    breaks = c("0.5","-0.5")) +
  theme(legend.position = "right")
# data.vpc
# 
# # Define Priors
# priors.vpc <- c(prior(normal(0, 1), class = b), # slope prior
#                 prior(normal(0, 1), class = b, coef = "Intercept"), # intercept prior
#                 prior(normal(0,1), class = b, dpar = sigma), # population variance slope
#                 #prior(normal(0,1), class = Intercept, dpar = sigma), # population variance intercept
#                 prior(cauchy(0,1), class = sd), # group variance
#                 prior(gamma(2, .1), class = nu) # degrees of freedom
# )
# 
# # # Fit Model
# # model.vpc <- brm(bf(VPC.s ~ 0 + Intercept + Index + Sex.c + Face + Side + (1|ID),
# #                     sigma ~ 0 + Intercept + Index),
# #                  data = df.s,
# #                  prior = priors.vpc,
# #                  sample_prior = "yes",
# #                  family = student(),
# #                  warmup = 1000, # burn-in
# #                  iter = 5000, # number of iterations
# #                  chains = 4,  # number of MCMC chains
# #                  control = list(adapt_delta = 0.95),# advanced MC settings
# #                  save_pars=save_pars(all=TRUE))
# # saveRDS(model.vpc,file = "BRMSFits/model.vpc.rds")
# model.vpc <- readRDS("BRMSFits/model.vpc.rds")
# 
# # Posterior Checks
# summary(model.vpc)
# plot(model.vpc, N = 9)
# pp_check(model.vpc, ndraws = 50)
# model.vpc.loo <- loo(model.vpc, moment_match = TRUE, save_psis = TRUE)
# model.vpc.loo
# model.vpc.w <- weights(model.vpc.loo$psis_object)
# ppc_loo_pit_overlay(df.s$VPC.s[!is.na(df.s$VPC.s)], 
#                     posterior_predict(model.vpc), 
#                     lw = model.vpc.w)
# 
# # Describe Posterior
# describe_posterior(model.vpc,
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.05, 0.05), rope_ci = 1,
#                    diagnostic = NULL, effects = "fixed", component = "all")
# 
# get_variables(model.vpc)
# 
# # Rescale Posterior
# rescale.vpc <- model.vpc %>% spread_draws(`b_.*`, `sd_.*`, nu, regex = TRUE) %>%
#   mutate(b_Intercept = (b_Intercept * scale.vpc) + center.vpc,
#          b_Index = b_Index * scale.vpc,
#          b_Sex.c = b_Sex.c * scale.vpc,
#          b_FaceL15 = b_FaceL15 * scale.vpc,
#          b_Side = b_Side * scale.vpc,
#          b_sigma_Intercept = exp(b_sigma_Intercept) * scale.vpc,
#          b_sigma_Index = exp(b_sigma_Index) * scale.vpc,
#          sd_ID__Intercept = sd_ID__Intercept * scale.vpc)
# describe_posterior(rescale.vpc[,-c(1:3)],
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.1*sd.vpc, 0.1*sd.vpc), rope_ci = 1,
#                    diagnostic = NULL)
# 
# # Plot Model
# model.vpc.dr <- fitted(model.vpc, newdata = nd2, re_formula = NA, probs = c(.055,.945)) %>%
#   data.frame() %>%
#   bind_cols(nd2) %>%
#   mutate(Estimate.r = (Estimate * scale.vpc) + center.vpc,
#          Q5.5.r = (Q5.5 * scale.vpc) + center.vpc,
#          Q94.5.r = (Q94.5 * scale.vpc) + center.vpc)
# 
# dose.vpc <- ggplot(model.vpc.dr, aes(x = Index, color = factor(Sex.c), fill = factor(Sex.c))) +
#   geom_smooth(aes(y = Estimate.r, ymin = Q5.5.r, ymax = Q94.5.r),
#               stat = "identity",
#               alpha = 0.2) +
#   facet_grid(Side ~ Face, labeller = labeller(Face = face.labs, Side = side.labs), scales = "free") +
#   theme_classic2() +
#   theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
#   scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq, expand = c(0,0)) +
#   scale_y_continuous(name = "IJV Pressure (mmHg)") +
#   scale_color_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                      labels = c("-0.5" = "Female", "0.5" = "Male"),
#                      breaks = c("0.5","-0.5")) +
#   scale_fill_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                     labels = c("-0.5" = "Female", "0.5" = "Male"),
#                     breaks = c("0.5","-0.5")) +
#   theme(legend.position = "bottom", panel.spacing = unit(1,"lines")) +
#   coord_cartesian(ylim = c(0,30), expand = FALSE)
# dose.vpc

#####################
### IJVF
#####################

# Plot Data
data.ijvf <- ggplot(data=subset(df.s, !is.na(IJVF)), aes(x = Index, color = as.factor(Sex), y = IJVF)) +
  geom_jitter(height = 0) +
  #stat_summary(aes(group = interaction(Sex,Face,Side)), fun=mean, geom="line", position = pd) +
  #stat_summary(aes(group = interaction(Sex,Face,Side)), fun.data = mean_se, geom="errorbar", width = 1, position = pd, linetype = "solid") +
  #stat_summary(aes(group = interaction(Sex,Face,Side)), fun=mean, geom="point", position = pd) +
  facet_grid(Side ~ Face, labeller = labeller(Face = face.labs, Side = side.labs), scales = "free") +
  theme_classic2() +
  theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
  scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq) +
  scale_y_discrete(name = "IJV Flow Pattern", breaks = c(1,2,3)) +
  scale_color_manual(name = "Sex:", values = c("Female" = "hotpink", "Male" = "royalblue"),
                     breaks = c("Male","Female")) +
  scale_linetype_manual(name = "Position:", values = c("L0" = "solid", "L15" = "dashed"),
                        labels = c("L0" = "0° Supine", "L15" = "15° HDT")) +
  scale_size_manual("Side:", values = c("0.5" = 0.75, "-0.5" = 0.25), labels = c("0.5" = "Right", "-0.5" = "Left"),
                    breaks = c("0.5","-0.5")) +
  theme(legend.position = "right") +
  coord_cartesian(ylim = c(1,3))
data.ijvf

# # Link Function Selection
# cumulativemodelfit <- function(formula, data, links=c("logit", "probit", "cloglog", "cauchit"),
#                                thresholds=c("flexible", "equidistant"), verbose=TRUE) {
#   names(links) <- links
#   names(thresholds) <- thresholds
#   llks <- outer(links, thresholds,
#                 Vectorize(function(link, threshold)
#                   # catch error for responses with 2 levels
#                   tryCatch(ordinal::clm(formula, data=data, link=link, threshold=threshold)$logLik,
#                            error = function(e) NA)))
#   print(llks)
#   if (verbose) {
#     bestfit <- which.max(llks)
#     cat("\nThe best link function is ", links[bestfit %% length(links)], " with a ",
#         thresholds[1 + bestfit %/% length(thresholds)], " threshold (logLik ", llks[bestfit],
#         ")\n", sep="")
#   }
#   invisible(llks)
# }
# 
# cumulativemodelfit(IJVF ~ 1, data=df.s)
# 
# # Define Priors
# priors.ijvf <- c(prior(normal(0, 1), class = b), # slope prior
#                 prior(normal(0, 10), class = Intercept), # intercept prior
#                 prior(cauchy(0,1), class = sd) # group variance
# )
# 
# # # Fit Model
# # model.ijvf <- brm(IJVF ~ Index + Sex.c + Face + Side + (1|ID),
# #                   data = df.s,
# #                   sample_prior = "yes",
# #                   family = cumulative(link = "logit", threshold = "flexible"),
# #                   prior = priors.ijvf,
# #                   warmup = 1000, # burn-in
# #                   iter = 5000, # number of iterations
# #                   chains = 4,  # number of MCMC chains
# #                   control = list(adapt_delta = 0.95),# advanced MC settings
# #                   save_pars=save_pars(all=TRUE))
# # saveRDS(model.ijvf,file = "BRMSFits/model.ijvf.rds")
model.ijvf <- readRDS("BRMSFits/model.ijvf.rds")
# 
# # Posterior Checks
# summary(model.ijvf)
# plot(model.ijvf, N = 9)
# pp_check(model.ijvf, ndraws = 50)
# model.ijvf.loo <- loo(model.ijvf, moment_match = TRUE, save_psis = TRUE)
# model.ijvf.loo
# 
# # Describe Posterior
# describe_posterior(model.ijvf,
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.1*pi/sqrt(3), 0.1*pi/sqrt(3)), rope_ci = 1,
#                    diagnostic = NULL, effects = "fixed", component = "all")
# 
# get_variables(model.ijvf)
# 
# 
# library(MASS)
# m <- polr(IJVF ~ Index*Sex.c*Face*Side, data = df.s, Hess=TRUE)
# summary(m)
# 
# 
# # Plot Model v1
# model.ijvf.dr <- cbind(nd2, fitted(model.ijvf, nd2, re_formula = NA, probs = c(.055,.945))) %>%
#   pivot_longer(cols = !Index:Face,
#                names_to = c(".value", "Pattern"),
#                names_sep = "P")
# 
# dose.ijvf <- ggplot(filter(model.ijvf.dr, Pattern == "(Y = 1)"), aes(x = Index, color = factor(Sex.c), fill = factor(Sex.c))) +
#   geom_smooth(aes(y = Estimate., ymin = Q5.5., ymax = Q94.5.),
#               stat = "identity",
#               alpha = 0.2) +
#   facet_grid(Side ~ Face, labeller = labeller(Face = face.labs, Side = side.labs), scales = "free") +
#   theme_classic2() +
#   theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
#   scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq, expand = c(0,0)) +
#   scale_y_continuous(name = "P(Grade 1 Flow)") +
#   scale_color_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                      labels = c("-0.5" = "Female", "0.5" = "Male"),
#                      breaks = c("0.5","-0.5")) +
#   scale_fill_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                     labels = c("-0.5" = "Female", "0.5" = "Male"),
#                     breaks = c("0.5","-0.5")) +
#   theme(legend.position = "bottom", panel.spacing = unit(1,"lines")) +
#   coord_cartesian(ylim = c(0,1), expand = FALSE)
# dose.ijvf
# 
# # Plot model v2
# model.ijvf.dr2 <- cbind(nd2, fitted(model.ijvf, nd2, re_formula = NA, probs = c(.055,.945))) %>%
#   mutate(estimate = 1 - `Estimate.P(Y = 1)`,
#          conf.high = 1 - `Q5.5.P(Y = 1)`,
#          conf.low = 1 - `Q94.5.P(Y = 1)`)
# 
# dose.ijvf2 <- ggplot(model.ijvf.dr2, aes(x = Index, color = factor(Sex.c), fill = factor(Sex.c))) +
#   geom_smooth(aes(y = estimate, ymin = conf.low, ymax = conf.high),
#               stat = "identity",
#               alpha = 0.2) +
#   facet_grid(Side ~ Face, labeller = labeller(Face = face.labs, Side = side.labs), scales = "free") +
#   theme_classic2() +
#   theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
#   scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq, expand = c(0,0)) +
#   scale_y_continuous(name = "1 - P(Grade 1 Flow)") +
#   scale_color_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                      labels = c("-0.5" = "Female", "0.5" = "Male"),
#                      breaks = c("0.5","-0.5")) +
#   scale_fill_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                     labels = c("-0.5" = "Female", "0.5" = "Male"),
#                     breaks = c("0.5","-0.5")) +
#   theme(legend.position = "right", panel.spacing = unit(1,"lines")) +
#   coord_cartesian(ylim = c(0,1), expand = FALSE)
# dose.ijvf2
# ggsave("FlowPattern.png", dose.ijvf2, width = 7, height = 4, units = "in")
# 
# # Plot model v3
# dose.ijvf3 <- ggplot(model.ijvf.dr, aes(x = Index, color = factor(Pattern), fill = factor(Pattern))) +
#   #geom_line(aes(y = Prob)) +
#   geom_smooth(aes(y = Estimate., ymin = Q5.5., ymax = Q94.5.),
#               stat = "identity",
#               alpha = 0.2) +
#   facet_grid(Side ~ interaction(Sex.c,Face), labeller = labeller(Side = side.labs, .cols = int.labs), scales = "free") +
#   theme_classic2() +
#   theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
#   scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq, expand = c(0,0)) +
#   scale_y_continuous(name = "Probability of Flow Pattern") +
#   scale_fill_brewer("Flow Pattern", palette = "Dark2") +
#   scale_color_brewer("Flow Pattern", palette = "Dark2") +
#   #scale_color_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#   #                   labels = c("-0.5" = "Female", "0.5" = "Male"),
#   #                   breaks = c("0.5","-0.5")) +
#   #scale_fill_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#   #                  labels = c("-0.5" = "Female", "0.5" = "Male"),
#   #                  breaks = c("0.5","-0.5")) +
#   theme(legend.position = "bottom", panel.spacing = unit(1,"lines")) +
#   coord_cartesian(ylim = c(0,1), expand = FALSE)
# dose.ijvf3

#####################
### A_CCA
#####################

# Find scales
scale.cca <- attributes(df.s$CCA.s)$`scaled:scale`
center.cca <- attributes(df.s$CCA.s)$`scaled:center`
sd.cca <- sd(df.s$CCA[!is.na(df.s$CCA)])

# Plot Data
data.cca <- ggplot(data = df.s, aes(x = Index, color = as.factor(Sex), linetype = Face, y = CCA, size = as.factor(Side))) +
  #geom_jitter(alpha = 0.2) +
  stat_summary(aes(group = interaction(Sex,Face,Side)), fun=mean, geom="line", position = pd) +
  stat_summary(aes(group = interaction(Sex,Face,Side)), fun.data = mean_se, geom="errorbar", width = 1, position = pd, linetype = "solid") +
  stat_summary(aes(group = interaction(Sex,Face,Side)), fun=mean, geom="point", position = pd) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
  scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq) +
  scale_y_continuous(name = expression(paste("CCA Area (",mm^2,")",sep=""))) +
  scale_color_manual(name = "Sex:", values = c("Female" = "hotpink", "Male" = "royalblue"),
                     breaks = c("Male","Female")) +
  scale_linetype_manual(name = "Position:", values = c("L0" = "solid", "L15" = "dashed"),
                        labels = c("L0" = "0° Supine", "L15" = "15° HDT")) +
  scale_size_manual("Side:", values = c("0.5" = 0.75, "-0.5" = 0.25), labels = c("0.5" = "Right", "-0.5" = "Left"),
                    breaks = c("0.5","-0.5")) +
  theme(legend.position = "right")
# data.cca
# 
# # Define Priors
# priors.cca <- c(prior(normal(0, 1), class = b), # slope prior
#                 prior(normal(0, 1), class = b, coef = "Intercept"), # intercept prior
#                 prior(cauchy(0,1), class = sigma), # population variance
#                 prior(cauchy(0,1), class = sd), # group variance
#                 prior(gamma(2, .1), class = nu) # degrees of freedom
# )
# 
# # # Fit Model
# # model.cca <- brm(CCA.s ~ 0 + Intercept + Index + Sex.c + Face + Side + (1|ID),
# #                  data = df.s,
# #                  prior = priors.cca,
# #                  sample_prior = "yes",
# #                  family = student(),
# #                  warmup = 1000, # burn-in
# #                  iter = 5000, # number of iterations
# #                  chains = 4,  # number of MCMC chains
# #                  control = list(adapt_delta = 0.95),# advanced MC settings
# #                  save_pars=save_pars(all=TRUE))
# # saveRDS(model.cca,file = "BRMSFits/model.cca.rds")
# model.cca <- readRDS("BRMSFits/model.cca.rds")
# 
# # Posterior Checks
# summary(model.cca)
# plot(model.cca, N = 9)
# pp_check(model.cca, ndraws = 50)
# model.cca.loo <- loo(model.cca, moment_match = TRUE, save_psis = TRUE)
# model.cca.loo
# model.cca.w <- weights(model.cca.loo$psis_object)
# ppc_loo_pit_overlay(df.s$CCA.s[!is.na(df.s$CCA.s)], 
#                     posterior_predict(model.cca), 
#                     lw = model.cca.w)
# 
# # Describe Posterior
# describe_posterior(model.cca,
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.05, 0.05), rope_ci = 1,
#                    diagnostic = NULL, effects = "fixed", component = "all")
# 
# get_variables(model.cca)
# 
# # Rescale Posterior
# rescale.cca <- model.cca %>% spread_draws(`b_.*`, `sigma`, `sd_.*`, nu, regex = TRUE) %>%
#   mutate(b_Intercept = (b_Intercept * scale.cca) + center.cca,
#          b_Index = b_Index * scale.cca,
#          b_Sex.c = b_Sex.c * scale.cca,
#          b_FaceL15 = b_FaceL15 * scale.cca,
#          b_Side = b_Side * scale.cca,
#          sigma = sigma * scale.cca,
#          sd_ID__Intercept = sd_ID__Intercept * scale.cca)
# describe_posterior(rescale.cca[,-c(1:3)],
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.1*sd.cca, 0.1*sd.cca), rope_ci = 1,
#                    diagnostic = NULL)
# 
# # Plot Model
# model.cca.dr <- fitted(model.cca, newdata = nd2, re_formula = NA, probs = c(.055,.945)) %>%
#   data.frame() %>%
#   bind_cols(nd2) %>%
#   mutate(Estimate.r = (Estimate * scale.cca) + center.cca,
#          Q5.5.r = (Q5.5 * scale.cca) + center.cca,
#          Q94.5.r = (Q94.5 * scale.cca) + center.cca)
# 
# dose.cca <- ggplot(model.cca.dr, aes(x = Index, color = factor(Sex.c), fill = factor(Sex.c))) +
#   geom_smooth(aes(y = Estimate.r, ymin = Q5.5.r, ymax = Q94.5.r),
#               stat = "identity",
#               alpha = 0.2) +
#   facet_grid(Side ~ Face, labeller = labeller(Face = face.labs, Side = side.labs), scales = "free") +
#   theme_classic2() +
#   theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
#   scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq, expand = c(0,0)) +
#   scale_y_continuous(name = expression(paste("CCA Area (",mm^2,")",sep=""))) +
#   scale_color_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                      labels = c("-0.5" = "Female", "0.5" = "Male"),
#                      breaks = c("0.5","-0.5")) +
#   scale_fill_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                     labels = c("-0.5" = "Female", "0.5" = "Male"),
#                     breaks = c("0.5","-0.5")) +
#   theme(legend.position = "bottom", panel.spacing = unit(1,"lines")) +
#   coord_cartesian(ylim = c(20,40), expand = FALSE)
# dose.cca

#####################
### A_IJV
#####################

# Find scales
scale.ijv <- attributes(df.s$IJV.s)$`scaled:scale`
center.ijv <- attributes(df.s$IJV.s)$`scaled:center`
sd.ijv <- sd(df.s$IJV[!is.na(df.s$IJV)])

# Plot Data
data.ijv <- ggplot(data = df.s, aes(x = Index, color = as.factor(Sex), linetype = Face, y = IJV, size = as.factor(Side))) +
  #geom_jitter(alpha = 0.2) +
  stat_summary(aes(group = interaction(Sex,Face,Side)), fun=mean, geom="line", position = pd) +
  stat_summary(aes(group = interaction(Sex,Face,Side)), fun.data = mean_se, geom="errorbar", width = 1, position = pd, linetype = "solid") +
  stat_summary(aes(group = interaction(Sex,Face,Side)), fun=mean, geom="point", position = pd) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
  scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq) +
  scale_y_continuous(name = expression(paste("IJV Area (",mm^2,")",sep=""))) +
  scale_color_manual(name = "Sex:", values = c("Female" = "hotpink", "Male" = "royalblue"),
                     breaks = c("Male","Female")) +
  scale_linetype_manual(name = "Position:", values = c("L0" = "solid", "L15" = "dashed"),
                        labels = c("L0" = "0° Supine", "L15" = "15° HDT")) +
  scale_size_manual("Side:", values = c("0.5" = 0.75, "-0.5" = 0.25), labels = c("0.5" = "Right", "-0.5" = "Left"),
                    breaks = c("0.5","-0.5")) +
  theme(legend.position = "right")
# data.ijv
# 
# # Define Priors
# priors.ijv <- c(prior(normal(0, 1), class = b), # slope prior
#                 prior(normal(0, 1), class = b, coef = "Intercept"), # intercept prior
#                 prior(cauchy(0,1), class = sigma), # population variance
#                 prior(cauchy(0,1), class = sd), # group variance
#                 prior(gamma(2, .1), class = nu) # degrees of freedom
# )
# 
# # # Fit Model
# # model.ijv <- brm(IJV.s ~ 0 + Intercept + Index + Sex.c + Face + Side + (1|ID),
# #                  data = df.s,
# #                  prior = priors.ijv,
# #                  sample_prior = "yes",
# #                  family = student(),
# #                  warmup = 1000, # burn-in
# #                  iter = 5000, # number of iterations
# #                  chains = 4,  # number of MCMC chains
# #                  control = list(adapt_delta = 0.95),# advanced MC settings
# #                  save_pars=save_pars(all=TRUE))
# # saveRDS(model.ijv,file = "BRMSFits/model.ijv.rds")
# model.ijv <- readRDS("BRMSFits/model.ijv.rds")
# 
# # Posterior Checks
# summary(model.ijv)
# plot(model.ijv, N = 9)
# pp_check(model.ijv, ndraws = 50)
# model.ijv.loo <- loo(model.ijv, moment_match = TRUE, save_psis = TRUE)
# model.ijv.loo
# model.ijv.w <- weights(model.ijv.loo$psis_object)
# ppc_loo_pit_overlay(df.s$IJV.s[!is.na(df.s$IJV.s)], 
#                     posterior_predict(model.ijv), 
#                     lw = model.ijv.w)
# 
# # Describe Posterior
# describe_posterior(model.ijv,
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.05, 0.05), rope_ci = 1,
#                    diagnostic = NULL, effects = "fixed", component = "all")
# 
# get_variables(model.ijv)
# 
# # Rescale Posterior
# rescale.ijv <- model.ijv %>% spread_draws(`b_.*`, `sigma`, `sd_.*`, nu, regex = TRUE) %>%
#   mutate(b_Intercept = (b_Intercept * scale.ijv) + center.ijv,
#          b_Index = b_Index * scale.ijv,
#          b_Sex.c = b_Sex.c * scale.ijv,
#          b_FaceL15 = b_FaceL15 * scale.ijv,
#          b_Side = b_Side * scale.ijv,
#          sigma = sigma * scale.ijv,
#          sd_ID__Intercept = sd_ID__Intercept * scale.ijv)
# describe_posterior(rescale.ijv[,-c(1:3)],
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.1*sd.ijv, 0.1*sd.ijv), rope_ci = 1,
#                    diagnostic = NULL)
# 
# # Plot Model
# model.ijv.dr <- fitted(model.ijv, newdata = nd2, re_formula = NA, probs = c(.055,.945)) %>%
#   data.frame() %>%
#   bind_cols(nd2) %>%
#   mutate(Estimate.r = (Estimate * scale.ijv) + center.ijv,
#          Q5.5.r = (Q5.5 * scale.ijv) + center.ijv,
#          Q94.5.r = (Q94.5 * scale.ijv) + center.ijv)
# 
# dose.ijv <- ggplot(model.ijv.dr, aes(x = Index, color = factor(Sex.c), fill = factor(Sex.c))) +
#   geom_smooth(aes(y = Estimate.r, ymin = Q5.5.r, ymax = Q94.5.r),
#               stat = "identity",
#               alpha = 0.2) +
#   facet_grid(Side ~ Face, labeller = labeller(Face = face.labs, Side = side.labs), scales = "free") +
#   theme_classic2() +
#   theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
#   scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq, expand = c(0,0)) +
#   scale_y_continuous(name = expression(paste("IJV Area (",mm^2,")",sep=""))) +
#   scale_color_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                      labels = c("-0.5" = "Female", "0.5" = "Male"),
#                      breaks = c("0.5","-0.5")) +
#   scale_fill_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                     labels = c("-0.5" = "Female", "0.5" = "Male"),
#                     breaks = c("0.5","-0.5")) +
#   theme(legend.position = "bottom", panel.spacing = unit(1,"lines")) +
#   coord_cartesian(ylim = c(0,120), expand = FALSE)
# dose.ijv

#####################
### CCA PSV
#####################

# Find scales
scale.psv <- attributes(df.s$PSV.s)$`scaled:scale`
center.psv <- attributes(df.s$PSV.s)$`scaled:center`
sd.psv <- sd(df.s$PSV[!is.na(df.s$PSV)])

# Plot Data
data.psv <- ggplot(data = df.s, aes(x = Index, color = as.factor(Sex), linetype = Face, y = PSV, size = as.factor(Side))) +
  #geom_jitter(alpha = 0.2) +
  stat_summary(aes(group = interaction(Sex,Face,Side)), fun=mean, geom="line", position = pd) +
  stat_summary(aes(group = interaction(Sex,Face,Side)), fun.data = mean_se, geom="errorbar", width = 1, position = pd, linetype = "solid") +
  stat_summary(aes(group = interaction(Sex,Face,Side)), fun=mean, geom="point", position = pd) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
  scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq) +
  scale_y_continuous(name = "Peak Systolic Velocity (cm/s)") +
  scale_color_manual(name = "Sex:", values = c("Female" = "hotpink", "Male" = "royalblue"),
                     breaks = c("Male","Female")) +
  scale_linetype_manual(name = "Position:", values = c("L0" = "solid", "L15" = "dashed"),
                        labels = c("L0" = "0° Supine", "L15" = "15° HDT")) +
  scale_size_manual("Side:", values = c("0.5" = 0.75, "-0.5" = 0.25), labels = c("0.5" = "Right", "-0.5" = "Left"),
                    breaks = c("0.5","-0.5")) +
  theme(legend.position = "right")
# data.psv
# 
# # Define Priors
# priors.psv <- c(prior(normal(0, 1), class = b), # slope prior
#                 prior(normal(0, 1), class = b, coef = "Intercept"), # intercept prior
#                 prior(cauchy(0,1), class = sigma), # population variance
#                 prior(cauchy(0,1), class = sd), # group variance
#                 prior(gamma(2, .1), class = nu) # degrees of freedom
# )
# 
# # # Fit Model
# # model.psv <- brm(PSV.s ~ 0 + Intercept + Index + Sex.c + Face + Side + (1|ID),
# #              data = df.s,
# #              prior = priors.psv,
# #              sample_prior = "yes",
# #              family = student(),
# #              warmup = 1000, # burn-in
# #              iter = 5000, # number of iterations
# #              chains = 4,  # number of MCMC chains
# #              control = list(adapt_delta = 0.95),# advanced MC settings
# #              save_pars=save_pars(all=TRUE))
# # saveRDS(model.psv,file = "BRMSFits/model.psv.rds")
# model.psv <- readRDS("BRMSFits/model.psv.rds")
# 
# # Posterior Checks
# summary(model.psv)
# plot(model.psv, N = 8)
# pp_check(model.psv, ndraws = 50)
# model.psv.loo <- loo(model.psv, moment_match = TRUE, save_psis = TRUE)
# model.psv.loo
# model.psv.w <- weights(model.psv.loo$psis_object)
# ppc_loo_pit_overlay(df.s$PSV.s[!is.na(df.s$PSV.s)], 
#                     posterior_predict(model.psv), 
#                     lw = model.psv.w)
# 
# # Describe Posterior
# describe_posterior(model.psv,
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.05, 0.05), rope_ci = 1,
#                    diagnostic = NULL, effects = "fixed", component = "all")
# 
# # Rescale Posterior
# rescale.psv <- model.psv %>% spread_draws(`b_.*`, sigma, `sd_.*`, nu, regex = TRUE) %>%
#   mutate(b_Intercept = (b_Intercept * scale.psv) + center.psv,
#          b_Index = b_Index * scale.psv,
#          b_Sex.c = b_Sex.c * scale.psv,
#          b_FaceL15 = b_FaceL15 * scale.psv,
#          b_Side = b_Side * scale.psv,
#          sigma = sigma * scale.psv,
#          sd_ID__Intercept = sd_ID__Intercept * scale.psv)
# describe_posterior(rescale.psv[,-c(1:3)],
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.1*sd.psv, 0.1*sd.psv), rope_ci = 1,
#                    diagnostic = NULL)
# 
# # Plot Model
# model.psv.dr <- fitted(model.psv, newdata = nd2, re_formula = NA, probs = c(.055,.945)) %>%
#   data.frame() %>%
#   bind_cols(nd2) %>%
#   mutate(Estimate.r = (Estimate * scale.psv) + center.psv,
#          Q5.5.r = (Q5.5 * scale.psv) + center.psv,
#          Q94.5.r = (Q94.5 * scale.psv) + center.psv)
# 
# dose.psv <- ggplot(model.psv.dr, aes(x = Index, color = factor(Sex.c), fill = factor(Sex.c))) +
#   geom_smooth(aes(y = Estimate.r, ymin = Q5.5.r, ymax = Q94.5.r),
#               stat = "identity",
#               alpha = 0.2) +
#   facet_grid(Side ~ Face, labeller = labeller(Face = face.labs, Side = side.labs), scales = "free") +
#   theme_classic2() +
#   theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
#   scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq, expand = c(0,0)) +
#   scale_y_continuous(name = "Peak Systolic Velocity (cm/s)") +
#   scale_color_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                      labels = c("-0.5" = "Female", "0.5" = "Male"),
#                      breaks = c("0.5","-0.5")) +
#   scale_fill_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                     labels = c("-0.5" = "Female", "0.5" = "Male"),
#                     breaks = c("0.5","-0.5")) +
#   theme(legend.position = "bottom", panel.spacing = unit(1,"lines")) +
#   coord_cartesian(ylim = c(60,120), expand = FALSE)
# dose.psv

#####################
### CCA EDV
#####################

# Find scales
scale.edv <- attributes(df.s$EDV.s)$`scaled:scale`
center.edv <- attributes(df.s$EDV.s)$`scaled:center`
sd.edv <- sd(df.s$EDV[!is.na(df.s$EDV)])

# Plot Data
data.edv <- ggplot(data = df.s, aes(x = Index, color = as.factor(Sex), linetype = Face, y = EDV, size = as.factor(Side))) +
  #geom_jitter(alpha = 0.2) +
  stat_summary(aes(group = interaction(Sex,Face,Side)), fun=mean, geom="line", position = pd) +
  stat_summary(aes(group = interaction(Sex,Face,Side)), fun.data = mean_se, geom="errorbar", width = 1, position = pd, linetype = "solid") +
  stat_summary(aes(group = interaction(Sex,Face,Side)), fun=mean, geom="point", position = pd) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
  scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq) +
  scale_y_continuous(name = "End Diastolic Velocity (cm/s)") +
  scale_color_manual(name = "Sex:", values = c("Female" = "hotpink", "Male" = "royalblue"),
                     breaks = c("Male","Female")) +
  scale_linetype_manual(name = "Position:", values = c("L0" = "solid", "L15" = "dashed"),
                        labels = c("L0" = "0° Supine", "L15" = "15° HDT")) +
  scale_size_manual("Side:", values = c("0.5" = 0.75, "-0.5" = 0.25), labels = c("0.5" = "Right", "-0.5" = "Left"),
                    breaks = c("0.5","-0.5")) +
  theme(legend.position = "right")
# data.edv
# 
# # Define Priors
# priors.edv <- c(prior(normal(0, 1), class = b), # slope prior
#                 prior(normal(0, 1), class = b, coef = "Intercept"), # intercept prior
#                 prior(cauchy(0,1), class = sigma), # population variance
#                 prior(cauchy(0,1), class = sd), # group variance
#                 prior(gamma(2, .1), class = nu) # degrees of freedom
# )
# 
# # # Fit Model
# # model.edv <- brm(EDV.s ~ 0 + Intercept + Index + Sex.c + Face + Side + 
# #                  Index:Sex.c + Sex.c:Face + Sex.c:Side + (1|ID),
# #                  data = df.s,
# #                  prior = priors.edv,
# #                  sample_prior = "yes",
# #                  family = student(),
# #                  warmup = 1000, # burn-in
# #                  iter = 5000, # number of iterations
# #                  chains = 4,  # number of MCMC chains
# #                  control = list(adapt_delta = 0.95),# advanced MC settings
# #                  save_pars=save_pars(all=TRUE))
# # saveRDS(model.edv,file = "BRMSFits/model.edv.rds")
# model.edv <- readRDS("BRMSFits/model.edv.rds")
# 
# # Posterior Checks
# summary(model.edv)
# plot(model.edv, N = 11)
# pp_check(model.edv, ndraws = 50)
# model.edv.loo <- loo(model.edv, moment_match = TRUE, save_psis = TRUE)
# model.edv.loo
# model.edv.w <- weights(model.edv.loo$psis_object)
# ppc_loo_pit_overlay(df.s$EDV.s[!is.na(df.s$EDV.s)], 
#                     posterior_predict(model.edv), 
#                     lw = model.edv.w)
# 
# # Describe Posterior
# describe_posterior(model.edv,
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.05, 0.05), rope_ci = 1,
#                    diagnostic = NULL, effects = "fixed", component = "all")
# 
# # Rescale Posterior
# rescale.edv <- model.edv %>% spread_draws(`b_.*`, sigma, `sd_.*`, nu, regex = TRUE) %>%
#   mutate(b_Intercept = (b_Intercept * scale.edv) + center.edv,
#          b_Index = b_Index * scale.edv,
#          b_Sex.c = b_Sex.c * scale.edv,
#          b_FaceL15 = b_FaceL15 * scale.edv,
#          b_Side = b_Side * scale.edv,
#          `b_Index:Sex.c` = `b_Index:Sex.c` * scale.edv,
#          `b_Sex.c:FaceL15` = `b_Sex.c:FaceL15` * scale.edv,
#          `b_Sex.c:Side` = `b_Sex.c:Side` * scale.edv,
#          sigma = sigma * scale.edv,
#          sd_ID__Intercept = sd_ID__Intercept * scale.edv)
# describe_posterior(rescale.edv[,-c(1:3)],
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.1*sd.edv, 0.1*sd.edv), rope_ci = 1,
#                    diagnostic = NULL)
# 
# # Plot Model
# model.edv.dr <- fitted(model.edv, newdata = nd2, re_formula = NA, probs = c(.055,.945)) %>%
#   data.frame() %>%
#   bind_cols(nd2) %>%
#   mutate(Estimate.r = (Estimate * scale.edv) + center.edv,
#          Q5.5.r = (Q5.5 * scale.edv) + center.edv,
#          Q94.5.r = (Q94.5 * scale.edv) + center.edv)
# 
# dose.edv <- ggplot(model.edv.dr, aes(x = Index, color = factor(Sex.c), fill = factor(Sex.c))) +
#   geom_smooth(aes(y = Estimate.r, ymin = Q5.5.r, ymax = Q94.5.r),
#               stat = "identity",
#               alpha = 0.2) +
#   facet_grid(Side ~ Face, labeller = labeller(Face = face.labs, Side = side.labs), scales = "free") +
#   theme_classic2() +
#   theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
#   scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq, expand = c(0,0)) +
#   scale_y_continuous(name = "End Diastolic Velocity (cm/s)") +
#   scale_color_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                      labels = c("-0.5" = "Female", "0.5" = "Male"),
#                      breaks = c("0.5","-0.5")) +
#   scale_fill_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                     labels = c("-0.5" = "Female", "0.5" = "Male"),
#                     breaks = c("0.5","-0.5")) +
#   theme(legend.position = "bottom", panel.spacing = unit(1,"lines")) +
#   coord_cartesian(ylim = c(25,40), expand = FALSE)
# dose.edv

#####################
### OPP
#####################

# Find scales
scale.opp <- attributes(df.s$OPP.s)$`scaled:scale`
center.opp <- attributes(df.s$OPP.s)$`scaled:center`
sd.opp <- sd(df.s$OPP[!is.na(df.s$OPP)])

# Plot Data
data.opp <- ggplot(data = df.s, aes(x = Index, color = as.factor(Sex), linetype = Face, y = OPP, size = as.factor(Side))) +
  #geom_jitter(alpha = 0.2) +
  stat_summary(aes(group = interaction(Sex,Face,Side)), fun=mean, geom="line", position = pd) +
  stat_summary(aes(group = interaction(Sex,Face,Side)), fun.data = mean_se, geom="errorbar", width = 1, position = pd, linetype = "solid") +
  stat_summary(aes(group = interaction(Sex,Face,Side)), fun=mean, geom="point", position = pd) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
  scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq) +
  scale_y_continuous(name = "Ocular Perfusion Pressure (mmHg)") +
  scale_color_manual(name = "Sex:", values = c("Female" = "hotpink", "Male" = "royalblue"),
                     breaks = c("Male","Female")) +
  scale_linetype_manual(name = "Position:", values = c("L0" = "solid", "L15" = "dashed"),
                        labels = c("L0" = "0° Supine", "L15" = "15° HDT")) +
  scale_size_manual("Side:", values = c("0.5" = 0.75, "-0.5" = 0.25), labels = c("0.5" = "Right", "-0.5" = "Left"),
                    breaks = c("0.5","-0.5")) +
  theme(legend.position = "right")
# data.opp
# 
# # Define Priors
# priors.opp <- c(prior(normal(0, 1), class = b), # slope prior
#                 prior(normal(0, 1), class = b, coef = "Intercept"), # intercept prior
#                 prior(cauchy(0,1), class = sigma), # population variance
#                 prior(cauchy(0,1), class = sd), # group variance
#                 prior(gamma(2, .1), class = nu) # degrees of freedom
# )
# 
# # # Fit Model
# # model.opp <- brm(OPP.s ~ 0 + Intercept + Index + Sex.c + Face + Side + (1|ID),
# #              data = df.s,
# #              prior = priors.opp,
# #              sample_prior = "yes",
# #              family = student(),
# #              warmup = 1000, # burn-in
# #              iter = 5000, # number of iterations
# #              chains = 4,  # number of MCMC chains
# #              control = list(adapt_delta = 0.95),# advanced MC settings
# #              save_pars=save_pars(all=TRUE))
# # saveRDS(model.opp,file = "BRMSFits/model.opp.rds")
# model.opp <- readRDS("BRMSFits/model.opp.rds")
# 
# # Posterior Checks
# summary(model.opp)
# plot(model.opp, N = 8)
# pp_check(model.opp, ndraws = 50)
# model.opp.loo <- loo(model.opp, moment_match = TRUE, save_psis = TRUE)
# model.opp.loo
# model.opp.w <- weights(model.opp.loo$psis_object)
# ppc_loo_pit_overlay(df.s$OPP.s[!is.na(df.s$OPP.s)], 
#                     posterior_predict(model.opp), 
#                     lw = model.opp.w)
# 
# # Describe Posterior
# describe_posterior(model.opp,
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.05, 0.05), rope_ci = 1,
#                    diagnostic = NULL, effects = "fixed", component = "all")
# 
# # Rescale Posterior
# rescale.opp <- model.opp %>% spread_draws(`b_.*`, sigma, `sd_.*`, nu, regex = TRUE) %>%
#   mutate(b_Intercept = (b_Intercept * scale.opp) + center.opp,
#          b_Index = b_Index * scale.opp,
#          b_Sex.c = b_Sex.c * scale.opp,
#          b_FaceL15 = b_FaceL15 * scale.opp,
#          b_Side = b_Side * scale.opp,
#          sigma = sigma * scale.opp,
#          sd_ID__Intercept = sd_ID__Intercept * scale.opp)
# describe_posterior(rescale.opp[,-c(1:3)],
#                    centrality = c("mean","MAP"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "rope"), rope_range = c(-0.1*sd.opp, 0.1*sd.opp), rope_ci = 1,
#                    diagnostic = NULL)
# 
# # Plot Model
# model.opp.dr <- fitted(model.opp, newdata = nd2, re_formula = NA, probs = c(.055,.945)) %>%
#   data.frame() %>%
#   bind_cols(nd2) %>%
#   mutate(Estimate.r = (Estimate * scale.opp) + center.opp,
#          Q5.5.r = (Q5.5 * scale.opp) + center.opp,
#          Q94.5.r = (Q94.5 * scale.opp) + center.opp)
# 
# dose.opp <- ggplot(model.opp.dr, aes(x = Index, color = factor(Sex.c), fill = factor(Sex.c))) +
#   geom_smooth(aes(y = Estimate.r, ymin = Q5.5.r, ymax = Q94.5.r),
#               stat = "identity",
#               alpha = 0.2) +
#   facet_grid(Side ~ Face, labeller = labeller(Face = face.labs, Side = side.labs), scales = "free") +
#   theme_classic2() +
#   theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
#   scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq, expand = c(0,0)) +
#   scale_y_continuous(name = "Ocular Perfusion Pressure (mmHg)") +
#   scale_color_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                      labels = c("-0.5" = "Female", "0.5" = "Male"),
#                      breaks = c("0.5","-0.5")) +
#   scale_fill_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                     labels = c("-0.5" = "Female", "0.5" = "Male"),
#                     breaks = c("0.5","-0.5")) +
#   theme(legend.position = "bottom", panel.spacing = unit(1,"lines")) +
#   coord_cartesian(ylim = c(65,85), expand = FALSE)
# dose.opp

##########################################
### Latent Variables
##########################################

# Find scales
scale.age <- attributes(df.s$Age.s)$`scaled:scale`
center.age <- attributes(df.s$Age.s)$`scaled:center`
sd.age <- sd(df.s$Age.y[!is.na(df.s$Age.y)])
scale.ht <- attributes(df.s$Height.s)$`scaled:scale`
center.ht <- attributes(df.s$Height.s)$`scaled:center`
sd.ht <- sd(df.s$Height[!is.na(df.s$Height)])
scale.wt <- attributes(df.s$Weight.s)$`scaled:scale`
center.wt <- attributes(df.s$Weight.s)$`scaled:center`
sd.wt <- sd(df.s$Weight[!is.na(df.s$Weight)])
scale.bmi <- attributes(df.s$BMI.s)$`scaled:scale`
center.bmi <- attributes(df.s$BMI.s)$`scaled:center`
sd.bmi <- sd(df.s$BMI[!is.na(df.s$BMI)])

df.s <- df.s %>% rowwise() %>%
  mutate(Age.s = Age.s + rnorm(1,0,1e-2),
         Height.s = Height.s + rnorm(1,0,1e-2),
         Weight.s = Weight.s + rnorm(1,0,1e-2),
         BMI.s = BMI.s + rnorm(1,0,1e-2))

#####################
### Age

# # Define Priors
# priors.age <- c(prior(normal(0,0.1), class = sigma), # population variance
#                 prior(normal(0.5,0.01), class = sd) # group variance
# )
# 
# # # Fit Model
# # model.age <- brm(Age.s ~ 0 + (1|ID),
# #              data = df.s,
# #              prior = priors.age,
# #              sample_prior = "yes",
# #              family = gaussian(),
# #              warmup = 1000, # burn-in
# #              iter = 5000, # number of iterations
# #              chains = 4,  # number of MCMC chains
# #              control = list(adapt_delta = 0.95),# advanced MC settings
# #              save_pars=save_pars(all=TRUE))
# # saveRDS(model.age,file = "BRMSFits/model.age.rds")
# model.age <- readRDS("BRMSFits/model.age.rds")
# 
# # Posterior Checks
# summary(model.age)
# plot(model.age, N = 2)
# pp_check(model.age, ndraws = 50)
# model.age.loo <- loo(model.age, moment_match = TRUE, save_psis = TRUE)
# model.age.loo
# model.age.w <- weights(model.age.loo$psis_object)
# ppc_loo_pit_overlay(df.s$Age.s[!is.na(df.s$Age.s)], 
#                     posterior_predict(model.age), 
#                     lw = model.age.w)
# 
# #####################
# ### Height
# 
# # Define Priors
# priors.ht <- c(prior(normal(0,0.1), class = sigma), # population variance
#                 prior(normal(0.5,0.01), class = sd) # group variance
# )
# 
# # # Fit Model
# # model.ht <- brm(Height.s ~ 0 + (1|ID),
# #                  data = df.s,
# #                  prior = priors.ht,
# #                  sample_prior = "yes",
# #                  family = gaussian(),
# #                  warmup = 1000, # burn-in
# #                  iter = 5000, # number of iterations
# #                  chains = 4,  # number of MCMC chains
# #                  control = list(adapt_delta = 0.95),# advanced MC settings
# #                  save_pars=save_pars(all=TRUE))
# # saveRDS(model.ht,file = "BRMSFits/model.ht.rds")
# model.ht <- readRDS("BRMSFits/model.ht.rds")
# 
# # Posterior Checks
# summary(model.ht)
# plot(model.ht, N = 2)
# pp_check(model.ht, ndraws = 50)
# model.ht.loo <- loo(model.ht, moment_match = TRUE, save_psis = TRUE)
# model.ht.loo
# model.ht.w <- weights(model.ht.loo$psis_object)
# ppc_loo_pit_overlay(df.s$Height.s[!is.na(df.s$Height.s)], 
#                     posterior_predict(model.ht), 
#                     lw = model.ht.w)
# 
# #####################
# ### Weight
# 
# # Define Priors
# priors.wt <- c(prior(normal(0,0.1), class = sigma), # population variance
#                prior(normal(0.5,0.01), class = sd) # group variance
# )
# 
# # # Fit Model
# # model.wt <- brm(Weight.s ~ 0 + (1|ID),
# #                  data = df.s,
# #                  prior = priors.wt,
# #                  sample_prior = "yes",
# #                  family = gaussian(),
# #                  warmup = 1000, # burn-in
# #                  iter = 5000, # number of iterations
# #                  chains = 4,  # number of MCMC chains
# #                  control = list(adapt_delta = 0.95),# advanced MC settings
# #                  save_pars=save_pars(all=TRUE))
# # saveRDS(model.wt,file = "BRMSFits/model.wt.rds")
# model.wt <- readRDS("BRMSFits/model.wt.rds")
# 
# # Posterior Checks
# summary(model.wt)
# plot(model.wt, N = 2)
# pp_check(model.wt, ndraws = 50)
# model.wt.loo <- loo(model.wt, moment_match = TRUE, save_psis = TRUE)
# model.wt.loo
# model.wt.w <- weights(model.wt.loo$psis_object)
# ppc_loo_pit_overlay(df.s$Weight.s[!is.na(df.s$Weight.s)], 
#                     posterior_predict(model.wt), 
#                     lw = model.wt.w)
# 
# #####################
# ### BMI
# 
# # Define Priors
# priors.bmi <- c(prior(normal(0,0.1), class = sigma), # population variance
#                prior(normal(0.5,0.01), class = sd) # group variance
# )
# 
# # # Fit Model
# # model.bmi <- brm(BMI.s ~ 0 + (1|ID),
# #                  data = df.s,
# #                  prior = priors.bmi,
# #                  sample_prior = "yes",
# #                  family = gaussian(),
# #                  warmup = 1000, # burn-in
# #                  iter = 5000, # number of iterations
# #                  chains = 4,  # number of MCMC chains
# #                  control = list(adapt_delta = 0.95),# advanced MC settings
# #                  save_pars=save_pars(all=TRUE))
# # saveRDS(model.bmi,file = "BRMSFits/model.bmi.rds")
# model.bmi <- readRDS("BRMSFits/model.bmi.rds")
# 
# # Posterior Checks
# summary(model.bmi)
# plot(model.bmi, N = 2)
# pp_check(model.bmi, ndraws = 50)
# model.bmi.loo <- loo(model.bmi, moment_match = TRUE, save_psis = TRUE)
# model.bmi.loo
# model.bmi.w <- weights(model.bmi.loo$psis_object)
# ppc_loo_pit_overlay(df.s$BMI.s[!is.na(df.s$BMI.s)], 
#                     posterior_predict(model.bmi), 
#                     lw = model.bmi.w)


##########################################
### Multivariate
##########################################

bf_hr <- bf(HR.s ~ 0 + Intercept + Index + Sex.c + Face + (1|p|ID)) + student()
bf_sv <- bf(SV.s ~ 0 + Intercept + Index + Sex.c + Face + Index:Sex.c + (1|p|ID)) + student()
bf_co <- bf(CO.s ~ 0 + Intercept + Index + Sex.c + Face + (1|p|ID)) + student()
bf_vo2 <- bf(VO2.s ~ 0 + Intercept + Index + Sex.c + Face + Sex.c:Face + (1|p|ID)) + student()
bf_sbp <- bf(SBP.s ~ 0 + Intercept + Index + Sex.c + Face + (1|p|ID)) + student()
bf_dbp <- bf(DBP.s ~ 0 + Intercept + Index + Sex.c + Face + (1|p|ID)) + student()
bf_rpp <- bf(RPP.s ~ 0 + Intercept + Index + Sex.c + Face + (1|p|ID)) + gaussian()
bf_mo <- bf(MO.s ~ 0 + Intercept + Index + Sex.c + Face + (1|p|ID)) + student()
bf_tpr <- bf(TPR.s ~ 0 + Intercept + Index + Sex.c + Face + (1|p|ID)) + student()
bf_ci <- bf(CI.s ~ 0 + Intercept + Index + Sex.c + Face + (1|p|ID)) + student()
bf_si <- bf(SI.s ~ 0 + Intercept + Index + Sex.c + Face + (1|p|ID)) + student()
bf_sdnn <- bf(SDNN.s ~ 0 + Intercept + Index + Sex.c + Face + (1|p|ID)) + student()
bf_hrvti <- bf(HRVTI.s ~ 0 + Intercept + Inde$x + Sex.c + Face + (1|p|ID)) + student()
bf_rmsdd <- bf(RMSDD.s ~ 0 + Intercept + Index + Sex.c + Face + (1|p|ID)) + student()
bf_brs <- bf(BRS.s ~ 0 + Intercept + Index + Sex.c + Face + (1|p|ID)) + student()
bf_lfn <- bf(LFNorm.s ~ 0 + Intercept + Index + Sex.c + Face + (1|p|ID)) + student()
bf_hfn <- bf(HFNorm.s ~ 0 + Intercept + Index + Sex.c + Face + (1|p|ID)) + student()
bf_lfhf <- bf(LFHF.s ~ 0 + Intercept + Index + Sex.c + Face + Index:Sex.c + (1|p|ID)) + student()
bf_iop <- bf(IOP.s ~ 0 + Intercept + Index + Sex.c + Face + Side + (1|p|ID)) + student()
bf_vpc <- bf(VPC.s ~ 0 + Intercept + Index + Sex.c + Face + Side + (1|p|ID), sigma ~ Index) + student()
bf_ijvf <- bf(IJVF ~ Index + Sex.c + Face + Side + (1|p|ID)) + cumulative(link = "logit", threshold = "flexible")
bf_opp <- bf(OPP.s ~ 0 + Intercept + Index + Sex.c + Face + Side + (1|p|ID)) + student()
bf_psv <- bf(PSV.s ~ 0 + Intercept + Index + Sex.c + Face + Side + (1|p|ID)) + student()
bf_edv <- bf(EDV.s ~ 0 + Intercept + Index + Sex.c + Face + Side + Index:Sex.c + Sex.c:Face + Sex.c:Side + (1|p|ID)) + student()
bf_cca <- bf(CCA.s ~ 0 + Intercept + Index + Sex.c + Face + Side + (1|p|ID)) + student()
bf_ijv <- bf(IJV.s ~ 0 + Intercept + Index + Sex.c + Face + Side + (1|p|ID)) + student()
bf_age <- bf(Age.s ~ 0 + (1|p|ID)) + gaussian()
bf_ht <- bf(Height.s ~ 0 + (1|p|ID)) + gaussian()
bf_wt <- bf(Weight.s ~ 0 + (1|p|ID)) + gaussian()
bf_bmi <- bf(BMI.s ~ 0 + (1|p|ID)) + gaussian()

priors.mult <- c(prior(normal(0, 1), class = b, resp = HRs), # slope prior
                 prior(normal(0, 1), class = b, resp = SVs), # slope prior
                 prior(normal(0, 1), class = b, resp = COs), # slope prior
                 prior(normal(0, 1), class = b, resp = VO2s), # slope prior
                 prior(normal(0, 1), class = b, resp = SBPs), # slope prior
                 prior(normal(0, 1), class = b, resp = DBPs), # slope prior
                 prior(normal(0, 1), class = b, resp = RPPs), # slope prior
                 prior(normal(0, 1), class = b, resp = MOs), # slope prior
                 prior(normal(0, 1), class = b, resp = TPRs), # slope prior
                 prior(normal(0, 1), class = b, resp = CIs), # slope prior
                 prior(normal(0, 1), class = b, resp = SIs), # slope prior
                 prior(normal(0, 1), class = b, resp = SDNNs), # slope prior
                 prior(normal(0, 1), class = b, resp = HRVTIs), # slope prior
                 prior(normal(0, 1), class = b, resp = RMSDDs), # slope prior
                 prior(normal(0, 1), class = b, resp = BRSs), # slope prior
                 prior(normal(0, 1), class = b, resp = LFNorms), # slope prior
                 prior(normal(0, 1), class = b, resp = HFNorms), # slope prior
                 prior(normal(0, 1), class = b, resp = LFHFs), # slope prior
                 prior(normal(0, 1), class = b, resp = IOPs), # slope prior
                 prior(normal(0, 1), class = b, resp = VPCs), # slope prior
                 prior(normal(0, 1), class = b, resp = IJVF), # slope prior
                 prior(normal(0, 1), class = b, resp = OPPs), # slope prior
                 prior(normal(0, 1), class = b, resp = PSVs), # slope prior
                 prior(normal(0, 1), class = b, resp = EDVs), # slope prior
                 prior(normal(0, 1), class = b, resp = CCAs), # slope prior
                 prior(normal(0, 1), class = b, resp = IJVs), # slope prior
                 prior(normal(0, 10), class = Intercept, resp = IJVF), # slope prior
                 prior(normal(0, 1), class = b, coef = Intercept, resp = HRs), # intercept prior
                 prior(normal(0, 1), class = b, coef = Intercept, resp = SVs), # intercept prior
                 prior(normal(0, 1), class = b, coef = Intercept, resp = COs), # intercept prior
                 prior(normal(0, 1), class = b, coef = Intercept, resp = VO2s), # intercept prior
                 prior(normal(0, 1), class = b, coef = Intercept, resp = SBPs), # intercept prior
                 prior(normal(0, 1), class = b, coef = Intercept, resp = DBPs), # intercept prior
                 prior(normal(0, 1), class = b, coef = Intercept, resp = RPPs), # intercept prior
                 prior(normal(0, 1), class = b, coef = Intercept, resp = MOs), # intercept prior
                 prior(normal(0, 1), class = b, coef = Intercept, resp = TPRs), # intercept prior
                 prior(normal(0, 1), class = b, coef = Intercept, resp = CIs), # intercept prior
                 prior(normal(0, 1), class = b, coef = Intercept, resp = SIs), # intercept prior
                 prior(normal(0, 1), class = b, coef = Intercept, resp = SDNNs), # intercept prior
                 prior(normal(0, 1), class = b, coef = Intercept, resp = HRVTIs), # intercept prior
                 prior(normal(0, 1), class = b, coef = Intercept, resp = RMSDDs), # intercept prior
                 prior(normal(0, 1), class = b, coef = Intercept, resp = BRSs), # intercept prior
                 prior(normal(0, 1), class = b, coef = Intercept, resp = LFNorms), # intercept prior
                 prior(normal(0, 1), class = b, coef = Intercept, resp = HFNorms), # intercept prior
                 prior(normal(0, 1), class = b, coef = Intercept, resp = LFHFs), # intercept prior
                 prior(normal(0, 1), class = b, coef = Intercept, resp = IOPs), # intercept prior
                 prior(normal(0, 1), class = b, coef = Intercept, resp = VPCs), # intercept prior
                 prior(normal(0, 1), class = b, coef = Intercept, resp = OPPs), # intercept prior
                 prior(normal(0, 1), class = b, coef = Intercept, resp = PSVs), # intercept prior
                 prior(normal(0, 1), class = b, coef = Intercept, resp = EDVs), # intercept prior
                 prior(normal(0, 1), class = b, coef = Intercept, resp = CCAs), # intercept prior
                 prior(normal(0, 1), class = b, coef = Intercept, resp = IJVs), # intercept prior
                 prior(cauchy(0, 1), class = sigma, resp = HRs), # population variance
                 prior(cauchy(0, 1), class = sigma, resp = SVs), # population variance
                 prior(cauchy(0, 1), class = sigma, resp = COs), # population variance
                 prior(cauchy(0, 1), class = sigma, resp = VO2s), # population variance
                 prior(cauchy(0, 1), class = sigma, resp = SBPs), # population variance
                 prior(cauchy(0, 1), class = sigma, resp = DBPs), # population variance
                 prior(cauchy(0, 1), class = sigma, resp = RPPs), # population variance
                 prior(cauchy(0, 1), class = sigma, resp = MOs), # population variance
                 prior(cauchy(0, 1), class = sigma, resp = TPRs), # population variance
                 prior(cauchy(0, 1), class = sigma, resp = CIs), # population variance
                 prior(cauchy(0, 1), class = sigma, resp = SIs), # population variance
                 prior(cauchy(0, 1), class = sigma, resp = SDNNs), # population variance
                 prior(cauchy(0, 1), class = sigma, resp = HRVTIs), # population variance
                 prior(cauchy(0, 1), class = sigma, resp = RMSDDs), # population variance
                 prior(cauchy(0, 1), class = sigma, resp = BRSs), # population variance
                 prior(cauchy(0, 1), class = sigma, resp = LFNorms), # population variance
                 prior(cauchy(0, 1), class = sigma, resp = HFNorms), # population variance
                 prior(cauchy(0, 1), class = sigma, resp = LFHFs), # population variance
                 prior(cauchy(0, 1), class = sigma, resp = IOPs), # population variance
                 prior(normal(0,1), class = b, dpar = sigma, resp = VPCs), # population variance slope
                 prior(cauchy(0, 1), class = sigma, resp = OPPs), # population variance
                 prior(cauchy(0, 1), class = sigma, resp = PSVs), # population variance
                 prior(cauchy(0, 1), class = sigma, resp = EDVs), # population variance
                 prior(cauchy(0, 1), class = sigma, resp = CCAs), # population variance
                 prior(cauchy(0, 1), class = sigma, resp = IJVs), # population variance
                 prior(cauchy(0, 1), class = sd, resp = HRs), # group variance
                 prior(cauchy(0, 1), class = sd, resp = SVs), # group variance
                 prior(cauchy(0, 1), class = sd, resp = COs), # group variance
                 prior(cauchy(0, 1), class = sd, resp = VO2s), # group variance
                 prior(cauchy(0, 1), class = sd, resp = SBPs), # group variance
                 prior(cauchy(0, 1), class = sd, resp = DBPs), # group variance
                 prior(cauchy(0, 1), class = sd, resp = RPPs), # group variance
                 prior(cauchy(0, 1), class = sd, resp = MOs), # group variance
                 prior(cauchy(0, 1), class = sd, resp = TPRs), # group variance
                 prior(cauchy(0, 1), class = sd, resp = CIs), # group variance
                 prior(cauchy(0, 1), class = sd, resp = SIs), # group variance
                 prior(cauchy(0, 1), class = sd, resp = SDNNs), # group variance
                 prior(cauchy(0, 1), class = sd, resp = HRVTIs), # group variance
                 prior(cauchy(0, 1), class = sd, resp = RMSDDs), # group variance
                 prior(cauchy(0, 1), class = sd, resp = BRSs), # group variance
                 prior(cauchy(0, 1), class = sd, resp = LFNorms), # group variance
                 prior(cauchy(0, 1), class = sd, resp = HFNorms), # group variance
                 prior(cauchy(0, 1), class = sd, resp = LFHFs), # group variance
                 prior(cauchy(0, 1), class = sd, resp = IOPs), # group variance
                 prior(cauchy(0, 1), class = sd, resp = VPCs), # group variance
                 prior(cauchy(0, 1), class = sd, resp = IJVF), # group variance
                 prior(cauchy(0, 1), class = sd, resp = OPPs), # group variance
                 prior(cauchy(0, 1), class = sd, resp = PSVs), # group variance
                 prior(cauchy(0, 1), class = sd, resp = EDVs), # group variance
                 prior(cauchy(0, 1), class = sd, resp = CCAs), # group variance
                 prior(cauchy(0, 1), class = sd, resp = IJVs), # group variance
                 prior(gamma(2, 0.1), class = nu, resp = HRs), # degrees of freedom
                 prior(gamma(2, 0.1), class = nu, resp = SVs), # degrees of freedom
                 prior(gamma(2, 0.1), class = nu, resp = COs), # degrees of freedom
                 prior(gamma(2, 0.1), class = nu, resp = VO2s), # degrees of freedom
                 prior(gamma(2, 0.1), class = nu, resp = SBPs), # degrees of freedom
                 prior(gamma(2, 0.1), class = nu, resp = DBPs), # degrees of freedom
                 prior(gamma(2, 0.1), class = nu, resp = MOs), # degrees of freedom
                 prior(gamma(2, 0.1), class = nu, resp = TPRs), # degrees of freedom
                 prior(gamma(2, 0.1), class = nu, resp = CIs), # degrees of freedom
                 prior(gamma(2, 0.1), class = nu, resp = SIs), # degrees of freedom
                 prior(gamma(2, 0.1), class = nu, resp = SDNNs), # degrees of freedom
                 prior(gamma(2, 0.1), class = nu, resp = HRVTIs), # degrees of freedom
                 prior(gamma(2, 0.1), class = nu, resp = RMSDDs), # degrees of freedom
                 prior(gamma(2, 0.1), class = nu, resp = BRSs), # degrees of freedom
                 prior(gamma(2, 0.1), class = nu, resp = LFNorms), # degrees of freedom
                 prior(gamma(2, 0.1), class = nu, resp = HFNorms), # degrees of freedom
                 prior(gamma(2, 0.1), class = nu, resp = LFHFs), # degrees of freedom
                 prior(gamma(2, 0.1), class = nu, resp = IOPs), # degrees of freedom
                 prior(gamma(2, 0.1), class = nu, resp = VPCs), # degrees of freedom
                 prior(gamma(2, 0.1), class = nu, resp = OPPs), # degrees of freedom
                 prior(gamma(2, 0.1), class = nu, resp = PSVs), # degrees of freedom
                 prior(gamma(2, 0.1), class = nu, resp = EDVs), # degrees of freedom
                 prior(gamma(2, 0.1), class = nu, resp = CCAs), # degrees of freedom
                 prior(gamma(2, 0.1), class = nu, resp = IJVs), # degrees of freedom
                 prior(normal(0,0.1), class = sigma, resp = Ages), # population variance
                 prior(normal(0,0.1), class = sigma, resp = Heights), # population variance
                 prior(normal(0,0.1), class = sigma, resp = Weights), # population variance
                 prior(normal(0,0.1), class = sigma, resp = BMIs), # population variance
                 prior(normal(0.5,0.01), class = sd, resp = Ages), # population variance
                 prior(normal(0.5,0.01), class = sd, resp = Heights), # group variance
                 prior(normal(0.5,0.01), class = sd, resp = Weights), # group variance
                 prior(normal(0.5,0.01), class = sd, resp = BMIs) # group variance
)

# # Fit Model
# model.mult <- brm(bf_hr + bf_sv + bf_co + bf_vo2 + bf_sbp + bf_dbp + bf_rpp + bf_mo + bf_tpr +
#                     bf_ci + bf_si +
#                     bf_sdnn + bf_hrvti + bf_rmsdd + bf_brs +
#                     bf_lfn + bf_hfn + bf_lfhf +
#                     bf_iop + bf_vpc + bf_ijvf + bf_opp +
#                     bf_psv + bf_edv + bf_cca + bf_ijv +
#                     bf_age + bf_ht + bf_wt + bf_bmi +
#                     set_rescor(FALSE),
#              data = df.s,
#              prior = priors.mult,
#              sample_prior = "yes",
#              warmup = 1000, # burn-in
#              iter = 20000, # number of iterations
#              chains = 4,  # number of MCMC chains
#              init_r = 0.1,
#              control = list(adapt_delta = 0.95, max_treedepth = 15),# advanced MC settings
#              save_pars=save_pars(all=TRUE))
# test <- stancode(model.mult)
# write.table(stancode(model.mult), "test.txt")
# saveRDS(model.mult,file = "BRMSFits/model.mult.rds")
model.mult <- readRDS("BRMSFits/model.mult.rds")
# summary(model.mult)
# test <- prior_summary(model.mult)

#####################
### Pressure Effect
#####################

fig.pres1 <- model.mult %>% gather_draws(`b_.*._Index`, regex = TRUE) %>%
  mutate(.variable = str_remove_all(.variable,"b_|s_Index")) %>%
  mutate(type = case_when(.variable %in% c("HR","SV","CO","VO2","SBP","DBP","RPP","MO","TPR","CI","SI") ~ "A",
                          .variable %in% c("SDNN","HRVTI","RMSDD","BRS","LFNorm","HFNorm","LFHF") ~ "B",
                          .variable %in% c("IOP","VPC","OPP","PSV","EDV","CCA","IJV") ~ "C",
                          TRUE ~ "D")) %>%
  filter(.variable != ("sigma_VPC")) %>% filter(.variable != ("IJVF_Index")) %>%
  ungroup() %>%
  mutate(.variable = reorder(.variable, .value)) %>%
  ggplot(aes(y = .variable, x = .value, fill = stat(abs(x) < .05))) +
  stat_halfeye(aes(color = type), .width = c(0.89,0.95), point_interval = "mode_hdi") +
  geom_vline(xintercept = c(-.05, .05), linetype = "dashed") +
  scale_fill_manual("ROPE:", values = c("gray80", "skyblue")) +
  scale_color_manual("Type:", values = c("A" = "lightgreen","B" = "tan1", "C" = "plum"), labels = c("A" = "Systemic","B" = "Autonomic", "C" = "Head/Neck")) +
  scale_x_continuous("Normalized Effect Size (with increasing negative pressure)") +
  scale_y_discrete("Parameter", labels = c("CI","SI","SV","CO","IJVP",expression(A[IJV]),"IOP","PSV","HFNorm","RMSDD","SBP",
                                           "VO2","BRS",expression(A[CCA]),"SDNN","HRVTi","DBP","OPP","EDV","LFNorm","HR",
                                           "RPP","LF/HF","MO","TPR")) +
  theme_classic() +
  theme(legend.position = "bottom", panel.grid.major.y = element_line(color = "gray80"))+
  labs(subtitle = "ROPE") + theme(plot.subtitle = element_text(hjust = 0.52))
fig.pres2 <- model.mult %>% gather_draws(`b_.*._Index`, regex = TRUE) %>%
  mutate(.variable = str_remove_all(.variable,"b_|s_Index")) %>%
  mutate(type = case_when(.variable %in% c("HR","SV","CO","VO2","SBP","DBP","RPP","MO","TPR","CI","SI") ~ "A",
                          .variable %in% c("SDNN","HRVTI","RMSDD","BRS","LFNorm","HFNorm","LFHF") ~ "B",
                          .variable %in% c("IOP","VPC","OPP","PSV","EDV","CCA","IJV","IJVF_Index") ~ "C",
                          TRUE ~ "D")) %>%
  filter(.variable != ("sigma_VPC")) %>% filter(.variable == ("IJVF_Index")) %>%
  mutate(.variable = "IJVF") %>%
  ungroup() %>%
  mutate(.variable = reorder(.variable, .value)) %>%
  ggplot(aes(y = .variable, x = .value, fill = stat(abs(x) < 0.1*pi/sqrt(3)))) +
  stat_halfeye(aes(color = type), .width = c(0.89,0.95), point_interval = "mode_hdi") +
  geom_vline(xintercept = c(-0.1*pi/sqrt(3), 0.1*pi/sqrt(3)), linetype = "dashed") +
  scale_fill_manual(name = expression(paste("ROPE (±",0.1*pi/sqrt(3),"):",sep="")), values = c("gray80", "skyblue")) +
  scale_color_manual("Type:", values = c("A" = "lightgreen","B" = "tan1", "C" = "plum"), labels = c("A" = "Systemic","B" = "Autonomic", "C" = "Head/Neck")) +
  scale_x_continuous("Normalized Effect Size (with increasing negative pressure): Logistic Regression") +
  scale_y_discrete(name = element_blank()) +
  theme_classic() +
  theme(legend.position = "bottom", panel.grid.major.y = element_line(color = "gray80"))+
  labs(subtitle = "ROPE") + theme(plot.subtitle = element_text(hjust = 0.92))

fig.pres <- ggarrange(fig.pres1,fig.pres2, ncol = 1, common.legend = TRUE, legend = "bottom", heights = c(10, 1.75), align = "v")
ggsave("LBNP-fig_dr_pres.png", fig.pres, width = 8, height = 7, units = "in")


#####################
### Sex Effect
#####################

fig.sex1 <- model.mult %>% gather_draws(`b_.*._Sex.c`, regex = TRUE) %>%
  mutate(.variable = str_remove_all(.variable,"b_|s_Sex\\.c")) %>%
  mutate(type = case_when(.variable %in% c("HR","SV","CO","VO2","SBP","DBP","RPP","MO","TPR","CI","SI") ~ "A",
                          .variable %in% c("SDNN","HRVTI","RMSDD","BRS","LFNorm","HFNorm","LFHF") ~ "B",
                          .variable %in% c("IOP","VPC","OPP","PSV","EDV","CCA","IJV") ~ "C",
                          TRUE ~ "D")) %>%
  filter(.variable != ("IJVF_Sex.c")) %>%
  ungroup() %>%
  mutate(.variable = reorder(.variable, .value)) %>%
  ggplot(aes(y = .variable, x = .value, fill = stat(abs(x) < .05))) +
  stat_halfeye(aes(color = type), .width = c(0.89,0.95), point_interval = "mode_hdi") +
  geom_vline(xintercept = c(-.05, .05), linetype = "dashed") +
  scale_fill_manual("ROPE:", values = c("gray80", "skyblue")) +
  scale_color_manual("Type:", values = c("A" = "lightgreen","B" = "tan1", "C" = "plum"), labels = c("A" = "Systemic","B" = "Autonomic", "C" = "Head/Neck")) +
  scale_x_continuous("Normalized Effect Size (from Female to Male)") +
  scale_y_discrete("Parameter", labels = c("HFNorm","TPR","HR","RPP","LF/HF","IOP","BRS","CI","MO","EDV","OPP",
                                           expression(A[CCA]),"IJVP","RMSDD","DBP",expression(A[IJV]),"SI","LFNorm","HRVTi","CO","SBP",
                                           "PSV","SDNN","VO2","SV")) +
  theme_classic() +
  theme(legend.position = "bottom", panel.grid.major.y = element_line(color = "gray80")) +
  coord_cartesian(xlim = c(-0.8,1))+
labs(subtitle = "ROPE") + theme(plot.subtitle = element_text(hjust = 0.45))

fig.sex2 <- model.mult %>% gather_draws(`b_.*._Sex.c`, regex = TRUE) %>%
  mutate(.variable = str_remove_all(.variable,"b_|s_FaceL15")) %>%
  mutate(type = case_when(.variable %in% c("HR","SV","CO","VO2","SBP","DBP","RPP","MO","TPR","CI","SI") ~ "A",
                          .variable %in% c("SDNN","HRVTI","RMSDD","BRS","LFNorm","HFNorm","LFHF") ~ "B",
                          .variable %in% c("IOP","VPC","OPP","PSV","EDV","CCA","IJV","IJVF_Sex.c") ~ "C",
                          TRUE ~ "D")) %>%
  filter(.variable != ("sigma_VPC")) %>% filter(.variable == ("IJVF_Sex.c")) %>%
  mutate(.variable = "IJVF") %>%
  ungroup() %>%
  mutate(.variable = reorder(.variable, .value)) %>%
  ggplot(aes(y = .variable, x = .value, fill = stat(abs(x) < 0.1*pi/sqrt(3)))) +
  stat_halfeye(aes(color = type), .width = c(0.89,0.95), point_interval = "mode_hdi") +
  geom_vline(xintercept = c(-0.1*pi/sqrt(3), 0.1*pi/sqrt(3)), linetype = "dashed") +
  scale_fill_manual(name = expression(paste("ROPE (±",0.1*pi/sqrt(3),"):",sep="")), values = c("gray80", "skyblue")) +
  scale_color_manual("Type:", values = c("A" = "lightgreen","B" = "tan1", "C" = "plum"), labels = c("A" = "Systemic","B" = "Autonomic", "C" = "Head/Neck")) +
  scale_x_continuous("Normalized Effect Size (from Female to Male): Logistic Regression") +
  scale_y_discrete(name = element_blank()) +
  theme_classic() +
  theme(legend.position = "bottom", panel.grid.major.y = element_line(color = "gray80"))+
  labs(subtitle = "ROPE") + theme(plot.subtitle = element_text(hjust = 0.48))

fig.sex <- ggarrange(fig.sex1,fig.sex2, ncol = 1, common.legend = TRUE, legend = "bottom", heights = c(10, 1.75), align = "v")
ggsave("LBNP-fig_dr_sex.png", fig.sex, width = 8, height = 7, units = "in")


#####################
### Position Effect
#####################

fig.pos1 <- model.mult %>% gather_draws(`b_.*._FaceL15`, regex = TRUE) %>%
  mutate(.variable = str_remove_all(.variable,"b_|s_FaceL15")) %>%
  mutate(type = case_when(.variable %in% c("HR","SV","CO","VO2","SBP","DBP","RPP","MO","TPR","CI","SI") ~ "A",
                          .variable %in% c("SDNN","HRVTI","RMSDD","BRS","LFNorm","HFNorm","LFHF") ~ "B",
                          .variable %in% c("IOP","VPC","OPP","PSV","EDV","CCA","IJV") ~ "C",
                          TRUE ~ "D")) %>%
  filter(.variable != ("IJVF_FaceL15")) %>%
  ungroup() %>%
  mutate(.variable = reorder(.variable, .value)) %>%
  ggplot(aes(y = .variable, x = .value, fill = stat(abs(x) < .05))) +
  stat_halfeye(aes(color = type), .width = c(0.89,0.95), point_interval = "mode_hdi") +
  geom_vline(xintercept = c(-.05, .05), linetype = "dashed") +
  scale_fill_manual("ROPE:", values = c("gray80", "skyblue")) +
  scale_color_manual("Type:", values = c("A" = "lightgreen","B" = "tan1", "C" = "plum"), labels = c("A" = "Systemic","B" = "Autonomic", "C" = "Head/Neck")) +
  scale_x_continuous("Normalized Effect Size (from 0° to 15° HDT)") +
  scale_y_discrete("Parameter", labels = c("TPR","HR","LF/HF","LFNorm","PSV","MO","RPP","EDV","SBP","DBP","HRVTi",
                                           "VO2","HFNorm","SDNN","OPP","RMSDD","CO","BRS","SV","CI","SI",
                                           expression(A[CCA]),"IOP","IJVP",expression(A[IJV]))) +
  theme_classic() +
  theme(legend.position = "bottom", panel.grid.major.y = element_line(color = "gray80"))+
  labs(subtitle = "ROPE") + theme(plot.subtitle = element_text(hjust = 0.35))

fig.pos2 <- model.mult %>% gather_draws(`b_.*._FaceL15`, regex = TRUE) %>%
  mutate(.variable = str_remove_all(.variable,"b_|s_FaceL15")) %>%
  mutate(type = case_when(.variable %in% c("HR","SV","CO","VO2","SBP","DBP","RPP","MO","TPR","CI","SI") ~ "A",
                          .variable %in% c("SDNN","HRVTI","RMSDD","BRS","LFNorm","HFNorm","LFHF") ~ "B",
                          .variable %in% c("IOP","VPC","OPP","PSV","EDV","CCA","IJV","IJVF_FaceL15") ~ "C",
                          TRUE ~ "D")) %>%
  filter(.variable != ("sigma_VPC")) %>% filter(.variable == ("IJVF_FaceL15")) %>%
  mutate(.variable = "IJVF") %>%
  ungroup() %>%
  mutate(.variable = reorder(.variable, .value)) %>%
  ggplot(aes(y = .variable, x = .value, fill = stat(abs(x) < 0.1*pi/sqrt(3)))) +
  stat_halfeye(aes(color = type), .width = c(0.89,0.95), point_interval = "mode_hdi") +
  geom_vline(xintercept = c(-0.1*pi/sqrt(3), 0.1*pi/sqrt(3)), linetype = "dashed") +
  scale_fill_manual(name = expression(paste("ROPE (±",0.1*pi/sqrt(3),"):",sep="")), values = c("gray80", "skyblue")) +
  scale_color_manual("Type:", values = c("A" = "lightgreen","B" = "tan1", "C" = "plum"), labels = c("A" = "Systemic","B" = "Autonomic", "C" = "Head/Neck")) +
  scale_x_continuous("Normalized Effect Size (from 0° to 15° HDT): Logistic Regression") +
  scale_y_discrete(name = element_blank()) +
  theme_classic() +
  theme(legend.position = "bottom", panel.grid.major.y = element_line(color = "gray80"))+
  labs(subtitle = "ROPE") + theme(plot.subtitle = element_text(hjust = 0.06))

fig.pos <- ggarrange(fig.pos1,fig.pos2, ncol = 1, common.legend = TRUE, legend = "bottom", heights = c(10, 1.75), align = "v")
ggsave("LBNP-fig_dr_pos.png", fig.pos, width = 8, height = 7, units = "in")


#####################
### Side Effect
#####################

fig.side1 <- model.mult %>% gather_draws(`b_.*._Side`, regex = TRUE) %>%
  mutate(.variable = str_remove_all(.variable,"b_|s_Side")) %>%
  mutate(type = case_when(.variable %in% c("HR","SV","CO","VO2","SBP","DBP","RPP","MO","TPR","CI","SI") ~ "A",
                          .variable %in% c("SDNN","HRVTI","RMSDD","BRS","LFNorm","HFNorm","LFHF") ~ "B",
                          .variable %in% c("IOP","VPC","OPP","PSV","EDV","CCA","IJV") ~ "C",
                          TRUE ~ "D")) %>%
  filter(.variable != ("IJVF_Side")) %>%
  ungroup() %>%
  mutate(.variable = reorder(.variable, .value)) %>%
  ggplot(aes(y = .variable, x = .value, fill = stat(abs(x) < .05))) +
  stat_halfeye(aes(color = type), .width = c(0.89,0.95), point_interval = "mode_hdi") +
  geom_vline(xintercept = c(-.05, .05), linetype = "dashed") +
  scale_fill_manual("ROPE:", values = c("gray80", "skyblue")) +
  scale_color_manual("Type:", values = c("A" = "lightgreen","B" = "tan1", "C" = "plum"), labels = c("A" = "Systemic","B" = "Autonomic", "C" = "Head/Neck")) +
  scale_x_continuous("Normalized Effect Size (from left to right)") +
  scale_y_discrete("Parameter", labels = c("OPP",expression(A[CCA]),"EDV","IOP",expression(A[IJV]),"IJVP","PSV")) +
  theme_classic() +
  theme(legend.position = "bottom", panel.grid.major.y = element_line(color = "gray80"))+
  labs(subtitle = "ROPE") + theme(plot.subtitle = element_text(hjust = 0.51))
#coord_cartesian(xlim = c(-0.25,0.25))

fig.side2 <- model.mult %>% gather_draws(`b_.*._Side`, regex = TRUE) %>%
  mutate(.variable = str_remove_all(.variable,"b_|s_FaceL15")) %>%
  mutate(type = case_when(.variable %in% c("HR","SV","CO","VO2","SBP","DBP","RPP","MO","TPR","CI","SI") ~ "A",
                          .variable %in% c("SDNN","HRVTI","RMSDD","BRS","LFNorm","HFNorm","LFHF") ~ "B",
                          .variable %in% c("IOP","VPC","OPP","PSV","EDV","CCA","IJV","IJVF_Side") ~ "C",
                          TRUE ~ "D")) %>%
  filter(.variable != ("sigma_VPC")) %>% filter(.variable == ("IJVF_Side")) %>%
  mutate(.variable = "IJVF") %>%
  ungroup() %>%
  mutate(.variable = reorder(.variable, .value)) %>%
  ggplot(aes(y = .variable, x = .value, fill = stat(abs(x) < 0.1*pi/sqrt(3)))) +
  stat_halfeye(aes(color = type), .width = c(0.89,0.95), point_interval = "mode_hdi") +
  geom_vline(xintercept = c(-0.1*pi/sqrt(3), 0.1*pi/sqrt(3)), linetype = "dashed") +
  scale_fill_manual(name = expression(paste("ROPE (±",0.1*pi/sqrt(3),"):",sep="")), values = c("gray80", "skyblue")) +
  scale_color_manual("Type:", values = c("A" = "lightgreen","B" = "tan1", "C" = "plum"), labels = c("A" = "Systemic","B" = "Autonomic", "C" = "Head/Neck")) +
  scale_x_continuous("Normalized Effect Size (from left to right): Logistic Regression") +
  scale_y_discrete(name = element_blank()) +
  theme_classic() +
  theme(legend.position = "bottom", panel.grid.major.y = element_line(color = "gray80"))+
  labs(subtitle = "ROPE") + theme(plot.subtitle = element_text(hjust = 0.475))

fig.side <- ggarrange(fig.side1,fig.side2, ncol = 1, common.legend = TRUE, legend = "bottom", heights = c(10, 3), align = "v")
ggsave("LBNP-fig_dr_side.png", fig.side, width = 8, height = 7, units = "in")


#####################
### Save Figures
#####################

fig.effect <- ggarrange(fig.pres,fig.sex,fig.pos, ncol = 3, labels = LETTERS[1:3],
                        common.legend = TRUE, legend = "bottom")
fig.effect

ggsave("LBNP.Effect.png",fig.effect, width = 15, height = 7, units = "in")

fig.effect2 <- ggarrange(fig.pres1+labs(title = "Pressure"),fig.pos1+labs(title = "Position"),
                        ncol = 2, nrow = 1, labels = LETTERS[1:2], align = "hv",
                        common.legend = TRUE, legend = "bottom")
ggsave("LBNP.Effect2.png",fig.effect2, width = 15, height = 7, units = "in")

ggsave("LBNP.Pres.png",fig.pres, width = 5, height = 8, units = "in")

ggarrange(data.iop,data.opp, nrow = 1, ncol = 2, align = "hv",
          common.legend = TRUE, legend = "bottom")


######################
### 15° HUT Hypothesis
#####################

fig.pos1
fig.hut <- fig.pos1 + scale_x_reverse("Normalized Effect Size (from 0° to 15° HUT)", breaks = c(0.6,0.4,0.2,0.0,-0.2),
                           labels = c("-0.6","-0.4","-0.2","0.0","0.2")) +
  labs(subtitle = "ROPE") + theme(plot.subtitle = element_text(hjust = 0.65))

fig.effect3 <- ggarrange(fig.pres1+labs(title = "Pressure"),fig.hut+labs(title = "Position"),
                         ncol = 2, nrow = 1, labels = LETTERS[1:2], align = "hv",
                         common.legend = TRUE, legend = "bottom")
ggsave("LBNP.Effect3.png",fig.effect3, width = 15, height = 7, units = "in")

#####################
### Graph Correlation Structure
#####################

# test <- VarCorr(model.mult)
#
# describe_posterior(model.mult,
#                    centrality = c("mean"), dispersion = TRUE,
#                    ci = 0.89, ci_method = "hdi",
#                    test = c("p_direction", "bf"), bf_prior, rope_ci = 1,
#                    diagnostic = NULL, effects = "random", component = "all")
#
# summary(model.mult)
# model.mult %>% gather_draws(`cor_ID__HRs_Intercept__SVs_Intercept`, regex = TRUE) %>% describe_posterior(test = "bf", bf_prior = uniform(0,1))
# h1 <- hypothesis(model.mult, "ID__DBPs_Intercept__TPRs_Intercept=0", class = "cor")
# 1 / h1$hypothesis$Evid.Ratio



cor.mult <- VarCorr(model.mult)$ID$cor[,1,]
rownames(cor.mult) <- c("HR", "SV", "CO", "VO2", "SBP", "DBP", "RPP", "MO", "TPR", "CI", "SI", "SDNN", "HRVTI",
                        "RMSDD", "BRS", "LF", "HF", "LF.HF", "IOP", "IJVP", "IJVF", "OPP", "PSV", "EDV", "CCA", "IJV", "Age", "Ht", "Wt", "BMI")
colnames(cor.mult) <- c("HR", "SV", "CO", "VO2", "SBP", "DBP", "RPP", "MO", "TPR", "CI", "SI", "SDNN", "HRVTI",
                        "RMSDD", "BRS", "LF", "HF", "LF.HF", "IOP", "IJVP", "IJVF", "OPP", "PSV", "EDV", "CCA", "IJV", "Age", "Ht", "Wt", "BMI")

corrplot(cor.mult, method = "color")

network <- graph_from_adjacency_matrix(cor.mult, weighted=T, mode="undirected", diag=F) %>%
  set_vertex_attr("type", value = c('A','A','A','A','A','A','A','A','A','A','A',
                                    'B','B','B','B','B','B','B',
                                    'C','C','C','C','C','C','C','C',
                                    'D','D','D','D'))
network.mod <- network %>%
  set_edge_attr("strength", value = edge_attr(network, "weight")) %>%
  set_edge_attr("weight", value = abs(edge_attr(network, "weight"))) %>%
  set_vertex_attr("title", value = c("HR", "SV", "CO", "VO2", "SBP", "DBP", "RPP", "MO", "TPR", "CI", "SI",
                                     "SDNN", "HRVTi", "RMSDD", "BRS", "LF", "HF", "LF/HF",
                                     "IOP", "IJVP", "IJVF", "OPP", "PSV", "EDV", "A[CCA]", "A[IJV]",
                                     "Age", "Height", "Weight", "BMI"))

set.seed(9)
cor.graph <- ggraph(network.mod, layout = 'fr') +
  geom_edge_link(aes(color = strength, width = weight)) +
  geom_node_point(aes(fill = factor(type)), size = 12, color = "black", shape = 21) +
  geom_node_text(aes(label = title), size = 2) +
  scale_edge_width_continuous(name = element_blank(), breaks = c(), limits = c(0.21,1)) +
  scale_edge_colour_gradient2("Correlation",
                              low = "darkred", mid = "#FFFFFF", high = "darkblue",
                              limits = c(-1,1)) +
  scale_fill_manual("Type", values = c("A" = "lightgreen", "B" = "tan1","C" = "plum", "D" = "skyblue"),
                    labels = c("A" = "Systemic", "B" = "Autonomic","C" = "Head/Neck", "D" = "Characteristics")) +
  scale_x_continuous(expand = expansion(mult = 0.1)) +
  scale_y_continuous(expand = expansion(mult = 0.1)) +
  theme_void() + theme_graph(base_family="sans")
cor.graph

ggsave("cor.graph.png",cor.graph, width = 8, height = 8, units = "in")

##########################################
### Unused
##########################################
legend <- get_legend(data.hr)
fig.data.hemo <- ggarrange(data.hr,data.sv,data.co,data.vo2,
                           data.sbp,data.dbp,data.rpp,data.mo,
                           data.tpr,data.si,data.ci,
                           ncol = 4, nrow = 3, common.legend = TRUE, legend = "none",
                           labels = LETTERS[1:11], align = "hv") + draw_grob(legend, 17/24, -0.03, 1/3, 0.5)
ggsave("LBNP-fig_hemo.png", fig.data.hemo, width = 18*0.85, height = 16*0.85, units = "in")
fig.data.auto <- ggarrange(data.sdnn,data.hrvti,data.rmsdd,data.brs,
                           data.lfn,data.hfn,data.lfhf,
                           ncol = 4, nrow = 2, common.legend = TRUE, legend = "none",
                           labels = LETTERS[1:7], align = "hv") + draw_grob(legend, 17/24, 0, 1/3, 0.5)
ggsave("LBNP-fig_auto.png", fig.data.auto, width = 18*0.85, height = 16*0.85*2/3, units = "in")
legend2 <- get_legend(data.iop)
fig.data.ceph <- ggarrange(data.iop,data.opp,data.ijv,data.vpc,
                           data.cca,data.psv,data.edv,
                           ncol = 4, nrow = 2, common.legend = TRUE, legend = "none",
                           labels = LETTERS[1:7], align = "hv") + draw_grob(legend2, 17/24, 0, 1/3, 0.5)
ggsave("LBNP-fig_ceph.png", fig.data.ceph, width = 18*0.85, height = 16*0.85*2/3, units = "in")
ggsave("LBNP-fig_ijvf.png", data.ijvf, width = 8, height = 6, units = "in")


# Plot Model v1
fitted <- fitted(model.mult, nd2, re_formula = NA, probs = c(.055,.945))

# model.ijvf.dr <- cbind(nd2, fitted) %>%
#   select(Index:Face,`Estimate.P(Y = 1)`:`Q94.5.P(Y = 3)`) %>%
#   pivot_longer(cols = !Index:Face,
#                names_to = c(".value", "Pattern"),
#                names_sep = "P")
# 
# dose.ijvf <- ggplot(filter(model.ijvf.dr, Pattern == "(Y = 1)"), aes(x = Index, color = factor(Sex.c), fill = factor(Sex.c))) +
#   geom_smooth(aes(y = Estimate., ymin = Q5.5., ymax = Q94.5.),
#               stat = "identity",
#               alpha = 0.2) +
#   facet_grid(Side ~ Face, labeller = labeller(Face = face.labs, Side = side.labs), scales = "free") +
#   theme_classic2() +
#   theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
#   scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq, expand = c(0,0)) +
#   scale_y_continuous(name = "P(Grade 1 Flow)") +
#   scale_color_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                      labels = c("-0.5" = "Female", "0.5" = "Male"),
#                      breaks = c("0.5","-0.5")) +
#   scale_fill_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
#                     labels = c("-0.5" = "Female", "0.5" = "Male"),
#                     breaks = c("0.5","-0.5")) +
#   theme(legend.position = "bottom", panel.spacing = unit(1,"lines")) +
#   coord_cartesian(ylim = c(0,1), expand = FALSE)
# dose.ijvf

model.ijvf.dr2 <- cbind(nd2, fitted) %>%
  select(Index:Face,`Estimate.P(Y = 1)`:`Q94.5.P(Y = 3)`) %>%
  mutate(estimate = 1 - `Estimate.P(Y = 1)`,
         conf.high = 1 - `Q5.5.P(Y = 1)`,
         conf.low = 1 - `Q94.5.P(Y = 1)`)

dose.ijvf2 <- ggplot(model.ijvf.dr2, aes(x = Index, color = factor(Sex.c), fill = factor(Sex.c))) +
  geom_smooth(aes(y = estimate, ymin = conf.low, ymax = conf.high),
              stat = "identity",
              alpha = 0.2) +
  facet_grid(Side ~ Face, labeller = labeller(Face = face.labs, Side = side.labs), scales = "free") +
  theme_classic2() +
  theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
  scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq, expand = c(0,0)) +
  scale_y_continuous(name = "1 - P(Grade 1 Flow)") +
  scale_color_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
                     labels = c("-0.5" = "Female", "0.5" = "Male"),
                     breaks = c("0.5","-0.5")) +
  scale_fill_manual(name = "MAP ± 89% CrI:", values = c("-0.5" = "hotpink", "0.5" = "royalblue"),
                    labels = c("-0.5" = "Female", "0.5" = "Male"),
                    breaks = c("0.5","-0.5")) +
  theme(legend.position = "right", panel.spacing = unit(1,"lines")) +
  coord_cartesian(ylim = c(0,1), expand = FALSE)
dose.ijvf2

nd3 <- tidyr::crossing(Index = seq(from = 0, to = 5, length.out = 51),
                       Sex.c = 0,
                       Side = 0,
                       Face = c("L0","L15"))
fitted2 <- fitted(model.mult, nd3, re_formula = NA, probs = c(.055,.945))

model.ijvf.dr.v2 <- cbind(nd3, fitted2) %>%
  select(Index:Face,`Estimate.P(Y = 1)`:`Q94.5.P(Y = 3)`) %>%
  mutate(estimate = 1 - `Estimate.P(Y = 1)`,
         conf.high = 1 - `Q5.5.P(Y = 1)`,
         conf.low = 1 - `Q94.5.P(Y = 1)`)

nd4 <- tidyr::crossing(Index = seq(from = 0, to = 5, length.out = 6),
                       Sex.c = 0,
                       Side = 0,
                       Face = c("L0","L15"))
fitted2 <- fitted(model.mult, nd4, re_formula = NA, probs = c(.055,.945))

model.ijvf.dr.v2 <- cbind(nd4, fitted2) %>%
  select(Index:Face,`Estimate.P(Y = 1)`:`Q94.5.P(Y = 3)`) %>%
  mutate(estimate = 1 - `Estimate.P(Y = 1)`,
         conf.high = 1 - `Q5.5.P(Y = 1)`,
         conf.low = 1 - `Q94.5.P(Y = 1)`) %>%
  select(Index, Face, estimate, conf.low, conf.high) %>% arrange(Face, Index) %>%
  mutate(estimate = round(100*estimate,1),
         conf.low = round(100*conf.low,1),
         conf.high = round(100*conf.high,1))

dose.ijvf2 <- ggplot(model.ijvf.dr.v2, aes(x = Index)) +
  # geom_smooth(aes(y = estimate, ymin = conf.low, ymax = conf.high),
  #             stat = "identity",
  #             alpha = 0.2, color = "red") +
  geom_line(aes(y = estimate), color = "green") +
  geom_ribbon(aes(y = estimate, ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "green") +
  facet_grid(~ Face, labeller = labeller(Face = face.labs), scales = "free") +
  theme_classic2() +
  theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
  scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq, expand = c(0,0)) +
  scale_y_continuous(name = "1 - P(Grade 1 Flow)") +
  theme(legend.position = "right", panel.spacing = unit(1,"lines")) +
  coord_cartesian(ylim = c(0,1), expand = FALSE)
dose.ijvf2
ggsave("LBNP-fig_dr_ijvf.png", dose.ijvf2, width = 6, height = 3, units = "in")








model.ijvp.dr <- cbind(nd4, fitted2) %>%
  select(Index:Face,Estimate.VPCs:Q94.5.VPCs) %>%
  mutate(estimate = (Estimate.VPCs*scale.vpc)+center.vpc,
         conf.high = (Q94.5.VPCs*scale.vpc)+center.vpc,
         conf.low = (Q5.5.VPCs*scale.vpc)+center.vpc) %>%
  select(Index, Face, estimate, conf.low, conf.high) %>% arrange(Face, Index)

dose.ijvp <- ggplot(model.ijvp.dr, aes(x = Index)) +
  # geom_smooth(aes(y = estimate, ymin = conf.low, ymax = conf.high),
  #             stat = "identity",
  #             alpha = 0.2, color = "red") +
  geom_line(aes(y = estimate), color = "green") +
  geom_ribbon(aes(y = estimate, ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "green") +
  facet_grid(~ Face, labeller = labeller(Face = face.labs), scales = "free") +
  theme_classic2() +
  theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
  scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq, expand = c(0,0)) +
  scale_y_continuous(name = "IJV Pressure (mmHg)") +
  theme(legend.position = "right", panel.spacing = unit(1,"lines")) +
  coord_cartesian(ylim = c(0,50), expand = FALSE)
dose.ijvp
ggsave("LBNP-fig_dr_ijvp.png", dose.ijvp, width = 6, height = 3, units = "in")

model.ijv.dr <- cbind(nd4, fitted2) %>%
  select(Index:Face,Estimate.IJVs:Q94.5.IJVs) %>%
  mutate(estimate = (Estimate.IJVs*scale.ijv)+center.ijv,
         conf.high = (Q94.5.IJVs*scale.ijv)+center.ijv,
         conf.low = (Q5.5.IJVs*scale.ijv)+center.ijv) %>%
  select(Index, Face, estimate, conf.low, conf.high) %>% arrange(Face, Index)

dose.ijv <- ggplot(model.ijv.dr, aes(x = Index)) +
  # geom_smooth(aes(y = estimate, ymin = conf.low, ymax = conf.high),
  #             stat = "identity",
  #             alpha = 0.2, color = "red") +
  geom_line(aes(y = estimate), color = "green") +
  geom_ribbon(aes(y = estimate, ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "green") +
  facet_grid(~ Face, labeller = labeller(Face = face.labs), scales = "free") +
  theme_classic2() +
  theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
  scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq, expand = c(0,0)) +
  scale_y_continuous(name = expression(paste("IJV Area (",mm^2,")",sep=""))) +
  theme(legend.position = "right", panel.spacing = unit(1,"lines")) +
  coord_cartesian(ylim = c(0,175), expand = FALSE)
dose.ijv
ggsave("LBNP-fig_dr_ijv.png", dose.ijv, width = 6, height = 3, units = "in")


##########################################
# Back on Original Scale
##########################################
get_variables(model.mult)
# Describe Posterior
describe_posterior(model.mult,
                   centrality = c("mean","MAP"), dispersion = TRUE,
                   ci = 0.89, ci_method = "hdi",
                   test = c("p_direction", "rope"), rope_range = c(-0.05, 0.05), rope_ci = 1,
                   diagnostic = NULL, effects = "fixed", component = "all")


# Rescale HR
rescale.hr <- model.mult %>% spread_draws(`b_HRs_.*`, sigma_HRs, sd_ID__HRs_Intercept, nu_HRs, regex = TRUE) %>%
  mutate(b_HRs_Intercept = (b_HRs_Intercept * scale.hr) + center.hr,
         b_HRs_Index = b_HRs_Index * scale.hr,
         b_HRs_Sex.c = b_HRs_Sex.c * scale.hr,
         b_HRs_FaceL15 = b_HRs_FaceL15 * scale.hr,
         sigma_HRs = sigma_HRs * scale.hr,
         sd_ID__HRs_Intercept = sd_ID__HRs_Intercept * scale.hr)
describe_posterior(rescale.hr[,-c(1:3)],
                   centrality = c("MAP"), dispersion = TRUE,
                   ci = 0.89, ci_method = "hdi",
                   test = c("p_direction", "rope"), rope_range = c(-0.1*sd.hr, 0.1*sd.hr), rope_ci = 1,
                   diagnostic = NULL)

# Rescale SV
rescale.sv <- model.mult %>% spread_draws(`b_SVs_.*`, sigma_SVs, sd_ID__SVs_Intercept, nu_SVs, regex = TRUE) %>%
  mutate(b_SVs_Intercept = (b_SVs_Intercept * scale.sv) + center.sv,
         b_SVs_Index = b_SVs_Index * scale.sv,
         b_SVs_Sex.c = b_SVs_Sex.c * scale.sv,
         b_SVs_FaceL15 = b_SVs_FaceL15 * scale.sv,
         `b_SVs_Index:Sex.c` = `b_SVs_Index:Sex.c` * scale.sv,
         sigma_SVs = sigma_SVs * scale.sv,
         sd_ID__SVs_Intercept = sd_ID__SVs_Intercept * scale.sv)
describe_posterior(rescale.sv[,-c(1:3)],
                   centrality = c("mean","MAP"), dispersion = TRUE,
                   ci = 0.89, ci_method = "hdi",
                   test = c("p_direction", "rope"), rope_range = c(-0.1*sd.sv, 0.1*sd.sv), rope_ci = 1,
                   diagnostic = NULL)

# Rescale CO
rescale.co <- model.mult %>% spread_draws(`b_COs_.*`, sigma_COs, sd_ID__COs_Intercept, nu_COs, regex = TRUE) %>%
  mutate(b_COs_Intercept = (b_COs_Intercept * scale.co) + center.co,
         b_COs_Index = b_COs_Index * scale.co,
         b_COs_Sex.c = b_COs_Sex.c * scale.co,
         b_COs_FaceL15 = b_COs_FaceL15 * scale.co,
         sigma_COs = sigma_COs * scale.co,
         sd_ID__COs_Intercept = sd_ID__COs_Intercept * scale.co)
describe_posterior(rescale.co[,-c(1:3)],
                   centrality = c("MAP"),
                   ci = 0.89, ci_method = "hdi",
                   test = c("p_direction", "rope"), rope_range = c(-0.1*sd.co, 0.1*sd.co), rope_ci = 1,
                   diagnostic = NULL)

# Rescale VO2
rescale.vo2 <- model.mult %>% spread_draws(`b_VO2s_.*`, sigma_VO2s, sd_ID__VO2s_Intercept, nu_VO2s, regex = TRUE) %>%
  mutate(b_VO2s_Intercept = (b_VO2s_Intercept * scale.vo2) + center.vo2,
         b_VO2s_Index = b_VO2s_Index * scale.vo2,
         b_VO2s_Sex.c = b_VO2s_Sex.c * scale.vo2,
         b_VO2s_FaceL15 = b_VO2s_FaceL15 * scale.vo2,
         `b_VO2s_Sex.c:FaceL15` = `b_VO2s_Sex.c:FaceL15` * scale.vo2,
         sigma_VO2s = sigma_VO2s * scale.vo2,
         sd_ID__VO2s_Intercept = sd_ID__VO2s_Intercept * scale.vo2)
describe_posterior(rescale.vo2[,-c(1:3)],
                   centrality = c("MAP"),
                   ci = 0.89, ci_method = "hdi",
                   test = c("p_direction", "rope"), rope_range = c(-0.1*sd.vo2, 0.1*sd.vo2), rope_ci = 1,
                   diagnostic = NULL)

# Rescale SBP
rescale.sbp <- model.mult %>% spread_draws(`b_SBPs_.*`, sigma_SBPs, sd_ID__SBPs_Intercept, nu_SBPs, regex = TRUE) %>%
  mutate(b_SBPs_Intercept = (b_SBPs_Intercept * scale.sbp) + center.sbp,
         b_SBPs_Index = b_SBPs_Index * scale.sbp,
         b_SBPs_Sex.c = b_SBPs_Sex.c * scale.sbp,
         b_SBPs_FaceL15 = b_SBPs_FaceL15 * scale.sbp,
         sigma_SBPs = sigma_SBPs * scale.sbp,
         sd_ID__SBPs_Intercept = sd_ID__SBPs_Intercept * scale.sbp)
describe_posterior(rescale.sbp[,-c(1:3)],
                   centrality = c("MAP"),
                   ci = 0.89, ci_method = "hdi",
                   test = c("p_direction", "rope"), rope_range = c(-0.1*sd.sbp, 0.1*sd.sbp), rope_ci = 1,
                   diagnostic = NULL)

# Rescale DBP
rescale.dbp <- model.mult %>% spread_draws(`b_DBPs_.*`, sigma_DBPs, sd_ID__DBPs_Intercept, nu_DBPs, regex = TRUE) %>%
  mutate(b_DBPs_Intercept = (b_DBPs_Intercept * scale.dbp) + center.dbp,
         b_DBPs_Index = b_DBPs_Index * scale.dbp,
         b_DBPs_Sex.c = b_DBPs_Sex.c * scale.dbp,
         b_DBPs_FaceL15 = b_DBPs_FaceL15 * scale.dbp,
         sigma_DBPs = sigma_DBPs * scale.dbp,
         sd_ID__DBPs_Intercept = sd_ID__DBPs_Intercept * scale.dbp)
describe_posterior(rescale.dbp[,-c(1:3)],
                   centrality = c("MAP"),
                   ci = 0.89, ci_method = "hdi",
                   test = c("p_direction", "rope"), rope_range = c(-0.1*sd.dbp, 0.1*sd.dbp), rope_ci = 1,
                   diagnostic = NULL)

# Rescale RPP
rescale.rpp <- model.mult %>% spread_draws(`b_RPPs_.*`, sigma_RPPs, sd_ID__RPPs_Intercept, regex = TRUE) %>%
  mutate(b_RPPs_Intercept = (b_RPPs_Intercept * scale.rpp) + center.rpp,
         b_RPPs_Index = b_RPPs_Index * scale.rpp,
         b_RPPs_Sex.c = b_RPPs_Sex.c * scale.rpp,
         b_RPPs_FaceL15 = b_RPPs_FaceL15 * scale.rpp,
         sigma_RPPs = sigma_RPPs * scale.rpp,
         sd_ID__RPPs_Intercept = sd_ID__RPPs_Intercept * scale.rpp)
describe_posterior(rescale.rpp[,-c(1:3)],
                   centrality = c("MAP"),
                   ci = 0.89, ci_method = "hdi",
                   test = c("p_direction", "rope"), rope_range = c(-0.1*sd.rpp, 0.1*sd.rpp), rope_ci = 1,
                   diagnostic = NULL)

# Rescale MO
rescale.mo <- model.mult %>% spread_draws(`b_MOs_.*`, sigma_MOs, sd_ID__MOs_Intercept, nu_MOs, regex = TRUE) %>%
  mutate(b_MOs_Intercept = (b_MOs_Intercept * scale.mo) + center.mo,
         b_MOs_Index = b_MOs_Index * scale.mo,
         b_MOs_Sex.c = b_MOs_Sex.c * scale.mo,
         b_MOs_FaceL15 = b_MOs_FaceL15 * scale.mo,
         sigma_MOs = sigma_MOs * scale.mo,
         sd_ID__MOs_Intercept = sd_ID__MOs_Intercept * scale.mo)
describe_posterior(rescale.mo[,-c(1:3)],
                   centrality = c("MAP"),
                   ci = 0.89, ci_method = "hdi",
                   test = c("p_direction", "rope"), rope_range = c(-0.1*sd.mo, 0.1*sd.mo), rope_ci = 1,
                   diagnostic = NULL)

# Rescale TPR
rescale.tpr <- model.mult %>% spread_draws(`b_TPRs_.*`, sigma_TPRs, sd_ID__TPRs_Intercept, nu_TPRs, regex = TRUE) %>%
  mutate(b_TPRs_Intercept = (b_TPRs_Intercept * scale.tpr) + center.tpr,
         b_TPRs_Index = b_TPRs_Index * scale.tpr,
         b_TPRs_Sex.c = b_TPRs_Sex.c * scale.tpr,
         b_TPRs_FaceL15 = b_TPRs_FaceL15 * scale.tpr,
         sigma_TPRs = sigma_TPRs * scale.tpr,
         sd_ID__TPRs_Intercept = sd_ID__TPRs_Intercept * scale.tpr)
describe_posterior(rescale.tpr[,-c(1:3)],
                   centrality = c("MAP"),
                   ci = 0.89, ci_method = "hdi",
                   test = c("p_direction", "rope"), rope_range = c(-0.1*sd.tpr, 0.1*sd.tpr), rope_ci = 1,
                   diagnostic = NULL)

# Rescale SI
rescale.si <- model.mult %>% spread_draws(`b_SIs_.*`, sigma_SIs, sd_ID__SIs_Intercept, nu_SIs, regex = TRUE) %>%
  mutate(b_SIs_Intercept = (b_SIs_Intercept * scale.si) + center.si,
         b_SIs_Index = b_SIs_Index * scale.si,
         b_SIs_Sex.c = b_SIs_Sex.c * scale.si,
         b_SIs_FaceL15 = b_SIs_FaceL15 * scale.si,
         sigma_SIs = sigma_SIs * scale.si,
         sd_ID__SIs_Intercept = sd_ID__SIs_Intercept * scale.si)
describe_posterior(rescale.si[,-c(1:3)],
                   centrality = c("MAP"),
                   ci = 0.89, ci_method = "hdi",
                   test = c("p_direction", "rope"), rope_range = c(-0.1*sd.si, 0.1*sd.si), rope_ci = 1,
                   diagnostic = NULL)

# Rescale CI
rescale.ci <- model.mult %>% spread_draws(`b_CIs_.*`, sigma_CIs, sd_ID__CIs_Intercept, nu_CIs, regex = TRUE) %>%
  mutate(b_CIs_Intercept = (b_CIs_Intercept * scale.ci) + center.ci,
         b_CIs_Index = b_CIs_Index * scale.ci,
         b_CIs_Sex.c = b_CIs_Sex.c * scale.ci,
         b_CIs_FaceL15 = b_CIs_FaceL15 * scale.ci,
         sigma_CIs = sigma_CIs * scale.ci,
         sd_ID__CIs_Intercept = sd_ID__CIs_Intercept * scale.ci)
describe_posterior(rescale.ci[,-c(1:3)],
                   centrality = c("MAP"),
                   ci = 0.89, ci_method = "hdi",
                   test = c("p_direction", "rope"), rope_range = c(-0.1*sd.ci, 0.1*sd.ci), rope_ci = 1,
                   diagnostic = NULL)

# Rescale SDNN
rescale.sdnn <- model.mult %>% spread_draws(`b_SDNNs_.*`, sigma_SDNNs, sd_ID__SDNNs_Intercept, nu_SDNNs, regex = TRUE) %>%
  mutate(b_SDNNs_Intercept = (b_SDNNs_Intercept * scale.sdnn) + center.sdnn,
         b_SDNNs_Index = b_SDNNs_Index * scale.sdnn,
         b_SDNNs_Sex.c = b_SDNNs_Sex.c * scale.sdnn,
         b_SDNNs_FaceL15 = b_SDNNs_FaceL15 * scale.sdnn,
         sigma_SDNNs = sigma_SDNNs * scale.sdnn,
         sd_ID__SDNNs_Intercept = sd_ID__SDNNs_Intercept * scale.sdnn)
describe_posterior(rescale.sdnn[,-c(1:3)],
                   centrality = c("MAP"),
                   ci = 0.89, ci_method = "hdi",
                   test = c("p_direction", "rope"), rope_range = c(-0.1*sd.sdnn, 0.1*sd.sdnn), rope_ci = 1,
                   diagnostic = NULL)

# Rescale HRVTI
rescale.hrvti <- model.mult %>% spread_draws(`b_HRVTIs_.*`, sigma_HRVTIs, sd_ID__HRVTIs_Intercept, nu_HRVTIs, regex = TRUE) %>%
  mutate(b_HRVTIs_Intercept = (b_HRVTIs_Intercept * scale.hrvti) + center.hrvti,
         b_HRVTIs_Index = b_HRVTIs_Index * scale.hrvti,
         b_HRVTIs_Sex.c = b_HRVTIs_Sex.c * scale.hrvti,
         b_HRVTIs_FaceL15 = b_HRVTIs_FaceL15 * scale.hrvti,
         sigma_HRVTIs = sigma_HRVTIs * scale.hrvti,
         sd_ID__HRVTIs_Intercept = sd_ID__HRVTIs_Intercept * scale.hrvti)
describe_posterior(rescale.hrvti[,-c(1:3)],
                   centrality = c("MAP"),
                   ci = 0.89, ci_method = "hdi",
                   test = c("p_direction", "rope"), rope_range = c(-0.1*sd.hrvti, 0.1*sd.hrvti), rope_ci = 1,
                   diagnostic = NULL)

# Rescale RMSDD
rescale.rmsdd <- model.mult %>% spread_draws(`b_RMSDDs_.*`, sigma_RMSDDs, sd_ID__RMSDDs_Intercept, nu_RMSDDs, regex = TRUE) %>%
  mutate(b_RMSDDs_Intercept = (b_RMSDDs_Intercept * scale.rmsdd) + center.rmsdd,
         b_RMSDDs_Index = b_RMSDDs_Index * scale.rmsdd,
         b_RMSDDs_Sex.c = b_RMSDDs_Sex.c * scale.rmsdd,
         b_RMSDDs_FaceL15 = b_RMSDDs_FaceL15 * scale.rmsdd,
         sigma_RMSDDs = sigma_RMSDDs * scale.rmsdd,
         sd_ID__RMSDDs_Intercept = sd_ID__RMSDDs_Intercept * scale.rmsdd)
describe_posterior(rescale.rmsdd[,-c(1:3)],
                   centrality = c("MAP"),
                   ci = 0.89, ci_method = "hdi",
                   test = c("p_direction", "rope"), rope_range = c(-0.1*sd.rmsdd, 0.1*sd.rmsdd), rope_ci = 1,
                   diagnostic = NULL)

# Rescale BRS
rescale.brs <- model.mult %>% spread_draws(`b_BRSs_.*`, sigma_BRSs, sd_ID__BRSs_Intercept, nu_BRSs, regex = TRUE) %>%
  mutate(b_BRSs_Intercept = (b_BRSs_Intercept * scale.brs) + center.brs,
         b_BRSs_Index = b_BRSs_Index * scale.brs,
         b_BRSs_Sex.c = b_BRSs_Sex.c * scale.brs,
         b_BRSs_FaceL15 = b_BRSs_FaceL15 * scale.brs,
         sigma_BRSs = sigma_BRSs * scale.brs,
         sd_ID__BRSs_Intercept = sd_ID__BRSs_Intercept * scale.brs)
describe_posterior(rescale.brs[,-c(1:3)],
                   centrality = c("MAP"),
                   ci = 0.89, ci_method = "hdi",
                   test = c("p_direction", "rope"), rope_range = c(-0.1*sd.brs, 0.1*sd.brs), rope_ci = 1,
                   diagnostic = NULL)

# Rescale LFNorm
rescale.lfn <- model.mult %>% spread_draws(`b_LFNorms_.*`, sigma_LFNorms, sd_ID__LFNorms_Intercept, nu_LFNorms, regex = TRUE) %>%
  mutate(b_LFNorms_Intercept = (b_LFNorms_Intercept * scale.lfn) + center.lfn,
         b_LFNorms_Index = b_LFNorms_Index * scale.lfn,
         b_LFNorms_Sex.c = b_LFNorms_Sex.c * scale.lfn,
         b_LFNorms_FaceL15 = b_LFNorms_FaceL15 * scale.lfn,
         sigma_LFNorms = sigma_LFNorms * scale.lfn,
         sd_ID__LFNorms_Intercept = sd_ID__LFNorms_Intercept * scale.lfn)
describe_posterior(rescale.lfn[,-c(1:3)],
                   centrality = c("MAP"),
                   ci = 0.89, ci_method = "hdi",
                   test = c("p_direction", "rope"), rope_range = c(-0.1*sd.lfn, 0.1*sd.lfn), rope_ci = 1,
                   diagnostic = NULL)

# Rescale HFNorm
rescale.hfn <- model.mult %>% spread_draws(`b_HFNorms_.*`, sigma_HFNorms, sd_ID__HFNorms_Intercept, nu_HFNorms, regex = TRUE) %>%
  mutate(b_HFNorms_Intercept = (b_HFNorms_Intercept * scale.hfn) + center.hfn,
         b_HFNorms_Index = b_HFNorms_Index * scale.hfn,
         b_HFNorms_Sex.c = b_HFNorms_Sex.c * scale.hfn,
         b_HFNorms_FaceL15 = b_HFNorms_FaceL15 * scale.hfn,
         sigma_HFNorms = sigma_HFNorms * scale.hfn,
         sd_ID__HFNorms_Intercept = sd_ID__HFNorms_Intercept * scale.hfn)
describe_posterior(rescale.hfn[,-c(1:3)],
                   centrality = c("MAP"),
                   ci = 0.89, ci_method = "hdi",
                   test = c("p_direction", "rope"), rope_range = c(-0.1*sd.hfn, 0.1*sd.hfn), rope_ci = 1,
                   diagnostic = NULL)

# Rescale LFHF
rescale.lfhf <- model.mult %>% spread_draws(`b_LFHFs_.*`, sigma_LFHFs, sd_ID__LFHFs_Intercept, nu_LFHFs, regex = TRUE) %>%
  mutate(b_LFHFs_Intercept = (b_LFHFs_Intercept * scale.lfhf) + center.lfhf,
         b_LFHFs_Index = b_LFHFs_Index * scale.lfhf,
         b_LFHFs_Sex.c = b_LFHFs_Sex.c * scale.lfhf,
         b_LFHFs_FaceL15 = b_LFHFs_FaceL15 * scale.lfhf,
         `b_LFHFs_Index:Sex.c` = `b_LFHFs_Index:Sex.c` * scale.lfhf,
         sigma_LFHFs = sigma_LFHFs * scale.lfhf,
         sd_ID__LFHFs_Intercept = sd_ID__LFHFs_Intercept * scale.lfhf)
describe_posterior(rescale.lfhf[,-c(1:3)],
                   centrality = c("MAP"),
                   ci = 0.89, ci_method = "hdi",
                   test = c("p_direction", "rope"), rope_range = c(-0.1*sd.lfhf, 0.1*sd.lfhf), rope_ci = 1,
                   diagnostic = NULL)

# Rescale IOP
rescale.iop <- model.mult %>% spread_draws(`b_IOPs_.*`, sigma_IOPs, sd_ID__IOPs_Intercept, nu_IOPs, regex = TRUE) %>%
  mutate(b_IOPs_Intercept = (b_IOPs_Intercept * scale.iop) + center.iop,
         b_IOPs_Index = b_IOPs_Index * scale.iop,
         b_IOPs_Sex.c = b_IOPs_Sex.c * scale.iop,
         b_IOPs_FaceL15 = b_IOPs_FaceL15 * scale.iop,
         b_IOPs_Side = b_IOPs_Side * scale.iop,
         sigma_IOPs = sigma_IOPs * scale.iop,
         sd_ID__IOPs_Intercept = sd_ID__IOPs_Intercept * scale.iop)
describe_posterior(rescale.iop[,-c(1:3)],
                   centrality = c("MAP"),
                   ci = 0.89, ci_method = "hdi",
                   test = c("p_direction", "rope"), rope_range = c(-0.1*sd.iop, 0.1*sd.iop), rope_ci = 1,
                   diagnostic = NULL)

# Rescale OPP
rescale.opp <- model.mult %>% spread_draws(`b_OPPs_.*`, sigma_OPPs, sd_ID__OPPs_Intercept, nu_OPPs, regex = TRUE) %>%
  mutate(b_OPPs_Intercept = (b_OPPs_Intercept * scale.opp) + center.opp,
         b_OPPs_Index = b_OPPs_Index * scale.opp,
         b_OPPs_Sex.c = b_OPPs_Sex.c * scale.opp,
         b_OPPs_FaceL15 = b_OPPs_FaceL15 * scale.opp,
         b_OPPs_Side = b_OPPs_Side * scale.opp,
         sigma_OPPs = sigma_OPPs * scale.opp,
         sd_ID__OPPs_Intercept = sd_ID__OPPs_Intercept * scale.opp)
describe_posterior(rescale.opp[,-c(1:3)],
                   centrality = c("MAP"),
                   ci = 0.89, ci_method = "hdi",
                   test = c("p_direction", "rope"), rope_range = c(-0.1*sd.opp, 0.1*sd.opp), rope_ci = 1,
                   diagnostic = NULL)

# Rescale CCA
rescale.cca <- model.mult %>% spread_draws(`b_CCAs_.*`, sigma_CCAs, sd_ID__CCAs_Intercept, nu_CCAs, regex = TRUE) %>%
  mutate(b_CCAs_Intercept = (b_CCAs_Intercept * scale.cca) + center.cca,
         b_CCAs_Index = b_CCAs_Index * scale.cca,
         b_CCAs_Sex.c = b_CCAs_Sex.c * scale.cca,
         b_CCAs_FaceL15 = b_CCAs_FaceL15 * scale.cca,
         b_CCAs_Side = b_CCAs_Side * scale.cca,
         sigma_CCAs = sigma_CCAs * scale.cca,
         sd_ID__CCAs_Intercept = sd_ID__CCAs_Intercept * scale.cca)
describe_posterior(rescale.cca[,-c(1:3)],
                   centrality = c("MAP"),
                   ci = 0.89, ci_method = "hdi",
                   test = c("p_direction", "rope"), rope_range = c(-0.1*sd.cca, 0.1*sd.cca), rope_ci = 1,
                   diagnostic = NULL)

# Rescale IJV
rescale.ijv <- model.mult %>% spread_draws(`b_IJVs_.*`, sigma_IJVs, sd_ID__IJVs_Intercept, nu_IJVs, regex = TRUE) %>%
  mutate(b_IJVs_Intercept = (b_IJVs_Intercept * scale.ijv) + center.ijv,
         b_IJVs_Index = b_IJVs_Index * scale.ijv,
         b_IJVs_Sex.c = b_IJVs_Sex.c * scale.ijv,
         b_IJVs_FaceL15 = b_IJVs_FaceL15 * scale.ijv,
         b_IJVs_Side = b_IJVs_Side * scale.ijv,
         sigma_IJVs = sigma_IJVs * scale.ijv,
         sd_ID__IJVs_Intercept = sd_ID__IJVs_Intercept * scale.ijv)
describe_posterior(rescale.ijv[,-c(1:3)],
                   centrality = c("MAP"),
                   ci = 0.89, ci_method = "hdi",
                   test = c("p_direction", "rope"), rope_range = c(-0.1*sd.ijv, 0.1*sd.ijv), rope_ci = 1,
                   diagnostic = NULL)

# Rescale PSV
rescale.psv <- model.mult %>% spread_draws(`b_PSVs_.*`, sigma_PSVs, sd_ID__PSVs_Intercept, nu_PSVs, regex = TRUE) %>%
  mutate(b_PSVs_Intercept = (b_PSVs_Intercept * scale.psv) + center.psv,
         b_PSVs_Index = b_PSVs_Index * scale.psv,
         b_PSVs_Sex.c = b_PSVs_Sex.c * scale.psv,
         b_PSVs_FaceL15 = b_PSVs_FaceL15 * scale.psv,
         b_PSVs_Side = b_PSVs_Side * scale.psv,
         sigma_PSVs = sigma_PSVs * scale.psv,
         sd_ID__PSVs_Intercept = sd_ID__PSVs_Intercept * scale.psv)
describe_posterior(rescale.psv[,-c(1:3)],
                   centrality = c("MAP"),
                   ci = 0.89, ci_method = "hdi",
                   test = c("p_direction", "rope"), rope_range = c(-0.1*sd.psv, 0.1*sd.psv), rope_ci = 1,
                   diagnostic = NULL)

# Rescale EDV
rescale.edv <- model.mult %>% spread_draws(`b_EDVs_.*`, sigma_EDVs, sd_ID__EDVs_Intercept, nu_EDVs, regex = TRUE) %>%
  mutate(b_EDVs_Intercept = (b_EDVs_Intercept * scale.edv) + center.edv,
         b_EDVs_Index = b_EDVs_Index * scale.edv,
         b_EDVs_Sex.c = b_EDVs_Sex.c * scale.edv,
         b_EDVs_FaceL15 = b_EDVs_FaceL15 * scale.edv,
         b_EDVs_Side = b_EDVs_Side * scale.edv,
         `b_EDVs_Index:Sex.c` = `b_EDVs_Index:Sex.c` * scale.edv,
         `b_EDVs_Sex.c:FaceL15` = `b_EDVs_Sex.c:FaceL15` * scale.edv,
         `b_EDVs_Sex.c:Side` = `b_EDVs_Sex.c:Side` * scale.edv,
         sigma_EDVs = sigma_EDVs * scale.edv,
         sd_ID__EDVs_Intercept = sd_ID__EDVs_Intercept * scale.edv)
describe_posterior(rescale.edv[,-c(1:3)],
                   centrality = c("MAP"),
                   ci = 0.89, ci_method = "hdi",
                   test = c("p_direction", "rope"), rope_range = c(-0.1*sd.edv, 0.1*sd.edv), rope_ci = 1,
                   diagnostic = NULL)

# Rescale VPC
rescale.vpc <- model.mult %>% spread_draws(`b_VPCs_.*`, b_sigma_VPCs_Index, Intercept_sigma_VPCs, sd_ID__VPCs_Intercept, nu_VPCs, regex = TRUE) %>%
  mutate(b_VPCs_Intercept = (b_VPCs_Intercept * scale.vpc) + center.vpc,
         b_VPCs_Index = b_VPCs_Index * scale.vpc,
         b_VPCs_Sex.c = b_VPCs_Sex.c * scale.vpc,
         b_VPCs_FaceL15 = b_VPCs_FaceL15 * scale.vpc,
         b_VPCs_Side = b_VPCs_Side * scale.vpc,
         sd_ID__VPCs_Intercept = sd_ID__VPCs_Intercept * scale.vpc,
         Intercept_sigma_VPCs = exp(Intercept_sigma_VPCs + log(scale.vpc)))
describe_posterior(rescale.vpc[,-c(1:3)],
                   centrality = c("MAP"),
                   ci = 0.89, ci_method = "hdi",
                   test = c("p_direction", "rope"), rope_range = c(-0.1*sd.vpc, 0.1*sd.vpc), rope_ci = 1,
                   diagnostic = NULL)

get_variables(model.mult)
rescale.ijvf <- model.mult %>% spread_draws(`b_IJVF_.*`, sd_ID__IJVF_Intercept, regex = TRUE)
describe_posterior(rescale.ijvf[,-c(1:3)],
                   centrality = c("MAP"),
                   ci = 0.89, ci_method = "hdi",
                   test = c("p_direction", "rope"), rope_range = c(-0.1*pi/sqrt(3), 0.1*pi/sqrt(3)), rope_ci = 1,
                   diagnostic = NULL)


model.ijvf <- readRDS("BRMSFits/model.ijvf.rds")
describe_posterior(model.ijvf,
                   centrality = c("mean","MAP"), dispersion = TRUE,
                   ci = 0.89, ci_method = "hdi",
                   test = c("p_direction", "rope"), rope_range = c(-0.05, 0.05), rope_ci = 1,
                   diagnostic = NULL, effects = "fixed", component = "all")

get_variables(model.ijvf)
summary(model.ijvf)


# Rescale Posterior
rescale.vpc <- model.vpc %>% spread_draws(`b_.*`, `sd_.*`, nu, regex = TRUE) %>%
  mutate(b_Intercept = (b_Intercept * scale.vpc) + center.vpc,
         b_Index = b_Index * scale.vpc,
         b_Sex.c = b_Sex.c * scale.vpc,
         b_FaceL15 = b_FaceL15 * scale.vpc,
         b_Side = b_Side * scale.vpc,
         b_sigma_Intercept = exp(b_sigma_Intercept) * scale.vpc,
         b_sigma_Index = exp(b_sigma_Index) * scale.vpc,
         sd_ID__Intercept = sd_ID__Intercept * scale.vpc)
describe_posterior(rescale.vpc[,-c(1:3)],
                   centrality = c("mean","MAP"), dispersion = TRUE,
                   ci = 0.89, ci_method = "hdi",
                   test = c("p_direction", "rope"), rope_range = c(-0.1*sd.vpc, 0.1*sd.vpc), rope_ci = 1,
                   diagnostic = NULL)

get_variables(model.mult)
test <- model.mult %>% spread_draws(cor_ID__IOPs_Intercept__Weights_Intercept)
describe_posterior(test,
                   centrality = c("mean","MAP"), dispersion = TRUE,
                   ci = 0.89, ci_method = "hdi",
                   test = c("p_direction", "rope"), rope_range = c(-0.1, 0.1), rope_ci = 1,
                   diagnostic = NULL)

##################################
# Pressure Effect

pres.eff <- model.mult %>% spread_draws(`b_.*._Index`, regex = TRUE) %>%
  select(-c(.chain,.iteration,.draw,b_sigma_VPCs_Index)) %>%
  describe_posterior(centrality = c("MAP"),
                     ci = 0.89, ci_method = "hdi",
                     diagnostic = NULL) %>%
  select(Parameter,MAP,CI_low,CI_high) %>% mutate(MAP = MAP/10,
                                        CI_low = CI_low/10,
                                        CI_high = CI_high/10)

pos.eff <- model.mult %>% spread_draws(`b_.*._FaceL15`, regex = TRUE) %>%
  select(-c(.chain,.iteration,.draw)) %>%
  describe_posterior(centrality = c("MAP"),
                     ci = 0.89, ci_method = "hdi",
                     diagnostic = NULL) %>%
  select(Parameter,MAP,CI_low,CI_high)

eff <- cbind(pres.eff,pos.eff)
names(eff) <- c("Parameter", "MAP.LBNP", "CI_low.LBNP", "CI_high.LBNP","Parameter2", "MAP.Pos", "CI_low.Pos", "CI_high.Pos")
eff2 <- eff %>% select(-Parameter2) %>% mutate(LBNP = -MAP.Pos/MAP.LBNP,
                                               LBNP_high = -CI_high.Pos/CI_low.LBNP,
                                               LBNP_low = -CI_low.Pos/CI_high.LBNP) %>%
  mutate_if(is.numeric, round, 4) %>% mutate(LBNP = round(LBNP,1))

scaling <- cbind(center.hr,scale.hr,center.sv,scale.sv,center.co,scale.co,
              center.sbp,scale.sbp,center.dbp,scale.dbp,center.tpr,scale.tpr,
              center.vpc,scale.vpc)
write.csv(scaling,"LBNP_Scaling_Diss.csv")




nd <- tidyr::crossing(Index = seq(from = 0, to = 5, length.out = 6),
                      Sex.c = c(0),
                      Side = 0,
                      Face = c("L0","L15"))
fitted <- as.data.frame(fitted(model.mult, nd, re_formula = NA, probs = c(.055,.945)))

scaling <- read.csv("LBNP_Scaling_Diss.csv")

fitted.iop <- fitted %>% select(Estimate.IOPs,Q5.5.IOPs,Q94.5.IOPs) %>%
  mutate_all(function(x) ((x*scale.iop)+center.iop)) %>%
  cbind(nd)

fitted.opp <- fitted %>% select(Estimate.OPPs,Q5.5.OPPs,Q94.5.OPPs) %>%
  mutate_all(function(x) ((x*scale.opp)+center.opp)) %>%
  cbind(nd)



Pressure_seq <- c("0 mmHg","-10 mmHg","-20 mmHg","-30 mmHg","-40 mmHg","-50 mmHg")

lbnp.1 <- ggplot(data = fitted.iop, aes(x = Index, fill = Face)) +
  geom_line(aes(y = Estimate.IOPs, color = Face)) +
  geom_ribbon(aes(y = Estimate.IOPs, ymax = Q94.5.IOPs, ymin = Q5.5.IOPs, fill = Face), alpha = 0.2) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
  scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq, expand = c(0,0)) +
  scale_y_continuous(name = "IOP (mmHg)") +
   scale_color_manual(name = "Position:", values = c("L0" = "blue", "L15" = "red"),
                      labels = c("L0" = "0° Supine", "L15" = "15° HDT"),
                      breaks = c("L0","L15")) +
   scale_fill_manual(name = "Position:", values = c("L0" = "blue", "L15" = "red"),
                     labels = c("L0" = "0° Supine", "L15" = "15° HDT"),
                     breaks = c("L0","L15")) +
  theme(legend.position = "bottom", panel.spacing = unit(1,"lines")) +
  coord_cartesian(expand = FALSE)
lbnp.1

lbnp.2 <- ggplot(data = fitted.opp, aes(x = Index, fill = Face)) +
  geom_line(aes(y = Estimate.OPPs, color = Face)) +
  geom_ribbon(aes(y = Estimate.OPPs, ymax = Q94.5.OPPs, ymin = Q5.5.OPPs, fill = Face), alpha = 0.2) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle=60, hjust = 1, vjust = 1), legend.key.size = grid::unit(2, "lines")) +
  scale_x_continuous(name = element_blank(), breaks = c(seq(0,length(Pressure_seq)-1)), labels = Pressure_seq, expand = c(0,0)) +
  scale_y_continuous(name = "OPP (mmHg)") +
  scale_color_manual(name = "Position:", values = c("L0" = "blue", "L15" = "red"),
                     labels = c("L0" = "0° Supine", "L15" = "15° HDT"),
                     breaks = c("L0","L15")) +
  scale_fill_manual(name = "Position:", values = c("L0" = "blue", "L15" = "red"),
                    labels = c("L0" = "0° Supine", "L15" = "15° HDT"),
                    breaks = c("L0","L15")) +
  theme(legend.position = "bottom", panel.spacing = unit(1,"lines")) +
  coord_cartesian(expand = FALSE)
lbnp.2

figsans <- ggarrange(lbnp.1,lbnp.2,nrow = 1, ncol = 2, common.legend = TRUE, legend = "bottom", align = "hv", labels = LETTERS[1:2])
ggsave("SANS2.png", figsans, width = 8, height = 4, units = "in")

data.iop2 <- data.iop + theme(legend.position = "bottom")
legsample <- get_legend(data.iop2)

figsample <- ggarrange(data.co,data.rmsdd,data.iop, nrow = 1, ncol = 3, common.legend = TRUE, legend = "bottom", legend.grob = legsample,
                       align = "h", labels = LETTERS[1:3])
figsample
ggsave("LBNPSample.png",figsample,height = 4, width = 8, units = "in")

figsample2 <- ggarrange(data.ci,data.rmsdd,data.iop, nrow = 1, ncol = 3, common.legend = TRUE, legend = "bottom", legend.grob = legsample,
                       align = "h", labels = LETTERS[1:3])
ggsave("LBNPSample2.png",figsample2,height = 4, width = 8, units = "in")

#####################################
## Body Weight and IOP
#####################################

df.iop <- df %>%
  select(ID,Sex,Face,Index,Weight,Height,IOP.D,IOP.S) %>%
  mutate(Weight = case_when(ID == "omni007" ~ 73,
                            ID == "omni012" ~ 75,
                            ID == "omni112" ~ 76,
                            TRUE ~ Weight)) %>%
  mutate(Height = case_when(ID == "omni007" ~ 183,
                            ID == "omni012" ~ 173,
                            ID == "omni112" ~ 172,
                            TRUE ~ Height)) %>%
  #filter(Index == 0) %>%
  pivot_longer(cols = IOP.D:IOP.S, names_to = "Side", values_to = "IOP", names_prefix = "IOP.") %>%
  pivot_wider(values_from = IOP, names_from = Face, names_sep = ".") %>%
  mutate(dIOP = L15 - L0) %>%
  select(ID,Sex,Index,Side,Height,Weight,dIOP) %>%
  mutate(BMI = Weight/(Height/100)^2) %>%
  na.omit()
  
fig.diop <- ggplot(df.iop, aes(x = Weight,y=dIOP)) +
  geom_point(aes(shape = Side, color = Sex, fill = Index)) +
  scale_shape_manual("Side:", values = c("D"=21,"S"=23),labels = c("D"="Right","S"="Left")) +
  scale_color_manual("Sex:", values = c("Male"="royalblue","Female"="hotpink"), breaks = c("Male","Female")) +
  scale_fill_gradient("Pressure (mmHg):", low = "#000000", high = "#FFFFFF", breaks = c(0,1,2,3,4,5),
                      labels = c("0","-10","-20","-30","-40","-50")) +
  scale_x_continuous("Weight (kg)") +
  scale_y_continuous(name = expression(Delta*"IOP"[pos]*" (mmHg)")) +
  theme_classic2()
fig.diop

x.diop <- cbind(df.iop$Weight,df.iop$dIOP)
data.diop = list(x=x.diop, N=nrow(x.diop))

cor.diop = stan(file="robust_correlation.stan", data=data.diop, 
                 iter=20000, warmup=500, chains=4, seed=230502)

stan_trace(cor.diop, pars=c("rho", "mu", "sigma", "nu"))
stan_dens(cor.diop, pars=c("rho", "mu", "sigma", "nu"))
stan_plot(cor.diop, pars=c("rho", "mu", "sigma", "nu"))
print(cor.diop)

diop.rand = data.frame(extract(cor.diop, c("x_rand"))[[1]])

fig.diop2 <- fig.diop +
  stat_ellipse(data = diop.rand, aes(x = X1, y = X2), geom = "polygon", alpha = 0.2, fill = "green", color = "green", linetype = "dashed", level = 0.89) +
  stat_ellipse(data = diop.rand, aes(x = X1, y = X2), geom = "polygon", alpha = 0.2, fill = "green", color = "green", level = 0.5) +
  stat_ellipse(data = diop.rand, aes(x = X1, y = X2), geom = "polygon", alpha = 0.2, fill = "green", level = 0.95) +
  annotate(geom = "text", label = "50%", x = 84, y = 5) +
  annotate(geom = "text", label = "89%", x = 84, y = 7) +
  annotate(geom = "text", label = "95%", x = 84, y = 9)
fig.diop2

rho.diop <- as.numeric(extract(cor.diop, "rho")[[1]])
length(rho.diop)  # number of MCMC samples
median(rho.diop)    # posterior mean
HPDinterval(as.mcmc(rho.diop), prob=0.99)  # 99% highest posterior density interval
HPDinterval(as.mcmc(rho.diop), prob=0.89)  # 99% highest posterior density interval
mean(rho.diop)
mean(rho.diop <= 0)
mean(rho.diop >= 0)
mean(rho.diop > -0.1 & rho.diop < 0.1)
1 - mean(rho.diop > -0.1 & rho.diop < 0.1)


rho.diop.test <- describe_posterior(rho.diop,
                   centrality = "MAP",
                   ci = 0.89, ci_method = "hdi",
                   test = c("p_direction", "rope"), rope_range = c(-0.1, 0.1), rope_ci = 1)

color_scheme_set(scheme = "red")
fig.rho <- mcmc_areas(cor.diop, pars = c("rho"), prob = 0.89, area_method = "equal height") +
  theme_classic2() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(title = expression(rho))
fig.rho

fig.diop.fin <- ggarrange(fig.diop2,fig.rho,nrow = 2, ncol = 1, align = "hv",
                          common.legend = FALSE, labels = LETTERS[1:2], heights = c(5,2))
fig.diop.fin
ggsave("diop.png", fig.diop.fin, width = 6, height = 7, units = "in")

#ggsave("rho.pdf", fig.rho, width = 5, height = 2, units = "in")
#ggsave("diop.pdf", fig.diop2, width = 6, height = 5, units = "in")

#####################################
## BMI and IOP
#####################################

x.diop.bmi <- cbind(df.iop$BMI,df.iop$dIOP)
data.diop.bmi = list(x=x.diop.bmi, N=nrow(x.diop.bmi))

cor.diop.bmi = stan(file="robust_correlation.stan", data=data.diop.bmi, 
                iter=20000, warmup=500, chains=4, seed=230502)

rho.diop.bmi <- as.numeric(extract(cor.diop.bmi, "rho")[[1]])
rho.diop.bmi.test <- describe_posterior(rho.diop.bmi,
                                    centrality = "MAP",
                                    ci = 0.89, ci_method = "hdi",
                                    test = c("p_direction", "rope"), rope_range = c(-0.1, 0.1), rope_ci = 1)

#####################################
## Body Weight and OPP
#####################################

df.opp <- df %>%
  select(ID,Sex,Face,Index,Weight,Height,OPP.D,OPP.S) %>%
  mutate(Weight = case_when(ID == "omni007" ~ 73,
                            ID == "omni012" ~ 75,
                            ID == "omni112" ~ 76,
                            TRUE ~ Weight)) %>%
  mutate(Height = case_when(ID == "omni007" ~ 183,
                            ID == "omni012" ~ 173,
                            ID == "omni112" ~ 172,
                            TRUE ~ Height)) %>%
  #filter(Index == 0) %>%
  pivot_longer(cols = OPP.D:OPP.S, names_to = "Side", values_to = "OPP", names_prefix = "OPP.") %>%
  pivot_wider(values_from = OPP, names_from = Face, names_sep = ".") %>%
  mutate(dOPP = L15 - L0) %>%
  select(ID,Sex,Index,Side,Height,Weight,dOPP) %>%
  mutate(BMI = Weight/(Height/100)^2) %>%
  na.omit()

x.dopp <- cbind(df.opp$Weight,df.opp$dOPP)
data.dopp = list(x=x.dopp, N=nrow(x.dopp))

cor.dopp = stan(file="robust_correlation.stan", data=data.dopp, 
                    iter=20000, warmup=500, chains=4, seed=230502)

rho.dopp <- as.numeric(extract(cor.dopp, "rho")[[1]])
rho.dopp.test <- describe_posterior(rho.dopp,
                                        centrality = "MAP",
                                        ci = 0.89, ci_method = "hdi",
                                        test = c("p_direction", "rope"), rope_range = c(-0.1, 0.1), rope_ci = 1)

#####################################
## BMI and IOP
#####################################

x.dopp.bmi <- cbind(df.opp$BMI,df.opp$dOPP)
data.dopp.bmi = list(x=x.dopp.bmi, N=nrow(x.dopp.bmi))

cor.dopp.bmi = stan(file="robust_correlation.stan", data=data.dopp.bmi, 
                    iter=20000, warmup=500, chains=4, seed=230502)

rho.dopp.bmi <- as.numeric(extract(cor.dopp.bmi, "rho")[[1]])
rho.dopp.bmi.test <- describe_posterior(rho.dopp.bmi,
                                        centrality = "MAP",
                                        ci = 0.89, ci_method = "hdi",
                                        test = c("p_direction", "rope"), rope_range = c(-0.1, 0.1), rope_ci = 1)


























x.diop.bmi <- cbind(df.iop$BMI,df.iop$dIOP)
data.diop.bmi = list(x=x.diop.bmi, N=nrow(x.diop.bmi))

cor.diop.bmi = stan(file="robust_correlation.stan", data=data.diop.bmi, 
                iter=20000, warmup=500, chains=4, seed=230502)

stan_trace(cor.diop.bmi, pars=c("rho", "mu", "sigma", "nu"))
stan_dens(cor.diop.bmi, pars=c("rho", "mu", "sigma", "nu"))
stan_plot(cor.diop.bmi, pars=c("rho", "mu", "sigma", "nu"))
print(cor.diop.bmi)

diop.rand = data.frame(extract(cor.diop, c("x_rand"))[[1]])



rho <- data.frame(x = rho.diop)
ggplot(rho, aes(x = x)) + stat_halfeye(.width = c(0.89,0.95), color = "firebrick", fill = "red")


library(rstan)    # to run the Bayesian model (stan)
library(coda)     # to obtain HPD intervals (HPDinterval)
library(mvtnorm)  # to generate random correlated data (rmvnorm)
library(car) 

sigma = c(20, 40)
rho = -0.95
cov.mat = matrix(c(sigma[1] ^ 2,
                   sigma[1] * sigma[2] * rho,
                   sigma[1] * sigma[2] * rho,
                   sigma[2] ^ 2),
                 nrow=2, byrow=T)

set.seed(210191)
x.clean = rmvnorm(n=40, sigma=cov.mat)
plot(x.clean, pch=16)

data.clean = list(x=x.clean, N=nrow(x.clean))

cor.clean = stan(file="robust_correlation.stan", data=data.clean, 
                 iter=2000, warmup=500, chains=4, seed=210191)

stan_trace(cor.clean, pars=c("rho", "mu", "sigma", "nu"))
stan_dens(cor.clean, pars=c("rho", "mu", "sigma", "nu"))
stan_plot(cor.clean, pars=c("rho", "mu", "sigma", "nu"))
print(cor.clean)

x.rand = extract(cor.clean, c("x_rand"))[[1]]
plot(x.clean, xlim=c(-60, 55), ylim=c(-120, 120), pch=16)
dataEllipse(x.rand, levels = c(0.5, 0.95),
            fill=T, plot.points = FALSE)

ggplot(df.iop, aes(x = Weight, y = dIOP, color = Sex)) + geom_point()

cor.test(df.iop$Weight,df.iop$dIOP)

x.diop <- cbind(df.iop$Weight,df.iop$dIOP)
plot(x.diop, pch=16)

data.diop = list(x=x.diop, N=nrow(x.diop))
cor.diop = stan(file="robust_correlation.stan", data=data.diop, 
                 iter=2000, warmup=500, chains=4, seed=230502)
stan_trace(cor.diop, pars=c("rho", "mu", "sigma", "nu"))
stan_dens(cor.diop, pars=c("rho", "mu", "sigma", "nu"))
stan_plot(cor.diop, pars=c("rho", "mu", "sigma", "nu"))
print(cor.diop)

diop.rand = extract(cor.diop, c("x_rand"))[[1]]
plot(x.diop, xlim=c(40, 130), ylim=c(-5, 10), pch=16)
dataEllipse(diop.rand, levels = c(0.5, 0.95),
            fill=T, plot.points = FALSE)

rho.diop <- as.numeric(extract(cor.diop, "rho")[[1]])
length(rho.diop)  # number of MCMC samples
median(rho.diop)    # posterior mean
HPDinterval(as.mcmc(rho.diop), prob=0.99)  # 99% highest posterior density interval
HPDinterval(as.mcmc(rho.diop), prob=0.89)  # 99% highest posterior density interval
mean(rho.diop)
mean(rho.diop <= 0)
mean(rho.diop >= 0)
mean(rho.diop > -0.1 & rho.diop < 0.1)

source("rob.cor.mcmc.R")
cor.noisy2 = rob.cor.mcmc(x.diop)

x.dopp <- cbind(df.iop$Weight,df.iop$dOPP)
cor.dopp = rob.cor.mcmc(x.dopp)

test.pres <- model.mult %>% gather_draws(`b_.*._Index`, regex = TRUE) %>%
  mutate(.variable = str_remove_all(.variable,"b_|s_Index")) %>%
  mutate(type = case_when(.variable %in% c("HR","SV","CO","VO2","SBP","DBP","RPP","MO","TPR","CI","SI") ~ "A",
                          .variable %in% c("SDNN","HRVTI","RMSDD","BRS","LFNorm","HFNorm","LFHF") ~ "B",
                          .variable %in% c("IOP","VPC","OPP","PSV","EDV","CCA","IJV") ~ "C",
                          TRUE ~ "D")) %>%
  filter(.variable != ("sigma_VPC")) %>% filter(.variable != ("IJVF_Index")) %>%
  ungroup() %>%
  mutate(.variable = reorder(.variable, .value))


#######################################
### Redeuce -6° to 0

##################################
# Pressure Effect

pres.eff <- model.mult %>% spread_draws(`b_.*._Index`, regex = TRUE) %>%
  select(-c(.chain,.iteration,.draw,b_sigma_VPCs_Index)) %>%
  describe_posterior(centrality = c("MAP"),
                     ci = 0.89, ci_method = "hdi",
                     diagnostic = NULL) %>%
  select(Parameter,MAP,CI_low,CI_high) %>% mutate(MAP = MAP/10,
                                                  CI_low = CI_low/10,
                                                  CI_high = CI_high/10)

pos.eff <- model.mult %>% spread_draws(`b_.*._FaceL15`, regex = TRUE) %>%
  select(-c(.chain,.iteration,.draw)) %>%
  describe_posterior(centrality = c("MAP"),
                     ci = 0.89, ci_method = "hdi",
                     diagnostic = NULL) %>%
  select(Parameter,MAP,CI_low,CI_high) %>% mutate(MAP = MAP*6/15,
                                              CI_low = CI_low*6/15,
                                              CI_high = CI_high*6/15)

eff <- cbind(pres.eff,pos.eff)
names(eff) <- c("Parameter", "MAP.LBNP", "CI_low.LBNP", "CI_high.LBNP","Parameter2", "MAP.Pos", "CI_low.Pos", "CI_high.Pos")
eff3 <- eff %>% select(-Parameter2) %>% mutate(LBNP = -MAP.Pos/MAP.LBNP,
                                               LBNP_high = -CI_high.Pos/CI_low.LBNP,
                                               LBNP_low = -CI_low.Pos/CI_high.LBNP) %>%
  mutate_if(is.numeric, round, 4) %>% mutate(LBNP = round(LBNP,1))

