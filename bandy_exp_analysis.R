# Load library
library(ggplot2) # use ggplot2 method of plotting
library(cowplot) # organise figures
library(dplyr) # for data manipulation like pipes %>%
library(tidyverse) 
library(colorspace)
library(brms)
library(rstan)
library(performance) # compare models
library(rptR) # repeatability estimation by glme models

# Functions
mytheme <- function() {
  theme_bw() + 
    theme(panel.border = element_rect(fill = NA, colour = "black", size = 0.8), # Border around plotting area.
          panel.grid.major = element_blank(), # Major grid lines blank
          panel.grid.minor = element_blank(), # Minor grid lines blank
          axis.line = element_blank(), # axis line size
          axis.ticks = element_line(colour = "black", size = 0.8),
          axis.text = element_text(size = 10, colour = "black"), # axis text size
          axis.title = element_text(size = 10), #axis title size
          panel.background = element_rect(fill = "transparent"), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          legend.background = element_rect(fill = "transparent"), # get rid of legend bg
          legend.box.background = element_rect(fill = "transparent", color = NA)) # get rid of legend panel bg)
} # set up plot theme

# Set directory
setwd('C:/Users/chwu6540/Dropbox (Personal)/Banded behaviour study') # on Desktop
setwd('C:/Users/nicho/Dropbox (Personal)/Banded behaviour study') # on Laptop

# Load data
adult_dat <- read.csv("bandy_exp_data.csv") # adult vs juvenile (within and between species)
str(adult_dat)

# Clean data
adult_dat <- adult_dat %>%
  dplyr::mutate(prey_type = factor(experiment, levels = c(1:6), 
                                   labels = c(
                                     "Normal hatchling", 
                                     "Patternless hatchling", 
                                     "Patternless hatchling (Painted)",
                                     "Patternless hatchling (Random painted)",
                                     "Cornsnake hatchling", 
                                     "Cornsnake hatchling (Painted)")))

# Proportion of predation attempt
aggregate(pred_attempt_bin ~ prey_type, adult_dat, mean)

## REPEATIBILITY ANALYSIS ## -----------------------------------------------------------------------------------------------------
# Estimate repeatability of predation attempt within juvenile prey choice via Binomial GLMM approach
rpt_model <- rpt(pred_attempt_bin ~ prey_type + (1|adult_ID), 
                 grname   = "adult_ID",
                 data     = adult_dat, 
                 datatype = "Binary",
                 nboot    = 1000, 
                 npermut  = 1000,
                 parallel = TRUE)

# Estimate repeatability of time to attack within juvenile prey choice via Poisson GLMM approach
rpt_model_2 <- rpt(time_to_attack_s ~ prey_type + (1|adult_ID), 
                   grname   = "adult_ID",
                   data     = adult_dat, 
                   datatype = "Poisson",
                   nboot    = 1000, npermut = 1000,
                   parallel = TRUE)

print(rpt_model)
print(rpt_model_2)

## PREDATION ATTEMPT ANALYSIS ## ------------------------------------------------------------------------------------------
# General STAN specs
rstan::rstan_options(auto_write = TRUE) # translate to STAN platform for running Bayesian model
options(mc.cores = parallel::detectCores()) # detects how many cores avaliable to use

# Account for variation in individual prey choice (1|juvenile_ID) and individual variation in behaviour of the test subjects (1|adult_ID)
# No covariates
predAttempt_m1 <- brms::brm(pred_attempt ~ prey_type + (1|juvenile_ID) + (1|adult_ID), data = adult_dat, family = 'bernoulli', prior = set_prior('normal(0, 3)'), iter = 5e3, warmup = 2.5e3, chains = 4, cores = 4, control = list(adapt_delta = .99, max_treedepth = 15))

# Prey size as covariate (prey size may influence predation attempt)
predAttempt_m2 <- brms::brm(pred_attempt ~ prey_type + juvenile_svl_cm + (1|juvenile_ID) + (1|adult_ID), data = adult_dat, family = 'bernoulli', prior = set_prior('normal(0, 3)'), iter = 5e3, warmup = 2.5e3, chains = 4, cores = 4, control = list(adapt_delta = .99, max_treedepth = 15))

# Sex as covariate (sex difference in behaviour?)
predAttempt_m3 <- brms::brm(pred_attempt ~ prey_type + adult_sex + (1|juvenile_ID) + (1|adult_ID), data = adult_dat, family = 'bernoulli', prior = set_prior('normal(0, 3)'), iter = 5e3, warmup = 2.5e3, chains = 4, cores = 4, control = list(adapt_delta = .99, max_treedepth = 15))

# Sex and size as covariate
predAttempt_m4 <- brms::brm(pred_attempt ~ prey_type + juvenile_svl_cm + adult_sex + (1|juvenile_ID) + (1|adult_ID), data = adult_dat, family = 'bernoulli', prior = set_prior('normal(0, 3)'), iter = 5e3, warmup = 2.5e3, chains = 4, cores = 4, control = list(adapt_delta = .99, max_treedepth = 15))

# Model diagnostics 
performance::compare_performance(predAttempt_m1, predAttempt_m2, predAttempt_m3, predAttempt_m4, rank = TRUE)
predAttempt_model <- predAttempt_m4 # m4 was the best model

summary(predAttempt_model)
brms::fixef(predAttempt_model)
brms::pp_check(predAttempt_model)

## PLOT FIG 2a ##
# Predict from brm model
pred_newdata               <- data.frame(prey_type       = levels(adult_dat$prey_type),
                                         juvenile_svl_cm = mean(adult_dat$juvenile_svl_cm),
                                         adult_sex       = mean(adult_dat$adult_sex))
pred_fit                   <- brms::fitted(predAttempt_model,
                                           pred_newdata = pred_newdata,
                                           re_formula   = NA, # Ignores random effects
                                           summary      = TRUE # mean and 95% CI
                                           ) * 100 # convert to %
predAttempt_pred           <- cbind(pred_newdata, pred_fit)
predAttempt_pred$prey_type <- factor(predAttempt_pred$prey_type, 
                                     levels = c(
                                       "Normal hatchling", 
                                       "Patternless hatchling", 
                                       "Patternless hatchling (Painted)",
                                       "Patternless hatchling (Random painted)",
                                       "Cornsnake hatchling", 
                                       "Cornsnake hatchling (Painted)")
                                     ) # correct factor levels

predAttempt_plot <- predAttempt_pred %>%
  ggplot(aes(x = prey_type, y = Estimate, colour = prey_type)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), size = 1, width = 0.1) + 
  colorspace::scale_colour_discrete_sequential(palette = "Viridis") +
  xlab(NULL) + ylab("Probability of predation attempt (%)") +
  ylim(0,100) +
  mytheme()

# Contrasts between prey types
predAttempt_fit <- as.data.frame(fitted(
  predAttempt_model,
  newdata    = data.frame(prey_type    = levels(adult_dat$prey_type),
                       juvenile_svl_cm = mean(adult_dat$juvenile_svl_cm),
                       adult_sex       = unique(adult_dat$adult_sex)),
  re_formula = NA,
  summary    = FALSE # extract the full MCMC
  ))

colnames(predAttempt_fit) <- newdata$prey_type
# Painted vs random painted
paint_vs_patternless <- predAttempt_fit$"Patternless hatchling (Painted)" - predAttempt_fit$"Patternless hatchling (Random painted)"
quantile(paint_vs_patternless, probs = c(.5, .025, .975))

# Corn snake vs painted corn snake
paint_vs_patternless <- predAttempt_fit$"Cornsnake hatchling" - predAttempt_fit$"Cornsnake hatchling (Painted)"
quantile(paint_vs_patternless, probs = c(.5, .025, .975))

# Patternless vs painted
paint_vs_patternless <- predAttempt_fit$"Patternless hatchling" - predAttempt_fit$"Patternless hatchling (Painted)"
quantile(paint_vs_patternless, probs = c(.5, .025, .975))

## TIME TO STRIKE ANALYSIS ## ----------------------------------------------------------------------------------------------------
# Convert time attack to proportion
adult_dat$time_to_attack_s_proportion <- adult_dat$time_to_attack_s / 900

# Same logic as predation attempt analysis
attack_m1 <- brms::brm(time_to_attack_s_proportion ~ prey_type + (1|juvenile_ID) + (1|adult_ID), data = adult_dat, family = zero_one_inflated_beta(), iter = 5e3, warmup = 2.5e3, chains = 4, cores = 4, control = list(adapt_delta = .99, max_treedepth = 15)) 
attack_m2 <- brms::brm(time_to_attack_s_proportion ~ prey_type + juvenile_svl_cm + (1|juvenile_ID) + (1|adult_ID), data = adult_dat, family = zero_one_inflated_beta(), iter = 5e3, warmup = 2.5e3, chains = 4, cores = 4, control = list(adapt_delta = .99, max_treedepth = 15)) 
attack_m3 <- brms::brm(time_to_attack_s_proportion ~ prey_type + adult_sex + (1|juvenile_ID) + (1|adult_ID), data = adult_dat, family = zero_one_inflated_beta(), iter = 5e3, warmup = 2.5e3, chains = 4, cores = 4, control = list(adapt_delta = .99, max_treedepth = 15)) 
attack_m4 <- brms::brm(time_to_attack_s_proportion ~ prey_type + juvenile_svl_cm + adult_sex + (1|juvenile_ID) + (1|adult_ID), data = adult_dat, family = zero_one_inflated_beta(), iter = 5e3, warmup = 2.5e3, chains = 4, cores = 4, control = list(adapt_delta = .99, max_treedepth = 15)) 

# Model diagnostics
performance::compare_performance(attack_m1, attack_m2, attack_m3, attack_m4, rank = TRUE)
attack_model <- attack_m3 # m3 was the best model
summary(attack_model)
brms::pp_check(attack_model)
brms::fixef(attack_model) # fixed effects

## PLOT FIG 2b ##
source('./geom_flat_violin.R')

attack_plot <- adult_dat %>%
  ggplot(aes(x = prey_type, y = time_to_attack_s, fill = prey_type)) +
  geom_jitter(aes(colour = prey_type), position = position_jitter(0.05), size = 2) +
  geom_flat_violin(alpha = 0.5, position = position_nudge(x = 0.1)) + # from geom_flat_violin.R
  mytheme() +
  geom_hline(yintercept = 900, linetype = "dashed") + 
  colorspace::scale_fill_discrete_sequential(palette = "Viridis") +
  colorspace::scale_colour_discrete_sequential(palette = "Viridis") +
  xlab(NULL) + ylab("Time to first strike (s)") +
  scale_y_reverse()

# CREATE FIG 2
plot_grid(predAttempt_plot + theme(legend.position = "none"), 
          attack_plot + theme(legend.position = "none"), 
          align = "h", axis = "bt", nrow = 2, labels = "auto")

