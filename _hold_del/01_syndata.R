library(here())

source(here('src', 'setup.R'), echo=F)

# Set seed for reproducibility
set.seed(42)

# Number of participants
n <- 1000

# COVARIATES ----

  # Generate age (normal distribution with mean=75, std=10)
    age <- round(rnorm(n, mean = 75, sd = 10), 1)
  
  # Generate gender (66% female, 34% male)
    female <- sample(c(1, 0), size = n, replace = TRUE, prob = c(0.66, 0.34))
  
  # Generate treatment status (20% treated, 90% control)
    treat <- sample(c(1, 0), size = n, replace = TRUE, prob = c(0.2, 0.9))

# SIMULATE FAILURE TIMES ----
  
  # Simulate failure times (exponential distribution for simplicity)
    hr_treatment = 0.8
    hr_gender = 0.9
    hr_age = 0.98
  
  # Calculate linear predictor (log-hazard)
    lp_outc <- log(hr_treatment) * treat +
      log(hr_gender) * female +
      log(hr_age) * age
  
  # Simulate baseline survival times (exponential distribution)
    baseline_hazard <- 0.05 
    baseline_survival <- rexp(n, rate = baseline_hazard)
    
  # Adjust survival times based on linear predictor
    t_outc <- round(baseline_survival * exp(-lp_outc))
    summary(t_outc)
  
  # TREATIME INITIATION TIMES
    lp_treat <- log(0.7) * female +
      log(1.02) * age
  
    baseline_hazard <- 0.1
    t_treat <- ifelse(treat==1, round(rexp(n, rate = baseline_hazard)* exp(-lp_treat)), Inf)
    
    summary(t_treat)

  # Create dataset
    data <- data.frame(
      X1 = age,
      X2 = female,
      treat = treat, #ever treated
      t_treat = t_treat, # time treat initiated
      t_outc = t_outc # time of outcome
    ) %>%
      mutate(id = row_number(),
             time = pmin(60, t_outc), # set 60 as end of follow-up
             event_outc = if_else(time==t_outc, 1L, 0L))

# quick test
survfit(Surv(time, event_outc) ~ treat, data = data)
  
saveRDS(data, here('dta', 'survdta.R'))