# setup -----
  library(here())
  
  source(here('src', 'setup.R'), echo=F)

  #increase size for bootstrapping procedure
  options(future.globals.maxSize = 4000*1024^2) 

  grid.draw.ggsurvplot = function(x) {
    survminer:::print.ggsurvplot(x, newpage=F)
  }
  
  #analysis options stored in list and saved
  d_output = list(runtime = Sys.time(),
                  #params = l_tte_params,
                  runplan = list(boots = 50,
                                 seed = as.integer(ymd('2024-11-16'))
                                 ))
  
  set.seed(d_output$runplan$seed)
  
  d_panel_treat =  readRDS(here('dta', 'survdta_treat_panel.R'))
  
  d_panel_outcome =  readRDS(here('dta', 'survdta_cloned_panel.R'))
  
# Estimate Weights ----
  # ADD BASELINE AND TIME_VARYING COVARIATES
  d_glm_wt = glm(outcome ~ poly(time, 2, raw=T) + poly(X1, 2) + X2, data=d_panel_treat, family=binomial())

  # estimate pr(treat==1)
  d_panel_treat$pr_treat = d_glm_wt$fitted.values
  
  # Calculative cumulative probability of non-treatment (treatment-free survival)
  setDT(d_panel_treat)
  d_panel_treat[, cumpr_notreat := cumprod(1-pr_treat), by = .(id)]
  
  # only need probabilities to join to outcome dataset
  d_panel_treat = select(d_panel_treat, id, time, cumpr_notreat)
  
  # join probability to outcome dataset
  d_panel_outcome = left_join(d_panel_outcome, d_panel_treat, by=c('id', 'time'))
  
  # calculate IPW
  setDT(d_panel_outcome)
  
  # PROJECT SPECIFIC!!! MUST CONSIDER CAREFULLY
  d_panel_outcome[, ipw := fcase(
    #cumpr_0; cumulative probability of no vaccination at time t
    assign==0, 1 / cumpr_notreat, 
    assign==1 & time < 12, 1, # trt - cant censor prior to grace
    assign==1 & time == 12 & treat_time < 12, 1, # assign=1 & treat time < grace end, cant censor at grace
    assign==1 & time==12 & treat_time==12, 1 / (1-cumpr_notreat), # treated at grace, then weight
    assign==1 & time==12 & treat_time>12, 0, # assign=1, and not treated before grace, censor
    assign==1 & time > 12, 1 # assign=1, after grace period, assign all weights as 1
  )]
  
  # cumprod of assign
  d_panel_outcome$ipw[d_panel_outcome$assign==1, ] = cumprod(d_panel_outcome$ipw[d_panel_outcome$assign==1, ])
  
  summary(d_panel_outcome$ipw)

# COVARIATE BALANCE ACROSS TIME ----
