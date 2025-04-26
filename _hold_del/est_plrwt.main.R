# setup -----
  source("P:/ORD_Mcconeghy_202404036D/rzv/src/r/set_start.R")
  
# Execute Run ----
  i_control = list(boots = 200,
                   seed = as.integer(ymd('2025-1-20')),
                   outcome = 'HZ',
                   ncores=10,
                   outcvars = c("age_6574", "age_7584", "age_85p", 
                                "bsl_Anticoag", "bsl_can", "bsl_Cancer", "bsl_CHF", "bsl_DMany", 
                                "bsl_hzhx", "bsl_hzvlive", "bsl_HIV", "bsl_Hosp_90", 
                                "bsl_Insulin", "bsl_LiverSevere", "bsl_LiverMild",
                                "bsl_Mets", 
                                "bsl_Pulmonary", "bsl_Rheumatic", 
                                "bsl_Smoker"),
                   censvars = c("age_6574", "age_7584", "age_85p", "clc_resident",
                                "bsl_Anticoag", "tv_Anticoag", 
                                "bsl_Cancer", "tv_Cancer", "bsl_can", "tv_can", "can_miss",
                                'bsl_CHF', 'tv_CHF', 
                                "bsl_DMany", "tv_DMany",
                                'bsl_hzhx', 'bsl_hzvlive', 
                                "bsl_HIV", "tv_HIV", 
                                "bsl_Homeless", "bsl_Hosp_90", "tv_Hosp_90", 
                                "bsl_Insulin", "tv_Insulin", 
                                "bsl_LiverSevere", 
                                'bsl_MI', 'tv_MI', "bsl_Mets", "tv_Mets", 
                                "bsl_Renal", 
                                "bsl_Smoker", "tv_Smoker", "bsl_Stroke", 
                                "bsl_Pulmonary", "tv_Pulmonary",
                                "Hosps_CLC_90"),
                   i_grace = l_tte_params$i_grace,
                   i_int = l_tte_params$i_int,
                   i_fup = l_tte_params$i_fup,
                   runtime = Sys.time())

# Run ----
  i_run = f_tstmp()
  i_runpath = here('out', 'plrwt', paste0(i_run))
  dir.create(i_runpath, recursive=T)
  
  set.seed(i_control$seed)
  message(paste0("follow-up: ", i_control$i_fup, " days"))
  message(paste0("grace period: ", i_control$i_grace, " days"))
  message(paste0("Interval of follow-up: ", i_control$i_int, " days"))
  
  source(here('src', 'r', 'ipw', 'est_plrwt_01_dta.R'))
  
  ## Outcome: Any HZ infection  ----
    d_res = list(params = i_control)
  
    # build survival data, estimate weights
      source(here('src', 'r', 'ipw', 'est_plrwt_02_survdta.R'))
      source(here('src', 'r', 'ipw', 'est_plrwt_03_ipw.R'))
  
    ### HZ Point estimates ----
      source(here('src', 'r', 'ipw', 'est_plrwt_04_outc.R'))
  
    d_res$pe %>%
      ggplot(data = ., aes(x = time, group = model)) +
      geom_line(aes(y = cir, color = model, linetype = model), linewidth = 0.8) +
      scale_x_continuous(limits = c(0, 1440), breaks = c(0, 360, 720, 1080, 1440)) +
      labs(x = 'Time (days)', y = 'CIR')
    
    ### HZ Bootstraps ----
      d_res$boots = list()

       for (i in 1:i_control$boots) {
        cat('.')
        set.seed(i_control$seed+i)
        
        # For Poisson bootstrap
        d_wt = mutate(distinct(d_cl_pnl, PatientICN), freqwt = rpois(n(), 1L))
        
        source(here('src', 'r', 'ipw', 'est_plrwt_02_survdta.R'))
        
        d_cens = left_join(d_cens, d_wt, by='PatientICN') 
        d_outc = left_join(d_outc, d_wt, by='PatientICN') 
        
        source(here('src', 'r', 'ipw', 'est_plrwt_03_ipw.R'))
        source(here('src', 'r', 'ipw', 'est_plrwt_04_outc.R'))
        
        d_res$boots[[paste0(i)]] = surv_est
      }
  
    ### Save ----
      saveRDS(d_res, file = paste0(i_runpath, '//', 'res_plrwt_', i_control$outcome, '.Rds'))
  
  ## Outcome: Dementia ----
    i_control$outcome = 'Dem'
    
    i_control$outcvars = c("bsl_Anticoag", "bsl_Cancer", "bsl_CHF", "bsl_DMany", 
                           "bsl_hzhx", "bsl_HIV", "bsl_Homeless", 
                           "bsl_Insulin", "bsl_Mets", 
                           "bsl_Pulmonary", "bsl_Rheumatic", 
                           "bsl_Stroke", "bsl_Smoker")
    
    i_control$censvars = c("age_6574", "age_7584", "age_85p", "clc_resident",
                           "bsl_Anticoag", "tv_Anticoag", 
                           "bsl_Cancer", "tv_Cancer", "bsl_can", "tv_can", "can_miss",
                           'bsl_CHF', 'tv_CHF', 
                           "bsl_DMany", "tv_DMany",
                           'bsl_hzhx', 'bsl_hzvlive', 
                           "bsl_HIV", "tv_HIV", 
                           "bsl_Homeless", "bsl_Hosp_90", "tv_Hosp_90", 
                           "bsl_Insulin", "tv_Insulin", 
                           'bsl_MI', 'tv_MI', "bsl_Mets", "tv_Mets", 
                           "bsl_Renal", "tv_Renal", 
                           "bsl_Smoker", "tv_Smoker", "bsl_Stroke", "tv_Stroke", 
                           "bsl_Pulmonary", "tv_Pulmonary",
                           "Hosps_CLC_90")
  
    d_res = list(params = i_control)
  
    # build survival data, estimate weights
      source(here('src', 'r', 'ipw', 'est_plrwt_02_survdta.R'))
      source(here('src', 'r', 'ipw', 'est_plrwt_03_ipw.R'))
    
    ### Dem PE ----
      source(here('src', 'r', 'ipw', 'est_plrwt_04_outc.R'))
  
    d_res$pe %>%
      ggplot(data = ., aes(x = time, group = model)) +
      geom_line(aes(y = cir, color = model, linetype = model), linewidth = 0.8) +
      scale_x_continuous(limits = c(0, 1440), breaks = c(0, 360, 720, 1080, 1440)) +
      labs(x = 'Time (days)', y = 'CIR')
    
    ### Dem Bootstraps ----
      d_res$boots = list()
      
      for (i in 1:i_control$boots) {
        cat('.')
        set.seed(i_control$seed+i)
        
        # For Poisson bootstrap
        d_wt = mutate(distinct(d_cl_pnl, PatientICN), freqwt = rpois(n(), 1L))
        
        source(here('src', 'r', 'ipw', 'est_plrwt_02_survdta.R'))
        
        d_cens = left_join(d_cens, d_wt, by='PatientICN') 
        d_outc = left_join(d_outc, d_wt, by='PatientICN') 
        
        source(here('src', 'r', 'ipw', 'est_plrwt_03_ipw.R'))
        source(here('src', 'r', 'ipw', 'est_plrwt_04_outc.R'))
        
        d_res$boots[[paste0(i)]] = surv_est
        }

    ### Save ----
      saveRDS(d_res, paste0(i_runpath, '//', 'res_plrwt_', i_control$outcome, '.Rds'))