est_f_plr = function(x, control) {

  d_return = list(params = control)

  # Data with outcome = vaccination
    dta_ipw = f_d_logit(x, censor=T, control)

  # estimate pr(vacc)
    d_return$mod_cens = f_logmod(dta_ipw[dta_ipw$assign==0, ], censor=T, control = control) 

  # Calculate cum. prob (vacc=0)
    dta_ipw = dta_ipw[dta_ipw$assign==0, ]
    setDT(dta_ipw)
    dta_ipw$pr_1 = d_return$mod_cens$pr_y
    dta_ipw[, cumpr_0 := cumprod(1-pr_1), by = .(assign, PatientICN)]
    dta_ipw = select(dta_ipw, PatientICN, time, cumpr_0)

  # Make outcome dataset
    dta = f_d_logit(d_cl_pnl, censor=F, control = control)

  # join probability to outcome dataset
    dta = left_join(dta, dta_ipw, by=c('PatientICN', 'time'))

  # calculate IPW
    setDT(dta)
    
    dta[, ipw := fcase(
      #censor_art = censor for nonadherence, i.e. vaccination status
      #cumpr_0; cumulative probability of no vaccination at time t
      assign==0 & censor_art==0, 1 / cumpr_0, 
      assign==0 & censor_art==1, 0, 
      assign==1 & cumdays < control$i_grace, 1, # trt - cant censor prior to grace
      # assign=treat & vacc prior to grace, cant censor
      assign==1 & cumdays == control$i_grace & t_hzv_rec_1 <= (control$i_grace - control$i_int), 1,
      # assign=treat, and treated on grace, 1 / 1 - cumpr_0
      assign==1 & cumdays == control$i_grace & 
        (t_hzv_rec_1 <= control$i_grace & t_hzv_rec_1 > (control$i_grace - control$i_int)), 1 / (1-cumpr_0),
      # assign=treat, and not treated on or before grace, censor
      assign==1 & cumdays == control$i_grace & (t_hzv_rec_1 > control$i_grace), 0,
      # assign=treat, after grace period, assign all weights as 1
      assign==1 & cumdays > control$i_grace, 1
    )]

  setDF(dta)

d_return$data = dta

# cycle through models; no adjustment, no weights, weights and outcome adjustments  
  control2 = control
  control2$outcvars = c()
  
  d_return$mod_outc_noipwnoadj = f_logmod(d_return$data, censor=F, ipw=F, control = control2)
  d_return$mod_outc_noipwnoadj$summ_est = f_summ_est(d_return$data, 
                                                     fitted.values = d_return$mod_outc_noipwnoadj$pr_y)
  
  d_return$mod_outc_noipwadj = f_logmod(d_return$data, censor=F, ipw=F, control = control)
  d_return$mod_outc_noipwadj$summ_est = f_summ_est(d_return$data, 
                                                   fitted.values = d_return$mod_outc_noipwadj$pr_y)
  
  d_return$mod_outc_ipwnoadj = f_logmod(d_return$data, censor=F, ipw=T, control = control2)
  d_return$mod_outc_ipwnoadj$summ_est = f_summ_est(d_return$data, 
                                                   fitted.values = d_return$mod_outc_ipwnoadj$pr_y)
  
  d_return$mod_outc_ipw = f_logmod(d_return$data, censor=F, ipw=T, control = control)
  d_return$mod_outc_ipw$summ_est = f_summ_est(d_return$data, 
                                              fitted.values = d_return$mod_outc_ipw$pr_y)
  
  if (!is.null(control$boots)) {
    message('beginning bootstrap procedure')
    
    summboots = map(.x = 1:control$boots, 
                    .f = \(boots) {
                      # For Poisson bootstrap
                      d_wt = mutate(distinct(x, PatientICN), freqwt = rpois(n(), 1L))
                      
                      dta_wt = left_join(x, d_wt, by='PatientICN') 
                      
                      # Data with outcome = vaccination
                      dta_ipw = f_d_logit(dta_wt, censor=T, control)
                      
                      # estimate pr(vacc)
                      mod_cens = f_logmod(dta_ipw[dta_ipw$assign==0, ], censor=T, control = control) 
                      
                      # Calculate cum. prob (vacc=0)
                      dta_ipw = dta_ipw[dta_ipw$assign==0, ]
                      setDT(dta_ipw)
                      dta_ipw$pr_1 = mod_cens$pr_y
                      dta_ipw[, cumpr_0 := cumprod(1-pr_1), by = .(assign, PatientICN)]
                      dta_ipw = select(dta_ipw, PatientICN, time, cumpr_0, freqwt)
                      
                      # Make outcome dataset
                      dta = f_d_logit(dta_wt, censor=F, control = control)
                      
                      # join probability to outcome dataset
                      dta = left_join(dta, dta_ipw, by=c('PatientICN', 'time'))
                      
                      # calculate IPW
                      setDT(dta)
                      
                      dta[, ipw := fcase(
                        #censor_art = censor for nonadherence, i.e. vaccination status
                        #cumpr_0; cumulative probability of no vaccination at time t
                        assign==0 & censor_art==0, 1 / cumpr_0, 
                        assign==0 & censor_art==1, 0, 
                        assign==1 & cumdays < control$i_grace, 1, # trt - cant censor prior to grace
                        # assign=treat & vacc prior to grace, cant censor
                        assign==1 & cumdays == control$i_grace & t_hzv_rec_1 <= (control$i_grace - control$i_int), 1,
                        # assign=treat, and treated on grace, 1 / 1 - cumpr_0
                        assign==1 & cumdays == control$i_grace & 
                          (t_hzv_rec_1 <= control$i_grace & t_hzv_rec_1 > (control$i_grace - control$i_int)), 1 / (1-cumpr_0),
                        # assign=treat, and not treated on or before grace, censor
                        assign==1 & cumdays == control$i_grace & (t_hzv_rec_1 > control$i_grace), 0,
                        # assign=treat, after grace period, assign all weights as 1
                        assign==1 & cumdays > control$i_grace, 1
                      )]
                      
                      setDF(dta)
                      
                      mod_outc_noipwnoadj = f_logmod(dta, censor=F, ipw=F, control = control2)
                      est_1 = f_summ_est(dta, fitted.values = mod_outc_noipwnoadj$pr_y) %>%
                        mutate(outcome = control$outcome, model = 'No IPW, No OM adj.')
                      
                      mod_outc_noipwadj = f_logmod(dta, censor=F, ipw=F, control = control)
                      est_2 = f_summ_est(dta, fitted.values = mod_outc_noipwadj$pr_y) %>%
                        mutate(outcome = control$outcome, model = 'No IPW, OM adj.')
                      
                      mod_outc_ipwnoadj = f_logmod(dta, censor=F, ipw=T, control = control2)
                      est_3 = f_summ_est(dta, fitted.values = mod_outc_ipwnoadj$pr_y) %>%
                        mutate(outcome = control$outcome, model = 'IPW, No OM adj.')
                      
                      mod_outc_ipw = f_logmod(dta, censor=F, ipw=T, control = control)
                      est_4 = f_summ_est(dta, fitted.values = mod_outc_ipw$pr_y) %>%
                        mutate(outcome = control$outcome, model = 'IPW, OM adj.')
                      
                      message('.')
                      
                      return(bind_rows(est_1, est_2, est_3, est_4))
                      
                    }
    ) %>%
      bind_rows(.id = 'boot') %>%
      mutate(boot = as.numeric(boot)) 
    
    d_return$boots = summboots
  }
  
return(d_return)

}






