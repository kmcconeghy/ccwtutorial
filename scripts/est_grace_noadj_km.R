
# Kaplan-Meier Estimator ----

# 1.  Naive estimator, K-M estimator, assignment is by clone (not person). `t_clone` is the follow-up time for each clone, taking into account artificial censoring.
# 2.  `estimate` from the model is the survival probability, so probability of event is 1 - estimate
# 3.  Data is in long form, so one colum per group with cumulative incidence/survival
# 4.  Estimands, `cir` is ratio analogous to relative risk, and `cid` is analogous to risk difference


  dta_c_person = readRDS(here('dta', 'dta_cloned_person.Rds'))
  
  d_mods = dta_c_person %>%
    nest(data = -model) %>%
    mutate(est_km = map(data, 
                        ~broom::tidy(
                          survfit(Surv(time, event) ~ assign, data=.)) %>% 
                          mutate(assign = str_extract(strata,  '(?<=assign=)\\d')
                          ) %>%
                          arrange(time) %>%
                          select(-strata) %>%
                          mutate(pr_ev = 1-estimate) %>% 
                          rename(pr_s = estimate) %>%
                          pivot_wider(id_cols = c(time),
                                      names_from = assign, 
                                      values_from = c(pr_ev, pr_s, 
                                                      n.risk, n.event)) %>% 
                          mutate(cir = pr_ev_1 / pr_ev_0, 
                                 cid = (pr_ev_1 - pr_ev_0)))) 
  
  d_gg = d_mods %>%
    mutate(km_gg = map(data, ~ggsurvplot(survfit(Surv(time, event) ~ assign, data=.),
                                         data=dta_c_person,
                                         fun = 'event',
                                         xlim = c(1, 60),
                                         xlab = 'Follow-up',
                                         break.time.by=6,
                                         palette = cbbPalette,
                                         linetype = 'strata',
                                         censor=F,
                                         cumevents=T,
                                         conf.int=F, pval=F,
                                         risk.table=T, risk.table.col = 'strata',
                                         risk.table.height=0.25, 
                                         ggtheme = theme_bw())))
  
  d_gg$km_gg[[1]] # A-> Y
  
  d_gg$km_gg[[2]] # X -> Y
  
  d_res = list(estimate = d_mods,
               graphs = list(ay = d_gg$km_gg[[1]], xy = d_gg$km_gg[[2]])
               )
  
  saveRDS(d_res, here('dta', 'est_noadj_km.Rds'))