f_plr_cuminc = function(x) {
  
  d_glm = glm(event==0 ~ poly(time, 2, raw=T)*assign, data=x, family=binomial()) # <1>
  
  d_plr_naive_est = crossing(assign = 1:0, time = 1:60)
  
  d_plr_naive_est$pr_surv = predict(d_glm, newdata = d_plr_naive_est, type = 'response') 
  
  d_plr_naive_est %>%
    group_by(assign) %>%
    mutate(pr_cumsurv = cumprod(pr_surv),
           pr_cumev = 1 - pr_cumsurv) %>% # <2>
    ungroup %>% # <3>
    pivot_wider(., id_cols =c('time'), # <3>
                names_from = assign, # <3>
                names_prefix = 'pr_ev_', # <3>
                values_from = pr_cumev # <3>
    ) %>%
    mutate(cid = pr_ev_1 - pr_ev_0, # <3>
           cir = pr_ev_1 / pr_ev_0)
}