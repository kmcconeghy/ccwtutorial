f_cuminc = function(dta) {
  dta %>%
    group_by(assign, id) %>%
    mutate(pr_cumsurv = cumprod(pr_surv),
           pr_cumev = 1 - pr_cumsurv) %>% 
    ungroup %>% 
    group_by(assign, time) %>%
    summarize(pr_cumsurv = mean(pr_cumsurv),
              pr_cumev = 1 - mean(pr_cumsurv), .groups = 'drop') %>% 
    ungroup %>%
    pivot_wider(., id_cols =c('time'), 
                names_from = assign,
                names_prefix = 'pr_ev_', 
                values_from = pr_cumev 
    ) %>%
    mutate(cid = pr_ev_1 - pr_ev_0, 
           cir = pr_ev_1 / pr_ev_0)
}


f_dtmod = function(dta, tfun = 'poly(time, 2, raw=T)', i_covars = c(''), modstrat=F, ...) {
  
  #Stratify models by assignment or using interaction
    if (modstrat) {i_assign = ' '} else {i_assign='*assign'}
  
  # Build model formula
    modform = as.formula(paste0('event==0 ~ ', tfun, i_assign, i_covars))
  
    if (modstrat) {
      
      d_glm_0 = glm(modform, data=dta[assign==0, ], model=F, ...) 
      d_glm_1 = glm(modform, data=dta[assign==1, ], model=F, ...) 
      
      d_glm = list(d_glm_0, d_glm_1)
    } else {
      d_glm = glm(modform, data=dta, model=F, ...) # <1>
    }
  
  return(d_glm)
  
}