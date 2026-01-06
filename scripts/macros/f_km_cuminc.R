f_km_cuminc = function(x) {
  broom::tidy(survfit(Surv(time, event) ~ assign, data=x)) %>% 
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
           cid = (pr_ev_1 - pr_ev_0))
}