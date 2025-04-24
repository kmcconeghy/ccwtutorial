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
                  runplan = list(boots = 50, #bootstraps
                                 workers = 8,
                                 seed = as.integer(ymd('2024-11-16'))
                                 ))

# Clone ----

  d_cloned = readRDS(here('dta', 'survdta_cloned.R'))
  d_ids = distinct(d_cloned, id) # list of unique person IDs (for bootstrap)
  
# Naive estimator ----

# Run a simple Kaplan-Meier, No adjustments made ----
  ## Kaplan-Meier function ----
  # one row per clone, not panel data, needs time, event, treat variables

# Point estimate ----

  d_surv_pe = broom::tidy(
    survfit(Surv(t_clone, event_outc) ~ assign, data=d_cloned)) %>%
    mutate(assign = str_extract(strata,  '(?<=assign=)\\d')
    ) %>%
    arrange(time) %>%
    select(-strata) %>%
    mutate(pr_e = 1-estimate) %>%
    rename(pr_s = estimate) %>%
    pivot_wider(id_cols = c(time), 
                names_from = assign,
                values_from = c(pr_e, pr_s,
                                n.risk, n.event)) %>%
    mutate(cir = pr_e_1 / pr_e_0,
           cid = (pr_e_1 - pr_e_0))

  d_survmin = ggsurvplot(survfit(Surv(t_clone, event_outc) ~ assign, data=d_cloned),
                            data=data_cloned,
                            fun = 'event',
                            xlim = c(0, 60),
                            xlab = 'Follow-up',
                            break.time.by=6,
                            palette = c('red', 'blue'),
                            linetype = 'strata',
                            censor=F,
                            cumevents=T,
                            conf.int=F, pval=F,
                            risk.table=T, risk.table.col = 'strata',
                            risk.table.height=0.25, 
                            ggtheme = theme_classic())
  
  ggsave( here('img', 
               paste0('survplot_survmin', '.jpeg')),
          d_survmin, width = 12, height=7, dpi=300)
  
# defined above
plan(multisession, workers = d_output$runplan$workers)
set.seed(d_output$runplan$seed)

# bootstrap 
d_fun_getkmres = function(dta, d_ids, ...) {
  
  # For Poisson bootstrap - cluster by person
  # Generate a frequency weight for each person
  # Wt ~ Poisson(Lambda=1)
  d_ids = mutate(d_ids, freqwt = rpois(n(), 1L))
  
  dta_2 = left_join(dta, d_ids, by='id') 
  
  fit_all = survfit(Surv(t_clone, event_outc) ~ assign, data=dta_2, weights = freqwt)
  
  # summarize results from KM estimator
  survdta = broom::tidy(fit_all) %>%
    mutate(assign = str_extract(strata,  '(?<=assign=)\\d')
    ) %>%
    arrange(time) %>%
    select(-strata) %>%
    mutate(pr_e = 1-estimate) %>% # cumulative incidence
    rename(pr_s = estimate) %>%
    pivot_wider(id_cols = c(time), 
                names_from = assign,
                values_from = c(pr_e, pr_s,
                                n.risk, n.event)) %>%
    mutate(cir = pr_e_1 / pr_e_0, # relative risk
           cid = (pr_e_1 - pr_e_0)) # risk difference
  
  return(survdta)
}

  d_bs = future_map(.x = 1:d_output$runplan$boots, 
                       .f = ~d_fun_getkmres(d_cloned, d_ids, .x),
                       .options = furrr_options(seed = T))

# Combine and summarize ----
    d_surv = d_bs %>%
      bind_rows(.id = 'boot') %>%
      mutate(boot = as.numeric(boot)) 
    
    d_summ_surv = d_surv %>%
      group_by(time) %>%
      summarize(pr_e_0_lc = quantile(pr_e_0, 0.025, na.rm=T),
                pr_e_1_lc = quantile(pr_e_1, 0.025, na.rm=T),
                cir_lc = quantile(cir, 0.025, na.rm=T),
                cid_lc = quantile(cid, 0.025, na.rm=T),
                pr_e_0_uc = quantile(pr_e_0, 0.975, na.rm=T),
                pr_e_1_uc = quantile(pr_e_1, 0.975, na.rm=T),
                cir_uc = quantile(cir, 0.975, na.rm=T),
                cid_uc = quantile(cid, 0.975, na.rm=T)) %>%
      bind_cols(select(d_surv_pe, pr_e_0, pr_e_1, cir, cid),
                .)

  ### Plot ----
    d_gg_ci = d_summ_surv %>%
      ggplot(aes(x=time)) +
      geom_line(aes(y = pr_e_0), color='red') +
      geom_line(aes(y = pr_e_1), color='blue') +
      geom_ribbon(aes(ymin = pr_e_0_lc, ymax = pr_e_0_uc), 
                  fill='red', alpha=0.2) + 
      geom_ribbon(aes(ymin = pr_e_1_lc, ymax = pr_e_1_uc), 
                  fill='blue', alpha=0.2) + 
      scale_x_continuous(breaks = seq(0, 60, 6),
                         limits = c(0, 60)) +
      theme_bw() +
      labs(x = 'Follow-up', y = 'Cumulative incidence')

    d_gg_rr = d_summ_surv %>%
      ggplot(aes(x=time)) +
      geom_line(aes(y = cir), color='green') +
      geom_hline(aes(yintercept=1), linetype=2) +
      geom_ribbon(aes(ymin = cir_lc, ymax = cir_uc), 
                  fill='green', alpha=0.2) + 
      scale_y_continuous(limits = c(0.5, 1.1)) +
      scale_x_continuous(breaks = seq(0, 60, 6),
                         limits = c(0, 60)) +
      theme_bw() +
      labs(x = 'Follow-up', y = 'Relative Risk')
    
    d_gg_1 = ggarrange(d_gg_ci, d_gg_rr,
                       nrow=1)
    
    ggsave( here('img', 
                 paste0('survplot_kmunadj', '.jpeg')),
            plot=d_gg_1, width = 12, height=8, dpi=300)

# Save KM est results ----

  write_csv(d_summ_surv,
    file = here('out', paste0('kmnoadjoutc', '.csv'
    ))
  )

d_output$results$summ = d_summ_surv
d_output$results$boot = d_surv

saveRDS(d_output, 
        file = here('out', paste0('kmnaive.Rds')))