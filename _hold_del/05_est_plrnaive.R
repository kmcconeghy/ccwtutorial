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
  
  d_panel =  readRDS(here('dta', 'survdta_cloned_panel.R'))
  setDT(d_panel)
  
# Naive estimator ----

# Run Plan ----
  
  # defined above
  set.seed(d_output$runplan$seed)
  
  d_glm = glm(event_outc ~ poly(time, 2, raw=T)*assign, data=d_panel, family=binomial())
  
  ## Survival probabilities ----
  d_panel$pr_ev = d_glm$fitted.values
  
  d_panel[, `:=`(pr_surv = cumprod(1 - pr_ev)), by=list(id, assign)] 
  
  d_res = d_panel %>%
    group_by(assign, time) %>%
    summarize(pr_ev = mean(1-pr_surv)) %>%
    ungroup %>%
    pivot_wider(., id_cols =c('time'), 
                names_from = assign, 
                names_prefix = 'pr_ev_',
                values_from = pr_ev
    ) %>%
    mutate(cid = pr_ev_1 - pr_ev_0,
           cir = pr_ev_1 / pr_ev_0)

  ### Plot ----
    d_gg_ci = d_res %>%
      ggplot(aes(x=time)) +
      geom_line(aes(y = pr_ev_0), color='red') +
      geom_line(aes(y = pr_ev_1), color='blue') + 
      scale_x_continuous(breaks = seq(0, 60, 6),
                         limits = c(0, 60)) +
      theme_bw() +
      labs(x = 'Follow-up (years)', y = 'Cumulative incidence')

    d_gg_rr = d_res %>%
      ggplot(aes(x=time)) +
      geom_line(aes(y = cir), color='green') + 
      scale_y_continuous(limits = c(0.5, 1.1)) + 
      scale_x_continuous(breaks = seq(0, 60, 6),
                         limits = c(0, 60)) +
      theme_bw() +
      labs(x = 'Follow-up (years)', y = 'Relative Risk')
    
    d_gg_1 = ggarrange(d_gg_ci, d_gg_rr,
                       nrow=1)
    
    ggsave( here('img', 
                 paste0('survplot_logmodunadj', '.jpeg')),
            plot=d_gg_1, width = 12, height=8, dpi=300)

# Save KM est results ----
  write_csv(d_res,
    file = here('out', paste0('plrnaiveoutc.csv'))
  )

d_output$results = d_res

saveRDS(d_output, 
        file = here('out', paste0('plrnaive.Rds')))