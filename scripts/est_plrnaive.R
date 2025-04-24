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
                  params = l_tte_params,
                  runplan = list(boots = 200,
                                 seed = as.integer(ymd('2024-11-16'))
                                 ))
  
  d_elig =  readRDS(here('dta', f_get_file(i_dtalist, 'elig')))
  
  d_trial = readRDS(here('dta', f_get_file(i_dtalist, 'trial')))

  source(here('src', 'r', 'chk_trial.R'))
  source(here('src', 'r', 'bld_pnl.R'))
  source(here('src', 'r', 'bld_clone.R'))

# Naive estimator ----

  message(paste0("follow-up: ", d_output$params$i_fup, " days"))
  message(paste0("grace period: ", d_output$params$i_grace, " days"))
  message(paste0("Interval of follow-up: ", d_output$params$i_int, " days"))

# Run Plan ----
  
  # defined above
  set.seed(d_output$runplan$seed)
  
# bootstrap 
  d_bs_hz = f_est_plr(d_cl_pnl, 
                      boots=d_output$runplan$boots, 
                      outcome = 'HZ', ncores=10, adjvars = c(''))
  
  d_bs_dem = f_est_plr(d_cl_pnl, 
                       boots=d_output$runplan$boots,
                       outcome = 'Dem', ncores=10, adjvars = c('')
                       )

# Combine and summarize ----
  ## Herpes Zoster ----
    
    d_summ_surv_hz = d_bs_hz$boots %>%
      group_by(time) %>%
      summarize(pr_e_0_lc = quantile(pr_e_0, 0.025, na.rm=T),
                pr_e_1_lc = quantile(pr_e_1, 0.025, na.rm=T),
                cir_lc = quantile(cir, 0.025, na.rm=T),
                cid_lc = quantile(cid, 0.025, na.rm=T),
                pr_e_0_uc = quantile(pr_e_0, 0.975, na.rm=T),
                pr_e_1_uc = quantile(pr_e_1, 0.975, na.rm=T),
                cir_uc = quantile(cir, 0.975, na.rm=T),
                cid_uc = quantile(cid, 0.975, na.rm=T)) %>%
      bind_cols(select(d_bs_hz$estimate, pr_e_0, pr_e_1, cir, cid),
                .)

  ### Plot ----
    d_gg_ci = d_summ_surv_hz %>%
      ggplot(aes(x=time)) +
      geom_line(aes(y = pr_e_0), color='red') +
      geom_line(aes(y = pr_e_1), color='blue') +
      geom_ribbon(aes(ymin = pr_e_0_lc, ymax = pr_e_0_uc), 
                  fill='red', alpha=0.2) + 
      geom_ribbon(aes(ymin = pr_e_1_lc, ymax = pr_e_1_uc), 
                  fill='blue', alpha=0.2) + 
      scale_x_continuous(breaks = seq(0, 365*4, 365),
                         limits = c(0, l_tte_params$i_fup),
                         labels = floor(seq(0, 365*4, 365)/365)) +
      theme_bw() +
      labs(x = 'Follow-up (years)', y = 'Cumulative incidence')

    d_gg_rr = d_summ_surv_hz %>%
      ggplot(aes(x=time)) +
      geom_line(aes(y = cir), color='green') +
      geom_ribbon(aes(ymin = cir_lc, ymax = cir_uc), 
                  fill='green', alpha=0.2) + 
      scale_y_continuous(limits = c(0.5, 1.1)) +
      scale_x_continuous(breaks = seq(0, 365*4, 365),
                         limits = c(0, l_tte_params$i_fup),
                         labels = floor(seq(0, 365*4, 365)/365)) +
      theme_bw() +
      labs(x = 'Follow-up (years)', y = 'Relative Risk')
    
    d_gg_1 = ggarrange(d_gg_ci, d_gg_rr,
                       nrow=1)
    
    ggsave( here('img', 
                 paste0(prj.specs$prj.prefix, '.', 'survplot_logmodunadj_hz', '.jpeg')),
            plot=d_gg_1, width = 12, height=8, dpi=300)

d_summ_surv_dem = d_bs_dem$boots %>%
  group_by(time) %>%
  summarize(pr_e_0_lc = quantile(pr_e_0, 0.025, na.rm=T),
            pr_e_1_lc = quantile(pr_e_1, 0.025, na.rm=T),
            cir_lc = quantile(cir, 0.025, na.rm=T),
            cid_lc = quantile(cid, 0.025, na.rm=T),
            pr_e_0_uc = quantile(pr_e_0, 0.975, na.rm=T),
            pr_e_1_uc = quantile(pr_e_1, 0.975, na.rm=T),
            cir_uc = quantile(cir, 0.975, na.rm=T),
            cid_uc = quantile(cid, 0.975, na.rm=T)) %>%
  bind_cols(select(d_bs_dem$estimate, pr_e_0, pr_e_1, cir, cid),
            .)

### Plot ----
d_gg_ci = d_summ_surv_dem %>%
  ggplot(aes(x=time)) +
  geom_line(aes(y = pr_e_0), color='red') +
  geom_line(aes(y = pr_e_1), color='blue') +
  geom_ribbon(aes(ymin = pr_e_0_lc, ymax = pr_e_0_uc), 
              fill='red', alpha=0.2) + 
  geom_ribbon(aes(ymin = pr_e_1_lc, ymax = pr_e_1_uc), 
              fill='blue', alpha=0.2) + 
  scale_y_continuous() +
  scale_x_continuous(breaks = seq(0, 365*4, 365),
                     limits = c(0, l_tte_params$i_fup),
                     labels = floor(seq(0, 365*4, 365)/365)) +
  theme_bw() +
  labs(x = 'Follow-up (years)', y = 'Cumulative incidence')

d_gg_rr = d_summ_surv_dem %>%
  ggplot(aes(x=time)) +
  geom_line(aes(y = cir), color='green', linewidth=1.2) +
  geom_ribbon(aes(ymin = cir_lc, ymax = cir_uc), 
              fill='green', alpha=0.2) + 
  scale_y_continuous(limits = c(0.75, 1.25)) +
  scale_x_continuous(breaks = seq(0, 365*4, 365),
                     limits = c(0, l_tte_params$i_fup),
                     labels = floor(seq(0, 365*4, 365)/365)) +
  theme_bw() +
  labs(x = 'Follow-up (years)', y = 'Relative Risk')

d_gg_1 = ggarrange(d_gg_ci, d_gg_rr,
                   nrow=1)

ggsave( here('img', paste0(prj.specs$prj.prefix, '.', 
                              'survplot_logmodunadj_dem', '.jpeg')),
        plot=d_gg_1,
        width = 12, height=8, dpi=300)

# Save KM est results ----
  d_summ_surv = bind_rows(
    mutate(d_summ_surv_hz, outcome = 'hz'),
    mutate(d_summ_surv_dem, outcome = 'dem'))

  write_csv(d_summ_surv,
    file = here('out', paste0(prj.specs$prj.prefix, '.', 
                                 'plrnaiveoutc', '.csv'
    ))
  )

d_output$results$hz = d_bs_hz 
d_output$results$dem = d_bs_dem 

saveRDS(d_output, 
        file = here('out', paste0('plrnaive', '.', f_tstmp(), '.Rds')))