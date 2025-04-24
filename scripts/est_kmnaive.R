# setup -----
  source("P:/ORD_Mcconeghy_202404036D/rzv/src/r/set_start.R")
  
  library(future)
  library(progressr)
  library(ggpubr)
  library(survival)
  library(survminer)
  library(furrr)

  #increase size for bootstrapping procedure
  options(future.globals.maxSize = 4000*1024^2) 

  grid.draw.ggsurvplot = function(x) {
    survminer:::print.ggsurvplot(x, newpage=F)
  }
  
  #analysis options stored in list and saved
  d_output = list(runtime = Sys.time(),
                  params = l_tte_params,
                  runplan = list(boots = 500,
                                 workers = 20,
                                 seed = as.integer(ymd('2024-11-16'))
                                 ))
  d_elig =  readRDS(here('dta', f_get_file(i_dtalist, 'elig')))
  
  d_trial = readRDS(here('dta', f_get_file(i_dtalist, 'trial')))

  source(here('src', 'r', 'chk_trial.R'))

# Clone ----

  d_cl = bind_rows(mutate(d_trial, assign=1),
                   mutate(d_trial, assign=0)) %>%
    mutate(dt_start = dt_clcadmit,
           dt_end = if_else(!is.na(DeathDate), 
                            DeathDate, 
                            pmin(dt_clcadmit + d_output$params$i_fup, 
                                 max(dt_clcadmit))
           ),
           t_end = as.integer((dt_end - dt_start)+1))
  
# Naive estimator ----

  message(paste0("follow-up: ", d_output$params$i_fup, " days"))
  message(paste0("grace period: ", d_output$params$i_grace, " days"))
  message(paste0("Interval of follow-up: ", d_output$params$i_int, " days"))

# Run a simple Kaplan-Meier, No adjustments made ----
  ## Kaplan-Meier function ----
  # one row per clone, not panel data, needs time, event, treat variables

  d_fun_getkmres = function(dta, d_ids, ...) {
    
    # For Poisson bootstrap
    d_ids = mutate(d_ids, freqwt = rpois(n(), 1L))
    
    dta_2 = left_join(dta, d_ids, by='PatientICN') 
    
    fit_all = survfit(Surv(time, event) ~ treat, 
                      data=dta_2, weights = freqwt)
    
    survdta = broom::tidy(fit_all) %>%
      mutate(treat = str_extract(strata,  '(?<=treat=)\\d')
      ) %>%
      arrange(time) %>%
      select(-strata) %>%
      mutate(pr_e = 1-estimate) %>%
      rename(pr_s = estimate) %>%
      pivot_wider(id_cols = c(time), 
                  names_from = treat,
                  values_from = c(pr_e, pr_s,
                                  n.risk, n.event)) %>%
      mutate(cir = pr_e_1 / pr_e_0,
             cid = (pr_e_1 - pr_e_0))
    
    return(survdta)
  }

  d_cl_hz = d_cl %>%
    mutate(time = case_when(
      assign == 1 & t_hzv_rec_1 > l_tte_params$i_grace ~ 
        pmin(t_hzdx, t_death, t_end, l_tte_params$i_grace),
      assign == 1 & t_hzv_rec_1 <= l_tte_params$i_grace ~ 
        pmin(t_hzdx, t_death, t_end),
      assign == 0 ~ pmin(t_hzv_rec_1, t_hzdx, 
                         t_death, t_end)),
      event = case_when(
        assign==1 & t_hzdx == time ~ 1L, 
        assign==1 & t_hzdx > time ~ 0L,
        assign==0 & t_hzdx == time ~ 1L,
        assign==0 & t_hzdx > time ~ 0L)
    ) %>%
    rename(treat = assign)

  d_cl_dem = d_cl %>%
    mutate(time = case_when(
      assign == 1 & t_hzv_rec_1 > l_tte_params$i_grace ~ pmin(t_demdx, t_death, t_end, l_tte_params$i_grace),
      assign == 1 & t_hzv_rec_1 <= l_tte_params$i_grace ~ pmin(t_demdx, t_death, t_end),
      assign == 0 ~ pmin(t_hzv_rec_1, t_demdx, 
                         t_death, t_end)),
      event = case_when(assign==1 & t_demdx == time ~ 1L, 
                        assign==1 ~ 0L,
                        assign==0 & t_demdx == time ~ 1L,
                        assign==0 ~ 0L)) %>%
    rename(treat = assign)

  d_ids = distinct(d_cl, PatientICN)


# Point estimate ----

  d_surv_pe_hz = broom::tidy(
    survfit(Surv(time, event) ~ treat, data=d_cl_hz)) %>%
    mutate(treat = str_extract(strata,  '(?<=treat=)\\d')
    ) %>%
    arrange(time) %>%
    select(-strata) %>%
    mutate(pr_e = 1-estimate) %>%
    rename(pr_s = estimate) %>%
    pivot_wider(id_cols = c(time), 
                names_from = treat,
                values_from = c(pr_e, pr_s,
                                n.risk, n.event)) %>%
    mutate(cir = pr_e_1 / pr_e_0,
           cid = (pr_e_1 - pr_e_0))

  d_survmin_hz = ggsurvplot(survfit(Surv(time, event) ~ treat, data=d_cl_hz),
                            data=d_cl_hz,
                            fun = 'event',
                            xlim = c(0, l_tte_params$i_fup),
                            xlab = 'Follow-up in days',
                            break.time.by=180,
                            palette = c('red', 'blue'),
                            linetype = 'strata',
                            censor=F,
                            cumevents=T,
                            conf.int=F, pval=F,
                            risk.table=T, risk.table.col = 'strata',
                            legend.labs  = c(l_tte_params$i_trt_labels[1], 
                                             l_tte_params$i_trt_labels[2]),
                            risk.table.height=0.25, 
                            ggtheme = theme_classic())
  
  ggsave( here('img', 
               paste0(prj.specs$prj.prefix, '.', 'survplot_survmin_hz', '.jpeg')),
          d_survmin_hz, width = 12, height=7, dpi=300)
  
  d_surv_pe_dem = broom::tidy(survfit(Surv(time, event) ~ treat, data=d_cl_dem)) %>%
    mutate(treat = str_extract(strata,  '(?<=treat=)\\d')
    ) %>%
    arrange(time) %>%
    select(-strata) %>%
    mutate(pr_e = 1-estimate) %>%
    rename(pr_s = estimate) %>%
    pivot_wider(id_cols = c(time), 
                names_from = treat,
                values_from = c(pr_e, pr_s,
                                n.risk, n.event)) %>%
    mutate(cir = pr_e_1 / pr_e_0,
           cid = (pr_e_1 - pr_e_0))

  d_survmin_dem = ggsurvplot(survfit(Surv(time, event) ~ treat, data=d_cl_dem),
                            data=d_cl_dem,
                            fun = 'event',
                            xlim = c(0, l_tte_params$i_fup),
                            xlab = 'Follow-up in days',
                            break.time.by=180,
                            palette = c('red', 'blue'),
                            linetype = 'strata',
                            censor=F,
                            cumevents=T,
                            conf.int=F, pval=F,
                            risk.table=T, risk.table.col = 'strata',
                            legend.labs  = c(l_tte_params$i_trt_labels[1], 
                                             l_tte_params$i_trt_labels[2]),
                            risk.table.height=0.25, 
                            ggtheme = theme_classic())
  
  ggsave( here('img', 
               paste0(prj.specs$prj.prefix, '.', 'survplot_survmin_dem', '.jpeg')),
          d_survmin_dem, width = 12, height=7, dpi=300)
  
# defined above
plan(multisession, workers = d_output$runplan$workers)
set.seed(d_output$runplan$seed)

# bootstrap 
  d_bs_hz = future_map(.x = 1:d_output$runplan$boots, 
                       .f = ~d_fun_getkmres(d_cl_hz, d_ids, .x),
                       .options = furrr_options(seed = T))
  
  d_bs_dem = future_map(.x = 1:d_output$runplan$boots, 
                        .f = ~d_fun_getkmres(d_cl_dem, d_ids, .x),
                        .options = furrr_options(seed = T))

# Combine and summarize ----
  ## Herpes Zoster ----
    d_surv_hz = d_bs_hz %>%
      bind_rows(.id = 'boot') %>%
      mutate(boot = as.numeric(boot)) 
    
    d_summ_surv_hz = d_surv_hz %>%
      group_by(time) %>%
      summarize(pr_e_0_lc = quantile(pr_e_0, 0.025, na.rm=T),
                pr_e_1_lc = quantile(pr_e_1, 0.025, na.rm=T),
                cir_lc = quantile(cir, 0.025, na.rm=T),
                cid_lc = quantile(cid, 0.025, na.rm=T),
                pr_e_0_uc = quantile(pr_e_0, 0.975, na.rm=T),
                pr_e_1_uc = quantile(pr_e_1, 0.975, na.rm=T),
                cir_uc = quantile(cir, 0.975, na.rm=T),
                cid_uc = quantile(cid, 0.975, na.rm=T)) %>%
      bind_cols(select(d_surv_pe_hz, pr_e_0, pr_e_1, cir, cid),
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
                 paste0(prj.specs$prj.prefix, '.', 'survplot_kmunadj_hz', '.jpeg')),
            plot=d_gg_1, width = 12, height=8, dpi=300)
    
d_surv_dem = d_bs_dem %>%
  bind_rows(.id = 'boot') %>%
  mutate(boot = as.numeric(boot)) 

d_summ_surv_dem = d_surv_dem %>%
  group_by(time) %>%
  summarize(pr_e_0_lc = quantile(pr_e_0, 0.025, na.rm=T),
            pr_e_1_lc = quantile(pr_e_1, 0.025, na.rm=T),
            cir_lc = quantile(cir, 0.025, na.rm=T),
            cid_lc = quantile(cid, 0.025, na.rm=T),
            pr_e_0_uc = quantile(pr_e_0, 0.975, na.rm=T),
            pr_e_1_uc = quantile(pr_e_1, 0.975, na.rm=T),
            cir_uc = quantile(cir, 0.975, na.rm=T),
            cid_uc = quantile(cid, 0.975, na.rm=T)) %>%
  bind_cols(select(d_surv_pe_dem, pr_e_0, pr_e_1, cir, cid),
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
                              'survplot_kmunadj_dem', '.jpeg')),
        plot=d_gg_1,
        width = 12, height=8, dpi=300)

# Save KM est results ----
  d_summ_surv = bind_rows(
    mutate(d_summ_surv_hz, outcome = 'hz'),
    mutate(d_summ_surv_dem, outcome = 'dem'))

  write_csv(d_summ_surv,
    file = here('out', paste0(prj.specs$prj.prefix, '.', 
                                 'kmnoadjoutc', '.csv'
    ))
  )

d_output$results$summ = d_summ_surv
d_output$results$boot = bind_rows(
  mutate(d_surv_hz, outcome = 'hz'),
  mutate(d_surv_dem, outcome = 'dem'))

saveRDS(d_output, 
        file = here('out', paste0('kmnaive', '.', f_tstmp(), '.Rds')))