# Define interventions ----

## Three a priori treatments: 
# Assign = 0 (No shingrix ever), 
# Assign = 1 (1+ doses by grace)
# Assign = 2 (2 doses by grace)

#### t_hzv_rec_1 = time to first vaccine, if none = INF
#### t_hzv_rec_2 = time to 2nd vaccine, if none = INF

## Modify panel ----
d_cl_pnl = bind_rows(# d_trt_0 - No shingrix ever ----
                     d_pnl %>%
                       mutate(assign = 0,
                              # censor if vaccinated in interval
                              censor_art = if_else(t_hzv_rec_1 < cumdays, 1L, 0L)
                       ) %>%
                       # keep uncensored intervals, and first interval with censor event
                       dplyr::filter(censor_art == 0 | (t_hzv_rec_1 >= cumdays-i_control$i_int)), 
                     # d_trt_1 - 1+ shingrix ----
                     d_pnl %>%
                       mutate(assign = 1,
                              censor_art = case_when(
                                # cannot censor prior to grace
                                cumdays < i_control$i_grace ~ 0L,
                                # censored if not vaccinated at last interval of grace
                                cumdays >= i_control$i_grace & t_hzv_rec_1 > i_control$i_grace ~ 1L,
                                cumdays >= i_control$i_grace & t_hzv_rec_1 <=i_control$i_grace ~ 0L 
                              )
                       ) %>%
                       # keep if vaccinated in grace
                       dplyr::filter(censor_art == 0 | cumdays == i_control$i_grace), 
                     # d_trt_2 - 2+ Shingrix ----
                     d_pnl %>%
                       mutate(assign = 2,
                              censor_art = case_when(
                                # cannot censor prior to grace
                                cumdays < i_control$i_grace ~ 0L,
                                # censored if not vaccinated at last interval of grace
                                cumdays>= i_control$i_grace & t_hzv_rec_2 > i_control$i_grace ~ 1L,
                                cumdays>= i_control$i_grace & t_hzv_rec_2 <=i_control$i_grace ~ 0L 
                              )
                       ) %>%
                       # keep if vaccinated in grace
                       dplyr::filter(censor_art == 0 | cumdays == i_control$i_grace)) %>%
  dplyr::filter(cumdays <= i_control$i_fup) %>%
  mutate(Strategy := factor(assign,
                            levels = c(0, 1, 2),
                            labels = l_tte_params$i_trt_labels)
         )