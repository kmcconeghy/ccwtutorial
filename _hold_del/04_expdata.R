library(here())

source(here('src', 'setup.R'), echo=F)

## a priori treatments: 
# Assign = 0 (No treatment ever), 
# Assign = 1 (treat by time 12); i.e. grace period of 12

# Outcome panel dataset ----

  # start - time 0 (1)
  # end - censor time (event, censor, end of follow-up etc.)
  d_cloned = readRDS(here('dta', 'survdta_cloned.R')) %>%
    mutate(start=1,
           end = t_clone)

  # Expand data so one row per unit of follow-up
  # data.table is faster  
    setDT(d_cloned)
    
    d_panel = d_cloned[rep(seq(.N), t_clone)]
                       
    d_panel[, exit := (seq_len(.N)), by = list(id, assign)]
    d_panel[, enter := exit-1]
    d_panel[, time := seq_len(.N), by = list(id, assign)]
    
    # Outcome is = 1 in row where event occurred
    d_panel[, event_outc := if_else( t_clone <= time, event_outc, 0L), by = list(id, assign)]
    
    setDF(d_panel)
    
    d_panel_2 = select(d_panel, id, time, event_outc, t_treat, assign, enter, exit, end) 

  saveRDS(d_panel_2, here('dta', 'survdta_cloned_panel.R'))

# TIME-VARYING FOR WEIGHTS ----
  d_treat = readRDS(here('dta', 'survdta_cloned.R')) %>%
    # keep one clone
    dplyr::filter(assign==0) %>%
    mutate(start=1,
           end = t_treat) %>%
    group_by(id) %>%
    mutate(
      time = t_clone, # censoring time is same as treatment time unless dead first
      outcome = if_else(time==t_treat, 0, 1)
    ) %>%
    ungroup

  # Expand data so one row per unit of follow-up
  # data.table is faster  
  setDT(d_treat)
  
  d_panel = d_treat[rep(seq(.N), time)]
  
  d_panel[, exit := (seq_len(.N)), by = list(id)]
  d_panel[, enter := exit-1]
  d_panel[, time := seq_len(.N), by = list(id)]
  
  # Outcome is = 1 in row where treat occurred
  d_panel[, event_treat := if_else(t_treat <= time, treat, 0L), by = list(id)]
  
  setDF(d_panel)
  
  d_panel_2 = select(d_panel, id, time, event_treat, t_treat, enter, exit, end, X1, X2) 

  saveRDS(d_panel_2, here('dta', 'survdta_treat_panel.R'))
  