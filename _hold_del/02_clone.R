# Define interventions ----

## Three a priori treatments: 
# Assign = 0 (No treatment ever), 
# Assign = 1 (treat by time 12); i.e. grace period of 12

## Modify panel ----
data_cloned = bind_rows(
                    # Assign = 0
                     data %>%
                       mutate(assign = 0,
                              # artificial censor if treated in interval
                              t_artcens = if_else(t_treat < time, t_treat, Inf)
                       ),
                     # Assign = 1
                     data %>%
                       mutate(assign = 1,
                              t_artcens = if_else(t_treat > 12, 12, Inf)
                              )
                       ) %>%
  # Generate new failure time, accounting for artificial censoring
  mutate(t_clone = pmin(time, t_artcens),
         # update event counts
         event_outc = if_else(t_clone<time, 0, event_outc))

head(data_cloned)

saveRDS(select(data_cloned, id, assign, t_clone, event_outc, time, t_artcens, t_treat, treat, X1, X2),
               here('dta', 'survdta_cloned.R'))