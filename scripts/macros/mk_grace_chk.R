f_clonedta = function(x, grace=12) {
  data_cloned = bind_rows( 
    x %>% 
      mutate(assign = 0, 
             censor = if_else(t_treat <= time, 1L, 0L), 
             event = if_else(censor==0 & t_outcome<=time, 1L, 0L), 
      ),
    x %>% 
      mutate(assign = 1, 
             censor = if_else(t_treat > grace & time>=grace, 1L, 0L), 
             event  = if_else(censor==0 & t_outcome<=time, 1L, 0L)
      )
  ) %>% 
    arrange(id, assign, time) %>%
    group_by(id, assign) %>%
      mutate(t_censor = min(time[censor==1])) %>%
      dplyr::filter(censor == 0 | time == min(time[censor==1])) %>% # <3>
    ungroup %>%
    select(id, time, t_outcome, t_treat, t_censor, censor, event, assign, 
           enter, exit, everything(), -any_of(c('fup', 'log_odds', 'p')))
}