  #source(here('r', 'onetime.R'))

  pkg_toload <- c('tidyverse', 'dplyr', 'tidyr', 'stringr',
                  'lubridate', 'here', 'knitr', 'kableExtra', 'quarto',
                  'survival', 'future', 'progressr', 
                  'ggpubr', 'survminer', 'furrr',
                  'data.table', 'parglm', 'gtsummary')
  
  hold_del <- sapply(pkg_toload, require, 
                     warn.conflicts=F, quietly=T,
                     character.only=T)
