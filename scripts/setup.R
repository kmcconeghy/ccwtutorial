  #source(here('r', 'onetime.R'))

  pkg_toload <- c('tidyverse', 
                  'lubridate', 'here', 'knitr', 'quarto',
                  'survival', 'future', 'progressr', 
                  'ggpubr', 'survminer', 'furrr',
                  'data.table', 'parglm')
  
  hold_del <- sapply(pkg_toload, require, 
                     warn.conflicts=F, quietly=T,
                     character.only=T)
  here()
