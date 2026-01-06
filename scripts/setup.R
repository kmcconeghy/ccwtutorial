  #source(here('r', 'onetime.R'))

  pkg_toload <- c('tidyverse', 'dplyr', 'tidyr', 'stringr',
                  'lubridate', 'here', 'knitr', 'kableExtra', 'quarto',
                  'survival', 'future', 'progressr', 
                  'ggpubr', 'survminer', 'purrr', 'furrr',
                  'data.table', 'parglm', 'gtsummary')
  
  hold_del <- sapply(pkg_toload, require, 
                     warn.conflicts=F, quietly=T,
                     character.only=T)
  
  cbbPalette = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  for (i in list.files(here('scripts', 'macros'))) {
    source(here('scripts', 'macros', i))
  }