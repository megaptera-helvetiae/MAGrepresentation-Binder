pkgs = c("tidyverse", "rmarkdown", "httr", "shinydashboard", "leaflet")
ncores = parallel::detectCores()
install.packages(pkgs, Ncpus = ncores)
