pkgs = c("tidyverse", "rmarkdown", "httr", "shinydashboard", "leaflet", "gplots", "plotly", "DT")
ncores = parallel::detectCores()
install.packages(pkgs, Ncpus = ncores)
