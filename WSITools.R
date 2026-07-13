# WSITools ####
# https://github.com/tkcaccia/wsiTools/blob/main/docs/installation.md

remotes::install_github(
  "tkcaccia/wsiTools",
  upgrade = "never",
  build_vignettes = FALSE
)

library(wsiTools)
wsi_backends()
wsi_setup_report()
wsi_diagnose(live_test = FALSE)