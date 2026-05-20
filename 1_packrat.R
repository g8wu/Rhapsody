# Packrat Docu: https://rstudio.github.io/packrat/commands.html ####
# Packrat initialization
install.packages("packrat")
library(packrat)

# Ensure package download repos include CRAN & BioC
options(repos = BiocManager::repositories())

# Initialize in desired directory (preferably project folder)
packrat::init()

# Save current state of session
packrat::snapshot()

# Restore most recent snapshot
packrat::restore()

# Remove unused packages
packrat::clean()

# Bundle for sharing
packrat::bundle()
packrat::unbundle()

# Unlink packrat entirely, can delete packrat folder too
unlink("packrat", recursive = T)