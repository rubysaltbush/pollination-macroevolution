# in this script I hope to run a BiSSE analysis on the tree with the vertebrate
# pollination state restricted to the Cenozoic (66 mya) to complement the main
# analyses of transition timing in text

# THOUGH I wonder if I do, will vertebrate pollination just be reconstructed as
# evolving basically at the edge of this boundary? Shall see.

#### set up environment ####
# remove all objects in global environment
rm(list=ls())
# run "garbage collection" to clear memory
gc()
# restart R ONLY WORKS IN RSTUDIO
.rs.restartR()

# set number of cores available on local computer (max - 1) for parallelisation
no_cores <- parallel::detectCores() - 1
# set up doParallel to run functions on multiple cores
cl <- parallel::makeCluster(no_cores)
doParallel::registerDoParallel(cl)
rm(cl)

# install (uncomment if necessary) and load diversitree package necessary for BiSSE analysis
# install.packages("diversitree")
library(diversitree)

# function to cache pre-prepared R data. If RDS already in cache will read data
source("scripts/functions/cache_RDS.R")

# for graphs assign fixed colours
source("scripts/my_colours.R")

#### prepare data ####

# read in consensus tree from Ramírez-Barahona et al (2020)
tree <- ape::read.tree("data_input/eFLOWER-1209_RC_complete_MCC.phy")

# read in prepped pollination data
source("scripts/pollination1209.R")




