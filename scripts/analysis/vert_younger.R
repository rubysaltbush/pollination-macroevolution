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

# cut tree down to only vert/insect tips
tree_vert_insect <- ape::drop.tip(tree, pollination1209$taxon_name[pollination1209$vert_insect == "?"])

# reformat vert/insect data as vector for diversitree analysis
vert_insect <- dplyr::na_if(pollination1209$vert_insect, "?")
names(vert_insect) <- pollination1209$taxon_name
# diversitree can't handle polymorphic data, convert these to NA
vert_insect <- dplyr::na_if(vert_insect, "6&5")
vert_insect <- gsub("5", "0", vert_insect)
vert_insect <- gsub("6", "1", vert_insect)

# run experimental analysis to reconstruct vert/insect pollination where
# vertebrate pollination forbidden before K-Pg mass extinction (66 mya)
# this function only possible with a binary character, thus can't run for
# whole insect-vertebrate-wind-water pollination scenario
diversitree::make.bisse.td(tree_vert_insect, states = vert_insect, n.epoch = 2,
              nt.extra=10, strict=TRUE, control=list())
# BUT I cannot work out looking through the documentation where to specify
# constant rates of speciation/extinction, or where to specify the time
# split as 66 mya, or where to specify that I want only insect pollination
# before that time split. And, is this even wise, or will it just give strange
# and spurious results?

