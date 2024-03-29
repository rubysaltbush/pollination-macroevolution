# uncomment below to install packages needed for analysis
# source("scripts/install_dependencies.R")

library(ape)
library(tidyverse)
library(sf)
library(raster)
library(CoordinateCleaner) # package to clean GBIF records
library(countrycode) # package to convert GBIF country codes for cleaning
library(corHMM) #  IMPORTANT INSTALL VERSION 2.8!!! other versions work differently, will mess up tests for correlated evolution in particular
library(phytools)
library(phylolm)
library(parallel)
library(doParallel)
library(plotrix)

# function to cache pre-prepared R data. If RDS already in cache will read data
source("scripts/functions/cache_RDS.R")

# read in consensus tree from Ramírez-Barahona et al (2020)
tree <- ape::read.tree("data_input/eFLOWER-1209_RC_complete_MCC.phy")

# set number of cores available on local computer (max - 1) for parallelisation
no_cores <- parallel::detectCores() - 1
# set up doParallel to run functions on multiple cores
cl <- parallel::makeCluster(no_cores)
doParallel::registerDoParallel(cl)
rm(cl)

# read in prepped pollination data
source("scripts/pollination1209.R")

# create figure of pollination systems vs syndromes across tree
source("scripts/analysis/sampling.R")

# run corHMM Ancestral State Reconstructions of pollination mode in angiosperms
# takes ~4 hours to run if not already cached
source("scripts/analysis/ASR.R")

# for graphs assign fixed colours
source("scripts/my_colours.R")

# stochastic character mapping of ASR results
# exports most major figures from paper
# takes ~30 minutes to run, will have to press enter twice
source("scripts/analysis/simmap.R")

## correlating evolution of pollination to environment ##

# first try phylogenetic logistic regression for continuous and categorical variables
source("scripts/analysis/PGLS.R")

# run correlation models for superbiome/pollination evolution (Supporting Information Notes S1)
source("scripts/analysis/biomes.R")

# stochastic character mapping on posterior trees to consider phylogenetic uncertainty
# needs to be run on machine with lots of memory (~3 TB) and many (~36) cores
# will take ~4 hours if so
source("scripts/analysis/multiPhylo_simmap.R")

# stochastic character mapping on alternative (younger) phylogeny
# to explore uncertainty in timing of transitions (Supporting Information Fig. S2)
source("scripts/analysis/simmap_younger.R")

