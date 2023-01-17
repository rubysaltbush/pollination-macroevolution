# run to install all packages needed for analysis
# R version for original analysis 4.1.3

install.packages("ape") # basic package for phylogenetic trees and analysis, run version 5.6-2
install.packages("tidyverse") # for data manipulation, run version 1.3.2
install.packages("sf") # for spatial data, run version 1.0-9
install.packages("raster") # for rasterised spatial data, run version 3.6-3
install.packages("CoordinateCleaner") # to clean GBIF records, run version 2.0-20
install.packages("rnaturalearthdata") # needed to run CoordinateCleaner tests, run version 0.1.0
install.packages("countrycode") # needed to run CoordinateCleaner tests, run version 1.4.0
# if installing corHMM on Linux need to install Rmpfr package in terminal first
# using command apt install   libmpfr-dev
install.packages("corHMM") # IMPORTANT INSTALL VERSION 2.8!!! other versions work very differently, will mess up tests for correlated evolution in particular
install.packages("phytools") # phylogenetic analysis, run version 1.0-3
install.packages("phylolm") # for phylogenetic logistic regression, run version 2.6.2
install.packages("doParallel") # for running tasks in parallel on multiple cores, run version 1.0.17
install.packages("plotrix") # for plotting circular phylogenies with labels, run version 3.8-2
