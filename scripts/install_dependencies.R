# run to install all packages needed for analysis

install.packages("ape") # basic package for phylogenetic trees and analysis
install.packages("tidyverse") # for data manipulation
install.packages("sf") # for spatial data
install.packages("raster") # for rasterised spatial data
install.packages("CoordinateCleaner") # to clean GBIF records
install.packages("rnaturalearthdata") # needed to run CoordinateCleaner tests
install.packages("countrycode") # needed to run CoordinateCleaner tests
install.packages("geodata") # downloads worldclim temp and precipitation data
# if installing corHMM on Linux need to install Rmpfr package in terminal first
# using command apt install   libmpfr-dev
install.packages("corHMM") # IMPORTANT INSTALL VERSION 2.8!!!
install.packages("phytools") # phylogenetic analysis
install.packages("phylolm") # for phylogenetic logistic regression
install.packages("doParallel") # for running tasks in parallel on multiple cores
install.packages("plotrix") # for plotting circular phylogenies with labels
