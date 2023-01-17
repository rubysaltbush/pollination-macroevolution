# pollination-macroevolution
Data and analysis for the paper *Insects pollinated flowering plants for most of angiosperm evolutionary history*

Authors: Ruby E. Stephens (1,2), Rachael V. Gallagher (1,3), Lily Dun (2,4), Will Cornwell (4), Hervé Sauquet (2,4)
+ Corresponding author: Ruby E. Stephens, stephenseruby@gmail.com

Author addresses:

1. School of Natural Sciences, Macquarie University, Sydney, Australia
2. National Herbarium of New South Wales (NSW), Royal Botanic Gardens and Domain Trust, Sydney, Australia
3. Hawkesbury Institute for the Environment, Western Sydney University, Sydney, Australia
4. Evolution and Ecology Research Centre, University of New South Wales, Sydney, Australia

## Re-running analysis

Using RStudio open [main.R](https://github.com/rubysaltbush/pollination-macroevolution/blob/main/main.R) 
and run scripts in order given in this main script.

Necessary packages can be installed by running [install_dependencies.R](https://github.com/rubysaltbush/pollination-macroevolution/blob/main/scripts/install_dependencies.R).
Check package versions mentioned in this script are consistent with your installed versions. Originally run in R version 4.1.3

Two data files are too large and not available from the GitHub version but can be downloaded from the Zenodo archive at **insert doi here**

These are:

1.`SA_thin_trees_corrected.tre`**insert link here** (155 mb, multiple angiosperm phylogenies from the posterior of Ramírez-Barahona et al. (2020)); and 

2.`filtered_obs_ruby.csv`**insert link here** (**X**mb, GBIF records for 1201 taxa downloaded and initially filtered by Will Cornwell)

The code caches several steps of the data cleaning and results so most analyses can run without needing to download the above files.
Running a completely fresh analysis will require downloading the above files, deleting all cached results and output in
[results](https://github.com/rubysaltbush/pollination-macroevolution/tree/main/results) and
[data_output](https://github.com/rubysaltbush/pollination-macroevolution/tree/main/data_output), and then running
[main.R](https://github.com/rubysaltbush/pollination-macroevolution/blob/main/main.R) from the start. 
To complete the posterior stochastic mapping in 
[multiPhylo_simmap.R](https://github.com/rubysaltbush/pollination-macroevolution/blob/main/scripts/analysis/multiPhylo_simmap.R) 
you will need access to a computer with a lot of CPUs and memory (I ran on a machine with 36 CPUs and 300GB memory).