#### simmap with multiphylo ####
# using RB2020 posterior trees to quantify impact of phylogenetic uncertainty,
# run stochastic mapping, with parallelisation to speed it up a bit

# TO SAVE MEMORY this code is run on each different simmapping separately
# with results exported in between, and R restarted in between each run
# to restart R, this code assumes you are using RStudio

# remove all objects in global environment
rm(list=ls())
# run "garbage collection" to clear memory
gc()
# restart R ONLY WORKS IN RSTUDIO
.rs.restartR()

#### WIND vs. ANIMAL ####

# set number of cores available on local computer (MAX for remote machine) for parallelisation
no_cores <- parallel::detectCores()
# set up doParallel to run functions on multiple cores
cl <- parallel::makeCluster(no_cores)
doParallel::registerDoParallel(cl)
rm(cl)

# read in pollination data
pollination1209 <- readr::read_csv("data_output/pollination1209.csv")

# read in posterior trees from Ramírez-Barahona et al (2020)
tree_multi <- ape::read.tree("data_input/SA_thin_trees_corrected.tre")
# remove first 25% of posterior trees in case these are "burn-in"
tree_multi <- tree_multi[738:2948]

# drop tips of multitree where binary characters missing
tree_multi_wind_animal <- phytools::drop.tip.multiPhylo(tree_multi, pollination1209$taxon_name[pollination1209$wind_animal == "?"])
# this also drops tips where pollination mode unknown

# prepare matrices with extra character states dropped
matrices <- list()

matrices$wind_animal <- pollination1209 %>%
  dplyr::select(taxon_name, value = wind_animal) %>%
  dplyr::filter(value != "?") %>%
  as.matrix()

# subsample 1000 trees from multiPhylo
multi_trees <- list()
multi_trees$wind_animal <- sample(tree_multi_wind_animal, size = 1000)
rm(tree_multi_wind_animal, tree_multi)
gc() # run "garbage collection" here to free up memory after removing trees

# run corHMM to estimate Q matrix and then simmapping
# foreach-looping through each posterior tree
multi_sims <- list()
start_time <- Sys.time()
for (name in names(matrices)) { # for each of my ONE binary characters
  data <- matrices[[name]] # use binary character data
  # foreach parallelisation here
  multi_sims[[name]] <- foreach(i = 1:length(multi_trees[[name]])) %dopar% { # then for each tree in multiPhylo
    phy <- multi_trees[[name]][[i]] # use that tree
    # run corHMM to get Q matrix
    ASR_for_Q <- corHMM::corHMM(
      phy = phy,
      data = data,
      model = "ARD",
      rate.cat = 1,
      nstarts = 1,
      n.cores = 1
    )
    model <- ASR_for_Q$solution
    model[is.na(model)] <- 0
    diag(model) <- -rowSums(model)
    # then run corHMM simmapping
    result <- corHMM::makeSimmap(
      tree = phy,
      data = data,
      model = model,
      rate.cat = 1,
      nSim = 100,
      nCores = 1
    )
    class(result) <- c("multiSimmap", "multiPhylo")
    result
  }
  
  # combine into one simmap
  multi_sims[[name]] <- do.call(c, multi_sims[[name]])
  class(multi_sims[[name]]) <- c("multiSimmap", "multiPhylo")
}
end_time <- Sys.time()
end_time - start_time
# Time difference of 24.28927 mins for 100 sims * 1000 trees on 36 cores
rm(data, name, matrices, start_time, end_time, multi_trees)
gc() # run "garbage collection" here to free up memory after removing things

# describe and calculate 95% confidence intervals etc for # of transitions in simmaps
# for SOME REASON rewriting below as foreach makes it take longer and use more memory
start_time <- Sys.time()
multisim_descriptions <- list()
multisim_density <- list()
for(name in names(multi_sims)){
  multisim_descriptions[[name]] <- phytools::describe.simmap(multi_sims[[name]])
  multisim_density[[name]] <- phytools::density.multiSimmap(multi_sims[[name]])
}
end_time <- Sys.time()
end_time - start_time
# Time difference of 1.137062 hours for 100 sims * 1000 trees on 36 cores
rm(name, start_time, end_time, multi_sims)
gc() # run "garbage collection" here to free up memory after removing things

# export print results from multi-simmapping to text file
sink(file = "results/multi_simmap_descriptions_windanimal.txt")
for(name in names(multisim_descriptions)){
  print(paste(name))
  print(multisim_descriptions[[name]])
  print(paste(name))
  print(multisim_density[[name]])
  print(paste(" "))
}
sink(file = NULL)

# remove all objects and clean memory before starting next round
# remove all objects in global environment
rm(list=ls())
# run "garbage collection" to clear memory
gc()
# restart R ONLY WORKS IN RSTUDIO
.rs.restartR()


#### INSECT vs VERTEBRATE ####


# set number of cores available on local computer (max) for parallelisation
no_cores <- parallel::detectCores()
# set up doParallel to run functions on multiple cores
cl <- parallel::makeCluster(no_cores)
doParallel::registerDoParallel(cl)
rm(cl)

# read in pollination data
pollination1209 <- readr::read_csv("data_output/pollination1209.csv")

# read in posterior trees from Ramírez-Barahona et al (2020)
tree_multi <- ape::read.tree("data_input/SA_thin_trees_corrected.tre")
# remove first 25% of posterior trees in case these are "burn-in"
tree_multi <- tree_multi[738:2948]

# drop tips of multitree where binary characters missing
tree_multi_vert_insect <- phytools::drop.tip.multiPhylo(tree_multi, pollination1209$taxon_name[pollination1209$vert_insect == "?"])
# this also drops tips where pollination mode unknown

# prepare matrices with extra character states dropped
matrices <- list()

matrices$vert_insect <- pollination1209 %>%
  dplyr::select(taxon_name, value = vert_insect) %>%
  dplyr::filter(value != "?") %>%
  as.matrix()

# subsample 1000 trees from multiPhylo
multi_trees <- list()
multi_trees$vert_insect <- sample(tree_multi_vert_insect, size = 1000)
rm(tree_multi_vert_insect, tree_multi)
gc() # run "garbage collection" here to free up memory after removing trees

# run corHMM to estimate Q matrix and then simmapping
# foreach-looping through each posterior tree
multi_sims <- list()
start_time <- Sys.time()
for (name in names(matrices)) { # for each of my ONE binary characters
  data <- matrices[[name]] # use binary character data
  # foreach parallelisation here
  multi_sims[[name]] <- foreach(i = 1:length(multi_trees[[name]])) %dopar% { # then for each tree in multiPhylo
    phy <- multi_trees[[name]][[i]] # use that tree
    # run corHMM to get Q matrix
    ASR_for_Q <- corHMM::corHMM(
      phy = phy,
      data = data,
      model = "ARD",
      rate.cat = 1,
      nstarts = 1,
      n.cores = 1
    )
    model <- ASR_for_Q$solution
    model[is.na(model)] <- 0
    diag(model) <- -rowSums(model)
    # then run corHMM simmapping
    result <- corHMM::makeSimmap(
      tree = phy,
      data = data,
      model = model,
      rate.cat = 1,
      nSim = 100,
      nCores = 1
    )
    class(result) <- c("multiSimmap", "multiPhylo")
    result
  }
  
  # combine into one simmap
  multi_sims[[name]] <- do.call(c, multi_sims[[name]])
  class(multi_sims[[name]]) <- c("multiSimmap", "multiPhylo")
}
end_time <- Sys.time()
end_time - start_time
# Time difference of 18.86719 mins for 100 sims * 1000 trees on 36 cores
rm(data, name, matrices, start_time, end_time, multi_trees)
gc() # run "garbage collection" here to free up memory after removing things

# describe and calculate 95% confidence intervals etc for # of transitions in simmaps
# for SOME REASON rewriting below as foreach makes it take longer and use more memory
start_time <- Sys.time()
multisim_descriptions <- list()
multisim_density <- list()
for(name in names(multi_sims)){
  multisim_descriptions[[name]] <- phytools::describe.simmap(multi_sims[[name]])
  multisim_density[[name]] <- phytools::density.multiSimmap(multi_sims[[name]])
}
end_time <- Sys.time()
end_time - start_time
# Time difference of 55.46987 mins for 100 sims * 1000 trees on 36 cores
rm(name, start_time, end_time, multi_sims)
gc() # run "garbage collection" here to free up memory after removing things

# export print results from multi-simmapping to text file
sink(file = "results/multi_simmap_descriptions_vertinsect.txt")
for(name in names(multisim_descriptions)){
  print(paste(name))
  print(multisim_descriptions[[name]])
  print(paste(name))
  print(multisim_density[[name]])
  print(paste(" "))
}
sink(file = NULL)

# remove all objects and clean memory before starting next round
# remove all objects in global environment
rm(list=ls())
# run "garbage collection" to clear memory
gc()
# restart R ONLY WORKS IN RSTUDIO
.rs.restartR()


#### ABIOTIC vs ANIMAL ####


# set number of cores available on local computer (max) for parallelisation
no_cores <- parallel::detectCores()
# set up doParallel to run functions on multiple cores
cl <- parallel::makeCluster(no_cores)
doParallel::registerDoParallel(cl)
rm(cl)

# read in pollination data
pollination1209 <- readr::read_csv("data_output/pollination1209.csv")

# read in posterior trees from Ramírez-Barahona et al (2020)
tree_multi <- ape::read.tree("data_input/SA_thin_trees_corrected.tre")
# remove first 25% of posterior trees in case these are "burn-in"
tree_multi <- tree_multi[738:2948]

# drop tips of multitree where binary characters missing
tree_multi_abiotic_animal <- phytools::drop.tip.multiPhylo(tree_multi, pollination1209$taxon_name[pollination1209$abiotic_animal == "?"])
# this also drops tips where pollination mode unknown

# prepare matrices with extra character states dropped
matrices <- list()

matrices$abiotic_animal <- pollination1209 %>%
  dplyr::select(taxon_name, value = abiotic_animal) %>%
  dplyr::filter(value != "?") %>%
  as.matrix()

# subsample 1000 trees from multiPhylo
multi_trees <- list()
multi_trees$abiotic_animal <- sample(tree_multi_abiotic_animal, size = 1000)
rm(tree_multi_abiotic_animal, tree_multi)
gc() # run "garbage collection" here to free up memory after removing trees

# run corHMM to estimate Q matrix and then simmapping
# foreach-looping through each posterior tree
multi_sims <- list()
start_time <- Sys.time()
for (name in names(matrices)) { # for each of my ONE binary characters
  data <- matrices[[name]] # use binary character data
  # foreach parallelisation here
  multi_sims[[name]] <- foreach(i = 1:length(multi_trees[[name]])) %dopar% { # then for each tree in multiPhylo
    phy <- multi_trees[[name]][[i]] # use that tree
    # run corHMM to get Q matrix
    ASR_for_Q <- corHMM::corHMM(
      phy = phy,
      data = data,
      model = "ARD",
      rate.cat = 1,
      nstarts = 1,
      n.cores = 1
    )
    model <- ASR_for_Q$solution
    model[is.na(model)] <- 0
    diag(model) <- -rowSums(model)
    # then run corHMM simmapping
    result <- corHMM::makeSimmap(
      tree = phy,
      data = data,
      model = model,
      rate.cat = 1,
      nSim = 100,
      nCores = 1
    )
    class(result) <- c("multiSimmap", "multiPhylo")
    result
  }
  
  # combine into one simmap
  multi_sims[[name]] <- do.call(c, multi_sims[[name]])
  class(multi_sims[[name]]) <- c("multiSimmap", "multiPhylo")
}
end_time <- Sys.time()
end_time - start_time
# Time difference of 26.27052 mins for nSim = 100 on 1000 posterior trees
# on 3TB memory and 36 cores: takes about 100GB of memory (yikes)
rm(data, name, matrices, start_time, end_time, multi_trees)
gc() # run "garbage collection" here to free up memory after removing things

# describe and calculate 95% confidence intervals etc for # of transitions in simmaps
# for SOME REASON rewriting below as foreach makes it take longer and use more memory
start_time <- Sys.time()
multisim_descriptions <- list()
multisim_density <- list()
for(name in names(multi_sims)){
  multisim_descriptions[[name]] <- phytools::describe.simmap(multi_sims[[name]])
  multisim_density[[name]] <- phytools::density.multiSimmap(multi_sims[[name]])
}
end_time <- Sys.time()
end_time - start_time
# Time difference of 1.223176 hours for 100 sims * 1000 trees on 36 cores!!
rm(name, start_time, end_time, multi_sims)
gc() # run "garbage collection" here to free up memory after removing things

# export print results from multi-simmapping to text file
sink(file = "results/multi_simmap_descriptions_abioticanimal.txt")
for(name in names(multisim_descriptions)){
  print(paste(name))
  print(multisim_descriptions[[name]])
  print(paste(name))
  print(multisim_density[[name]])
  print(paste(" "))
}
sink(file = NULL)

# remove all objects and clean memory before starting next round
# remove all objects in global environment
rm(list=ls())
# run "garbage collection" to clear memory
gc()
# restart R ONLY WORKS IN RSTUDIO
.rs.restartR()


#### WIND WATER VERT INSECT ####


# set number of cores available on local computer (max) for parallelisation
no_cores <- parallel::detectCores()
# set up doParallel to run functions on multiple cores
cl <- parallel::makeCluster(no_cores)
doParallel::registerDoParallel(cl)
rm(cl)

# read in pollination data
pollination1209 <- readr::read_csv("data_output/pollination1209.csv")

# read in posterior trees from Ramírez-Barahona et al (2020)
tree_multi <- ape::read.tree("data_input/SA_thin_trees_corrected.tre")
# remove first 25% of posterior trees in case these are "burn-in"
tree_multi <- tree_multi[738:2948]

# drop tips of multitree where binary characters missing
tree_multi_wind_water_vert_insect <- phytools::drop.tip.multiPhylo(tree_multi, pollination1209$taxon_name[pollination1209$wind_water_vert_insect == "?"])
# this also drops tips where pollination mode unknown

# prepare matrices with extra character states dropped
matrices <- list()

matrices$wind_water_vert_insect <- pollination1209 %>%
  dplyr::select(taxon_name, value = wind_water_vert_insect) %>%
  dplyr::filter(value != "?") %>%
  as.matrix()

# subsample 1000 trees from multiPhylo
multi_trees <- list()
multi_trees$wind_water_vert_insect <- sample(tree_multi_wind_water_vert_insect, size = 1000)
rm(tree_multi_wind_water_vert_insect, tree_multi)
gc() # run "garbage collection" here to free up memory after removing trees

# run corHMM to estimate Q matrix and then simmapping
# foreach-looping through each posterior tree
multi_sims <- list()
start_time <- Sys.time()
for (name in names(matrices)) { # for each of my ONE binary characters
  data <- matrices[[name]] # use binary character data
  # foreach parallelisation here
  multi_sims[[name]] <- foreach(i = 1:length(multi_trees[[name]])) %dopar% { # then for each tree in multiPhylo
    phy <- multi_trees[[name]][[i]] # use that tree
    # run corHMM to get Q matrix
    ASR_for_Q <- corHMM::corHMM(
      phy = phy,
      data = data,
      model = "ARD",
      rate.cat = 1,
      nstarts = 1,
      n.cores = 1
    )
    model <- ASR_for_Q$solution
    model[is.na(model)] <- 0
    diag(model) <- -rowSums(model)
    # then run corHMM simmapping
    result <- corHMM::makeSimmap(
      tree = phy,
      data = data,
      model = model,
      rate.cat = 1,
      nSim = 100,
      nCores = 1
    )
    class(result) <- c("multiSimmap", "multiPhylo")
    result
  }
  
  # combine into one simmap
  multi_sims[[name]] <- do.call(c, multi_sims[[name]])
  class(multi_sims[[name]]) <- c("multiSimmap", "multiPhylo")
}
end_time <- Sys.time()
end_time - start_time
# Time difference of 2.805684 hours for 100 sims * 1000 trees on 36 cores
rm(data, name, matrices, start_time, end_time, multi_trees)
gc() # run "garbage collection" here to free up memory after removing things

# describe and calculate 95% confidence intervals etc for # of transitions in simmaps
start_time <- Sys.time()
multisim_descriptions <- list()
multisim_density <- list()
for(name in names(multi_sims)){
  multisim_descriptions[[name]] <- phytools::describe.simmap(multi_sims[[name]])
  multisim_density[[name]] <- phytools::density.multiSimmap(multi_sims[[name]])
}
end_time <- Sys.time()
end_time - start_time
# Time difference of 8.134669 hours for 100 sims * 1000 trees on 36 cores
rm(name, start_time, end_time, multi_sims)
gc() # run "garbage collection" here to free up memory after removing things

# export print results from multi-simmapping to text file
sink(file = "results/multi_simmap_descriptions_windwatervertinsect.txt")
for(name in names(multisim_descriptions)){
  print(paste(name))
  print(multisim_descriptions[[name]])
  print(paste(name))
  print(multisim_density[[name]])
  print(paste(" "))
}
sink(file = NULL)

# remove all objects and clean memory before starting next round
# remove all objects in global environment
rm(list=ls())
# run "garbage collection" to clear memory
gc()
# restart R ONLY WORKS IN RSTUDIO
.rs.restartR()

# for now just exporting these results, don't think consensus tree or density
# maps useful for this tiny bit of extra analysis



