#####correlated evolution of biomes and pollination mode#####

# biomes shows biome by cleaned GBIF records of contemporary species distribution, 
# with different thresholds for inclusion in a particular biome. 
# e.g. SB50 = >50% of species records are in a particular superbiome 
# (temperate = 0, tropical = 1, arid = 2).
# 7 species with arctic/antarctic distributions are excluded.

# pollination binary trait of interest = wind/animal pollination

# convert data to morphological df for analysis of correlated evolution
pollination_biome <- pollination1209 %>%
  dplyr::select(taxon_name, wind_water_animal, SB50) %>%
  as.data.frame()

#### prepare tree ####
# drop water pollinated tips from tree (leaving missing data as is)
to_drop <- pollination_biome %>%
  dplyr::filter(wind_water_animal %in% c("3", "2&3", "3&4"))
tree_nowater <- ape::drop.tip(tree, to_drop$taxon_name)
# and remove water pollinated taxa from morphological data
pollination_biome <- pollination_biome %>%
  dplyr::filter(wind_water_animal %in% c("?", "2", "4", "2&4")) %>%
  dplyr::rename(wind_animal = wind_water_animal)
rm(to_drop)

#### prepare morphological data ####

pollbiome <- list()

# for analyses with biomes as binary trait, can separate by tree type
# as binary trait 3 ways to split biomes - trop/extra-trop, arid/non-arid, temp/non-temp

## tropical vs non-tropical ##
pollbiome$trop <- pollination_biome
# convert superbiome to tropical/extra-tropical by replacing "arid" (2) with 0 (non-tropical)
pollbiome$trop$tropSB50 <- gsub("2", "0", pollbiome$trop$SB50)
# reduces df to 2 binary traits:
# 0 = non-tropical
# 1 = tropical
# 2 = wind pollinated
# 4 = animal pollinated
pollbiome$trop <- pollbiome$trop %>%
  dplyr::select(taxon_name, wind_animal, tropSB50)
table(pollbiome$trop$wind_animal)
table(pollbiome$trop$tropSB50)
pollbiome$trop %>% dplyr::group_by(tropSB50, wind_animal) %>% dplyr::summarise(n())

## arid vs non-arid ##
pollbiome$arid <- pollination_biome
# convert superbiome to arid/extra-arid by replacing "tropical" (1) with 0 (non-arid)
pollbiome$arid$aridSB50 <- gsub("1", "0", pollbiome$arid$SB50)
pollbiome$arid$aridSB50 <- gsub("2", "1", pollbiome$arid$aridSB50)
# reduces df to 2 binary traits:
# 0 = non-arid
# 1 = arid
# 2 = wind pollinated
# 4 = animal pollinated
pollbiome$arid <- pollbiome$arid %>%
  dplyr::select(taxon_name, wind_animal, aridSB50)
table(pollbiome$arid$wind_animal)
table(pollbiome$arid$aridSB50)
pollbiome$arid %>% dplyr::group_by(aridSB50, wind_animal) %>% dplyr::summarise(n())
# only 24 wind pollinated arid species, not very many

## temperate vs non-temperate ##
pollbiome$temp <- pollination_biome
# convert superbiome to temperate/non-temperate
pollbiome$temp$tempSB50 <- gsub("0", "3", pollbiome$temp$SB50)
pollbiome$temp$tempSB50 <- gsub("2", "0", pollbiome$temp$tempSB50)
pollbiome$temp$tempSB50 <- gsub("1", "0", pollbiome$temp$tempSB50)
pollbiome$temp$tempSB50 <- gsub("3", "1", pollbiome$temp$tempSB50)
# reduces df to 2 binary traits:
# 0 = non-temperate
# 1 = temperate
# 2 = wind pollinated
# 4 = animal pollinated
pollbiome$temp <- pollbiome$temp %>%
  dplyr::select(taxon_name, wind_animal, tempSB50)
table(pollbiome$temp$wind_animal)
table(pollbiome$temp$tempSB50)
pollbiome$temp %>% dplyr::group_by(tempSB50, wind_animal) %>% dplyr::summarise(n())

#### prepare rate matrices ####

# will need 4 different rate matrices - store all in list
rate_matrix <- list()

# all biomes
# correlated model
rate_matrix$allbiomes_corr <- corHMM::getStateMat4Dat(pollination_biome, model = "ARD", dual = FALSE)$rate.mat
rate_matrix$allbiomes_corr
# uncorrelated model
rate_matrix$allbiomes_uncorr <- corHMM::getStateMat4Dat(pollination_biome, model = "ARD", dual = FALSE)$rate.mat
# equate certain transitions for uncorrelated model MAKE SURE THIS GOES AS EXPECTED
rate_matrix$allbiomes_uncorr <- corHMM::equateStateMatPars(rate_matrix$allbiomes_uncorr, c(1,11))
rate_matrix$allbiomes_uncorr
rate_matrix$allbiomes_uncorr <- corHMM::equateStateMatPars(rate_matrix$allbiomes_uncorr, c(1,10))
rate_matrix$allbiomes_uncorr
rate_matrix$allbiomes_uncorr <- corHMM::equateStateMatPars(rate_matrix$allbiomes_uncorr, c(1,4,7))
rate_matrix$allbiomes_uncorr
rate_matrix$allbiomes_uncorr <- corHMM::equateStateMatPars(rate_matrix$allbiomes_uncorr, c(1,7))
rate_matrix$allbiomes_uncorr
rate_matrix$allbiomes_uncorr <- corHMM::equateStateMatPars(rate_matrix$allbiomes_uncorr, c(1,6))
rate_matrix$allbiomes_uncorr
rate_matrix$allbiomes_uncorr <- corHMM::equateStateMatPars(rate_matrix$allbiomes_uncorr, c(1,6))
rate_matrix$allbiomes_uncorr
rate_matrix$allbiomes_uncorr <- corHMM::equateStateMatPars(rate_matrix$allbiomes_uncorr, c(1,5))
rate_matrix$allbiomes_uncorr
rate_matrix$allbiomes_uncorr <- corHMM::equateStateMatPars(rate_matrix$allbiomes_uncorr, c(1,2,3))
rate_matrix$allbiomes_uncorr

# binary biomes (trop/non-trop etc)
# CORRELATED MODEL, i.e. ARD (works with dual trait data automatically yay!)
rate_matrix$binary_corr <- corHMM::getStateMat4Dat(pollbiome$trop, model = "ARD", dual = FALSE)$rate.mat
rate_matrix$binary_corr

#UNCORRELATED MODEL, i.e. ARD
rate_matrix$binary_uncorr <- corHMM::getStateMat4Dat(pollbiome$trop, model = "ARD", dual = FALSE)$rate.mat
rate_matrix$binary_uncorr

#equates certain transitions for uncorrelated model MAKE SURE AS EXPECTED
rate_matrix$binary_uncorr <- corHMM::equateStateMatPars(rate_matrix$binary_uncorr, c(1,6))
rate_matrix$binary_uncorr
rate_matrix$binary_uncorr <- corHMM::equateStateMatPars(rate_matrix$binary_uncorr, c(1,3))
rate_matrix$binary_uncorr
rate_matrix$binary_uncorr <- corHMM::equateStateMatPars(rate_matrix$binary_uncorr, c(1,4))
rate_matrix$binary_uncorr
rate_matrix$binary_uncorr <- corHMM::equateStateMatPars(rate_matrix$binary_uncorr, c(1,2))
rate_matrix$binary_uncorr

# corHMM analyses take ~1.5 hours to run on 3 cores (16 models total)
corHMM <- cache_RDS("results/corHMM_corr.rds", function(){
  #### run corHMM ####
  
  # run models of correlated and uncorrelated evolution, compare AICs to 
  # determine if evolution is correlated or not
  
  # store results from all corHMM analyses in a list
  corHMM <- list()
  
  # first run correlated and uncorrelated models with all 3 biomes, using
  # equal rate prior which as good as any other by my tests
  
  # run corHMM with no hidden states on two characters, will model correlated evolution
  # below runs ASR on wind/animal and 3 biomes (6 states)
  start_time <- Sys.time()
  corHMM$allbiomes_corr <- corHMM::corHMM(phy = tree_nowater, 
                                                   data = pollination_biome,
                                                   rate.cat = 1, 
                                                   rate.mat = rate_matrix$allbiomes_corr,
                                                   model = "ARD", 
                                                   nstarts = 10, n.cores = no_cores)
  corHMM$allbiomes_uncorr <- corHMM::corHMM(phy = tree_nowater, 
                                                   data = pollination_biome,
                                                   rate.cat = 1, 
                                                   rate.mat = rate_matrix$allbiomes_uncorr,
                                                   model = "ARD", 
                                                   nstarts = 10, n.cores = no_cores)
  end_time <- Sys.time()
  end_time - start_time
  # Time difference of 36.85346 mins
  
  # then run models on three binary biome definitions with 2 different root priors
  start_time <- Sys.time()
  for(name in names(pollbiome)){
    corHMM[[paste(name, "_corr_eqwt", sep = "")]] <- corHMM(phy = tree_nowater, 
                                                               data = pollbiome[[name]], 
                                                               rate.cat = 1, 
                                                               rate.mat = rate_matrix$binary_corr, 
                                                               node.states = "marginal",
                                                               root.p = "NULL",
                                                               nstarts = 10,
                                                               n.cores = no_cores)
    corHMM[[paste(name, "_uncorr_eqwt", sep = "")]] <- corHMM(phy = tree_nowater, 
                                                                 data = pollbiome[[name]], 
                                                                 rate.cat = 1, 
                                                                 rate.mat = rate_matrix$binary_uncorr, 
                                                                 node.states = "marginal",
                                                                 root.p = "NULL",
                                                                 nstarts = 10,
                                                                 n.cores = no_cores)
  }
  end_time <- Sys.time()
  end_time - start_time
  # Time difference of 55.42604 mins
  
  # save RDS of model output so I don't have to re-run these all the time!
  saveRDS(corHMM, file = "results/corHMM_corr.rds")
  
})

corHMM_results <- cache_RDS("results/corHMM_corr_results.csv", 
                            read_function = readr::read_csv,
                            save_function = write_csv,
                            function(){
  # export corHMM results by loop through list and export csv of states, rates and AIC
  corHMM_results <- data.frame()
  for (name in names(corHMM)){
    # export states as csv
    write.csv(x = corHMM[[name]]$states, 
              file = paste("results/biomes/", name, "_states.csv", sep = ""))
    # export rates as csv
    write.csv(x = corHMM[[name]]$solution,
              file = paste("results/biomes/", name, "_rates.csv", sep = ""))
    # assemble model fit data into df to export later
    results_row <- data.frame(model = paste(name, sep = ""), 
                              loglik = corHMM[[name]]$loglik,
                              AIC = corHMM[[name]]$AIC,
                              AICc = corHMM[[name]]$AICc)
    corHMM_results <- rbind(corHMM_results, results_row)
  }
  
  readr::write_csv(corHMM_results, "results/corHMM_corr_results.csv")
  rm(name, results_row)
})

# export print results from correlation models
sink(file = "results/pagels_biome_summaries.txt")
for(name in names(corHMM)){
  print(paste(name))
  print(corHMM[[name]])
  print(paste(" "))
}
sink(file = NULL)

# table of pollination mode vs superbiome occupancy
pollination_biome %>% dplyr::group_by(wind_animal, SB50) %>% dplyr::summarise(n())

rm(name, pollbiome, pollination_biome, rate_matrix, corHMM, corHMM_results, tree_nowater)

