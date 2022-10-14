# run corHMM Ancestral State Reconstructions of pollination modes in angiosperms
# this will take some time to run depending on number of computer cores


#####Ancestral State Reconstructions#####
ASR <- cache_RDS("results/ASR.rds", function(){
  
  # make matrices for ASR
  matrices <- list()
  # corHMM::rayDISC needs traits as matrix
  matrices$abiotic_animal <- pollination1209 %>%
    dplyr::select(taxon_name, value = abiotic_animal) %>%
    as.matrix()
  
  matrices$abiotic_animal_nopoly <- pollination1209 %>%
    dplyr::select(taxon_name, value = abiotic_animal_nopoly) %>%
    as.matrix()
  
  matrices$wind_water_animal <- pollination1209 %>%
    dplyr::select(taxon_name, value = wind_water_animal) %>%
    as.matrix()
  
  matrices$wind_water_animal_nopoly <- pollination1209 %>%
    dplyr::select(taxon_name, value = wind_water_animal_nopoly) %>%
    as.matrix()
  
  matrices$abiotic_vert_insect <- pollination1209 %>%
    dplyr::select(taxon_name, value = abiotic_vert_insect) %>%
    as.matrix()
  
  matrices$abiotic_vert_insect_nopoly <- pollination1209 %>%
    dplyr::select(taxon_name, value = abiotic_vert_insect_nopoly) %>%
    as.matrix()
  
  matrices$wind_water_vert_insect <- pollination1209 %>%
    dplyr::select(taxon_name, value = wind_water_vert_insect) %>%
    as.matrix()
  
  # not running ASR on wind_water_vert_insect with no polymorphisms for now as 
  # this overfits models (way too many parameters)
  
  # store results of ancestral state reconstructions in a list
  ASR <- list()
  
  # run Equal Rates and All Rates Different ASR
  # corHMM default root prior is yang
  start_time <- Sys.time()
  # run in parallel on as many cores as available
  # first ER models
  ASR <- foreach(name = names(matrices)) %dopar% {
    result <- corHMM::corHMM(tree, matrices[[name]], 
                         model = "ER", 
                         rate.cat = 1, 
                         nstarts = 10, 
                         n.cores = 1)
    result
  }
  names(ASR) <- paste(names(matrices), "_ER", sep = "") # name ER results
  # then ARD models
  ASR_ARD <- foreach(name = names(matrices)) %dopar% {
    result <- corHMM::corHMM(tree, matrices[[name]], 
                             model = "ARD", 
                             rate.cat = 1, 
                             nstarts = 10, 
                             n.cores = 1)
    result
  }
  names(ASR_ARD) <- paste(names(matrices), "_ARD", sep = "") # name ARD results
  # combine models list of all results
  ASR <- c(ASR, ASR_ARD)
  rm(ASR_ARD)
  end_time <- Sys.time()
  end_time - start_time
  # time elapsed = 
  
  # output printout of different models to text file
  sink(file = "results/ASR_descriptions.txt")
  for(name in names(ASR)){
    print(paste(name))
    print(ASR[[name]])
    print(paste("states"))
    print(head(ASR[[name]]$states))
    print(paste(" "))
  }
  sink(file = NULL)
  rm(name)
  
  # look at results
  ASR$abiotic_animal_ER
  head(ASR$abiotic_animal_ER$states)
  # look at results
  ASR$abiotic_animal_ARD
  head(ASR$abiotic_animal_ARD$states)
  # with all rates different it looks like transition rates from abiotic to animal
  # are slightly higher, AICc slightly higher too though
  
  # look at results
  ASR$abiotic_animal_nopoly_ER
  head(ASR$abiotic_animal_nopoly_ER$states)
  # look at results
  ASR$abiotic_animal_nopoly_ARD
  head(ASR$abiotic_animal_nopoly_ARD$states)
  # ARD has lower AICc, suggests that polymorphic state is more 
  # frequently transitional for animal->wind pollination? 
  
  # look at results
  ASR$wind_water_animal_ER
  head(ASR$wind_water_animal_ER$states)
  # look at results
  ASR$wind_water_animal_ARD
  head(ASR$wind_water_animal_ARD$states)
  # look at results
  ASR$wind_water_animal_nopoly_ER
  head(ASR$wind_water_animal_nopoly_ER$states)
  # look at results
  ASR$wind_water_animal_nopoly_ARD
  head(ASR$wind_water_animal_nopoly_ARD$states)
  # look at results
  ASR$abiotic_vert_insect_ER
  head(ASR$abiotic_vert_insect_ER$states)
  # look at results
  ASR$abiotic_vert_insect_ARD
  head(ASR$abiotic_vert_insect_ARD$states)
  # look at results
  ASR$abiotic_vert_insect_nopoly_ER
  head(ASR$abiotic_vert_insect_nopoly_ER$states)
  # look at results
  ASR$abiotic_vert_insect_nopoly_ARD
  head(ASR$abiotic_vert_insect_nopoly_ARD$states)
  # look at results
  ASR$wind_water_vert_insect_ER
  head(ASR$wind_water_vert_insect_ER$states)
  # look at results
  ASR$wind_water_vert_insect_ARD
  head(ASR$wind_water_vert_insect_ARD$states)
  
  
  # export ASR results by loop through list and export csv of states, rates and AIC
  ASR_results <- data.frame()
  for (name in names(ASR)){
    # export states as csv
    write.csv(x = ASR[[name]]$states, 
              file = paste("results/", name, "_states.csv", sep = ""))
    # export rates as csv
    write.csv(x = ASR[[name]]$solution,
              file = paste("results/", name, "_rates.csv", sep = ""))
    # assemble model fit data into df to export later
    results_row <- data.frame(model = paste(name, sep = ""), 
                              loglik = ASR[[name]]$loglik,
                              AIC = ASR[[name]]$AIC,
                              AICc = ASR[[name]]$AICc,
                              root_prior = ASR[[name]]$root.p)
    ASR_results <- rbind(ASR_results, results_row)
  }
  
  readr::write_csv(ASR_results, "results/ASR_results.csv")
  rm(name, results_row, matrices, ASR_results)
  
  # save RDS of model output so I don't have to re-run these all the time!
  saveRDS(ASR, file = "results/ASR.rds")
})
