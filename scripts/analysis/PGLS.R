# run phylogenetic logistic regression on pollination mode (wind/animal) vs.
# absolute mean latitude and LAI as test of relatedness

#### prepare data and tree ####

# subset data to variables of interest for phylogenetic logistic regression
poll_pgls <- pollination1209 %>%
  dplyr::select(taxon_name, wind_animal, abs_meanlat, meanLAI) %>%
  as.data.frame()

# drop water pollinated, polymorphic and missing data tips from tree
# lose 180 tips
to_drop <- poll_pgls %>%
  dplyr::filter(wind_animal %in% c("?", "2&4")|is.na(meanLAI)|is.na(abs_meanlat))
tree_nowaterpolymissing <- ape::drop.tip(tree, to_drop$taxon_name)
rm(to_drop)
# and remove water pollinated & missing data taxa from morphological data
poll_pgls <- poll_pgls %>%
  dplyr::filter(wind_animal %in% c("2", "4") & !is.na(abs_meanlat) & !is.na(meanLAI))

# taxon_name to row names
rownames(poll_pgls) <- poll_pgls[,1]
poll_pgls[,1] <- NULL

# redefine 2 (wind poll) and 4 (animal poll) as 0 and 1
poll_pgls$wind_animal <- gsub("2", "0", poll_pgls$wind_animal)
poll_pgls$wind_animal <- gsub("4", "1", poll_pgls$wind_animal)
table(poll_pgls$wind_animal)
# 114 wind pollinated taxa to 908 animal pollinated taxa

# double check distribution of continuous variables
plot(poll_pgls$abs_meanlat) # some outlying high values
plot(poll_pgls$meanLAI) # few outlying high values

hist(poll_pgls$abs_meanlat) # few outlying high values
hist(poll_pgls$meanLAI) # few outlying high values

# looks okay, don't think transformations necessary

#### run model ####
# below adapted from Joly and Schoen (2021)
# works without polymorphic or missing data
# Model fit with Ives and Garlan optimisation

PGLS <- list()

PGLS$abs_lat <- phylolm::phyloglm(wind_animal ~ abs_meanlat, 
                              data = poll_pgls, 
                              phy = tree_nowaterpolymissing,
                              method = "logistic_IG10", 
                              boot = 100)

summary(PGLS$abs_lat)

PGLS$LAI <- phylolm::phyloglm(wind_animal ~ meanLAI, 
                              data = poll_pgls, 
                              phy = tree_nowaterpolymissing,
                              method = "logistic_IG10", 
                              boot = 100)
summary(PGLS$LAI)

# export summary of model output
sink(file = "results/phylologis_summaries.txt")
for(name in names(PGLS)){
  print(paste(name))
  print(PGLS[[name]])
  print(paste(" "))
}
sink(file = NULL)

# export summary of results to table
PGLS_results <- data.frame()
for(name in names(PGLS)){
  results_row <- data.frame(coef(summary(PGLS[[name]])))
  PGLS_results <- rbind(PGLS_results, results_row)
}
PGLS_results$name <- rownames(PGLS_results)
readr::write_csv(PGLS_results, "results/PGLS_results.csv")
rm(name, results_row)

#### visualise ####

# plot LAI vs absolute latitude coloured by pollination mode
pdf(file="figures/meanLAI_meanlat_windanimal.pdf", width = 8, height = 5)
par(las = 1, bty = "l") # remove plot outline
palette(c("#F6BA65", "#BD8DAA"))# set colour palette
plot(meanLAI ~ abs_meanlat, 
     data = poll_pgls, 
     pch = 16, cex = 1, 
     col = factor(wind_animal),
     xlab = "Species mean latitude (absolute)", 
     ylab = "Species mean Leaf Area Index")
legend("topright",
       legend = c("animal pollinated", "wind pollinated"),
       col = c(my_colours$wind_water_animal[3], my_colours$wind_water_animal[1]),
       pch = 16, bty = "n")
dev.off()

# mean absolute latitude
pdf(file="figures/windanimal_absmeanlat.pdf", width = 8, height = 5)
par(las = 1, bty = "l") # remove plot outline
plot(wind_animal ~ abs_meanlat, 
     data = poll_pgls,
     pch = 1, cex = 1, 
     col = factor(wind_animal),
     xlab = "Mean latitude (absolute)", ylab = "Probability of animal pollination",
     xlim = c(0, 65))
cc <- coef(PGLS$abs_lat)
curve(plogis(cc[1] + cc[2] * x), col = "red", add = TRUE)
legend("right",
       legend = c("animal pollinated", "wind pollinated"),
       col = c(my_colours$wind_water_animal[3], my_colours$wind_water_animal[1]),
       pch = 1, bty = "n")
dev.off()
rm(cc)

# mean LAI
pdf(file="figures/windanimal_meanLAI.pdf", width = 8, height = 5)
par(las = 1, bty = "l") # remove plot outline
plot(wind_animal ~ meanLAI, 
     data = poll_pgls,
     pch = 1, cex = 1, 
     col = factor(wind_animal),
     xlab = "Mean Leaf Area Index", ylab = "Probability of animal pollination",
     xlim = c(0, 8))
cc <- coef(PGLS$LAI)
curve(plogis(cc[1] + cc[2] * x), col = "red", add = TRUE)
dev.off()
rm(cc)
# line shows likelihood of animal pollination increases with increasing LAI

palette("default") # set colour palette back to default
rm(PGLS, PGLS_results, poll_pgls, tree_nowaterpolymissing)


#### phylogenetic logistic regression for categorical biomes ####
# using superbiome classification from R-B 2020

## prepare data
pollbiome_pgls <- pollination1209 %>%
  dplyr::select(taxon_name, wind_animal, SB50) %>%
  as.data.frame()

# drop water pollinated, polymorphic and missing data tips from tree
# lose 189 tips
to_drop <- pollbiome_pgls %>%
  dplyr::filter(wind_animal %in% c("?", "2&4")|SB50 %in% ("?"))
tree_nowaterpolymissing2 <- ape::drop.tip(tree, to_drop$taxon_name)
rm(to_drop)
# and remove same taxa from morphological data
pollbiome_pgls <- pollbiome_pgls %>%
  dplyr::filter(wind_animal %in% c("2", "4") & SB50 %in% c("0", "1", "2"))

# taxon_name to row names
rownames(pollbiome_pgls) <- pollbiome_pgls[,1]
pollbiome_pgls[,1] <- NULL
# redefine 2 (wind poll) and 4 (animal poll) as 0 and 1
pollbiome_pgls$wind_animal <- gsub("2", "0", pollbiome_pgls$wind_animal)
pollbiome_pgls$wind_animal <- gsub("4", "1", pollbiome_pgls$wind_animal)
table(pollbiome_pgls$wind_animal)
# 108 wind, 904 animal

# order superbiomes according to hypothesis, that animal pollination most likely
# in tropical superbiome and least likely in arid superbiome
# thus tropical = 0, temperate = 1 and arid = 2
pollbiome_pgls$SB50 <- gsub("1", "3", pollbiome_pgls$SB50)
pollbiome_pgls$SB50 <- gsub("0", "1", pollbiome_pgls$SB50)
pollbiome_pgls$SB50 <- gsub("3", "0", pollbiome_pgls$SB50)

# add arid/non-arid, temperate/non-temp and tropical/non-trop columns to data
# convert superbiome to tropical/extra-tropical by replacing "arid" (2) with 0 (non-tropical)
pollbiome_pgls$SB50_trop <- gsub("0", "3", pollbiome_pgls$SB50)
pollbiome_pgls$SB50_trop <- gsub("2", "0", pollbiome_pgls$SB50_trop)
pollbiome_pgls$SB50_trop <- gsub("1", "0", pollbiome_pgls$SB50_trop)
pollbiome_pgls$SB50_trop <- gsub("3", "1", pollbiome_pgls$SB50_trop)
table(pollbiome_pgls$SB50_trop)

# convert superbiome to arid/extra-arid by replacing "temperate" (1) with 0 (non-arid)
pollbiome_pgls$SB50_arid <- gsub("1", "0", pollbiome_pgls$SB50)
pollbiome_pgls$SB50_arid <- gsub("2", "1", pollbiome_pgls$SB50_arid)
table(pollbiome_pgls$SB50_arid)

# convert superbiome to temperate/non-temperate
pollbiome_pgls$SB50_temp <- gsub("2", "1", pollbiome_pgls$SB50)
table(pollbiome_pgls$SB50_temp)

## run models
PGLS_biome <- list()
PGLS_biome$biome_all <- phylolm::phyloglm(wind_animal ~ SB50, 
                          data = pollbiome_pgls, 
                          phy = tree_nowaterpolymissing2,
                          method = "logistic_IG10", 
                          boot = 100)
summary(PGLS_biome$biome_all)

PGLS_biome$biome_trop <- phylolm::phyloglm(wind_animal ~ SB50_trop, 
                                    data = pollbiome_pgls, 
                                    phy = tree_nowaterpolymissing2,
                                    method = "logistic_IG10", 
                                    boot = 100)
summary(PGLS_biome$biome_trop)

PGLS_biome$biome_arid <- phylolm::phyloglm(wind_animal ~ SB50_arid, 
                                     data = pollbiome_pgls, 
                                     phy = tree_nowaterpolymissing2,
                                     method = "logistic_IG10", 
                                     boot = 100)
summary(PGLS_biome$biome_arid)

PGLS_biome$biome_temp <- phylolm::phyloglm(wind_animal ~ SB50_temp, 
                                     data = pollbiome_pgls, 
                                     phy = tree_nowaterpolymissing2,
                                     method = "logistic_IG10", 
                                     boot = 100)
summary(PGLS_biome$biome_temp)

# export summary of model output
sink(file = "results/phylologis_biome_summaries.txt")
for(name in names(PGLS_biome)){
  print(paste(name))
  print(PGLS_biome[[name]])
  print(paste(" "))
}
sink(file = NULL)

# export summary of results to table
PGLS_biome_results <- data.frame()
for(name in names(PGLS_biome)){
  results_row <- data.frame(coef(summary(PGLS_biome[[name]])))
  PGLS_biome_results <- rbind(PGLS_biome_results, results_row)
}
PGLS_biome_results$name <- rownames(PGLS_biome_results)
readr::write_csv(PGLS_biome_results, "results/PGLS_biome_results.csv")
rm(name, results_row)

# none of these significant! ach well.

rm(pollbiome_pgls, tree_nowaterpolymissing2, PGLS_biome, PGLS_biome_results)
