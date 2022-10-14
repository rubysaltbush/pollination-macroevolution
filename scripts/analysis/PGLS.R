# run phylogenetic logistic regression on pollination mode (wind/animal) and
# mean latitude, temperature, rainfall and LAI as general tests of relatedness

#### prepare data and tree ####

# subset data to variables of interest for phylogenetic logistic regression
poll_pgls <- pollination1209 %>%
  dplyr::select(taxon_name, wind_animal, meanlat, abs_meanlat, meanMAT, meanTAP, meanLAI) %>%
  as.data.frame()

# drop water pollinated, polymorphic and missing data tips from tree
# lose 180 tips
to_drop <- poll_pgls %>%
  dplyr::filter(wind_animal %in% c("?", "2&4")|is.na(meanlat)|is.na(meanMAT)|
                  is.na(meanTAP)|is.na(meanLAI)|is.na(abs_meanlat))
tree_nowaterpolymissing <- ape::drop.tip(tree, to_drop$taxon_name)
rm(to_drop)
# and remove water pollinated & missing data taxa from morphological data
poll_pgls <- poll_pgls %>%
  dplyr::filter(wind_animal %in% c("2", "4") & !is.na(meanlat) & !is.na(meanMAT) 
                & !is.na(meanTAP) & !is.na(meanLAI))

# taxon_name to row names
rownames(poll_pgls) <- poll_pgls[,1]
poll_pgls[,1] <- NULL

# redefine 2 (wind poll) and 4 (animal poll) as 0 and 1
poll_pgls$wind_animal <- gsub("2", "0", poll_pgls$wind_animal)
poll_pgls$wind_animal <- gsub("4", "1", poll_pgls$wind_animal)
table(poll_pgls$wind_animal)
# 114 wind pollinated taxa to 907 animal pollinated taxa

# double check distribution of continuous variables
plot(poll_pgls$meanlat)
plot(poll_pgls$abs_meanlat) # more outlying high values??
plot(poll_pgls$meanMAT) # few outlying low values
plot(poll_pgls$meanTAP) # few outlying high values
plot(poll_pgls$meanLAI) # few outlying high values

hist(poll_pgls$meanlat)
hist(poll_pgls$abs_meanlat) # few outlying high values
hist(poll_pgls$meanMAT) # few outlying low values
hist(poll_pgls$meanTAP) # few outlying high values
hist(poll_pgls$meanLAI) # few outlying high values

# looks okay, don't think transformations necessary

# plot MAT vs latitude, TAP vs LAI 
# (as reason for why not reporting MAT/TAP in main text)
plot(poll_pgls$meanlat ~ poll_pgls$meanMAT)

# first run linear models
lm <- list()
lm$lat_MAT <- lm(poll_pgls$meanMAT ~ poll_pgls$abs_meanlat)
lm$LAI_TAP <- lm(poll_pgls$meanTAP ~ poll_pgls$meanLAI)

pdf(file = "figures/latitude_vs_MAT.pdf", width = 8, height = 5)
par(las = 1, bty = "l") # remove plot outline
plot(poll_pgls$meanMAT ~ poll_pgls$abs_meanlat,
     xlab = "Species mean absolute latitude",
     ylab = "Species mean Mean Annual Temperature (ºC)")
# add line from regression
abline(lm$lat_MAT, col = "red") 
# add R2 and p-value from regression as plot title
title(main = paste("R² = ", signif(summary(lm$lat_MAT)$r.squared, 2),
                   "    P = ", 
                   format.pval(summary(lm$lat_MAT)$coef[2,4], 
                               eps = .001, digits = 2))) 
dev.off()

pdf(file = "figures/LAI_vs_TAP.pdf", width = 8, height = 5)
par(las = 1, bty = "l") # remove plot outline
plot(poll_pgls$meanTAP ~ poll_pgls$meanLAI,
     xlab = "Species mean Leaf Area Index",
     ylab = "Species mean Total Annual Precipitation (mm)")
abline(lm$LAI_TAP, col = "red") # add line from regression
# add R2 and p-value from regression as plot title
title(main = paste("R² = ", signif(summary(lm$LAI_TAP)$r.squared, 2),
                   "    P = ", 
                   format.pval(summary(lm$LAI_TAP)$coef[2,4], 
                               eps = .001, digits = 2))) 
dev.off()

rm(lm)

#### run model ####
# below adapted from Joly and Schoen (2021)
# works without polymorphic or missing data
# Model fit with Ives and Garlan optimisation

PGLS <- list()
PGLS$lat <- phylolm::phyloglm(wind_animal ~ meanlat, 
                          data = poll_pgls, 
                          phy = tree_nowaterpolymissing,
                          method = "logistic_IG10", 
                          boot = 100)

summary(PGLS$lat)
# NOT SIGNIFICANT!

PGLS$abs_lat <- phylolm::phyloglm(wind_animal ~ abs_meanlat, 
                              data = poll_pgls, 
                              phy = tree_nowaterpolymissing,
                              method = "logistic_IG10", 
                              boot = 100)

summary(PGLS$abs_lat)
# significant!!! makes sense if absolute lat values more linear

PGLS$MAT <- phylolm::phyloglm(wind_animal ~ meanMAT, 
                              data = poll_pgls, 
                              phy = tree_nowaterpolymissing,
                              method = "logistic_IG10", 
                              boot = 100)

summary(PGLS$MAT)
# NOT SIGNIFICANT! strange given linear relationship with abs latitude values

PGLS$TAP <- phylolm::phyloglm(wind_animal ~ meanTAP, 
                              data = poll_pgls, 
                              phy = tree_nowaterpolymissing,
                              method = "logistic_IG10", 
                              boot = 100)
summary(PGLS$TAP)
# JUST SIGNIFICANT
# maybe this makes sense with my almost-support for correlated evolution in arid/non-arid biome

PGLS$LAI <- phylolm::phyloglm(wind_animal ~ meanLAI, 
                              data = poll_pgls, 
                              phy = tree_nowaterpolymissing,
                              method = "logistic_IG10", 
                              boot = 100)
summary(PGLS$LAI)
# JUST SIGNIFICANT

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
     pch = 2, cex = 1, 
     col = factor(wind_animal),
     xlab = "Species mean latitude (absolute)", 
     ylab = "Species mean Leaf Area Index")
legend("topright",
       legend = c("animal pollinated (1)", "wind pollinated (0)"),
       col = c(my_colours$wind_water_animal[3], my_colours$wind_water_animal[1]),
       pch = 2, bty = "n")
dev.off()


# mean latitude (no significant relationship, don't plot line)
pdf(file="figures/windanimal_meanlat.pdf", width = 8, height = 5)
par(las = 1, bty = "l") # remove plot outline
plot(wind_animal ~ meanlat, 
     data = poll_pgls,
     pch = 2, cex = 1, 
     col = factor(wind_animal),
     xlab = "Mean latitude", ylab = "Pollination mode",
     xlim = c(-65, 65))
dev.off()

# mean absolute latitude
pdf(file="figures/windanimal_absmeanlat.pdf", width = 8, height = 5)
par(las = 1, bty = "l") # remove plot outline
plot(wind_animal ~ abs_meanlat, 
     data = poll_pgls,
     pch = 2, cex = 1, 
     col = factor(wind_animal),
     xlab = "Mean latitude (absolute)", ylab = "Pollination mode",
     xlim = c(0, 65))
cc <- coef(PGLS$abs_lat)
curve(plogis(cc[1] + cc[2] * x), col = "red", add = TRUE)
dev.off()
rm(cc)

# mean MAT (no significant relationship, don't plot line)
pdf(file="figures/windanimal_meanMAT.pdf", width = 8, height = 5)
par(las = 1, bty = "l")  # remove plot outline
plot(wind_animal ~ meanMAT, 
     data = poll_pgls,
     pch = 2, cex = 1, 
     col = factor(wind_animal),
     xlab = "Average Mean Annual Temperature (ºC)", ylab = "Pollination mode",
     xlim = c(-6, 30))
dev.off()

# mean TAP
pdf(file="figures/windanimal_meanTAP.pdf", width = 8, height = 5)
par(las = 1, bty = "l")  # remove plot outline
plot(wind_animal ~ meanTAP, 
     data = poll_pgls,
     pch = 2, cex = 1, 
     col = factor(wind_animal),
     xlab = "Mean Total Annual Precipitation (mm)", ylab = "Pollination mode",
     xlim = c(15, 4000))
cc <- coef(PGLS$TAP)
curve(plogis(cc[1] + cc[2] * x), col = "red", add = TRUE)
dev.off()
rm(cc)
# line shows not a hugely strong relationship, but likelihood of animal pollination
# increases with increasing TAP

# mean LAI
pdf(file="figures/windanimal_meanLAI.pdf", width = 8, height = 5)
par(las = 1, bty = "l") # remove plot outline
plot(wind_animal ~ meanLAI, 
     data = poll_pgls,
     pch = 2, cex = 1, 
     col = factor(wind_animal),
     xlab = "Mean Leaf Area Index", ylab = "Pollination mode",
     xlim = c(0, 80))
cc <- coef(PGLS$LAI)
curve(plogis(cc[1] + cc[2] * x), col = "red", add = TRUE)
dev.off()
rm(cc)
# line shows not a hugely strong relationship, but likelihood of animal pollination
# increases with increasing LAI

palette("default") # set colour palette back to default
rm(PGLS, PGLS_results, poll_pgls, tree_nowaterpolymissing)

#### phylogenetic logistic regression for binary LAI/latitude categories ####

## prepare data

pollbinary_pgls <- pollination1209 %>%
  dplyr::select(taxon_name, wind_animal, trop_nontrop, open_closed) %>%
  as.data.frame()

# drop water pollinated, polymorphic and missing data tips from tree
# lose 180 tips
to_drop <- pollbinary_pgls %>%
  dplyr::filter(wind_animal %in% c("?", "2&4")|is.na(trop_nontrop)|is.na(open_closed))
tree_nowaterpolymissingbinary <- ape::drop.tip(tree, to_drop$taxon_name)
rm(to_drop)
# and remove same taxa from morphological data
pollbinary_pgls <- pollbinary_pgls %>%
  dplyr::filter(wind_animal %in% c("2", "4") & trop_nontrop %in% c("tropical", "nontropical") & open_closed %in% c(1,2))

# taxon_name to row names
rownames(pollbinary_pgls) <- pollbinary_pgls[,1]
pollbinary_pgls[,1] <- NULL
# redefine 2 (wind poll) and 4 (animal poll) as 0 and 1
pollbinary_pgls$wind_animal <- gsub("2", "0", pollbinary_pgls$wind_animal)
pollbinary_pgls$wind_animal <- gsub("4", "1", pollbinary_pgls$wind_animal)
table(pollbinary_pgls$wind_animal)
# 114 wind, 907 animal

# convert trop_nontrop to binary integers - 1 is tropical, 0 nontropical
pollbinary_pgls$trop_nontrop <- gsub("nontropical", "0", pollbinary_pgls$trop_nontrop)
pollbinary_pgls$trop_nontrop <- gsub("tropical", "1", pollbinary_pgls$trop_nontrop)
table(pollbinary_pgls$trop_nontrop)

# convert open_closed to binary integers - now 0 is open, 1 is closed
pollbinary_pgls$open_closed <- gsub("2", "0", pollbinary_pgls$open_closed)
table(pollbinary_pgls$open_closed)

## run models

PGLS$trop_nontrop <- phylolm::phyloglm(wind_animal ~ trop_nontrop, 
                                    data = pollbinary_pgls, 
                                    phy = tree_nowaterpolymissingbinary,
                                    method = "logistic_IG10", 
                                    boot = 100)
summary(PGLS$trop_nontrop)

PGLS$open_closed <- phylolm::phyloglm(wind_animal ~ open_closed, 
                                     data = pollbinary_pgls, 
                                     phy = tree_nowaterpolymissingbinary,
                                     method = "logistic_IG10", 
                                     boot = 100)
summary(PGLS$open_closed)

# none of these significant! suspect we just don't have enough data points for
# broad categories like this, thus why continuous significant but cat not

rm(pollbinary_pgls, tree_nowaterpolymissingbinary)

#### phylogenetic logistic regression for biomes ####

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

# add arid/non-arid, temperate/non-temp and tropical/non-trop columns to data
# convert superbiome to tropical/extra-tropical by replacing "arid" (2) with 0 (non-tropical)
pollbiome_pgls$SB50_trop <- gsub("2", "0", pollbiome_pgls$SB50)
table(pollbiome_pgls$SB50_trop)

# convert superbiome to arid/extra-arid by replacing "tropical" (1) with 0 (non-arid)
pollbiome_pgls$SB50_arid <- gsub("1", "0", pollbiome_pgls$SB50)
pollbiome_pgls$SB50_arid <- gsub("2", "1", pollbiome_pgls$SB50_arid)
table(pollbiome_pgls$SB50_arid)

# convert superbiome to temperate/non-temperate
pollbiome_pgls$SB50_temp <- gsub("0", "3", pollbiome_pgls$SB50)
pollbiome_pgls$SB50_temp <- gsub("2", "0", pollbiome_pgls$SB50_temp)
pollbiome_pgls$SB50_temp <- gsub("1", "0", pollbiome_pgls$SB50_temp)
pollbiome_pgls$SB50_temp <- gsub("3", "1", pollbiome_pgls$SB50_temp)
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
