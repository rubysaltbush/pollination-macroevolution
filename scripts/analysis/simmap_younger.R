# this script repeats the wind-animal and insect-vertebrate stochastic
# character mapping on an alternative (younger) phylogeny from Ramírez-Barahona
# et al. (2020) to consider alternative transition times

#####Stochastic character mapping#####

### set up environment ###
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

# function to cache pre-prepared R data. If RDS already in cache will read data
source("scripts/functions/cache_RDS.R")

# for graphs assign fixed colours
source("scripts/my_colours.R")

#### prepare data ####

# read in prepped pollination data
source("scripts/pollination1209.R")

# read in alternative angiosperm phylogeny, with all the same tips as the main
# phylogeny used in analysis but with root (crown node) dated at 154 mya instead
# of 196 mya

tree_younger <- ape::read.tree("data_input/eFLOWER-1209_CC_complete_MCC.phy")
# double check root age of tree
max(nodeHeights(tree_younger))
# 154 mya, 40 million years younger than tree in main analyses

# prepare Ancestral State Reconstructions on younger tree
# drop tips of tree where binary characters missing
tree_wind_animal <- ape::drop.tip(tree_younger, pollination1209$taxon_name[pollination1209$wind_animal == "?"])
tree_vert_insect <- ape::drop.tip(tree_younger, pollination1209$taxon_name[pollination1209$vert_insect == "?"])
rm(tree_younger)

# prepare matrices with extra character states dropped
matrices <- list()

matrices$wind_animal <- pollination1209 %>%
  dplyr::select(taxon_name, value = wind_animal) %>%
  dplyr::filter(value != "?") %>%
  as.matrix()

matrices$vert_insect <- pollination1209 %>%
  dplyr::select(taxon_name, value = vert_insect) %>%
  dplyr::filter(value != "?") %>%
  as.matrix()

# run ARD ASR for these two binary categorisations using modified trees
ASR_forsimmap_young <- list()
ASR_forsimmap_young$wind_animal_ARD <- corHMM::corHMM(tree_wind_animal, 
                                                matrices$wind_animal, 
                                                model = "ARD", 
                                                rate.cat = 1, 
                                                nstarts = 10,
                                                n.cores = no_cores)
ASR_forsimmap_young$wind_animal_ER <- corHMM::corHMM(tree_wind_animal, 
                                               matrices$wind_animal, 
                                               model = "ER", 
                                               rate.cat = 1, 
                                               nstarts = 10,
                                               n.cores = no_cores)
ASR_forsimmap_young$vert_insect_ARD <- corHMM::corHMM(tree_vert_insect, 
                                                matrices$vert_insect, 
                                                model = "ARD", 
                                                rate.cat = 1, 
                                                nstarts = 10,
                                                n.cores = no_cores)
ASR_forsimmap_young$vert_insect_ER <- corHMM::corHMM(tree_vert_insect, 
                                               matrices$vert_insect, 
                                               model = "ER", 
                                               rate.cat = 1, 
                                               nstarts = 10,
                                               n.cores = no_cores)
# select most supported model from ER or ARD for stochastic mapping?
ASR_forsimmap_young$wind_animal_ARD
ASR_forsimmap_young$wind_animal_ER
# <2 (i.e. no) difference in AICc, ARD transition rate from wind to animal higher
ASR_forsimmap_young$vert_insect_ARD
ASR_forsimmap_young$vert_insect_ER
# 30 difference in AICc, ARD lower. Transition rate from vert to insect higher
# interesting that these results identical to those on main tree! makes sense
# given only dating changed

# remove ER models from list
ASR_forsimmap_young <- ASR_forsimmap_young[-c(2,4)]
# remove items no longer needed
rm(matrices, pollination1209, cache_RDS)

#### run simmap analyses ####

# for loop over ASRs with binary states, run simmapping for each model using corHMMM

# list to store results of simmaps
simmaps <- list()

# run simmaps in parallel
start_time <- Sys.time()
simmaps <- foreach(name = names(ASR_forsimmap_young)) %dopar% {
  phy <- ASR_forsimmap_young[[name]]$phy
  data <- ASR_forsimmap_young[[name]]$data
  data[is.na(data)] <- "?"
  model <- ASR_forsimmap_young[[name]]$solution
  model[is.na(model)] <- 0
  diag(model) <- -rowSums(model)
  result <- corHMM::makeSimmap(
    tree = phy,
    data = data,
    model = model,
    rate.cat = 1,
    nSim = 1000,
    nCores = floor(no_cores/length(ASR_forsimmap_young))
  )
  class(result) <- c("multiSimmap", "multiPhylo")
  result
}
names(simmaps) <- names(ASR_forsimmap_young)
end_time <- Sys.time()
end_time - start_time
# Time difference of 3 mins for 2*1000 sims on 3 cores
rm(end_time, start_time, ASR_forsimmap_young)

# describe and calculate 95% confidence intervals etc for # of transitions in simmaps
# for SOME REASON rewriting below as foreach makes it take longer and use more memory
start_time <- Sys.time()
simmap_descriptions <- list()
simmap_density <- list()
for(name in names(simmaps)){
  simmap_descriptions[[name]] <- phytools::describe.simmap(simmaps[[name]])
  simmap_density[[name]] <- phytools::density.multiSimmap(simmaps[[name]])
}
end_time <- Sys.time()
end_time - start_time
# Time difference of 1 min
rm(name, start_time, end_time)

# check out results!
simmap_descriptions$wind_animal_ARD
simmap_density$wind_animal_ARD
simmap_descriptions$vert_insect_ARD
simmap_density$vert_insect_ARD

# still more transitions on average between insect and vertebrate pollination (89)
# than between wind and animal (53)

# export print results from simmapping to text file
sink(file = "results/simmap_descriptions_younger.txt")
for(name in names(simmap_descriptions)){
  print(paste(name))
  print(simmap_descriptions[[name]])
  print(paste(" "))
  print(simmap_density[[name]])
  print(paste(" "))
}
sink(file = NULL)
rm(name, simmap_density)

#### graphing simmaps ####

# density maps for each binary multisimmap - each one takes >10 minutes to build
density_maps <- list()
# wind and animal
# to change colours in density map, first build density map as object TAKES AGES
density_maps$wind_animal_ARD <- phytools::densityMap(simmaps$wind_animal_ARD, plot = FALSE)
# vertebrate and insect
# to change colours in density map, first build density map as object TAKES AGES
density_maps$vert_insect_ARD <- phytools::densityMap(simmaps$vert_insect_ARD, plot = FALSE)

# wind and animal
# assign new colours to density map (taken from phytools blog)
n <- length(density_maps$wind_animal_ARD$cols)
density_maps$wind_animal_ARD$cols[1:n] <- colorRampPalette(c(my_colours$wind_water_animal[1], my_colours$wind_water_animal[3]), space="Lab")(n)

# then plot it and export to pdf
pdf(file = "figures/densitymap_wind_animal_ARD_1000_tall_younger.pdf", width = 12.6, height = 20, useDingbats = FALSE)
plot(density_maps$wind_animal_ARD, ftype = "off", lwd = 4, xlim = c(-5, 205))

# draw time axis 
# adapted from http://blog.phytools.org/2018/02/another-technique-for-including-time.html
T <- max(nodeHeights(tree_wind_animal))
obj <- axis(1, pos = -2, at = seq(T, -10, by = -40), cex.axis = 3,
            labels = FALSE, lwd = 2)
axis(side = 1, pos = -2, at = seq(T, -10, by = -40), cex.axis = 3,
     labels = FALSE, lwd = 2)
text(x = obj, y = rep(-26, length(obj)), T - obj, cex = 2)
text(x = mean(obj), y = -53, "time (mya)", cex = 2.1)
rm(T, obj)

# label clades
#ANA
segments(x0 = 158, y0 = 1, x1 = 158, y1 = 10, lwd = 10, col = "#bdbdbd", lend = "butt")
text(x = 158+5, y = 5.5, "ANA", srt = 0, cex = 2, col = "#bdbdbd", adj = 0)
#Magnoliids
segments(x0 = 158, y0 = 11, x1 = 158, y1 = 63, lwd = 10, col = "#636363", lend = "butt")
text(x = 158+5, y = 37, "Magnoliids", srt = 0, cex = 2, col = "#636363", adj = 0)
#Monocots
segments(x0 = 158, y0 = 66, x1 = 158, y1 = 296, lwd = 10, col = "#bdbdbd", lend = "butt")
text(x = 158+5, y = 181, "Monocots", srt = 0, cex = 2, col = "#bdbdbd", adj = 0)
#Commelinids
segments(x0 = 158+2, y0 = 193, x1 = 158+2, y1 = 296, lwd = 10, col = "#636363", lend = "butt")
text(x = 158+5, y = 244.5, "Commelinids", srt = 0, cex = 2, col = "#636363", adj = 0)
#Eudicots
segments(x0 = 158, y0 = 297, x1 = 158, y1 = 1145, lwd = 10, col = "black", lend = "butt")
text(x = 158+5, y = 721, "Eudicots", srt = 0, cex = 2, col = "black", adj = 0)
#Rosids
segments(x0 = 158+2, y0 = 368, x1 = 158+2, y1 = 680, lwd = 10, col = "#bdbdbd", lend = "butt")
text(x = 158+5, y = 524, "Rosids", srt = 0, cex = 2, col = "#bdbdbd", adj = 0)
#Asterids
segments(x0 = 158+2, y0 = 806, x1 = 158+2, y1 = 1145, lwd = 10, col = "#bdbdbd", lend = "butt")
text(x = 158+5, y = 975.5, "Asterids", srt = 0, cex = 2, col = "#bdbdbd", adj = 0)
#Chloranthales
segments(x0 = 158, y0 = 64, x1 = 158, y1 = 65, lwd = 10, col = "black", lend = "butt")
text(x = 158+5, y = 64.5, "Chloranthales", srt = 0, cex = 2, col = "black", adj = 0)

dev.off()

# vertebrate and insect
# assign new colours to density map (taken from phytools blog)
n <- length(density_maps$vert_insect_ARD$cols)
density_maps$vert_insect_ARD$cols[1:n] <- colorRampPalette(c(my_colours$abiotic_vert_insect[2], my_colours$abiotic_vert_insect[3]), space="Lab")(n)
# then plot it and export to pdf
pdf(file = "figures/densitymap_vert_insect_ARD_1000_tall_younger.pdf", width = 12.6, height = 20, useDingbats = FALSE)

plot(density_maps$vert_insect_ARD, ftype = "off", lwd = 4, xlim = c(-15, 198))

# draw time axis 
# adapted from http://blog.phytools.org/2018/02/another-technique-for-including-time.html
T <- max(nodeHeights(tree_vert_insect))
obj <- axis(1, pos = -2, at = seq(T, -10, by = -40), cex.axis = 3,
            labels = FALSE, lwd = 2)
axis(side = 1, pos = -2, at = seq(T, -10, by = -40), cex.axis = 3,
     labels = FALSE, lwd = 2)
text(x = obj, y = rep(-26, length(obj)), T - obj, cex = 2)
text(x = mean(obj), y = -53, "time (mya)", cex = 2.1)
rm(T, obj)

# label clades
#ANA
segments(x0 = 158, y0 = 1, x1 = 158, y1 = 5, lwd = 10, col = "#bdbdbd", lend = "butt")
text(x = 158+5, y = 2.5, "ANA", srt = 0, cex = 2, col = "#bdbdbd", adj = 0)
#Magnoliids
segments(x0 = 158, y0 = 6, x1 = 158, y1 = 50, lwd = 10, col = "#636363", lend = "butt")
text(x = 158+5, y = 28, "Magnoliids", srt = 0, cex = 2, col = "#636363", adj = 0)
#Monocots
segments(x0 = 158, y0 = 52, x1 = 158, y1 = 230, lwd = 10, col = "#bdbdbd", lend = "butt")
text(x = 158+5, y = 141, "Monocots", srt = 0, cex = 2, col = "#bdbdbd", adj = 0)
#Commelinids
segments(x0 = 158+2, y0 = 171, x1 = 158+2, y1 = 230, lwd = 10, col = "#636363", lend = "butt")
text(x = 158+5, y = 200.5, "Commelinids", srt = 0, cex = 2, col = "#636363", adj = 0)
#Eudicots
segments(x0 = 158, y0 = 231, x1 = 158, y1 = 962, lwd = 10, col = "black", lend = "butt")
text(x = 158+5, y = 596.5, "Eudicots", srt = 0, cex = 2, col = "black", adj = 0)
#Rosids
segments(x0 = 158+2, y0 = 281, x1 = 158+2, y1 = 533, lwd = 10, col = "#bdbdbd", lend = "butt")
text(x = 158+5, y = 407, "Rosids", srt = 0, cex = 2, col = "#bdbdbd", adj = 0)
#Asterids
segments(x0 = 158+2, y0 = 635, x1 = 158+2, y1 = 962, lwd = 10, col = "#bdbdbd", lend = "butt")
text(x = 158+5, y = 798.5, "Asterids", srt = 0, cex = 2, col = "#bdbdbd", adj = 0)
#Chloranthales
segments(x0 = 158, y0 = 50.5, x1 = 158, y1 = 51.5, lwd = 10, col = "black", lend = "butt")
text(x = 158+5, y = 51, "Chloranthales", srt = 0, cex = 2, col = "black", adj = 0)

dev.off()
rm(n)

rm(tree_vert_insect, tree_wind_animal, density_maps, simmap_descriptions)

#### time of transitions ####

# plot time of transitions across simmaps

# function to extract transition times from a tree
source("scripts/functions/transition_times.R")

### WIND ANIMAL ###
# apply function across list of multiple simulations
wa_transitions <- data.frame()
for(i in 1:length(simmaps$wind_animal_ARD)){
  temp <- cbind(i, transition_times(simmaps$wind_animal_ARD[[i]]))
  wa_transitions <- rbind(wa_transitions, temp)
}
rm(temp, i)

table(wa_transitions$transition)

# build new data frame with cumulative number of transitions
wa_trans_cumul <- data.frame()
for(n in 1:1000){
  trans <- wa_transitions %>%
    dplyr::filter(i == n) %>%
    dplyr::group_by(transition) %>%
    dplyr::mutate(trans_no = row_number(i))
  
  wa_trans_cumul <- rbind(wa_trans_cumul, trans)
}
rm(n, trans)

# replace 2->4 with wind to animal
wa_trans_cumul$transition <- gsub("2->4", "wind to animal", wa_trans_cumul$transition)
wa_trans_cumul$transition <- gsub("4->2", "animal to wind", wa_trans_cumul$transition)

# TAKE AVERAGE OF ALL TRANSITION TIMES so this can be graphed

# first need to know average number of transitions, rounded
trans_avg_length <- wa_trans_cumul %>%
  dplyr::group_by(transition, i) %>%
  dplyr::mutate(no_trans = max(trans_no)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(trans_no == no_trans) %>%
  dplyr::group_by(transition) %>%
  dplyr::mutate(avg_length = round(mean(no_trans), digits = 0)) %>%
  dplyr::select(transition, avg_length) %>%
  dplyr::ungroup() %>%
  dplyr::distinct()

# now knowing this, can rearrange data and average times across rows
avg_trans_times <- wa_trans_cumul %>%
  dplyr::group_by(transition, trans_no) %>%
  dplyr::mutate(avg_time = mean(time)) %>%
  dplyr::mutate(SE_time = sqrt(var(time) / length(time))) %>%
  dplyr::select(transition, trans_no, avg_time, SE_time) %>%
  dplyr::distinct()

# reduce avg_trans_times to trans_avg_length
avg_trans_times_w2a <- avg_trans_times %>%
  dplyr::filter(transition == "wind to animal") %>%
  dplyr::filter(trans_no <= trans_avg_length[1,2])
avg_trans_times_a2w <- avg_trans_times %>%
  dplyr::filter(transition == "animal to wind") %>%
  dplyr::filter(trans_no <= trans_avg_length[2,2])
avg_trans_times <- rbind(avg_trans_times_a2w, avg_trans_times_w2a)
rm(avg_trans_times_a2w, avg_trans_times_w2a)

# export these results to csv in case I need them
readr::write_csv(avg_trans_times, "results/mean_transition_times_wind_animal_younger.csv")

## figures ##

# histogram of wind and animal pollination transition timing
pdf("figures/wind_animal_transition_times_mean_hist_younger.pdf", height = 2.8, width = 6)
min <- min(avg_trans_times$avg_time)
max <- max(avg_trans_times$avg_time)
ax <- pretty(min:max, n = 20)
wind_to_animal <- avg_trans_times %>%
  dplyr::filter(transition == "wind to animal") %>%
  dplyr::mutate(wind_to_animal = avg_time) %>%
  dplyr::ungroup() %>%
  dplyr::select(wind_to_animal)
animal_to_wind <- avg_trans_times %>%
  dplyr::filter(transition == "animal to wind") %>%
  dplyr::mutate(animal_to_wind = avg_time) %>%
  dplyr::ungroup() %>%
  dplyr::select(animal_to_wind)
wb <- hist(wind_to_animal$wind_to_animal, breaks = ax, plot = FALSE)
bw <- hist(animal_to_wind$animal_to_wind, breaks = ax, plot = FALSE)
plot (bw, col = my_colours$wind_water_animal[1], xlab = "Time of transitions (mya)",  
      main = "", ylab = "number of transitions", 
      ylim = c(0, 8), xlim = c(160,0)) # alter if x values change!
plot (wb, col = alpha(my_colours$wind_water_animal[3], 0.9), add = TRUE)
dev.off()
rm(min, max, ax, wb, bw, wind_to_animal, animal_to_wind)

# density plot of transitions from wind to animal pollination and back
pdf("figures/density_windanimal_younger.pdf", height = 3, width = 6)
# data for density graph
wind_to_animal <- wa_transitions %>%
  dplyr::filter(transition == "2->4") %>%
  dplyr::mutate(wind_to_animal = time) %>%
  dplyr::select(wind_to_animal, simulation = i)
animal_to_wind <- wa_transitions %>%
  dplyr::filter(transition == "4->2") %>%
  dplyr::mutate(animal_to_wind = time) %>%
  dplyr::select(animal_to_wind, simulation = i)
# calculate density curve
density_wb <- density(wind_to_animal$wind_to_animal)
density_bw <- density(animal_to_wind$animal_to_wind)
# plot the density
plot(density_wb, lwd = 2, col = my_colours$wind_water_animal[3], 
     xlim = c(160,0), xlab = "Time of transitions (mya)", bty = "l",
     cex.lab = 1.4, cex.axis = 1.4, main = NULL, sub = NULL, title = NULL)
lines(density_bw, lwd = 2, col = my_colours$wind_water_animal[1], xlim = c(160,0))
# add data-points with noise in the X-axis
rug(jitter(animal_to_wind$animal_to_wind), col = alpha(my_colours$wind_water_animal[1], 0.5))
rug(jitter(wind_to_animal$wind_to_animal), col = alpha(my_colours$wind_water_animal[3], 0.5))
dev.off()
rm(density_wb, density_bw, wind_to_animal, animal_to_wind)

# graph times all together (lines and their average)
pdf("figures/wind_animal_transitions_plus_mean_younger.pdf", height = 8, width = 10)

par(mar = c(5.1, 5.1, 4.1, 2.1))    # increase margins

animal_to_wind <- wa_trans_cumul %>%
  dplyr::filter(i == 1, transition == "4->2")

# first set up basic plot parameters
plot(animal_to_wind$trans_no ~ animal_to_wind$time,
     type = "p", bty = "l", xlim = c(160,0), ylim = c(0,60),
     col = alpha(my_colours$wind_water_animal[1], 0.1), pch = 15,
     xlab = "Time of transitions (mya)", 
     ylab = "Cumulative number of transitions",
     cex.lab = 1.8, cex.axis = 1.8)

# then loop through all data and add all points and lines to plot for 1000 simulations
for(n in 1:1000){
  test <- wa_transitions %>%
    dplyr::filter(i == n) %>%
    dplyr::group_by(transition) %>%
    dplyr::mutate(trans_no = row_number(i))
  
  wind_to_animal <- test %>%
    dplyr::filter(transition == "2->4")
  animal_to_wind <- test %>%
    dplyr::filter(transition == "4->2")
  
  points(animal_to_wind$trans_no ~ animal_to_wind$time, 
         col = alpha(my_colours$wind_water_animal[1], 0.5), pch = 15)
  lines(animal_to_wind$trans_no ~ animal_to_wind$time, 
        col = alpha(my_colours$wind_water_animal[1], 0.5))
  points(wind_to_animal$trans_no ~ wind_to_animal$time, 
         col = alpha(my_colours$wind_water_animal[3], 0.5), pch = 17)
  lines(wind_to_animal$trans_no ~ wind_to_animal$time, 
        col = alpha(my_colours$wind_water_animal[3], 0.5))
}

# then add average points and line in blue
# first prep data
wind_to_animal <- avg_trans_times %>%
  dplyr::filter(transition == "wind to animal") %>%
  dplyr::filter(trans_no <= trans_avg_length[1,2])
animal_to_wind <- avg_trans_times %>%
  dplyr::filter(transition == "animal to wind") %>%
  dplyr::filter(trans_no <= trans_avg_length[2,2])
# then add to plot
points(animal_to_wind$trans_no ~ animal_to_wind$avg_time,
       col = "blue", pch = 15)
lines(animal_to_wind$trans_no ~ animal_to_wind$avg_time, 
      col = "blue")
points(wind_to_animal$trans_no ~ wind_to_animal$avg_time,
       col = "blue", pch = 17)
lines(wind_to_animal$trans_no ~ wind_to_animal$avg_time, 
      col = "blue")
# add legend
legend("topleft",
       legend = c("animal to wind", "wind to animal"),
       col = c(my_colours$wind_water_animal[1], my_colours$wind_water_animal[3]),
       pch = c(15, 17), pt.lwd = 0.001, bty = "n", cex = 1.8)
dev.off()
rm(n, test, animal_to_wind, wind_to_animal)

rm(wa_trans_cumul, trans_avg_length, avg_trans_times, wa_transitions)

### VERT INSECT ###
# apply function across list of multiple simulations
vi_transitions <- data.frame()
for(i in 1:length(simmaps$vert_insect_ARD)){
  temp <- cbind(i, transition_times(simmaps$vert_insect_ARD[[i]]))
  vi_transitions <- rbind(vi_transitions, temp)
}
rm(temp, i)

table(vi_transitions$transition)

# build new data frame with cumulative number of transitions
vi_trans_cumul <- data.frame()
for(n in 1:1000){
  trans <- vi_transitions %>%
    dplyr::filter(i == n) %>%
    dplyr::group_by(transition) %>%
    dplyr::mutate(trans_no = row_number(i))
  
  vi_trans_cumul <- rbind(vi_trans_cumul, trans)
}
rm(n, trans)

vi_trans_cumul$transition <- gsub("5->6", "insect to vertebrate", vi_trans_cumul$transition)
vi_trans_cumul$transition <- gsub("6->5", "vertebrate to insect", vi_trans_cumul$transition)

# TAKE AVERAGE OF ALL TRANSITION TIMES so this can be graphed

# first need to know average number of transitions, rounded
trans_avg_length <- vi_trans_cumul %>%
  dplyr::group_by(transition, i) %>%
  dplyr::mutate(no_trans = max(trans_no)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(trans_no == no_trans) %>%
  dplyr::group_by(transition) %>%
  dplyr::mutate(avg_length = round(mean(no_trans), digits = 0)) %>%
  dplyr::select(transition, avg_length) %>%
  dplyr::ungroup() %>%
  dplyr::distinct()

# now knowing this, can rearrange data and average times across rows
avg_trans_times <- vi_trans_cumul %>%
  dplyr::group_by(transition, trans_no) %>%
  dplyr::mutate(avg_time = mean(time)) %>%
  dplyr::mutate(SE_time = sqrt(var(time) / length(time))) %>%
  dplyr::select(transition, trans_no, avg_time, SE_time) %>%
  dplyr::distinct()

# reduce avg_trans_times to trans_avg_length
avg_trans_times_i2v <- avg_trans_times %>%
  dplyr::filter(transition == "insect to vertebrate") %>%
  dplyr::filter(trans_no <= trans_avg_length[1,2])
avg_trans_times_v2i <- avg_trans_times %>%
  dplyr::filter(transition == "vertebrate to insect") %>%
  dplyr::filter(trans_no <= trans_avg_length[2,2])
avg_trans_times <- rbind(avg_trans_times_i2v, avg_trans_times_v2i)
rm(avg_trans_times_i2v, avg_trans_times_v2i)

# export these results to csv in case I need them
readr::write_csv(avg_trans_times, "results/mean_transition_times_vert_insect_younger.csv")

## time figures ##

# histogram of vertebrate and insect pollination transition timing average
pdf("figures/vert_insect_transition_times_mean_hist_younger.pdf", height = 2.8, width = 6)
min <- min(avg_trans_times$avg_time)
max <- max(avg_trans_times$avg_time)
ax <- pretty(min:max, n = 20)
insect_to_vert <- avg_trans_times %>%
  dplyr::filter(transition == "insect to vertebrate") %>%
  dplyr::mutate(insect_to_vert = avg_time) %>%
  dplyr::ungroup() %>%
  dplyr::select(insect_to_vert)
vert_to_insect <- avg_trans_times %>%
  dplyr::filter(transition == "vertebrate to insect") %>%
  dplyr::mutate(vert_to_insect = avg_time) %>%
  dplyr::ungroup() %>%
  dplyr::select(vert_to_insect)
wb <- hist(insect_to_vert$insect_to_vert, breaks = ax, plot = FALSE)
bw <- hist(vert_to_insect$vert_to_insect, breaks = ax, plot = FALSE)
plot (bw, col = my_colours$abiotic_vert_insect[2], xlab = "Time of transitions (mya)",  
      main = "", ylab = "number of transitions", 
      ylim = c(0, 8), xlim = c(160,0)) # alter if values change!
plot (wb, col = alpha(my_colours$abiotic_vert_insect[3], 0.7), add = TRUE)
dev.off()
rm(min, max, ax, wb, bw, insect_to_vert, vert_to_insect)

# density plot of transitions from insect to vert pollination and back
pdf("figures/density_vertinsect_younger.pdf", height = 3, width = 6)
# prep data for figure
insect_to_vert <- vi_transitions %>%
  dplyr::filter(transition == "5->6") %>%
  dplyr::mutate(insect_to_vert = time) %>%
  dplyr::select(insect_to_vert, simulation = i)
vert_to_insect <- vi_transitions %>%
  dplyr::filter(transition == "6->5") %>%
  dplyr::mutate(vert_to_insect = time) %>%
  dplyr::select(vert_to_insect, simulation = i)
# calculate density curve
density_vi <- density(vert_to_insect$vert_to_insect)
density_iv <- density(insect_to_vert$insect_to_vert)
# plot the density
plot(density_vi, lwd = 2, col = my_colours$abiotic_vert_insect[2], 
     xlim = c(160,0), xlab = "Time of transitions (mya)", bty = "l",
     cex.lab = 1.4, cex.axis = 1.4, main = NULL, sub = NULL, title = NULL)
lines(density_iv, lwd = 2, col = my_colours$abiotic_vert_insect[3], xlim = c(160,0))
# add data-points with noise in the X-axis
rug(jitter(insect_to_vert$insect_to_vert), col = alpha(my_colours$abiotic_vert_insect[3], 0.5))
rug(jitter(vert_to_insect$vert_to_insect), col = alpha(my_colours$abiotic_vert_insect[2], 0.5))
dev.off()
rm(density_vi, density_iv, insect_to_vert, vert_to_insect)

# graph all together (lines and their average)
pdf("figures/vert_insect_transitions_plus_mean_younger.pdf", height = 8, width = 10)

par(mar = c(5.1, 5.1, 4.1, 2.1))    # increase margins

vert_to_insect <- vi_trans_cumul %>%
  dplyr::filter(i == 1, transition == "6->5")

# first set up basic plot parameters
plot(vert_to_insect$trans_no ~ vert_to_insect$time,
     type = "p", bty = "l", xlim = c(160,0), ylim = c(0,75),
     col = alpha(my_colours$abiotic_vert_insect[2], 0.1), pch = 17,
     xlab = "Time of transitions (mya)", 
     ylab = "Cumulative number of transitions",
     cex.lab = 1.8, cex.axis = 1.8)

# then loop through all data and add all points and lines to plot for 1000 simulations
for(n in 1:1000){
  test <- vi_transitions %>%
    dplyr::filter(i == n) %>%
    dplyr::group_by(transition) %>%
    dplyr::mutate(trans_no = row_number(i))
  
  insect_to_vert <- test %>%
    dplyr::filter(transition == "5->6")
  vert_to_insect <- test %>%
    dplyr::filter(transition == "6->5")
  
  points(vert_to_insect$trans_no ~ vert_to_insect$time, 
         col = alpha(my_colours$abiotic_vert_insect[2], 0.5), pch = 17)
  lines(vert_to_insect$trans_no ~ vert_to_insect$time, 
        col = alpha(my_colours$abiotic_vert_insect[2], 0.5))
  points(insect_to_vert$trans_no ~ insect_to_vert$time, 
         col = alpha(my_colours$abiotic_vert_insect[3], 0.5), pch = 15)
  lines(insect_to_vert$trans_no ~ insect_to_vert$time, 
        col = alpha(my_colours$abiotic_vert_insect[3], 0.5))
}

# then add average points and line in blue
# first prep data
insect_to_vert <- avg_trans_times %>%
  dplyr::filter(transition == "insect to vertebrate") %>%
  dplyr::filter(trans_no <= trans_avg_length[1,2])
vert_to_insect <- avg_trans_times %>%
  dplyr::filter(transition == "vertebrate to insect") %>%
  dplyr::filter(trans_no <= trans_avg_length[2,2])
# then add to plot
points(vert_to_insect$trans_no ~ vert_to_insect$avg_time,
       col = "blue", pch = 17)
lines(vert_to_insect$trans_no ~ vert_to_insect$avg_time, 
      col = "blue")
points(insect_to_vert$trans_no ~ insect_to_vert$avg_time,
       col = "blue", pch = 15)
lines(insect_to_vert$trans_no ~ insect_to_vert$avg_time, 
      col = "blue")
# add legend
legend("topleft",
       legend = c("insect to vertebrate", "vertebrate to insect"),
       col = c(my_colours$abiotic_vert_insect[3], my_colours$abiotic_vert_insect[2]),
       pch = c(15, 17), pt.lwd = 0.001, bty = "n", cex = 1.8)
dev.off()
rm(n, test, vert_to_insect, insect_to_vert)

rm(transition_times, vi_trans_cumul, trans_avg_length, avg_trans_times, vi_transitions)

rm(simmaps, my_colours, no_cores)
