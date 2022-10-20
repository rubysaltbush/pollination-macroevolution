#####Stochastic character mapping#####

#### prepare data ####

# drop tips of tree where binary characters missing
tree_wind_animal <- ape::drop.tip(tree, pollination1209$taxon_name[pollination1209$wind_animal == "?"])
tree_vert_insect <- ape::drop.tip(tree, pollination1209$taxon_name[pollination1209$vert_insect == "?"])

# remove ER models from list
ASR_forsimmap <- ASR_forsimmap[-c(2,4)]

# and add full 4 state model to list
ASR_forsimmap$wind_water_vert_insect_ARD <- ASR_forsimmap$wind_water_vert_insect_ARD

# remove other ASR, not needed now
rm(ASR)

#### run simmap analyses ####

# for loop over ASRs with binary states, run simmapping for each model using corHMMM

# list to store results of simmaps
simmaps <- list()

# run simmaps in parallel
start_time <- Sys.time()
simmaps <- foreach(name = names(ASR_forsimmap)) %dopar% {
  phy <- ASR_forsimmap[[name]]$phy
  data <- ASR_forsimmap[[name]]$data
  data[is.na(data)] <- "?"
  model <- ASR_forsimmap[[name]]$solution
  model[is.na(model)] <- 0
  diag(model) <- -rowSums(model)
  result <- corHMM::makeSimmap(
    tree = phy,
    data = data,
    model = model,
    rate.cat = 1,
    nSim = 1000,
    nCores = floor(no_cores/length(ASR_forsimmap))
  )
  class(result) <- c("multiSimmap", "multiPhylo")
  result
}
names(simmaps) <- names(ASR_forsimmap)
end_time <- Sys.time()
end_time - start_time
# Time difference of 3 mins for 3*1000 sims on 3 cores
rm(end_time, start_time)

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
# Time difference of 2 mins
rm(name, start_time, end_time)

# check out results!
simmap_descriptions$wind_animal_ARD
simmap_density$wind_animal_ARD
simmap_descriptions$vert_insect_ARD
simmap_density$vert_insect_ARD
simmap_descriptions$wind_water_vert_insect_ARD
simmap_density$wind_water_vert_insect_ARD

# more transitions on average between insect and vertebrate pollination (73)
# than between wind and animal (53)
# 130 changes between pollination modes overall

# export print results from simmapping to text file
sink(file = "results/simmap_descriptions.txt")
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

### CIRCULAR PLOT ###
# PUBLICATION WORTHY FIGURE showing all pollination modes
# based on single random simmap but with ancestral nodes from ASR
pdf(file="figures/corHMM_ARD_wind_water_vert_insect_circular.pdf", 
    width = 30, height = 30, useDingbats = FALSE)

# first prep tip labels for polymorphic tips
# wind, water, vert and insect tips
wwvi_tips <- as.data.frame(ASR_forsimmap$wind_water_vert_insect_ARD$tip.states) %>%
  dplyr::mutate(statesum = dplyr::select(. , V1:V4) %>% rowSums) %>%
  dplyr::mutate(tip_num = 1:nrow(.)) %>% # add number to match tip number in tree
  dplyr::filter(statesum < 4) %>% # filter out taxa with missing tip data
  dplyr::mutate(across(c(1:4), .fns = ~./statesum)) %>% # divide states by sum of states
  tibble::remove_rownames() %>%
  tibble::column_to_rownames(var = "tip_num") %>% # rename rows with tip numbers
  dplyr::select(V1, V2, V3, V4) %>%
  as.matrix()

# then calculate TIPS of important clades to label - try orders for now
for_cladelabels <- pollination1209 %>%
  dplyr::select(position, ParOrdTax) %>%
  dplyr::group_by(ParOrdTax) %>%
  tidyr::nest() %>%
  dplyr::ungroup() %>%
  dplyr::mutate(mintip = purrr::map(data, ~{min(.x$position)})) %>%
  dplyr::mutate(maxtip = purrr::map(data, ~{max(.x$position)})) %>%
  tidyr::unnest(cols = c(data, mintip, maxtip)) %>%
  dplyr::select(ParOrdTax, mintip, maxtip) %>%
  dplyr::distinct() %>%
  dplyr::mutate(tiprange = maxtip - mintip)

# for nodes, calculate MRCA of all these tips to plot nodes for all orders
for_nodes <- for_cladelabels %>%
  dplyr::rowwise() %>%
  dplyr::mutate(MRCAn = ape::getMRCA(tree, c(mintip, maxtip)))

# then filter out smaller orders for clade labelling, leaves 26 of 63 orders 
for_cladelabels <- for_cladelabels %>%
  dplyr::filter(tiprange > 15|ParOrdTax == "Fagales") #filter
# and replace hard-to-label clades with number codes
for_cladelabels$ParOrdTax <- gsub("Ranunculales", "Ranun.", for_cladelabels$ParOrdTax)
for_cladelabels$ParOrdTax <- gsub("Fagales", "Faga.", for_cladelabels$ParOrdTax)
for_cladelabels$ParOrdTax <- gsub("Cucurbitales", "Cucu.", for_cladelabels$ParOrdTax)
for_cladelabels$ParOrdTax <- gsub("Solanales", "Solan.", for_cladelabels$ParOrdTax)
for_cladelabels$ParOrdTax <- gsub("Boraginales", "Borag.", for_cladelabels$ParOrdTax)

# and plot it!
cols <- my_colours$wind_water_vert_insect
names(cols) <- c(2, 3, 5, 6)
phytools::plotSimmap(simmaps$wind_water_vert_insect_ARD[1], 
                     colors = cols, 
                     ftype = "off", 
                     type = "fan",
                     xlim = c(-220, 220), ylim = c(-220, 240), lwd = 4,
                     tips = seq(5, 1195, by = 1190/1200), maxY = 1201 # make space between two ends of phylogeny for time axis
                     #part = 0.99 #try plotting in just 99% of space so room for time axis
                     )


# first label Cretaceous-Palaeogene boundary at 66 mya
plotrix::draw.circle(0, 0, radius = max(nodeHeights(tree)) - 66, 
                     col = "#dadada", lty = 0)

# then label Jurassic-Cretaceous boundary at 145 mya
plotrix::draw.circle(0, 0, radius = max(nodeHeights(tree)) - 145, 
                     col = "#f8f6f7", lty = 0)

# then add text to boundary points
text(x = max(nodeHeights(tree)) - 66, y = -4, "66 mya", cex = 2, col = "#636363")
text(x = max(nodeHeights(tree)) - 145, y = -6, "145 mya", cex = 2, col = "#636363")

# repeat plot call to draw plot over top of time labels
phytools::plotSimmap(simmaps$wind_water_vert_insect_ARD[1], 
                     colors = cols, 
                     ftype = "off", 
                     type = "fan",
                     xlim = c(-220, 220), ylim = c(-220, 240), lwd = 4,
                     tips = seq(5, 1195, by = 1190/1200), maxY = 1201, # make space between two ends of phylogeny for time axis
                     add = TRUE
)

# plot nodes to of MRCA of all orders with pollination state pie from ASR
ape::nodelabels(node = for_nodes$MRCAn,
                pie = ASR_forsimmap$wind_water_vert_insect_ARD$states[for_nodes$MRCAn-1201,], 
                piecol = my_colours$wind_water_vert_insect, cex = 0.5)
# plot tip labels with pie charts for polymorphies
names(cols) <- c("V1", "V2", "V3", "V4")
ape::tiplabels(tip = as.numeric(rownames(wwvi_tips)), pie = wwvi_tips,
                piecol = cols, cex = 0.15)

# source custom arc labelling function
source("scripts/functions/arclabel.R")

# # loop through and draw labels on phylogeny for larger orders
for(i in 1:length(for_cladelabels$ParOrdTax)) {
 arclabel(text = for_cladelabels$ParOrdTax[i],
          tips = c(for_cladelabels$mintip[i], for_cladelabels$maxtip[i]),
          cex = 2)
}
rm(i, for_cladelabels)

# label large clades as per RB2020 labelling
source("scripts/functions/RB2020_cladelabels_fan.R")
RB2020_cladelabels_fan(offset = 1.015)

# insert legend
legend(x = -220, y = 220, legend = names(my_colours$wind_water_vert_insect), bg = "white",
       fill = my_colours$wind_water_vert_insect, cex = 3, pt.lwd = 0.001, bty = "n",
       title = "Pollination modes")
dev.off()
rm(cols, RB2020_cladelabels_fan, arclabel, for_nodes)

### TALL PLOT ###
# tall plot same as above but with more detail for Supplementary Information
pdf(file="figures/corHMM_ARD_wind_water_vert_insect_tall.pdf", 
    width = 8, height = 80, useDingbats = FALSE)

# set colours
cols <- my_colours$wind_water_vert_insect
names(cols) <- c(2, 3, 5, 6)
# and plot it!
phytools::plotSimmap(simmaps$wind_water_vert_insect_ARD[1], 
                     colors = cols, lwd = 3,
                     ftype = "off", xlim = c(0, 310))

# plot all nodes with pie charts from ASR
ape::nodelabels(pie = ASR_forsimmap$wind_water_vert_insect_ARD$states, 
                piecol = my_colours$wind_water_vert_insect, cex = 0.5)

# plot tip labels with pie charts for polymorphies
names(cols) <- c("V1", "V2", "V3", "V4")
ape::tiplabels(tip = as.numeric(rownames(wwvi_tips)), pie = wwvi_tips,
               piecol = cols, cex = 0.25)

# label clades, orders and families
source("scripts/functions/RB2020_cladelabels.R")
source("scripts/functions/family_and_order_labels.R")
RB2020_cladelabels(xpos = 300)
order_labels(xpos = 250)
family_labels(xpos = 205)

# insert legend
legend(x = 0, y = 1100, legend = names(my_colours$wind_water_vert_insect), bg = "white",
       fill = my_colours$wind_water_vert_insect, cex = 1, pt.lwd = 0.001, bty = "n",
       title = "Pollination modes")

# manually draw time axis adapted from http://blog.phytools.org/2018/02/another-technique-for-including-time.html
# probably more complicated than it needs to be
T <- max(nodeHeights(tree))
tick.spacing <- 40
min.tick <- -10
obj <- axis(1, pos = 0, at = seq(T, min.tick, by = -tick.spacing), cex.axis = 0.5,
          labels = FALSE)
axis(side = 1, pos = 0, at = seq(T, min.tick, by = -tick.spacing), cex.axis = 0.5,
     labels = FALSE)
text(obj, rep(-5, length(obj)), T - obj, cex = 0.6)
text(mean(obj), -10 , "time (mya)", cex = 0.8)
rm(T, tick.spacing, min.tick, obj)
dev.off()
rm(wwvi_tips, RB2020_cladelabels, order_labels, family_labels, cols)


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
pdf(file = "figures/densitymap_wind_animal_ARD_1000_tall.pdf", width = 12.6, height = 20, useDingbats = FALSE)
plot(density_maps$wind_animal_ARD, ftype = "off", lwd = 4, xlim = c(-5, 245))

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
segments(x0 = 198, y0 = 1, x1 = 198, y1 = 10, lwd = 10, col = "#bdbdbd", lend = "butt")
text(x = 198+5, y = 5.5, "ANA", srt = 0, cex = 2, col = "#bdbdbd", adj = 0)
#Magnoliids
segments(x0 = 198, y0 = 11, x1 = 198, y1 = 63, lwd = 10, col = "#636363", lend = "butt")
text(x = 198+5, y = 37, "Magnoliids", srt = 0, cex = 2, col = "#636363", adj = 0)
#Monocots
segments(x0 = 198, y0 = 66, x1 = 198, y1 = 296, lwd = 10, col = "#bdbdbd", lend = "butt")
text(x = 198+5, y = 181, "Monocots", srt = 0, cex = 2, col = "#bdbdbd", adj = 0)
#Commelinids
segments(x0 = 198+2, y0 = 193, x1 = 198+2, y1 = 296, lwd = 10, col = "#636363", lend = "butt")
text(x = 198+5, y = 244.5, "Commelinids", srt = 0, cex = 2, col = "#636363", adj = 0)
#Eudicots
segments(x0 = 198, y0 = 297, x1 = 198, y1 = 1145, lwd = 10, col = "black", lend = "butt")
text(x = 198+5, y = 721, "Eudicots", srt = 0, cex = 2, col = "black", adj = 0)
#Rosids
segments(x0 = 198+2, y0 = 368, x1 = 198+2, y1 = 680, lwd = 10, col = "#bdbdbd", lend = "butt")
text(x = 198+5, y = 524, "Rosids", srt = 0, cex = 2, col = "#bdbdbd", adj = 0)
#Asterids
segments(x0 = 198+2, y0 = 806, x1 = 198+2, y1 = 1145, lwd = 10, col = "#bdbdbd", lend = "butt")
text(x = 198+5, y = 975.5, "Asterids", srt = 0, cex = 2, col = "#bdbdbd", adj = 0)
#Chloranthales
segments(x0 = 198, y0 = 64, x1 = 198, y1 = 65, lwd = 10, col = "black", lend = "butt")
text(x = 198+5, y = 64.5, "Chloranthales", srt = 0, cex = 2, col = "black", adj = 0)

dev.off()

# vertebrate and insect
# assign new colours to density map (taken from phytools blog)
n <- length(density_maps$vert_insect_ARD$cols)
density_maps$vert_insect_ARD$cols[1:n] <- colorRampPalette(c(my_colours$abiotic_vert_insect[2], my_colours$abiotic_vert_insect[3]), space="Lab")(n)
# then plot it and export to pdf
pdf(file = "figures/densitymap_vert_insect_ARD_1000_tall.pdf", width = 12.6, height = 20, useDingbats = FALSE)

plot(density_maps$vert_insect_ARD, ftype = "off", lwd = 4, xlim = c(-15, 233))

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
segments(x0 = 190, y0 = 1, x1 = 190, y1 = 5, lwd = 10, col = "#bdbdbd", lend = "butt")
text(x = 190+5, y = 2.5, "ANA", srt = 0, cex = 2, col = "#bdbdbd", adj = 0)
#Magnoliids
segments(x0 = 190, y0 = 6, x1 = 190, y1 = 50, lwd = 10, col = "#636363", lend = "butt")
text(x = 190+5, y = 28, "Magnoliids", srt = 0, cex = 2, col = "#636363", adj = 0)
#Monocots
segments(x0 = 190, y0 = 52, x1 = 190, y1 = 230, lwd = 10, col = "#bdbdbd", lend = "butt")
text(x = 190+5, y = 141, "Monocots", srt = 0, cex = 2, col = "#bdbdbd", adj = 0)
#Commelinids
segments(x0 = 190+2, y0 = 171, x1 = 190+2, y1 = 230, lwd = 10, col = "#636363", lend = "butt")
text(x = 190+5, y = 200.5, "Commelinids", srt = 0, cex = 2, col = "#636363", adj = 0)
#Eudicots
segments(x0 = 190, y0 = 231, x1 = 190, y1 = 962, lwd = 10, col = "black", lend = "butt")
text(x = 190+5, y = 596.5, "Eudicots", srt = 0, cex = 2, col = "black", adj = 0)
#Rosids
segments(x0 = 190+2, y0 = 281, x1 = 190+2, y1 = 533, lwd = 10, col = "#bdbdbd", lend = "butt")
text(x = 190+5, y = 407, "Rosids", srt = 0, cex = 2, col = "#bdbdbd", adj = 0)
#Asterids
segments(x0 = 190+2, y0 = 635, x1 = 190+2, y1 = 962, lwd = 10, col = "#bdbdbd", lend = "butt")
text(x = 190+5, y = 798.5, "Asterids", srt = 0, cex = 2, col = "#bdbdbd", adj = 0)
#Chloranthales
segments(x0 = 190, y0 = 50.5, x1 = 190, y1 = 51.5, lwd = 10, col = "black", lend = "butt")
text(x = 190+5, y = 51, "Chloranthales", srt = 0, cex = 2, col = "black", adj = 0)

dev.off()
rm(n)

rm(tree_vert_insect, tree_wind_animal, density_maps)

# Jakub's method to make histogram of # of transitions from simmapping

# wind and animal
pdf("figures/wind_animal_ARD_1000_transitions.pdf", height = 5, width = 5)
bwb <- simmap_descriptions$wind_animal_ARD$count[,2:3]
min <- min (bwb)
max <- max (bwb)
ax <- pretty(min:max, n = 20)
bw <- hist(bwb[,1], breaks = ax, plot = FALSE)
wb <- hist(bwb[,2], breaks = ax, plot = FALSE)
plot (bw, col = my_colours$wind_water_animal[3], 
      xlim = c(0,60), ylim = c(0,400), 
      xlab = "total number of transitions",  main = "", 
      ylab = "frequency across 1000 simulations")
plot (wb, col = my_colours$wind_water_animal[1], add = TRUE)
legend("topleft", legend = c("animal to wind", "wind to animal"),
       fill = c(my_colours$wind_water_animal[1], my_colours$wind_water_animal[3]),
       pt.lwd = 0.001, bty = "n")
dev.off()
rm(bwb, min, max, ax, wb, bw)

# vertebrate and insect
pdf("figures/vert_insect_ARD_1000_transitions.pdf", height = 5, width = 5)
ivi <- simmap_descriptions$vert_insect_ARD$count[,2:3]
min <- min (ivi)
max <- max (ivi)
ax <- pretty(min:max, n = 20)
vi <- hist(ivi[,1], breaks = ax, plot = FALSE)
iv <- hist(ivi[,2], breaks = ax, plot = FALSE)
plot (vi, col = my_colours$abiotic_vert_insect[3], 
      xlim = c(10,70), ylim = c(0, 400), 
      xlab = "total number of transitions",  main = "", 
      ylab = "frequency across 1000 simulations")
plot (iv, col = alpha(my_colours$abiotic_vert_insect[2], 0.7), add = TRUE)
legend("topleft", legend = c("insect to vertebrate", "vertebrate to insect"),
       fill = c(my_colours$abiotic_vert_insect[3], my_colours$abiotic_vert_insect[2]),
       pt.lwd = 0.001, bty = "n")
dev.off()
rm(ivi, min, max, ax, iv, vi)

# wind, water, vert, insect
# pooled so that displays transitions TO 4 states (not from)
pdf("figures/wwvi_ARD_1000_transitions.pdf", height = 5, width = 5)
wwvi <- data.frame(simmap_descriptions$wind_water_vert_insect_ARD$count[,2:13]) # get # of state changes
wwvi$to_wind <- wwvi[,4] + wwvi[,7] + wwvi[,10]
wwvi$to_water <- wwvi[,1] + wwvi[,8] + wwvi[,11]
wwvi$to_insect <- wwvi[,2] + wwvi[,5] + wwvi[,12]
wwvi$to_vert <- wwvi[,3] + wwvi[,6] + wwvi[,9]
wwvi <- as.matrix(wwvi[,13:16])
min <- min (wwvi)
max <- max (wwvi)
ax <- pretty(min:max, n = 20)
to_wind <- hist(wwvi[,1], breaks = ax, plot = FALSE)
to_water <- hist(wwvi[,2], breaks = ax, plot = FALSE)
to_insect <- hist(wwvi[,3], breaks = ax, plot = FALSE)
to_vert <- hist(wwvi[,4], breaks = ax, plot = FALSE)
plot (to_wind, col = my_colours$wind_water_vert_insect[1], 
      xlim = c(0,80), ylim = c(0,1000), 
      xlab = "total number of transitions",  main = "", 
      ylab = "frequency across 1000 simulations")
plot (to_water, col = my_colours$wind_water_vert_insect[2], 
      add = TRUE)
plot (to_vert, col = alpha(my_colours$wind_water_vert_insect[4], 0.7), 
      add = TRUE)
plot (to_insect, col = alpha(my_colours$wind_water_vert_insect[3], 0.7), 
      add = TRUE)
# insert legend
legend("topright", legend = c("to wind", "to water", "to insect", "to vertebrate"), 
       fill = my_colours$wind_water_vert_insect, pt.lwd = 0.001, bty = "n")
dev.off()
rm(wwvi, min, max, ax, to_wind, to_water, to_insect, to_vert, simmap_descriptions)

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
for(n in 1:500){
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

# export these results to csv in case I need them
readr::write_csv(avg_trans_times, "results/mean_transition_times_wind_animal.csv")

## figures ##

# histogram of wind and animal pollination transition timing
pdf("figures/wind_animal_transition_times_mean_hist.pdf", height = 3, width = 5)
min <- min(avg_trans_times$avg_time)
max <- max(avg_trans_times$avg_time)
ax <- pretty(min:max, n = 20)
wind_to_animal <- avg_trans_times %>%
  dplyr::ungroup() %>%
  dplyr::filter(transition == "wind to animal") %>%
  dplyr::mutate(wind_to_animal = avg_time) %>%
  dplyr::select(wind_to_animal)
animal_to_wind <- avg_trans_times %>%
  dplyr::ungroup() %>%
  dplyr::filter(transition == "animal to wind") %>%
  dplyr::mutate(animal_to_wind = avg_time) %>%
  dplyr::select(animal_to_wind)
wb <- hist(wind_to_animal$wind_to_animal, breaks = ax, plot = FALSE)
bw <- hist(animal_to_wind$animal_to_wind, breaks = ax, plot = FALSE)
plot (bw, col = my_colours$wind_water_animal[1], xlab = "time of transitions (mya)",  
      main = "", ylab = "number of transitions", 
      ylim = c(0, 6), xlim = c(200,0)) # alter if x values change!
plot (wb, col = alpha(my_colours$wind_water_animal[3], 0.9), add = TRUE)
dev.off()
rm(min, max, ax, wb, bw, wind_to_animal, animal_to_wind)

# density plot of transitions from wind to animal pollination and back
pdf("figures/density_windanimal.pdf", height = 3, width = 6)
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
     xlim = c(200,0), xlab = "Time of transitions (mya)", bty = "l",
     cex.lab = 1.4, cex.axis = 1.4, main = NULL, sub = NULL, title = NULL)
lines(density_bw, lwd = 2, col = my_colours$wind_water_animal[1], xlim = c(200,0))
# add data-points with noise in the X-axis
rug(jitter(animal_to_wind$animal_to_wind), col = alpha(my_colours$wind_water_animal[1], 0.5))
rug(jitter(wind_to_animal$wind_to_animal), col = alpha(my_colours$wind_water_animal[3], 0.5))
dev.off()
rm(density_wb, density_bw, wind_to_animal, animal_to_wind)

# graph times all together (lines and their average)
pdf("figures/wind_animal_transitions_plus_mean.pdf", height = 8, width = 10)

par(mar = c(5.1, 5.1, 4.1, 2.1))    # increase margins

animal_to_wind <- wa_trans_cumul %>%
  dplyr::filter(i == 1, transition == "4->2")

# first set up basic plot parameters
plot(animal_to_wind$trans_no ~ animal_to_wind$time,
     type = "p", bty = "l", xlim = c(200,0), ylim = c(0,60),
     col = alpha(my_colours$wind_water_animal[1], 0.1), pch = 15,
     xlab = "Time of transitions (mya)", 
     ylab = "Cumulative number of transitions",
     cex.lab = 1.8, cex.axis = 1.8)

# then loop through all data and add all points and lines to plot for 500 simulations
for(n in 1:500){
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
for(n in 1:500){
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

# export these results to csv in case I need them
readr::write_csv(avg_trans_times, "results/mean_transition_times_vert_insect.csv")

## time figures ##

# histogram of wind and animal pollination transition timing average
pdf("figures/vert_insect_transition_times_mean_hist.pdf", height = 3, width = 5)
min <- min(avg_trans_times$avg_time)
max <- max(avg_trans_times$avg_time)
ax <- pretty(min:max, n = 20)
insect_to_vert <- avg_trans_times %>%
  dplyr::ungroup() %>%
  dplyr::filter(transition == "insect to vertebrate") %>%
  dplyr::mutate(insect_to_vert = avg_time) %>%
  dplyr::select(insect_to_vert)
vert_to_insect <- avg_trans_times %>%
  dplyr::ungroup() %>%
  dplyr::filter(transition == "vertebrate to insect") %>%
  dplyr::mutate(vert_to_insect = avg_time) %>%
  dplyr::select(vert_to_insect)
wb <- hist(insect_to_vert$insect_to_vert, breaks = ax, plot = FALSE)
bw <- hist(vert_to_insect$vert_to_insect, breaks = ax, plot = FALSE)
plot (bw, col = my_colours$abiotic_vert_insect[2], xlab = "Time of transitions (mya)",  
      main = "", ylab = "number of transitions", 
      ylim = c(0, 12), xlim = c(200,0)) # alter if x values change!
plot (wb, col = alpha(my_colours$abiotic_vert_insect[3], 0.8), add = TRUE)
dev.off()
rm(min, max, ax, wb, bw, insect_to_vert, vert_to_insect)

# density plot of transitions from insect to vert pollination and back
pdf("figures/density_vertinsect.pdf", height = 3, width = 6)
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
     xlim = c(200,0), xlab = "Time of transitions (mya)", bty = "l",
     cex.lab = 1.4, cex.axis = 1.4, main = NULL, sub = NULL, title = NULL)
lines(density_iv, lwd = 2, col = my_colours$abiotic_vert_insect[3], xlim = c(200,0))
# add data-points with noise in the X-axis
rug(jitter(insect_to_vert$insect_to_vert), col = alpha(my_colours$abiotic_vert_insect[3], 0.5))
rug(jitter(vert_to_insect$vert_to_insect), col = alpha(my_colours$abiotic_vert_insect[2], 0.5))
dev.off()
rm(density_vi, density_iv, insect_to_vert, vert_to_insect)

# graph all together (lines and their average)
pdf("figures/vert_insect_transitions_plus_mean.pdf", height = 8, width = 10)

par(mar = c(5.1, 5.1, 4.1, 2.1))    # increase margins

vert_to_insect <- vi_trans_cumul %>%
  dplyr::filter(i == 1, transition == "6->5")

# first set up basic plot parameters
plot(vert_to_insect$trans_no ~ vert_to_insect$time,
     type = "p", bty = "l", xlim = c(200,0), ylim = c(0,65),
     col = alpha(my_colours$abiotic_vert_insect[2], 0.1), pch = 17,
     xlab = "Time of transitions (mya)", 
     ylab = "Cumulative number of transitions",
     cex.lab = 1.8, cex.axis = 1.8)

# then loop through all data and add all points and lines to plot for 500 simulations
for(n in 1:500){
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

rm(simmaps, ASR_forsimmap)
