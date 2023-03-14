# script to visualise distribution of pollination system vs syndrome data across
# angiosperm phylogeny

# first subset data and prep
pollination_sampling <- pollination1209 %>%
  dplyr::select(taxon_name:ParCladeTax, conf_score,
                wind_pollination_explicitly_tested:syndrome_or_system, 
                wind_water_vert_insect, position)

# now map system/syndrome onto phylogeny and have a look

### CIRCULAR PLOT ###
pdf(file = "figures/pollination_sampling.pdf", width = 30, height = 30, useDingbats = FALSE)

# calculate TIPS of important clades to label - try orders for now
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

# filter out smaller orders for clade labelling, leaves 26 of 63 orders 
for_cladelabels <- for_cladelabels %>%
  dplyr::filter(tiprange > 15|ParOrdTax == "Fagales") #filter
# and replace hard-to-label clades with abbreviations
for_cladelabels$ParOrdTax <- gsub("Ranunculales", "Ranun.", for_cladelabels$ParOrdTax)
for_cladelabels$ParOrdTax <- gsub("Fagales", "Faga.", for_cladelabels$ParOrdTax)
for_cladelabels$ParOrdTax <- gsub("Cucurbitales", "Cucu.", for_cladelabels$ParOrdTax)
for_cladelabels$ParOrdTax <- gsub("Solanales", "Solan.", for_cladelabels$ParOrdTax)
for_cladelabels$ParOrdTax <- gsub("Boraginales", "Borag.", for_cladelabels$ParOrdTax)

# plotting tip labels needs traits as data frame with species row names
# for pollination systems vs syndromes
system_syndrome <- pollination_sampling %>%
  dplyr::select(taxon_name, value = syndrome_or_system)
system_syndrome$value <- gsub("system", "1", system_syndrome$value)
system_syndrome$value <- gsub("syndrome", "0", system_syndrome$value)
system_syndrome <- as.data.frame(system_syndrome)
rownames(system_syndrome) <- system_syndrome[,1]
system_syndrome[,1] <- NULL

# and plot it!
cols <- c("#b2abd2", "#e66101")
names(cols) <- c(0, 1)
ape::plot.phylo(tree, type = "fan", show.tip.label = FALSE,
                     x.lim = c(-220, 220), y.lim = c(-220, 240), lwd = 1)
tiplabels(pch = 15, col = cols[as.numeric(system_syndrome[tree$tip.label,])+1])

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
legend(x = -220, y = 220, legend = c("syndrome", "system"), bg = "white",
       fill = c("#b2abd2", "#e66101"), cex = 3, pt.lwd = 0.001, bty = "n",
       title = "Pollination data source")

dev.off()
rm(cols, RB2020_cladelabels_fan, arclabel, system_syndrome)
