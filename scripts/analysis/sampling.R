# script to assess distribution of pollination system vs syndrome data across
# angiosperm phylogeny, and quality (confidence score) of said data

# first subset data and prep

pollination_sampling <- pollination1209 %>%
  dplyr::select(taxon_name:ParCladeTax, conf_score,
                wind_pollination_explicitly_tested:syndrome_or_system, 
                wind_water_vert_insect, position)

# now what? just map system/syndrome onto phylogeny and have a look?

# plotting tip labels needs traits as data frame with species row names
# for pollination systems vs syndromes
system_syndrome <- pollination_sampling %>%
  dplyr::select(taxon_name, value = syndrome_or_system)
system_syndrome$value <- gsub("system", "1", system_syndrome$value)
system_syndrome$value <- gsub("syndrome", "0", system_syndrome$value)
system_syndrome <- as.data.frame(system_syndrome)
rownames(system_syndrome) <- system_syndrome[,1]
system_syndrome[,1] <- NULL
# and for pollination data confidence scores
conf_score <- pollination_sampling %>%
  dplyr::select(taxon_name, value = conf_score)
conf_score <- as.data.frame(conf_score)
rownames(conf_score) <- conf_score[,1]
conf_score[,1] <- NULL

co <- c("#fc8d62", "#66c2a5")
col2 <- c("#ffeda0", "#feb24c", "#f03b20", "#e0f3db", "#a8ddb5", "#43a2ca")
plot(tree, cex = 0.5)
tiplabels(pch = 19, col = co[as.numeric(system_syndrome[tree$tip.label,])+1], 
          cex = 0.5, ftype = "off")
tiplabels(pch = 15, col = col2[as.numeric(conf_score[tree$tip.label,])+1], 
          cex = 0.5, ftype = "off", offset = 5)

# looks okay, now just need to do this but with circular phylogeny and labels?
# can I also add tiplabels for 

### CIRCULAR PLOT ###
pdf(file = "figures/pollination_sampling.pdf", 
    width = 50, height = 30, useDingbats = FALSE)

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

# and plot it!
cols1 <- c("#fc8d62", "#66c2a5")
names(cols1) <- c(0, 1)
cols2 <- c("#ffeda0", "#feb24c", "#f03b20", "#e0f3db", "#a8ddb5", "#43a2ca")
names(cols2) <- c(1, 2, 3, 4, 5, 6)
ape::plot.phylo(tree, type = "fan", show.tip.label = FALSE,
                     x.lim = c(-220, 220), y.lim = c(-220, 300), lwd = 1#,
                     #tips = seq(5, 1195, by = 1190/1200), maxY = 1201 # make space between two ends of phylogeny for time axis
                     #part = 0.99 #try plotting in just 99% of space so room for time axis
)
tiplabels(pch = 19, col = cols1[as.numeric(system_syndrome[tree$tip.label,])+1])
tiplabels(pch = 15, col = cols2[as.numeric(conf_score[tree$tip.label,])+1], 
         offset = 2)

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
       fill = c("#fc8d62", "#66c2a5"), cex = 3, pt.lwd = 0.001, bty = "n",
       title = "Pollination data source")
legend(x = 220, y = 220, 
       legend = c(stringr::str_wrap("6 - Pollination system, wind pollination explicitly tested for, pollination rewards clear, pollination efficiency of different pollinators tested"),
                  stringr::str_wrap("5 - Pollination system but wind pollination not explicitly tested for, pollination rewards clear, pollination efficiency of different pollinators not explicitly tested but pollinators observed transporting pollen and effecting pollination"),
                  stringr::str_wrap("4 - Pollination system based on general observations of pollinator visitation and informed by syndrome"),
                  stringr::str_wrap("3 - Pollination syndrome but specific adaptations for insect, vertebrate, wind or water pollination clearly described"),
                  stringr::str_wrap("2 - Pollination syndrome based on interpretations of images or illustrations and text, and pollination in similar closely related species"),
                  stringr::str_wrap("1 - Pollination syndrome based on interpretations of images or illustrations and text but least confident")), bg = "white",
       fill = c("#43a2ca", "#a8ddb5", "#e0f3db", "#f03b20", "#feb24c", "#ffeda0"), cex = 2, pt.lwd = 0.001, bty = "n",
       title = "Pollination data confidence score")
dev.off()
rm(cols1, cols2, RB2020_cladelabels_fan, arclabel)
