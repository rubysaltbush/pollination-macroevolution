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

# for some reason my data do not work with the phytools::paintBranches function
# but I have fixed it in below version copied from phytools github
paintBranchesFixed <- function(tree, edge, state, anc.state = "1"){
  if(!inherits(tree, "phylo")) stop("tree should be an object of class \"phylo\".")
  if(is.null(tree$maps)) maps<-lapply(tree$edge.length,function(x) setNames(x,anc.state))
  else maps<-tree$maps
  ii<-sapply(edge,function(x,y) which(y==x),y=tree$edge[,2])
  ii<-unlist(ii) # fix added to prevent "Error in tree$edge.length[[ii[i]]] : invalid subscript type 'list'"
  for(i in 1:length(ii)) maps[[ii[i]]]<-setNames(tree$edge.length[[ii[i]]],state)
  ## build mapped.edge matrix
  s<-vector()
  for(i in 1:nrow(tree$edge)) s<-c(s,names(maps[[i]]))
  s<-unique(s)
  mapped.edge<-matrix(0,length(tree$edge.length),length(s),dimnames=list(edge=apply(tree$edge,1,function(x) paste(x,collapse=",")),state=s))
  for(i in 1:length(maps)) for(j in 1:length(maps[[i]])) mapped.edge[i,names(maps[[i]])[j]]<-mapped.edge[i,names(maps[[i]])[j]]+maps[[i]][j]
  ## add attributes to the tree
  tree$mapped.edge<-mapped.edge
  tree$maps<-maps
  class(tree)<-c("simmap",setdiff(class(tree),"simmap"))
  tree
}

# prep data for paintBranches plotting
sys_syn <- as.factor(system_syndrome$value)
names(sys_syn) <- rownames(system_syndrome)

syn <- names(sys_syn)[sys_syn == "0"]
sys <- names(sys_syn)[sys_syn == "1"]

syntree <- paintBranchesFixed(tree, edge = sapply(syn, match, tree$tip.label), 
                              state = "0", anc.state = "3")
syntree <- paintBranchesFixed(syntree, edge = sapply(sys, match, tree$tip.label), 
                              state = "1")

# prepare colours for branches and tips
cols <- c("#b2abd2", "#e66101", "black")
names(cols) <- c("0", "1", "3")
# and plot tree
phytools::plotSimmap(syntree, type = "fan", ftype = "off", colors = cols,
                     xlim = c(-220, 220), ylim = c(-220, 240), lwd = 2)
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
rm(cols, RB2020_cladelabels_fan, arclabel, system_syndrome, sys, syn, sys_syn,
   paintBranchesFixed, syntree, pollination_sampling)
