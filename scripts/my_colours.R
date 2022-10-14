# assign fixed colour scales to pollination modes
my_colours <- list()

# try to assign fixed colour to each possible state
my_colours$all <- c("#92BFB1", "#F6BA65", "#306BAC", "#BD8DAA", "#A891BC", 
                    "#FB8072", "#F1BFA0", "#F1BFA0", "#BFC1C6", "#F1BFA0",
                    "#C5A777", "#F1BFA0", "#C5A777", "#BFC1C6")
names(my_colours$all) <- c("abiotic", "wind", "water", "animal", "insect", 
                           "vertebrate", "abiotic&animal", "wind&animal",
                           "water&animal", "abiotic&insect",
                           "abiotic&vertebrate", "wind&insect",
                           "wind&vertebrate", "water&insect")
my_colours$abiotic_animal <- c("#92BFB1","#BD8DAA")
names(my_colours$abiotic_animal) <- c("abiotic", "animal")
my_colours$abiotic_animal_nopoly <- c("#92BFB1","#BD8DAA", "#F1BFA0")
names(my_colours$abiotic_animal_nopoly) <- c("abiotic", "animal", "abiotic&animal")
my_colours$wind_water_animal <- c("#F6BA65", "#306BAC", "#BD8DAA")
names(my_colours$wind_water_animal) <- c("wind", "water", "animal")
my_colours$wind_water_animal_nopoly <- c("#92BFB1", "#F6BA65", "#306BAC",
                                         "#BD8DAA", "#F1BFA0", "#BFC1C6")
names(my_colours$wind_water_animal_nopoly) <- c("wind&water", "wind", "water", 
                                                "animal", "wind&animal", 
                                                "water&animal")
my_colours$abiotic_vert_insect <- c("#92BFB1", "#A891BC", "#FB8072")
names(my_colours$abiotic_vert_insect) <- c("abiotic", "insect", "vertebrate")
my_colours$abiotic_vert_insect_nopoly <- c("#92BFB1", "#A891BC", "#FB8072", "#F1BFA0",
                                           "#C5A777", "#BD8DAA")
names(my_colours$abiotic_vert_insect_nopoly) <- c("abiotic", "insect", "vertebrate",
                                                  "abiotic&insect", "abiotic&vertebrate",
                                                  "vertebrate&insect")
my_colours$wind_water_vert_insect <- c("#F6BA65", "#306BAC", "#A891BC", "#FB8072")
names(my_colours$wind_water_vert_insect) <- c("wind", "water", "insect", 
                                              "vertebrate")
my_colours$wind_water_vert_insect_nopoly <- c("#F6BA65", "#306BAC", "#A891BC",
                                              "#FB8072", "#F1BFA0", "#92BFB1",
                                              "#C5A777", "#BFC1C6", "#BD8DAA")
names(my_colours$wind_water_vert_insect_nopoly) <- c("wind", "water", "insect", 
                                                     "vertebrate",
                                                     "wind&insect",
                                                     "wind&water",
                                                     "wind&vertebrate",
                                                     "water&insect",
                                                     "insect&vertebrate")

# colour scales for transitions between pollination modes
#assign fixed colour scale to biomes

temp <- c("#F6BA65", "#BD8DAA")
names(temp) <- c("animal to wind", "wind to animal")
my_colours$wind_animal <- scale_colour_manual(name = "transition", values = temp)
rm(temp)

# to add - colours for biome and pollination mode?




