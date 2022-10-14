# greyscale clade labelling as per Ramirez-Barahona et al (2020) for circular (fan) tree
# using custom arclabel function

RB2020_cladelabels_fan <- function(offset = 1){ # use offset argument to move labels closer (<1) or further away (>1) from tree
  source("scripts/arclabel.R") # get arclabel function
  arclabel(text = "ANA", tips = c(1, 12), 
           lwd = 40, cex = 3.2, col = "#bdbdbd",
           ln.offset = offset + .07, lab.offset = offset + .11,
           lend = "butt")
  arclabel(text = "Magnoliids", tips = c(13, 66), 
           lwd = 40, cex = 3.2, col = "#636363",
           ln.offset = offset + .07, lab.offset = offset + .11,
           lend = "butt")
  arclabel(text = "Chlor.", tips = c(67, 68),
           lwd = 80, cex = 2, col = "black",
           ln.offset = offset + .05, lab.offset = offset + .095,
           orientation = "perpendicular",
           lend = "butt")
  arclabel(text = "Monocots", tips = c(69, 314), 
           lwd = 40, cex = 3.2, col = "#bdbdbd",
           ln.offset = offset + .07, lab.offset = offset + .11,
           lend = "butt")
  arclabel(text = "Commelinids", tips = c(209, 314), 
           lwd = 30, cex = 3.2, col = "#636363",
           ln.offset = offset + .06, lab.offset = offset + .11,
           lend = "butt")
  arclabel(text = "Cerat.", tips = c(315, 316),
           lwd = 80, cex = 2, col = "black",
           ln.offset = offset + .05, lab.offset = offset + .095,
           orientation = "perpendicular",
           lend = "butt")
  arclabel(text = "Eudicots", tips = c(317, 1201),
           lwd = 40, cex = 3.2, col = "#636363",
           ln.offset = offset + .07, lab.offset = offset + .11,
           lend = "butt")
  arclabel(text = "Rosids", tips = c(390, 712),
           lwd = 30, cex = 3.2, col = "#bdbdbd",
           ln.offset = offset + .06, lab.offset = offset + .11,
           lend = "butt")
  arclabel(text = "Asterids", tips = c(845, 1201), 
           lwd = 30, cex = 3.2, col = "#bdbdbd",
           ln.offset = offset + .06, lab.offset = offset + .11,
           lend = "butt")
}
