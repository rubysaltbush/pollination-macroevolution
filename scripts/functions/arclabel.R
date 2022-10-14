## function to create curved clade labels for a fan tree THAT ACTUALLY WORKS
## adapted from arc.cladelabels by Liam J. Revell 2017, 2022 by Gregory McIntyre

# TO USE within plotting call e.g. arclabel(text = "Clade1", tips = c(1, 12))
# orientation can be curved, perpendicular or horizontal

arclabel <- function(text, tips, ln.offset = 1.02, lab.offset = 1.035, 
                     orientation = "curved", lwd = 1, col = "black", cex = 1,
                     lend = "round") {
  # get plot from environment
  obj <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  # calculate distance from centre of circle to tip
  h <- max(sqrt(obj$xx^2 + obj$yy^2))
  # calculate angle (in degrees) from centre to point to label
  deg <- atan2(obj$yy[tips], obj$xx[tips]) * 180/pi
  # ensure deg2 is bigger than deg1, so when we calculate the midpoint using mean it will work fine
  if (deg[2] < deg[1]) {
    deg[2] <- deg[2] + 360
  }
  # get midpoint of arc to centre text here
  a <- mean(deg * pi/180)
  # draw line
  plotrix::draw.arc(x = 0, y = 0, radius = h * ln.offset, deg1 = deg[1], 
                    deg2 = deg[2], lwd = lwd, col = col, lend = lend)
  # draw text at different possible orientations
  if(orientation == "curved"){
    plotrix::arctext(text, radius = h * lab.offset, middle = a, 
                     cex = cex, col = col)
  }
  else if(orientation == "perpendicular"){
    plotrix::radialtext(text, start = h * lab.offset, angle = a, 
                        cex = cex, col = col)
  }
  else if(orientation == "horizontal"){
    x0 <- lab.offset * cos(median(deg) * pi/180) * h
    y0 <- lab.offset * sin(median(deg) * pi/180) * h
    text(x = x0, y = y0, label = text,
         adj = c(if(x0 >= 0) 0 else 1, if(y0 >= 0) 0 else 1),
         offset = 0, cex = cex, col = col)
  }
}

