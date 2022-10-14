# clade labelling as per Ramirez-Barahona at al (2020) for a straight tree

RB2020_cladelabels <- function(xpos = 245){
  # clade labelling as per Ramirez-Barahona at al (2020)
  segments(x0 = xpos, y0 = 1, x1 = xpos, y1 = 12, lwd = 3, col = "#BBCDE9")
  text(x = xpos+3, y = 6.5, "ANA", srt = 270, cex = 0.5, col = "#BBCDE9")
  segments(x0 = xpos, y0 = 13, x1 = xpos, y1 = 66, lwd = 3, col = "#47A1D1")
  text(x = xpos+3, y = 39.5, "Magnoliids", srt = 270, cex = 0.5, col = "#47A1D1")
  segments(x0 = xpos-10, y0 = 67, x1 = xpos-10, y1 = 68, lwd = 3)
  text(x = xpos-7, y = 67.5, "Chloranthales", adj = 0, srt = 0, cex = 0.5)
  segments(x0 = xpos, y0 = 69, x1 = xpos, y1 = 314, lwd = 3, col = "#59BE1C")
  text(x = xpos+3, y = 191.5, "Monocots", srt = 270, cex = 0.5, col = "#59BE1C")
  segments(x0 = xpos-10, y0 = 209, x1 = xpos-10, y1 = 314, lwd = 3, col = "#0C9934")
  text(x = xpos-7, y = 261.5, "Commelinids", srt = 270, cex = 0.5, col = "#0C9934")
  segments(x0 = xpos-10, y0 = 315, x1 = xpos-10, y1 = 316, lwd = 3, col = "#E95208")
  text(x = xpos-7, y = 315.5, "Ceratophyllales", adj = 0, srt = 0, cex = 0.5, col = "#E95208")
  segments(x0 = xpos, y0 = 317, x1 = xpos, y1 = 1201, lwd = 3, col = "#F0D01B")
  text(x = xpos+3, y = 759, "Eudicots", srt = 270, cex = 0.5, col = "#F0D01B")
  segments(x0 = xpos-10, y0 = 390, x1 = xpos-10, y1 = 712, lwd = 3, col = "#F2B211")
  text(x = xpos-7, y = 551, "Rosids", srt = 270, cex = 0.5, col = "#F2B211")
  segments(x0 = xpos-10, y0 = 845, x1 = xpos-10, y1 = 1201, lwd = 3, col = "#FCB98E")
  text(x = xpos-7, y = 1023, "Asterids", srt = 270, cex = 0.5, col = "#FCB98E")
}