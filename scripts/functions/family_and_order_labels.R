# Two functions to label straight trees

# function to work out highest and lowest tip numbers for each family,
# then loop through these to draw segment and text labels for each family
# labelling families not orders as some orders in df paraphyletic
family_labels <- function(xpos = 200){
  tip_pos <- pollination1209 %>%
    dplyr::group_by(ParFamTax) %>%
    dplyr::summarise(y0 = min(position), y1 = max(position))
  for (i in 1:nrow(tip_pos)){
    segments(x0 = xpos, y0 = tip_pos$y0[[i]], x1 = xpos, y1 = tip_pos$y1[[i]], lwd = 2)
    text(x = xpos+3, y = tip_pos$y0[[i]] + ((tip_pos$y1[[i]] - tip_pos$y0[[i]])/2), 
         paste(tip_pos$ParFamTax[[i]], sep = ""), adj = 0, srt = 0, cex = 0.4)
  }
  rm(tip_pos)
}
# will have to change some of these labels manually as two families paraphyletic

# function to work out highest and lowest tip numbers for each order,
# then loop through these to draw segment and text labels for each order
order_labels <- function(xpos = 200){
  tip_pos <- pollination1209 %>%
    dplyr::group_by(ParOrdTax) %>%
    dplyr::summarise(y0 = min(position), y1 = max(position))
  for (i in 1:nrow(tip_pos)){
    segments(x0 = xpos, y0 = tip_pos$y0[[i]], x1 = xpos, y1 = tip_pos$y1[[i]], lwd = 2)
    text(x = xpos + 3, y = tip_pos$y0[[i]] + ((tip_pos$y1[[i]] - tip_pos$y0[[i]])/2), 
         paste(tip_pos$ParOrdTax[[i]], sep = ""), adj = 0, srt = 0, cex = 0.4)
  }
  rm(tip_pos)
}