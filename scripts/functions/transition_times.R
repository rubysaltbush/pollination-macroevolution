# function to extract transition times from simulation mapping

transition_times <- function(simmap){
  # below copied from phytools blog 
  # http://blog.phytools.org/2015/08/getting-timing-of-trait-changes-from.html
  # extracts raw transition times from a simmap (collapses multiple transitions
  # down into single transition events e.g. 4->2->4 becomes ?????)
  # get tips and their states
  x <- phytools::getStates(simmap,"tips")
  # get unique states
  states <- sort(unique(x))
  # get length of states
  m <- length(states)
  # below makes a little matrix describing transitions
  ct <- sapply(states, 
               function(x,y) sapply(y, function(y,x) paste(x,"->", y, sep=""), 
                                    x = x), y = states)
  rm(x, states)
  # then a matrix to invalidate self->self transitions
  ii <- matrix(TRUE, m, m)
  diag(ii) <- rep(FALSE, m)
  # then a list to store results in
  changes <- vector(mode="list", length = m*(m - 1))
  rm(m)
  # named by types of transitions
  names(changes) <- as.vector(ct[ii])
  rm(ct, ii)
  # then singling out maps where transitions happen (where there is more than 1 state)
  nc <- sapply(simmap$maps, length) - 1
  ind <- which(nc > 0)
  nc <- nc[ind]
  
  # It would be nice to return the position where these changes are happening
  # and maybe this could be gotten here, but too late/tired to think of it now

  # getting the node heights (measure of time) across the tree
  H <- phytools::nodeHeights(simmap)
  maps <- simmap$maps[ind]
  # then looping through and calculating the node heights (i.e. time) of each
  # transition
  for(i in 1:length(maps)){
    for(j in 1:nc[i]){
      sc <- paste(names(maps[[i]])[j:(j + 1)], collapse = "->")
      h <- H[ind[i], 1] + cumsum(maps[[i]])[j]
      changes[[sc]] <- c(changes[[sc]], as.numeric(h))
    }
  }
  rm(nc, ind, h, H, i, j, sc, maps)
  # removing any nulls from list of changes and sorting small to large
  changes <- changes[!sapply(changes, is.null)]
  changes <- lapply(changes, sort, decreasing = TRUE)
  
  # now convert this changes list into nice data frame output RUBY WROTE THIS BIT
  output <- data.frame()
  for(i in 1:length(changes)){
    df <- dplyr::bind_cols(changes[i])
    df <- df %>%
      mutate(transition = colnames(df)) %>%
      rename(time = 1)
    output <- rbind(output, df)
  }
  
  # and return the output! to graph etc.
  output
}

# BUT THE NUMBERS IN "TIME" COLUMN ARE NOT ACTUALLY TIMES - THEY ARE "MAPS"!!! HOW TO CONVERT???
