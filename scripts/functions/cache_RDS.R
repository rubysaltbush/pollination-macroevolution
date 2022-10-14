# function to cache pre-prepared R data. If RDS already in cache will read data
# from cache, if not will run defined script to produce and then cache data

cache_RDS <- function(path, fetch_function, read_function = readRDS, 
                      save_function = saveRDS) {
  if (file.exists(path)) {
    read_function(path)
  } else {
    data <- fetch_function()
    save_function(data, path)
    data
  }
}
