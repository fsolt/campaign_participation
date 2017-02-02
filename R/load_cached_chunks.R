setwd("paper")
text <- readLines("campaign_participation.Rnw")
chunks <- grep("^<<.*>>=", text, value = T)
extractchunks <- gsub("^<<label=([^,>]*)[,>]*.*", "\\1", chunks)[-1]

envir <- parent.frame()
oldls <- ls(envir, all = TRUE)

# For each chunk...
for(ch in extractchunks){   
  # Detect the file name of the database...
  pat <- paste0("^", ch, "_.*\\.rdb")
  val <- gsub(".rdb", "", dir("cache", pattern = pat))
  # Lazy load the database
  lazyLoad(file.path("cache", val), envir = envir)
}
setwd("..")
