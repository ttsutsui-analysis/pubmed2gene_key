library(shiny)
library(pbapply)
library(future.apply)
library(progressr)

# count overlaps
count.sapply <- function(all, ids){
  x <- pbsapply(all, function(x){
    incProgress(1/length(all))
    sum(x %in% ids)
  })
  return(x)
}

# parallel
count.para <- function(all, pmid, cores){
  cores <- as.integer(cores)
  plan(multisession, workers=cores)
  p <- progressr::progressor(100)
  count <- future_sapply(seq_along(all), function(x){
    y <- sum(pmid %in% all[[x]])
    if(x %% floor(length(all)/100)==0){p()}
    return(y)
  })
  plan(sequential)
  return(count)
}
