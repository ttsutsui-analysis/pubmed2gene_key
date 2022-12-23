library(shiny)
library(pbapply)
library(future.apply)
library(progressr)

# efetch using 18.2
fetch <- function(term){
  options(warn=2)
  com <- paste0("esearch -db pubmed -query '", term, "' | efetch -format uid")
  ids <- system(com, intern = T)
  options(warn=1)
  return(ids)
}

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

# pickup overlapped pmid for each cut-offed gene
pick.pmid <- function(all, pmid, cutoff){
  all.cut <- all[cutoff]
  count <- pbsapply(all.cut, function(x){
    tmp <- x[x %in% pmid]
    incProgress(1/length(all.cut))
    tmp <- paste0(tmp, collapse = ", ")
    return(tmp)
  })
  return(count)
}
