library(data.table)
source("R/globals.R") ## global variables and functions
 
#' @title Obtain a data set according to a study name
#' @description \code{getDataset} requests a particular dataset
#' according to user specified study
#'
#' @param study [vector] vector of studies from which to obtain data
#'
#' @param summary [boolean] optional parameter of whether or not to return summary data 
#' if summary=TRUE, summary data is returned
#' if summary=FALSE, dose response data is returned
#' default is summary=FALSE
#'
#' @param stats [vector] vector of summary statistics values to return
#' default is NULL because default for summary is FALSE
#' 
#' @return list of data.frame containing data if found or NA if data does not exist
#' 
#' @examples
#' getDataset("CCLE") ## get the CCLE dose-response data
#' getDataset(c("CGP", "GNE")) ## get CGP and GNE dose-response data
#' getDataset("CTRP", summary=TRUE) ## get all summary statistics for CTRP data
#' getDataset("CGP", stats = c("IC50_published")) ## get IC50 values in CGP 
#' @export

getDataset <- function(study,summary=FALSE,stats=NULL) {

  ## Construct URLs
  if (!is.null(stats)) {
    summary = TRUE
  }

  if (summary) {
    data <- sapply(study, function (x) paste0(host, "/data/study/summary/", x, "/"))
    if (!is.null(stats)) {
      data <- sapply(data, function (val) sapply(stats, function (x) paste0(val, x, "/")))
    }
  }
  else {
    data <- sapply(study, function (x) paste0(host, "/data/study/", x, "/"))
  }

  ## Query URLs and then parse data
  dfs <- lapply(data, function (x) {
    if (url.exists(x)) {
      json <- getURL(x)
      data <- processJSONData(json)
      return(data)
    }
    else {
    warning(paste("Data at", x, "does not exist"))
    return(NA)
    }
  })

  ## Quick fix for names of summary data
  if (!is.null(stats)) {
    names(dfs) <- stats
  }

  return(dfs)
}

## Helper function to change json strings to data.frames
## json: json.string containing dose-response data or summary data
processJSONData <- function(json) {
  data <- processJSON(json, "data")
  data <- lapply(data, data.frame)
  completeCases <- max(sapply(data, ncol))  ## Get complete cases
  data <-  lapply(data, function (x) if (ncol(x) == completeCases) { x }) ## Remove incomplete observations
  data <- rbindlist(data)
  return(data)
}

#' @title Obtain a data set according to a study name
#' @description \code{getExpData} requests a particular dataset
#' according to user specified study
#'
#' @param cellline [vector] vector of cell lines for which to obtain data
#' if cellline=NULL and drug != NULL, then data for all celllines 
#' tested on that vector of drugs is returned.
#'
#' @param drug [vector] vector of drugs for which to obtain data
#' if drug=NULL and cellline != NULL, then data for all
#' drugs tested on the cellline vector returned
#'
#' if drug!=NULL and cellline !=NULL, then data for all
#' experiments in the set cellline, drug pairs for which
#' there is data is returned
#'
#' @param summary [boolean] whether or not to return summary data 
#' if summary=TRUE, summary data is returned
#' if summary=FALSE, dose response data is returned
#' default is summary=FALSE
#'
#' @param stats [vector] vector of summary statistics values to return
#' default is stats=NULL because summary=FALSE
#' 
#' @return list of data.frame containing data if found or NA if data does not exist
#' 
#' @examples
#' getExpData(cellline="HCC70") ## get all dose-response curves tested on HCC70
#' ## Get the published IC50 values for all cell lines tested on Erlotinib or 17-AAG
#' getExpData(drug="Erlotinib", stats = "IC50_Published")
#' ## Get all summary statistics for experiments tested on MCF7 and 1321N1
#' ## and with erlotinib or 17-AAG
#' getExpData(cellline = c("MCF7", "1321N1"), drug = c("erlotinib", "17-AAG"), summary = TRUE)
#'
#' @export
getExpData <- function(cellline=NULL, drug=NULL, summary=FALSE, stats=NULL) {

  ## Check parameters used properly
  if (is.null(cellline) && is.null(drug)) {
    stop("No data requested!")
  }

  if (!is.null(stats)) {
    summary=TRUE
  }

  ## Construct the URLs 
  if (is.null(cellline)) {
    url <- constructDrugURL(drug, summary)
  }
  if (is.null(drug)) {
    url <- constructCellURL(cellline, summary)
  }
  if (!is.null(drug) & !is.null(cellline)) {
    url <- constructExperimentURL(cellline, drug, summary)
  }

  ## Query URLs and parse
  dfs <- lapply(url, function (x) {
    if (url.exists(x)) {
      json <- getURL(x)
      data <- processJSONData(json)
      return(data)
    }
    else {
    warning(paste("Data at", x, "does not exist"))
    return(NA)
    }
  })

  ## Filter if stats != NULL
  if (!is.null(stats)) {
    dfs <- filterByType(dfs, stats)
  }


  return(dfs)
}

## Helper functions to construct drug URLs
## drug: vector of drugs 
## summary: whether or not summary statistics are to be returned
constructDrugURL <- function(drug, summary) {
  url <- sapply(drug, function (x) paste0(host, "/data/drug/", x, "/"))
  if (summary) {
    url <- sapply(url, function (x) paste0(x, "summary/"))
  }
  return(url)
}

## Helper function to construct cell URLs
## cell: vector of cellline
## summary: whether or not summary statistics are to be returned
constructCellURL <- function(cellline, summary) {
  url <- sapply(cellline, function (x) paste0(host, "/data/cell/", x, "/"))
  if (summary) {
    url <- sapply(url, function (x) paste0(x, "summary/"))
  }
  return(url)
}


## Helper function to construct experiment URLs
## cell: vector of cellline
## drug: vector of drugs 
## summary: whether or not summary statistics are to be returned
constructExperimentURL <- function(cellline, drug, summary) {
  exp <- sapply(cellline, function (x) sapply(drug, function (y) paste0(x, "/", y)))
  url <- sapply(cellline, function (x)
                           sapply(drug, function (y) paste0(host, "/data/", x, "/", y, "/")))
  if (summary) {
    url <- sapply(url, function (x) paste0(x, "summary/"))
  }
  names(url) <- exp 
  return(url)
}

## filter obtained data by statistic types
## frames: list of curation
## stats: statistics to keep in curation
filterByType <- function(frames, stats) {
  dfs <- lapply(frames, function (frame) {
    if (is.null(nrow(frame))) {
      return(NA)
    }
    else {
    frame[frame$stat_name %in% stats, ]
    }})
  return(dfs)
}


