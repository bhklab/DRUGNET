source("R/globals.R")
library(data.table)

#' @title Obtain all cell lines tested on a set of drugs
#' @description Same as title
#'
#' @param drug [vector] vector of drugs to search
#' @return list of data.frame or NA containing search results
#' @examples
#' getAllCells(c("Erlotinib", "17-AAG"))
#' @export

getAllCells <- function(drug) {
  urls <- lapply(drug, function (x) paste0(host, "/test/cell/", x, "/"))
  vals <- lapply(urls, function (x) {
  	if (url.exists(x)) {
  		json <- getURL(x)
  		frame <- rbindlist(lapply(fromJSON(json)$cells, function (frame) data.frame("cellline_id" = frame["cellline_id"], "study_name" = frame["study_name"])))
  		return(frame)
  	}
  	return(NA)
  })
  return(vals)
}

#' @title Obtain all drugs tested on a set of cell lines
#' @description Same as title
#'
#' @param cell [vector] vector of cell lines to search
#' @return list of data.frame or NA containing search results
#' @examples
#' getAllDrugs(c("MCF7", "1321N1"))
#' @export
getAllDrugs <- function(cell) {
  urls <- lapply(cell, function (x) paste0(host, "/test/drug/", x, "/"))
  vals <- lapply(urls, function (x) {
    if (url.exists(x)) {
      json <- getURL(x)
      frame <- rbindlist(lapply(fromJSON(json)$drugs, function (frame) data.frame("drug_id" = frame["drug_id"], "study_name" = frame["study_name"])))
      return(frame)
    }
    return(NA)
  })
  return(vals)
}

#' @title Get list of studies in the database
#' @description Same as title
#'
#' @return vector of studies
#' @export
getStudies <- function() {
  json <- getURL(paste0(host, "/studies/"))
  return(fromJSON(json)$studies)
}


#' @title Get list of studies with dose-response data/summary data in the database
#' @description Same as title
#'
#' @return vector of studies
#' @export
getDataStudies <- function() {
  json <- getURL(paste0(host, "/data/studies/"))
  return(fromJSON(json)$studies)
}

#' @title Get list of summary statistics available in a study
#' @description Same as title
#' @param [vector] vector of studies
#' @return list of vector of summary statistics
#' @examples
#' listSummaryStats(c("CCLE", "CGP"))
#' @export
listSummaryStats <- function(study) {
  vals <- lapply(study, function (x){
    url <- paste0(host, "/data/summary/", x, "/")
    if (url.exists(url)) {
      json <- getURL(url)
      return(fromJSON(json)$stats)
    }
    else { NA }
    })
  return(vals)
}

