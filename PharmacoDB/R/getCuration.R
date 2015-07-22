library(data.table)
source("R/globals.R") ## global variables and functions
 
#' @title Obtain a curation
#' @description \code{getCuration} requests a particular curation 
#'  according to user specified parameters from PharmacoDB web service.
#' @param curation the curation to be requested: 
#' 
#' "cell" is listing of cell lines in each study and corresponding synonyms 
#' 
#' "drug" is listing of drugs in each study and corresponding synonyms
#' "tissue" is listing of cell lines in each study, their corresponding tissue types, 
#' and matching of tissue types in each study to COSMIC tissue types
#'
#' "histology" is listing of cell lines in each study, their corresponding histological subtypes,  
#'  matching between subtypes based on manual curation by BHK lab
#'
#' "gene_expression" is listing of gene expression metadata by study
#' 
#' @param study optional vector of studies to subset curation
#' '
#' @param subtype optional vector of subtype so only curation for that subtype is displayed
#'        only valid for "tissue", "histology" and "gene expression" curations
#'
#' @return list containing curations if found or NA if curation not found
#' 
#' @examples
#' getCuration("cell") 
#' ## get the complete cell line curation
#' getCuration("drug", study=c("CGP", "CCLE")) 
#' ## get all drugs in CGP and in CCLE
#' getCuration("tissue", study=c("CCLE", "GRAY", "GNE"), subtype=c("breast", "testis")) 
#' ## get all breast and testis cancer cell lines tested in CCLE and GRAY
#' getCuration("histology", subtype = "carcinoma") 
#' ## get all cell lines which are carcinomas
#' getCuration("gene_expression", "MEIS")
#' ## returns warning since "MEIS" study does not exist
#' @export

getCuration <- function(curation = c("cell", "drug", "tissue", "histology", "gene_expression"), study = NULL, subtype = NULL) {

  curation <- match.arg(curation)
  
  ## Check that subtype parameter used properly
  if (!is.null(subtype) & curation %in% c("cell", "drug")) {
    stop("Subtype or profile should not be specified for cell or drug curation.")
  }
  
  ## Case with no subtypes - specified curation and optional study parameter
  if (is.null(subtype)) {
    curationFrames <- getCurationByStudy(curation, study)
  }

  ## Case with a subtype 
  if (!is.null(subtype)) {
    curationFrames <- getCurationByType(curation, subtype)
    if (!is.null(study)) { ## Filter by study here
      curationFrames <- filterByStudy(curationFrames, study)
    }
  }
  return(curationFrames)
}

## Helper function to process a JSON curation and convert it into a dataframe
## json: json string
processJSONCuration <- function(json) {
  df <- processJSON(json, "curation")
  df <- do.call("rbind", df)
  df <- data.frame(df)
  return(df)
}

## Obtain a curation by study 
## type: type of curation
## study: studies in which curation is requested
getCurationByStudy <- function(type, study) {
  if (is.null(study)) {
    url <- paste0(host, "/curation/", type, "/")
    json <- getURL(url)
    df <- processJSONCuration(json)
  }
  else {
    url <- sapply(study, function (x) paste0(host, "/curation/", type, "/", x, "/"))
    df <- lapply(url, function (x) {
        if (url.exists(x)) {
          json <- getURL(x)
          frame <- processJSONCuration(json)
          return(frame)
        }
        else {
          warning(paste("Curation at", x, "does not exist"))
          NA
        }
      })
  }
  return(df)
}

## Obtain a curation by type
## type: type of curation
## subtype: subtypes requested
getCurationByType <- function(type, subtype) {
  url <- sapply(subtype, function (x) paste0(host, "/curation/", type, "/type/", x, "/"))
  df <- lapply(url, function (x) {
        if (url.exists(x)) {
          json <- getURL(x)
          frame <- processJSONCuration(json)
          return(frame)
        }
        else {
          warning(paste("Curation at", x, "does not exist"))
          NA
        }
      })
  return(df)
}

## filter obtained curations by type 
## frames: list of curation
## study: studies to keep in curation
filterByStudy <- function(frames, study) {
  curationFrames <- lapply(frames, function (frame) {
    if (is.null(nrow(frame))) {
      return(NA)
    }
    else {
    frame[frame$study_name %in% study, ]
    }})
  return(curationFrames)
}

