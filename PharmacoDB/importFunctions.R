## Functions to process raw csv and import them into PHarmacoDB
## Author: Adrian She
## Last Updated: June 11, 2015

library(RMySQL)

## Mordor connDB() 
rmysql.settingsfile <- "/mnt/work1/users/home2/bhcoop4/.my.cnf"

## Local connection
connDB <- function() {
 #  con <- dbConnect(MySQL(), username = "root", password = "", host = "127.0.0.1", db = "PharmacoDB")
  con <- dbConnect(MySQL(), group = "pharmacodb")
  return(con)
}

## === STUDIES FUNCTIONS === ##

## Write a study vec into the STUDIES table
doStudiesQuery  <- function (study) {
  connection <- connDB()
  dbGetQuery(connection, paste0("INSERT INTO STUDIES (study_name) VALUES (\'", study, "\')")) 
  dbDisconnect(connection)
}

## function to get a study ID from name of study
getStudyID <- function(study) {
  connection <- connDB()
  val <- dbGetQuery(connection, paste0("SELECT study_id FROM STUDIES WHERE study_name = \"", study, "\""))
  if (nrow(val) == 0) {
    stop(paste(study, "does not exist in the database!"))
  }
  dbDisconnect(connection)
  return(as.numeric(val))
}

## ==== CELL TABLE FUNCTIONS ==== ###

## Write a cell metadata item into cell line metadata table
doCellMetadataQuery  <- function (item) {
  connection <- connDB()
  dbGetQuery(connection, paste0("INSERT INTO CELL_METADATA_IDX (metadata_name) VALUES (\'", item, "\')")) 
  dbDisconnect(connection)
}

## function to get a metadata ID from name of metadata
getCellMetadata <- function(metadata) {
  connection <- connDB()
  val <- dbGetQuery(connection, paste0("SELECT metadata_id FROM CELL_METADATA_IDX WHERE metadata_name = \"", metadata, "\""))
  if (nrow(val) == 0) {
    stop(paste(metadata, "does not exist in the database!"))
  }
  return(as.numeric(val))
  dbDisconnect(connection)
}

## Extract cell line names and synonyms for a particular study as a data frame
## data: read csv file for "cell_line_annotation.csv"
## study: name of the study
## study_id: integer identifying the study name
processCell <- function(data, study) { 
  message(paste("Processing cell names for", study))
  frame <- data.frame("cellline_id" =  data[,"unique.cellid"], 
                      "cell_study_id" = getStudyID(study), 
                      "cellline_name" = data[,paste0(toupper(study), ".cellid")]) # construct frame of cell line info

  frame <- frame[complete.cases(frame),] # remove NAs
  rownames(frame) <- NULL  # remove row.names caused by subsetting in complete.cases
  message(paste("Importing", study, "to CELL_LINE_CURATION"))
  connection <- connDB()
 dbWriteTable(connDB(), "CELL_LINE_CURATION", frame, append = TRUE, row.names = FALSE)
 dbDisconnect(connection)
 }

## Extract tissue types for cell lines in a particular study as a data frame
## data: data frame which stores "cell_line_annotation_all_ATTENTION.csv"
## study: name of the study
## study_id: integer identifying the study name
processCellTissues <- function(data, study, tissueCor) { 
  message(paste("Processing tissues type for", study))
  frame <- data.frame("cellline_id" =  data[,"unique.cellid"], 
                      "study_id" = getStudyID(study), 
                      "tissue_name" = data[,paste0(toupper(study), ".tissueid")]) # construct frames for cell line info

  frame <- frame[complete.cases(frame),] ## remove NAs
  # rename according to colnames of table in the schema
  rownames(frame) <- NULL # remove row.names caused by subesetting 
  frame <-  matchTissue(tissueCor, frame, study)
  message(paste("Importing", study, "tissues to TISSUE_CURATION"))
  connection <- connDB()
  dbWriteTable(connection, "TISSUE_CURATION", frame, append = TRUE, row.names = FALSE)
  dbDisconnect(connection)
 }

## Add matching tissues
## tissue_match : data.frame of matching tissue
## frame: the data frame which contains non-unique tissue types
## study: name of study
matchTissue <- function(tissue_match, frame, study) {
  if (paste0(study, ".tissueid") %in% colnames(tissue_match)) {
    col <- paste0(study, ".tissueid")
    tissues <- as.vector(t(frame$tissue_name))
    frame$unique_tissue_name <- 
      sapply(tissues, function (x) {
        paste(tissue_match[which(toupper(tissue_match[, col]) == toupper(x)), "COSMIC.tissueid"], collapse = '/')
      })
  }
  return(frame)
}

## == DRUG CURATION FUNCTIONS ==== 

## Extract drugs used in a particular study  as a data frame
## drugs: read CSV file for "drug_annotation_STUDYNAME.csv"
## study_id: integer identifying the study
processDrugs <- function(drugs, study) {
  message(paste("Processing drug names for", study))
  
  ## Construct frame
  frame <- data.frame("drug_id" = drugs[,"unique.drugid"], 
                          "drug_study_id" = getStudyID(study),
                          "drug_name" = drugs[, paste0(toupper(study), ".drugid")])
  frame <- frame[complete.cases(frame),]
	rownames(frame) <- NULL
  message(paste("Importing", study, "drug to DRUG_CURATION"))
  connection <- connDB()
  dbWriteTable(connection, "DRUG_CURATION", frame, append = TRUE, row.names = FALSE)
  dbDisconnect(connection)
}

## Write a drug metadata item into cell line metadata table
doKeyQuery  <- function (item) {
  connection <- connDB()
  dbGetQuery(connection, paste0("INSERT INTO DRUG_METADATA_IDX (key_name) VALUES (\'", item, "\')")) 
  dbDisconnect(connection)
}

## function to get a metadata ID from name of metadata
getDrugMetadataID <- function(metadata) {
  connection <- connDB()
  val <- dbGetQuery(connDB(), paste0("SELECT key_id FROM DRUG_METADATA_IDX WHERE key_name = \"", metadata, "\""))
  if (nrow(val) == 0) {
    stop(paste(metadata, "does not exist in the database!"))
  }
  dbDisconnect(connection)
  return(as.numeric(val))
}

## Write drug metadata into drug metadata table
processDrugsMetadata <- function(data) {
  frame <- data.frame("drug_id" = data[,"unique.drugid"],
                      "smiles" = data[,"smiles"],
                      "inchikey" = data[, "inchikey"],
                        "cid" = data[,"cid"])
  message("Importing Drug Metadata")
  connection <- connDB()
  dbWriteTable(connection, "DRUG_METADATA", frame, overwrite = TRUE, row.names = FALSE)
  dbDisconnect(connection)
}

## ==== SENSITIVITY DATA AND SUMMARY DATA ===== ##


## Convert read CSV file of dose-response sensitivity data 
## into format which can be inputted into PharmacoDB as a data frame
## data: read .csv frame with original data 
## study_id: integer which identifes the study
processDoseResponseData <- function(data, study) {
  ## get number of observations
  firstDoseCol <- grep("log10_dose_1$", colnames(data)) 
  firstViabilityCol <- grep("viability_1$", colnames(data)) 

  study_id <- getStudyID(study)
  doseCols <- data[,seq(firstDoseCol, firstViabilityCol - 1, 1)]
  viabilityCols <- data[, seq(firstViabilityCol, ncol(data), 1)]
  
  frame <- data.frame()
  for (i in 1:ncol(doseCols)) {
    message(paste("Processing Observation", i))
    newFrame <- data.frame(
      "cellline_id" = data$cellline,
      "drug_id" = data$drug,
      "study_id" = study_id,
      "log10dose" = doseCols[,i],
      "viability" = viabilityCols[,i])
    frame <- rbind(frame, newFrame)
  }
  
  rownames(frame) <- NULL
  
  ## Write Frame
  message(paste("Importing", study, "to DOSE_RESPONSE_DATA"))
  connection <- connDB()
  dbWriteTable(connection, "DOSE_RESPONSE_DATA", frame, append = TRUE,
               row.names = FALSE)
  dbDisconnect(connection)
}

## function to write an element of studies vec into the table
doStatQuery  <- function (item) {
  connection <- connDB()
  dbGetQuery(connection, paste0("INSERT INTO SUMMARY_STATISTICS (stat_name) VALUES (\'", item, "\')"))
  dbDisconnect(connection)
}

## Get stat ids from a list of 
getStatIDs <- function (stats) {
  connection <- connDB()
  vals <- sapply(stats, function (x) {
    val <- dbGetQuery(connection,
               paste0("SELECT stat_id FROM SUMMARY_STATISTICS WHERE stat_name = \"", x, "\""))
  })
  dbDisconnect(connection)
  return(as.numeric(vals))
}


## Convert read CSV file of summary statistics
## into format which can be inputted into PharmacoDB as a data frame
## data: read .csv file with the original data
## study_id: integer which identifies the current study
## stat_id: integer vector which identifies summary statistics
## stat_name: names of the summary statistics being considered
processSummaryData <- function(data, study) {
  
  names <- colnames(data)[-1:-3]
  ids <- getStatIDs(names)
  study_id <- getStudyID(study)
  ## create a empty data frame for processed data 
  frame <- data.frame()

  ## Iterate through each statistic
  for (i in 1:length(ids)) {
        message(paste("Processing", names[i], "data"))
        frame <- rbind(frame, 
                       data.frame(
                         "cellline_id" = data$cellline,
                         "drug_id" = data$drug,
                         "stat_id" = ids[i],
                         "study_id" = study_id,
                         "stat_value" = data[,names[i]]))
    }
   ## return frame
   rownames(frame) <- NULL
   message(paste("Importing", study, "to SUMMARY_STATISTICS_DATA"))
   connection <- connDB()
   dbWriteTable(connection, "SUMMARY_STATISTICS_DATA", frame, append = TRUE,
                 row.names = FALSE)
   dbDisconnect(connection)
}
  




