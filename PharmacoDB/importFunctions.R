library(RMySQL)
source("processCSV.R")

## Connect to PharmacoDB
drv <- dbDriver("MySQL")

## Mordor connection 
# rmysql.settingsfile <- "/mnt/work1/users/home2/bhcoop4/.my.cnf" # TODO: Change dir of R config fil
# connection <- dbConnect(drv, group = "pharmacodb")

## Local connection
connection <- dbConnect(drv, username = "root", password = "", host = "127.0.0.1", db = "PharmacoDB")

## Write a study vec into the STUDIES table
doStudiesQuery  <- function (study) {
 query <- paste0("INSERT INTO STUDIES (study_name) VALUES (\'", study, "\')")
 dbGetQuery(connection, query) 
}

## function to get a study ID from name of study
getStudyID <- function(study) {
	val <- dbGetQuery(connection, paste0("SELECT study_id FROM STUDIES WHERE study_name = \"", study, "\""))
	if (nrow(val) == 0) {
		stop(paste(study, "does not exist in the database!"))
	}
	return(val)
}

## Write a metadata item into METADATA table
doMetadataQuery  <- function (item) {
 query <- paste0("INSERT INTO METADATA (metadata_name) VALUES (\'", item, "\')")
 dbGetQuery(connection, query) 
}

## Write to cell tables from cell_annotation_all.csv
writeToCellTables <- function(data, study, toAppend = TRUE) {
	studyID <- getStudyID(study)
	frame <- processCell(data, study, studyID)
	message(paste("Importing", study, "to CELL_LINE_CURATION"))
	dbWriteTable(connection, "CELL_LINE_CURATION", frame, append = TRUE, row.names = FALSE)
}

## Write To tissues tables
writeToCellTissues <- function(data, study) {
	studyID <- getStudyID(study)
	metadataID <- dbGetQuery(connection, paste0("SELECT metadata_id FROM METADATA WHERE metadata_name = \'tissue\'"))
	frame <-  processCellTissues(data, study, studyID, metadataID)
	message(paste("Importing", study, "tissues to CELL_LINE_METADATA"))
	dbWriteTable(connection, "CELL_LINE_METADATA", frame, append = TRUE, row.names = FALSE)
}

## Write to drug tables
writeToDrugTables <- function(drugs, study) {
	## Use individual drug annotations
	## fileName <- paste0("drug_annotation_", study, ".csv")
	## if (!file.exists(fileName)) {
	##	stop(paste0("drug_annotation_", study, "does not exist!"))
	## }
	## drugs <- read.csv(fileName, na.strings = "") ## read file
	studyID <- getStudyID(study)
	frame <- processDrugs(drugs, study, studyID)
	message(paste("Importing", study, "to DRUG_CURATION"))
	dbWriteTable(connection, "DRUG_CURATION", frame, append = TRUE, row.names = FALSE)
}

## Write to tissue tables
writeToTissueTables <- function(data, study) {
    studyID <- getStudyID(study)
    frame <- processTissue(data, study, studyID)
    message(paste("Importing", study, "to TISSUE_TYPE"))
    dbWriteTable(connection, "TISSUE_TYPE", frame, append = TRUE, row.names = FALSE)
}

## Read data and import
writeToSensitivityTables <- function(study) {
	fileName <- paste0(tolower(study), "_sensitivity_detail.csv")
	if (!file.exists(fileName)) {
		stop(paste0(tolower(study), "_sensitivity_detail.csv does not exist!"))
	}
	sensitivityCSV <- read.csv(fileName, stringsAsFactors = FALSE)
	studyID <- getStudyID(study)

	if (study == "GRAY") { ## GRAY has experimental replicates in the dataset and 
		## thus must be processed differently
		frame <- processDoseResponseByRow(sensitivityCSV, studyID)
	}
	else {
		frame <- processDoseResponseData(sensitivityCSV, studyID)
	}
    message(paste("Importing", study, "to SENSITIVITY_DATA"))
	dbWriteTable(connection, "SENSITIVITY_DATA", frame, append = TRUE, row.names = FALSE)
}

## function to write an element of studies vec into the table
doStatQuery  <- function (item) {
    query <- paste0("INSERT INTO SUMMARY_STATISTICS (stat_name) VALUES (\'", item, "\')")
    dbGetQuery(connection, query)
}

## Get stat ids from a list of 
getStatIDs <- function (stats) {
	statIDs <- c()
   for (i in 1:length(stats)) {
	val <- dbGetQuery(connection, paste0("SELECT stat_id FROM SUMMARY_STATISTICS WHERE stat_name = \"", stats[i], "\""))
	if (nrow(val) == 0) {
    		stop(paste(stat, "does not exist in SUMMARY_STATISTICS!"))
    }
    statIDs[i] <- val
	}
    return(statIDs)
}

## function to read drug summary data for a study
writeToSummaryTables <- function(study) {
	fileName <- paste0(tolower(study), "_sensitivity.csv")
	if (!file.exists(fileName)) {
		stop(paste0(tolower(study), "_sensitivity.csv does not exist!"))
	}

	summaryCSV <- read.csv(fileName, stringsAsFactors = FALSE) ## read file
    studyID <- getStudyID(study)

    if ("X" %in% colnames(summaryCSV)) { ## remove extraneous row of row names
        summaryCSV <- summaryCSV[,-1]
    }

    stats <- colnames(summaryCSV)[-1:-2] ## get all names of stats in CSV w/o cellline or drugs
    statIDs <- getStatIDs(stats) ## get corresponding IDs
       
	if (study == "GRAY") { ## GRAY has experimental replicates in the dataset and 
		## thus must be processed differently
		frame <- processSummaryDataByRow(summaryCSV, studyID, statIDs, stats)
	}
	else {
		frame <- processSummaryData(summaryCSV, studyID, statIDs, stats)
	}
    message(paste("Importing", study, "to SUMMARY_STATISTICS_DATA"))
    dbWriteTable(connection, "SUMMARY_STATISTICS_DATA", frame, append = TRUE,
    row.names = FALSE)
}
