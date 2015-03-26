library(RMySQL)
source("processCSV.R")

## Connect to PharmacoDB
# rmysql.settingsfile <- "/mnt/work1/users/home2/bhcoop4/.my.cnf" # TODO: Change dir of R config file
drv <- dbDriver("MySQL")
# connection <- dbConnect(drv, group = "pharmacodb")
connection <- dbConnect(drv, username = "root", password = "", host = "127.0.0.1", db = "PharmacoDB")

## Construct STUDIES Table
message("Writing to STUDIES table")
studies <- c("CCLE", "CGP")

## function to write an element of studies vec into the table
 doStudiesQuery  <- function (study) {
 query <- paste0("INSERT INTO STUDIES (study_name) VALUES (\'", study, "\')")
 dbGetQuery(connection, query) 
}

## import studies vec into STUDIES table
sapply(studies, doStudiesQuery)

## Construct METADATA table
message("Writing to METADATA table")
metadata <- c("tissue") 

## function to write an element of metadata into the table
doMetadataQuery  <- function (item) {
 query <- paste0("INSERT INTO METADATA (metadata_name) VALUES (\'", item, "\')")
 dbGetQuery(connection, query) 
}

## import metadata vec into METADATA table
sapply(metadata, doMetadataQuery)

## Construct cell line curation and cell line metadata tables
message("Writing to CELL_LINE_CURATION and CELL_LINE_METADATA table")

## Read data from "cell_line_annotation_all_ATTENTION.csv"
if (!file.exists("cell_line_annotation_all_ATTENTION.csv")) {
	stop("Cell line curations do not exist.")
}
cells <- read.csv("cell_line_annotation_all_ATTENTION.csv", na.strings = "")
writeToCellTables <- function(study) {
	message(paste("Processing", study))
	studyID <- dbGetQuery(connection, paste0("SELECT study_id FROM STUDIES WHERE study_name = \"", study, "\""))
	frame <- processCell(cells, study, studyID)
	message(paste("Importing", study, "to CELL_LINE_CURATION"))
	dbWriteTable(connection, "CELL_LINE_CURATION", frame, append = TRUE, row.names = FALSE)
	metadataID <- dbGetQuery(connection, paste0("SELECT metadata_id FROM METADATA WHERE metadata_name = \'tissue\'"))
	frame <-  processCellTissues(cells, study, studyID, metadataID)
	dbWriteTable(connection, "CELL_LINE_METADATA", frame, append = TRUE, row.names = FALSE)
}

## write to the cell tables
sapply(studies, writeToCellTables)

## Construct DRUG_CURATION TABLE
message("Writing to DRUG_CURATION table")

## Read data and import
writeToDrugTables <- function(study) {
	message(paste("Processing", study))
	fileName <- paste0("drug_annotation_", study, ".csv")
	if (!file.exists(fileName)) { 
		stop(paste0("drug_annotation_", study, "does not exist!")) 
	}
	drugCSV <- read.csv(fileName)
	studyID <- dbGetQuery(connection, paste0("SELECT study_id FROM STUDIES WHERE study_name = \"", study, "\""))
	frame <- processDrugs(drugCSV, studyID)
	dbWriteTable(connection, "DRUG_CURATION", frame, append = TRUE, row.names = FALSE)
}

## Import data into DRUG_CURAITON table
sapply(studies, writeToDrugTables)

## Construct SENSITIVITY_DATA TABLE
message("Writing to SENSITIVITY_DATA")

## Read data and import
writeToSensitivityTables <- function(study) {
	message(paste("Processing", study))
	fileName <- paste0(tolower(study), "_sensitivity_detail.csv")
	if (!file.exists(fileName)) {
		stop(paste0(tolower(study), "_sensitivity_detail.csv does not exist!"))
	}
	sensitivityCSV <- read.csv(fileName)
	studyID <- dbGetQuery(connection, paste0("SELECT study_id FROM STUDIES WHERE study_name = \"", study, "\""))
	frame <- processDoseResponseData(sensitivityCSV, studyID)
	dbWriteTable(connection, "SENSITIVITY_DATA", frame, append = TRUE, row.names = FALSE)
}

## Import data into table
sapply(studies, writeToSensitivityTables)

## Construct SUMMARY_STATISTICS
message("Writing to SUMMARY_STATISTICS")
stats <- c("AUC", "AUC.", "slope", "slope.")

## function to write an element of studies vec into the table
doStatQuery  <- function (item) {
    query <- paste0("INSERT INTO SUMMARY_STATISTICS (stat_name) VALUES (\'", item, "\')")
    dbGetQuery(connection, query)
}

## import each element of studies into the table
sapply(stats, doStatQuery)

## Construct SUMMARY_STATISTICS_DATA
message("Writing to SUMMARY_STATISTICS_DATA")

## function to read drug summary data for a study
writeToSummaryTables <- function(study) {
	message(paste("Processing", study))
	fileName <- paste0(tolower(study), "_sensitivity.csv")
	if (!file.exists(fileName)) {
		stop(paste0(tolower(study), "_sensitivity.csv does not exist!"))
	}
	summaryCSV <- read.csv(fileName)
    studyID <- dbGetQuery(connection, paste0("SELECT study_id FROM STUDIES WHERE study_name = \"", study, "\""))
    if ("X" %in% colnames(summaryCSV)) { ## remove extraneous row of row names
        summaryCSV <- summaryCSV[,-1]
    }
    stats <- colnames(summaryCSV)[-1:-2] ## remove cellline/drugname columns of original CSV
       statIDs <- c()
    ## Query SUMMARY_STATISTICS table to obtain  ids
    for (i in 1:length(stats)) {
    	val <- dbGetQuery(connection, paste0("SELECT stat_id FROM SUMMARY_STATISTICS WHERE stat_name = \"", stats[i], "\""))
    	if (length(val)== 0) {
    		stop(paste(stats[i], "does not exist in SUMMARY_STATISTICS!"))
    	}
    	statIDs <- c(statIDs, val)
    }
    frame <- processSummaryData(summaryCSV, studyID, statIDs, stats)

    dbWriteTable(connection, "SUMMARY_STATISTICS_DATA", frame, append = TRUE,
    row.names = FALSE)
}

## Import data into table
sapply(studies, writeToSummaryTables)