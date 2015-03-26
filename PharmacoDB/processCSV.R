## Functions which process existing CSV files 
## which hold data for PharmacoDB

## Author: Adrian She
## Last Updated: March 23, 2015

library(webchem)

## Extract cell line names and synonyms for a particular study as a data frame
## data: read csv file for "cell_line_annotation_all_ATTENTION.csv"
## study: name of the study
## study_id: integer identifying the study name
processCell <- function(data, study, study_id) { 
  message(paste("Processing cell names for", study))
  study <- toupper(study) # cell line columns are upper case
  uniqueIDs <- data["unique.cellid"]
  studyIDs <- data[paste0(study, ".cellid")]
  frame <- cbind(uniqueIDs, study_id, studyIDs) # construct frame of cell line info
  frame <- frame[complete.cases(frame),] # remove NAs
  cellCurationNames <- c("cellline_id", "study_id", "cellline_name") 
  colnames(frame) <- cellCurationNames
  # rename according to colnames of table in the schema
  rownames(frame) <- NULL  # remove row.names caused by subsetting in complete.cases
  return(frame)
 }

## Extract tissue types for cell lines in a particular study as a data frame
## data: data frame which stores "cell_line_annotation_all_ATTENTION.csv"
## study: name of the study
## study_id: integer identifying the study name
processCellTissues <- function(data, study, study_id, metadata_id) { 
  message(paste("Processing tissues type for", study))
  study <- toupper(study) # cell line columns are upper case
  uniqueIDs <- data["unique.cellid"]
  studyIDs <- data[paste0(study, ".tissueid")]
  frame <- cbind(uniqueIDs, metadata_id, study_id, studyIDs) # construct frames for cell line info
  frame <- frame[complete.cases(frame),] ## remove NAs
  # rename according to colnames of table in the schema
  cellMetadataNames <- c("cellline_id", "metadata_id", "study_id", "metadata_value")
  colnames(frame) <-  cellMetadataNames
  rownames(frame) <- NULL # remove row.names caused by subesetting 
  return(frame)
 }

## Extract drugs used in a particular study  as a data frame
## drugs: read CSV file for "drug_annotation_STUDYNAME.csv"
## study_id: integer identifying the study
processDrugs <- function(drugs, study_id) {
	drugs <- drugs["drugid"] # extract name of the drug
	drugsVec <- as.vector(drugs[complete.cases(drugs),]) # get necessary PubChem IDs
	# construct frame
	cidVec <- convertToCID(drugsVec)
  drugFrame <- data.frame(drugsVec,  study_id, cidVec)
	# rename according to colnames of table in the schema
  drugCurationNames <- c("drug_name", "study_id", "drug_id")
	colnames(drugFrame) <- drugCurationNames
	rownames(drugFrame) <- NULL
	return(drugFrame)
}

## Convert a vector of chemical names to PubChem CID
convertToCID <- function(drugsVec) {
  cidVec <- c()
  for (i in 1:length(drugsVec)) {
    cidVec <- c(cidVec, 
      as.integer(cts_convert(query = drugsVec[i], from = "Chemical Name", to = "PubChem CID"))[1])
  }
  return(cidVec)
}

## Convert read CSV file of dose-response sensitivity data 
## into format which can be inputted into PharmacoDB as a data frame
## data: read .csv file with the original data
## study_id: integer which identifes the study
## TODO: deal with issues which may arise as a result of loss of precision
processDoseResponseData <- function(data, study_id) {
  ## create a empty data frame for processed data 
  doseResponseFrame <- data.frame()
  ## get list of cell lines and drugs used in each study 
  cells <- data["cellline"]
  drugs <- data["drug"]
  ## get number of observations
  colBeforeFirstDose <- grep("log10_dose_1", colnames(data)) - 1
  colBeforeFirstViability <- grep("viability_1", colnames(data)) - 1
  observations <- colBeforeFirstViability - colBeforeFirstDose 
  sensitivityDataNames <- c("cellline_id", "drug_id", "study_id", "log10dose", "viability")
  ## iterate through each column of dosage and viability  
	for (i in 1:observations) {
        message(paste("Processing Observation", i))
        ## get corresponding dose and viability for that observation
        doseVec <- data[,colBeforeFirstDose + i]
        viabilityVec <- data[, colBeforeFirstViability + i]
        ## create frame containing dose-respone data for one dosage
	    newFrame <- cbind(cells, drugs, study_id, doseVec, viabilityVec)
        newFrame <- setNames(newFrame, sensitivityDataNames)
        ## append current frame to data processed so far
        doseResponseFrame <- rbind(doseResponseFrame, newFrame)
	}
  ## return frame
  rownames(doseResponseFrame) <- NULL
  return(doseResponseFrame)
}

## Convert read CSV file of summary statistics
## into format which can be inputted into PharmacoDB as a data frame
## data: read .csv file with the original data
## study_id: integer which identifies the current study
## stat_id: integer vector which identifies summary statistics
## stat_name: names of the summary statistics being considered
## TODO: deal with issues which may arise as a result of loss of precision
processSummaryData <- function(data, study_id, stat_id, stat_name) {
  ## create a empty data frame for processed data 
  summaryFrame <- data.frame()
  ## get list of cell lines and drugs used in each study 
  cells <- data["cellline"]
  drugs <- data["drug"]
  ## iterate through each statistic
  if (length(stat_id) != length(stat_name)) {
  	stop("Number of ids do not correspond to number of statistics!")
  }
  summaryDataNames <- c("cellline_id", "drug_id", "study_id", "stat_id", "stat_value")
  for (i in 1:length(stat_id)) {
        message(paste("Processing", stat_name[i], "data"))
        curStat <- data[stat_name[i]]
        ## create frame containing data for current statistic
        curStatFrame <- cbind(cells, drugs, study_id, stat_id[i], curStat)
		curStatFrame <- setNames(curStatFrame, summaryDataNames)
	    ## append current frame to data processed so far
		summaryFrame <- rbind(summaryFrame, curStatFrame)
    }
   ## return frame
   rownames(summaryFrame) <- NULL
   return(summaryFrame)
 }

