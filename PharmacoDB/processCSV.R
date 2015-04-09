## Functions which process existing CSV files 
## which hold data for PharmacoDB

## Author: Adrian She
## Last Updated: April 1, 2015

library(webchem)

## Extract cell line names and synonyms for a particular study as a data frame
## data: read csv file for "cell_line_annotation_all_ATTENTION.csv"
## study: name of the study
## study_id: integer identifying the study name
processCell <- function(data, study, study_id) { 
  message(paste("Processing cell names for", study))
  study <- toupper(study) # cell line columns are upper case
  uniqueIDs <- data["unique.cellid"] ## get unique names
  cells <- data[paste0(study, ".cellid")] ## get names in that study
  frame <- cbind(uniqueIDs, study_id, cells) # construct frame of cell line info
  frame <- frame[complete.cases(frame),] # remove NAs
  cellCurationNames <- c("cellline_id", "study_id", "cellline_name") 
  colnames(frame) <- cellCurationNames # rename according to db colnames
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
  tissues <- data[paste0(study, ".tissueid")]
  frame <- cbind(uniqueIDs, metadata_id, study_id, tissues) # construct frames for cell line info
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
processDrugsByStudy <- function(drugs, study, study_id) {
	if (study == "CGP") { ## remove problematic drug from CGP frame
		drugs <- drugs[-40,] 
	}
	
  drugs <- cbind(drugs["drugid"], drugs["drug.name"])
  drugs <- drugs[complete.cases(drugs),] # get drugs only in current study
  cidVec <- convertToCID(drugs)

  ## Construct frame
  drugFrame <- data.frame(cbind(drugs["drugid"],  study_id, drugs["drug.name"], cidVec))
	# rename according to colnames of table in the schema
  drugCurationNames <- c("drug_id", "study_id", "drug_name", "pubchem_cid")
	colnames(drugFrame) <- drugCurationNames
	rownames(drugFrame) <- NULL
	return(drugFrame)
}

## Extract drugs used in a particular study  as a data frame
## drugs: read CSV file for "drug_annotation_STUDYNAME.csv"
## study_id: integer identifying the study
processDrugs <- function(drugs, study, study_id) {
    study <- toupper(study) ## study has to be upper case
	drugs <- drugs[-40,] ## remove problematic drug 681640 from data
    drugs <- cbind(drugs["unique.drugid"], drugs[paste0(study, ".drugid")]) ## Subset required rows
    drugs <- drugs[complete.cases(drugs),] # get drugs only in current study
    cidVec <- convertToCID(drugs[,2])
    
    ## Construct frame
    drugFrame <- data.frame(cbind(drugs[,1],  study_id, drugs[,2], cidVec))
	# rename according to colnames of table in the schema
    drugCurationNames <- c("drug_id", "study_id", "drug_name", "pubchem_cid")
	colnames(drugFrame) <- drugCurationNames
	rownames(drugFrame) <- NULL
	return(drugFrame)
}

## Convert a vector of chemical names to PubChem CID
convertToCID <- function(drugs) {
  ## convert to current drugs to vector
  drugsVec <- as.vector(t(drugs))
  ## remove white spaces from all drug names, otherwise webchem will crash
  drugsVec <- gsub("\\s+", "", drugsVec)
  # get necessary PubChem IDs
  cidVec <- c()
  for (i in 1:length(drugsVec)) {
    cidVec <- c(cidVec, 
      as.integer(cts_convert(query = drugsVec[i], from = "Chemical Name", to = "PubChem CID"))[1])
  }
  return(cidVec)
}

## Process the file matching_tissue.csv for import into the database 
## data: Read csv file for matching_tissue.csv
## study: Name of study
## study_id: Integer which identifies the study
processTissue <- function(data, study, study_id) { 
  message(paste("Processing tissue names for", study))
  study <- toupper(study) # cell line columns are upper case
  uniqueIDs <- data["COSMIC.tissueid"]
  tissues <- data[paste0(study, ".tissueid")]
  frame <- cbind(uniqueIDs, study_id, tissues) # construct frame of cell line info
  frame <- frame[complete.cases(frame),] # remove NAs
  tissueCurationNames <- c("unique_tissue_name", "study_id", "tissue_name") 
  colnames(frame) <- tissueCurationNames
  # rename according to colnames of table in the schema
  rownames(frame) <- NULL  # remove row.names caused by subsetting in complete.cases
  return(frame)
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

  expid <- 1
  sensitivityDataNames <- c("cellline_id", "drug_id", "study_id", "log10dose", "exp_id", "viability")
  ## iterate through each column of dosage and viability  
	for (i in 1:observations) {
        message(paste("Processing Observation", i))
        ## get corresponding dose and viability for that observation
        doseVec <- data[,colBeforeFirstDose + i]
        viabilityVec <- data[, colBeforeFirstViability + i]
        ## create frame containing dose-respone data for one dosage
	      newFrame <- cbind(cells, drugs, study_id, doseVec, expid, viabilityVec)
        newFrame <- setNames(newFrame, sensitivityDataNames)
        ## append current frame to data processed so far
        doseResponseFrame <- rbind(doseResponseFrame, newFrame)
	}
  ## return frame
  rownames(doseResponseFrame) <- NULL
  return(doseResponseFrame)
}


## When there are experimental replicates, need to process by row
processDoseResponseByRow <- function(data, study_id) {
  doseResponseFrame <- data.frame()
  
  ## Get all unique experiments
  experiments <- unique(cbind(data$cellline, data$drug))

  ## Get number of observations done in the experiment
  colBeforeFirstDose <- grep("log10_dose_1", colnames(data)) - 1
  colBeforeFirstViability <- grep("viability_1", colnames(data)) - 1
  observations <- colBeforeFirstViability - colBeforeFirstDose 
  sensitivityDataNames <- c("cellline_id", "drug_id", "study_id", "log10dose", "exp_id", "viability")

  ## Process by experiment
  for (i in 1:nrow(experiments)) {
    curFrame <- data.frame()
    ## Get all data
    ## corresponding to current experiment
    cell <- experiments[i, 1]
    drug <- experiments[i, 2]
    exps <- data[which(data$cellline == cell & data$drug == drug),] 

    ## Iterate over no. of observations
    for (j in 1:observations) {
        doseFrame <- cbind(exps["cellline"], exps["drug"], study_id, 
          exps[,colBeforeFirstDose + j], 1:nrow(exps), exps[, colBeforeFirstViability + j])
        doseFrame <- setNames(doseFrame, sensitivityDataNames)
        curFrame <- rbind(curFrame, doseFrame) ## bind current observation to current frame
    }
    doseResponseFrame <- rbind(doseResponseFrame, curFrame) ## bind current experiment
  }

  ## Return the frame
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
  expid <- 1
  summaryDataNames <- c("cellline_id", "drug_id", "study_id", "stat_id", "exp_id", "stat_value")

  ## Iterate through each statistic
  for (i in 1:length(stat_id)) {
        message(paste("Processing", stat_name[i], "data"))
        curStat <- data[stat_name[i]]
        ## create frame containing data for current statistic
        curStatFrame <- cbind(cells, drugs, study_id, stat_id[i], expid, curStat)
		    curStatFrame <- setNames(curStatFrame, summaryDataNames)
	      ## append current frame to data processed so far
		    summaryFrame <- rbind(summaryFrame, curStatFrame)
    }
   ## return frame
   rownames(summaryFrame) <- NULL
   return(summaryFrame)
 }

 ## When there are experimental replicates, need to process by row
processSummaryDataByRow <- function(data, study_id, stat_id, stat_name) {
  summaryFrame <- data.frame()
  
  ## Get all unique experiments
  experiments <- unique(cbind(data$cellline, data$drug))

  ## Process by experiment
  summaryDataNames <- c("cellline_id", "drug_id", "study_id", "stat_id", "exp_id", "stat_value")
  for (i in 1:nrow(experiments)) {
    statFrame <- data.frame()

    ## Get all data corresponding to current experiment
    cell <- experiments[i, 1]
    drug <- experiments[i, 2]
    exps <- data[which(data$cellline == cell & data$drug == drug),]

    ## Iterate over no. of summary statistics
    for (j in 1:length(stat_id)) {
        curStat <- exps[stat_name[j]] ## get current statistic
        ## create frame containing data for current statistic
        curFrame <- cbind(exps["cellline"], exps["drug"], study_id, stat_id[j], 1:nrow(exps), curStat)
        curFrame <- setNames(curFrame, summaryDataNames)
        ## append current frame to data processed so far
        statFrame <- rbind(statFrame, curFrame)
    }
    summaryFrame <- rbind(summaryFrame, statFrame)
  }
  ## Return the frame
  rownames(summaryFrame) <- NULL
  return(summaryFrame)
}



