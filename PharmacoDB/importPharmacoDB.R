library(RMySQL)
source("importFunctions.R")
source("convertCGPSensitivity.R")

## Construct STUDIES Table
message("Writing to STUDIES table")
studies <- c("CCLE", "CGP", "GRAY", "GNE")

## import studies vec into STUDIES table
sapply(studies, doStudiesQuery)

## Construct METADATA table
message("Writing to METADATA table")
metadata <- c("tissue")

## import metadata vec into METADATA table
sapply(metadata, doMetadataQuery)

## Construct cell line curation and cell line metadata tables
message("Writing to CELL_LINE_CURATION and CELL_LINE_METADATA table")

## Read data from "cell_line_annotation_all_ATTENTION.csv"
if (!file.exists("cell_annotation_all.csv")) {
	stop("Cell line curations do not exist.")
}
cells <- read.csv("cell_annotation_all.csv", na.strings = "")

## write to the cell tables
sapply(studies, function (study) writeToCellTables (cells, study))
sapply(studies, function (study) writeToCellTissues (cells, study))

## Construct DRUG_CURATION TABLE
message("Writing to DRUG_CURATION table")

## Import data into DRUG_CURAITON table
sapply(studies[1:3], writeToDrugTables) # no drug info for GNE

## write to TISSUE tables
message("Writing to TISSUE_TYPE")

## Read data and import
tissues <- read.csv("matching_tissue.csv", na.strings = "")

## Import data into TISSUE_TYPE table
sapply(studies, function (study) writeToTissueTables (tissues, study))

## Construct SENSITIVITY_DATA TABLE
message("Writing to SENSITIVITY_DATA")

## Import data into table
sapply(studies, writeToSensitivityTables)

## Construct SUMMARY_STATISTICS
message("Writing to SUMMARY_STATISTICS")

## Read all summary statistics names from sample file "ccLe_sensitivity.csv"
## Remove row of row names
stats <- colnames(read.csv("ccle_sensitivity.csv")) ## Use ccle since it contains
## all other summary statistics
stats <- stats[4:length(stats)] 

## import each element of studies into the table
sapply(stats, doStatQuery)

## Construct SUMMARY_STATISTICS_DATA
message("Writing to SUMMARY_STATISTICS_DATA")

## Import data into table
sapply(studies, writeToSummaryTables)

## remove workspace objs after import
rm(list = ls())