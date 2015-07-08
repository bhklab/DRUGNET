library(RMySQL)
source("importFunctions.R")

## Construct STUDIES Table
message("Writing to STUDIES table")
studies <- c("CCLE", "CGP", "GRAY", "GNE", "GSK", "CTRP", "NCI60")

## import studies vec into STUDIES table
#sapply(studies, doStudiesQuery)


## Construct cell line curation and cell line metadata tables
message("Writing to CELL_LINE_CURATION table")

## Read data from "cell_line_annotation_all_ATTENTION.csv"
setwd("curation")
if (!file.exists("cell_annotation_all.csv")) {
	stop("Cell line curations do not exist.")
}
#cells <- read.csv("cell_annotation_all.csv", na.strings = "")

## write to the cell tables
tissue_match <- read.csv("matching_tissue.csv", stringsAsFactor = FALSE, na.strings = "")
#lapply(studies, function (study) processCell (cells, study))
#lapply(studies, function (study) processCellTissues (cells, study))

#source("importHistologies.R")

## Construct DRUG_CURATION TABLE
message("Writing to DRUG_CURATION table")

drugs <- read.csv("drugs_with_ids.csv")
## Import data into DRUG_CURAITON table
## TODO: Process duplicates here
#processDrugsMetadata(drugs)
#sapply(studies, function (study) processDrugs(drugs, study))

##  Run script to populate gene expression metadata.
setwd("~/Desktop/pharmacodb_dat/gene expression")
#setwd("/mnt/work1/users/home2/bhcoop4/pharmacodb_dat")
source("processGenes.R")

## Construct SENSITIVITY_DATA TABLE
message("Writing to SENSITIVITY_DATA")
#setwd("/mnt/work1/users/home2/bhcoop4/pharmacodb_dat")
## Import data into table
setwd("~/Desktop/pharmacodb_dat/sensitivity")
#setwd("~/mnt/work1/users/bhklab/users/adrian/sensitivity")
#source("convertCGPSensitivity.R")
studies_with_dat <- c("ccle", "cgp", "gne", "gray")
files <- sapply(studies_with_dat, function (x) paste0(x, "_sensitivity_detail.csv"))
files_summary <-  sapply(studies_with_dat, function (x) paste0(x, "_sensitivity.csv"))
mapply(processDoseResponseData, lapply(files, read.csv), studies_with_dat)
processDoseResponseData(read.delim("ctrp_sensitivity_detail.txt"), "ctrp")

## Construct SUMMARY_STATISTICS
message("Writing to SUMMARY_STATISTICS")

## Import all summary statistics into table
stats <- c("applied.concentration.no", "AUC", "slope", "slope0", "AUC_Published", "AUC_Published_raw", "IC50_Published")
sapply(stats, doStatQuery)

## Construct SUMMARY_STATISTICS_DATA
message("Writing to SUMMARY_STATISTICS_DATA")

## Import data into table
mapply(processSummaryData, lapply(files_summary, read.csv), studies_with_dat)
processSummaryData(read.delim("ctrp_sensitivity.txt"), "ctrp")

## remove workspace objs after import and disconnect any remaining connections
lapply(dbListConnections(MySQL()), dbDisconnect)
rm(list = ls())
