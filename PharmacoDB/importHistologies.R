frame <- data.frame()

matching <- read.csv("diagnosis_type_match.csv", stringsAsFactors = F)

## vec: vector of histology names
## study: study
matchHistology <- function(hists, study) {
    if (paste0("type.", study) %in% colnames(matching)) {
      col <- paste0("type.", study)
      histology <- as.vector(t(hists))
      vals <- sapply(histology, function (x) {
          match <- paste(matching[which(toupper(matching[, col]) == toupper(x)), "unique.cancer.type"], collapse = '/')
          if (match == "") { return(NA)}
          return(match)
        })
      
      return(vals)
    }
}

## name: name of study
## colname: Name of column in CSV which stores histology information
processHistology <- function(name, colname) {
  dat <- read.csv(paste0("cell_line_annotation_", toupper(name), ".csv"), stringsAsFactors = FALSE)
  frame <- data.frame("cellline_id" = dat$cellid,
                      "study_id" = getStudyID(name),
                      "histology_name" = dat[, colname])
  frame <- frame[complete.cases(frame),]
  frame$unique_histology_type <- matchHistology(frame$histology_name, tolower(name))
  frame <- as.data.frame(frame, stringsAsFactors = FALSE)
  message(paste0("Processing histologies for", name))
  return(frame)
}

histology_studies <- c("CCLE", "GRAY", "GNE", "GSK", "CTRP")
colname <- c("Hist.Subtype1", "Transcriptional.subtype", "Tissue_Diagnosis", "Characteristics.DiseaseState.", "cell_line_sublineage")
histologies <- mapply(processHistology, histology_studies, colname, SIMPLIFY = FALSE)

## GRAY is solely a breast cancer data set
message("Matching histologies for GRAY, CTRP, GSK")
histologies[["GRAY"]]$unique_histology_type <- histologies[["GRAY"]]$histology_name
histologies[["CTRP"]]$unique_histology_type <- tolower(histologies[["CTRP"]]$histology_name)
histologies[["GSK"]]$unique_histology_type <- sapply(histologies[["GSK"]]$histology_name, function (x) gsub(" ", "_", tolower(x)))
## Quick fix for GSK until matching_tissue_updated

## Process CGP here. 
library(openxlsx)
cgp <- read.xlsx("Supplementary_data_final_Apr6.xlsx", 2)
cgp_info <- data.frame("cellline_id" = cgp$Cell.Line, "study_id" = getStudyID("CGP"), "histology_name" = cgp[,"Cancer-type"])
cgp_info$unique_histology_type <- matchHistology(cgp_info$histology_name, "cgp")
histologies[["CGP"]] <- cgp_info
message("Processed CGP")

# Write final frame to DB
toImport <- do.call("rbind", histologies)
toImport$id <- 1:nrow(toImport)
connection <- connDB()
dbWriteTable(connection, "HISTOLOGY_CURATION", toImport, append = TRUE, row.names = FALSE)
dbDisconnect(connection)
