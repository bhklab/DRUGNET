frame <- data.frame()

matching <- read.csv("diagnosis_type_match.csv", stringsAsFactors = F, na.strings = "")

matchHistology <- function(frame, study) {
    if (paste0("type.", study) %in% colnames(matching)) {
      col <- paste0("type.", study)
      histology <- as.vector(t(frame$histology_name))
      frame$unique_histology_type <- 
        sapply(histology, function (x) {
          paste(matching[which(toupper(matching[, col]) == toupper(x)), "unique.cancer.type"], collapse = '/')
        })
    }
    else {
      frame$unique_histology_type <- NA
    }
    return(frame)
}

## name: name of study
## colname: Name of column in CSV which stores histology information
processHistology <- function(name, colname) {
  dat <- read.csv(paste0("cell_line_annotation_", toupper(name), ".csv"), stringsAsFactors = FALSE)
  frame <- data.frame("cellline_id" = dat$cellid,
                      "study_id" = getStudyID(name),
                      "histology_name" = dat[, colname])
  frame <- frame[complete.cases(frame),]
  frame <- matchHistology(frame, tolower(name))
  frame <- as.data.frame(frame)
  message(paste0("Processing histologies for", name))
  return(frame)
}

histology_studies <- c("CCLE", "GRAY", "GNE", "GSK", "CTRP")
colname <- c("Hist.Subtype1", "Transcriptional.subtype", "Tissue_Diagnosis", "Characteristics.DiseaseState.", "cell_line_sublineage")
histologies <- mapply(processHistology, histology_studies, colname, SIMPLIFY = FALSE)


## Process CGP here. 
library(openxlsx)
cgp <- read.xlsx("Supplementary_data_final_Apr6.xlsx", 2)
cgp_info <- data.frame("cellline_id" = cgp$Cell.Line, "study_id" = getStudyID("CGP"), "histology_name" = cgp[,"Cancer-type"])
cgp_info <- matchHistology(cgp_info, "cgp")
histologies[["CGP"]] <- cgp_info
message("Processed CGP")

# Write final frame to DB
toImport <- do.call("rbind", histologies)
toImport$id <- 1:nrow(toImport)
connection <- connDB()
dbWriteTable(connection, "HISTOLOGY_CURATION", toImport, append = TRUE, row.names = FALSE)
dbDisconnect(connection)