library(webchem)
library(RJSONIO)
library(RCurl)

## Global variables
chemblURL <- "https://www.ebi.ac.uk/chembl/api/data/molecule/"

## Convert chemical name into a SMILES string based on record in ChEMBL database
## string: a chemical name
getSMILESfromName <- function(string) {
  chembl <- getChEMBLfromName(string)
  url <- paste0(chemblURL, chembl, ".json")
  vals <- getURL(url)
  parsed <- fromJSON(vals)
  return(parsed[["molecule_structures"]]["canonical_smiles"])
}

## Convert chemical name into a ChEMBL identifier based on CTS
## string: a chemical name
getChEMBLfromName <- function(string) {
  tryCatch(cts_convert(string, 'Chemical Name', 'CHEMBL', first = TRUE), error = function (e) NA)
}



