library(RMySQL)
library(rcdk)
source("chemicalTransUtils.R")

## Script to initialize application objects by
## getting drug curation, drug fingerprints, and
## chembl ids

## Get table for drug similarity application
con <- dbConnect(MySQL(), dbname = "PharmacoDB", hostname = "localhost", username = "root", password = "")
drugs <- dbGetQuery(con, "SELECT * FROM  drug_similarity_app")
dbDisconnect(con)

## Remove na values for smiles to get fingerprints
drugs <- drugs[!is.na(drugs$smiles), ]
target.mols <- parse.smiles(as.character(drugs[,"smiles"]))
target.fps <- lapply(target.mols, get.fingerprint, type="extended")
drugs$fps <- target.fps
save(drugs, file = "drug_dat_with_fps.RData")

## Get CHEMBL IDS
chembls <- lapply(unique(drugs$drug_id), getChEMBLfromName)
names(chembls) <- unique(drugs$drug_id)
chembls <- Filter(function (x) (length(x) > 0 && !is.na(x)), chembls) ## Filter ChEMBL ids which were not found
save(chembls, file = "chembls.RData")

