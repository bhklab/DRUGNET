library(RMySQL)
library(rcdk)

## Get table for drug similarity application
con <- dbConnect(MySQL(), dbname = "PharmacoDB", hostname = "localhost", username = "root", password = "")
frame <- dbGetQuery(con, "SELECT * FROM drug_similiarity_app_view")
dbDisconnect(con)

## Remove na values for smiles to get fingerprints
frame <- frame[!is.na(frame$smiles), ]
target.mols <- parse.smiles(as.character(frame[,"smiles"]))
target.fps <- lapply(target.mols, get.fingerprint, type="extended")

save(frame, file = "data_obj.RData")
save(target.fps, file = "fingerprints.RData")