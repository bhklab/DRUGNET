drugs <- read.csv("drug_annotation_all.csv", na.strings = "", stringsAsFactors = FALSE)
frame <- read.csv("all_ids.csv")

frame <- frame[!duplicated(frame$drug_id),]
drug_names <-  as.character(frame$drug_id)
smiles <- as.character(frame$smiles)
inchikey <- as.character(frame$inchikey)
cid <- as.character(frame$cid)
names(smiles) <- drug_names
names(inchikey) <- drug_names
names(cid) <- drug_names

appendCol <- function (entry, id) {
  entry <- unlist(entry)
  return(id[entry["unique.drugid"]])
}

print("Getting SMILES")
drugs$smiles <- apply(drugs, 1, function (x) appendCol(x, smiles))
print("Getting InChIKey")
drugs$inchikey <- apply(drugs, 1, function (x) appendCol(x, inchikey))
print("Getting PubChem CID")
drugs$cid <- apply(drugs, 1, function (x) appendCol(x, cid))

postProcess <- function(x) {
  unlist(lapply(x, function (x) ifelse(is.null(x), NA, x)))
}

drugs$smiles <- postProcess(drugs$smiles)
drugs$inchikey <- postProcess(drugs$inchikey)
drugs$cid <- postProcess(drugs$cid)
write.csv(drugs, file = "drugs_with_ids.csv")

