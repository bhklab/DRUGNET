
################## Code for NCI60 compounds, around 8000 ###########

######### install CTSgetR
library(devtools) # install.packages("devtools") if missing
library(jsonlite) # install.packages("jsonlite")
library(CTSgetR)

#### Install webchem package

#install.packages("devtools")
library("devtools")
#install_github("edild/webchem")
library(webchem)
library(gdata)
library(pbapply)

######## just a compound or vector of NSCs
drug.nci <- read.csv(file.path(path.out, "drug_annotation_NCI60.csv"), stringsAsFactor=F)
id <- drug.nci$NSC..

####### Get compound name from cid, InChIKey and PubChem CID
message("Getting Names from NSC")
from <- "DTP/NCI"
to <- c("Chemical Name", "InChIKey", "PubChem CID")
for (i in seq(1,length(id), 1000)) {
    message(paste("Processing batch", i))
    toTranslate <- id[i:(i+999)]
    toTranslate <- toTranslate[!is.na(toTranslate)]
    results <- pbapply::pblapply(toTranslate, function (x) tryCatch(multi.CTSgetR(x,from,to, progress=FALSE), error = function (e) { NA }))
    save(results, file = paste0("batch_", i, ".RData"))
}

## Rbind each batch into one data frame
nci60_complete_withnas <- data.frame()
for (i in seq(1, 49001, 1000)) {
    load(paste0("batch_", i, ".RData"))
	nas <- grep("NA", results)
    results <- results[!is.na(results)]
	results <- data.table::rbindlist(results)
	nci60_complete_withnas <- rbind(nci60_complete_withnas, results)
}
nci60_complete_withnas <- unknownToNA(nci60_complete_withnas, "NA")
save(nci60_complete_withnas, file = "NCI60.RData")


## Get SMILES from data frame
getAllSMILES <- function(frame) {
	cids <- as.character(frame$"PubChem CID")
	smiles <- pblapply(cids,
    function (x) {
        if (!is.na(x)) {
            cid_compinfo(x)$CanonicalSmiles
        }
        else {
            NA
        }
    })
	return(smiles)
}
message("Getting SMILES")
for (i in seq(1, nrow(nci60_complete_withnas), 1000)) {
    if ((i + 999) > nrow(nci60_complete_withnas)) {
        smiles <- getAllSMILES(nci60_complete_withnas[i:nrow(nci60_complete_withnas),])
    }
    else {
        smiles <- getAllSMILES(nci60_complete_withnas[i:(i+999),])
    }
    save(smiles, file = paste0("smiles_batch", i, ".RData"))
}
    
all_smiles <- list()
for (i in seq(1, 49001, 1000)) {
    load(paste0("smiles_batch", i, ".RData"))
    all_smiles <- c(all_smiles, smiles)
}
    
nci60_complete_withnas$smiles <- as.character(all_smiles)
nci60_complete_withnas$smiles <- unknownToNA(nci60_complete_withnas$smiles, "NA")
nci60_complete_withnas$smiles <- as.factor(nci60_complete_withnas$smiles)
nci60_with_smiles <- nci60_complete_withnas
write.csv(nci60_with_smiles, file.path(path.out, "drug_annotation_NCI60.csv"))
