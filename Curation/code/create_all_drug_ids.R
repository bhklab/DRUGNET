## Load original annotation without studies
load("drugs_with_smiles.RData")
frame  <- frame[frame$study_name != "NCI60", -6]
frame <- frame[!duplicated(frame$drug_id),-2]

## Load NCI60 drugs with the substituted names.
nci60 <- read.csv("drug_annotation_NCI60_new_clean.csv", stringsAsFactors = FALSE)
nci60 <- nci60[!is.na(nci60$drug.name),]
nci60 <- nci60[,3:6]

## Load CMAP drugs
load("cmap_with_cbid_smiles.RData")
cmap <- as.data.frame(cmap)
cmap <- cmap[,c(1,2,3,6)] ## Get name, inchikey, pubchem cid, cbid_smiles

## Bring together all curation and write a csv for it
colnames(nci60) <- c("drug_id", "inchikey", "cid", "smiles")
colnames(cmap) <-c("drug_id", "inchikey", "cid", "smiles")
all <- rbind(frame, nci60, cmap)
write.csv(all, file = "all_ids.csv")

