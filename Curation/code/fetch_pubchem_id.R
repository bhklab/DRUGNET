####################################
path.out <- file.path ("output")
drugall <- read.csv("drug_annotation_all.csv")

library(devtools) # install.packages("devtools") if missing
library(jsonlite) # install.packages("jsonlite")
install_github(repo = "CTSgetR", username = "dgrapov")
library(CTSgetR)

id<- drugall[,1]
from<-"Chemical Name"
to<- "PubChem CID"
allcids <- multi.CTSgetR(id,from,to,progress=TRUE)

save(allcids, file = file.path(path.out, "pubchem_id.RData"))
################################


path.out <- "output"
################## Code for NCI60 compounds, around 8000 ###########

require(devtools) || install.packages("devtools")


######### install CTSgetR
require(jsonlite) || install.packages("jsonlite")
require(CTSgetR) || devtools::install_github(repo = "CTSgetR", username = "dgrapov")

#### Install webchem package
require(webchem) || devtools::install_github("edild/webchem")

############## Convert NSC's to PubChem CID

######## just a compound or vector of NSCs
drug.nci <- read.csv(file.path(path.out, "drug_annotation_NCI60.csv"), stringsAsFactor=F, na.strings = c("", " ", "NA"))
drug.nci2 <- drug.nci[!is.na(drug.nci[,2]),]
id <- drug.nci2[,1]
from <- "DTP/NCI"
to <- c("PubChem CID")
results <- multi.CTSgetR(id,from,to,progress=TRUE)
save(results, file="nsc_to_pubchem_nci60.RData")

###### Use this function to convert Pubchem CIDs from NCI60 to smiles and inchikeys

results.nci <- list()
for(i in 1:nrow(results)){
  testn <- cid_compinfo(results[i,2], verbose = TRUE)
  testn <- c(testn$CID,testn$CanonicalSmiles, testn$InChIKey )
  results.nci[[i]] <- testn
}

save(results.nci, file="smiles_inchikey_nci60.RData")

results.nci <- lapply(results.nci, function(list.element) {
  if(length(list.element) > 0)
    return(list.element)
  else
    return(NA)
})

df.nci <- do.call(rbind, results.nci)

df.nci.identifiers <- data.frame("nsc"=drug.nci2[,1], "name"=drug.nci2[,2], "cid"=df.nci[,1],"smiles"=df.nci[,2], "inchikey"=df.nci[,3])

#### Remove NAs and duplicates to parse smiles 
df.nci.identifiers <- na.omit(df.nci.identifiers)
df.nci.identifiers <- df.nci.identifiers[!duplicated(df.nci.identifiers[,"name"]),]

#### Case study: similarity search using the rcdk package
require(rcdk)
library(rcdk)
#query.mol <- parse.smiles("CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)OC(C)(C)C)O)O)OC(=O)C6=CC=CC=C6)(CO4)OC(=O)C)O)C)O")[[1]] # docetaxel
#query.mol <- parse.smiles("CC(C)C[C@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)c2cnccn2)B(O)O")[[1]] # bortezomib
#query.mol <- parse.smiles("COCCOC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC=CC(=C3)C#C)OCCOC")[[1]] # erlotinib
#query.mol <- parse.smiles("CC12C(CC(O1)N3C4=CC=CC=C4C5=C6C(=C7C8=CC=CC=C8N2C7=C53)CNC6=O)(CO)O")[[1]] # lestaurtinib
#query.mol <- parse.smiles("CCN(CC)C(=S)SSC(=S)N(CC)CC")[[1]] # disulfiram
#query.mol <- parse.smiles("COC1=C(C=C2C(=C1)N=CN=C2NC3=CC(=C(C=C3)F)Cl)OCCCN4CCOCC4")[[1]] # gefitinib
#query.mol <- parse.smiles("[NH2-].[NH2-].Cl[Pt+2]Cl")[[1]]

target.mols <- parse.smiles(as.character(df.nci.identifiers[,"smiles"]))

query.fp <- get.fingerprint(query.mol, type="extended")

target.fps <- lapply(target.mols, get.fingerprint, type="extended")

sims <- unlist(lapply(target.fps, distance,fp2=query.fp, method="tanimoto"))

#### Get similar compounds (you can use any threshold, I like 0.3-1 for biological significance!)
sims2 <- data.frame(sims)
rownames(sims2) <- df.nci.identifiers[,"name"]
sims2 <- sims2[sims2[,1] >= 0.3,,drop=F]
sims2 <- sims2[order(sims2[,1],decreasing=TRUE),,drop=F]

