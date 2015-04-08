## Original CGP sensitivity detail includes drug.ids only
## need to change the names for import into PharmacoDB
cgpSensitivity <- read.csv("cgp_sensitivity_detail.csv") ## CGP Sensitivity file
cgpSumm <- read.csv("cgp_sensitivity.csv") ## Summary statistics
cgpDrugDetail <- read.csv("drug_annotation_CGP.csv", stringsAsFactors = FALSE) ## CGP Drug Annotation File

cells <- cgpSensitivity["cellline"]
odrug <- as.vector(t(cgpSensitivity["drug"]))
vals <- cgpSensitivity[,4:ncol(cgpSensitivity)]

## Convert each drug.id into a drug Name. This list is the same
## in both sensitivity/summary because each experiment has a
## sumary value
idNameCorrespondence <- cbind(cgpDrugDetail["drug.id"], cgpDrugDetail["Name"])
drug <- c()
for (i in 1:length(odrug)) {
  drug[i] <- idNameCorrespondence[which(idNameCorrespondence$drug.id ==  odrug[i]), "Name"]
}

write.csv(cbind(cells, drug, vals), file = "cgp_sensitivity_detail.csv")
write.csv(cbind(cgpSumm["cellline"], drug, cgpSumm[,4:ncol(cgpSumm)]),
file = "cgp_sensitivity.csv")