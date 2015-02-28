## path to files
path.data <- file.path("data", "NCI60")
path.drug <- file.path(path.data, "drug")
path.cell <- file.path(path.data, "celline")
## create directories
if(!file.exists(path.data)) { dir.create(path.data, showWarnings=FALSE, recursive=TRUE) }
if(!file.exists(path.drug)) { dir.create(path.drug, showWarnings=FALSE, recursive=TRUE) }
if(!file.exists(path.cell)) { dir.create(path.cell, showWarnings=FALSE, recursive=TRUE) }

## drug pheno
message("Download drug sensitivity measurements")
myfn <- file.path(path.drug, "tmp", "nci60_drug_pheno_file.txt")
if (!file.exists(myfn)) {
  if(!file.exists(file.path(path.drug, "tmp"))) { dir.create(file.path(path.drug, "tmp"), showWarnings=FALSE, recursive=TRUE) }
  ## drug sensitivity data from the addendum in Nature
  dwl.status <- download.file(url="http://discover.nci.nih.gov/cellminerdata/rawdata/DTP_NCI60_RAW.txt.zip", destfile=file.path(path.drug, "tmp", "DTP_NCI60_RAW.txt.zip"), quiet=TRUE)
  if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
  fff <- unzip(zipfile=file.path(path.drug, "tmp", "DTP_NCI60_RAW.txt.zip"), list=TRUE)
  file.copy(from=file.path(path.drug, "tmp", "DTP_NCI60_RAW_1_5.txt"), to=myfn)
}

## read drug pheno
drugpheno <- read.csv(file.path(path.drug, "tmp", "nci60_drug_pheno_file.txt"), sep="\t")
drugpheno[!is.na(drugpheno) & drugpheno == ""] <- NA

sampleinfo <- data.frame("cellid"=colnames(drugpheno)[-c(1:6)], "tissueid"=sapply(strsplit(colnames(drugpheno)[-c(1:6)], split="[.]"), function (x) { return(x[[1]]) }))
sampleinfo[ , "cellid"] <- apply(sampleinfo, 1, function (x) {
  x[1] <- gsub(sprintf("^%s.", x[2]), sprintf("%s:", x[2]), x[1])
})
druginfo <- drugpheno[ , c(1:4)]
druginfo <- data.frame("drugid"=as.character(druginfo[ , "NSC.."]), "drug.name"=as.character(druginfo[ , "Drug.Name"]), druginfo)
druginfo <- druginfo[!duplicated(druginfo[ , "drugid"]), , drop=FALSE]
myx <- match(c("drugid", "drug.name"), colnames(druginfo))
druginfo <- druginfo[ , c(myx, setdiff(1:ncol(druginfo), myx)), drop=FALSE]

write.csv(sampleinfo, file=file.path(path.out, "cell_line_annotation_NCI60.csv"), row.names=FALSE)
write.csv(druginfo, file=file.path(path.out, "drug_annotation_NCI60.csv"), row.names=FALSE)