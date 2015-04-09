## path to files
path.data <- file.path("data", "GNE")
path.drug <- file.path(path.data, "drug")
path.cell <- file.path(path.data, "celline")
## create directories
if(!file.exists(path.data)) { dir.create(path.data, showWarnings=FALSE, recursive=TRUE) }
if(!file.exists(path.drug)) { dir.create(path.drug, showWarnings=FALSE, recursive=TRUE) }
if(!file.exists(path.cell)) { dir.create(path.cell, showWarnings=FALSE, recursive=TRUE) }


## sample information
message("Download sample information")
myfn <- file.path(path.cell, "tmp", "gne_ge_sampleinfo.map")
if (!file.exists(myfn)) {
	if(!file.exists(file.path(path.cell, "tmp"))) { dir.create(file.path(path.cell, "tmp"), showWarnings=FALSE, recursive=TRUE) }
  	message(sprintf("Please download Metadata for series EGAD00001000725 from https://www.ebi.ac.uk/ega/datasets/EGAD00001000725 and copy %s to %s", "Run_Sample_meta_info.map", dirname(myfn)))
  file.copy(from=file.path(path.cell, "tmp", "Run_Sample_meta_info.map"), to=myfn)
}

## drug sensitivity data is not available online and recieved by email from Chris Klijn
myfn <- file.path(path.drug, "tmp", "KlijnEtAl_drugDosingMeanViab.txt")
druginfo <- read.table(myfn, stringsAsFactors = FALSE, sep = "\t")

colnames(druginfo) <- druginfo[1,]
druginfo <- druginfo[-1,]
druginfo$"Cell line" <- gsub(",", "_", druginfo$"Cell line")



## Additional drug sensitivity just for one drug which is downloaded from internet
 message("Download drug sensitivity measurements")
 myfn <- file.path(path.drug, "tmp", "gne_drug_pheno.xls")
 if (!file.exists(myfn)) {
   if(!file.exists(file.path(path.drug, "tmp"))) { dir.create(file.path(path.drug, "tmp"), showWarnings=FALSE, recursive=TRUE) }
   dwl.status <- download.file(url="http://www.nature.com/nbt/journal/v33/n3/extref/nbt.3080-S9.xls", destfile=file.path(path.drug, "tmp", "nbt.3080-S9.xls"), quiet=TRUE)
   if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
   file.copy(from=file.path(path.drug, "tmp", "nbt.3080-S9.xls"), to=myfn)
 }

require(gdata) || stop("Library gdata is not available!")
temp <- gdata::read.xls(xls=myfn, sheet=1,stringsAsFactors = FALSE)
colnames(temp) <- temp[1,]
temp <- temp[-1,]
temp <- cbind(temp[,1:11], 1, 1/3^8, temp[,12])
colnames(temp) <- colnames(druginfo)
druginfo <-rbind(druginfo,temp)

druginfo <- data.frame("drugid"=toupper(gsub(badchars, "", unique(druginfo$Molecule))), "drug.name"=unique(druginfo$Molecule))


## read sample info
sampleinfo <- read.csv(file.path(path.cell, "tmp", "gne_ge_sampleinfo.map"), sep="\t", header=FALSE)
colnames(sampleinfo) <- c("SAMPLE_ID", "SAMPLE_ALIAS", "BIOSAMPLE_ID", "SAMPLE_TITLE", "ATTRIBUTES")
sampleinfo[!is.na(sampleinfo) & sampleinfo == ""] <- NA
## split the attributes
rr <- apply(sampleinfo[ , "ATTRIBUTES", drop=FALSE], 1, function (x) {
	x <- gsub("90%; DMSO", "90%, DMSO", x)
	x <- unlist(strsplit(x, ";"))
	x <- x[grep("=", x)]
	nn <- sapply(strsplit(x, "="), function (x) { return (x[1]) })
	x <- sapply(strsplit(x, "="), function (x) { return (x[2]) })
	names(x) <- nn
	return (x)
})
cn <- unique(unlist(sapply(rr, names)))
dd <- t(sapply(rr, function (x, y) {
	xx <- array(NA, dim=length(y), dimnames=list(y))
	xx[names(x)] <- x
	return (xx)
}, y=cn))
sampleinfo <- data.frame(sampleinfo[ , !is.element(colnames(sampleinfo), "ATTRIBUTES"), drop=FALSE], dd)
sampleinfo <- data.frame("xpid"=sampleinfo[ , "SAMPLE_ID"], "cellid"=sampleinfo[ , "Cell_line"], "tissueid"=sampleinfo[ , "Primary_Tissue"], sampleinfo)

cellineinfo <- sampleinfo

write.csv(cellineinfo, file=file.path(path.out, "cell_line_annotation_GNE.csv"), row.names=FALSE)
write.csv(druginfo, file=file.path(path.out, "drug_annotation_GNE.csv"), row.names=FALSE)

