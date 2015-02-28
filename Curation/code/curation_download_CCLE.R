## path to files
path.data <- file.path("data", "CCLE")
path.drug <- file.path(path.data, "drug")
path.cell <- file.path(path.data, "celline")
## create directories
if(!file.exists(path.data)) { dir.create(path.data, showWarnings=FALSE, recursive=TRUE) }
if(!file.exists(path.drug)) { dir.create(path.drug, showWarnings=FALSE, recursive=TRUE) }
if(!file.exists(path.cell)) { dir.create(path.cell, showWarnings=FALSE, recursive=TRUE) }

## cell line annotations
message("Download cell line annotation")
dir.create(file.path(path.cell, "tmp"), showWarnings=FALSE)
myfn <- file.path(path.cell, "tmp", "ccle_sample_info_file.txt")
if (!file.exists(myfn)) {
  if(!file.exists(file.path(path.cell, "tmp"))) { dir.create(file.path(path.cell, "tmp"), showWarnings=FALSE, recursive=TRUE) }
  dwl.status <- download.file(url="http://www.broadinstitute.org/ccle/downloadFile/DefaultSystemRoot/exp_10/ds_22/CCLE_sample_info_file_2012-10-18.txt?downloadff=true&fileId=6801", destfile=file.path(path.cell, "tmp", "CCLE_sample_info_file_2012-10-18.txt"), quiet=TRUE)
  if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
  file.copy(from=file.path(path.cell, "tmp", "CCLE_sample_info_file_2012-10-18.txt"), to=myfn)
}

## drug info
message("Download drug information")
dir.create(file.path(path.drug, "tmp"), showWarnings=FALSE)
myfn <- file.path(path.drug, "tmp", "ccle_drug_info_file.csv")
if (!file.exists(myfn)) {
  if(!file.exists(file.path(path.drug, "tmp"))) { dir.create(file.path(path.drug, "tmp"), showWarnings=FALSE, recursive=TRUE) }
  dwl.status <- download.file(url="http://www.broadinstitute.org/ccle/downloadFile/DefaultSystemRoot/exp_10/ds_27/CCLE_NP24.2009_profiling_2012.02.20.csv?downloadff=true&fileId=3422", destfile=file.path(path.drug, "tmp", "CCLE_NP24.2009_profiling_2012.02.20.csv"), quiet=TRUE)
  if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
  file.copy(from=file.path(path.drug, "tmp", "CCLE_NP24.2009_profiling_2012.02.20.csv"), to=myfn)
}

## drug pheno
message("Download drug sensitivity measurements")
myfn <- file.path(path.drug, "tmp", "ccle_drug_pheno_file.xls")
if (!file.exists(myfn)) {
  if(!file.exists(file.path(path.drug, "tmp"))) { dir.create(file.path(path.drug, "tmp"), showWarnings=FALSE, recursive=TRUE) }
  ## drug sensitivity data from the addendum in Nature
  dwl.status <- download.file(url="http://www.nature.com/nature/journal/v492/n7428/extref/nature11735-s2.xls", destfile=file.path(path.drug, "tmp", "nature11735-s2.xls"), quiet=TRUE)
  if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
  file.copy(from=file.path(path.drug, "tmp", "nature11735-s2.xls"), to=myfn)
}


message("Read drug information")
druginfo <- read.csv(file.path(path.drug, "tmp", "ccle_drug_info_file.csv"))
druginfo <- apply(druginfo, 2, function (x) { return (gsub("\xa0", "", x)) })
druginfo[druginfo == "" | druginfo == " "] <- NA
druginfo <- data.frame(druginfo, "drugid"=paste("drugid", toupper(gsub(pattern=badchars, replacement="", x=toupper(druginfo[ ,"Compound..code.or.generic.name."]))), sep="_"))
rownames(druginfo) <- druginfo[ ,"drugid"]


message("Read drug sensitivity measurements")
## drug sensitivity data from the addendum in Nature
myfn2 <- file.path(path.drug, "tmp", "ccle_drug_pheno_file.RData")
if(!file.exists(myfn2)) {
  drugpheno <- gdata::read.xls(xls=file.path(path.drug, "tmp", "ccle_drug_pheno_file.xls"), sheet=12)
  drugpheno[drugpheno == "" | drugpheno == " "] <- NA
  save(list="drugpheno", compress=TRUE, file=myfn2)
} else { load(myfn2) }
nn <- apply(drugpheno[1, ], 2, as.character)
nn <- gsub(pattern=badchars, replacement="_", x=nn)
drugpheno <- drugpheno[-c(1),!is.na(nn), drop=FALSE]
nn <- nn[!is.na(nn)]
nn[match(c("Primary_Cell_Line_Name", "IC50_(_M)_(norm)", "EC50_(_M)_(norm)", "Amax_(norm)" ,"ActArea_(norm)", "ActArea_(raw)", "Doses_(uM)", "Activity_Data_(normalized_data)", "Activity_Error_(SD)_(normalized_data)", "FitType_(norm)"), nn)] <- c("Primary.Cell.Line.Name", "IC50..uM.norm",  "EC50..uM.norm", "Amax.norm", "ActArea.norm", "ActArea.raw", "Doses..uM.", "Activity.Data.norm", "Activity.SD.norm", "FitType.norm")
colnames(drugpheno) <- nn
drugpheno <- setcolclass.df(df=drugpheno, colclass=rep("character", ncol(drugpheno)), factor.levels=sapply(drugpheno, levels))
drugpheno <- setcolclass.df(df=drugpheno, colclass=c(rep("character", 8), rep("numeric", 4), "character", "numeric", "numeric"))
drugpheno[ ,"Compound"] <- toupper(gsub(pattern=badchars, replacement="", drugpheno[ ,"Compound"]))

drugpheno <- data.frame(drugpheno, "drugid"=paste("drugid", toupper(gsub(pattern=badchars, replacement="", x=as.character(drugpheno[ ,"Compound"]))), sep="_"), "cellid"=as.character(drugpheno[ ,"Primary.Cell.Line.Name"]))
## drug screening concentrations
ll <- sapply(strsplit(drugpheno[ , "Doses..uM."], ","), length)
drugconc <- lapply(strsplit(drugpheno[ , "Doses..uM."], ","), function(x) {
  xx <- as.numeric(x)
  xx2 <- c(0.0025, 0.0080, 0.0250, 0.0800, 0.2500, 0.8000, 2.5300, 8.0000)
  if(any(!is.element(xx, xx2))) { stop("Unexpected drug screening concentrations!") }
  xx3 <- rep(NA, length(xx2))
  names(xx3) <- xx2
  xx3[match(xx, names(xx3))] <- xx
  return(xx3)
  })
drugconc <- do.call(rbind, drugconc)
drugconc <- data.frame("cellid"=as.character(drugpheno[ , "cellid"]), "drugid"=as.character(drugpheno[ , "drugid"]), "nbr.conc.tested"=ll, drugconc)
dimnames(drugconc) <- list(paste(drugconc[ , "drugid"], drugconc[ , "cellid"], sep="..."), c("cellid", "drugid", "nbr.conc.tested", sprintf("Dose%i.uM", 1:(ncol(drugconc) - 3))))


## info about each experiment
message("Read sample information")
sampleinfo <- read.csv(file.path(path.cell, "tmp", "ccle_sample_info_file.txt"), sep="\t")
sampleinfo[sampleinfo == "" | sampleinfo == " "] <- NA
sampleinfo <- data.frame("cellid"=as.character(sampleinfo[ ,"Cell.line.primary.name"]), sampleinfo)
## remove duplicated cell line hybridization
## only the most recent experiment (as determine by hyridization date or CEL file timestamp) will be kept for each cell line
sampleinfo <- sampleinfo[!duplicated(sampleinfo[ ,"cellid"]), ,drop=FALSE]
sampleinfo[ , "cellid"] <- as.character(sampleinfo[ , "cellid"])
rownames(sampleinfo) <- as.character(sampleinfo[ , "cellid"])

## cell lines and drugs
cellnall <- sort(unique(c(as.character(sampleinfo[ ,"cellid"]), as.character(drugpheno[ , "cellid"]), as.character(drugconc[ , "cellid"]))))
drugnall <- sort(unique(c(rownames(druginfo), as.character(drugconc[ , "drugid"]), as.character(drugpheno[ , "drugid"]))))

## update sampleinfo
dd <- data.frame(matrix(NA, nrow=length(cellnall), ncol=ncol(sampleinfo), dimnames=list(cellnall, colnames(sampleinfo))))
newlev <- sapply(sampleinfo, levels)
newlev$cellid <- cellnall
dd <- genefu::setcolclass.df(df=dd, colclass=sapply(sampleinfo, class), factor.levels=newlev)
dd[rownames(sampleinfo), colnames(sampleinfo)] <- sampleinfo
dd[ , "cellid"] <- rownames(dd)
dd <- data.frame(dd, "tissueid"=as.character(dd[ , "Site.Primary"]))
sampleinfo <- dd

## update druginfo
druginfo <- data.frame("drug.name"=druginfo[ , "Compound..code.or.generic.name."], druginfo)
dd <- data.frame(matrix(NA, nrow=length(drugnall), ncol=ncol(druginfo), dimnames=list(drugnall, colnames(druginfo))))
dd <- genefu::setcolclass.df(df=dd, colclass=sapply(druginfo, class), factor.levels=newlev)
dd[match(rownames(druginfo), drugnall), colnames(druginfo)] <- druginfo
myx <- is.na(dd[ , "Compound..code.or.generic.name."])
dd[myx, "Compound..code.or.generic.name."] <- dd[myx, "drugid"]
dd[ , "drugid"] <- gsub("drugid_", "", drugnall)
druginfo <- dd
myx <- match(c("drugid", "drug.name"), colnames(druginfo))
druginfo <- druginfo[ , c(myx, setdiff(1:ncol(druginfo), myx)), drop=FALSE]


write.csv(sampleinfo, file=file.path(path.out, "cell_line_annotation_CCLE.csv"), row.names=FALSE)
write.csv(druginfo, file=file.path(path.out, "drug_annotation_CCLE.csv"), row.names=FALSE)