## path to files
path.data <- file.path("data", "SU2C")
path.drug <- file.path(path.data, "drug")
path.cell <- file.path(path.data, "celline")
## create directories
if(!file.exists(path.data)) { dir.create(path.data, showWarnings=FALSE, recursive=TRUE) }
if(!file.exists(path.drug)) { dir.create(path.drug, showWarnings=FALSE, recursive=TRUE) }
if(!file.exists(path.cell)) { dir.create(path.cell, showWarnings=FALSE, recursive=TRUE) }


## sample information
message("Download sample information")
myfn <- file.path(path.cell, "tmp", "su2c_ge_sampleinfo.txt")
if (!file.exists(myfn)) {
  if(!file.exists(file.path(path.cell, "tmp"))) { dir.create(file.path(path.cell, "tmp"), showWarnings=FALSE, recursive=TRUE) }
  ## drug sensitivity data from the addendum in Nature
  dwl.status <- download.file(url="http://www.ebi.ac.uk/arrayexpress/files/E-MTAB-181/E-MTAB-181.sdrf.txt", destfile=file.path(path.cell, "tmp", "E-MTAB-181.sdrf.txt"), quiet=TRUE)
  if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
  file.copy(from=file.path(path.cell, "tmp", "E-MTAB-181.sdrf.txt"), to=myfn)
}

## cell line information
message("Download cell line information")
myfn <- file.path(path.cell, "tmp", "su2c_cellineinfo.xls")
if (!file.exists(myfn)) {
  if(!file.exists(file.path(path.cell, "tmp"))) { dir.create(file.path(path.cell, "tmp"), showWarnings=FALSE, recursive=TRUE) }
  ## drug sensitivity data from the addendum in Nature
  dwl.status <- download.file(url="http://www.pnas.org/content/suppl/2011/10/14/1018854108.DCSupplemental/sd01.xls", destfile=file.path(path.cell, "tmp", "sd01.xls"), quiet=TRUE)
  if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
  file.copy(from=file.path(path.cell, "tmp", "sd01.xls"), to=myfn)
}

## drug pheno
message("Download drug sensitivity measurements")
myfn <- file.path(path.drug, "tmp", "su2c_drug_pheno.xlsx")
if (!file.exists(myfn)) {
  if(!file.exists(file.path(path.drug, "tmp"))) { dir.create(file.path(path.drug, "tmp"), showWarnings=FALSE, recursive=TRUE) }
  ## drug sensitivity data from the addendum in Nature
  dwl.status <- download.file(url="http://www.pnas.org/content/suppl/2011/10/14/1018854108.DCSupplemental/sd02.xlsx", destfile=file.path(path.drug, "tmp", "sd02.xlsx"), quiet=TRUE)
  if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
  file.copy(from=file.path(path.drug, "tmp", "sd02.xlsx"), to=myfn)
}


## read sample info
sampleinfo <- read.csv(file.path(path.cell, "tmp", "su2c_ge_sampleinfo.txt"), sep="\t")
sampleinfo[!is.na(sampleinfo) & sampleinfo == ""] <- NA
sampleinfo <- data.frame("cellid"=sampleinfo[ , "Characteristics..CellLine."], sampleinfo)
rownames(sampleinfo) <- gsub("[.CEL]", "", sampleinfo[ , "Array.Data.File"])

## read cell line info
cellineinfo <- gdata::read.xls(xls=file.path(path.cell, "tmp", "su2c_cellineinfo.xls"), sheet=1)
cellineinfo[!is.na(cellineinfo) & cellineinfo == ""] <- NA
cellineinfo <- cellineinfo[1:(min(which(cellineinfo[ ,1] == "a")) - 1), , drop=FALSE]
rn <- cellineinfo[-1, 1]
cn <- cellineinfo[1, -1]
cn <- gsub(badchars, ".", cn)
cellineinfo <- cellineinfo[-1, -1]
dimnames(cellineinfo) <- list(rn, cn)
cellineinfo <- data.frame("cellid"=rn, "tissueid"="breast", cellineinfo)

## read drug pheno
drugpheno <- gdata::read.xls(file.path(path.drug, "tmp", "su2c_drug_pheno.xlsx"), sheet=1)
drugpheno[!is.na(drugpheno) & drugpheno == ""] <- NA
nn <- drugpheno[-1, 1]
druginfo <- data.frame("drugid"=toupper(gsub(badchars, "", unlist(drugpheno[1, -1, drop=T]))), "drug.name"=unlist(drugpheno[1, -1, drop=T]))
myx <- match(c("drugid", "drug.name"), colnames(druginfo))
druginfo <- druginfo[ , c(myx, setdiff(1:ncol(druginfo), myx)), drop=FALSE]

drugpheno <- drugpheno[-1, -1]
dimnames(drugpheno) <- list(nn, druginfo[ , "drugid"])

## merge all cell lines
cellnall <- fold(union, sampleinfo[ , "cellid"], cellineinfo[ , "cellid"], rownames(drugpheno))

## update cell line info
dd <- data.frame(matrix(NA, nrow=length(cellnall), ncol=ncol(cellineinfo), dimnames=list(cellnall, colnames(cellineinfo))))
dd[rownames(cellineinfo), colnames(cellineinfo)] <- cellineinfo
dd[ , "cellid"] <- cellnall
dd["tissueid"] <- "breast"
cellineinfo <- dd

write.csv(cellineinfo, file=file.path(path.out, "cell_line_annotation_SU2C.csv"), row.names=FALSE)
write.csv(druginfo, file=file.path(path.out, "drug_annotation_SU2C.csv"), row.names=FALSE)

