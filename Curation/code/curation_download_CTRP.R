## path to files
path.data <- file.path("data", "CTRP")
path.drug <- file.path(path.data, "drug")
path.cell <- file.path(path.data, "cell")

## create directories
if(!file.exists(path.data)) { dir.create(path.data, showWarnings=FALSE, recursive=TRUE) }
if(!file.exists(path.drug)) { dir.create(path.drug, showWarnings=FALSE, recursive=TRUE) }
if(!file.exists(path.cell)) { dir.create(path.cell, showWarnings=FALSE, recursive=TRUE) }



ftpdir <- "ftp://caftpd.nci.nih.gov/pub/dcc_ctd2/Broad/Cancer_Cell_Line_Profiling/"

myfn <- file.path(path.cell, "tmp", "Broad.CTD2.M2.cell_line_info.txt")
if(!file.exists(file.path(path.cell, "tmp"))){dir.create(file.path(path.cell, "tmp"), showWarnings=FALSE)}
if (!file.exists(myfn)) {
  dwl.status <- download.file(url=sprintf("%s/Broad.CTD2.CellResource.data_package.zip", ftpdir), destfile=file.path(path.cell, "tmp", "Broad.CTD2.CellResource.data_package.zip"), quiet=TRUE)
  if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
  unzip(zipfile=file.path(path.cell, "tmp", "Broad.CTD2.CellResource.data_package.zip"), exdir=file.path(path.cell, "tmp"))
}

cellineinfo <- read.csv(file = myfn,header = TRUE, sep = "\t",stringsAsFactors = FALSE, na.strings = c(""," "))
cellineinfo <- data.frame("cellid"=cellineinfo$cell_line_name, "tissueid"=cellineinfo$cell_line_lineage, cellineinfo)
rownames(cellineinfo) <- cellineinfo$cell_line_name

myfn <- file.path(path.drug, "tmp", "Broad.CTD2.M1.informer_set.txt")
if(!file.exists(file.path(path.drug, "tmp"))){dir.create(file.path(path.drug, "tmp"), showWarnings=FALSE)}
if (!file.exists(myfn)) {
  dwl.status <- download.file(url=sprintf("%s/Broad.CTD2.CellResource.data_package.zip", ftpdir), destfile=file.path(path.drug, "tmp", "Broad.CTD2.CellResource.data_package.zip"), quiet=TRUE)
  if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
  unzip(zipfile=file.path(path.drug, "tmp", "Broad.CTD2.CellResource.data_package.zip"), exdir=file.path(path.drug, "tmp"))
}

druginfo <- read.csv(file = myfn,header = TRUE, sep = "\t",stringsAsFactors = FALSE, na.strings = c(""," "))
druginfo <- data.frame("drugid"=toupper(gsub(badchars, "", unique(druginfo$compound_name))), "drug.name"=unique(sampleinfo$compound_name), druginfo)
rownames(druginfo) <- druginfo$compound_name

myfn <- file.path(path.drug, "tmp", "Broad.CTD2.D2.avg_pct_viability_data.txt")
if(!file.exists(file.path(path.drug, "tmp"))){dir.create(file.path(path.drug, "tmp"), showWarnings=FALSE)}
if (!file.exists(myfn)) {
  dwl.status <- download.file(url=sprintf("%s/Broad.CTD2.CellResource.data_package.zip", ftpdir), destfile=file.path(path.drug, "tmp", "Broad.CTD2.CellResource.data_package.zip"), quiet=TRUE)
  if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
  unzip(zipfile=file.path(path.drug, "tmp", "Broad.CTD2.CellResource.data_package.zip"), exdir=file.path(path.drug, "tmp"))
}

sampleinfo <- read.csv(file = myfn,header = TRUE, sep = "\t",stringsAsFactors = FALSE, na.strings = c(""," "))
sampleinfo <- data.frame("cellid"=sampleinfo$cell_line_name, "drug.name"=sampleinfo$compound_name, sampleinfo)


## merge all cell lines
cellnall <- fold(union, sampleinfo[ , "cellid"], cellineinfo[ , "cellid"])

## update cell line info
dd <- data.frame(matrix(NA, nrow=length(cellnall), ncol=ncol(cellineinfo), dimnames=list(cellnall, colnames(cellineinfo))))
dd[rownames(cellineinfo), colnames(cellineinfo)] <- cellineinfo
dd[ , "cellid"] <- cellnall
cellineinfo <- dd


## merge all drugs
drugnall <- fold(union, sampleinfo[ , "drug.name"], druginfo[ , "drug.name"])

## update cell line info
dd <- data.frame(matrix(NA, nrow=length(drugnall), ncol=ncol(druginfo), dimnames=list(drugnall, colnames(druginfo))))
dd[rownames(druginfo), colnames(druginfo)] <- druginfo
dd[ , "drug.name"] <- drugnall
druginfo <- dd

write.csv(cellineinfo, file=file.path(path.out, "cell_line_annotation_CTRP.csv"), row.names=FALSE)
write.csv(druginfo, file=file.path(path.out, "drug_annotation_CTRP.csv"), row.names=FALSE)

