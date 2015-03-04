## path to files
path.data <- file.path("data", "GRAY")
path.drug <- file.path(path.data, "drug")
path.cell <- file.path(path.data, "celline")
## create directories
if(!file.exists(path.data)) { dir.create(path.data, showWarnings=FALSE, recursive=TRUE) }
if(!file.exists(path.drug)) { dir.create(path.drug, showWarnings=FALSE, recursive=TRUE) }
if(!file.exists(path.cell)) { dir.create(path.cell, showWarnings=FALSE, recursive=TRUE) }


## sample information
message("Download sample information")
archivn <- "gb-2013-14-10-r110-s9.txt"
myfn <- file.path(path.cell, "tmp", archivn)

if (!file.exists(myfn)) {
  if(!file.exists(file.path(path.cell, "tmp"))) { dir.create(file.path(path.cell, "tmp"), showWarnings=FALSE, recursive=TRUE) }
  ## drug sensitivity data from the addendum in Nature
  dwl.status <- download.file(url=sprintf("http://genomebiology.com/content/supplementary/%s", archivn), destfile=myfn, quiet=TRUE)
  if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
}

## read sample info
sampleinfo <- read.csv(file = myfn,header = TRUE, sep = "\t",stringsAsFactors = FALSE)
sampleinfo[!is.na(sampleinfo) & sampleinfo == ""] <- NA
sampleinfo <- data.frame("cellid"=sampleinfo$cellline, sampleinfo)

## read cell line info
message("Download cell line information")

archivn <- "gb-2013-14-10-r110-s1.xlsx"
myfn <- file.path(path.cell, "tmp", archivn)

if (!file.exists(myfn)) {
  if(!file.exists(file.path(path.cell, "tmp"))) { dir.create(file.path(path.cell, "tmp"), showWarnings=FALSE, recursive=TRUE) }
  ## drug sensitivity data from the addendum in Nature
  dwl.status <- download.file(url=sprintf("http://genomebiology.com/content/supplementary/%s", archivn), destfile=myfn, quiet=TRUE)
  if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
}
cellineinfo <- gdata::read.xls(xls=myfn, sheet=1)
cellineinfo[!is.na(cellineinfo) & cellineinfo == ""] <- NA
cellineinfo <- cellineinfo[1:(min(grep("For", cellineinfo[,1])) - 1), , drop=FALSE]
rn <- cellineinfo[-1, 1]
cn <- t(cellineinfo[1, -1])
cn <- gsub(badchars, ".", cn)
cellineinfo <- cellineinfo[-1, -1]
dimnames(cellineinfo) <- list(rn, cn)
cellineinfo <- data.frame("cellid"=rn, "tissueid"="breast", cellineinfo[,1:10])

druginfo <- data.frame("drugid"=toupper(gsub(badchars, "", unique(sampleinfo$compound))), "drug.name"=unique(sampleinfo$compound))

## merge all cell lines
cellnall <- fold(union, sampleinfo[ , "cellid"], cellineinfo[ , "cellid"])

## update cell line info
dd <- data.frame(matrix(NA, nrow=length(cellnall), ncol=ncol(cellineinfo), dimnames=list(cellnall, colnames(cellineinfo))))
dd[rownames(cellineinfo), colnames(cellineinfo)] <- cellineinfo
dd[ , "cellid"] <- cellnall
dd["tissueid"] <- "breast"
cellineinfo <- dd

write.csv(cellineinfo, file=file.path(path.out, "cell_line_annotation_GRAY.csv"), row.names=FALSE)
write.csv(druginfo, file=file.path(path.out, "drug_annotation_GRAY.csv"), row.names=FALSE)

