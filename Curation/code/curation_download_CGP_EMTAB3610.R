## path to files
path.data <- file.path("data", "CGP_EMTAB3610")
path.cell <- file.path(path.data, "cell")

## create directories
if(!file.exists(path.data)) { dir.create(path.data, showWarnings=FALSE, recursive=TRUE) }
if(!file.exists(path.cell)) { dir.create(path.cell, showWarnings=FALSE, recursive=TRUE) }



ftpdir <- "https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-3610"

myfn <- file.path(path.cell, "tmp", "E-MTAB-3610.sdrf.txt")
if(!file.exists(file.path(path.cell, "tmp"))){dir.create(file.path(path.cell, "tmp"), showWarnings=FALSE)}
if (!file.exists(myfn)) {
  dwl.status <- download.file(url=sprintf("%s/E-MTAB-3610.sdrf.txt", ftpdir), destfile=file.path(path.cell, "tmp", "E-MTAB-3610.sdrf.txt"), quiet=TRUE)
  if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
}

cellineinfo <- read.csv(file = myfn,header = TRUE, sep = "\t",stringsAsFactors = FALSE, na.strings = c(""," "))
tt <- strsplit(cellineinfo$Source.Name, split="_")
tissue <- sapply(tt, function(x) {
  if(length(x) > 5) {
    y <- paste(x[4:(length(x)-1)], collapse="_")
  }else{
    y <- x[4]
  }
}, simplify=TRUE)
cellineinfo <- data.frame("cellid"=cellineinfo$Characteristics.cell.line., "tissueid"=tissue, cellineinfo)
cellineinfo <- cellineinfo[which(!duplicated(cellineinfo$cellid)), , drop=FALSE]
rownames(cellineinfo) <- cellineinfo$Source.Name

write.csv(cellineinfo, file=file.path(path.out, "cell_line_annotation_CGP_EMTAB3610.csv"), row.names=FALSE)

