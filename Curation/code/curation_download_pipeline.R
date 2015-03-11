## choose your favorite mirror
chooseCRANmirror(graphics=FALSE, ind=15)

## set path to local directory if it is not properly set up
.libPaths(c("/mnt/work1/users/bhklab/Rlib", .libPaths()))

## set method for downloading
options(download.file.method="auto")
# options(download.file.method="wget")
## change to curl, wget or internal depending on your system

## prevent strings to be converted into factors
options(stringsAsFactors=FALSE)

require(affy) || stop("Library affy is not available!")
require(R.utils) || stop("Library R.utils is not available!")
require(stringdist) || stop("Library stringdist is not available!")
require(PharmacoGx) || stop("Library PharmacoGx is not available!")

source(file.path("code", "curation_foo.R"))
source(file.path("code", "lincsAPIQuery.R"))

path.out <- file.path ("output")
if(!file.exists(path.out)) { dir.create(path.out, showWarnings=FALSE, recursive=TRUE) }

badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[\xa0]|[ ]"

mylog <- "curation_download_log.txt"
if(!file.exists(mylog)) {
    steps <- c("curation_download_CGP", "curation_download_CCLE", "curation_download_GSK", "curation_download_NCI60", "curation_download_GRAY", "curation_download_GNE", "curation_download_LINCS", "curation_download_LINCS_HMS")
    progress.log <- cbind(steps, rep("...", length(steps)))
    dimnames(progress.log) <- list(paste("step", 1:length(steps), sep="."), c("script", "progress"))
    write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=mylog, quote=FALSE)
} else {
    progress.log <- read.table(file=mylog, sep="\t", header=TRUE, stringsAsFactor=FALSE)
}

########################
## curation, annotation and normalization of CGP data

message("\n----------------------------------\n| Download information for CGP   |\n----------------------------------")
if (progress.log["step.1", "progress"] != "done") {
    progress.log["step.1", "progress"] <- "in progress"
    write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=mylog, quote=FALSE)
    source(file.path("code", "curation_download_CGP.R"))
    progress.log["step.1", "progress"] <- "done"
    write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=mylog, quote=FALSE)
}
message("\t-> DONE")

#######################
## curation, annotation and normalization of CCLE data

message("\n----------------------------------\n| Download information for CCLE  |\n----------------------------------")
if (progress.log["step.2", "progress"] != "done") {
    progress.log["step.2", "progress"] <- "in progress"
    write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=mylog, quote=FALSE)
    source(file.path("code", "curation_download_CCLE.R"))
    progress.log["step.2", "progress"] <- "done"
    write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=mylog, quote=FALSE)
}
message("\t-> DONE")

#######################
## curation, annotation and normalization of GSK data

message("\n----------------------------------\n| Download information for GSK   |\n----------------------------------")
if (progress.log["step.3", "progress"] != "done") {
    progress.log["step.3", "progress"] <- "in progress"
    write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=mylog, quote=FALSE)
    source(file.path("code", "curation_download_GSK.R"))
    progress.log["step.3", "progress"] <- "done"
    write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=mylog, quote=FALSE)
}
message("\t-> DONE")

#######################
## curation, annotation and normalization of NCI60 data
message("\n----------------------------------\n| Download information for NCI60 |\n----------------------------------")
if (progress.log["step.4", "progress"] != "done") {
    progress.log["step.4", "progress"] <- "in progress"
    write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=mylog, quote=FALSE)
    source(file.path("code", "curation_download_NCI60.R"))
    progress.log["step.4", "progress"] <- "done"
    write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=mylog, quote=FALSE)
}
message("\t-> DONE")

#######################
## curation, annotation and normalization of GRAY data
message("\n----------------------------------\n| Download information for GRAY  |\n----------------------------------")
if (progress.log["step.5", "progress"] != "done") {
    progress.log["step.5", "progress"] <- "in progress"
    write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=mylog, quote=FALSE)
    source(file.path("code", "curation_download_GRAY.R"))
    progress.log["step.5", "progress"] <- "done"
    write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=mylog, quote=FALSE)
}
message("\t-> DONE")

#######################
## curation, annotation and normalization of GNE data
message("\n----------------------------------\n| Download information for GNE  |\n----------------------------------")
if (progress.log["step.6", "progress"] != "done") {
    progress.log["step.6", "progress"] <- "in progress"
    write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=mylog, quote=FALSE)
    source(file.path("code", "curation_download_GNE.R"))
    progress.log["step.6", "progress"] <- "done"
    write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=mylog, quote=FALSE)
}
message("\t-> DONE")

#######################
## curation, annotation and normalization of LINCS data
message("\n----------------------------------\n| Download information for LINCS  |\n----------------------------------")
if (progress.log["step.7", "progress"] != "done") {
    progress.log["step.7", "progress"] <- "in progress"
    write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=mylog, quote=FALSE)
    source(file.path("code", "curation_download_LINCS.R"))
    progress.log["step.7", "progress"] <- "done"
    write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=mylog, quote=FALSE)
}

message("\t-> DONE")

#######################
## curation, annotation and normalization of LINCS-HMS data
message("\n----------------------------------\n| Download information for LINCS-HMS  |\n----------------------------------")
if (progress.log["step.8", "progress"] != "done") {
    progress.log["step.8", "progress"] <- "in progress"
    write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=mylog, quote=FALSE)
    source(file.path("code", "curation_download_LINCS_HMS.R"))
    progress.log["step.8", "progress"] <- "done"
    write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=mylog, quote=FALSE)
}

message("\t-> DONE")

#######################
## merge the annotations

## cell lines
cell.annot <- NULL
cura <- read.csv("matching_cell.csv")
cura[!is.na(cura) & cura == ""] <- NA
rownames(cura) <- cura[ , "unique.cellid"]
## CGP
tt <- read.csv(file.path(path.out, "cell_annotation_CGP.csv"))
cell.annot <- cbind(cell.annot, "CGP.cellid"=tt[ , "cellid"], "CGP.tissueid"=tt[ , "tissueid"])
rownames(cell.annot) <- tt[ , "cellid"]
## CCLE
tt <- read.csv(file.path(path.out, "cell_annotation_CCLE.csv"))
rownames(tt) <- tt[ , "cellid"]
curat <- cura[!is.na(cura[ , "CCLE.cellid"]), , drop=FALSE]
rownames(tt)[match(curat[ , "CCLE.cellid"], rownames(tt))] <- rownames(curat)
myx2 <- setdiff(rownames(tt), rownames(cell.annot))
cell.annot <- rbind(cell.annot, matrix(NA, nrow=length(myx2), ncol=ncol(cell.annot), dimnames=list(myx2, colnames(cell.annot))))
cell.annot <- cbind(cell.annot, "CCLE.cellid"=NA, "CCLE.tissueid"=NA)
cell.annot[rownames(tt), "CCLE.cellid"] <- tt[ , "cellid"]
cell.annot[rownames(tt), "CCLE.tissueid"] <- tt[ , "tissueid"]
## GSK
tt <- read.csv(file.path(path.out, "cell_annotation_GSK.csv"))
rownames(tt) <- tt[ , "cellid"]
curat <- cura[!is.na(cura[ , "GSK.cellid"]), , drop=FALSE]
rownames(tt)[match(curat[ , "GSK.cellid"], rownames(tt))] <- rownames(curat)
myx2 <- setdiff(rownames(tt), rownames(cell.annot))
cell.annot <- rbind(cell.annot, matrix(NA, nrow=length(myx2), ncol=ncol(cell.annot), dimnames=list(myx2, colnames(cell.annot))))
cell.annot <- cbind(cell.annot, "GSK.cellid"=NA, "GSK.tissueid"=NA)
cell.annot[rownames(tt), "GSK.cellid"] <- tt[ , "cellid"]
cell.annot[rownames(tt), "GSK.tissueid"] <- tt[ , "tissueid"]
## NCI60
tt <- read.csv(file.path(path.out, "cell_annotation_NCI60.csv"))
rownames(tt) <- tt[ , "cellid"]
curat <- cura[!is.na(cura[ , "NCI60.cellid"]), , drop=FALSE]
rownames(tt)[match(curat[ , "NCI60.cellid"], rownames(tt))] <- rownames(curat)
myx2 <- setdiff(rownames(tt), rownames(cell.annot))
cell.annot <- rbind(cell.annot, matrix(NA, nrow=length(myx2), ncol=ncol(cell.annot), dimnames=list(myx2, colnames(cell.annot))))
cell.annot <- cbind(cell.annot, "NCI60.cellid"=NA, "NCI60.tissueid"=NA)
cell.annot[rownames(tt), "NCI60.cellid"] <- tt[ , "cellid"]
cell.annot[rownames(tt), "NCI60.tissueid"] <- tt[ , "tissueid"]
## GRAY
tt <- read.csv(file.path(path.out, "cell_annotation_GRAY.csv"))
rownames(tt) <- tt[ , "cellid"]
curat <- cura[!is.na(cura[ , "GRAY.cellid"]), , drop=FALSE]
rownames(tt)[match(curat[ , "GRAY.cellid"], rownames(tt))] <- rownames(curat)
myx2 <- setdiff(rownames(tt), rownames(cell.annot))
cell.annot <- rbind(cell.annot, matrix(NA, nrow=length(myx2), ncol=ncol(cell.annot), dimnames=list(myx2, colnames(cell.annot))))
cell.annot <- cbind(cell.annot, "GRAY.cellid"=NA, "GRAY.tissueid"=NA)
cell.annot[rownames(tt), "GRAY.cellid"] <- tt[ , "cellid"]
cell.annot[rownames(tt), "GRAY.tissueid"] <- tt[ , "tissueid"]
## GNE
tt <- read.csv(file.path(path.out, "cell_annotation_GNE.csv"))
rownames(tt) <- tt[ , "cellid"]
curat <- cura[!is.na(cura[ , "GNE.cellid"]), , drop=FALSE]
rownames(tt)[match(curat[ , "GNE.cellid"], rownames(tt))] <- rownames(curat)
myx2 <- setdiff(rownames(tt), rownames(cell.annot))
cell.annot <- rbind(cell.annot, matrix(NA, nrow=length(myx2), ncol=ncol(cell.annot), dimnames=list(myx2, colnames(cell.annot))))
cell.annot <- cbind(cell.annot, "GNE.cellid"=NA, "GNE.tissueid"=NA)
cell.annot[rownames(tt), "GNE.cellid"] <- tt[ , "cellid"]
cell.annot[rownames(tt), "GNE.tissueid"] <- tt[ , "tissueid"]
## LINCS
tt <- read.csv(file.path(path.out, "cell_annotation_LINCS.csv"))
rownames(tt) <- as.character(tt[ , "cellid"])
curat <- cura[!is.na(cura[ , "LINCS.cellid"]), , drop=FALSE]
rownames(tt)[match(as.character(curat[ , "LINCS.cellid"]), as.character(rownames(tt)))] <- rownames(curat)
myx2 <- setdiff(rownames(tt), rownames(cell.annot))
cell.annot <- rbind(cell.annot, matrix(NA, nrow=length(myx2), ncol=ncol(cell.annot), dimnames=list(myx2, colnames(cell.annot))))
cell.annot <- cbind(cell.annot, "LINCS.cellid"=NA, "LINCS.tissueid"=NA)
cell.annot[rownames(tt), "LINCS.cellid"] <- tt[ , "cellid"]
cell.annot[rownames(tt), "LINCS.tissueid"] <- tt[ , "tissueid"]
## LINCS-HMS
tt <- read.csv(file.path(path.out, "cell_annotation_LINCS_HMS.csv"))
rownames(tt) <- as.character(tt[ , "cellid"])
rownames(tt) [match("5637.0",rownames(tt))] <- c("5637")##Remove the extra 0s
rownames(tt) [match("697.0",rownames(tt))] <- c("697")
curat <- cura[!is.na(cura[ , "LINCS_HMS.cellid"]), , drop=FALSE]
rownames(tt)[match((curat[ , "LINCS_HMS.cellid"]), rownames(tt))] <- rownames(curat)
myx2 <- setdiff(rownames(tt), rownames(cell.annot))
cell.annot <- rbind(cell.annot, matrix(NA, nrow=length(myx2), ncol=ncol(cell.annot), dimnames=list(myx2, colnames(cell.annot))))
cell.annot <- cbind(cell.annot, "LINCS_HMS.cellid"=NA, "LINCS_HMS.tissueid"=NA)
cell.annot[rownames(tt), "LINCS_HMS.cellid"] <- tt[ , "cellid"]
cell.annot[rownames(tt), "LINCS_HMS.tissueid"] <- tt[ , "tissueid"]


## reorderiung of the columns
cell.annot <- cell.annot[ , c(grep("cellid", colnames(cell.annot)), grep("tissueid", colnames(cell.annot)))]
## add COSMIC tissue type
cosmic <- read.csv(file.path(path.out, "cell_annotation_COSMIC.csv"))
rownames(cosmic) <- cosmic[ , "cellid"]
cell.annot <- cbind(cell.annot, "COSMIC.tissueid"=NA)
myx <- intersect(rownames(cell.annot), cosmic[ , "cellid"])
cell.annot[myx, "COSMIC.tissueid"] <- cosmic[myx, "tissueid"]
## find the closest match in cosmic
myx <- setdiff(rownames(cell.annot), rownames(cosmic))
ss <- sapply(myx, function (x, y) {
    dd <- stringdist::stringdist(a=x, b=y)
    dd <- paste(y[!is.na(dd) & dd == min(dd, na.rm=TRUE)], collapse="|")
    if (length(dd) == 0) { dd <- "no_match"}
    return (dd)
}, y=rownames(cosmic))
cell.annot <- cbind(cell.annot, "COSMIC.bestmatch"=NA)
cell.annot[myx, "COSMIC.bestmatch"] <- ss
## create unique tissueid
cell.annot <- cbind(cell.annot, "unique.tissueid"=cell.annot[ , "COSMIC.tissueid"])
## add CCLE tissue type when not found in COSMIC
myx <- is.na(cell.annot[ , "unique.tissueid"]) & !is.na(cell.annot[ , "CCLE.tissueid"])
cell.annot[myx, "unique.tissueid"] <- cell.annot[myx, "CCLE.tissueid"]
## save consolidated annotation file
tt <- cbind("unique.cellid"=rownames(cell.annot), cell.annot)
tt[is.na(tt)] <- ""
write.csv(tt, row.names=FALSE, file=file.path(path.out, "cell_annotation_all_new.csv"))


# myx <- intersect(rownames(cura), cosmic[ , "cellid"])
# cura[myx, "COSMIC.tissueid"] <- cosmic[myx, "tissueid"]
# cura[is.na(cura)] <- ""
# write.csv(cura, "matching_cell.csv")