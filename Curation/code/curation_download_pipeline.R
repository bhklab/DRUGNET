## choose your favorite mirror
chooseCRANmirror(graphics=FALSE, ind=15)

## set path to local directory if it is not properly set up
## on Mordor
# .libPaths(c("/mnt/work1/users/bhklab/Rlib", .libPaths()))
## on Guillimin (module add R)
.libPaths(c("/sb/project/afb-431/Rlib", .libPaths()))

## set method for downloading
options(download.file.method="auto")
# options(download.file.method="wget")
## change to curl, wget or internal depending on your system

## prevent strings to be converted into factors
options(stringsAsFactors=FALSE)

#require(affy) || stop("Library affy is not available!")
require(affyio) || stop("Library affyio is not available!")
require(R.utils) || stop("Library R.utils is not available!")
require(stringdist) || stop("Library stringdist is not available!")
require(RCurl) || stop("Library RCurl is not available!")
require(genefu)|| stop("Library genefu is not available!")
require(XML)|| stop("Library XML is not available!")
require(Hmisc)|| stop("Library Hmisc is not available!")
require(gdata)|| stop("Library gdata is not available!")
# require(PharmacoGx) || stop("Library PharmacoGx is not available!")

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
cura <- read.csv("matching_cell.csv", na.strings=c("NA","NaN", " ",""))
#cura[!is.na(cura) & cura == ""] <- NA
rownames(cura) <- cura[ , "unique.cellid"]
## CGP
tt <- read.csv(file.path(path.out, "cell_line_annotation_CGP.csv"))
tt <- subset(tt, !is.na(tt[ , "cellid"])&!duplicated(tt[ , "cellid"]))
rownames(tt) <- tt[ , "cellid"]
curat <- cura[!is.na(cura[ , "CGP.cellid"]), , drop=FALSE]
rownames(tt)[match(curat[ , "CGP.cellid"], rownames(tt))] <- rownames(curat)
cell.annot <- cbind(cell.annot, "CGP.cellid"=tt[ , "cellid"], "CGP.tissueid"=tt[ , "tissueid"])
rownames(cell.annot) <- rownames(tt)
## CCLE
tt <- read.csv(file.path(path.out, "cell_line_annotation_CCLE.csv"))
rownames(tt) <- tt[ , "cellid"]
curat <- cura[!is.na(cura[ , "CCLE.cellid"]), , drop=FALSE]
rownames(tt)[match(curat[ , "CCLE.cellid"], rownames(tt))] <- rownames(curat)
myx2 <- setdiff(rownames(tt), rownames(cell.annot))
cell.annot <- rbind(cell.annot, matrix(NA, nrow=length(myx2), ncol=ncol(cell.annot), dimnames=list(myx2, colnames(cell.annot))))
cell.annot <- cbind(cell.annot, "CCLE.cellid"=NA, "CCLE.tissueid"=NA)
cell.annot[rownames(tt), "CCLE.cellid"] <- tt[ , "cellid"]
cell.annot[rownames(tt), "CCLE.tissueid"] <- tt[ , "tissueid"]
## GSK
tt <- read.csv(file.path(path.out, "cell_line_annotation_GSK.csv"))
rownames(tt) <- tt[ , "cellid"]
curat <- cura[!is.na(cura[ , "GSK.cellid"]), , drop=FALSE]
rownames(tt)[match(curat[ , "GSK.cellid"], rownames(tt))] <- rownames(curat)
myx2 <- setdiff(rownames(tt), rownames(cell.annot))
cell.annot <- rbind(cell.annot, matrix(NA, nrow=length(myx2), ncol=ncol(cell.annot), dimnames=list(myx2, colnames(cell.annot))))
cell.annot <- cbind(cell.annot, "GSK.cellid"=NA, "GSK.tissueid"=NA)
cell.annot[rownames(tt), "GSK.cellid"] <- tt[ , "cellid"]
cell.annot[rownames(tt), "GSK.tissueid"] <- tt[ , "tissueid"]
## NCI60
tt <- read.csv(file.path(path.out, "cell_line_annotation_NCI60.csv"))
rownames(tt) <- tt[ , "cellid"]
curat <- cura[!is.na(cura[ , "NCI60.cellid"]), , drop=FALSE]
rownames(tt)[match(curat[ , "NCI60.cellid"], rownames(tt))] <- rownames(curat)
myx2 <- setdiff(rownames(tt), rownames(cell.annot))
cell.annot <- rbind(cell.annot, matrix(NA, nrow=length(myx2), ncol=ncol(cell.annot), dimnames=list(myx2, colnames(cell.annot))))
cell.annot <- cbind(cell.annot, "NCI60.cellid"=NA, "NCI60.tissueid"=NA)
cell.annot[rownames(tt), "NCI60.cellid"] <- tt[ , "cellid"]
cell.annot[rownames(tt), "NCI60.tissueid"] <- tt[ , "tissueid"]
## GRAY
tt <- read.csv(file.path(path.out, "cell_line_annotation_GRAY.csv"))
rownames(tt) <- tt[ , "cellid"]
curat <- cura[!is.na(cura[ , "GRAY.cellid"]), , drop=FALSE]
rownames(tt)[match(curat[ , "GRAY.cellid"], rownames(tt))] <- rownames(curat)
myx2 <- setdiff(rownames(tt), rownames(cell.annot))
cell.annot <- rbind(cell.annot, matrix(NA, nrow=length(myx2), ncol=ncol(cell.annot), dimnames=list(myx2, colnames(cell.annot))))
cell.annot <- cbind(cell.annot, "GRAY.cellid"=NA, "GRAY.tissueid"=NA)
cell.annot[rownames(tt), "GRAY.cellid"] <- tt[ , "cellid"]
cell.annot[rownames(tt), "GRAY.tissueid"] <- tt[ , "tissueid"]
## GNE
tt <- read.csv(file.path(path.out, "cell_line_annotation_GNE.csv"))
tt <- subset(tt, !is.na(tt[ , "cellid"]))
rownames(tt) <- tt[ , "cellid"]
curat <- cura[!is.na(cura[ , "GNE.cellid"]), , drop=FALSE]
rownames(tt)[match(curat[ , "GNE.cellid"], rownames(tt))] <- rownames(curat)
myx2 <- setdiff(rownames(tt), rownames(cell.annot))
cell.annot <- rbind(cell.annot, matrix(NA, nrow=length(myx2), ncol=ncol(cell.annot), dimnames=list(myx2, colnames(cell.annot))))
cell.annot <- cbind(cell.annot, "GNE.cellid"=NA, "GNE.tissueid"=NA)
cell.annot[rownames(tt), "GNE.cellid"] <- tt[ , "cellid"]
cell.annot[rownames(tt), "GNE.tissueid"] <- tt[ , "tissueid"]
## CTRP
tt <- read.csv(file.path(path.out, "cell_line_annotation_CTRP.csv"))
tt <- subset(tt, !is.na(tt[ , "cellid"]))
rownames(tt) <- tt[ , "cellid"]
curat <- cura[!is.na(cura[ , "CTRP.cellid"]), , drop=FALSE]
rownames(tt)[match(curat[ , "CTRP.cellid"], rownames(tt))] <- rownames(curat)
myx2 <- setdiff(rownames(tt), rownames(cell.annot))
cell.annot <- rbind(cell.annot, matrix(NA, nrow=length(myx2), ncol=ncol(cell.annot), dimnames=list(myx2, colnames(cell.annot))))
cell.annot <- cbind(cell.annot, "CTRP.cellid"=NA, "CTRP.tissueid"=NA)
cell.annot[rownames(tt), "CTRP.cellid"] <- tt[ , "cellid"]
cell.annot[rownames(tt), "CTRP.tissueid"] <- tt[ , "tissueid"]
## CGP EMTAB3610
tt <- read.csv(file.path(path.out, "cell_line_annotation_CGP_EMTAB3610.csv"))
tt <- subset(tt, !is.na(tt[ , "cellid"]))
rownames(tt) <- tt[ , "cellid"]
curat <- cura[!is.na(cura[ , "CGP_EMTAB3610.cellid"]), , drop=FALSE]
rownames(tt)[match(curat[ , "CGP_EMTAB3610.cellid"],rownames(tt))] <- rownames(curat)
myx2 <- setdiff(rownames(tt), rownames(cell.annot))
cell.annot <- rbind(cell.annot, matrix(NA, nrow=length(myx2), ncol=ncol(cell.annot), dimnames=list(myx2, colnames(cell.annot))))
cell.annot <- cbind(cell.annot, "CGP_EMTAB3610.cellid"=NA, "CGP_EMTAB3610.tissueid"=NA)
cell.annot[rownames(tt), "CGP_EMTAB3610.cellid"] <- tt[ , "cellid"]
cell.annot[rownames(tt), "CGP_EMTAB3610.tissueid"] <- tt[ , "tissueid"]
## CCLE rnaseq
tt <- read.csv(file.path(path.out, "cell_line_annotation_CCLE_rnaseq.csv"))
tt <- subset(tt, !is.na(tt[ , "cellid"]))
rownames(tt) <- tt[ , "cellid"]
curat <- cura[!is.na(cura[ , "CCLE_rnaseq.cellid"]), , drop=FALSE]
rownames(tt)[match(curat[ , "CCLE_rnaseq.cellid"],rownames(tt))] <- rownames(curat)
myx2 <- setdiff(rownames(tt), rownames(cell.annot))
cell.annot <- rbind(cell.annot, matrix(NA, nrow=length(myx2), ncol=ncol(cell.annot), dimnames=list(myx2, colnames(cell.annot))))
cell.annot <- cbind(cell.annot, "CCLE_rnaseq.cellid"=NA, "CCLE_rnaseq.tissueid"=NA)
cell.annot[rownames(tt), "CCLE_rnaseq.cellid"] <- tt[ , "cellid"]
cell.annot[rownames(tt), "CCLE_rnaseq.tissueid"] <- tt[ , "tissueid"]
## GDSC SNP
tt <- read.csv(file.path(path.out, "cell_line_annotation_GDSC_SNP.csv"))
tt <- subset(tt, !is.na(tt[ , "cellid"]))
rownames(tt) <- tt[ , "cellid"]
curat <- cura[!is.na(cura[ , "GDSC.SNP.cellid"]), , drop=FALSE]
rownames(tt)[match(curat[ , "GDSC.SNP.cellid"],rownames(tt))] <- rownames(curat)
myx2 <- setdiff(rownames(tt), rownames(cell.annot))
cell.annot <- rbind(cell.annot, matrix(NA, nrow=length(myx2), ncol=ncol(cell.annot), dimnames=list(myx2, colnames(cell.annot))))
cell.annot <- cbind(cell.annot, "GDSC.SNP.cellid"=NA, "GDSC.SNP.tissueid"=NA)
cell.annot[rownames(tt), "GDSC.SNP.cellid"] <- tt[ , "cellid"]
cell.annot[rownames(tt), "GDSC.SNP.tissueid"] <- tt[ , "tissueid"]
## Ben Neel
tt <- read.csv(file.path(path.out, "cell_line_annotation_BenNeel.csv"))
tt <- subset(tt, !is.na(tt[ , "cellid"]))
rownames(tt) <- tt[ , "cellid"]
curat <- cura[!is.na(cura[ , "Ben_Neel.cellid"]), , drop=FALSE]
rownames(tt)[match(curat[ , "Ben_Neel.cellid"],rownames(tt))] <- rownames(curat)
myx2 <- setdiff(rownames(tt), rownames(cell.annot))
cell.annot <- rbind(cell.annot, matrix(NA, nrow=length(myx2), ncol=ncol(cell.annot), dimnames=list(myx2, colnames(cell.annot))))
cell.annot <- cbind(cell.annot, "Ben_Neel.cellid"=NA, "Ben_Neel.tissueid"=NA)
cell.annot[rownames(tt), "Ben_Neel.cellid"] <- tt[ , "cellid"]
cell.annot[rownames(tt), "Ben_Neel.tissueid"] <- tt[ , "tissueid"]

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
cosmic <- read.csv(file.path(path.out, "cell_line_annotation_COSMIC.csv"))
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

myx <- is.na(cell.annot[ , "unique.tissueid"]) & !is.na(cell.annot[ , "GSK.tissueid"])
tt <- cell.annot[myx, "GSK.tissueid"] 
tt[which(tt == "hematopoietic_and_lymphatic_system")] <- "haematopoietic_and_lymphoid_tissue"
tt[which(tt == "cervix_uteri")] <- "cervix"
cell.annot[myx, "unique.tissueid"] <- tt

myx <- is.na(cell.annot[ , "unique.tissueid"]) & !is.na(cell.annot[ , "NCI60.tissueid"])
cell.annot[myx, "unique.tissueid"] <- ifelse(cell.annot[myx, "NCI60.tissueid"] == "ME", "skin", cell.annot[myx, "NCI60.tissueid"])

myx <- is.na(cell.annot[ , "unique.tissueid"]) & !is.na(cell.annot[ , "GRAY.tissueid"])
cell.annot[myx, "unique.tissueid"] <- cell.annot[myx, "GRAY.tissueid"]

myx <- is.na(cell.annot[ , "unique.tissueid"]) & !is.na(cell.annot[ , "GNE.tissueid"])
tt <- tolower(cell.annot[myx, "GNE.tissueid"])
tt[which(tt == "adrenal")] <- "adrenal_gland"
tt[which(tt == "blood")] <- "haematopoietic_and_lymphoid_tissue"
tt[which(tt == "lymph node")] <- "haematopoietic_and_lymphoid_tissue"
tt[which(tt == "skeletal muscle")] <- "skeletal_muscle"
tt[which(tt == "oral cavity")] <- "oral_cavity"
cell.annot[myx, "unique.tissueid"] <- tt

myx <- is.na(cell.annot[ , "unique.tissueid"]) & !is.na(cell.annot[ , "CTRP.tissueid"])
cell.annot[myx, "unique.tissueid"] <- cell.annot[myx, "CTRP.tissueid"]

myx <- is.na(cell.annot[ , "unique.tissueid"]) & !is.na(cell.annot[ , "CGP_EMTAB3610.tissueid"])
tt <- tolower(cell.annot[myx, "CGP_EMTAB3610.tissueid"])
tt[which(tt == "0")] <- "other"
tt[which(tt == "blood")] <- "haematopoietic_and_lymphoid_tissue"
tt[which(tt == "head & neck")] <- "head_and_neck"
cell.annot[myx, "unique.tissueid"] <- tt

myx <- is.na(cell.annot[ , "unique.tissueid"]) & !is.na(cell.annot[ , "CCLE_rnaseq.tissueid"])
cell.annot[myx, "unique.tissueid"] <- cell.annot[myx, "CCLE_rnaseq.tissueid"]

myx <- is.na(cell.annot[ , "unique.tissueid"]) & !is.na(cell.annot[ , "GDSC.SNP.tissueid"])
cell.annot[myx, "unique.tissueid"] <- cell.annot[myx, "GDSC.SNP.tissueid"]

myx <- is.na(cell.annot[ , "unique.tissueid"]) & !is.na(cell.annot[ , "Ben_Neel.tissueid"])
cell.annot[myx, "unique.tissueid"] <- cell.annot[myx, "Ben_Neel.tissueid"]

myx <- is.na(cell.annot[ , "unique.tissueid"]) & !is.na(cell.annot[ , "CGP.tissueid"])
cell.annot[myx, "unique.tissueid"] <- cell.annot[myx, "CGP.tissueid"]

##manual tissue assignment for those without any tissue type
cell.annot["CRO-AP3", "unique.tissueid"] <- "haematopoietic_and_lymphoid_tissue"
cell.annot["NTERA-2_cl.D1", "unique.tissueid"] <- "testis"
cell.annot["PC-3_JPC-3", "unique.tissueid"] <- "lung"
cell.annot["PL4", "unique.tissueid"] <- "pancreas"
cell.annot["Sarc9371", "unique.tissueid"] <- "bone"
cell.annot["U-CH2", "unique.tissueid"] <- "placenta"
cell.annot["UDSCC2", "unique.tissueid"] <- "upper_aerodigestive_tract"
cell.annot["WM793B", "unique.tissueid"] <- "skin"

## save consolidated annotation file
tt <- cbind("unique.cellid"=rownames(cell.annot), cell.annot)
tt[is.na(tt)] <- ""
write.csv(tt, row.names=FALSE, file=file.path(path.out, "cell_annotation_all.csv"))


# myx <- intersect(rownames(cura), cosmic[ , "cellid"])
# cura[myx, "COSMIC.tissueid"] <- cosmic[myx, "tissueid"]
# cura[is.na(cura)] <- ""
# write.csv(cura, "matching_cell.csv")

## drugs
drug.annot <- NULL
cura <- read.csv("matching_drug.csv", na.strings=c("NA","NaN", " ",""))
#cura[!is.na(cura) & cura == ""] <- NA
rownames(cura) <- cura[ , "unique.drugid"]
## CGP
tt <- read.csv(file.path(path.out, "drug_annotation_CGP.csv"))
tt <- subset(tt, !is.na(tt[ , "drug.name"])&!duplicated(tt[ , "drug.name"]))
rownames(tt) <- tt[ , "drug.name"]
curat <- cura[!is.na(cura[ , "CGP.drugid"]), , drop=FALSE]
rownames(tt)[match(curat[ , "CGP.drugid"], rownames(tt))] <- rownames(curat)
drug.annot <- cbind(drug.annot, "CGP.drugid"=tt[ , "drug.name"])
rownames(drug.annot) <- rownames(tt)
## CCLE
tt <- read.csv(file.path(path.out, "drug_annotation_CCLE.csv"))
tt <- subset(tt, !is.na(tt[ , "drug.name"]))
rownames(tt) <- tt[ , "drug.name"]
curat <- cura[!is.na(cura[ , "CCLE.drugid"]), , drop=FALSE]
rownames(tt)[match(curat[ , "CCLE.drugid"], rownames(tt))] <- rownames(curat)
myx2 <- setdiff(rownames(tt), rownames(drug.annot))
drug.annot <- rbind(drug.annot, matrix(NA, nrow=length(myx2), ncol=ncol(drug.annot), dimnames=list(myx2, colnames(drug.annot))))
drug.annot <- cbind(drug.annot, "CCLE.drugid"=NA)
drug.annot[rownames(tt), "CCLE.drugid"] <- tt[ , "drug.name"]
## GSK
tt <- read.csv(file.path(path.out, "drug_annotation_GSK.csv"))
#tt <- subset(tt, !is.na(tt[ , "drug.name"]))
rownames(tt) <- tt[ , "drug.name"]
curat <- cura[!is.na(cura[ , "GSK.drugid"]), , drop=FALSE]
rownames(tt)[match(curat[ , "GSK.drugid"], rownames(tt))] <- rownames(curat)
myx2 <- setdiff(rownames(tt), rownames(drug.annot))
drug.annot <- rbind(drug.annot, matrix(NA, nrow=length(myx2), ncol=ncol(drug.annot), dimnames=list(myx2, colnames(drug.annot))))
drug.annot <- cbind(drug.annot, "GSK.drugid"=NA)
drug.annot[rownames(tt), "GSK.drugid"] <- tt[ , "drug.name"]
## NCI60
tt <- read.csv(file.path(path.out, "drug_annotation_NCI60.csv"))
#tt <- subset(tt, !is.na(tt[ ,"drug.name"]) & !duplicated(tt[,"drug.name"]))
#rownames(tt) <- tt[ , "drug.name"]
rownames(tt) <- tt[,"drugid"]
curat <- cura[!is.na(cura[ , "NCI60.drugid"]), , drop=FALSE]
rownames(tt)[match(curat[ , "NCI60.drugid"], rownames(tt))] <- rownames(curat)
myx2 <- setdiff(rownames(tt), rownames(drug.annot))
drug.annot <- rbind(drug.annot, matrix(NA, nrow=length(myx2), ncol=ncol(drug.annot), dimnames=list(myx2, colnames(drug.annot))))
drug.annot <- cbind(drug.annot, "NCI60.drugid"=NA)
drug.annot[rownames(tt), "NCI60.drugid"] <- tt[ , "drug.name"]
## GRAY
tt <- read.csv(file.path(path.out, "drug_annotation_GRAY.csv"))
#tt <- subset(tt, !is.na(tt[ , "drug.name"]))
rownames(tt) <- tt[ , "drug.name"]
curat <- cura[!is.na(cura[ , "GRAY.drugid"]), , drop=FALSE]
rownames(tt)[match(curat[ , "GRAY.drugid"], rownames(tt))] <- rownames(curat)
myx2 <- setdiff(rownames(tt), rownames(drug.annot))
drug.annot <- rbind(drug.annot, matrix(NA, nrow=length(myx2), ncol=ncol(drug.annot), dimnames=list(myx2, colnames(drug.annot))))
drug.annot <- cbind(drug.annot, "GRAY.drugid"=NA)
drug.annot[rownames(tt), "GRAY.drugid"] <- tt[ , "drug.name"]
## GNE
tt <- read.csv(file.path(path.out, "drug_annotation_GNE.csv"))
#tt <- subset(tt, !is.na(tt[ , "drug.name"]))
rownames(tt) <- tt[ , "drug.name"]
curat <- cura[!is.na(cura[ , "GNE.drugid"]), , drop=FALSE]
rownames(tt)[match(curat[ , "GNE.drugid"], rownames(tt))] <- rownames(curat)
myx2 <- setdiff(rownames(tt), rownames(drug.annot))
drug.annot <- rbind(drug.annot, matrix(NA, nrow=length(myx2), ncol=ncol(drug.annot), dimnames=list(myx2, colnames(drug.annot))))
drug.annot <- cbind(drug.annot, "GNE.drugid"=NA)
drug.annot[rownames(tt), "GNE.drugid"] <- tt[ , "drug.name"]
## CTRP
tt <- read.csv(file.path(path.out, "drug_annotation_CTRP.csv"))
#tt <- subset(tt, !is.na(tt[ , "drug.name"]))
rownames(tt) <- tt[ , "drug.name"]
curat <- cura[!is.na(cura[ , "CTRP.drugid"]), , drop=FALSE]
rownames(tt)[match(curat[ , "CTRP.drugid"], rownames(tt))] <- rownames(curat)
myx2 <- setdiff(rownames(tt), rownames(drug.annot))
drug.annot <- rbind(drug.annot, matrix(NA, nrow=length(myx2), ncol=ncol(drug.annot), dimnames=list(myx2, colnames(drug.annot))))
drug.annot <- cbind(drug.annot, "CTRP.drugid"=NA)
drug.annot[rownames(tt), "CTRP.drugid"] <- tt[ , "drug.name"]

## CMAP
tt <- read.csv(file.path(path.out, "cmap_CBID_Lamb.csv"))
#tt <- subset(tt, !is.na(tt[ , "drug.name"]))
rownames(tt) <- tt[ , "cmap_name"]
curat <- cura[!is.na(cura[ , "CMAP.drugid"]), , drop=FALSE]
rownames(tt)[match(curat[ , "CMAP.drugid"], rownames(tt))] <- rownames(curat)
myx2 <- setdiff(rownames(tt), rownames(drug.annot))
drug.annot <- rbind(drug.annot, matrix(NA, nrow=length(myx2), ncol=ncol(drug.annot), dimnames=list(myx2, colnames(drug.annot))))
drug.annot <- cbind(drug.annot, "CMAP.drugid"=NA)
drug.annot[rownames(tt), "CMAP.drugid"] <- tt[ , "cmap_name"]



tt <- cbind("unique.drugid"=rownames(drug.annot), drug.annot)
tt[is.na(tt)] <- ""
write.csv(tt, row.names=FALSE, file=file.path(path.out, "drug_annotation_all.csv"))



##tissues
drug.annot <- NULL
cura <- read.csv("matching_tissue.csv", na.strings=c("NA","NaN", " ",""))
#cura[!is.na(cura) & cura == ""] <- NA
rownames(cura) <- cura[ , "COSMIC.tissueid"]
## CGP
tt <- read.csv(file.path(path.out, "cell_line_annotation_CGP.csv"))
drug.annot <- cbind(drug.annot, "CGP.drugid"=tt[ , "drug.name"])
rownames(drug.annot) <- tt[ , "drug.name"]





