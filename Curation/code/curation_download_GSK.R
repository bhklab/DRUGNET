## path to files
path.data <- file.path("data", "GSK")
path.drug <- file.path(path.data, "drug")
path.cell <- file.path(path.data, "celline")
## create directories
if(!file.exists(path.data)) { dir.create(path.data, showWarnings=FALSE, recursive=TRUE) }
if(!file.exists(path.drug)) { dir.create(path.drug, showWarnings=FALSE, recursive=TRUE) }
if(!file.exists(path.cell)) { dir.create(path.cell, showWarnings=FALSE, recursive=TRUE) }


ftpdir <- "ftp://caftpd.nci.nih.gov/pub/caARRAY/transcript_profiling"

## download sample information
message("Download sample information")
myfn <- file.path(path.cell, "tmp", "gsk_ge_sampleinfo.txt")
if (!file.exists(myfn)) {
  dir.create(file.path(path.cell, "tmp"), showWarnings=FALSE, recursive=TRUE)
  dwl.status <- download.file(url=sprintf("%s/GSK_RNA.sdrf", ftpdir), destfile=file.path(path.cell, "tmp", "GSK_RNA.sdrf"), quiet=TRUE)
  if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
  file.copy(from=file.path(path.cell, "tmp", "GSK_RNA.sdrf"), to=myfn)
}
  
## download drug sensitivity
message("Download drug sensitivity measurements")
myfn <- file.path(path.drug, "tmp", "gsk_drug_sensitivity.xls")
if (!file.exists(myfn)) {
  dir.create(file.path(path.drug, "tmp"), showWarnings=FALSE)
  dwl.status <- download.file(url="http://cancerres.aacrjournals.org/content/suppl/2010/04/19/0008-5472.CAN-09-3788.DC1/stab_2.xls", destfile=file.path(path.drug, "tmp", "stab_2.xls"), quiet=TRUE)
  if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
  file.copy(from=file.path(path.drug, "tmp", "stab_2.xls"), to=myfn)
}

## info about each experiment
message("Read sample information")
sampleinfo <- read.csv(file.path(path.cell, "tmp", "gsk_ge_sampleinfo.txt"), sep="\t", comment.char="#")
sampleinfo[sampleinfo == "" | sampleinfo == " "] <- NA
## curate cell line names
rownames(sampleinfo) <- gsub(" - Replicate ", "_rep", sampleinfo[ , "Source.Name"])
sampleinfo <- data.frame("samplename"=rownames(sampleinfo), "cellid"=sampleinfo[ , "Characteristics.Cell.Line.Name."], "filename"=sprintf("%s.gz", sampleinfo[ , "Array.Data.File"]), "tissueid"=tolower(gsub("_$", "", gsub(badchars, "_", sampleinfo[ , "Characteristics.OrganismPart."]))), sampleinfo)
## manual curation of cell line names
sampleinfo[!is.na(sampleinfo[ , "cellid"]) & sampleinfo[ , "cellid"] == "Wi38", "cellid"] <- "WI38"
sampleinfo[!is.na(sampleinfo[ , "cellid"]) & sampleinfo[ , "cellid"] == "DoTc2", "cellid"] <- "DOTC24510"
sampleinfo[!is.na(sampleinfo[ , "cellid"]) & sampleinfo[ , "cellid"] == "MV4II", "cellid"] <- "MV411"
sampleinfo[!is.na(sampleinfo[ , "cellid"]) & sampleinfo[ , "cellid"] == "SKNFI", "cellid"] <- "SKNF1"

## remove replicates
sampleinfo <- sampleinfo[!duplicated(sampleinfo[ , "cellid"]), , drop=FALSE]


## read drug phenotypes
drugpheno <- gdata::read.xls(xls=file.path(path.drug, "tmp", "gsk_drug_sensitivity.xls"), sheet=1)
drugpheno[drugpheno == "" | drugpheno == "*"] <- NA
cn <- gsub(badchars, ".", drugpheno[6, ])
cn <- gsub("CL.ID", "CL_ID", cn)
drugpheno <- drugpheno[-(1:6), , drop=FALSE]
drugpheno <- drugpheno[!is.na(drugpheno[ , 2]), , drop=FALSE]
dimnames(drugpheno) <- list(drugpheno[ , 2], cn)
celline.info <- drugpheno[ , c("CL_ID", "Cell.Line", "Site", "Dx"), drop=FALSE]
celline.info <- data.frame(celline.info, "cellid"=as.character(celline.info[ , "Cell.Line"]), "tissueid"=as.character(celline.info[ , "Site"]))
drugpheno <- data.matrix(drugpheno[ , grep("IC50", cn), drop=FALSE])
## IC50 in nano molar

## cell line identifiers are different between gene expression experiments and sensitivity screening
matchix <- match(toupper(gsub(badchars, "", sampleinfo[ , "cellid"])), toupper(gsub(badchars, "", celline.info[ , "cellid"])))
ix0 <- which(!is.na(matchix))
ix1 <- which(is.na(matchix))
ix2 <- setdiff(1:nrow(celline.info), matchix[!is.na(matchix)])
dd <- data.frame(matrix(NA, ncol=4, nrow=nrow(sampleinfo) + length(ix2), dimnames=list(c(toupper(gsub(badchars, "", sampleinfo[ , "cellid"])), toupper(gsub(badchars, "", celline.info[ix2, "cellid"]))), c("sampleinfo.cellid", "sampleinfo.tissueid", "drugpheno.cellid", "drugpheno.tissueid"))))
dd[toupper(gsub(badchars, "", sampleinfo[ , "cellid"])), c("sampleinfo.cellid", "sampleinfo.tissueid")] <- sampleinfo[ , c("cellid", "tissueid")]
dd[ix0, c("drugpheno.cellid", "drugpheno.tissueid")] <- celline.info[matchix[!is.na(matchix)], c("cellid", "tissueid")]
dd[toupper(gsub(badchars, "", celline.info[ix2, "cellid"])), c("drugpheno.cellid", "drugpheno.tissueid")] <- celline.info[ix2, c("cellid", "tissueid")]
write.csv(dd, file=file.path(path.cell, "tmp", "gsk_cell_line_matching.csv"))

rownames(drugpheno) <- toupper(gsub(badchars, "", celline.info[ , "cellid"]))
rownames(celline.info) <- toupper(gsub(badchars, "", celline.info[ , "cellid"]))
rownames(sampleinfo) <- toupper(gsub(badchars, "", sampleinfo[ , "cellid"]))


## druginfo
xx <- gsub("[(]gIC50..nM[)]", "", colnames(drugpheno))
xx <- gsub("[.]*$", "", xx)
druginfo <- data.frame("drugid"=toupper(xx), "drug.name"=xx)
rownames(druginfo) <- toupper(xx)
myx <- match(c("drugid", "drug.name"), colnames(druginfo))
druginfo <- druginfo[ , c(myx, setdiff(1:ncol(druginfo), myx)), drop=FALSE]

## cell lines and drugs
cellnall <- sort(unique(c(rownames(sampleinfo), rownames(drugpheno))))
drugnall <- sort(unique(druginfo[ , "drugid"]))

## update sampleinfo
dd <- data.frame(matrix(NA, nrow=length(cellnall), ncol=ncol(sampleinfo), dimnames=list(cellnall, colnames(sampleinfo))))
newlev <- sapply(sampleinfo, levels)
newlev$cellid <- cellnall
dd <- genefu::setcolclass.df(df=dd, colclass=sapply(sampleinfo, class), factor.levels=newlev)
dd[rownames(sampleinfo), colnames(sampleinfo)] <- sampleinfo
dd[ , "cellid"] <- rownames(dd)
dd[is.na(dd[ , "tissueid"]), "tissueid"] <- tolower(gsub(badchars, "_", celline.info[rownames(dd)[is.na(dd[ , "tissueid"])], "tissueid"]))
sampleinfo <- dd


write.csv(sampleinfo, file=file.path(path.out, "cell_line_annotation_GSK.csv"), row.names=FALSE)
write.csv(druginfo, file=file.path(path.out, "drug_annotation_GSK.csv"), row.names=FALSE)

