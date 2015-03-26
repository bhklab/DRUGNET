## path to files
path.data <- file.path("data", "CGP")
path.drug <- file.path(path.data, "drug")
path.cell <- file.path(path.data, "celline")
## create directories
if(!file.exists(path.data)) { dir.create(path.data, showWarnings=TRUE, recursive=TRUE) }
if(!file.exists(path.drug)) { dir.create(path.drug, showWarnings=TRUE, recursive=TRUE) }
if(!file.exists(path.cell)) { dir.create(path.cell, showWarnings=TRUE, recursive=TRUE) }

########################

## download sample information
message("Download sample information")
ftpdir <- "ftp://ftp.ebi.ac.uk//pub/databases/microarray/data/experiment/MTAB/E-MTAB-783/"
myfn <- file.path(path.cell, "tmp", "cgp_ge_sampleinfo.txt")
if (!file.exists(myfn)) {
  dir.create(file.path(path.cell, "tmp"), showWarnings=FALSE, recursive=TRUE)
  dwl.status <- download.file(url=sprintf("%s/E-MTAB-783.sdrf.txt", ftpdir), destfile=file.path(path.cell, "tmp", "E-MTAB-783.sdrf.txt"), quiet=TRUE)
  if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
  file.copy(from=file.path(path.cell, "tmp", "E-MTAB-783.sdrf.txt"), to=myfn)
}
  
## download drug sensitivity
message("Download drug sensitivity measurements")
myfn <- file.path(path.drug, "tmp", "cgp_drug_sensitivity.csv")
if (!file.exists(myfn)) {
  dir.create(file.path(path.drug, "tmp"), showWarnings=FALSE, recursive=TRUE)
  dwl.status <- download.file(url="ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-5.0/gdsc_manova_input_w5.csv", destfile=file.path(path.drug, "tmp", "gdsc_manova_input_w5.csv"), quiet=TRUE)
  if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
  file.copy(from=file.path(path.drug, "tmp", "gdsc_manova_input_w5.csv"), to=myfn)
}

## download drug concentration
message("Download screening drug concentrations")
myfn <- file.path(path.drug, "tmp", "cgp_drug_concentration.csv")
if (!file.exists(myfn)) {
  dir.create(file.path(path.drug, "tmp"), showWarnings=FALSE, recursive=TRUE)
  dwl.status <- download.file(url="ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-5.0/gdsc_compounds_conc_w5.csv", destfile=file.path(path.drug, "tmp", "gdsc_compounds_conc_w5.csv"), quiet=TRUE)
  if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
  file.copy(from=file.path(path.drug, "tmp", "gdsc_compounds_conc_w5.csv"), to=myfn)
}

## download cell line annotations and COSMIC IDs
## annotations from COSMIC cell line project
myfn <- file.path(path.cell, "tmp", "cosmic_annotations.RData")
if(!file.exists(myfn)) {
  message("Download COSMIC annotations for cell lines")
  myfn2 <- file.path(path.cell, "tmp", "cosmic_cell_line_collection.txt")
  if(!file.exists(myfn2)) {
    dir.create(file.path(path.cell, "tmp"), showWarnings=FALSE, recursive=TRUE)
    dwl.status <- getCosmic(em="bhk.labgroup@gmail.com", passw="pharmacogenomics", directory=file.path(path.cell, "tmp"))
    # dwl.status <- download.file(url=sprintf("http://cancer.sanger.ac.uk/files/cosmic/current_release/CosmicCompleteExport.tsv.gz"), destfile=file.path(path.cell, "tmp", sprintf("CosmicCompleteExport.tsv.gz")), quiet=TRUE)
    if(dwl.status != 0) { stop("Download failed, please rerun the pipeline") }
    ## untar
    res <- R.utils::gunzip(filename=file.path(path.cell, "tmp", sprintf("CosmicCompleteExport.tsv.gz")), overwrite=TRUE)
    file.copy(from=file.path(path.cell, "tmp", "CosmicCompleteExport.tsv"), to=myfn2)
  }
  message("Process COSMIC annotations")
  cosmic.celline <- read.csv(file=file.path(path.cell, "tmp", "cosmic_cell_line_collection.txt"), sep="\t")
  # cosmic.celline <- cosmic.celline[- c(grep("row selected", cosmic.celline[ ,1]), grep("rows selected", cosmic.celline[ ,1])), , drop=FALSE]
  cosmic.celline <- cosmic.celline[complete.cases(cosmic.celline[ , c("Sample.name", "Sample.source")]) & cosmic.celline[ , "Sample.source"] == "cell-line", , drop=FALSE]
  cosmic.celline[cosmic.celline == "NS" | cosmic.celline == "" | cosmic.celline == " " | cosmic.celline == "  "] <- NA
  ## merge the gene targets
  dupln <- sort(unique(cosmic.celline[ , "Sample.name"][duplicated(cosmic.celline[ , "Sample.name"])]))
  tt <- cosmic.celline
  ## select unique cell lines
  iix.rm <- NULL
  for(i in 1:length(dupln)) {
    duplix <- cosmic.celline[ ,"Sample.name"] == dupln[i]
    iix <- sort((which(duplix)), decreasing=FALSE)[1]
    iix.rm <- c(iix.rm, setdiff(which(duplix), iix))
    ## get the most frequent tissue type
    tissuet <- table(cosmic.celline[duplix, "Primary.site"])
    if (length(tissuet) == 0) {
      tt[iix, "Primary.site"] <- NA
    } else {
      tt[iix, "Primary.site"] <- names(sort(tissuet, decreasing=TRUE))[1]
    }
    # tt[iix, "Gene.name"] <- paste(cosmic.celline[duplix, "Gene.name"], collapse="///")
    # tt[iix, "UniProt.ID"] <- paste(cosmic.celline[duplix, "UniProt.ID"], collapse="///")
    # tt[iix, "Zygosity"] <- paste(cosmic.celline[duplix, "Zygosity"], collapse="///")
    # tt[iix, "CDS_MUT_SYNTAX"] <- paste(cosmic.celline[duplix, "CDS_MUT_SYNTAX"], collapse="///")
    # tt[iix, "AA_MUT_SYNTAX"] <- paste(cosmic.celline[duplix, "AA_MUT_SYNTAX"], collapse="///")
    # tt[iix, "NCBI36.genome.position"] <- paste(cosmic.celline[duplix, "NCBI36.genome.position"], collapse="///")
    # tt[iix, "GRCh37.genome.position"] <- paste(cosmic.celline[duplix, "GRCh37.genome.position"], collapse="///")
  }
  tt <- tt[-iix.rm, , drop=FALSE]
  tt <- tt[!is.na(tt[ , "Sample.name"]), , drop=FALSE]
  rownames(tt) <- tt[ , "Sample.name"]
  ## remove unnecessary annotations
  tt <- tt[ , c("Sample.name", "ID_sample", "ID_tumour", "Primary.site", "Site.subtype", "Primary.histology", "Histology.subtype", "Sample.source", "Tumour.origin", "Comments"), drop=FALSE]
  cosmic.celline <- tt
  save(list=c("cosmic.celline"), compress=TRUE, file=myfn)
} else { load(myfn) }

## annotations from GDSC (Genomics of Drug Sensitivity in Cancer)
myfn <- file.path(path.cell, "tmp", "gdsc_annotations.RData")
if(!file.exists(myfn)) {
  message("Download GDSC annotations for cell liness")
  myfn2 <- file.path(path.cell, "tmp", "cgp_celline_collection.csv")
  if(!file.exists(myfn2)) {
    dir.create(file.path(path.cell, "tmp"), showWarnings=FALSE, recursive=TRUE)
    dwl.status <- download.file(url="ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-5.0/gdsc_cell_lines_w5.csv", destfile=file.path(path.cell, "tmp", "gdsc_cell_lines_w5.csv"), quiet=TRUE)
    if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
    file.copy(from=file.path(path.cell, "tmp", "gdsc_cell_lines_w5.csv"), to=myfn2)
  }
  gdsc.celline <- read.csv(file=file.path(path.cell, "tmp", "cgp_celline_collection.csv"))
  gdsc.celline[gdsc.celline == "" | gdsc.celline == " " | gdsc.celline == "  "] <- NA
  gdsc.celline <- gdsc.celline[!is.na(gdsc.celline[ , "CELL_LINE_NAME"]), , drop=FALSE]
  dupln <- unique(gdsc.celline[ , "CELL_LINE_NAME"][duplicated(gdsc.celline[ , "CELL_LINE_NAME"])])
  gdsc.celline <- gdsc.celline[!duplicated(gdsc.celline[ , "CELL_LINE_NAME"]), , drop=FALSE]
  rownames(gdsc.celline) <- gdsc.celline[ , "CELL_LINE_NAME"]
  save(list=c("gdsc.celline"), compress=TRUE, file=myfn)
} else { load(myfn) }

## merge GDSC and COSMIC annotations through COSMIC_ID
message("Merge COSMIC and GDSC annotations for cell liness")
iix <- which(complete.cases(gdsc.celline[ , c("CELL_LINE_NAME", "COSMIC_ID")]) & !is.element(gdsc.celline[ , "COSMIC_ID"], cosmic.celline[ , "ID_sample"]) & !is.element(gdsc.celline[ , "CELL_LINE_NAME"], cosmic.celline[ , "Sample.name"]))
tt <- data.frame(matrix(NA, nrow=nrow(cosmic.celline) + length(iix), ncol=ncol(cosmic.celline), dimnames=list(c(rownames(cosmic.celline), rownames(gdsc.celline)[iix]), colnames(cosmic.celline))))
tt[rownames(cosmic.celline), ] <- cosmic.celline
tt[rownames(gdsc.celline)[iix], "Sample.name"] <- gdsc.celline[iix, "CELL_LINE_NAME"]
tt[rownames(gdsc.celline)[iix], "ID_sample"] <- gdsc.celline[iix, "COSMIC_ID"]
celline.cgp <- tt

## download drug information
message("Download drug information")
myfn <- file.path(path.drug, "tmp",  "cgp_drug_information.csv")
if (!file.exists(myfn)) {
  dir.create(file.path(path.drug, "tmp"), showWarnings=FALSE, recursive=TRUE)
  # dwl.status <- download.file(url="http://www.cancerrxgene.org/action/ExportJsonTable/CSV", destfile=file.path(path.drug, "tmp", "export-Automatically_generated_table_data.csv"), quiet=TRUE)
  # if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }  
  tables <- XML::readHTMLTable("http://www.cancerrxgene.org/translation/Drug")
  drugs <- tables[1][[1]]
  write.csv(drugs, row.names=FALSE, file=file.path(path.drug, "tmp", "export.csv"))
  file.copy(from=file.path(path.drug, "tmp", "export.csv"), to=myfn)
}
myfn <- file.path(path.drug, "tmp", "nature_supplementary_information.xls")
if (!file.exists(myfn)) {
  dir.create(file.path(path.drug, "tmp"), showWarnings=FALSE, recursive=TRUE)
  dwl.status <- download.file(url="http://www.nature.com/nature/journal/v483/n7391/extref/nature11005-s2.zip", destfile=file.path(path.drug, "tmp", "nature11005-s2.zip"), quiet=TRUE)
  ff <- as.character(unzip(zipfile=file.path(path.drug, "tmp", "nature11005-s2.zip"), list=TRUE)[1, 1])
  unzip(zipfile=file.path(path.drug, "tmp", "nature11005-s2.zip"), exdir=file.path(path.drug, "tmp"))
  file.copy(from=file.path(path.drug, "tmp", ff), to=myfn)
}


message("Read drug sensitivity measurements")
myfn2 <- file.path(path.drug, "tmp", "cgp_drug_sensitivity.RData")
if(!file.exists(myfn2)) {
  drugpheno <- read.csv(file.path(path.drug, "tmp", "cgp_drug_sensitivity.csv"))
  drugpheno[drugpheno == "" | drugpheno == " "] <- NA
  save(list="drugpheno", compress=TRUE, file=myfn2)
} else { load(myfn2) }
## format column names
coln2 <- unlist(drugpheno[1, ,drop=TRUE])
coln2[coln2 == ""] <- NA
drugpheno <- drugpheno[!is.na(drugpheno[ , "Cell.Line"]), ,drop=FALSE]
coln <- colnames(drugpheno)
coln2[is.na(coln2)] <- coln[is.na(coln2)]
coln2 <- genefu::rename.duplicate(x=coln2, sep="_dupl")$new.x
myx <- sapply(sapply(strsplit(coln2, "_"), function(x) { return(x[[1]]) }), Hmisc::all.is.numeric)
coln2[myx] <- paste("drugid", gsub(pattern=badchars, replacement="_", x=toupper(coln2[myx])), sep="_")
colnames(drugpheno) <- coln2
## drug identifiers and names
dn <- toupper(gsub(badchars, "", sapply(strsplit(coln, "_"), function(x) { return(x[[1]]) })))
## manual curation for drug names starting with a figure
dn[!is.na(dn) & dn == "X17AAG"] <- "17AAG"
dn[!is.na(dn) & dn == "X681640"] <- "681640"
did <- sapply(strsplit(coln2, "_"), function(x) { if(x[[1]] == "drugid") { return(x[[2]]) } else { return(NA) } })
drugnid <- cbind("drug.name2"=dn, "drug.id"=did)[!is.na(did) & !duplicated(did), ]
rownames(drugnid) <- paste("drugid", drugnid[ , "drug.id"], sep="_")

## cell line identifiers
dupln <- duplicated(drugpheno[ ,"Cell.Line"])
if(sum(dupln) > 1) { warning("some cell lines are duplicated, only the first instance is kept") }
drugpheno <- drugpheno[!dupln, , drop=FALSE]
if(any(!is.element(drugpheno[ ,"Cell.Line"], celline.cgp[ , "Sample.name"]))) { stop("Some cell line names are not included in the COSMIC database") }
celln <- drugpheno[ ,"Cell.Line"]
drugpheno <- data.frame("cellid"=celln, drugpheno)
rownames(drugpheno) <- celln


## info about each experiment
message("Read sample information")
sampleinfo <- read.csv(file.path(path.cell, "tmp", "cgp_ge_sampleinfo.txt"), sep="\t")
sampleinfo[sampleinfo == "" | sampleinfo == " "] <- NA
## curate cell line names
sampleinfo[sampleinfo[ , "Source.Name"] == "MZ2-MEL.", "Source.Name"] <- "MZ2-MEL"
iix <- which(!duplicated(sampleinfo[ , "Source.Name"]) & !is.element(sampleinfo[ , "Source.Name"], celline.cgp[ , "Sample.name"]))
if(length(iix) > 0) {
  ## enrich the list of cell lines
  tt <- matrix(NA, nrow=length(iix), ncol=ncol(celline.cgp), dimnames=list(sampleinfo[iix, "Source.Name"], colnames(celline.cgp)))
  tt[ , "Sample.name"] <- sampleinfo[iix, "Source.Name"]
  celline.cgp <- rbind(celline.cgp, tt)
}
fn <- gsub(patter="[.]CEL", replacement="", x=sampleinfo[ ,"Array.Data.File"])
rownames(sampleinfo) <- fn
sampleinfo <- data.frame("cellid"=sampleinfo[ , "Source.Name"], sampleinfo)
## remove duplcated cell line hybridization
dupln <- duplicated(sampleinfo[ ,"cellid"])
if(sum(dupln) > 1) { warning("some cell lines have been profiled for gene expression multiple times, only the first instance is kept") }
sampleinfo <- sampleinfo[!dupln, , drop=FALSE]
rownames(sampleinfo) <- sampleinfo[ ,"cellid"]

## update of cgp cell line collection
celline.cgp <- data.frame("cellid"=as.character(celline.cgp[ , "Sample.name"]), celline.cgp)
celline.cgp[ , "cellid"] <- as.character(celline.cgp[ , "cellid"])
## add url based on COSMIC IDs
uurl <- paste("http://cancer.sanger.ac.uk/cell_lines/sample/overview?id=", celline.cgp[ , "ID_sample"], sep="")
uurl[is.na(celline.cgp[ , "ID_sample"])] <- NA
celline.cgp <- data.frame("cellid"=celline.cgp[ , "cellid"], "tissueid"=celline.cgp[ , "Primary.site"], "link"=uurl, celline.cgp[ , !is.element(colnames(celline.cgp), "cellid")])

## drug information
message("Read drug information")
druginfo <- read.csv(file.path(path.drug, "tmp", "cgp_drug_information.csv"))
druginfo[!is.na(druginfo) & (druginfo == " " | druginfo == " ")] <- NA
druginfo <- data.frame("drug.name2"=toupper(gsub(badchars, "", druginfo[ ,"Name"])), "drug.name"=druginfo[ ,"Name"], druginfo)
myx <- match(druginfo[ , "drug.name2"], drugnid[ , "drug.name2"])
if (any(is.na(myx))) { stop ("Some drugs have missing annotations") }
## correct ambiguity for AZD6482: drugid_156 corresponds to the first occurence of AZD6482 while drugid_1066 corresponds to the second
## table(!is.na(drugpheno[ , "drugid_156_AUC"]))
## table(!is.na(drugpheno[ , "drugid_1066_AUC"]))
myx[druginfo[ , "drug.name2"] == "AZD6482"][2] <- which(drugnid[ , "drug.name2"] == "AZD6482")[2]
druginfo <- data.frame("drugid"=rownames(drugnid)[myx], drugnid[myx, , drop=FALSE], druginfo)
rownames(druginfo) <- as.character(druginfo[ ,"drugid"])
## complement drug infomration with the supplementary infomration from the Nature website
myfn2 <- file.path(path.drug, "tmp", "nature_supplinfo_druginfo_cgp.RData")
if(!file.exists(myfn2)) {
  druginfo.nature <- gdata::read.xls(xls=file.path(path.drug, "tmp", "nature_supplementary_information.xls"), sheet=4)
  druginfo.nature[druginfo.nature == "" | druginfo.nature == " "] <- NA
  save(list="druginfo.nature", compress=TRUE, file=myfn2)
} else { load(myfn2) }
rownames(druginfo.nature) <- paste("drugid", druginfo.nature[ , "Drug.ID"], sep="_")
druginfo <- data.frame(druginfo, druginfo.nature[rownames(druginfo), c("Brand.name", "Site.of.screening", "Drug.type", "Drug.class.I", "Drug.class.II", "Target.family", "Effector.pathway.biological.process", "Clinical.trials", "Source")])

## drug concentration
message("Read drug concentration")
drugconc <- read.csv(file.path(path.drug, "tmp", "cgp_drug_concentration.csv"))
drugconc[!is.na(drugconc) & (drugconc == "" | drugconc == " ")] <- NA
drugconc <- data.frame("drug.name2"=toupper(gsub(badchars, "", drugconc[ ,"Compound.Name"])), drugconc)
if(all(!is.element(drugconc[ , "drug.name2"], drugnid[ , "drug.name2"]))) { stop("Screening concentration for drugs without identifiers!") }
myx <- match(drugconc[ , "drug.name2"], drugnid[ , "drug.name2"])
## correct ambiguity for AZD6482: drugid_156 corresponds to the first occurence of AZD6482 while drugid_1066 corresponds to the second
myx[drugconc[ , "drug.name2"] == "AZD6482"][2] <- which(drugnid[ , "drug.name2"] == "AZD6482")[2]
rownames(drugconc) <- rownames(drugnid)[myx]
drugconc <- data.frame("drugid"=rownames(drugconc), drugconc)

## combne all cell lines
cellnall <- sort(unique(c(as.character(sampleinfo[ ,"cellid"]), as.character(drugpheno[ ,"cellid"]))))

## combine all drugs
dix <- sort(unique(c(rownames(druginfo), rownames(drugconc), paste("drugid", sapply(strsplit(colnames(drugpheno)[grep("^drugid_", colnames(drugpheno))], "_"), function(x) { return(x[[2]]) }), sep="_"))))
## update druginfo
druginfo2 <- data.frame(matrix(NA, nrow=length(dix), ncol=ncol(druginfo), dimnames=list(dix, colnames(druginfo))))
newlev <- sapply(druginfo, levels)
newlev$drugid <- sapply(strsplit(dix, split="_"), function(x) { return(x[2]) })
druginfo2 <- genefu::setcolclass.df(df=druginfo2, colclass=sapply(druginfo, class), factor.levels=newlev)
druginfo2[match(rownames(druginfo), dix), colnames(druginfo)] <- druginfo
druginfo2[ , "drugid"] <- newlev$drugid
druginfo <- druginfo2
## update drugconc
drugconc2 <- data.frame(matrix(NA, nrow=length(dix), ncol=ncol(drugconc), dimnames=list(dix, colnames(drugconc))))
newlev <- sapply(drugconc, levels)
newlev$drugid <- sapply(strsplit(dix, split="_"), function(x) { return(x[2]) })
drugconc2 <- genefu::setcolclass.df(df=drugconc2, colclass=sapply(drugconc, class), factor.levels=newlev)
drugconc2[match(rownames(drugconc), dix), colnames(drugconc)] <- drugconc
drugconc2[ , "drugid"] <- newlev$drugid
drugconc <- drugconc2

## report concentrations per cell line and per drug
drugconc2 <- data.frame(matrix(NA, nrow=nrow(drugconc) * length(cellnall), ncol=6, dimnames=list(paste(rep(rownames(drugconc), times=length(cellnall)), rep(cellnall, each=nrow(drugconc)), sep="..."), c("cellid", "drugid", "drug.name2", "nbr.conc.tested", "min.Dose.uM", "max.Dose.uM"))))
drugconc2[ , "cellid"] <- rep(cellnall, times=nrow(drugconc))
drugconc2[ , "drugid"] <- rep(rownames(drugconc), each=length(cellnall))
drugconc2[ , "drug.name2"] <- rep(as.character(drugconc[ ,"drug.name2"]), each=length(cellnall))
## as mentioned in the supplementary information of Garnett et al., a single cell line is used on each plate and treated with 28 different drugs over a 9-pt, 256-fold concentration range
drugconc2[ , "nbr.conc.tested"] <- 9
drugconc2[ , "min.Dose.uM"] <- rep(drugconc[ , "Min.Concentration.micromolar."], each=length(cellnall))
drugconc2[ , "max.Dose.uM"] <- rep(drugconc[ , "Max.Concentration.micromolar."], each=length(cellnall))
drugconc <- drugconc2

## combine al cell lines and drugs
cellnall <- sort(unique(c(as.character(sampleinfo[ ,"cellid"]), as.character(drugpheno[ ,"cellid"]))))
drugnall <- sort(unique(c(rownames(druginfo), as.character(drugconc[ , "drugid"]), paste("drugid", sapply(strsplit(colnames(drugpheno)[grep("^drugid_", colnames(drugpheno))], "_"), function(x) { return(x[[2]]) }), sep="_"))))

## update drugpheno
dd <- data.frame(matrix(NA, ncol=ncol(drugpheno), nrow=length(cellnall), dimnames=list(cellnall, colnames(drugpheno))))
newlev <- sapply(drugpheno, levels)
newlev$cellid <- cellnall
dd <- genefu::setcolclass.df(df=dd, colclass=sapply(drugpheno, class), factor.levels=newlev)
dd[rownames(drugpheno),colnames(drugpheno)] <- drugpheno
dd[ ,"cellid"] <- cellnall
drugpheno <- dd

## update sampleinfo
dd <- data.frame(matrix(NA, nrow=length(cellnall), ncol=ncol(sampleinfo), dimnames=list(cellnall, colnames(sampleinfo))))
newlev <- sapply(sampleinfo, levels)
newlev$cellid <- cellnall
dd <- genefu::setcolclass.df(df=dd, colclass=sapply(sampleinfo, class), factor.levels=newlev)
dd[rownames(sampleinfo), colnames(sampleinfo)] <- sampleinfo
dd[ , "cellid"] <- rownames(dd)
dd <- data.frame(dd, "tissueid"=celline.cgp[match(dd[ , "cellid"], celline.cgp[ , "cellid"]), "tissueid"])
sampleinfo<- dd

## update druginfo
dd <- data.frame(matrix(NA, nrow=length(drugnall), ncol=ncol(druginfo), dimnames=list(drugnall, colnames(druginfo))))
newlev <- sapply(druginfo, levels)
newlev$drugid <- sapply(strsplit(drugnall, split="_"), function(x) { return(x[2]) })
dd <- genefu::setcolclass.df(df=dd, colclass=sapply(druginfo, class), factor.levels=newlev)
dd[match(rownames(druginfo), drugnall), colnames(druginfo)] <- druginfo
dd[ , "drugid"] <- dd[ , "drug.name2"]
druginfo <- dd
myx <- match(c("drugid", "drug.name"), colnames(druginfo))
druginfo <- druginfo[ , c(myx, setdiff(1:ncol(druginfo), myx)), drop=FALSE]

write.csv(celline.cgp, file=file.path(path.out, "cell_line_annotation_COSMIC.csv"), row.names=FALSE)
write.csv(sampleinfo, file=file.path(path.out, "cell_line_annotation_CGP.csv"), row.names=FALSE)
write.csv(druginfo, file=file.path(path.out, "drug_annotation_CGP.csv"), row.names=FALSE)