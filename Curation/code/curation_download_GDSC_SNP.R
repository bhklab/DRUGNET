
load("~/Documents/PharmacoGx_test/working/data/GDSC/cnv/gdsc.cnv.eset.RData")
write.csv(cbind("cellid"=rownames(pData(gdsc.eset)), "tissueid"=as.character(pData(gdsc.eset)[,"Tissue"])), file="output/cell_line_annotation_GDSC_SNP.csv", row.names=F)






##add to matching
cell_line_annotation_GDSC_SNP <- read.csv("output/cell_line_annotation_GDSC_SNP.csv", header=T)
cell_all <- read.csv("output/cell_annotation_all.csv", na.strings=c("", " ", "NA"))
myx <- cell_line_annotation_GDSC_SNP$cellid[which(is.na(match(cell_line_annotation_GDSC_SNP$cellid, cell_all$unique.cellid)))]

matching_cell <- read.csv("matching_cell.csv", header=T, na.strings=c("", " ", "NA"))
ss <- sapply(myx, function (x, y) {
  dd <- stringdist::stringdist(a=x, b=y)
  dd <- paste(y[!is.na(dd) & dd == min(dd, na.rm=TRUE)], collapse="|")
  if (length(dd) == 0) { dd <- "no_match"}
  return (dd)
}, y=matching_cell$unique.cellid)

#matching_cell[is.na(matching_cell)]<-""
#write.csv(matching_cell[order(matching_cell$unique.cellid),], file="matching_cell.csv", row.names=F)
