

ben_neel <- read.csv("data/BenNeel/bc_cellines_neel_subtypes.csv", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)
ben_neel$tissueid <- "breast"
write.csv(ben_neel, file="output/cell_line_annotation_BenNeel.csv", row.names=F)


##add to matching
cell_line_annotation_BenNeel <- read.csv("output/cell_line_annotation_BenNeel.csv", header=T)
cell_all <- read.csv("output/cell_annotation_all.csv", na.strings=c("", " ", "NA"))
myx <- cell_line_annotation_BenNeel$cellid[which(is.na(match(cell_line_annotation_BenNeel$cellid, cell_all$unique.cellid)))]

matching_cell <- read.csv("matching_cell.csv", header=T, na.strings=c("", " ", "NA"))
ss <- sapply(myx, function (x, y) {
  dd <- stringdist::stringdist(a=x, b=y)
  dd <- paste(y[!is.na(dd) & dd == min(dd, na.rm=TRUE)], collapse="|")
  if (length(dd) == 0) { dd <- "no_match"}
  return (dd)
}, y=matching_cell$unique.cellid)

###cleaning matching 
# matching_cell[is.na(matching_cell)]<-""
# write.csv(matching_cell[order(matching_cell$unique.cellid),], file="matching_cell.csv", row.names=F)
# cell_all[is.na(cell_all)]<-""
# write.csv(cell_all[order(cell_all$unique.cellid),], file="output/cell_annotation_all.csv", row.names=F)

# cell_all <- read.csv("output/cell_annotation_all.csv", na.strings=c("", " ", "NA"))
# matching <- read.csv("matching_cell_copy.csv", na.strings=c("", " ", "NA"))
# 
# tt <- unlist(apply(matching, MARGIN=1, function(x){table(is.na(x))["TRUE"]}))
# matching <- matching[-which(tt == (ncol(matching)-3)),]
# matching[is.na(matching)]<-""
# write.csv(matching[order(matching$unique.cellid),], file="matching_cell.csv", row.names=F)
# 
# matching <- read.csv("matching_cell.csv", na.strings=c("", " ", "NA"))
# 
# matching <- matching[-which(is.na(match(matching$unique.cellid,cell_all$unique.cellid)))[1],]
# matching[which(is.na(match(matching$unique.cellid,cell_all$unique.cellid)))[1], "unique.cellid"] <- "MCF10F"
# cc <- cell_all[match(matching$unique.cellid,cell_all$unique.cellid),match(colnames(matching)[c(1:6, 8:9,13:14)], colnames(cell_all))]
# matching[, colnames(matching)[c(1:6, 8:9,13:14)]] <- cc 
# matching[is.na(matching)]<-""
# write.csv(matching[order(matching$unique.cellid),], file="matching_cell.csv", row.names=F)
# 
# matching <- read.csv("matching_cell.csv", na.strings=c("", " ", "NA"))
# 
# tt <- unlist(apply(matching, MARGIN=1, function(x){y<-x["unique.cellid"];table(x==y)["FALSE"]}))
# matching <- matching[-which(tt == 1),]
# matching[is.na(matching)]<-""
# write.csv(matching[order(matching$unique.cellid),], file="matching_cell.csv", row.names=F)
# 
# matching <- read.csv("matching_cell.csv", na.strings=c("", " ", "NA"))
