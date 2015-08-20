path.out <- "output"
drug_all <- read.csv(file = file.path(path.out, "drug_annotation_all.csv"), na.strings = c(""," "), stringsAsFactor = F)
subset(drug_all, !is.na(drug_all$GRAY.drugid) & !is.na(drug_all$CTRP.drugid))
subset(drug_all, !is.na(drug_all$GNE.drugid) & !is.na(drug_all$CTRP.drugid))


cell_all <- read.csv(file = file.path(path.out, "cell_annotation_all.csv"), na.strings = c(""," "), stringsAsFactor = F)
rownames(cell_all) <- cell_all$unique.cellid
cell_all_unique.cellid <- cell_all$unique.cellid
names(cell_all_unique.cellid) <- toupper(gsub(badchars, "",cell_all$unique.cellid))


match_cell <- read.csv(file = "matching_cell.csv", stringsAsFactor = F, na.strings = c(""," "))
match_unique.cellid <- match_cell$unique.cellid
names(match_unique.cellid) <- toupper(gsub(badchars, "",match_cell$unique.cellid))


###CTRP test
ctrp_cell <- read.csv(file = file.path(path.out, "cell_line_annotation_CTRP.csv"), na.strings = c(""," "), stringsAsFactor = F)
ctrp <- ctrp_cell$cellid
names(ctrp) <- toupper(gsub(badchars, "",ctrp))
length(ctrp)

match_ctrp <- intersect(names(cell_all_unique.cellid),names(ctrp))
length(match_ctrp)

match_ctrp.unique.cellid <- cell_all_unique.cellid[names(cell_all_unique.cellid) %in% match_ctrp]
mm <- which(match_cell$unique.cellid %in% match_ctrp.unique.cellid)
match_cell[mm,"CTRP.cellid"] <- ctrp[toupper(gsub(badchars, "",match_cell[mm,"unique.cellid"]))]
match_cell[is.na(match_cell)]<-""
write.csv(match_cell, file = "matching_cell.csv",row.names=FALSE)


reminder <- match_ctrp.unique.cellid[which(!(match_ctrp.unique.cellid %in% match_cell$unique.cellid))]
reminder <- data.frame(cbind("unique"= reminder[names(reminder)], "ctrp"= ctrp[names(reminder)]))
reminder <- subset(reminder, reminder$unique != reminder$ctrp)
write.csv(reminder,file = "reminder.csv")

match_cell <- read.csv(file = "matching_cell.csv", stringsAsFactor = F, na.strings = c(""," "))

tt<- cell_all[reminder$unique,c("unique.cellid","CGP.cellid","CCLE.cellid","GSK.cellid","NCI60.cellid","GRAY.cellid","GNE.cellid")]
tt <-cbind(tt,"CTRP.cellid" = reminder$ctrp,
              "PGP.cellid" = NA,
              "SU2C.cellid" = NA, 
              "LINCS.cellid" = NA, 
              "LINCS_HMS.cellid" = NA, 
              "COSMIC.tissueid" = cell_all[reminder$unique,"COSMIC.tissueid"],
              "Comment" = NA)

match_cell <- rbind(match_cell, tt)
match_cell[is.na(match_cell)]<-""
write.csv(match_cell[order(match_cell$unique.cellid),], file = "matching_cell.csv",row.names=FALSE)

####biosample test

all.celllines <- read.csv("output/cell_annotation_all.csv", stringsAsFactors = F, na.strings = c("", " "))
rownames(all.celllines) <-all.celllines$unique.cellid
all.celllines.upper<- as.data.frame(apply(all.celllines, MARGIN = 2, function(x){toupper(gsub(pattern = badchars, replacement = "", x))}))

sampleiftt <- read.csv("data/biosample/cellline/biosample.csv", stringsAsFactors = F, na.strings = c("", " "))
sampleiftt <- subset(sampleiftt, !duplicated(sampleiftt$cell.line))
rownames(sampleiftt) <- sampleiftt$cell.line

cellids <- grep("cellid",colnames(all.celllines))
tt <- NULL
ss <- NULL
for(cellid in cellids)
{
  tt <- union(tt, rownames(subset(all.celllines, all.celllines[,cellid] %in%sampleiftt$"cell.line")))
  ss <- union(ss, rownames(subset(sampleiftt, sampleiftt$"cell.line" %in% all.celllines[,cellid])))
}
length(tt)
length(ss)
for(cellid in cellids)
{
  tt <- union(tt, rownames(subset(all.celllines.upper, all.celllines.upper[,cellid] %in% toupper(gsub(pattern = badchars, replacement = "", sampleiftt$cell.line)))))
  ss <- union(ss, rownames(subset(sampleiftt, toupper(gsub(pattern = badchars, replacement = "", sampleiftt$cell.line)) %in% all.celllines[,cellid])))
}
length(tt)
length(ss)

reminder <- matrix(NA, ncol = 2, nrow= max(length(all.celllines[which(!(all.celllines$unique.cellid %in% tt)),"unique.cellid"]), 
                                           length(sampleiftt[which(!(sampleiftt$cell.line %in% ss)),"cell.line"])))
colnames(reminder) = c("unique.cellid", "biosample.cellid")
reminder[1:length(all.celllines[which(!(all.celllines$unique.cellid %in% tt)),"unique.cellid"]), "unique.cellid"] <- all.celllines[which(!(all.celllines$unique.cellid %in% tt)),"unique.cellid"]
reminder[,"biosample.cellid"] <- sampleiftt[which(!(sampleiftt$cell.line %in% ss)),"cell.line"]
write.csv(reminder, file = "biosample.reminder.csv")

reminder <- read.csv(file = "biosample.reminder.csv", stringsAsFactors = F)
nn <- vector(length = length(gsub(".cellid","", colnames(all.celllines)[grep("cellid",colnames(all.celllines))])))
names(nn) <- gsub(".cellid","", colnames(all.celllines)[grep("cellid",colnames(all.celllines))])
nn[names(nn)] <- 0
for(i in 1:length(reminder$unique.cellid))
{
  tt<- all.celllines[reminder$unique.cellid[i],grep("cellid",colnames(all.celllines))]
  aa <- gsub(".cellid","", names(tt[which(!is.na(tt))]))
  nn[aa] <- nn[aa] + 1
}
bb <- all.celllines[which(all.celllines$unique.cellid %in% reminder$unique.cellid),]
pdf(file = "unmatched_cell.pdf", width = 5, height =  5)
pie(nn[2:length(nn)], col=mycol, labels = sprintf("%s(%s)",names(nn)[2:length(nn)], nn[2:length(nn)]), main="distribution of 330 unmatched cell lines across studies")
dev.off()

library("gdata")
ccle.biosample <- gdata::read.xls(xls = "~/Downloads/nature14397-s1/Supplementary Tables 1-9 and 11-14.xlsx", sheet = 9)
bb.ccle <- bb[which(!is.na(bb$CCLE.cellid)), "CCLE.cellid"]
ccle.unmatched <- ccle.biosample[which(ccle.biosample$Cell.line.primary.name %in% bb.ccle),]
setdiff(bb.ccle, ccle.biosample$Cell.line.primary.name)
nn<- table(ccle.unmatched$Source)
  pdf(file = "unmatched_ccle_source.pdf", width = 8, height =  8)
  pie(nn, col=mycol, labels = sprintf("%s(%s)",names(nn), nn), main="source of 180 unmatched ccle cell lines")
  dev.off()

supp.sheet2 <- gdata::read.xls(xls = "~/Downloads/nature14397-s1/Supplementary Tables 1-9 and 11-14.xlsx", sheet = 2)
##read.xls bug 
##This note form Note column has considered as a cell line name: both of which are compatible with HPB-ALL only.\
supp.sheet2<- supp.sheet2[which(regexpr("both", supp.sheet2$Cell.Line.Name) < 0),]

nrow(supp.sheet2[which(duplicated(supp.sheet2$Cell.Line.Name)),]) ##6  ##These cell lines are duplicated but have different information in the other columns
nrow(supp.sheet2[which(duplicated(supp.sheet2$Canonical.Name)),]) ##260

cellids <- grep("cellid",colnames(all.celllines))
tt <- NULL
ss <- NULL
for(cellid in cellids)
{
  tt <- union(tt, rownames(subset(all.celllines, all.celllines[,cellid] %in%supp.sheet2$"Cell.Line.Name")))
  ss <- union(ss, rownames(subset(supp.sheet2, supp.sheet2$"Cell.Line.Name" %in% all.celllines[,cellid])))
}
length(tt)
length(ss)
for(cellid in cellids)
{
  tt <- union(tt, rownames(subset(all.celllines.upper, all.celllines.upper[,cellid] %in% toupper(gsub(pattern = badchars, replacement = "", supp.sheet2$"Cell.Line.Name")))))
  ss <- union(ss, rownames(subset(supp.sheet2, toupper(gsub(pattern = badchars, replacement = "", supp.sheet2$"Cell.Line.Name")) %in% all.celllines[,cellid])))
}
length(tt)
length(ss)


reminder <- matrix(NA, ncol = 2, nrow= max(length(all.celllines[which(!(all.celllines$unique.cellid %in% tt)),"unique.cellid"]), 
                                           length(supp.sheet2[which(!(supp.sheet2$Cell.Line.Name %in% ss)),"Cell.Line.Name"])))
colnames(reminder) = c("unique.cellid", "supp.sheet2.cellid")
reminder[1:length(all.celllines[which(!(all.celllines$unique.cellid %in% tt)),"unique.cellid"]), "unique.cellid"] <- all.celllines[which(!(all.celllines$unique.cellid %in% tt)),"unique.cellid"]
reminder[,"supp.sheet2.cellid"] <- supp.sheet2[which(!(supp.sheet2$Cell.Line.Name %in% ss)),"Cell.Line.Name"]
write.csv(reminder, file = "supp.sheet2.reminder.csv")

reminder <- read.csv(file = "supp.sheet2.reminder.csv", stringsAsFactors = F)


nn <- vector(length = length(gsub(".cellid","", colnames(all.celllines)[grep("cellid",colnames(all.celllines))])))
names(nn) <- gsub(".cellid","", colnames(all.celllines)[grep("cellid",colnames(all.celllines))])
nn[names(nn)] <- 0
for(i in 1:length(reminder$unique.cellid))
{
  tt<- all.celllines[reminder$unique.cellid[i],grep("cellid",colnames(all.celllines))]
  aa <- gsub(".cellid","", names(tt[which(!is.na(tt))]))
  nn[aa] <- nn[aa] + 1
}
bb <- all.celllines[which(all.celllines$unique.cellid %in% reminder$unique.cellid),]
pdf(file = "unmatched_cell_supp2.pdf", width = 7, height =  7)
pie(nn[2:length(nn)], col=mycol, labels = sprintf("%s(%s)",names(nn)[2:length(nn)], nn[2:length(nn)]), main="distribution of 223 unmatched cell lines across studies")
dev.off()
bb[which(!is.na(bb$GNE.cellid)), "GNE.cellid"]
bb[which(!is.na(bb$NCI60.cellid)), "NCI60.cellid"]

library("gdata")
ccle.biosample <- gdata::read.xls(xls = "~/Downloads/nature14397-s1/Supplementary Tables 1-9 and 11-14.xlsx", sheet = 9)
bb.ccle <- bb[which(!is.na(bb$CCLE.cellid)), "CCLE.cellid"]
ccle.unmatched <- ccle.biosample[which(ccle.biosample$Cell.line.primary.name %in% bb.ccle),]
setdiff(bb.ccle, ccle.biosample$Cell.line.primary.name)
nn<- table(ccle.unmatched$Source)
pdf(file = "unmatched_ccle_source_supp2.pdf", width = 8, height =  8)
pie(nn, col=mycol, labels = sprintf("%s(%s)",names(nn), nn), main="source of 109(out of 111) unmatched ccle cell lines")
dev.off()


ccle.unmatched[which(ccle.unmatched$Source == "ATCC"), "Cell.line.primary.name"]
ccle.unmatched[which(ccle.unmatched$Source == "ICLC"), "Cell.line.primary.name"]


cgp.biosample <- gdata::read.xls(xls = "~/Downloads/nature14397-s1/Supplementary Tables 1-9 and 11-14.xlsx", sheet = 10)
bb.cgp <- bb[which(!is.na(bb$CGP.cellid)), "CGP.cellid"]
cgp.unmatched <- cgp.biosample[which(cgp.biosample$Sample.name %in% bb.cgp),]
setdiff(bb.cgp, cgp.biosample$Sample.name)

temp <- cgp.unmatched
temp$Source <- ""
temp[which(regexpr("ATCC",temp$Availability..Institute.Address.Catalogue.Number.) > 0),"Source"] = "ATCC"
temp[which(regexpr("Interlab Cell Line Collection",temp$Availability..Institute.Address.Catalogue.Number.)> 0),"Source"] = "ICLC"
temp[which(regexpr("RIKEN",temp$Availability..Institute.Address.Catalogue.Number.)> 0),"Source"] = "JCRB"
temp[which(regexpr("German",temp$Availability..Institute.Address.Catalogue.Number.)> 0),"Source"] = "German collection"
temp[which(regexpr("Frederick",temp$Availability..Institute.Address.Catalogue.Number.)> 0),"Source"] = "Frederick"
temp[which(regexpr("Rinku",temp$Availability..Institute.Address.Catalogue.Number.)> 0),"Source"] = "Rinku Japan"
cgp.unmatched <- temp

nn<- table(cgp.unmatched$Source)
pdf(file = "unmatched_cgp_source_supp2.pdf", width = 8, height =  8)
pie(nn, col=mycol, labels = sprintf("%s(%s)",names(nn), nn), main="source of 75(out of 91) unmatched cgp cell lines")
dev.off()

cgp.unmatched[which(cgp.unmatched$Source == "ATCC"), "Sample.name"]

##########
View(ccle.unmatched)
bio.siff <- setdiff(sampleiftt$"cell line", cell)
dnet.siff <- setdiff(all.celllines$unique.cellid, cell)
names(dnet.siff) <- toupper(gsub(pattern = badchars, replacement = "", dnet.siff))
names(bio.siff) <- toupper(gsub(pattern = badchars, replacement = "", bio.siff))
cell.2 <- intersect(names(dnet.siff),names(bio.siff))

