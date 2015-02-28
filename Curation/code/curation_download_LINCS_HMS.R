## path to files
path.data <- file.path("data", "LINCS_HMS")
path.drug <- file.path(path.data, "drug")
path.cell <- file.path(path.data, "celline")
## create directories
if(!file.exists(path.data)) { dir.create(path.data, showWarnings=FALSE, recursive=TRUE) }
if(!file.exists(path.drug)) { dir.create(path.drug, showWarnings=FALSE, recursive=TRUE) }
if(!file.exists(path.cell)) { dir.create(path.cell, showWarnings=FALSE, recursive=TRUE) }

   ## cell line information
  message("Download cell line information")
  myfn <- file.path(path.cell, "tmp", "LINCS_HMS.csv")
  if (!file.exists(myfn)) {
    if(!file.exists(file.path(path.cell, "tmp"))) { dir.create(file.path(path.cell, "tmp"), showWarnings=FALSE, recursive=TRUE) }
    dwl.status <- download.file(url="http://lincs.hms.harvard.edu/db/cells/?search=&output_type=.csv", destfile=myfn, quiet=TRUE)
    if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
    cell_info<-read.csv(myfn)
  } else {
      cell_info<-read.csv(myfn)
  }
  


  cell_annot_info <- data.frame(cellid=cell_info[,"Cell.Name"], tissueid=cell_info[,"Organ"])
  
  cell_info <- cbind(cell_info, cell_annot_info)
  
  write.csv(cell_info, file=file.path(path.out, "cell_line_annotation_LINCS_HMS.csv"))
