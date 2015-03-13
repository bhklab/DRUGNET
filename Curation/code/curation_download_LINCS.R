## path to files
path.data <- file.path("data", "LINCS")
path.drug <- file.path(path.data, "drug")
path.cell <- file.path(path.data, "celline")

## read amazon and secret keys
kk <- read.csv(file.path("code", "LINCS_keys.csv"))
amazon_key <- kk["amazon_key", "key"]
secret_key <- kk["secret_key", "key"]
lincs_key <- kk["lincs_key", "key"]

## create directories
if(!file.exists(path.data)) { dir.create(path.data, showWarnings=FALSE, recursive=TRUE) }
if(!file.exists(path.drug)) { dir.create(path.drug, showWarnings=FALSE, recursive=TRUE) }
if(!file.exists(path.cell)) { dir.create(path.cell, showWarnings=FALSE, recursive=TRUE) }

  ## sample information
  message("Download sample information")
  myfn <- file.path(path.cell, "tmp", "inst.info")
  
  
  
  if (!file.exists(myfn)){
    if (!file.exists(file.path(path.data, "aws", "bin", "aws"))){
      
      system(command=paste("cd ", path.data, " && ", "aws/install -i $PWD/", "aws", sep=""))
      
    }
    if(!file.exists(file.path(path.cell, "tmp"))) { dir.create(file.path(path.cell, "tmp"), showWarnings=FALSE, recursive=TRUE) }
    system(command=paste(path.data, "/aws/bin/aws configure set aws_access_key_id ", amazon_key, sep=""))
    system(command=paste(path.data, "/aws/bin/aws configure set aws_secret_access_key ", secret_key , sep=""))
    
    system(command=paste(path.data, "/aws/bin/aws s3 cp s3://data.lincscloud.org/l1000/metadata/inst.info ", myfn, sep=""))
    sample_info <- read.delim(myfn, stringsAsFactors=F)
    
  } else {
    
    sample_info <- read.delim(myfn, stringsAsFactors=F)
    
  }
  
  ## cell line information
  message("Download cell line information")
  myfn <- file.path(path.cell, "tmp", "cell_info.csv")
  if (!file.exists(myfn)) {
    if(!file.exists(file.path(path.cell, "tmp"))) { dir.create(file.path(path.cell, "tmp"), showWarnings=FALSE, recursive=TRUE) }
  
    cell_info<-lincsAPIQuery(infoType="cellinfo", user_key=lincs_key)  
    write.csv(cell_info, file=myfn)
  } else {
      cell_info<-read.csv(myfn)
  }
  
  message("Download drug information")
  myfn <- file.path(path.drug, "tmp", "pert_info.csv")
  if (!file.exists(myfn)) {
    if(!file.exists(file.path(path.drug, "tmp"))) { dir.create(file.path(path.drug, "tmp"), showWarnings=FALSE, recursive=TRUE) }
    
    pert_info <- NULL
    
    #pert_list <- lincsAPIQuery("pertinfo", q="pert_type:trt_cp", d="pert_id", user_key="afe4848514df2b491be1c49b49475cab")
        
    success <- FALSE
    
    while (!success){
      success <- TRUE
      tryCatch(pert_info <- rbind(pert_info, lincsAPIQuery(infoType="pertinfo", query="pert_type:trt_cp", user_key=lincs_key)), error = function(e) {success <- FALSE})
    }
    
    write.csv(pert_info, file=myfn)
  } else {
      
      pert_info <- read.csv(myfn)
  }
  
  
  cell_tested <- cell_info[ cell_info[,"cell_id"] %in% unique(sample_info$cell_id),]
  
  cell_annot_info <- data.frame(cellid=cell_tested[,"cell_id"], tissueid=cell_tested[,"cell_lineage"])
  
  cell_tested <- cbind(cell_tested, cell_annot_info)
  
  write.csv(cell_tested, file=file.path(path.out, "cell_line_annotation_LINCS.csv"))
  write.csv(pert_info, file=file.path(path.out, "drug_annotation_LINCS.csv"))
