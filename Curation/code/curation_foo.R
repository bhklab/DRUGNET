########################
## Benjamin Haibe-Kains
## Code under License Artistic-2.0
## September 15, 2013
########################

########################
## functions


## Pareto distribution
## http://stats.stackexchange.com/questions/78168/how-to-know-if-my-data-fits-pareto-distribution

getCosmic <- function(em, passw, directory="tmp") {
  if (missing(em)) { stop ("Email must be provided") }
  if (missing(passw)) { stop ("Password must be provided") } 
  myagent <- "Firefox/23.0" 
  mycurl <- RCurl::getCurlHandle()
  RCurl::curlSetOpt(cookiejar=file.path(directory, "cookies.txt"), useragent=myagent, followlocation=TRUE, autoreferer=TRUE, curl=mycurl)
  params <- list("email"=em, "password"=passw)
  resp <- RCurl::postForm("https://cancer.sanger.ac.uk/cosmic/login", .params=params, curl=mycurl, style="POST")
  RCurl::curlSetOpt(cookiejar=file.path(directory, "cookies.txt"), useragent=myagent, curl=mycurl)
  resp <- RCurl::getBinaryURL("http://cancer.sanger.ac.uk/files/cosmic/current_release/CosmicCompleteExport.tsv.gz", curl=mycurl)
  if (!file.exists(file.path(directory))) { dir.create(file.path(directory), showWarnings=FALSE, recursive=TRUE) }
  myfl <- file(file.path(directory, "CosmicCompleteExport.tsv.gz"), "wb")
  writeBin(resp, myfl)
  close(myfl)
  unlink(file.path(directory, "cookies.txt"))
  return (0)
}

## or use gPdtest::gpd.test
pareto.mle <- function (x) {
  xm <- min(x)
  alpha <- length(x)/(sum(log(x))-length(x)*log(xm))
  return(list(xm = xm, alpha = alpha))
}

pareto.test <- function (x, B=1000) {
  require(gPdtest)
  # distribution, cdf, quantile and random functions for Pareto distributions
  dpareto <- function(x, xm, alpha) {
    ifelse(x > xm , alpha*xm**alpha/(x**(alpha+1)), 0)
  }
  ppareto <- function(q, xm, alpha) {
    ifelse(q > xm , 1 - (xm/q)**alpha, 0 )
  }
  qpareto <- function(p, xm, alpha) {
    ifelse(p < 0 | p > 1, NaN, xm*(1-p)**(-1/alpha))
  }
  rpareto <- function(n, xm, alpha) {
    qpareto(runif(n), xm, alpha)
  }
  
  # a <- pareto.mle(x)
  a <- gPdtest::gpd.fit(x=x, method="amle")
  a <- list(xm=a[2,1], alpha=a[1,1])

  # KS statistic
  D <- ks.test(x, function (q) { return (ppareto(q, a$xm, a$alpha)) }, exact=FALSE)$statistic

  # estimating p value with parametric bootstrap
  n <- length(x)
  emp.D <- numeric(B)
  for(b in 1:B) {
    xx <- rpareto(n, a$xm, a$alpha);
    # aa <- pareto.mle(xx)
    aa <- gPdtest::gpd.fit(x=x, method="amle")
    aa <- list(xm=aa[2,1], alpha=aa[1,1])
    emp.D[b] <- ks.test(xx, function(q) { return (ppareto(q, aa$xm, aa$alpha)) }, exact=FALSE)$statistic
  }
  pp <- sum(emp.D > D)/B
  if (pp == 0) { pp <- 1 / (B + 1) }
  return(list("xm" = a$xm, "alpha"=a$alpha, "D"=D, "p"= pp))
}

## compute overlap for Venn diagrams such as library VennDiagram
overlap <- function(l) {
  results <- lapply(l, unique)
  # combinations of m elements of list l
  for (m in seq(along=l)[-1]) {
    # generate and iterate through combinations of length m
    for (indices in combn(seq(length(l)), m, simplify=FALSE)) {
      # make name by concatenating the names of the elements
      # of l that we are intersecting
      name1 <- paste(names(l)[indices[-m]], collapse="_")
      name2 <- names(l)[indices[m]]
      name <- paste(name1, name2, sep="_")
      results[[name]] <- intersect(results[[name1]], results[[name2]])
    }
  }
  return(results)
}

## courtesy of Matthew McCall
celfileDateHour <- function(filename) {
	require(affyio)
	h <- affyio::read.celfile.header(filename, info="full")
	#ddate <- grep("/", strsplit(h$DatHeader, " ")[[1]], value=TRUE)
	#ddate <- strsplit(ddate, split="/")[[1]]
	#CC <- ifelse(substr(ddate[3],1,1)=="9", "19", "20")
	if(length(h$ScanDate) > 0) {
	    h$ScanDate <- gsub(pattern="T", replacement=" ", x=h$ScanDate)
	    ddate <- strsplit(h$ScanDate, " ")[[1]]
    } else { ddate <- rep(NA, 2)}
    names(ddate) <- c("day", "hour")
	return(ddate)
}

celfileChip <- function(filename) {
	require(affyio)
	h <- affyio::read.celfile.header(filename, info="full")
	return(as.character(h$cdfName))
}


## sample size calculation for correlation coefficients (Pearson, kendall and SPearman)
## Bonett, D. G., & Wright, T. A. (2000). Sample size requirements for estimating pearson, kendall and spearman correlations. Psychometrika, 65(1), 23–28. doi:10.1007/BF02294183
## example: power.cor(rho=0.5, w=0.1, alpha=0.05, method="spearman")
power.cor <- function(rho, w, alpha=0.05, method=c("pearson", "kendall", "spearman")) {
  method <- match.arg(method)
  bb <- c(3, 4, 3)
  cc <- c(1, sqrt(0.437), sqrt(1 + (rho^2 / 2)))
  names(bb) <- names(cc) <- c("pearson", "kendall", "spearman")
  bb <- bb[method]
  cc <- cc[method]
  nn0 <- 4 * cc^2 * (1 - rho^2)^2 * (qnorm(p=alpha/2, lower.tail=FALSE) / w) + bb
  if(nn0 < 10) { nn0t <- 10 } else { nn0t <- ceiling(nn0) }
  w0w <- sqrt(nn0t - bb) / sqrt(nn0 - bb)
  nn <- ceiling((nn0 - bb) * w0w^2 + bb)
  return(nn)
}

## intersection of more than 2 sets
fold <- function(f, x, y, ...){
    if (missing(...)) { f(x, y) } else { f(x, fold(f, y, ...)) }
}

#A function that is sometimes useful in determining the 
#coordinate(i.e. row and column number) of a matrix position
#(and vice-versa). 
#Either a vector of positions ("pos") 
#OR a 2 column matrix of matrix coordinates, ("coord", i.e. cbind(row,col)), 
#AND the matrix dimentions must be supplied (dim.mat, i.e. c(nrow,ncol)).
pos2coord <- function(pos, coord, dim.mat=NULL) {
  if(missing(pos) && missing(coord) || missing(dim.mat)){
    stop("must supply either 'pos' or 'coord', and 'dim.mat'")
  }
  if(missing(pos) && !missing(coord)) {
    if(!is.matrix(coord)) { coord <- t(coord) }
    if(any(coord < 0) || any(coord[ , 1] > dim.mat[1]) || any(coord[ , 2] > dim.mat[2])) { stop("wrong coordinates") }
    pos <- ((coord[ , 2] - 1) * dim.mat[1]) + coord[ , 1] 
    return(pos)
  } else {
    coord <- matrix(NA, nrow=length(pos), ncol=2)
    coord[ , 1] <- ((pos - 1) %% dim.mat[1]) + 1
    coord[ , 2] <- ((pos - 1) %/% dim.mat[1]) + 1
    return(coord)
  }
}


blueyellow <- function(N) {
  x <- (1:N)-1
  rgb(x, x, rev(x), maxColorValue=N)
}

myBluePalette <- colorRampPalette(colors=blues9[-(1:3)])

myScatterPlot <- function(x, y, method=c("plain", "transparent", "smooth"), transparency=0.10, smooth.pch=".", pch=16, minp=50, col=blues9[7], smooth.col=c("white", blues9), ...) {
  require(grDevices) || stop("Library grDevices is not available!")
  method <- match.arg(method)
  if (length(col) != length(x)) {
    col <- rep(col, length.out=length(x))
  }
  ccix <- complete.cases(x, y)
  x <- x[ccix]
  y <- y[ccix]
  col <- col[ccix]

  if (sum(ccix) < minp) {
    ## too few points, no transparency, no smoothing
    if (sum(ccix) > 0) { rr <- plot(x=x, y=y, col=col, pch=pch, ...) } else { rr <- plot(x=x, y=y, col=col, pch=pch, ...) }
  } else {
    ## enough data points
    switch(method,
      "plain"={
        rr <- plot(x=x, y=y, col=col, pch=pch, ...)
      },
      "transparent"={
        myrgb <- sapply(col, grDevices::col2rgb, alpha=FALSE) / 255
        myrgb <- apply(myrgb, 2, function (x, transparency) {
          return (rgb(red=x[1], green=x[2], blue=x[3], alpha=transparency, maxColorValue=1))
          }, transparency=transparency)
        rr <- plot(x=x, y=y, pch=pch, col=myrgb, ...)
      },
      "smooth"={
        rr <- smoothScatter(x=x, y=y, col="lightgray", colramp=colorRampPalette(smooth.col), pch=smooth.pch, ...)
      }
    )
  }
  
  invisible(rr)
}

gene.drug.assocs <- function(x, y, z, method=c("lm", "cindex")) {
  method <- match.arg(method)
  switch(method, 
    "lm" = {
      ## standardized coefficient using linear regression controlled for tissue type
      xx <- as.numeric(x)
      names(xx) <- names(x)
      ## avoid underestimation of sd due to truncated or placeholders values
      xxsd <- sd(xx[!duplicated(xx)], na.rm=TRUE)
      x <- xx / xxsd
      if(!is.factor(y)) {
        yy <- as.numeric(y)
        names(yy) <- names(y)
        ## avoid underestimation of sd due to truncated or placeholders values
        yysd <- sd(yy[!duplicated(yy)], na.rm=TRUE)
        y <- yy / yysd
      } else {
        if(sum(table(y) > 0) != 2) { stop("y must have 2 levels!") }
      }
      ff <- sprintf("y ~ x")
      nn <- sum(complete.cases(x, y))
      if(sum(table(z) > 0) > 1) {
        ff <- sprintf("%s + z", ff)
        nn <- sum(complete.cases(x, y, z))
      }
      if(nn >= 3) {
        cc <- summary(glm(formula(ff), family=ifelse(is.factor(y), "binomial", "gaussian")))
        if(is.element("x", rownames(cc$coefficients))) {
          res <- c("estimate"=cc$coefficients["x", 1], "se"=cc$coefficients["x", 2], "n"=nn, "stat"=cc$coefficients["x", 3], "pvalue"=cc$coefficients["x", 4])
          } else { res <- c("estimate"=NA, "se"=NA, "n"=0, "stat"=NA, "pvalue"=NA) }
      } else { res <- c("estimate"=NA, "se"=NA, "n"=nn, "stat"=NA, "pvalue"=NA) }
    },
    "cindex" = {
      ## Somers' D index stratified by tissue type
      requite(survcomp)
      cc <- survcomp::concordance.index(x=-x, cl=y, strat=z, outx=TRUE, method="noether", na.rm=TRUE)
      ss <- (cc$c.index - 0.5) / cc$se
      res <- c("estimate"=(cc$c.index - 0.5) * 2, "se"=cc$se * 2, "n"=cc$n, "stat"=ss, "pvalue"=cc$p.value)
    })
    return(res)
}


gene.drug.call.assocs <- function(x, y, z) {

  ## standardized coefficient using linear regression controlled for tissue type
  xx <- as.numeric(x)
  names(xx) <- names(x)
  ## avoid underestimation of sd due to truncated or placeholders values
  xxsd <- sd(xx[!duplicated(xx)], na.rm=TRUE)
  x <- xx / xxsd
  if(sum(table(y) > 0) != 2) { stop("y must have 2 levels!") }
  ff <- sprintf("y ~ x")
  nn <- sum(complete.cases(x, y))
  if(sum(table(z) > 0) > 1) {
    ff <- sprintf("%s + z", ff)
    nn <- sum(complete.cases(x, y, z))
  }
  if(nn >= 3) {
    ## remove intermediate sensitivity calls
    levels(y)[-c(1, length(levels(y)))] <- NA
    cc <- summary(glm(formula(ff), family="binomial"))
    if(is.element("x", rownames(cc$coefficients))) {
      res <- c("estimate"=cc$coefficients["x", 1], "se"=cc$coefficients["x", 2], "n"=nn, "stat"=cc$coefficients["x", 3], "pvalue"=cc$coefficients["x", 4])
      } else { res <- c("estimate"=NA, "se"=NA, "n"=0, "stat"=NA, "pvalue"=NA) }
  } else { res <- c("estimate"=NA, "se"=NA, "n"=nn, "stat"=NA, "pvalue"=NA) }
    return(res)
}


########################
## Benjamin Haibe-Kains
## All rights Reserved
## October 23, 2013
########################

## Drug sensitivity calling using waterfall plots
## Method:
## 1. Sensitivity calls were made using one of IC50, ActArea or Amax
## 2. Sort log IC50s (or ActArea or Amax) of the cell lines to generate a “waterfall distribution”
## 3. Identify cutoff:
##  3.1 If the waterfall distribution is non-linear (pearson cc to the linear fit <=0.95), estimate the major inflection point of the log IC50 curve as the point on the curve with the maximal distance to a line drawn between the start and end points of the distribution.
##  3.2 If the waterfall distribution appears linear (pearson cc to the linear fit > 0.95), then use the median IC50 instead.
## 4. Cell lines within a 4-fold IC50 (or within a 1.2-fold ActArea or 20% Amax difference) difference centered around this inflection point are classified as being “intermediate”,  cell lines with lower IC50s (or ActArea/Amax values) than this range are defined as sensitive, and those with IC50s (or ActArea/Amax) higher than this range are called “insensitive”.
## 5. Require at least x sensitive and x insensitive cell lines after applying these criteria (x=5 in our case).

## Input:
##  ic50: IC50 values in micro molar (positive values)
##  actarea: Activity Area, that is area under the drug activity curve (positive values)
##  amax: Activity at max concentration (positive values)
##  intermediate.fold: vector of fold changes used to define the intermediate sensitivities for ic50, actarea and amax respectively
callingWaterfallAUC <- function (x, drug.class=c("Unspecified", "Targeted", "Cytotoxic"), drug.ks=0.5, intermediate.fold=c("Targeted"=1.2, "Cytotoxic"=2), cor.min.linear=0.95, name="Drug", plot=FALSE) {
  
  if (is.null(names(x))) { names(x) <- paste("X", 1:length(x), sep=".") }
  xx <- x[complete.cases(x)]
  ## not enough data
  if (length(xx) < 3) {
    tt <- array(NA, dim=length(x), dimnames=list(names(x)))
    if (any(intermediate.fold == 0)) {
      tt <- factor(tt, levels=c("resistant", "sensitive"))
    } else {
      tt <- factor(tt, levels=c("resistant", "intermediate", "sensitive"))
    }
    return (tt)
  }
  mad.stat <- mad(xx)
  ## determine drug class if unspecified
  drug.class <- match.arg(drug.class)
  if (drug.class == "Unspecified") {
    ## using the KS statistics on the pareto distribution of AUC values
    ks.stat <- pareto.test(x=xx, B=1000)$D    
    if (ks.stat >= drug.ks) {
      drug.class <- "Cytotoxic"
    } else {
      drug.class <- "Targeted"
    }
  }
  interfold <- intermediate.fold[drug.class]
  ylabel <- "AUC"
  
  oo <- order(xx, decreasing=TRUE)
  ## test linearity with Pearson correlation
  cc <- cor.test(-xx[oo], 1:length(oo), method="pearson")
  ## line between the two extreme sensitivity values
  dd <- cbind("y"=xx[oo][c(1, length(oo))], "x"=c(1, length(oo)))
  rr <- lm(y ~ x, data=data.frame(dd))
  ## compute distance from sensitivity values and the line between the two extreme sensitivity values
  ddi <- apply(cbind(1:length(oo), xx[oo]), 1, function(x, slope, intercept) {
    return(distancePointLine(x=x[1], y=x[2], slope=slope, intercept=intercept))
  }, slope=rr$coefficients[2], intercept=rr$coefficients[1])
  if(cc$estimate > cor.min.linear){
    ## approximately linear waterfall
    cutoff <- which.min(abs(xx[oo] - median(xx[oo])))
    cutoffn <- names(cutoff)[1]
  } else {
    ## non linear waterfall
    ## identify cutoff as the maximum distance
    cutoff <- which.max(abs(ddi))
    cutoffn <- names(ddi)[cutoff]
  }
  ## identify intermediate sensitivities
  if (interfold == 0) {
    rang <- c(xx[oo][cutoff], xx[oo][cutoff])
  } else {
    rang <- c(xx[oo][cutoff] / interfold, xx[oo][cutoff] * interfold)
  }
  
  ## outlier detection for targeted drugs
  if (drug.class == "Targeted") {
    require(extremevalues)
    xx2 <- xx
    xx2[xx == 0] <- 1e-16
    auc.out <- extremevalues::getOutliersI(y=xx2, rho=c(1,ceiling(length(xx) * 0.10)), FLim=c(0.01,0.99), distribution="pareto")$iRight
    if (length(auc.out) > 0) {
      rang[2] <- min(xx[names(auc.out)])
    } else {
      rang[2] <- max(xx)
    }
  }
  
  ## compute calls
  calls <- array(NA, dim=length(xx), dimnames=list(names(xx)))
  calls[xx < rang[1]] <- "resistant"
  calls[xx >= rang[2]] <- "sensitive"
  calls[xx >= rang[1] & xx < rang[2]] <- "intermediate"
  
  if (plot) {
    par(mfrow=c(2, 1))
    ccols <- rainbow(4)
    mycol <- rep("grey", length(xx))
    names(mycol) <- names(xx)
    mycol[calls == "sensitive"] <- ccols[2]
    mycol[calls == "intermediate"] <- ccols[3]
    mycol[calls == "resistant"] <- ccols[4]
    mycol[cutoffn] <- ccols[1]
    mypch <- rep(16, length(xx))
    names(mypch) <- names(xx)
    mypch[cutoffn] <- 19
    plot(xx[oo], col=mycol[oo], pch=mypch[oo], ylab=ylabel, main=sprintf("%s\nWaterfall", name))
    points(x=cutoff, y=xx[cutoffn], pch=mypch[cutoffn], col=mycol[cutoffn])
    abline(a=rr$coefficients[1], b=rr$coefficients[2], lwd=2, col="darkgrey")
    lines(x=c(cutoff, cutoff), y=c(par("usr")[3], xx[cutoffn]), col="red")
    lines(x=c(par("usr")[1], cutoff), y=c(xx[cutoffn], xx[cutoffn]), col="red")
    legend("topright", legend=c(sprintf("resistant (n=%i)", sum(!is.na(calls) & calls == "resistant")), sprintf("intermediate (n=%i)", sum(!is.na(calls) & calls == "intermediate")), sprintf("sensitive (n=%i)", sum(!is.na(calls) & calls == "sensitive")), "cutoff", sprintf("R=%.3g", cc$estimate), sprintf("KS=%.3g", ks.stat), sprintf("MAD=%.3g", mad.stat)), col=c(rev(ccols), NA, NA, NA), pch=c(16, 16, 16, 19, NA, NA, NA), bty="n")
        
    plot(ddi, pch=mypch[oo], col=mycol[oo], ylab="Distance", main=sprintf("%s\n%s", name, "Distance from min--max line"))
    points(x=cutoff, y=ddi[cutoffn], pch=mypch[cutoffn], col=mycol[cutoffn])
    legend("topright", legend=c("resistant", "intermediate", "sensitive", "cutoff"), col=rev(ccols), pch=c(16, 16, 16, 19), bty="n")
  } 
  
  tt <- rep(NA, length(x))
  names(tt) <- names(x)
  tt[names(calls)] <- calls
  if (interfold == 0) {
    tt <- factor(tt, levels=c("resistant", "sensitive"))
  } else {
    tt <- factor(tt, levels=c("resistant", "intermediate", "sensitive"))
  }
  return(tt)  
}


## improvement of the heatmap.2 function  to include more bars
## "https://gist.github.com/nachocab/3853004"
#
# EXAMPLE USAGE
# # example of colsidecolors rowsidecolors (single column, single row)
# mat <- matrix(1:100, byrow=T, nrow=10)
# column_annotation <- sample(c("red", "blue", "green"), 10, replace=T)
# column_annotation <- as.matrix(column_annotation)
# colnames(column_annotation) <- c("Variable X")
# row_annotation <- sample(c("red", "blue", "green"), 10, replace=T)
# row_annotation <- as.matrix(t(row_annotation))
# rownames(row_annotation) <- c("Variable Y")
# heatmap.3(mat, RowSideColors=row_annotation, ColSideColors=column_annotation)
#  
# # multiple column and row
# mat <- matrix(1:100, byrow=T, nrow=10)
# column_annotation <- matrix(sample(c("red", "blue", "green"), 20, replace=T), ncol=2)
# colnames(column_annotation) <- c("Variable X1", "Variable X2")
# row_annotation <- matrix(sample(c("red", "blue", "green"), 20, replace=T), nrow=2)
# rownames(row_annotation) <- c("Variable Y1", "Variable Y2")
# heatmap.3(mat, RowSideColors=row_annotation, ColSideColors=column_annotation)
heatmap.3 <- function(x,
                      Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                      distfun = dist,
                      hclustfun = hclust,
                      dendrogram = c("both","row", "column", "none"),
                      symm = FALSE,
                      scale = c("none","row", "column"),
                      na.rm = TRUE,
                      revC = identical(Colv,"Rowv"),
                      add.expr,
                      breaks,
                      symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                      col = "heat.colors",
                      colsep,
                      rowsep,
                      sepcolor = "white",
                      sepwidth = c(0.05, 0.05),
                      cellnote,
                      notecex = 1,
                      notecol = "cyan",
                      na.color = par("bg"),
                      trace = c("none", "column","row", "both"),
                      tracecol = "cyan",
                      hline = median(breaks),
                      vline = median(breaks),
                      linecol = tracecol,
                      margins = c(5,5),
                      ColSideColors,
                      RowSideColors,
                      side.height.fraction=0.3,
                      cexRow = 0.2 + 1/log10(nr),
                      cexCol = 0.2 + 1/log10(nc),
                      labRow = NULL,
                      labCol = NULL,
                      colabRow = "black",
                      colabCol = "black",
                      key = TRUE,
                      keysize = 1.5,
                      density.info = c("none", "histogram", "density"),
                      denscol = tracecol,
                      symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                      densadj = 0.25,
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      lmat = NULL,
                      lhei = NULL,
                      lwid = NULL,
                      NumColSideColors = 1,
                      NumRowSideColors = 1,
                      KeyValueName="Value",
                      keyTitle="Color Key",...){
 
    invalid <- function (x) {
      if (missing(x) || is.null(x) || length(x) == 0)
          return(TRUE)
      if (is.list(x))
          return(all(sapply(x, invalid)))
      else if (is.vector(x))
          return(all(is.na(x)))
      else return(FALSE)
    }
 
    x <- as.matrix(x)
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }
    retval <- list()
    scale <- if (symm && missing(scale))
        "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col))
        col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none"))
        warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv))
        Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv))
        Colv <- FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv))
        Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x))
        stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1)
        stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2)
        stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote))
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
            c("both", "row"))) {
            if (is.logical(Colv) && (Colv))
                dendrogram <- "column"
            else dedrogram <- "none"
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting row dendogram.")
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
            c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv))
                dendrogram <- "row"
            else dendrogram <- "none"
            warning("Discrepancy: Colv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting column dendogram.")
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else {
        rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc)
            stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
        hcc <- hclustfun(distfun(if (symm)
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        hcc <- hclustfun(distfun(if (symm)
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else {
        colInd <- 1:nc
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow))
        labRow <- if (is.null(rownames(x)))
            (1:nr)[rowInd]
        else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol))
        labCol <- if (is.null(colnames(x)))
            (1:nc)[colInd]
        else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
        retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
        x <- sweep(x, 1, rm)
        retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
        x <- sweep(x, 2, rm)
        retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
        if (missing(col) || is.function(col))
            breaks <- 16
        else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
        if (!symbreaks)
            breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                length = breaks)
        else {
            extreme <- max(abs(x), na.rm = TRUE)
            breaks <- seq(-extreme, extreme, length = breaks)
        }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function")
        col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei))
        lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid))
        lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(4:3, 2:1)
 
        if (!missing(ColSideColors)) {
           #if (!is.matrix(ColSideColors))
           #stop("'ColSideColors' must be a matrix")
            if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
                stop("'ColSideColors' must be a matrix of nrow(x) rows")
            lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
            #lhei <- c(lhei[1], 0.2, lhei[2])
             lhei=c(lhei[1], side.height.fraction*NumColSideColors, lhei[2])
        }
 
        if (!missing(RowSideColors)) {
            #if (!is.matrix(RowSideColors))
            #stop("'RowSideColors' must be a matrix")
            if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
                stop("'RowSideColors' must be a matrix of ncol(x) columns")
            lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
            #lwid <- c(lwid[1], 0.2, lwid[2])
            lwid <- c(lwid[1], side.height.fraction*NumRowSideColors, lwid[2])
        }
        lmat[is.na(lmat)] <- 0
    }
 
    if (length(lhei) != nrow(lmat))
        stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat))
        stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
 
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
 
    if (!missing(RowSideColors)) {
        if (!is.matrix(RowSideColors)){
          par(mar = c(margins[1], 0, 0, 0.5))
          image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)      
        } else {
            par(mar = c(margins[1], 0, 0, 0.5))
            rsc = t(RowSideColors[,rowInd, drop=F])
            rsc.colors = matrix()
            rsc.names = names(table(rsc))
            rsc.i = 1
            for (rsc.name in rsc.names) {
                rsc.colors[rsc.i] = rsc.name
                rsc[rsc == rsc.name] = rsc.i
                rsc.i = rsc.i + 1
            }
            rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
            image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
            if (length(colnames(RowSideColors)) > 0) {
                axis(1, 0:(dim(rsc)[2] - 1)/(dim(rsc)[2] - 1), colnames(RowSideColors), las = 2, tick = FALSE)
            }
        }
        ## separators should be displayed
        # if (!missing(rowsep)) {
        #   for (rsep in rowsep) {
        #     rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
        #   }
        # }
    }
 
    if (!missing(ColSideColors)) {
 
        if (!is.matrix(ColSideColors)){
            par(mar = c(0.5, 0, 0, margins[2]))
            image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
        } else {
            par(mar = c(0.5, 0, 0, margins[2]))
            csc = ColSideColors[colInd, , drop=F]
            csc.colors = matrix()
            csc.names = names(table(csc))
            csc.i = 1
            for (csc.name in csc.names) {
                csc.colors[csc.i] = csc.name
                csc[csc == csc.name] = csc.i
                csc.i = csc.i + 1
            }
            csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
            image(csc, col = as.vector(csc.colors), axes = FALSE)
            if (length(colnames(ColSideColors)) > 0) {
                axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
            }
        }
    }
 
    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
        iy <- nr:1
        if (exists("ddr"))
            ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr"))
        retval$rowDendrogram <- ddr
    if (exists("ddc"))
        retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
        mmat <- ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
            col = na.color, add = TRUE)
    }
    # axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, cex.axis = cexCol)
    text(x=1 + (1:nc) + par("usr")[1], y=par("usr")[3] - par("usr")[4] * 0.02, pos=2, labels=labCol, srt=90, xpd=TRUE, font=2, col=colabCol, cex=cexCol)
    
    if (!is.null(xlab))
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    # axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, cex.axis = cexRow)
    text(x=par("usr")[2], y=iy, pos=4, labels=labRow, xpd=TRUE, font=2, col=colabRow, cex=cexRow)
    
    if (!is.null(ylab))
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr))
        eval(substitute(add.expr))
    if (!missing(colsep))
        for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!missing(rowsep))
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
        retval$vline <- vline
        vline.vals <- scale01(vline, min.scale, max.scale)
        for (i in colInd) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol,
                  lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace %in% c("both", "row")) {
        retval$hline <- hline
        hline.vals <- scale01(hline, min.scale, max.scale)
        for (i in rowInd) {
            if (!is.null(hline)) {
                abline(h = i + hline, col = linecol, lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (!missing(cellnote))
        text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
            col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else plot.new()
    if (!is.null(main))
        title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
        par(mar = c(5, 4, 2, 1), cex = 0.75)
        tmpbreaks <- breaks
        if (symkey) {
            max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
            min.raw <- -max.raw
            tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
            tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
        }
        else {
            min.raw <- min(x, na.rm = TRUE)
            max.raw <- max(x, na.rm = TRUE)
        }
 
        z <- seq(min.raw, max.raw, length = length(col))
        image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
            xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        lv <- pretty(breaks)
        xv <- scale01(as.numeric(lv), min.raw, max.raw)
        axis(1, at = xv, labels = lv)
        if (scale == "row")
            mtext(side = 1, "Row Z-Score", line = 2)
        else if (scale == "column")
            mtext(side = 1, "Column Z-Score", line = 2)
        else mtext(side = 1, KeyValueName, line = 2)
        if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[-omit]
            dens$y <- dens$y[-omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
                lwd = 1)
            axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
            title(keyTitle)
            par(cex = 0.5)
            mtext(side = 2, "Density", line = 2)
        }
        else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
                col = denscol)
            axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
            title(keyTitle)
            par(cex = 0.5)
            mtext(side = 2, "Count", line = 2)
        }
        else title(keyTitle)
    }
    else plot.new()
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
        high = retval$breaks[-1], color = retval$col)
    invisible(retval)
}





