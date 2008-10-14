image.splitSet <- function(x, filter.fdr=1, main="", max.label.length=50, full.names=TRUE, 
                           xlab=NULL, sample.labels=FALSE, 
                           col=c("yellow","red"), invert=FALSE,
                           outfile=NULL, res=72, pointsize=7, ...) {

  # filter splitSet
  if (!is.null(x$qvalues)) {
    select <- x$qvalues < filter.fdr
    x$cuts <- x$cuts[select,]
    x$scores <- x$scores[select]
    x$pvalues <- x$pvalues[select]
    x$qvalues <- x$qvalues[select]
  }

  # determine left labels
  ids <- rownames(x$cuts)
  if (is.null(x$pvalues)) left.labels <- rep("", nrow(x$cuts))
  else left.labels <- format.pval(x$pvalues, digits=2, eps=1e-4)

  # determine right labels if wanted
  right.labels <- NULL
  if (full.names) {
    for (i in 1:length(ids)) {
      rowid <- ids[i]
      if (length(grep("^GO:", rowid)) > 0) { # GO
        tmp <- GOTERM[[rowid]]
        desc <- attr(tmp, "Term")
        onto <- attr(tmp, "Ontology")
      } else { # KEGG
        id <- sub("^KEGG:", "", rowid)
        desc <- KEGGPATHID2NAME[[id]]
        onto <- "KEGG"
      }
      desc <- paste(onto, ":", desc)
      if (nchar(desc) > max.label.length) 
        desc <- paste(substr(desc, 1, max.label.length-3), "...")
      right.labels <- c(right.labels, desc)
    }
  }

  # determine cluster orderings
  dx <- agnes(x$cuts)
  dy <- agnes(t(x$cuts))
  inv.dy <- rep(0, length(dy$order))
  for (i in 1:length(dy$order)) inv.dy[dy$order[i]] <- i
  inv.dx <- rep(0, length(dx$order))
  for (i in 1:length(dx$order)) inv.dx[dx$order[i]] <- i

  # determine widths of labels to optimize layout
  if (is.null(outfile)) X11(pointsize=pointsize)
  else postscript(pointsize=pointsize)
  plot(c(0), type="n") # plot nothing to make strwidth happy
  left.label.widths  <- strwidth( left.labels, units="inches")
  right.label.widths <- strwidth(right.labels, units="inches")
  if (sample.labels) sample.label.widths <- strwidth(colnames(x$cuts), units="inches")
  else sample.label.widths <- 0
  left.add <- 0.5
  margins <- c(0.5 + max(sample.label.widths), max(left.label.widths) + left.add, 
               0.5, max(right.label.widths) + 0.5)
  height <- 1.0 + max(sample.label.widths) + 1.3 * nrow(x$cuts) * pointsize/72
  if (sample.labels) 
    width <- max(left.label.widths) + left.add + max(right.label.widths) + 0.5 + 
             ncol(x$cuts) * pointsize/72
  else width <- 11
  dev.off()
  leftsplit  <- 0.2
  rightsplit <- 0.9
  l <- nrow(x$cuts)
  label.pos <- (seq(l)-1)/(l-1)

  # draw graph
  if (is.null(outfile)) X11(width=width, height=height, pointsize=pointsize)
  else { 
    postscript(outfile, width=width, height=height, pointsize=pointsize, 
               horizontal=TRUE, paper="special")
  }
  if (invert) par(bg="black", fg="white", col.axis="white", col.main="white") 

  # draw the image
  par(las=1, mai=margins)
  image(t(x$cuts[dx$order, dy$order]), axes=FALSE, zlim=c(0,1), main=main, col=col, ...)
  grid(ncol(x$cuts), nrow(x$cuts), lty=1, lwd=0.5, col="white")
  if (!is.null(x$pvalues)) {
    axis(side=2, at=label.pos, labels=left.labels[dx$order])
    mtext("p-values", side=3, line=1, at=0, adj=1)
    if (!is.null(x$qvalues)) {
      mtext(paste("FDR: ", round(max(x$qvalues)*100, 1), "%", sep=""), 
            side=1, line=1, at=0, adj=1)
    }
  }
  if (!is.null(right.labels)) axis(side=4, at=label.pos, labels=right.labels[dx$order])
  box()
  l <- ncol(x$cuts)
  label.pos <- (seq(l)-1)/(l-1)
  if (sample.labels) axis(side=1, at=label.pos, labels=colnames(x$cuts)[dy$order], las=2)
  else {
    if (is.null(xlab)) {
      xlab <- paste(nrow(x$cuts), "annotation driven clusterings on",
                    ncol(x$cuts), "patient samples")
    }
    mtext(xlab, side=1, line=2.5)
  }

  if (!is.null(outfile)) dev.off()
  return(invisible(NULL))
}

computeOverlap <- function(cuts, chip.name, expr.mat=NULL) {
  # determine genes associated to cuts
  require(package=paste(chip.name,".db",sep=""), character.only=TRUE)
  genes <- list()
  env <- eval(as.symbol(paste(chip.name, "GO2ALLPROBES", sep="")))
  for (i in grep("GO:", rownames(cuts))) {
    id <- rownames(cuts)[i]
    genes[[i]] <- as.character(get(id, env))
  }
  env <- eval(as.symbol(paste(chip.name, "PATH2PROBE", sep="")))
  for (i in grep("KEGG:", rownames(cuts))) {
    id <- sub("^KEGG:", "", rownames(cuts)[i])
    genes[[i]] <- as.character(get(id, env))
  }
  if (!is.null(expr.mat)) {
    for (i in 1:length(genes)) genes[[i]] <- intersect(genes[[i]], rownames(expr.mat))
  }

  res <- matrix(0, nrow=nrow(cuts), ncol=nrow(cuts))
  for (i in 1:(nrow(cuts))) {
    for (j in i:nrow(cuts)) {
      ninter <- length(intersect(genes[[i]], genes[[j]]))
#      nunion <- length(union(genes[[i]], genes[[j]]))
      res[i,j] <- ninter
      res[j,i] <- res[i,j]
    }
  }
  colnames(res) <- rownames(cuts)
  rownames(res) <- rownames(cuts)
  return(res)
}

computeHamming <- function(cuts) {
  res <- matrix(0, nrow=nrow(cuts), ncol=nrow(cuts))
  for (i in 1:(nrow(cuts)-1)) {
    for (j in (i+1):nrow(cuts)) {
      dist <- sum(cuts[i,] != cuts[j,])
      res[i,j] <- min(dist, ncol(cuts) - dist)
      res[j,i] <- res[i,j]
    }
  }
  colnames(res) <- rownames(cuts)
  rownames(res) <- rownames(cuts)
  return(res)
}

scatterpairs <- function(res, chip.name="hgu133a", expr.mat=NULL) {
  if (nrow(res$cuts) <=1) stop("No pairs available")
  overlap <- as.vector(computeOverlap(res$cuts, chip.name, expr.mat))
  hamming <- as.vector(computeHamming(res$cuts))
  plot(overlap, hamming)
}
