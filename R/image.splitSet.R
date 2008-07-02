image.splitSet <- function(x, filter.fdr=1, main="", show.graph=FALSE, 
                           max.label.length=50, full.names=TRUE, 
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
  left.add <- ifelse(show.graph, 0.2, 0.5)
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
  if (show.graph) {
    par(mai=c(margins[1], 0.5, margins[3], 0.0), fig=c(0.0, leftsplit, 0.0, 1.0))
    root.dist <<- rep(0, length(ids))
    names(root.dist) <<- ids
    graph.lines <- collect.lines("GO:0003673", "", 0, grep("GO:", ids, value=TRUE))
    tmp <- 0.5/(ncol(x$cuts)-1)
    image(matrix(0,2,2), col="transparent", xlab="", ylab="", axes=FALSE,
         xlim=c(min(root.dist)-1, max(root.dist)+1), ylim=c(-tmp,1.0+tmp))
    points(root.dist, label.pos[inv.dx], pch=19)
#    for (i in 1:length(root.dist)) {
#      lines(c(root.dist[i], max(root.dist)), rep(label.pos[inv.dx[i]], 2), lty=1,lwd=0.1)
#    }
    for (i in 1:nrow(graph.lines)) {
      pos <- graph.lines[i,]
      lines(root.dist[pos], label.pos[inv.dx[pos]])
    }
    mtext("GO-", side=1, line=0.5)
    mtext("Structure", side=1, line=1.5)
  }

  # draw the image
  if(show.graph) {
    par(las=1, mai=margins, fig=c(leftsplit,1.0,0.0,1.0), new=TRUE)
  } else {
    par(las=1, mai=margins)
  }
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

figure_ccimage <- function() {
  load("../results/brain_pomeroy/filtered.rdat")
  cc_image(res, pointsize=8, outfile="../plos/allsplits_pomeroy.eps")
#  system("epstopdf ../plos/allsplits_pomeroy.eps")

  # scatter plot for pairs
#  postscript("../plos/scatterpairs.eps", pointsize=8, width=7.0, height=6.0,
#             paper="special", horizontal=FALSE)
#  scatterpairs(res, "hgu133a")
#  dev.off()
#  system("epstopdf ../plos/scatterpairs.eps")
}

# Figure
# for induced graph for clinical correlation
figure_induced <- function() {
  load("../results/evalCorrelation.RData")
  goids <- as.character(rossResult$frame[1:7,1])
  gcc <- oneGOGraph(goids[1], GOCCPARENTS)
  gbp <- GOGraph(goids[2:7], GOBPPARENTS)
  g <- join(gbp, gcc)

  # write a .dot-file since I don't manage to use Rgraphviz
  cat("digraph {\ncenter=TRUE;\nrankdir=LR;\nsize=\"9,6\";\n",
      "node [style=filled,fontname=\"Helvetica\",fontsize=\"10\"];\n",
      sep="", file="../plos/induced.dot")
  for (n in nodes(g)) {
    curColor <- ifelse(n %in% goids, "yellow", "white")
    cat('"', n, '" [fillcolor=', curColor, "];\n",
        sep="", file="../plos/induced.dot", append=TRUE)
  }
  e <- edges(g)
  for (i in 1:length(e)) {
    for (n in e[[i]]) {
      cat('"', n, '" -> "', names(e)[i], "\"\n",
          sep="", file="../plos/induced.dot", append=TRUE)
    }
  }
  cat("}\n", file="../plos/induced.dot", append=TRUE)

  system("dot -Tps < ../plos/induced.dot > ../plos/induced.eps")

# replaced the following stuff with using dot instead of Rgraphviz
#  postscript("../plos/induced.eps", pointsize=8, width=6.0, height=6.0,
#             paper="special", horizontal=FALSE)
#  seedColors <- rep("yellow", length(goids))
#  names(seedColors) <- goids
#  plot(g,
#       attrs=list(graph=list(rankdir="LR"),
#                  node=list(shape="ellipse", fontsize="8"),
#                  edges=list(arrowhead="none", dir="back")),
#       nodeAttrs=list(fillcolor=seedColors))
#  dev.off()

  system("epstopdf ../plos/induced.eps")
  system("cp ../plos/induced.pdf /home/web/lottaz/tmp")
  system("chmod 644 /home/web/lottaz/tmp/induced.pdf")

  # write and compile a tex-file including the GO terms for g
  setwd("../plos")
  nds <- sort(nodes(g))
  nbucket <- round(length(nds)/3, 0)
  rest <- length(nds) %% nbucket  # to be added to the first minipage
  cat('\\documentclass{article}',
      '\\usepackage{graphicx}',
      '\\pagestyle{empty}',
      '\\textwidth=17cm',
      '\\begin{document}',
      '  \\begin{center}',
      '    \\includegraphics[width=12cm]{induced}',
      '    \\begin{tabular}{l@{}l@{}l}',
      '      \\begin{minipage}[b]{5cm} \\tiny',
      sep="\n", file="induced_go.tex")
  for (i in 1:(nbucket+rest)) {
    cat("       ", nds[i], gsub("_"," ",attr(get(nds[i], GOTERM), "Term")), '\\\\\n',
        sep=" ", file="induced_go.tex", append=TRUE)
  }
  cat('      \\end{minipage} &',
      '      \\begin{minipage}[b]{5cm} \\tiny ',

      sep="\n", file="induced_go.tex", append=TRUE)
  for (i in (nbucket+rest+1):(2*nbucket+rest)) {
    cat("       ", nds[i], gsub("_"," ",attr(get(nds[i], GOTERM), "Term")), '\\\\\n',
        sep=" ", file="induced_go.tex", append=TRUE)
  }
  cat('      \\end{minipage} &',
      '      \\begin{minipage}[b]{5cm} \\tiny ',
      sep="\n", file="induced_go.tex", append=TRUE)
  for (i in (2*nbucket+rest+1):length(nds)) {
    cat("       ", nds[i], gsub("_"," ",attr(get(nds[i], GOTERM), "Term")), '\\\\\n',
        sep=" ", file="induced_go.tex", append=TRUE)
  }
  cat('      \\end{minipage} \\\\',
      '    \\end{tabular}',
      '  \\end{center}',
      '\\end{document}',
      sep="\n", file="induced_go.tex", append=TRUE)
  system("latex induced_go")
  system("dvips induced_go")
  readline(paste("Use Ghostview interactively to convert induced_go.ps to induced_go.eps",
                 "then hit return...", sep="\n"))
  system("epstopdf induced_go.eps")

  setwd("../src")
}
