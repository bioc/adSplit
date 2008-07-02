adSplit <- function(mydata, annotation.ids, chip.name, min.probes=20, max.probes=NULL,
                    B=NULL, min.group.size=5, ngenes=50, ignore.genes=5) {

  # initialize
  if (is(mydata,"ExpressionSet")) mydata <- exprs(mydata)
  nsamples <- ncol(mydata)
  nprobesets <- nrow(mydata)
  if (is.null(max.probes)) max.probes <- round(nrow(mydata)/10, 0)
  collected.res <- list(cuts=NULL, scores=NULL, pvalues=NULL)
  class(collected.res) <- "splitSet"

  # get annotation ids if not provided directly
  require(package=paste(chip.name,".db",sep=""), character.only=TRUE)
  GOenv   <- eval(as.symbol(paste(chip.name, "GO2ALLPROBES", sep="")))
  KEGGenv <- eval(as.symbol(paste(chip.name, "PATH2PROBE", sep="")))
  chipProebeSets <- length(eval(as.symbol(paste(chip.name, "ACCNUM",sep=""))))
  if (chipProebeSets != nprobesets) 
    stop("Expected ", chipProebeSets, " probe-sets for ", chip.name, ", received ", nprobesets)
  if (length(annotation.ids) == 1) {
    if (annotation.ids == "GO") {
      annotation.ids <- ls(GOenv)
    } else {
      if (annotation.ids == "KEGG") {
        annotation.ids <- paste("KEGG:", Rkeys(KEGGenv), sep="")
      }  else {
        if (annotation.ids == "all")  
          annotation.ids <- c(ls(GOenv), paste("KEGG:", ls(KEGGenv), sep=""))
      }
    }
  }  
  for (id in annotation.ids) {
    # find probesets
    if (id == "all") next
    if (length(grep("GO:", id)) == 1) {
      probes <- GOenv[[id]]
    } else {
      if (length(grep("KEGG:", id)) == 1) {
        tmp.id <- sub("KEGG:","",id)
        probes <- KEGGenv[[tmp.id]]
      } else {
        stop("'", id, "' is not recognized as GO or KEGG identifier")
      }
    }
    probes <- probes[probes %in% rownames(mydata)]
    nprobes <- length(probes)

    # compute score for current id
    cat("Evaluating identifier", id, "with", nprobes, "probesets...\n")
    if (nprobes < min.probes) {
      cat(" -> skipped, too few probes associated (", nprobes, ")\n", sep="")
      next
    }
    if (nprobes > max.probes) {
      cat(" -> skipped, too many probes associated (", nprobes, ")\n", sep="")
      next
    }
    
    id.data <- mydata[probes,,drop=FALSE]
    res <- diana2means(id.data, min.group.size, ngenes, ignore.genes, return.cut=TRUE)

    if (!is.na(res$score)) {
      if (!is.null(B)) {
        # compute random scores and empirical p-value
        rand.scores <- randomDiana2means(nprobes, mydata, chip.name, 
                                         B, ngenes, ignore.genes)
        res$p <- sum(rand.scores > res$score)/length(rand.scores)
      }

      # collect results
      collected.res$cuts   <- rbind(collected.res$cuts, res$cut)
      collected.res$scores <- c(collected.res$scores, res$score)
      names(collected.res$scores)[length(collected.res$scores)] <- id
      if (!is.null(B)) {
        collected.res$pvalues <- c(collected.res$pvalues, res$p)
        names(collected.res$pvalues)[length(collected.res$pvalues)] <- id
      }
      rownames(collected.res$cuts)[nrow(collected.res$cuts)] <- id
      colnames(collected.res$cuts) <- colnames(mydata)
    }
  }

  # add q-values
  if (length(collected.res$pvalues) > 1) {
    q <- mt.rawp2adjp(collected.res$pvalues, "BH")
    collected.res$qvalues <- q$adjp[,2]
    collected.res$cuts <- collected.res$cuts[q$index,]
    collected.res$scores <- collected.res$scores[q$index]
    collected.res$pvalues <- collected.res$pvalues[q$index]
    names(collected.res$qvalues) <- names(collected.res$pvalues)
  }
  return(collected.res)
}

print.splitSet <- function(x, ...) {
  if (is.null(x$cuts)) {
    cat("Empty split set\n")
    return()
  }
  if (nrow(x$cuts) > 1) { 
    cat("Annotation-driven split set \n",
        " holds", nrow(x$cuts), "splits on", ncol(x$cuts), "elements\n",
        " scores range is:", range(x$scores), "\n")
    if (!is.null(x$pvalues)) {
      cat("  empirical p-values range is:", range(x$pvalues), "\n")
      cat("  q-value range is:", range(x$qvalues), "\n")
    } else {
      cat("  no empirical p-values computed\n")
    }
  } else { # just a single split
    id <- rownames(x$cuts)
    if (length(grep("GO:", id)) == 1) {
      term <- attr(GOTERM[[id]],"Term")
    } else {
      if (length(grep("KEGG:", id)) == 1) {
        tmp.id <- sub("KEGG:","",id)
        term <- KEGGPATHID2NAME[[tmp.id]]
      } else {
        stop("'", id, "' is not recognized as GO or KEGG identifier")
      }
    }
    term <- paste(term, " (", id, ")", sep="")
    cat("Annotation-driven split set \n",
        " holds 1 split on", ncol(x$cuts), "elements\n",
        " associated to:", term, "\n",
        " object distribution is:", table(x$cuts[1,]), "\n",
        " score is:", x$scores, "\n")
    if (!is.null(x$pvalues)) {
      cat("  empirical p-value is:", x$pvalues, "\n")
    } else {
      cat("  no empirical p-values computed\n")
    }
  }
}

hist.splitSet <- function(x, main="Distribution of p-Values", xlab="p-values", 
                          col="grey", xlim=c(0,1), ...) {
  mymai <- par("mai")
  mymai[4] <- mymai[1]
  par(mai=mymai)
  h <- hist(x$pvalues, main=main, xlab=xlab, col=col, xlim=xlim, ...)
  m <- max(h$counts)
  axis(side=4, at=m*seq(0.0, 1.0, 0.2), labels=seq(0.0, 1.0, 0.2))
  mtext("q-value", side=4, at=0.5*m, line=3, cex=1.0)
  lines(x$pvalues, m*x$qvalues)
}

makeEID2PROBESenv <- function(EIDenv) {
  # expects a hash with ENTREZIDs as values and Affy-probes as keys
  # generates a hash with ENTREZIDs as keys and Affy-probes as values
  env <- new.env(hash=TRUE)
  for (affyid in Lkeys(EIDenv)) {
    lid <- as.character(EIDenv[[affyid]])
    if (!is.null(lid)) assign(lid, affyid, env)
    else assign(lid, c(get(lid, env), affyid), env)
  }
  return(env)
}

drawRandomPS <- function(nps, EID2PSenv, allEIDs){
  # nps: number of probe set to draw at random
  # chip: name of the chip
  drawnEIDs <- sample(allEIDs,nps,replace=FALSE)
  drawnpsids <- unlist(mget(as.character(drawnEIDs), EID2PSenv))
  return(drawnpsids[1:nps])
} #drawRandomPS

randomDiana2means <- function(nprobes,data,chip,ndraws=10000, ngenes=50,ignore.genes=5){
  #nprobes: number of probesets to draw at random
  #data: rows=genes, columns=samples
  #ndraws: number of random draws
  #ngenes: number of genes to compute score on
  #ignore.genes: number of genes to disregard for score calculation
  cat("  determining", ndraws, "random DLD-scores with",nprobes,
      "probe sets each (wait for", round(ndraws/100,0), "dots)\n  ")
  require(package=paste(chip,".db",sep=""), character.only=TRUE)
  randscores <- numeric(ndraws) # initialize
  EIDenv <- eval(as.symbol(paste(chip, "ENTREZID", sep="")))
  EID2PSenv <- makeEID2PROBESenv(EIDenv)
  allEIDs <- ls(EID2PSenv)
  for (i in 1:ndraws){
    if(i%%100==0) cat(".")
    randps <- drawRandomPS(nprobes, EID2PSenv, allEIDs)
#    randps <- sample(1:nrow(data), nprobes, replace=FALSE)
    randdata <- data[randps,,drop=FALSE]
    randscores[i] <- diana2means(randdata, 1, ngenes,
                                 ignore.genes, return.cut=FALSE)
  } #for
  cat("\n")
  randscores[is.na(randscores)] <- 0
  return(randscores)
}#randomDiana2means
