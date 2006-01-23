# determine split with diana and kmeans
diana2means <- function(mydata, mingroupsize=5, ngenes=50, ignore.genes=5, 
                        return.cut=FALSE) {

  if (class(mydata) == "exprSet") mydata <- exprs(mydata)
  nsamples <- ncol(mydata)
  na.cut <- rep(NA, nsamples)

  # 1. hierarchical clustering for initial cluster centroids
  mydiana <- diana(t(mydata),diss=FALSE,keep.diss=FALSE)
  mydianaclusters <- cutree(as.hclust(mydiana),k=2)
  mycent1 <- rowMeans(mydata[,mydianaclusters==1,drop=FALSE])
  mycent2 <- rowMeans(mydata[,mydianaclusters==2,drop=FALSE])
  mycents <- rbind(t(mycent1),t(mycent2))

  # Now use 2-means:
  mykmeans <- try(kmeans(t(mydata),mycents,100),silent=TRUE)
  if (class(mykmeans)=="try-error") {
    mycut <- na.cut
    myscore <- NA
  } else {
    mycut <- matrix(mykmeans$cluster-1,nrow=1)
    # check if cut is valid:
    if ((sum(mycut)<mingroupsize)|(sum(mycut)>(nsamples-mingroupsize))){
      mycut <- na.cut
      myscore <- NA
    } else {
      myscore <- try(tscore(mydata,mycut,p=ngenes,p.offs=ignore.genes),silent=TRUE)
      if (class(myscore)=="try-error") {mycut <- na.cut;myscore <- NA}
    } # else
  } # else
  if (return.cut) {
    tmp <- list(cut=mycut, score=myscore)
    class(tmp) <- "split"
    return(tmp)
  } else return(myscore)
} #diana2means

print.split <- function(x, ...) {
  cat("2-means split holding ", length(x$cut), " elements\n",
      "  distribution is: ", table(x$cut), "\n",
      "  DLD-score: ", x$score, "\n")
}

tscore <- function(dat, split, p = 50, p.offs = 0) {
  if (is.vector(split)) { ncs <- length(split); nrs <- 1          }  
  if (is.matrix(split)) { ncs <- ncol(split);   nrs <- nrow(split)}
  if (ncs != ncol(dat))
    stop(paste("Fatal error in tscore: number of columns of second argument \"split\" is", ncs,
               "but first argument \"dat\" has", ncol(dat), "columns.",
               "They should be the same.\n"))
  if(any(is.na(dat)))
     stop("Fatal error in tscore: first argument \"dat\" contains NA values.\n")
  p <- min(p, nrow(dat))
  
  means <- matrix(rep(     apply(dat, 1, mean), ncs), nrow=nrow(dat), ncol=ncs)
  stds  <- matrix(rep(sqrt(apply(dat, 1, var)), ncs), nrow=nrow(dat), ncol=ncs)
    
  res <- .C("isis", "tscore", 
            as.double(t((dat-means)/stds)),       
            as.integer(nrow(dat)),
            as.integer(ncol(dat)),
            as.integer(split),    
            as.integer(nrs),
            as.integer(c(p, p.offs, 0)),
            t = numeric(nrs),
            returncode = as.integer(0),       
            DUP = TRUE)
  
  if (res$returncode) stop(paste("Fatal error in tscore: caught exception", res$returncode))
  return (res$t)
}

