

#' @title manhattan plot
#' @export
manhattan.plot <- function(x, p.threshold, ylim, ...) {
  
  c.names <- colnames(x)
  if (!("chr" %in%  c.names & "pos" %in%  c.names & "p.value" %in%  c.names)) 
    stop("Make sure your data frame contains columns 'chr', 'pos' and 'p.value'")   
  x <- x[na.omit(order(x$chr, x$pos)), ]
  if (any(x$p.value < 0 & x$p.value > 1))
    stop("P-values out of rage")
  if (any(x$p.value == 0)){
    warning("P-values of zero changed to 4.940656e-324")
    x$p.value[x$p.value == 0] <- 4.940656e-324
  }
  x$nlog.p <- -log10(x$p.value)
  if (missing(ylim))
    ylim <- c(0, max(x$nlog.p)+1)
  nChr <- max(x$chr)
  start.pos <- rep(NA_integer_, nChr)
  end.pos <- rep(NA_integer_, nChr)
  for(i in 1:nChr) {
    start.pos[i] <- min(x$pos[x$chr == i])
    end.pos[i] <- max(x$pos[x$chr == i])
  }
  chr.pos <- end.pos-start.pos
  layout(matrix(1:nChr, 1), widths=c(chr.pos/sum(chr.pos)))
  def.par <- par(no.readonly=TRUE, mai=c(.1, .1, .1, .1), oma=c(5, 6, 2, 2))
  for (i in 1:nChr) {
    plot(x$pos[x$chr == i], x$nlog.p[x$chr == i], bty='n',
         ylim=ylim, yaxs="i", yaxt="n", ylab="",
         xaxs="i", xaxt="n", xlab="",
         ...)
    axis(1, line=1, labels=FALSE, at=c(start.pos[i], end.pos[i]), 
         col="gray50", lwd=2)
    mtext(paste("Chr.", i), 1, line=3, col="gray50")
    if (i == 1L) {
       axis(2, line=1, col.axis="gray50", col="gray50", lwd=2, cex.axis=1.3)
       mtext(expression(-log[10](italic(p))), 2, line=4, col="gray50")
     }
    if (!missing(p.threshold)) {
      abline(h=-log10(p.threshold), lwd=1, col="gray50", xpd=TRUE)
    }
  }
  par(def.par)
  layout(1)
}



#' @title qq plot
#' @export
qq.plot <- function(x, ...) {
  if (!is.numeric(x$p.value)) 
    stop("P value vector is not numeric.")
  o <- -log10(sort(x$p.value[!is.na(x$p.value) & x$p.value < 1 & x$p.value > 0], 
                   decreasing=F))
  #e <- -log10(1:length(o)/length(o))
  e <- -log10(ppoints(length(o)))
  plot(e, o, pch=19, cex=1, xlab=expression(Expected~~-log[10](italic(p))), 
       ylab=expression(Observed~~-log[10](italic(p))), xlim=c(0,max(e)),
       ylim=c(0,max(o)), ...)
  abline(0, 1, lty=2)
}
